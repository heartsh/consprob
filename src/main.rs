extern crate rnafamprob;
extern crate scoped_threadpool;
extern crate bio;
extern crate num_cpus;
extern crate itertools;

use rfamprob::*;
use self::scoped_threadpool::Pool;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::fs::create_dir;
use std::thread;
use itertools::multizip;

type NumOfThreads = u32;

const DEFAULT_OPENING_GAP_PENALTY: StaFreeEnergy = 0.;
const DEFAULT_EXTENDING_GAP_PENALTY: StaFreeEnergy = 0.;
const DEFAULT_STA_FE_SCALE_PARAM: LogProb = 1. / INVERSE_TEMPERATURE;
const DEFAULT_MIN_BPP: Prob = 0.005;
const DEFAULT_GAP_NUM: usize = 0;
const BPP_MAT_ON_SS_FILE_NAME: &'static str = "bpp_mats_on_ss.dat";
const BAP_MAT_FILE_NAME: &'static str = "bap_mats.dat";
const BPAP_MAT_FILE_NAME: &'static str = "bpap_mats.dat";
const BPP_MAT_ON_STA_FILE_NAME: &'static str = "bpp_mats_on_sta.dat";
const UPP_MAT_FILE_NAME: &'static str = "upp_mats.dat";
const VERSION: &'static str = "0.1.0";

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "opening_gap_penalty", &format!("An opening-gap penalty (Uses {} by default)", DEFAULT_OPENING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "extending_gap_penalty", &format!("An extending-gap penalty (Uses {} by default)", DEFAULT_EXTENDING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "struct_align_free_energy_scale_param", &format!("A structural-alignment free-energy scale parameter (Uses {} by default)", DEFAULT_STA_FE_SCALE_PARAM), "FLOAT");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "gap_num", &format!("A gap number for setting the maximum number of the gaps in a structural alignment; This gap number plus the absolute value of the difference of 2 RNA sequence lengths is this maximum number (Uses {} by default)", DEFAULT_GAP_NUM), "UINT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").expect("Failed to get the path to an input FASTA file containing RNA sequences from command arguments.");
  let input_file_path = Path::new(&input_file_path);
  let output_dir_path = matches.opt_str("o").expect("Failed to get the path to an output directory from command arguments.");
  let output_dir_path = Path::new(&output_dir_path);
  let opening_gap_penalty = if matches.opt_present("opening_gap_penalty") {
    matches.opt_str("opening_gap_penalty").expect("Failed to get an opening-gap penalty from command arguments.").parse().expect("Failed to parse an opening-gap penalty.")
  } else {
    DEFAULT_OPENING_GAP_PENALTY
  };
  let extending_gap_penalty = if matches.opt_present("extending_gap_penalty") {
    matches.opt_str("extending_gap_penalty").expect("Failed to get an extending-gap penalty from command arguments.").parse().expect("Failed to parse an extending-gap penalty.")
  } else {
    DEFAULT_EXTENDING_GAP_PENALTY
  };
  let sta_fe_scale_param = if matches.opt_present("struct_align_free_energy_scale_param") {
    matches.opt_str("struct_align_free_energy_scale_param").expect("Failed to get a structural-alignment free-energy scale parameter from command arguments.").parse().expect("Failed to parse a structural-alignment free-energy scale parameter.")
  } else {
    DEFAULT_STA_FE_SCALE_PARAM
  };
  let min_bpp = if matches.opt_present("min_base_pair_prob") {
    matches.opt_str("min_base_pair_prob").expect("Failed to get a minimum base-pairing-probability from command arguments.").parse().expect("Failed to parse a minimum base-pairing-probability.")
  } else {
    DEFAULT_MIN_BPP
  };
  let gap_num = if matches.opt_present("gap_num") {
    matches.opt_str("gap_num").expect("Failed to get a gap number for setting the maximum number of the gaps in a structural alignment from command arguments.").parse().expect("Failed to parse a gap number for setting the maximum number of the gaps in a structural alignment.")
  } else {
    DEFAULT_GAP_NUM
  };
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").expect("Failed to get the number of threads in multithreading from command arguments.").parse().expect("Failed to parse the number of threads in multithreading.")
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).expect("Failed to set a FASTA file reader.");
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.expect("Failed to read a FASTA record.");
    let mut seq: Seq = unsafe {from_utf8_unchecked(fasta_record.seq()).to_uppercase().as_bytes().iter().filter(|&&base| {is_rna_base(base)}).map(|&base| {base}).collect()};
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let num_of_fasta_records = fasta_records.len();
  let mut thread_pool = Pool::new(num_of_threads);
  let mut sparse_bpp_mats_1 = vec![SparseProbMat::default(); num_of_fasta_records];
  let mut sparse_bpp_mats_2 = sparse_bpp_mats_1.clone();
  let mut upp_mats = vec![Probs::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (sparse_bpp_mat_1, sparse_bpp_mat_2, upp_mat, fasta_record) in multizip((sparse_bpp_mats_1.iter_mut(), sparse_bpp_mats_2.iter_mut(), upp_mats.iter_mut(), fasta_records.iter_mut())) {
      let seq_len = fasta_record.seq.len();
      scope.execute(move || {
        let prob_mat_pair = get_bpp_and_unpair_prob_mats(&fasta_record.seq[1 .. seq_len - 1]);
        *sparse_bpp_mat_1 = remove_zeros_from_bpp_mat(&prob_mat_pair.0, seq_len);
        *sparse_bpp_mat_2 = remove_small_bpps_from_bpp_mat(&sparse_bpp_mat_1, min_bpp);
        *upp_mat = prob_mat_pair.1;
        upp_mat.insert(0, 0.);
        upp_mat.push(0.);
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let mut buf_4_writer_2_bpp_mat_on_ss_file = format!("; The version {} of the RNAfamProb program.\n; The path to the input file for computing the Base-Pairing Probability Matrices (= BPPMs) on secondary structure (= SS) in this file = \"{}\".\n; The values of the parameters used for computing these matrices are as follows.\n; \"num_of_threads\" = {}.", VERSION, input_file_path.display(), num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the BPPM of this sequence on SS.";
  let bpp_mat_on_ss_file_path = output_dir_path.join(BPP_MAT_ON_SS_FILE_NAME);
  let mut writer_2_bpp_mat_on_ss_file = BufWriter::new(File::create(bpp_mat_on_ss_file_path).expect("Failed to create an output file."));
  for (rna_id, bpp_mat) in sparse_bpp_mats_1.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in bpp_mat.iter() {
      if i == 0 {continue;}
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
    }
    buf_4_writer_2_bpp_mat_on_ss_file.push_str(&buf_4_rna_id);
  }
  let thread_1 = thread::spawn(move || {
    let _ = writer_2_bpp_mat_on_ss_file.write_all(buf_4_writer_2_bpp_mat_on_ss_file.as_bytes());
  });
  let mut sta_fe_param_sets_with_rna_id_pairs = StaFeParamSetsWithRnaIdPairs::default();
  let mut log_prob_mat_pairs_on_sta_with_rna_id_pairs = LogProbMatPairsOnStaWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      sta_fe_param_sets_with_rna_id_pairs.insert(rna_id_pair, StaFeParams::origin());
      log_prob_mat_pairs_on_sta_with_rna_id_pairs.insert(rna_id_pair, LogProbMatPairOnSta::new());
    }
  }
  let mut log_prob_mat_pairs_on_sta_with_rna_id_pairs_from_pct = log_prob_mat_pairs_on_sta_with_rna_id_pairs.clone();
  thread_pool.scoped(|scope| {
    for (rna_id_pair, sta_fe_params) in sta_fe_param_sets_with_rna_id_pairs.iter_mut() {
      let seq_len_pair = (fasta_records[rna_id_pair.0].seq.len(), fasta_records[rna_id_pair.1].seq.len());
      let max_gap_num = gap_num + get_seq_len_diff(&seq_len_pair);
      let ref ref_2_fasta_records = fasta_records;
      let ref ref_2_bpp_mats = sparse_bpp_mats_2;
      scope.execute(move || {
        *sta_fe_params = StaFeParams::new(rna_id_pair, ref_2_fasta_records, max_gap_num, ref_2_bpp_mats, sta_fe_scale_param, opening_gap_penalty, extending_gap_penalty);
      });
    }
  });
  thread_pool.scoped(|scope| {
    for (rna_id_pair, log_prob_mat_pair_on_sta) in log_prob_mat_pairs_on_sta_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
      let max_gap_num = gap_num + get_seq_len_diff(&seq_len_pair);
      let ref sta_fe_params = sta_fe_param_sets_with_rna_id_pairs[&rna_id_pair];
      scope.execute(move || {
        *log_prob_mat_pair_on_sta = io_algo_4_log_prob_mat_pair_on_rna_sta(&seq_pair, &seq_len_pair, sta_fe_params, max_gap_num);
      });
    }
  });
  thread_pool.scoped(|scope| {
    for (rna_id_pair, log_prob_mat_pair_on_sta) in log_prob_mat_pairs_on_sta_with_rna_id_pairs_from_pct.iter_mut() {
      let ref ref_2_log_prob_mat_pairs_on_sta_with_rna_id_pairs = log_prob_mat_pairs_on_sta_with_rna_id_pairs;
      scope.execute(move || {
        *log_prob_mat_pair_on_sta = prob_cons_transformation_of_log_prob_mat_pair_on_sta(ref_2_log_prob_mat_pairs_on_sta_with_rna_id_pairs, rna_id_pair, num_of_fasta_records);
      });
    }
  });
  thread_pool.scoped(|scope| {
    for (rna_id, bpp_mat, upp_mat) in multizip((0 .. num_of_fasta_records, sparse_bpp_mats_1.iter_mut(), upp_mats.iter_mut())) {
      let ref ref_2_log_prob_mat_pairs_on_sta_with_rna_id_pairs = log_prob_mat_pairs_on_sta_with_rna_id_pairs;
      scope.execute(move || {
        let prob_mat_pair = pct_of_bpp_and_upp_mat(ref_2_log_prob_mat_pairs_on_sta_with_rna_id_pairs, rna_id, num_of_fasta_records, bpp_mat, upp_mat);
        *bpp_mat = prob_mat_pair.0;
        *upp_mat = prob_mat_pair.1;
      });
    }
  });
  let prob_mat_pairs_on_sta_with_rna_id_pairs = get_prob_mat_pairs_on_sta_with_rna_id_pairs(&log_prob_mat_pairs_on_sta_with_rna_id_pairs_from_pct);
  let output_file_header = format!(" in this file = \"{}\".\n; The values of the parameters used for computing these matrices are as follows.\n; \"opening_gap_penalty\" = {}, \"extending_gap_penalty\" = {}, \"struct_align_free_energy_scale_param\" = {}, \"min_base_pair_prob\" = {}, \"gap_num\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), opening_gap_penalty, extending_gap_penalty, sta_fe_scale_param, min_bpp, gap_num, num_of_threads);
  let bap_mat_file_path = output_dir_path.join(BAP_MAT_FILE_NAME);
  let bpap_mat_file_path = output_dir_path.join(BPAP_MAT_FILE_NAME);
  let bpp_mat_on_sta_file_path = output_dir_path.join(BPP_MAT_ON_STA_FILE_NAME);
  let upp_mat_file_path = output_dir_path.join(UPP_MAT_FILE_NAME);
  let mut writer_2_bap_mat_file = BufWriter::new(File::create(bap_mat_file_path).expect("Failed to create an output file."));
  let mut writer_2_bpap_mat_file = BufWriter::new(File::create(bpap_mat_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_bap_mat_file = format!("; The version {} of the RNAfamProb program.\n; The path to the input file for computing base alignment probability matrices (= BAPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the BAPM between these 2 sequences.";
  let mut buf_4_writer_2_bpap_mat_file = format!("; The version {} of the RNAfamProb program.\n; The path to the input file for computing base pair alignment probability matrices (= BPAPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the BPAPM between these 2 sequences.";
  for (rna_id_pair, prob_mat_pair_on_sta) in &prob_mat_pairs_on_sta_with_rna_id_pairs {
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j), bap) in prob_mat_pair_on_sta.base_align_prob_mat.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, bap));
    }
    buf_4_writer_2_bap_mat_file.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k, l), &bpap) in prob_mat_pair_on_sta.base_pair_align_prob_mat.iter() {
      if i == 0 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{},{} ", i - 1, j - 1, k - 1, l - 1, bpap));
    }
    buf_4_writer_2_bpap_mat_file.push_str(&buf_4_rna_id_pair);
  }
  let thread_2 = thread::spawn(move || {
    let _ = writer_2_bap_mat_file.write_all(buf_4_writer_2_bap_mat_file.as_bytes());
  });
  let thread_3 = thread::spawn(move || {
    let _ = writer_2_bpap_mat_file.write_all(buf_4_writer_2_bpap_mat_file.as_bytes());
  });
  let mut writer_2_bpp_mat_on_sta_file = BufWriter::new(File::create(bpp_mat_on_sta_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_bpp_mat_on_sta_file = format!("; The version {} of the RNAfamProb program.\n; The path to the input file for computing the Base-Pairing Probability Matrices (= BPPMs) on STructural Alignment (= STA) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the BPPM of this sequence on STA.";
  for (rna_id, bpp_mat) in sparse_bpp_mats_1.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in bpp_mat.iter() {
      if i == 0 {continue;}
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
    }
    buf_4_writer_2_bpp_mat_on_sta_file.push_str(&buf_4_rna_id);
  }
  let thread_4 = thread::spawn(move || {
    let _ = writer_2_bpp_mat_on_sta_file.write_all(buf_4_writer_2_bpp_mat_on_sta_file.as_bytes());
  });
  let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_upp_mat_file = format!("; The version {} of the RNAfamProb program.\n; The path to the input file for computing the unpairing probability matrices in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the BPPM of this sequence on STA.";
  for (rna_id, upp_mat) in upp_mats.iter().enumerate() {
    let seq_len = fasta_records[rna_id].seq.len();
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (i, &upp) in upp_mat.iter().enumerate() {
      if i == 0 || i == seq_len - 1 {continue;}
      buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
    }
    buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
  let _ = thread_1.join();
  let _ = thread_2.join();
  let _ = thread_3.join();
  let _ = thread_4.join();
}
