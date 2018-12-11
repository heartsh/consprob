extern crate stem;
extern crate scoped_threadpool;
extern crate bio;
extern crate num_cpus;

use stem::*;
use self::scoped_threadpool::Pool;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::fs::create_dir;
use std::thread;

type NumOfThreads = u32;

const DEFAULT_STA_FE_SCALE_PARAM: LogProb = 1. / INVERSE_TEMPERATURE;
const DEFAULT_MIN_BPP: Prob = 0.01;
const DEFAULT_GAP_NUM: usize = 0;
const DEFAULT_NUM_OF_TIMES_OF_EXPECT_MAX_ITER: usize = 1;
const BPP_MAT_ON_SS_FILE_NAME: &'static str = "bpp_mats_on_ss.dat";
const BAP_MAT_FILE_NAME: &'static str = "bap_mats.dat";
const OGP_MAT_FILE_NAME_1: &'static str = "ogp_mats_1.dat";
const OGP_MAT_FILE_NAME_2: &'static str = "ogp_mats_2.dat";
const EGP_MAT_FILE_NAME_1: &'static str = "egp_mats_1.dat";
const EGP_MAT_FILE_NAME_2: &'static str = "egp_mats_2.dat";
const BPAP_MAT_FILE_NAME: &'static str = "bpap_mats.dat";
const BPP_MAT_ON_STA_FILE_NAME: &'static str = "bpp_mats_on_sta.dat";
const VERSION: &'static str = "0.1.0";

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "struct_align_free_energy_scale_param", &format!("A structural-alignment free-energy scale parameter (Uses {} by default)", DEFAULT_STA_FE_SCALE_PARAM), "FLOAT");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "gap_num", &format!("A gap number for setting the maximum number of the gaps in a structural alignment; This gap number plus the absolute value of the difference of 2 RNA sequence lengths is this maximum number (Uses {} by default)", DEFAULT_GAP_NUM), "UINT");
  opts.optopt("", "num_of_times_of_expect_max_iter", &format!("The number of times of an expectation-maximization iteration (Uses {} by default)", DEFAULT_NUM_OF_TIMES_OF_EXPECT_MAX_ITER), "UINT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_file_path = opts.opt_str("i").expect("Failed to get the path to an input FASTA file containing RNA sequences from command arguments.");
  let input_file_path = Path::new(&input_file_path);
  let output_dir_path = opts.opt_str("o").expect("Failed to get the path to an output directory from command arguments.");
  let output_dir_path = Path::new(&output_dir_path);
  let sta_fe_scale_param = if opts.opt_present("struct_align_free_energy_scale_param") {
    opts.opt_str("struct_align_free_energy_scale_param").expect("Failed to get a structural-alignment free-energy scale parameter from command arguments.").parse().expect("Failed to parse a structural-alignment free-energy scale parameter.")
  } else {
    DEFAULT_STA_FE_SCALE_PARAM
  };
  let min_bpp = if opts.opt_present("min_base_pair_prob") {
    opts.opt_str("min_base_pair_prob").expect("Failed to get a minimum base-pairing-probability from command arguments.").parse().expect("Failed to parse a minimum base-pairing-probability.")
  } else {
    DEFAULT_MIN_BPP
  };
  let min_lbpp = min_bpp.ln();
  let gap_num = if opts.opt_present("gap_num") {
    opts.opt_str("gap_num").expect("Failed to get a gap number for setting the maximum number of the gaps in a structural alignment from command arguments.").parse().expect("Failed to parse a gap number for setting the maximum number of the gaps in a structural alignment.")
  } else {
    DEFAULT_GAP_NUM
  };
  let num_of_times_of_expect_max_iter = if opts.opt_present("num_of_times_of_expect_max_iter") {
    opts.opt_str("num_of_times_of_expect_max_iter").expect("Failed to get the number of times of an expectation-maximization iteration from command arguments.").parse().expect("Failed to parse the number of times of an expectation-maximization iteration.")
  } else {
    DEFAULT_NUM_OF_TIMES_OF_EXPECT_MAX_ITER
  };
  let num_of_threads = if opts.opt_present("t") {
    opts.opt_str("t").expect("Failed to get the number of threads in multithreading from command arguments.").parse().expect("Failed to parse the number of threads in multithreading.")
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
  let mut lbpp_mats = vec![SparseLogProbMat::default(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (lbpp_mat, fasta_record) in lbpp_mats.iter_mut().zip(fasta_records.iter_mut()) {
      let seq_len = fasta_record.seq.len();
      scope.execute(move || {
        *lbpp_mat = remove_little_lbpps_from_lbpp_mat(&get_log_bpp_mat(&fasta_record.seq[1 .. seq_len - 1]), min_lbpp);
        lbpp_mat.insert((0, seq_len - 1), 0.);
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let mut buf_4_writer_2_bpp_mat_on_ss_file = format!("; The version {} of the STEM program.\n; The path to the input file for computing the Base-Pairing Probability Matrices (= BPPMs) on secondary structure (= SS) in this file = \"{}\".\n; The values of the parameters used for computing these matrices are as follows.\n; \"min_base_pair_prob\" = {}, \"num_of_threads\" = {}.", VERSION, input_file_path.display(), min_bpp, num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the BPPM of this sequence on SS.";
  let bpp_mat_on_ss_file_path = output_dir_path.join(BPP_MAT_ON_SS_FILE_NAME);
  let mut writer_2_bpp_mat_on_ss_file = BufWriter::new(File::create(bpp_mat_on_ss_file_path).expect("Failed to create an output file."));
  for (rna_id, lbpp_mat) in lbpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &lbpp) in lbpp_mat.iter() {
      if i == 0 {continue;}
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, lbpp.exp()));
    }
    buf_4_writer_2_bpp_mat_on_ss_file.push_str(&buf_4_rna_id);
  }
  thread::spawn(move || {
    let _ = writer_2_bpp_mat_on_ss_file.write_all(buf_4_writer_2_bpp_mat_on_ss_file.as_bytes());
  });
  let mut sta_fe_param_sets_with_rna_id_pairs = StaFeParamSetsWithRnaIdPairs::default();
  let mut lstapmts_with_rna_id_pairs = LstapmtsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      sta_fe_param_sets_with_rna_id_pairs.insert(rna_id_pair, StaFeParams::origin());
      lstapmts_with_rna_id_pairs.insert(rna_id_pair, LogStapmt::origin());
    }
  }
  let mut lstapmts_with_rna_id_pairs_from_pct = lstapmts_with_rna_id_pairs.clone();
  thread_pool.scoped(|scope| {
    for (rna_id_pair, sta_fe_params) in sta_fe_param_sets_with_rna_id_pairs.iter_mut() {
      let seq_len_pair = (fasta_records[rna_id_pair.0].seq.len(), fasta_records[rna_id_pair.1].seq.len());
      let max_gap_num = gap_num + get_seq_len_diff(&seq_len_pair);
      let ref ref_2_fasta_records = fasta_records;
      let ref ref_2_lbpp_mats = lbpp_mats;
      scope.execute(move || {
        *sta_fe_params = StaFeParams::new(rna_id_pair, ref_2_fasta_records, max_gap_num, ref_2_lbpp_mats);
      });
    }
  });
  for i in 0 .. num_of_times_of_expect_max_iter {
    thread_pool.scoped(|scope| {
      for (rna_id_pair, lstapmt) in lstapmts_with_rna_id_pairs.iter_mut() {
        let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
        let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
        let max_gap_num = gap_num + get_seq_len_diff(&seq_len_pair);
        let ref sta_fe_params = sta_fe_param_sets_with_rna_id_pairs[&rna_id_pair];
        scope.execute(move || {
          *lstapmt = io_algo_4_rna_lstapmt(&seq_pair, &seq_len_pair, sta_fe_params, sta_fe_scale_param, max_gap_num);
        });
      }
    });
    thread_pool.scoped(|scope| {
      for (rna_id_pair, lstapmt) in lstapmts_with_rna_id_pairs_from_pct.iter_mut() {
        let ref ref_2_lstapmts_with_rna_id_pairs = lstapmts_with_rna_id_pairs;
        scope.execute(move || {
          *lstapmt = prob_cons_transformation_of_lstapmt(ref_2_lstapmts_with_rna_id_pairs, rna_id_pair, num_of_fasta_records);
        });
      }
    });
    if i == num_of_times_of_expect_max_iter - 1 {continue;}
    for (rna_id_pair, lstapmt) in lstapmts_with_rna_id_pairs.iter_mut() {
      let new_lstapmt = lstapmts_with_rna_id_pairs_from_pct[rna_id_pair].clone();
      *lstapmt = new_lstapmt.clone();
      sta_fe_param_sets_with_rna_id_pairs.get_mut(rna_id_pair).expect("Failed to get an element from a hash map.").lstapmt = new_lstapmt;
    }
  }
  thread_pool.scoped(|scope| {
    for (rna_id, lbpp_mat) in lbpp_mats.iter_mut().enumerate() {
      let ref ref_2_lstapmts_with_rna_id_pairs = lstapmts_with_rna_id_pairs;
      scope.execute(move || {
        *lbpp_mat = pct_of_lbpp_mat(ref_2_lstapmts_with_rna_id_pairs, rna_id, num_of_fasta_records, lbpp_mat);
      });
    }
  });
  let stapmts_with_rna_id_pairs = get_stapmts_with_rna_id_pairs(&lstapmts_with_rna_id_pairs_from_pct);
  let output_file_header = format!(" in this file = \"{}\".\n; The values of the parameters used for computing these matrices are as follows.\n; \"struct_align_free_energy_scale_param\" = {}, \"min_base_pair_prob\" = {}, \"gap_num\" = {}, \"num_of_times_of_expect_max_iter\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), sta_fe_scale_param, min_bpp, gap_num, num_of_times_of_expect_max_iter, num_of_threads);
  let bap_mat_file_path = output_dir_path.join(BAP_MAT_FILE_NAME);
  let ogp_mat_file_path_1 = output_dir_path.join(OGP_MAT_FILE_NAME_1);
  let ogp_mat_file_path_2 = output_dir_path.join(OGP_MAT_FILE_NAME_2);
  let egp_mat_file_path_1 = output_dir_path.join(EGP_MAT_FILE_NAME_1);
  let egp_mat_file_path_2 = output_dir_path.join(EGP_MAT_FILE_NAME_2);
  let bpap_mat_file_path = output_dir_path.join(BPAP_MAT_FILE_NAME);
  let bpp_mat_on_sta_file_path = output_dir_path.join(BPP_MAT_ON_STA_FILE_NAME);
  let mut writer_2_bap_mat_file = BufWriter::new(File::create(bap_mat_file_path).expect("Failed to create an output file."));
  let mut writer_2_ogp_mat_file_1 = BufWriter::new(File::create(ogp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_ogp_mat_file_2 = BufWriter::new(File::create(ogp_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_egp_mat_file_1 = BufWriter::new(File::create(egp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_egp_mat_file_2 = BufWriter::new(File::create(egp_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_bpap_mat_file = BufWriter::new(File::create(bpap_mat_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_bap_mat_file = format!("; The version {} of the STEM program.\n; The path to the input file for computing base alignment probability matrices (= BAPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the BAPM between these 2 sequences.";
  let mut buf_4_writer_2_ogp_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing opening gap probability matrices (= OGPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the OGPM between these 2 sequences.";
  let mut buf_4_writer_2_ogp_mat_file_2 = buf_4_writer_2_ogp_mat_file_1.clone();
  let mut buf_4_writer_2_egp_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing extending gap probability matrices (= EGPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the EGPM between these 2 sequences.";
  let mut buf_4_writer_2_egp_mat_file_2 = buf_4_writer_2_egp_mat_file_1.clone();
  let mut buf_4_writer_2_bpap_mat_file = format!("; The version {} of the STEM program.\n; The path to the input file for computing base pair alignment probability matrices (= BPAPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the BPAPM between these 2 sequences.";
  for (rna_id_pair, stapmt) in &stapmts_with_rna_id_pairs {
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    let seq_len_pair = (fasta_records[rna_id_pair.0].seq.len(), fasta_records[rna_id_pair.1].seq.len());
    for (&(i, j), bap) in stapmt.base_align_prob_mat.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, bap));
    }
    buf_4_writer_2_bap_mat_file.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 1 .. seq_len_pair.0 - 1 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i - 1, stapmt.opening_gap_prob_mat_1[i]));
    }
    buf_4_writer_2_ogp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 1 .. seq_len_pair.1 - 1 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i - 1, stapmt.ogp_mat_2[i]));
    }
    buf_4_writer_2_ogp_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 1 .. seq_len_pair.0 - 1 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i - 1, stapmt.extending_gap_prob_mat_1[i]));
    }
    buf_4_writer_2_egp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 1 .. seq_len_pair.1 - 1 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i - 1, stapmt.egp_mat_2[i]));
    }
    buf_4_writer_2_egp_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k, l), &bpap) in stapmt.base_pair_align_prob_mat.iter() {
      if i == 0 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{},{} ", i - 1, j - 1, k - 1, l - 1, bpap));
    }
    buf_4_writer_2_bpap_mat_file.push_str(&buf_4_rna_id_pair);
  }
  thread::spawn(move || {
    let _ = writer_2_bap_mat_file.write_all(buf_4_writer_2_bap_mat_file.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_ogp_mat_file_1.write_all(buf_4_writer_2_ogp_mat_file_1.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_ogp_mat_file_2.write_all(buf_4_writer_2_ogp_mat_file_2.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_egp_mat_file_1.write_all(buf_4_writer_2_egp_mat_file_1.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_egp_mat_file_2.write_all(buf_4_writer_2_egp_mat_file_2.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_bpap_mat_file.write_all(buf_4_writer_2_bpap_mat_file.as_bytes());
  });
  let mut buf_4_writer_2_bpp_mat_on_sta_file = format!("; The version {} of the STEM program.\n; The path to the input file for computing the Base-Pairing Probability Matrices (= BPPMs) on STructural Alignment (= STA) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the BPPM of this sequence on STA.";
  let mut writer_2_bpp_mat_on_sta_file = BufWriter::new(File::create(bpp_mat_on_sta_file_path).expect("Failed to create an output file."));
  for (rna_id, lbpp_mat) in lbpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &lbpp) in lbpp_mat.iter() {
      if i == 0 {continue;}
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, lbpp.exp()));
    }
    buf_4_writer_2_bpp_mat_on_sta_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_on_sta_file.write_all(buf_4_writer_2_bpp_mat_on_sta_file.as_bytes());
}
