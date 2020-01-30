extern crate phyloprob;
extern crate scoped_threadpool;
extern crate bio;
extern crate num_cpus;
extern crate itertools;

use phyloprob::*;
use self::scoped_threadpool::Pool;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::fs::create_dir;
use itertools::multizip;

type NumOfThreads = u32;

const DEFAULT_OPENING_GAP_PENALTY: FreeEnergy = 0.;
const DEFAULT_EXTENDING_GAP_PENALTY: FreeEnergy = 0.;
const DEFAULT_MIN_BPP: Prob = 0.005;
const DEFAULT_MAX_POS_DIST_4_IL: usize = 5;
const DEFAULT_OFFSET_4_MAX_POS_DIST_ON_EL: usize = 0;
const DEFAULT_MAX_SUBSTR_DIST: usize = 5;
const BPP_MAT_ON_SS_FILE_NAME: &'static str = "bpp_mats_on_ss.dat";
const BPP_MAT_ON_STA_FILE_NAME: &'static str = "bpp_mats_on_sta.dat";
const UPP_MAT_ON_SS_FILE_NAME: &'static str = "upp_mats_on_ss.dat";
const UPP_MAT_ON_STA_FILE_NAME: &'static str = "upp_mats_on_sta.dat";
const VERSION: &'static str = "0.1.0";

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "opening_gap_penalty", &format!("An opening-gap penalty (Uses {} by default)", DEFAULT_OPENING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "extending_gap_penalty", &format!("An extending-gap penalty (Uses {} by default)", DEFAULT_EXTENDING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "max_pos_dist_4_int_loop", &format!("A maximum position distance for internal loop (Uses {} by default)", DEFAULT_MAX_POS_DIST_4_IL), "UINT");
  opts.optopt("", "offset_4_max_pos_dist_on_ext_loop", &format!("An offset for maximum position distances for external loop (Uses {} by default)", DEFAULT_OFFSET_4_MAX_POS_DIST_ON_EL), "UINT");
  opts.optopt("", "max_substr_dist", &format!("A maximum substring distance (Uses {} by default)", DEFAULT_MAX_SUBSTR_DIST), "UINT");
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
  }.exp();
  let exp_opening_gap_penalty = opening_gap_penalty.exp();
  let extending_gap_penalty = if matches.opt_present("extending_gap_penalty") {
    matches.opt_str("extending_gap_penalty").expect("Failed to get an extending-gap penalty from command arguments.").parse().expect("Failed to parse an extending-gap penalty.")
  } else {
    DEFAULT_EXTENDING_GAP_PENALTY
  };
  let exp_extending_gap_penalty = extending_gap_penalty.exp();
  let min_bpp = if matches.opt_present("min_base_pair_prob") {
    matches.opt_str("min_base_pair_prob").expect("Failed to get a minimum base-pairing-probability from command arguments.").parse().expect("Failed to parse a minimum base-pairing-probability.")
  } else {
    DEFAULT_MIN_BPP
  };
  let max_pos_dist_4_il = if matches.opt_present("max_pos_dist_4_int_loop") {
    matches.opt_str("max_pos_dist_4_int_loop").expect("Failed to get a maximum position distance for internal loop from command arguments.").parse().expect("Failed to parse a maximum position distance for internal loop.")
  } else {
    DEFAULT_MAX_POS_DIST_4_IL
  };
  let offset_4_max_pos_dist_on_el = if matches.opt_present("offset_4_max_pos_dist_on_ext_loop") {
    matches.opt_str("offset_4_max_pos_dist_on_ext_loop").expect("Failed to get an offset for maximum position distances for external loop from command arguments.").parse().expect("Failed to parse an offset for maximum position distances for external loop.")
  } else {
    DEFAULT_OFFSET_4_MAX_POS_DIST_ON_EL
  };
  let max_substr_dist = if matches.opt_present("max_substr_dist") {
    matches.opt_str("max_substr_dist").expect("Failed to get a maximum substring distance from command arguments.").parse().expect("Failed to parse a maximum substring distance.")
  } else {
    DEFAULT_MAX_SUBSTR_DIST
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
  let mut bpp_mats = vec![SparseProbMat::default(); num_of_fasta_records];
  let mut sparse_bpp_mats = bpp_mats.clone();
  let mut upp_mats = vec![Probs::new(); num_of_fasta_records];
  let mut max_bp_spans = vec![0; num_of_fasta_records];
  let mut invert_exp_max_free_energies = vec![0.; num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (bpp_mat, sparse_bpp_mat, upp_mat, max_bp_span, fasta_record, invert_exp_max_free_energy) in multizip((bpp_mats.iter_mut(), sparse_bpp_mats.iter_mut(), upp_mats.iter_mut(), max_bp_spans.iter_mut(), fasta_records.iter_mut(), invert_exp_max_free_energies.iter_mut())) {
      let seq_len = fasta_record.seq.len();
      scope.execute(move || {
        let (obtained_bpp_mat, obtained_upp_mat, max_free_energy) = get_bpp_and_unpair_prob_mats(&fasta_record.seq[1 .. seq_len - 1]);
        *bpp_mat = obtained_bpp_mat;
        *sparse_bpp_mat = remove_small_bpps_from_bpp_mat(&bpp_mat, min_bpp);
        *upp_mat = obtained_upp_mat;
        *invert_exp_max_free_energy = 1. / max_free_energy.exp();
        *max_bp_span = get_max_bp_span(sparse_bpp_mat);
        upp_mat.insert(0, 0.);
        upp_mat.push(0.);
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let mut buf_4_writer_2_bpp_mat_on_ss_file = format!("; The version {} of the Zprob program.\n; The path to the input file in order to compute the base-pairing probability matrices on secondary structure in this file = \"{}\".\n; The values of the parameters used in order to compute the matrices are as follows.\n; \"num_of_threads\" = {}.", VERSION, input_file_path.display(), num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the base-pairing probability marix on secondary structure of the sequence.";
  let bpp_mat_on_ss_file_path = output_dir_path.join(BPP_MAT_ON_SS_FILE_NAME);
  let mut writer_2_bpp_mat_on_ss_file = BufWriter::new(File::create(bpp_mat_on_ss_file_path).expect("Failed to create an output file."));
  for (rna_id, bpp_mat) in bpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in bpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
    }
    buf_4_writer_2_bpp_mat_on_ss_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_on_ss_file.write_all(buf_4_writer_2_bpp_mat_on_ss_file.as_bytes());
  let upp_mat_on_ss_file_path = output_dir_path.join(UPP_MAT_ON_SS_FILE_NAME);
  let mut writer_2_upp_mat_on_ss_file = BufWriter::new(File::create(upp_mat_on_ss_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_upp_mat_on_ss_file = format!("; The version {} of the Zprob program.\n; The path to the input file in order to compute the unpairing probability matrices on secondary structure in this file = \"{}\".\n; The values of the parameters used in order to compute the matrices are as follows.\n; \"num_of_threads\" = {}.", VERSION, input_file_path.display(), num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the unpairing probability marix on secondary structure of the sequence.";
  for (rna_id, upp_mat) in upp_mats.iter().enumerate() {
    let seq_len = fasta_records[rna_id].seq.len();
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (i, &upp) in upp_mat.iter().enumerate() {
      if i == 0 || i == seq_len - 1 {continue;}
      buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
    }
    buf_4_writer_2_upp_mat_on_ss_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_upp_mat_on_ss_file.write_all(buf_4_writer_2_upp_mat_on_ss_file.as_bytes());
  let mut sta_fe_param_sets_with_rna_id_pairs = StaFeParamSetsWithRnaIdPairs::default();
  let mut bpap_mats_with_rna_id_pairs = Prob4dMatsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      sta_fe_param_sets_with_rna_id_pairs.insert(rna_id_pair, StaFeParams::origin());
      bpap_mats_with_rna_id_pairs.insert(rna_id_pair, Prob4dMat::default());
    }
  }
  for (rna_id_pair, sta_fe_params) in sta_fe_param_sets_with_rna_id_pairs.iter_mut() {
    let seq_len_pair = (fasta_records[rna_id_pair.0].seq.len(), fasta_records[rna_id_pair.1].seq.len());
    let max_pos_dist_4_el = offset_4_max_pos_dist_on_el + get_seq_len_diff(&seq_len_pair);
    let max_bp_span_pair = (max_bp_spans[rna_id_pair.0], max_bp_spans[rna_id_pair.1]);
    let ref ref_2_fasta_records = fasta_records;
    let ref ref_2_invert_exp_max_free_energies = invert_exp_max_free_energies;
    let ref ref_2_bpp_mats = sparse_bpp_mats;
    *sta_fe_params = StaFeParams::new(rna_id_pair, ref_2_fasta_records, ref_2_invert_exp_max_free_energies, &max_bp_span_pair, max_pos_dist_4_il, max_pos_dist_4_el, max_substr_dist, ref_2_bpp_mats, exp_opening_gap_penalty, exp_extending_gap_penalty);
  }
  thread_pool.scoped(|scope| {
    for (rna_id_pair, bpap_mat) in bpap_mats_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
      let max_bp_span_pair = (max_bp_spans[rna_id_pair.0], max_bp_spans[rna_id_pair.1]);
      let max_pos_dist_4_el = offset_4_max_pos_dist_on_el + get_seq_len_diff(&seq_len_pair);
      let ref sta_fe_params = sta_fe_param_sets_with_rna_id_pairs[&rna_id_pair];
      scope.execute(move || {
        *bpap_mat = io_algo_4_bpap_mat(&seq_pair, &seq_len_pair, sta_fe_params, &max_bp_span_pair, max_pos_dist_4_il, max_pos_dist_4_el, max_substr_dist);
      });
    }
  });
  thread_pool.scoped(|scope| {
    for (rna_id, bpp_mat, upp_mat) in multizip((0 .. num_of_fasta_records, bpp_mats.iter_mut(), upp_mats.iter_mut())) {
      let ref ref_2_bpap_mats_with_rna_id_pairs = bpap_mats_with_rna_id_pairs;
      scope.execute(move || {
        let prob_mat_pair = pct_of_bpp_and_upp_mat(ref_2_bpap_mats_with_rna_id_pairs, rna_id, num_of_fasta_records, bpp_mat, upp_mat.len());
        *bpp_mat = prob_mat_pair.0;
        *upp_mat = prob_mat_pair.1;
      });
    }
  });
  let output_file_header = format!(" in this file = \"{}\".\n; The values of the parameters used to the matrices are as follows.\n; \"opening_gap_penalty\" = {}, \"extending_gap_penalty\" = {}, \"min_bpp\" = {}, \"max_pos_dist_4_il\" = {}, \"offset_4_max_pos_dist_on_el\" = {}, \"max_substr_dist\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), opening_gap_penalty, extending_gap_penalty, min_bpp, max_pos_dist_4_il, offset_4_max_pos_dist_on_el, max_substr_dist, num_of_threads);
  let bpp_mat_on_sta_file_path = output_dir_path.join(BPP_MAT_ON_STA_FILE_NAME);
  let upp_mat_on_sta_file_path = output_dir_path.join(UPP_MAT_ON_STA_FILE_NAME);
  let mut writer_2_bpp_mat_on_sta_file = BufWriter::new(File::create(bpp_mat_on_sta_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_bpp_mat_on_sta_file = format!("; The version {} of the Zprob program.\n; The path to the input file in order to compute the average base-pairing probability matrices on structural alignment in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of the input RNA sequence correpsonding to the row. The next row to the row is with the average base-pairing probability matrix on structural alignment of the sequence.";
  for (rna_id, bpp_mat) in bpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in bpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
    }
    buf_4_writer_2_bpp_mat_on_sta_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_on_sta_file.write_all(buf_4_writer_2_bpp_mat_on_sta_file.as_bytes());
  let mut writer_2_upp_mat_on_sta_file = BufWriter::new(File::create(upp_mat_on_sta_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_upp_mat_on_sta_file = format!("; The version {} of the Zprob program.\n; The path to the input file in order to compute the average unpairing probability matrices on structural alignment in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The next row to the row is with the average unpairing probability matrix on structural alignment of the sequence.";
  for (rna_id, upp_mat) in upp_mats.iter().enumerate() {
    let seq_len = fasta_records[rna_id].seq.len();
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (i, &upp) in upp_mat.iter().enumerate() {
      if i == 0 || i == seq_len - 1 {continue;}
      buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
    }
    buf_4_writer_2_upp_mat_on_sta_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_upp_mat_on_sta_file.write_all(buf_4_writer_2_upp_mat_on_sta_file.as_bytes());
}
