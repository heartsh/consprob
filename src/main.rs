extern crate stem;
extern crate scoped_threadpool;
extern crate bio;
extern crate itertools;
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
use itertools::multizip;

type NumOfThreads = u32;

const DEFAULT_MIN_BPP: Prob = 0.001;
const DEFAULT_MAX_GAP_NUM: usize = 3;
const DEFAULT_MAX_BP_SPAN: usize = 0;
const DEFAULT_NUM_OF_TIMES_OF_EXPECT_MAX_ITER: usize = 1;
const BPP_MAT_ON_SS_FILE_NAME: &'static str = "bpp_mats_on_ss.dat";
const BAP_MAT_FILE_NAME: &'static str = "bap_mats.dat";
const OGP_MAT_FILE_NAME_1: &'static str = "ogp_mats_1.dat";
const OGP_MAT_FILE_NAME_2: &'static str = "ogp_mats_2.dat";
const EGP_MAT_FILE_NAME_1: &'static str = "egp_mats_1.dat";
const EGP_MAT_FILE_NAME_2: &'static str = "egp_mats_2.dat";
const BPAP_MAT_FILE_NAME_1: &'static str = "bpap_mats_1.dat";
const BPAP_MAT_FILE_NAME_2: &'static str = "bpap_mats_2.dat";
const BPAP_MAT_FILE_NAME_3: &'static str = "bpap_mats_3.dat";
const OGPP_MAT_FILE_NAME_1: &'static str = "ogpp_mats_1.dat";
const OGPP_MAT_FILE_NAME_2: &'static str = "ogpp_mats_2.dat";
const EGPP_MAT_FILE_NAME_1: &'static str = "egpp_mats_1.dat";
const EGPP_MAT_FILE_NAME_2: &'static str = "egpp_mats_2.dat";
const LGP_MAT_FILE_NAME_1: &'static str = "lgp_mats_1.dat";
const LGP_MAT_FILE_NAME_2: &'static str = "lgp_mats_2.dat";
const RGP_MAT_FILE_NAME_1: &'static str = "rgp_mats_1.dat";
const RGP_MAT_FILE_NAME_2: &'static str = "rgp_mats_2.dat";
const BPP_MAT_ON_STA_FILE_NAME: &'static str = "bpp_mats_on_sta.dat";
const NBPP_MAT_FILE_NAME: &'static str = "nbpp_mats.dat";
const VERSION: &'static str = "0.1.0";

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "max_gap_num", &format!("The maximum number of the gaps in a structural alignment (Uses {} by default)", DEFAULT_MAX_GAP_NUM), "UINT");
  opts.optopt("", "max_base_pairing_span", &format!("A maximum base-pairing span (Uses {} by default)", DEFAULT_MAX_BP_SPAN), "FLOAT");
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
  let min_bpp = if opts.opt_present("min_base_pair_prob") {
    opts.opt_str("min_base_pair_prob").expect("Failed to get a minimum base-pairing-probability from command arguments.").parse().expect("Failed to parse a minimum base-pairing-probability.")
  } else {
    DEFAULT_MIN_BPP
  };
  let min_lbpp = min_bpp.ln();
  let max_gap_num = if opts.opt_present("max_gap_num") {
    opts.opt_str("gap_num").expect("Failed to get the maximum number of the gaps in a structural alignment from command arguments.").parse().expect("Failed to parse the maximum number of the gaps in a structural alignment.")
  } else {
    DEFAULT_MAX_GAP_NUM
  };
  let max_bp_span = if opts.opt_present("max_base_pairing_span") {
    opts.opt_str("max_base_pairing_span").expect("Failed to get a maximum base-pairing span from command arguments.").parse().expect("Failed to parse a max base-pairing span.")
  } else {
    DEFAULT_MAX_BP_SPAN
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
  let mut original_fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.expect("Failed to read a FASTA record.");
    let seq = unsafe {from_utf8_unchecked(fasta_record.seq()).to_uppercase().as_bytes().iter().filter(|&&base| {is_rna_base(base)}).map(|&base| {base}).collect()};
    original_fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let num_of_fasta_records = original_fasta_records.len();
  let mut thread_pool = Pool::new(num_of_threads);
  let mut lbpp_mats = vec![SparseLogProbMat::default(); num_of_fasta_records];
  let mut lnbpp_mats = vec![LogProbs::new(); num_of_fasta_records];
  let mut fasta_records = original_fasta_records.clone();
  let max_seq_len = original_fasta_records.iter().map(|fasta_record| fasta_record.seq.len()).max().expect("Failed to get the maximum among sequence lengths.");
  thread_pool.scoped(|scope| {
    for (lbpp_mat, lnbpp_mat, fasta_record) in multizip((lbpp_mats.iter_mut(), lnbpp_mats.iter_mut(), fasta_records.iter_mut())) {
      let seq_len = fasta_record.seq.len();
      let num_of_pseudo_bases = max_seq_len - seq_len;
      scope.execute(move || {
        *lbpp_mat = get_lbpp_mat(&fasta_record.seq[..], max_bp_span);
        *lnbpp_mat = get_lnbpp_mat(lbpp_mat, seq_len);
        *lbpp_mat = remove_little_lbpps_from_lbpp_mat(lbpp_mat, min_lbpp);
        append_pseudo_bases(&mut fasta_record.seq, num_of_pseudo_bases);
        lnbpp_mat.append(&mut vec![0.; num_of_pseudo_bases]);
        lnbpp_mat.insert(0, NEG_INFINITY);
        lnbpp_mat.push(NEG_INFINITY);
        lbpp_mat.insert((0, max_seq_len + 1), 0.);
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let mut buf_4_writer_2_bpp_mat_on_ss_file = format!("; The version {} of the STEM program.\n; The path to the input file to compute the Base-Pairing Probability Matrices (= BPPMs) on secondary structure (= SS) in this file = \"{}\".\n; The values of the parameters used to compute these matrices are as follows.\n; \"min_base_pair_prob\" = {}, \"max_base_pairing_span\" = {}, \"num_of_threads\" = {}.", VERSION, input_file_path.display(), min_bpp, max_bp_span, num_of_threads) + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the BPPM of this sequence on SS.";
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
      let ref ref_2_fasta_records = fasta_records;
      let ref ref_2_lbpp_mats = lbpp_mats;
      let ref ref_2_lnbpp_mats = lnbpp_mats;
      scope.execute(move || {
        *sta_fe_params = StaFeParams::new(rna_id_pair, ref_2_fasta_records, max_gap_num, ref_2_lbpp_mats, ref_2_lnbpp_mats);
      });
    }
  });
  for _ in 0 .. num_of_times_of_expect_max_iter {
    thread_pool.scoped(|scope| {
      for (rna_id_pair, lstapmt) in lstapmts_with_rna_id_pairs.iter_mut() {
        let seq_len_pair = (fasta_records[rna_id_pair.0].seq.len(), fasta_records[rna_id_pair.1].seq.len());
        let ref sta_fe_params = sta_fe_param_sets_with_rna_id_pairs[&rna_id_pair];
        scope.execute(move || {
          *lstapmt = io_algo_4_rna_lstapmt(&seq_len_pair, sta_fe_params, max_gap_num);
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
    thread_pool.scoped(|scope| {
      for (rna_id, lbpp_mat, lnbpp_mat) in multizip((0 .. num_of_fasta_records, lbpp_mats.iter_mut(), lnbpp_mats.iter_mut())) {
        let ref ref_2_lstapmts_with_rna_id_pairs = lstapmts_with_rna_id_pairs;
        let seq_len = fasta_records[rna_id].seq.len();
        scope.execute(move || {
          *lbpp_mat = pct_of_lbpp_mat(ref_2_lstapmts_with_rna_id_pairs, rna_id, num_of_fasta_records, lbpp_mat);
          *lnbpp_mat = pct_of_lnbpp_mat(ref_2_lstapmts_with_rna_id_pairs, rna_id, num_of_fasta_records, lnbpp_mat, seq_len);
        });
      }
    });
    for (rna_id_pair, lstapmt) in lstapmts_with_rna_id_pairs.iter_mut() {
      *lstapmt = lstapmts_with_rna_id_pairs_from_pct[rna_id_pair].clone();
    }
    update_lstapmts_with_rna_id_pairs(&mut lstapmts_with_rna_id_pairs, &lbpp_mats, &lnbpp_mats);
  }
  let stapmts_with_rna_id_pairs = get_stapmts_with_rna_id_pairs(&lstapmts_with_rna_id_pairs);
  let output_file_header = format!(" in this file = \"{}\".\n; The values of the parameters used to compute these matrices are as follows.\n; \"min_base_pair_prob\" = {}, \"max_gap_num\" = {}, \"max_base_pairing_span\" = {}, \"num_of_times_of_expect_max_iter\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), min_bpp, max_gap_num, max_bp_span, num_of_times_of_expect_max_iter, num_of_threads);
  let bap_mat_file_path = output_dir_path.join(BAP_MAT_FILE_NAME);
  let ogp_mat_file_path_1 = output_dir_path.join(OGP_MAT_FILE_NAME_1);
  let ogp_mat_file_path_2 = output_dir_path.join(OGP_MAT_FILE_NAME_2);
  let egp_mat_file_path_1 = output_dir_path.join(EGP_MAT_FILE_NAME_1);
  let egp_mat_file_path_2 = output_dir_path.join(EGP_MAT_FILE_NAME_2);
  let bpap_mat_file_path_1 = output_dir_path.join(BPAP_MAT_FILE_NAME_1);
  let bpap_mat_file_path_2 = output_dir_path.join(BPAP_MAT_FILE_NAME_2);
  let bpap_mat_file_path_3 = output_dir_path.join(BPAP_MAT_FILE_NAME_3);
  let ogpp_mat_file_path_1 = output_dir_path.join(OGPP_MAT_FILE_NAME_1);
  let ogpp_mat_file_path_2 = output_dir_path.join(OGPP_MAT_FILE_NAME_2);
  let egpp_mat_file_path_1 = output_dir_path.join(EGPP_MAT_FILE_NAME_1);
  let egpp_mat_file_path_2 = output_dir_path.join(EGPP_MAT_FILE_NAME_2);
  let lgp_mat_file_path_1 = output_dir_path.join(LGP_MAT_FILE_NAME_1);
  let lgp_mat_file_path_2 = output_dir_path.join(LGP_MAT_FILE_NAME_2);
  let rgp_mat_file_path_1 = output_dir_path.join(RGP_MAT_FILE_NAME_1);
  let rgp_mat_file_path_2 = output_dir_path.join(RGP_MAT_FILE_NAME_2);
  let bpp_mat_on_sta_file_path = output_dir_path.join(BPP_MAT_ON_STA_FILE_NAME);
  let nbpp_mat_file_path = output_dir_path.join(NBPP_MAT_FILE_NAME);
  let mut writer_2_bap_mat_file = BufWriter::new(File::create(bap_mat_file_path).expect("Failed to create an output file."));
  let mut writer_2_ogp_mat_file_1 = BufWriter::new(File::create(ogp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_ogp_mat_file_2 = BufWriter::new(File::create(ogp_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_egp_mat_file_1 = BufWriter::new(File::create(egp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_egp_mat_file_2 = BufWriter::new(File::create(egp_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_bpap_mat_file_1 = BufWriter::new(File::create(bpap_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_bpap_mat_file_2 = BufWriter::new(File::create(bpap_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_bpap_mat_file_3 = BufWriter::new(File::create(bpap_mat_file_path_3).expect("Failed to create an output file."));
  let mut writer_2_ogpp_mat_file_1 = BufWriter::new(File::create(ogpp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_ogpp_mat_file_2 = BufWriter::new(File::create(ogpp_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_egpp_mat_file_1 = BufWriter::new(File::create(egpp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_egpp_mat_file_2 = BufWriter::new(File::create(egpp_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_lgp_mat_file_1 = BufWriter::new(File::create(lgp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_lgp_mat_file_2 = BufWriter::new(File::create(lgp_mat_file_path_2).expect("Failed to create an output file."));
  let mut writer_2_rgp_mat_file_1 = BufWriter::new(File::create(rgp_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_rgp_mat_file_2 = BufWriter::new(File::create(rgp_mat_file_path_2).expect("Failed to create an output file."));
  let mut buf_4_writer_2_bap_mat_file = format!("; The version {} of the STEM program.\n; The path to the input file for computing base alignment probability matrices (= BAPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the BAPM between these 2 sequences.";
  let mut buf_4_writer_2_ogp_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing opening gap probability matrices (= OGPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the OGPM between these 2 sequences.";
  let mut buf_4_writer_2_ogp_mat_file_2 = buf_4_writer_2_ogp_mat_file_1.clone();
  let mut buf_4_writer_2_egp_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing extending gap probability matrices (= EGPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the EGPM between these 2 sequences.";
  let mut buf_4_writer_2_egp_mat_file_2 = buf_4_writer_2_egp_mat_file_1.clone();
  let mut buf_4_writer_2_bpap_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing base pair alignment probability matrices (= BPAPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the BPAPM between these 2 sequences.";
  let mut buf_4_writer_2_bpap_mat_file_2 = buf_4_writer_2_bpap_mat_file_1.clone();
  let mut buf_4_writer_2_bpap_mat_file_3 = buf_4_writer_2_bpap_mat_file_1.clone();
  let mut buf_4_writer_2_ogpp_mat_file_1 = format!("; The version {} of the the STEM program.\n; The path to the input file for computing opening gap pair probability matrices (= OGPPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the OGPPM between these 2 sequences.";
  let mut buf_4_writer_2_ogpp_mat_file_2 = buf_4_writer_2_ogpp_mat_file_1.clone();
  let mut buf_4_writer_2_egpp_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing extending gap pair probability matrices (= EGPPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the EGPPM between these 2 sequences.";
  let mut buf_4_writer_2_egpp_mat_file_2 = buf_4_writer_2_egpp_mat_file_1.clone();
  let mut buf_4_writer_2_lgp_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing left gap probability matrices (= LGPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the LGPM between these 2 sequences.";
  let mut buf_4_writer_2_lgp_mat_file_2 = buf_4_writer_2_lgp_mat_file_1.clone();
  let mut buf_4_writer_2_rgp_mat_file_1 = format!("; The version {} of the STEM program.\n; The path to the input file for computing right gap probability matrices (= RGPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of each of 2 RNA sequences. The row next to this row is with the RGPM between these 2 sequences.";
  let mut buf_4_writer_2_rgp_mat_file_2 = buf_4_writer_2_rgp_mat_file_1.clone();
  for (rna_id_pair, stapmt) in &stapmts_with_rna_id_pairs {
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    let seq_len_pair = (original_fasta_records[rna_id_pair.0].seq.len(), original_fasta_records[rna_id_pair.1].seq.len());
    for (&(i, j), bap) in stapmt.base_align_prob_mat.iter() {
      if i > seq_len_pair.0 || j > seq_len_pair.1 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, bap));
    }
    buf_4_writer_2_bap_mat_file.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 0 .. seq_len_pair.0 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i, stapmt.opening_gap_prob_mat_1[i + 1]));
    }
    buf_4_writer_2_ogp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 0 .. seq_len_pair.1 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i, stapmt.ogp_mat_2[i + 1]));
    }
    buf_4_writer_2_ogp_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 0 .. seq_len_pair.0 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i, stapmt.extending_gap_prob_mat_1[i + 1]));
    }
    buf_4_writer_2_egp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for i in 0 .. seq_len_pair.1 {
      buf_4_rna_id_pair.push_str(&format!("{},{} ", i, stapmt.egp_mat_2[i + 1]));
    }
    buf_4_writer_2_egp_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k, l), &bpap) in stapmt.base_pair_align_prob_mat_1.iter() {
      if i == 0 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{},{} ", i - 1, j - 1, k - 1, l - 1, bpap));
    }
    buf_4_writer_2_bpap_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k, l), &bpap) in stapmt.bpap_mat_2.iter() {
      if i > seq_len_pair.0 || j > seq_len_pair.0 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{},{} ", i - 1, j - 1, k - 1, l - 1, bpap));
    }
    buf_4_writer_2_bpap_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k, l), &bpap) in stapmt.bpap_mat_3.iter() {
      if k > seq_len_pair.1 || l > seq_len_pair.1 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{},{} ", i - 1, j - 1, k - 1, l - 1, bpap));
    }
    buf_4_writer_2_bpap_mat_file_3.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j), &ogpp) in stapmt.opening_gap_pair_prob_mat_1.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, ogpp));
    }
    buf_4_writer_2_ogpp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j), &ogpp) in stapmt.ogpp_mat_2.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, ogpp));
    }
    buf_4_writer_2_ogpp_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j), &egpp) in stapmt.extending_gap_pair_prob_mat_1.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, egpp));
    }
    buf_4_writer_2_egpp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j), &egpp) in stapmt.egpp_mat_2.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, egpp));
    }
    buf_4_writer_2_egpp_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k), &lgp) in stapmt.left_gap_prob_mat_1.iter() {
      if k > seq_len_pair.1 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{} ", i - 1, j - 1, k - 1, lgp));
    }
    buf_4_writer_2_lgp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k), &lgp) in stapmt.lgp_mat_2.iter() {
      if i > seq_len_pair.0 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{} ", i - 1, j - 1, k - 1, lgp));
    }
    buf_4_writer_2_lgp_mat_file_2.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k), &rgp) in stapmt.right_gap_prob_mat_1.iter() {
      if k > seq_len_pair.1 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{} ", i - 1, j - 1, k - 1, rgp));
    }
    buf_4_writer_2_rgp_mat_file_1.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k), &rgp) in stapmt.rgp_mat_2.iter() {
      if i > seq_len_pair.0 {continue;}
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{} ", i - 1, j - 1, k - 1, rgp));
    }
    buf_4_writer_2_rgp_mat_file_2.push_str(&buf_4_rna_id_pair);
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
    let _ = writer_2_bpap_mat_file_1.write_all(buf_4_writer_2_bpap_mat_file_1.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_bpap_mat_file_2.write_all(buf_4_writer_2_bpap_mat_file_2.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_bpap_mat_file_3.write_all(buf_4_writer_2_bpap_mat_file_3.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_ogpp_mat_file_1.write_all(buf_4_writer_2_ogpp_mat_file_1.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_ogpp_mat_file_2.write_all(buf_4_writer_2_ogpp_mat_file_2.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_egpp_mat_file_1.write_all(buf_4_writer_2_egpp_mat_file_1.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_egpp_mat_file_2.write_all(buf_4_writer_2_egpp_mat_file_2.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_lgp_mat_file_1.write_all(buf_4_writer_2_lgp_mat_file_1.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_lgp_mat_file_2.write_all(buf_4_writer_2_lgp_mat_file_2.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_rgp_mat_file_1.write_all(buf_4_writer_2_rgp_mat_file_1.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_rgp_mat_file_2.write_all(buf_4_writer_2_rgp_mat_file_2.as_bytes());
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
  thread::spawn(move || {
    let _ = writer_2_bpp_mat_on_sta_file.write_all(buf_4_writer_2_bpp_mat_on_sta_file.as_bytes());
  });
  let mut buf_4_writer_2_nbpp_mat_file = format!("; The version {} of the STEM program.\n; The path to the input file for computing the Not Base-Pairing Probability Matrices (= NBPPMs) in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The row next to this row is with the NBPPM of this sequence.";
  let mut writer_2_nbpp_mat_file = BufWriter::new(File::create(nbpp_mat_file_path).expect("Failed to create an output file."));
  for (rna_id, lnbpp_mat) in lnbpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    let seq_len = fasta_records[rna_id].seq.len();
    for i in 0 .. seq_len {
      buf_4_rna_id.push_str(&format!("{},{} ", i, lnbpp_mat[i + 1].exp()));
    }
    buf_4_writer_2_nbpp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_nbpp_mat_file.write_all(buf_4_writer_2_nbpp_mat_file.as_bytes());
}
