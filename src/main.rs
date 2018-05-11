extern crate io_algo_4_rna_stapmq;
extern crate bio_seq_algos;
#[macro_use]
extern crate lazy_static;
extern crate getopts;
extern crate scoped_threadpool;
extern crate bio;
extern crate itertools;
extern crate num_cpus;

use io_algo_4_rna_stapmq::*;
use bio_seq_algos::durbin_algo::*;
use getopts::Options;
use self::scoped_threadpool::Pool;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::fs::create_dir;
use std::process::{Command, Output};
use std::str::from_utf8_unchecked;
use std::thread;

type Arg = String;
type NumOfThreads = u32;
type FastaId = String;
type FastaRecord = (FastaId, Seq);
type FastaRecords = Vec<FastaRecord>;
type LstapmqsWithRnaIdPairs = HashMap<RnaIdPair, LogStapmq, Hasher>;

const DEFAULT_BMS: SaScore = 0.1;
const DEFAULT_BMMS: SaScore = -0.1;
const DEFAULT_BOGP_ON_SA: SaScore = -1.;
const DEFAULT_BEGP_ON_SA: SaScore = -0.1;
lazy_static! {
  static ref DEFAULT_LOG_NH_BPP: StaScore = {(0.01 as StaScore).log2()};
  static ref DEFAULT_LOG_NH_BAP: StaScore = *DEFAULT_LOG_NH_BPP;
  static ref DEFAULT_SCALE_PARAM_4_BPA_SCORE: StaScore = - *DEFAULT_LOG_NH_BPP * 2.;
  static ref DEFAULT_BOGP_ON_STA: StaScore = {(0.00_001 as StaScore).log2()};
  static ref DEFAULT_BEGP_ON_STA: StaScore = *DEFAULT_LOG_NH_BPP;
  static ref DEFAULT_LOGP: StaScore = *DEFAULT_BOGP_ON_STA;
  static ref DEFAULT_LEGP: StaScore = *DEFAULT_BEGP_ON_STA;
}
const DEFAULT_OFFSET_BPA_SCORE: StaScore = 0.5;
const DEFAULT_MIN_BPP: Prob = 0.00_1;
const DEFAULT_MIN_BAP: Prob = 0.03;
const DEFAULT_MAX_BP_SPAN: usize = 200;
const DEFAULT_NUM_OF_TIMES_OF_IMPROVEMENTS_OF_STAPMQS: usize = 5;
const PARASOR_COMMAND: &'static str = "ParasoR";
const BPP_MAT_ON_SS_FILE_NAME: &'static str = "bpp_mats_on_ss.dat";
const BAP_MAT_ON_SA_FILE_NAME: &'static str = "bap_mats_on_sa.dat";
const BPAP_MAT_FILE_NAME: &'static str = "bpap_mats.dat";
const BAP_MAT_ON_STA_FILE_NAME: &'static str = "bap_mats_on_sta.dat";
const BPIP_MAT_FILE_NAME_1: &'static str = "bpip_mats_1.dat";
const BPIP_MAT_FILE_NAME_2: &'static str = "bpip_mats_2.dat";
const BPP_MAT_ON_STA_FILE_NAME: &'static str = "bpp_mats_on_sta.dat";

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "base_match_score", &format!("A base match score (Uses {} by default)", DEFAULT_BMS), "FLOAT");
  opts.optopt("", "base_mismatch_score", &format!("A base mismatch score (Uses {} by default)", DEFAULT_BMMS), "FLOAT");
  opts.optopt("", "log_null_hypothesis_base_pairing_prob", &format!("The logarithm Base-Pairing Probability (= BPP) on a null hypothesis (Uses {} by default)", *DEFAULT_LOG_NH_BPP), "FLOAT");
  opts.optopt("", "log_null_hypothesis_base_align_prob", &format!("The logarithm base alignment probability (= BAP) on a null hypothesis (Uses {} by default)", *DEFAULT_LOG_NH_BAP), "FLOAT");
  opts.optopt("", "scale_param_4_base_pair_align_score", &format!("A scale parameter for a base pair alignment score (Uses {} by default)", *DEFAULT_SCALE_PARAM_4_BPA_SCORE), "FLOAT");
  opts.optopt("", "offset_base_pair_align_score", &format!("An offset base pair alignment score (Uses {} by default)", DEFAULT_OFFSET_BPA_SCORE), "FLOAT");
  opts.optopt("", "base_opening_gap_penalty_on_sa", &format!("A base opening gap penalty on sequence alignment (= SA) (Uses {} by default)", DEFAULT_BOGP_ON_SA), "FLOAT");
  opts.optopt("", "base_extending_gap_penalty_on_sa", &format!("A base extending gap penalty on SA (Uses {} by default)", DEFAULT_BEGP_ON_SA), "FLOAT");
  opts.optopt("", "base_opening_gap_penalty_on_sta", &format!("A base opening gap penalty on STructural Alignment (= STA) (Uses {} by default)", *DEFAULT_BOGP_ON_STA), "FLOAT");
  opts.optopt("", "base_extending_gap_penalty_on_sta", &format!("A base extending gap penalty on STA (Uses {} by default)", *DEFAULT_BEGP_ON_STA), "FLOAT");
  opts.optopt("", "loop_opening_gap_penalty", &format!("A loop opening gap penalty (Uses {} by default)", *DEFAULT_LOGP), "FLOAT");
  opts.optopt("", "loop_extending_gap_penalty", &format!("A loop extending gap penalty (Uses {} by default)", *DEFAULT_LEGP), "FLOAT");
  opts.optopt("", "min_base_pairing_prob", &format!("A minimum BPP (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "min_base_align_prob", &format!("A minimum BAP (Uses {} by default)", DEFAULT_MIN_BAP), "FLOAT");
  opts.optopt("", "max_base_pairing_span", &format!("A maximum base-pairing span (Uses {} by default)", DEFAULT_MAX_BP_SPAN), "FLOAT");
  opts.optopt("", "num_of_times_of_improvements_of_stapmqs", &format!("The number of times of the improvements of structural-alignment probability quadruples (Uses {} by default)", DEFAULT_NUM_OF_TIMES_OF_IMPROVEMENTS_OF_STAPMQS), "FLOAT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_file_path = opts.opt_str("i").expect("Failed to get the path to an input FASTA file containing RNA sequences from command arguments.");
  let input_file_path = Path::new(&input_file_path);
  let output_dir_path = opts.opt_str("o").expect("Failed to get the path to an output directory from command arguments.");
  let output_dir_path = Path::new(&output_dir_path);
  let bms = if opts.opt_present("base_match_score") {
    opts.opt_str("base_match_score").expect("Failed to get a base match score from command arguments.").parse().expect("Failed to parse a base match score.")
  } else {
    DEFAULT_BMS
  };
  let bmms = if opts.opt_present("base_mismatch_score") {
    opts.opt_str("base_mismatch_score").expect("Failed to get a base mismatch score from command arguments.").parse().expect("Failed to parse a base mismatch score.")
  } else {
    DEFAULT_BMMS
  };
  let bogp_on_sa = if opts.opt_present("base_opening_gap_penalty_on_sa") {
    opts.opt_str("base_opening_gap_penalty_on_sa").expect("Failed to get a base opening gap penalty on sequence alignment from command arguments.").parse().expect("Failed to parse a base opening gap penalty on sequence alignment.")
  } else {
    DEFAULT_BOGP_ON_SA
  };
  let begp_on_sa = if opts.opt_present("base_extending_gap_penalty_on_sa") {
    opts.opt_str("base_extending_gap_penalty_on_sa").expect("Failed to get a base extending gap penalty on sequence alignment from command arguments.").parse().expect("Failed to parse a base extending gap penalty on sequence alignment.")
  } else {
    DEFAULT_BEGP_ON_SA
  };
  let log_nh_bpp = if opts.opt_present("log_null_hypothesis_base_pairing_prob") {
    opts.opt_str("log_null_hypothesis_base_pairing_prob").expect("Failed to get the logarithm base-pairing probability on a null hypothesis from command arguments.").parse().expect("Failed to parse the logarithm base-pairing probability on a null hypothesis.")
  } else {
    *DEFAULT_LOG_NH_BPP
  };
  let log_nh_bap = if opts.opt_present("log_null_hypothesis_base_align_prob") {
    opts.opt_str("log_null_hypothesis_base_align_prob").expect("Failed to get the logarithm base alignment probability on a null hypothesis from command arguments.").parse().expect("Failed to parse the logarithm base alignment probability on a null hypothesis.")
  } else {
    *DEFAULT_LOG_NH_BAP
  };
  let scale_param_4_bpa_score = if opts.opt_present("scale_param_4_base_pair_align_score") {
    opts.opt_str("scale_param_4_base_pair_align_score").expect("Failed to get a scale parameter for a base pair alignment score from command arguments.").parse().expect("Failed to parse a scale parameter for a base pair alignment score.")
  } else {
    *DEFAULT_SCALE_PARAM_4_BPA_SCORE
  };
  let offset_bpa_score = if opts.opt_present("offset_base_pair_align_score") {
    opts.opt_str("offset_base_pair_align_score").expect("Failed to get an offset base pair alignment score from command arguments.").parse().expect("Failed to parse an offset base pair alignment score.")
  } else {
    DEFAULT_OFFSET_BPA_SCORE
  };
  let bogp_on_sta = if opts.opt_present("base_opening_gap_penalty_on_sta") {
    opts.opt_str("base_opening_gap_penalty_on_sta").expect("Failed to get a base opening gap penalty on structural alignment from command arguments.").parse().expect("Failed to parse a base opening gap penalty on structural alignment.")
  } else {
    *DEFAULT_BOGP_ON_STA
  };
  let begp_on_sta = if opts.opt_present("base_opening_gap_penalty_on_sta") {
    opts.opt_str("base_extending_gap_penalty_on_sta").expect("Failed to get a base extending gap penalty on structural alignment from command arguments.").parse().expect("Failed to parse a base extending gap penalty on structural alignment.")
  } else {
    *DEFAULT_BEGP_ON_STA
  };
  let logp = if opts.opt_present("loop_opening_gap_penalty") {
    opts.opt_str("loop_opening_gap_penalty").expect("Failed to get a loop opening gap penalty from command arguments.").parse().expect("Failed to parse a loop opening gap penalty.")
  } else {
    *DEFAULT_LOGP
  };
  let legp = if opts.opt_present("loop_extending_gap_penalty") {
    opts.opt_str("loop_extending_gap_penalty").expect("Failed to get a loop extending gap penalty from command arguments.").parse().expect("Failed to parse a loop extending gap penalty.")
  } else {
    *DEFAULT_LEGP
  };
  let min_bpp = if opts.opt_present("min_base_pairing_prob") {
    opts.opt_str("min_base_pairing_prob").expect("Failed to get a minimum base pairing probability from command arguments.").parse().expect("Failed to parse a minimum base pairing probability.")
  } else {
    DEFAULT_MIN_BPP
  };
  let min_lbpp = min_bpp.log2();
  let min_bap = if opts.opt_present("min_base_align_prob") {
    opts.opt_str("min_base_align_prob").expect("Failed to get a minimum base alignment probability from command arguments.").parse().expect("Failed to parse a minimum base alignment probability.")
  } else {
    DEFAULT_MIN_BAP
  };
  let min_lbap = min_bap.log2();
  let max_bp_span = if opts.opt_present("max_base_pairing_span") {
    opts.opt_str("max_base_pairing_span").expect("Failed to get a maximum base-pairing span from command arguments.").parse().expect("Failed to parse a max base-pairing span.")
  } else {
    DEFAULT_MAX_BP_SPAN
  };
  let num_of_times_of_improvements_of_stapmqs = if opts.opt_present("max_base_pairing_span") {
    opts.opt_str("num_of_times_of_improvements").expect("Failed to get the number of times of the improvements of structural-alignment probability quadruples from command arguments.").parse().expect("Failed to parse the number of times of the improvements of structural-alignment probability quadruples.")
  } else {
    DEFAULT_NUM_OF_TIMES_OF_IMPROVEMENTS_OF_STAPMQS
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
    let seq = unsafe {from_utf8_unchecked(fasta_record.seq()).to_uppercase().as_bytes().iter().filter(|&&base| {is_rna_base(base)}).map(|&base| {base}).collect()};
    fasta_records.push((String::from(fasta_record.id().expect("Failed to get the ID of a FASTA record.")), seq));
  }
  let num_of_fasta_records = fasta_records.len();
  let mut lbpp_mats = vec![SparseLogProbMat::default(); num_of_fasta_records];
  let mut thread_pool = Pool::new(num_of_threads);
  thread_pool.scoped(|scope| {
    for (lbpp_mat, fasta_record) in lbpp_mats.iter_mut().zip(fasta_records.iter()) {
      scope.execute(move || {
        *lbpp_mat = get_lbpp_mat_from_parasor(&fasta_record.1[..], max_bp_span, min_bpp);
      });
    }
  });
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let mut buf_4_writer_2_bpp_mat_on_ss_file = format!("; The path to the input file to compute the base-pairing matrices on secondary structure in this file = \"{}\".\n; The values of the parameters used to compute these matrices are as follows.\n; \"min_base_pairing_prob\" = {}, \"max_base_pairing_span\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), min_bpp, max_bp_span, num_of_threads);
  let bpp_mat_on_ss_file_path = output_dir_path.join(BPP_MAT_ON_SS_FILE_NAME);
  let mut writer_2_bpp_mat_on_ss_file = BufWriter::new(File::create(bpp_mat_on_ss_file_path).expect("Failed to create an output file."));
  for (rna_id, lbpp_mat) in lbpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &lbpp) in lbpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, lbpp.exp2()));
    }
    buf_4_writer_2_bpp_mat_on_ss_file.push_str(&buf_4_rna_id);
  }
  thread::spawn(move || {
    let _ = writer_2_bpp_mat_on_ss_file.write_all(buf_4_writer_2_bpp_mat_on_ss_file.as_bytes());
  });
  let mut lbap_mats_with_rna_id_pairs = LogProbMatsWithRnaIdPairs::default();
  let mut lstapmqs_with_rna_id_pairs = LstapmqsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      lbap_mats_with_rna_id_pairs.insert((rna_id_1, rna_id_2), SparseLogProbMat::default());
      lstapmqs_with_rna_id_pairs.insert((rna_id_1, rna_id_2), LogStapmq::new());
    }
  }
  let mut ca_score_mat = CaScoreMat::default();
  let alphabet = b"AUGC";
  for &base_1 in alphabet {
    for &base_2 in alphabet {
      ca_score_mat.insert((base_1, base_2), if base_1 == base_2 {bms} else {bmms});
    }
  }
  let sa_scoring_params = SaScoringParams::new(&ca_score_mat, bogp_on_sa, begp_on_sa);
  thread_pool.scoped(|scope| {
    for (rna_id_pair, lbap_mat) in lbap_mats_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].1[..], &fasta_records[rna_id_pair.1].1[..]);
      let ref ref_2_sa_scoring_params = sa_scoring_params;
      scope.execute(move || {
        *lbap_mat = get_sparse_lbap_mat(&get_log_cap_mat(&seq_pair, ref_2_sa_scoring_params), min_lbap);
      });
    }
  });
  let mut lbap_mats_with_rna_id_pairs_from_pct = LogProbMatsWithRnaIdPairs::default();
  let mut lstapmqs_with_rna_id_pairs_from_pct = LstapmqsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      lbap_mats_with_rna_id_pairs_from_pct.insert((rna_id_1, rna_id_2), SparseLogProbMat::default());
      lstapmqs_with_rna_id_pairs_from_pct.insert((rna_id_1, rna_id_2), LogStapmq::new());
    }
  }
  thread_pool.scoped(|scope| {
    for (rna_id_pair, lbap_mat) in lbap_mats_with_rna_id_pairs_from_pct.iter_mut() {
      let ref ref_2_lbap_mats_with_rna_id_pairs = lbap_mats_with_rna_id_pairs;
      scope.execute(move || {
        *lbap_mat = remove_little_lbaps_from_sparse_lbap_mat(&pct_of_lbap_mat(ref_2_lbap_mats_with_rna_id_pairs, rna_id_pair, num_of_fasta_records), min_lbap);
      });
    }
  });
  let mut lbap_mats_with_rna_id_pairs = lbap_mats_with_rna_id_pairs_from_pct;
  let mut buf_4_writer_2_bap_mat_on_sa_file = format!("; The path to the input file to compute the base alignment probability matrices on sequence alignment in this file = \"{}\".\n; The values of the parameters used to compute these matrices are as follows.\n; \"base_match_score\" = {}, \"base_mismatch_score\" = {}, \"base_opening_gap_penalty_on_sa\" = {}, \"base_extending_gap_penalty_on_sa\" = {}, \"min_base_align_prob\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), bms, bmms, bogp_on_sa, begp_on_sa, min_bap, num_of_threads);
  let bap_mat_on_sa_file_path = output_dir_path.join(BAP_MAT_ON_SA_FILE_NAME);
  let mut writer_2_bap_mat_on_sa_file = BufWriter::new(File::create(bap_mat_on_sa_file_path).expect("Failed to create an output file."));
  for (rna_id_pair, lbap_mat) in &lbap_mats_with_rna_id_pairs {
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j), &lbap) in lbap_mat.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, lbap.exp2()));
    }
    buf_4_writer_2_bap_mat_on_sa_file.push_str(&buf_4_rna_id_pair);
  }
  thread::spawn(move || {
    let _ = writer_2_bap_mat_on_sa_file.write_all(buf_4_writer_2_bap_mat_on_sa_file.as_bytes());
  });
  let sta_scoring_params = StaScoringParams::new(log_nh_bpp, log_nh_bap, scale_param_4_bpa_score, offset_bpa_score, bogp_on_sta, begp_on_sta, logp, legp);
  for i in 0 .. num_of_times_of_improvements_of_stapmqs + 1 {
    thread_pool.scoped(|scope| {
      for (rna_id_pair, lstapmq) in lstapmqs_with_rna_id_pairs.iter_mut() {
        let seq_len_pair = (fasta_records[rna_id_pair.0].1.len(), fasta_records[rna_id_pair.1].1.len());
        let lbpp_mat_pair = (&lbpp_mats[rna_id_pair.0], &lbpp_mats[rna_id_pair.1]);
        let ref lbap_mat = lbap_mats_with_rna_id_pairs[rna_id_pair];
        let ref ref_2_sta_scoring_params = sta_scoring_params;
        scope.execute(move || {
          *lstapmq = io_algo_4_rna_lstapmq(&seq_len_pair, &lbpp_mat_pair, lbap_mat, ref_2_sta_scoring_params, min_lbpp);
        });
      }
    });
    thread_pool.scoped(|scope| {
      for (rna_id_pair, lstapmq) in lstapmqs_with_rna_id_pairs_from_pct.iter_mut() {
        let ref ref_2_lstapmqs_with_rna_id_pairs = lstapmqs_with_rna_id_pairs;
        scope.execute(move || {
          *lstapmq = prob_cons_transformation_of_lstapmq(ref_2_lstapmqs_with_rna_id_pairs, rna_id_pair, num_of_fasta_records);
        });
      }
    });
    let mut lbpp_mat_pairs_with_rna_id_pairs = LogProbMatPairsWithRnaIdPairs::default();
    for rna_id_1 in 0 .. num_of_fasta_records {
      for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
        lbpp_mat_pairs_with_rna_id_pairs.insert((rna_id_1, rna_id_2), (SparseLogProbMat::default(), SparseLogProbMat::default()));
      }
    }
    thread_pool.scoped(|scope| {
      for (rna_id_pair, lbpp_mat_pair) in lbpp_mat_pairs_with_rna_id_pairs.iter_mut() {
        let ref ref_2_lstapmqs_with_rna_id_pairs = lstapmqs_with_rna_id_pairs_from_pct;
        scope.execute(move || {
          *lbpp_mat_pair = get_lbpp_mat_pair(ref_2_lstapmqs_with_rna_id_pairs, rna_id_pair);
        });
      }
    });
    thread_pool.scoped(|scope| {
      for (rna_id, lbpp_mat) in lbpp_mats.iter_mut().enumerate() {
        let ref ref_2_lbpp_mat_pairs_with_rna_id_pairs = lbpp_mat_pairs_with_rna_id_pairs;
        scope.execute(move || {
          *lbpp_mat = pct_of_lbpp_mat(ref_2_lbpp_mat_pairs_with_rna_id_pairs, rna_id, num_of_fasta_records);
          if i < num_of_times_of_improvements_of_stapmqs {
            *lbpp_mat = remove_little_lbpps_from_sparse_lbpp_mat(lbpp_mat, min_lbpp);
          }
        });
      }
      if i < num_of_times_of_improvements_of_stapmqs {
        for (rna_id_pair, lbap_mat) in lbap_mats_with_rna_id_pairs.iter_mut() {
          let ref ref_2_lbap_mat = lstapmqs_with_rna_id_pairs_from_pct[rna_id_pair].log_bap_mat;
          scope.execute(move || {
            *lbap_mat = remove_little_lbaps_from_sparse_lbap_mat(ref_2_lbap_mat, min_lbap);
          });
        }
      }
    });
  }
  let stapmqs_with_rna_id_pairs = get_stapmqs_with_rna_id_pairs(&lstapmqs_with_rna_id_pairs_from_pct);
  let output_file_header = format!(" in this file = \"{}\".\n; The values of the parameters used to compute these matrices are as follows.\n; \"base_match_score\" = {}, \"base_mismatch_score\" = {}, \"base_opening_gap_penalty_on_sa\" = {}, \"base_extending_gap_penalty_on_sa\" = {}, \"log_null_hypothesis_base_pairing_prob\" = {}, \"log_null_hypothesis_base_align_prob\" = {}, \"scale_param_4_base_pair_align_score\" = {}, \"offset_base_pair_align_score\" = {}, \"base_opening_gap_penalty_on_sta\" = {}, \"base_extending_gap_penalty_on_sta\" = {}, \"loop_opening_gap_penalty\" = {}, \"loop_extending_gap_penalty\" = {}, \"min_base_pairing_prob\" = {}, \"min_base_align_prob\" = {}, \"max_base_pairing_span\" = {}, \"num_of_times_of_improvements_of_stapmqs\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), bms, bmms, bogp_on_sa, begp_on_sa, log_nh_bpp, log_nh_bap, scale_param_4_bpa_score, offset_bpa_score, bogp_on_sta, begp_on_sta, logp, legp, min_bpp, min_bap, max_bp_span, num_of_times_of_improvements_of_stapmqs, num_of_threads);
  let bpap_mat_file_path = output_dir_path.join(BPAP_MAT_FILE_NAME);
  let bap_mat_on_sta_file_path = output_dir_path.join(BAP_MAT_ON_STA_FILE_NAME);
  let bpip_mat_file_path_1 = output_dir_path.join(BPIP_MAT_FILE_NAME_1);
  let bpip_mat_file_path_2 = output_dir_path.join(BPIP_MAT_FILE_NAME_2);
  let bpp_mat_on_sta_file_path = output_dir_path.join(BPP_MAT_ON_STA_FILE_NAME);
  let mut buf_4_writer_2_bpp_mat_on_sta_file = String::from("; The path to the input file to compute the base-pairing matrices on structural alignment") + &output_file_header;
  let mut writer_2_bpp_mat_on_sta_file = BufWriter::new(File::create(bpp_mat_on_sta_file_path).expect("Failed to create an output file."));
  for (rna_id, lbpp_mat) in lbpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &lbpp) in lbpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, lbpp.exp2()));
    }
    buf_4_writer_2_bpp_mat_on_sta_file.push_str(&buf_4_rna_id);
  }
  thread::spawn(move || {
    let _ = writer_2_bpp_mat_on_sta_file.write_all(buf_4_writer_2_bpp_mat_on_sta_file.as_bytes());
  });
  let mut writer_2_bpap_mat_file = BufWriter::new(File::create(bpap_mat_file_path).expect("Failed to create an output file."));
  let mut writer_2_bap_mat_on_sta_file = BufWriter::new(File::create(bap_mat_on_sta_file_path).expect("Failed to create an output file."));
  let mut writer_2_bpip_mat_file_1 = BufWriter::new(File::create(bpip_mat_file_path_1).expect("Failed to create an output file."));
  let mut writer_2_bpip_mat_file_2 = BufWriter::new(File::create(bpip_mat_file_path_2).expect("Failed to create an output file."));
  let mut buf_4_writer_2_bpap_mat_file = String::from("; The path to the input file to compute the base pair alignment probability matrices") + &output_file_header;
  let mut buf_4_writer_2_bap_mat_on_sta_file = String::from("; The path to the input file to compute the base alignment probability matrices on structural alignment") + &output_file_header;
  let mut buf_4_writer_2_bpip_mat_file_1 = String::from("; The path to the input file to compute the base pair indel matrices") + &output_file_header;
  let mut buf_4_writer_2_bpip_mat_file_2 = buf_4_writer_2_bpip_mat_file_1.clone();
  for (rna_id_pair, stapmq) in &stapmqs_with_rna_id_pairs {
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j, k, l), &bpap) in stapmq.base_pair_align_prob_mat.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{},{},{} ", i - 1, j - 1, k - 1, l - 1, bpap));
    }
    buf_4_writer_2_bpap_mat_file.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
    for (&(i, j), &bap) in stapmq.base_align_prob_mat.iter() {
      buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - 1, j - 1, bap));
    }
    buf_4_writer_2_bap_mat_on_sta_file.push_str(&buf_4_rna_id_pair);
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id_pair.0);
    for (&(i, j), &bpip) in stapmq.base_pair_indel_prob_mat_1.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpip));
    }
    buf_4_writer_2_bpip_mat_file_1.push_str(&buf_4_rna_id);
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id_pair.1);
    for (&(i, j), &bpip) in stapmq.bpip_mat_2.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpip));
    }
    buf_4_writer_2_bpip_mat_file_2.push_str(&buf_4_rna_id);
  }
  thread::spawn(move || {
    let _ = writer_2_bpap_mat_file.write_all(buf_4_writer_2_bpap_mat_file.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_bap_mat_on_sta_file.write_all(buf_4_writer_2_bap_mat_on_sta_file.as_bytes());
  });
  thread::spawn(move || {
    let _ = writer_2_bpip_mat_file_1.write_all(buf_4_writer_2_bpip_mat_file_1.as_bytes());
  });
  let _ = writer_2_bpip_mat_file_2.write_all(buf_4_writer_2_bpip_mat_file_2.as_bytes());
}

#[inline]
fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
}

#[inline]
fn get_lbpp_mat_from_parasor(seq: SeqSlice, max_bp_span: usize, min_bpp: Prob) -> SparseProbMat {
  let arg_1 = format!("--bpp={}", min_bpp);
  let arg_2 = max_bp_span.to_string();
  let args = unsafe {vec!["-f", from_utf8_unchecked(seq), "--pre", &arg_1, "--constraint", &arg_2]};
  let parasor_output = unsafe {String::from_utf8_unchecked(run_command(PARASOR_COMMAND, &args, "Failed to run ParasoR.").stdout)};
  parasor_output.lines().filter(|line| {!line.starts_with("#")}).map(|line| {
    let strings = line.trim().splitn(4, '\t').collect::<Vec<&str>>();
    (
      (
        strings[1].parse().expect("Failed to parse a string."),
        strings[2].parse().expect("Failed to parse a string."),
      ),
      strings[3].parse::<Prob>().expect("Failed to parse a string.").log2()
    )
  }).collect::<SparseLogProbMat>()
}

fn run_command(command: &str, args: &[&str], expect: &str) -> Output {
  Command::new(command).args(args).output().expect(expect)
}
