extern crate io_algo_4_rna_stap_matrix_pair;
extern crate bio_seq_algos;
extern crate rna_algos;
#[macro_use]
extern crate lazy_static;

use io_algo_4_rna_stap_matrix_pair::*;
use bio_seq_algos::utils::*;
use bio_seq_algos::durbin_algo::*;
use rna_algos::mccaskill_algo::*;
use std::f64::consts::LOG2_E;

type SeqPair = (Seq, Seq);

lazy_static! {
  static ref TEST_SEQ_PAIR: SeqPair = {
    (
      String::from("AUGCAAGGGGGCUUUAACAC").into_bytes(),
      String::from("GAUGCAUGCAAGGGCGCUUUGACA").into_bytes(),
    )
  };
  static ref TS_LEN_PAIR: (usize, usize) = {
    (TEST_SEQ_PAIR.0.len(), TEST_SEQ_PAIR.1.len())
  };
}

#[test]
fn test_stap_mp() {
  let mut ca_sm = CaScoreMatrix::default();
  let alphabet = b"AUGC";
  for &base_1 in alphabet {
    for &base_2 in alphabet {
      ca_sm.insert((base_1, base_2), if base_1 == base_2 {1.} else {-1.});
    }
  }
  let sa_sps = SaScoringParams::new(&ca_sm, -4., -1.);
  let log_bap_matrix = get_log_cap_matrix(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]), &sa_sps);
  let log_bpp_mp = (get_log_bpp_matrix(&TEST_SEQ_PAIR.0[..]), get_log_bpp_matrix(&TEST_SEQ_PAIR.1[..]));
  let sta_sps = StaScoringParams::new((0.00_000_5 as f64).ln(), (0.00_000_5 as f64).ln(), -(0.00_000_5 as f64).ln() / LOG2_E.sqrt(), -11., -1., -11., -1.);
  let bpp_mp = (get_bpp_matrix(&log_bpp_mp.0), get_bpp_matrix(&log_bpp_mp.1));
  let nbpp_mp = (
    bpp_mp.0.iter().map(|xs| xs.iter().map(|&x| 1. - x).collect()).collect::<LogProbMatrix>(),
    bpp_mp.1.iter().map(|xs| xs.iter().map(|&x| 1. - x).collect()).collect::<LogProbMatrix>(),
  );
  let log_nbpp_mp = (
    nbpp_mp.0.iter().map(|xs| xs.iter().map(|&x| fast_ln(x)).collect()).collect::<ProbMatrix>(),
    nbpp_mp.1.iter().map(|xs| xs.iter().map(|&x| fast_ln(x)).collect()).collect::<ProbMatrix>(),
  );
  let stap_mp = io_algo_4_rna_stap_mp(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]), &log_bpp_mp, &log_bap_matrix, &bpp_mp, &nbpp_mp, &log_nbpp_mp, &sta_sps);
  println!("The BPAP matrix pair for the seq. pair \"{}\" and \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ_PAIR.0[..]), String::from_utf8_lossy(&TEST_SEQ_PAIR.1[..]), &stap_mp.bpap_matrix);
  for sub_bpap_matrix in &stap_mp.bpap_matrix {
    for sub_sub_bpap_matrix in sub_bpap_matrix {
      for bpaps in sub_sub_bpap_matrix {
        for &bpap in bpaps {assert!((0. <= bpap && bpap <= 1.));}
      }
    }
  }
  println!("The BAP matrix pair for the seq. pair \"{}\" and \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ_PAIR.0[..]), String::from_utf8_lossy(&TEST_SEQ_PAIR.1[..]), &stap_mp.bap_matrix);
  for baps in &stap_mp.bap_matrix {
    for &bap in baps {assert!((0. <= bap && bap <= 1.));}
  }
}
