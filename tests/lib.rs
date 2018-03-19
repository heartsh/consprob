extern crate io_algo_4_rna_stapmp;
extern crate bio_seq_algos;
extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;

use io_algo_4_rna_stapmp::*;
use bio_seq_algos::durbin_algo::*;
use rna_algos::mccaskill_algo::*;
use time::precise_time_s;
use std::f64::consts::LN_2;

type SeqPair = (Seq, Seq);

lazy_static! {
  static ref TEST_SEQ_PAIR: SeqPair = {
    (
      String::from("AUGCAAGGGGGCUUUAACAC").into_bytes(),
      String::from("GAUGCAUGCAAGGGCGCUUUGACA").into_bytes(),
    )
  };
}

#[test]
fn test_stapmp() {
  let slp = (TEST_SEQ_PAIR.0.len(), TEST_SEQ_PAIR.1.len());
  let mut ca_sm = CaScoreMatrix::default();
  let alphabet = b"AUGC";
  for &base_1 in alphabet {
    for &base_2 in alphabet {
      ca_sm.insert((base_1, base_2), if base_1 == base_2 {0.1} else {-0.1});
    }
  }
  let sa_sps = SaScoringParams::new(&ca_sm, -1., -0.1);
  let log_bap_matrix = get_log_cap_matrix(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]), &sa_sps);
  let log_prob_1 = (0.01 as StaScore).ln();
  let log_prob_2 = (0.00_001 as StaScore).ln();
  let sta_sps = StaScoringParams::new(log_prob_1, log_prob_1, -log_prob_1 / LN_2.sqrt(), log_prob_2, log_prob_1, log_prob_2, log_prob_1);
  let bppt_4_sta_ss_1 = 0.;
  let bpp_mp = (mccaskill_algo(&TEST_SEQ_PAIR.0[..]), mccaskill_algo(&TEST_SEQ_PAIR.1[..]));
  let begin = precise_time_s();
  let stapmp = io_algo_4_rna_stapmp(&slp, &bpp_mp, &log_bap_matrix, &sta_sps, bppt_4_sta_ss_1);
  check_stapmp(&stapmp);
  let elapsed_time = precise_time_s() - begin;
  let bppt_4_sta_ss_2 = 0.00_1;
  let begin = precise_time_s();
  let stapmp = io_algo_4_rna_stapmp(&slp, &bpp_mp, &log_bap_matrix, &sta_sps, bppt_4_sta_ss_2);
  check_stapmp(&stapmp);
  let acceleration = elapsed_time / (precise_time_s() - begin);
  println!("The acceleration with the Base-Pairing Probability Threshold (= BPPT) = {} is {}-fold compared with the BPPT = {}.", bppt_4_sta_ss_2, acceleration, bppt_4_sta_ss_1);
}

fn check_stapmp(stapmp: &Stapmp) {
  for sub_bpap_matrix in &stapmp.bpap_matrix {
    for sub_sub_bpap_matrix in sub_bpap_matrix {
      for bpaps in sub_sub_bpap_matrix {
        for &bpap in bpaps {assert!((0. <= bpap && bpap <= 1.));}
      }
    }
  }
  for baps in &stapmp.bap_matrix {
    for &bap in baps {assert!((0. <= bap && bap <= 1.));}
  }
}
