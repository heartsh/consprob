extern crate io_algo_4_rna_stapmq;
extern crate bio_seq_algos;
extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;

use io_algo_4_rna_stapmq::*;
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
fn test_stapmq() {
  let seq_len_pair = (TEST_SEQ_PAIR.0.len(), TEST_SEQ_PAIR.1.len());
  let mut ca_score_mat = CaScoreMat::default();
  let alphabet = b"AUGC";
  for &base_1 in alphabet {
    for &base_2 in alphabet {
      ca_score_mat.insert((base_1, base_2), if base_1 == base_2 {0.1} else {-0.1});
    }
  }
  let sa_scoring_params = SaScoringParams::new(&ca_score_mat, -1., -0.1);
  let mut lbap_mat = get_log_cap_mat(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]), &sa_scoring_params);
  for lbaps in &mut lbap_mat {
    for lbap in lbaps {
      *lbap = *lbap - LN_2;
    }
  }
  let bpp_mat_pair = (mccaskill_algo(&TEST_SEQ_PAIR.0[..]), mccaskill_algo(&TEST_SEQ_PAIR.1[..]));
  let log_prob_1 = (0.01 as StaScore).log2();
  let log_prob_2 = (0.00_001 as StaScore).log2();
  let sta_scoring_params = StaScoringParams::new(log_prob_1, log_prob_1, -log_prob_1 * 2., 0.5, log_prob_2, log_prob_1, log_prob_2, log_prob_1);
  let min_bpp_1 = 0. as Prob;
  let min_lbap_1 = min_bpp_1.log2();
  let begin = precise_time_s();
  let stapmq = io_algo_4_rna_stapmq(&seq_len_pair, &bpp_mat_pair, &lbap_mat, &sta_scoring_params, min_bpp_1, min_lbap_1);
  let elapsed_time = precise_time_s() - begin;
  check_stapmq(&stapmq);
  let min_bpp_2 = 0.01 as Prob;
  let min_lbap_2 = (0.05 as Prob).log2();
  let begin = precise_time_s();
  let stapmq = io_algo_4_rna_stapmq(&seq_len_pair, &bpp_mat_pair, &lbap_mat, &sta_scoring_params, min_bpp_2, min_lbap_2);
  let acceleration = elapsed_time / (precise_time_s() - begin);
  check_stapmq(&stapmq);
  println!("The acceleration with the minimum Base-Pairing-Probability (= BPP) {} and minimum Base-Alignment-Probability (= BAP) {} is {}-fold compared with the minimum BPP {} and minimum BAP {}.", min_bpp_1, min_lbap_1.exp2(), acceleration, min_bpp_2, min_lbap_2.exp2());
}

fn check_stapmq(stapmq: &Stapmq) {
  for &bpap in stapmq.base_pair_align_prob_mat.values() {
    assert!((0. <= bpap && bpap <= 1.));
  }
  for &bap in stapmq.base_align_prob_mat.values() {
    assert!(0. <= bap && bap <= 1.);
  }
  for &bpip in stapmq.base_pair_indel_prob_mat_1.values() {
    assert!(0. <= bpip && bpip <= 1.);
  }
  for &bpip in stapmq.bpip_mat_2.values() {
    assert!(0. <= bpip && bpip <= 1.);
  }
}
