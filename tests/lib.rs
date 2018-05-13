extern crate strap;
extern crate bio_seq_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;

use strap::*;
use bio_seq_algos::durbin_algo::*;
use time::precise_time_s;

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
  let lbap_mat = get_log_cap_mat(&(&TEST_SEQ_PAIR.0[..], &TEST_SEQ_PAIR.1[..]), &sa_scoring_params);
  let log_prob_1 = (0.01 as StaScore).log2();
  let log_prob_2 = (0.00_001 as StaScore).log2();
  let sta_scoring_params = StaScoringParams::new(log_prob_1, log_prob_1, -log_prob_1 * 2., 0.5, log_prob_2, log_prob_1, log_prob_2, log_prob_1);
  let min_bpp_1 = 0.;
  let min_lbap_1 = (min_bpp_1 as Prob).log2();
  let max_bp_span = 200;
  let sparse_lbpp_mat_pair = (&get_lbpp_mat_from_parasor(&TEST_SEQ_PAIR.0[..], max_bp_span, min_bpp_1), &get_lbpp_mat_from_parasor(&TEST_SEQ_PAIR.1[..], max_bp_span, min_bpp_1));
  let sparse_lbap_mat = get_sparse_lbap_mat(&lbap_mat, min_lbap_1);
  let begin = precise_time_s();
  let stapmq = io_algo_4_rna_stapmq(&seq_len_pair, &sparse_lbpp_mat_pair, &sparse_lbap_mat, &sta_scoring_params);
  let elapsed_time = precise_time_s() - begin;
  check_stapmq(&stapmq);
  let min_bpp_2 = 0.01;
  let min_lbap_2 = (0.01 as Prob).log2();
  let sparse_lbpp_mat_pair = (&get_lbpp_mat_from_parasor(&TEST_SEQ_PAIR.0[..], max_bp_span, min_bpp_2), &get_lbpp_mat_from_parasor(&TEST_SEQ_PAIR.1[..], max_bp_span, min_bpp_2));
  let sparse_lbap_mat = get_sparse_lbap_mat(&lbap_mat, min_lbap_2);
  let begin = precise_time_s();
  let stapmq = io_algo_4_rna_stapmq(&seq_len_pair, &sparse_lbpp_mat_pair, &sparse_lbap_mat, &sta_scoring_params);
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
