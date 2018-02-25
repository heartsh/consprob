extern crate io_algo_4_rna_stap_mp;
extern crate bio_seq_algos;
extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate time;
// extern crate num;

use io_algo_4_rna_stap_mp::*;
use bio_seq_algos::utils::*;
use bio_seq_algos::durbin_algo::*;
use rna_algos::mccaskill_algo::*;
use std::f32::consts::LOG2_E;
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
  let sta_sps = StaScoringParams::new((0.00_000_5 as StaScore).ln(), (0.00_000_5 as StaScore).ln(), -(0.00_000_5 as StaScore).ln() / LOG2_E.sqrt(), -11., -1., -11., -1.);
  let bppt_4_sta_ss = 0.25;
  let bpp_mp = (mccaskill_algo(&TEST_SEQ_PAIR.0[..]), mccaskill_algo(&TEST_SEQ_PAIR.1[..]));
  let slp = (TEST_SEQ_PAIR.0.len(), TEST_SEQ_PAIR.1.len());
  let begin = precise_time_s();
  let stap_mp = io_algo_4_rna_stap_mp(&slp, &bpp_mp, &log_bap_matrix, &sta_sps, bppt_4_sta_ss);
  let elapsed_time = precise_time_s() - begin;
  println!("The elapsed time = {}[s].", elapsed_time);
  // println!("The BPAP matrix pair for the seq. pair \"{}\" and \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ_PAIR.0[..]), String::from_utf8_lossy(&TEST_SEQ_PAIR.1[..]), &stap_mp.bpap_matrix);
  for sub_bpap_matrix in &stap_mp.bpap_matrix {
    for sub_sub_bpap_matrix in sub_bpap_matrix {
      for bpaps in sub_sub_bpap_matrix {
        for &bpap in bpaps {assert!((0. <= bpap && bpap <= 1.));}
      }
    }
  }
  // println!("The BAP matrix pair for the seq. pair \"{}\" and \"{}\" = \"{:?}\".", String::from_utf8_lossy(&TEST_SEQ_PAIR.0[..]), String::from_utf8_lossy(&TEST_SEQ_PAIR.1[..]), &stap_mp.bap_matrix);
  /* for baps in &stap_mp.bap_matrix {
    for &bap in baps {assert!((0. <= bap && bap <= 1.));}
  } */
}
