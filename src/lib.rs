extern crate bio_seq_algos;
extern crate rna_algos;
extern crate itertools;

pub use rna_algos::utils::*;
use itertools::multizip;
use std::f64::consts::LN_2;
use std::process::{Command, Output};
pub use std::str::from_utf8_unchecked;

pub type PosQuadruple = (Pos, Pos, Pos, Pos);
pub type Prob4dMat = HashMap<PosQuadruple, Prob, Hasher>;
pub type SparseProbMat = HashMap<PosPair, Prob, Hasher>;
pub struct Stapmq {
  pub base_pair_align_prob_mat: Prob4dMat, 
  pub base_align_prob_mat: SparseProbMat,
  pub base_pair_indel_prob_mat_1: SparseProbMat,
  pub bpip_mat_2: SparseProbMat,
}
pub type LogProb4dMat = HashMap<PosQuadruple, LogProb, Hasher>;
pub type SparseLogProbMat = HashMap<PosPair, LogProb, Hasher>;
#[derive(Clone)]
pub struct LogStapmq {
  pub log_bpap_mat: LogProb4dMat,
  pub log_bap_mat: SparseLogProbMat,
  pub log_bpip_mat_1: SparseLogProbMat,
  pub log_bpip_mat_2: SparseLogProbMat,
}
type LogPpf4dMat = HashMap<PosQuadruple, LogPf, Hasher>;
type LogPpf4dMatPair = (LogPpf4dMat, LogPpf4dMat);
pub type StaScore = LogProb;
pub struct StaScoringParams {
  pub log_null_hypothesis_bpp: LogProb,
  pub log_nh_bap: LogProb,
  pub scale_param_4_bpa_score: StaScore,
  pub offset_bpa_score: StaScore,
  pub base_opening_gap_penalty: StaScore,
  pub base_extending_gap_penalty: StaScore,
  pub loop_opening_gap_penalty: StaScore,
  pub loop_extending_gap_penalty: StaScore,
}
type ProbDist = Probs;
type LogProbDist = LogProbs;
type ProbDistSlice<'a> = &'a[Prob];
type LogProbDistSlice<'a> = &'a[LogProb];
type ProbDistDiv = LogProb;
type ProbDistDist = ProbDistDiv;
type PosSpan = usize;
pub type Poss = Vec<Pos>;
type PosSeqsWithPosSpans = HashMap<PosSpan, Poss, Hasher>;
type PosSpanPair = (PosSpan, PosSpan);
pub type PosPairs = Vec<PosPair>;
type PosPairSeqsWithPosSpanPairs = HashMap<PosSpanPair, PosPairs, Hasher>;
type SparseProbMatPair = (SparseProbMat, SparseProbMat);
type LogProbMatPair = (SparseLogProbMat, SparseLogProbMat);
pub type PosSeqsWithPoss = HashMap<Pos, Poss, Hasher>;
type PosSeqSetPairWithPoss = (PosSeqsWithPoss, PosSeqsWithPoss);
pub type PosPairSeqsWithPosPairs = HashMap<PosPair, PosPairs, Hasher>;
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
type StapmqsWithRnaIdPairs = HashMap<RnaIdPair, Stapmq, Hasher>;
type LstapmqsWithRnaIdPairs = HashMap<RnaIdPair, LogStapmq, Hasher>;
pub type LogProbMatsWithRnaIdPairs = HashMap<RnaIdPair, SparseLogProbMat, Hasher>;
pub type Ref2LogProbMatPair<'a> = (&'a SparseLogProbMat, &'a SparseLogProbMat);
type SeqsOfEpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, ExpPartOfTerm4LogProb, Hasher>;
pub type LogProbMatPairsWithRnaIdPairs = HashMap<RnaIdPair, LogProbMatPair, Hasher>;
type SeqsOfEpsOfTerms4LogProbsWithPosQuadruples = HashMap<PosQuadruple, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosQuadruples = HashMap<PosQuadruple, ExpPartOfTerm4LogProb, Hasher>;

const PARASOR_COMMAND: &'static str = "ParasoR";

impl LogStapmq {
  pub fn new() -> LogStapmq {
    LogStapmq {
      log_bpap_mat: LogProb4dMat::default(),
      log_bap_mat: SparseLogProbMat::default(), 
      log_bpip_mat_1: SparseLogProbMat::default(), 
      log_bpip_mat_2: SparseLogProbMat::default(), 
    }
  }
}

impl StaScoringParams {
  pub fn new(lnhbpp: LogProb, lnhbap: LogProb, sp4bpas: StaScore, obpas: StaScore, bogp: StaScore, begp: StaScore, logp: StaScore, legp: StaScore) -> StaScoringParams {
    StaScoringParams {
      log_null_hypothesis_bpp: lnhbpp,
      log_nh_bap: lnhbap,
      scale_param_4_bpa_score: sp4bpas,
      offset_bpa_score: obpas,
      base_opening_gap_penalty: bogp,
      base_extending_gap_penalty: begp,
      loop_opening_gap_penalty: logp,
      loop_extending_gap_penalty: legp,
    }
  }
}

#[inline]
pub fn io_algo_4_rna_stapmq(seq_len_pair: &(usize, usize), lbpp_mat_pair: &Ref2LogProbMatPair, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams) -> Stapmq {
  let lstapmq = io_algo_4_rna_lstapmq(seq_len_pair, lbpp_mat_pair, lbap_mat, sta_scoring_params);
  get_stapmq(&lstapmq)
}

#[inline]
pub fn io_algo_4_rna_lstapmq(seq_len_pair: &(usize, usize), lbpp_mat_pair: &Ref2LogProbMatPair, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams) -> LogStapmq {
  let mut pos_seq_set_pair_with_pos_spans = (PosSeqsWithPosSpans::default(), PosSeqsWithPosSpans::default());
  let mut bpp_mat_pair = (SparseProbMat::default(), SparseProbMat::default());
  let mut not_bpp_mat_pair = (SparseProbMat::default(), SparseProbMat::default());
  let mut log_nbpp_mat_pair = (SparseLogProbMat::default(), SparseLogProbMat::default());
  for substr_len in 2 .. seq_len_pair.0 + 1 {
    let mut poss = Poss::new();
    for i in 1 .. seq_len_pair.0 - substr_len + 2 {
      let j = i + substr_len - 1;
      match lbpp_mat_pair.0.get(&(i, j)) {
        Some(&lbpp) => {
          let bpp = lbpp.exp2();
          let nbpp = 1. - bpp;
          let lnbpp = nbpp.log2();
          poss.push(i);
          let pos_pair = (i, j);
          bpp_mat_pair.0.insert(pos_pair, bpp);
          not_bpp_mat_pair.0.insert(pos_pair, nbpp);
          log_nbpp_mat_pair.0.insert(pos_pair, lnbpp);
        }, None => {},
      }
    }
    if poss.len() > 0 {
      pos_seq_set_pair_with_pos_spans.0.insert(substr_len, poss);
    }
  }
  for substr_len in 2 .. seq_len_pair.1 + 1 {
    let mut poss = Poss::new();
    for i in 1 .. seq_len_pair.1 - substr_len + 2 {
      let j = i + substr_len - 1;
      match lbpp_mat_pair.1.get(&(i, j)) {
        Some(&lbpp) => {
          let bpp = lbpp.exp2();
          let nbpp = 1. - bpp;
          let lnbpp = nbpp.log2();
          poss.push(i);
          let pos_pair = (i, j);
          bpp_mat_pair.1.insert(pos_pair, bpp);
          not_bpp_mat_pair.1.insert(pos_pair, nbpp);
          log_nbpp_mat_pair.1.insert(pos_pair, lnbpp);
        }, None => {},
      }
    }
    if poss.len() > 0 {
      pos_seq_set_pair_with_pos_spans.1.insert(substr_len, poss);
    }
  }
  let mut pos_pair_seqs_with_pos_span_pairs = PosPairSeqsWithPosSpanPairs::default();
  for substr_len_1 in 2 .. seq_len_pair.0 + 1 {
    match pos_seq_set_pair_with_pos_spans.0.get(&substr_len_1) {
      Some(poss_1) => {
        for substr_len_2 in 2 .. seq_len_pair.1 + 1 {
          match pos_seq_set_pair_with_pos_spans.1.get(&substr_len_2) {
            Some(poss_2) => {
              let mut pos_pairs = PosPairs::new();
              for &i in poss_1 {
                for &k in poss_2 {
                  match lbap_mat.get(&(i, k)) {
                    Some(_) => {
                      let j = i + substr_len_1 - 1;
                      let l = k + substr_len_2 - 1;
                      match lbap_mat.get(&(j, l)) {
                        Some(_) => {
                          pos_pairs.push((i, k));
                        }, None => {},
                      }
                    }, None => {},
                  }
                }
              }
              if pos_pairs.len() > 0 {
                pos_pair_seqs_with_pos_span_pairs.insert((substr_len_1, substr_len_2), pos_pairs);
              }
            }, None => {},
          }
        }
      }, None => {},
    }
  }
  pos_pair_seqs_with_pos_span_pairs.insert((seq_len_pair.0 + 2, seq_len_pair.1 + 2), vec![(0, 0)]);
  let mut pos_seq_set_pair_with_poss_4_forward_bps = (PosSeqsWithPoss::default(), PosSeqsWithPoss::default());
  let mut pos_seq_set_pair_with_poss_4_backward_bps = (PosSeqsWithPoss::default(), PosSeqsWithPoss::default());
  for &(i, j) in lbpp_mat_pair.0.keys() {
    let poss_with_j_exists = match pos_seq_set_pair_with_poss_4_forward_bps.0.get(&j) {
      Some(_) => {true},
      None => {false},
    };
    if poss_with_j_exists {
      pos_seq_set_pair_with_poss_4_forward_bps.0.get_mut(&j).expect("Failed to get an element from a hash map.").push(i);
    } else {
      pos_seq_set_pair_with_poss_4_forward_bps.0.insert(j, vec![i]);
    }
    let poss_with_i_exists = match pos_seq_set_pair_with_poss_4_backward_bps.0.get(&i) {
      Some(_) => {true},
      None => {false},
    };
    if poss_with_i_exists {
      pos_seq_set_pair_with_poss_4_backward_bps.0.get_mut(&i).expect("Failed to get an element from a hash map.").push(j);
    } else {
      pos_seq_set_pair_with_poss_4_backward_bps.0.insert(i, vec![j]);
    }
  }
  for &(i, j) in lbpp_mat_pair.1.keys() {
    let poss_with_j_exists = match pos_seq_set_pair_with_poss_4_forward_bps.1.get(&j) {
      Some(_) => {true},
      None => {false},
    };
    if poss_with_j_exists {
      pos_seq_set_pair_with_poss_4_forward_bps.1.get_mut(&j).expect("Failed to get an element from a hash map.").push(i);
    } else {
      pos_seq_set_pair_with_poss_4_forward_bps.1.insert(j, vec![i]);
    }
    let poss_with_i_exists = match pos_seq_set_pair_with_poss_4_backward_bps.1.get(&i) {
      Some(_) => {true},
      None => {false},
    };
    if poss_with_i_exists {
      pos_seq_set_pair_with_poss_4_backward_bps.1.get_mut(&i).expect("Failed to get an element from a hash map.").push(j);
    } else {
      pos_seq_set_pair_with_poss_4_backward_bps.1.insert(i, vec![j]);
    }
  }
  let mut pos_seqs_with_poss_4_bas = PosSeqsWithPoss::default();
  for &(i, j) in lbap_mat.keys() {
    let poss_with_i_exists = match pos_seqs_with_poss_4_bas.get(&i) {
      Some(_) => {true},
      None => {false},
    };
    if poss_with_i_exists {
      pos_seqs_with_poss_4_bas.get_mut(&i).expect("Failed to get an element from a hash map.").push(j);
    } else {
      pos_seqs_with_poss_4_bas.insert(i, vec![j]);
    }
  }
  let mut pos_pair_seqs_with_pos_pairs_4_forward_bpas = PosPairSeqsWithPosPairs::default();
  let mut pos_pair_seqs_with_pos_pairs_4_backward_bpas = PosPairSeqsWithPosPairs::default();
  for i in 1 .. seq_len_pair.0 + 1 {
    match pos_seqs_with_poss_4_bas.get(&i) {
      Some(poss) => {
        for &k in poss {
          let mut forward_pos_pairs = PosPairs::new();
          match pos_seq_set_pair_with_poss_4_forward_bps.0.get(&i) {
            Some(poss_1) => {
              match pos_seq_set_pair_with_poss_4_forward_bps.1.get(&k) {
                Some(poss_2) => {
                  for &j in poss_1 {
                    for &l in poss_2 {
                      match lbap_mat.get(&(j, l)) {
                        Some(_) => {
                          forward_pos_pairs.push((j, l));
                        }, None => {},
                      }
                    }
                  }
                }, None => {},
              }
            }, None => {},
          }
          if forward_pos_pairs.len() > 0 {
            pos_pair_seqs_with_pos_pairs_4_forward_bpas.insert((i, k), forward_pos_pairs);
          }
          let mut backward_pos_pairs = PosPairs::new();
          match pos_seq_set_pair_with_poss_4_backward_bps.0.get(&i) {
            Some(poss_1) => {
              match pos_seq_set_pair_with_poss_4_backward_bps.1.get(&k) {
                Some(poss_2) => {
                  for &j in poss_1 {
                    for &l in poss_2 {
                      match lbap_mat.get(&(j, l)) {
                        Some(_) => {
                          backward_pos_pairs.push((j, l));
                        }, None => {},
                      }
                    }
                  }
                }, None => {},
              }
            }, None => {},
          }
          if backward_pos_pairs.len() > 0 {
            pos_pair_seqs_with_pos_pairs_4_backward_bpas.insert((i, k), backward_pos_pairs);
          }
        }
      }, None => {},
    }
  }
  let (log_sta_ppf_mat_4_bpas_1, log_sta_ppf_mat_4_bpas_2) = get_log_sta_ppf_mat_pair_4_bpas(&seq_len_pair, &lbpp_mat_pair, &lbap_mat, sta_scoring_params, &bpp_mat_pair, &not_bpp_mat_pair, &log_nbpp_mat_pair, &pos_pair_seqs_with_pos_span_pairs, &pos_seq_set_pair_with_poss_4_forward_bps, &pos_pair_seqs_with_pos_pairs_4_forward_bpas);
  get_lstapmq(&log_sta_ppf_mat_4_bpas_1, &log_sta_ppf_mat_4_bpas_2, &seq_len_pair, &lbpp_mat_pair, &lbap_mat, sta_scoring_params, &bpp_mat_pair, &not_bpp_mat_pair, &log_nbpp_mat_pair, &pos_seqs_with_poss_4_bas, &pos_pair_seqs_with_pos_span_pairs, &pos_pair_seqs_with_pos_pairs_4_forward_bpas, &pos_pair_seqs_with_pos_pairs_4_backward_bpas)
}

#[inline]
fn get_stapmq(lstapmq: &LogStapmq) -> Stapmq {
  Stapmq {
    base_pair_align_prob_mat: lstapmq.log_bpap_mat.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect(),
    base_align_prob_mat: lstapmq.log_bap_mat.iter().map(|(pos_pair, &lbap)| (*pos_pair, lbap.exp())).collect(),
    base_pair_indel_prob_mat_1: lstapmq.log_bpip_mat_1.iter().map(|(pos_pair, &lbpip)| (*pos_pair, lbpip.exp())).collect(),
    bpip_mat_2: lstapmq.log_bpip_mat_2.iter().map(|(pos_pair, &lbpip)| (*pos_pair, lbpip.exp())).collect(),
  }
}

#[inline]
pub fn get_log_sta_ppf_mat_pair_4_bpas(seq_len_pair: &(usize, usize), lbpp_mat_pair: &Ref2LogProbMatPair, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams, bpp_mat_pair: &SparseProbMatPair, not_bpp_mat_pair: &SparseProbMatPair, log_nbpp_mat_pair: &LogProbMatPair, pos_pair_seqs_with_pos_span_pairs: &PosPairSeqsWithPosSpanPairs, pos_seq_set_pair_with_poss_4_forward_bps: &PosSeqSetPairWithPoss, pos_pair_seqs_with_pos_pairs_4_forward_bpas: &PosPairSeqsWithPosPairs) -> LogPpf4dMatPair {
  let mut log_sta_ppf_mat_4_bpas_1 = LogPpf4dMat::default();
  let mut log_sta_ppf_mat_4_bpas_2 = log_sta_ppf_mat_4_bpas_1.clone();
  for substr_len_1 in 2 .. seq_len_pair.0 + 3 {
    for substr_len_2 in 2 .. seq_len_pair.1 + 3 {
      match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let j = i + substr_len_1 - 1;
            let l = k + substr_len_2 - 1;
            let log_sta_ppf_mat_4_2_loops = get_log_sta_ppf_mat_4_2_loops_1(&(i, j, k, l), lbap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_forward_bpas, &log_sta_ppf_mat_4_bpas_1);
            let mut log_sta_ppfs_4_2_loop_deletions_1 = vec![NEG_INFINITY; substr_len_1 - 1];
            let mut log_sta_ppfs_4_2_lds_2 = vec![NEG_INFINITY; substr_len_2 - 1];
            for n in i + 2 .. j {
              let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
              let mut max_ep_of_term_4_log_pf = log_sta_ppfs_4_2_loop_deletions_1[n - i - 1] + sta_scoring_params.loop_extending_gap_penalty;
              if max_ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
              }
              match pos_seq_set_pair_with_poss_4_forward_bps.0.get(&n) {
                Some(poss) => {
                  for &m in poss {
                    if m <= i {continue;}
                    match log_sta_ppf_mat_4_bpas_2.get(&(m, n, k, l)) {
                      Some (log_sta_pf_4_bpa) => {
                        let ep_of_term_4_log_pf = get_legp(&(i + 1, m - 1), sta_scoring_params) + get_logp(&(m, n), &lbpp_mat_pair.0, sta_scoring_params) + log_sta_pf_4_bpa;
                        if ep_of_term_4_log_pf.is_finite() {
                          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                          }
                        }
                      }, None => {},
                    }
                  }
                  if eps_of_terms_4_log_pf.len() > 0 {
                    log_sta_ppfs_4_2_loop_deletions_1[n - i] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
                  }
                }, None => {},
              }
            }
            for p in k + 2 .. l {
              let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
              let mut max_ep_of_term_4_log_pf = log_sta_ppfs_4_2_lds_2[p - k - 1] + sta_scoring_params.loop_extending_gap_penalty;
              if max_ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
              }
              match pos_seq_set_pair_with_poss_4_forward_bps.1.get(&p) {
                Some(poss) => {
                  for &o in poss {
                    if o <= k {continue;}
                    match log_sta_ppf_mat_4_bpas_2.get(&(i, j, o, p)) {
                      Some (log_sta_pf_4_bpa) => {
                        let ep_of_term_4_log_pf = get_legp(&(k + 1, o - 1), sta_scoring_params) + get_logp(&(o, p), &lbpp_mat_pair.1, sta_scoring_params) + log_sta_pf_4_bpa;
                        if ep_of_term_4_log_pf.is_finite() {
                          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                          }
                        }
                      }, None => {},
                    }
                  }
                  if eps_of_terms_4_log_pf.len() > 0 {
                    log_sta_ppfs_4_2_lds_2[p - k] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
                  }
                }, None => {},
              }
            }
            let eps_of_terms_4_log_pf = [
              log_sta_ppf_mat_4_2_loops[j - i - 1][l - k - 1],
              log_sta_ppfs_4_2_loop_deletions_1[j - i - 1],
              log_sta_ppfs_4_2_lds_2[l - k - 1],
            ];
            let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
            for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_4_bpas_2.insert((i, j, k, l), if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {NEG_INFINITY});
            log_sta_ppf_mat_4_bpas_1.insert((i, j, k, l), get_bpa_score(&(i, j, k, l), seq_len_pair, lbpp_mat_pair, lbap_mat, sta_scoring_params, bpp_mat_pair, not_bpp_mat_pair, log_nbpp_mat_pair) + log_sta_ppf_mat_4_bpas_2[&(i, j, k, l)]);
          }
        }, None => {},
      }
    }
  }
  (log_sta_ppf_mat_4_bpas_1, log_sta_ppf_mat_4_bpas_2)
}

#[inline]
fn get_log_sta_ppf_mat_4_2_loops_1(pos_quadruple: &PosQuadruple, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams, pos_pair_seqs_with_pos_pairs_4_forward_bpas: &PosPairSeqsWithPosPairs, log_sta_ppf_mat_4_bpas_1: &LogPpf4dMat) -> LogPpfMat {
  let &(i, j, k, l) = pos_quadruple;
  let substr_len_1 = j - i + 1;
  let substr_len_2 = l - k + 1;
  let mut log_sta_ppf_mat_4_2_loops = vec![vec![NEG_INFINITY; substr_len_2 - 1]; substr_len_1 - 1];
  let mut log_sta_ppf_mat_4_2_loops_and_bas = log_sta_ppf_mat_4_2_loops.clone();
  let mut log_sta_ppf_mat_4_2_loops_and_base_gaps_1 = log_sta_ppf_mat_4_2_loops.clone();
  let mut log_sta_ppf_mat_4_2_loops_and_bgs_2 = log_sta_ppf_mat_4_2_loops.clone();
  for n in i + 1 .. j {
    let bg_penalty = sta_scoring_params.base_opening_gap_penalty + get_begp(&(i + 2, n), sta_scoring_params);
    log_sta_ppf_mat_4_2_loops_and_base_gaps_1[n - i][0] = bg_penalty;
    log_sta_ppf_mat_4_2_loops[n - i][0] = log_sta_ppf_mat_4_2_loops_and_base_gaps_1[n - i][0];
  }
  for p in k + 1 .. l {
    let bgp = sta_scoring_params.base_opening_gap_penalty + get_begp(&(k + 2, p), sta_scoring_params);
    log_sta_ppf_mat_4_2_loops_and_bgs_2[0][p - k] = bgp;
    log_sta_ppf_mat_4_2_loops[0][p - k] = log_sta_ppf_mat_4_2_loops_and_bgs_2[0][p - k];
  }
  for n in i + 1 .. j {
    for p in k + 1 .. l {
      match lbap_mat.get(&(n, p)) {
        Some(_) => {
          let ba_score = get_ba_score(&(n, p), lbap_mat, sta_scoring_params);
          if n == i + 1 && p == k + 1 {
            log_sta_ppf_mat_4_2_loops_and_bas[n - i][p - k] = ba_score;
          } else {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops[n - i - 1][p - k - 1] + ba_score;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            match pos_pair_seqs_with_pos_pairs_4_forward_bpas.get(&(n, p)) {
              Some(pos_pairs) => {
                for &(m, o) in pos_pairs {
                  if m <= i || o <= k {continue;}
                  let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops[m - i - 1][o - k - 1] + log_sta_ppf_mat_4_bpas_1[&(m, n, o, p)];
                  if ep_of_term_4_log_pf.is_finite() {
                    eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                      max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    }
                  }
                }
              }, None => {},
            }
            if eps_of_terms_4_log_pf.len() > 0 {
              log_sta_ppf_mat_4_2_loops_and_bas[n - i][p - k] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
            }
          }
        }, None => {},
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bas[n - i - 1][p - k] + sta_scoring_params.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bgs_2[n - i - 1][p - k] + sta_scoring_params.base_opening_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_base_gaps_1[n - i - 1][p - k] + sta_scoring_params.base_extending_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_sta_ppf_mat_4_2_loops_and_base_gaps_1[n - i][p - k] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bas[n - i][p - k - 1] + sta_scoring_params.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_base_gaps_1[n - i][p - k - 1] + sta_scoring_params.base_opening_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bgs_2[n - i][p - k - 1] + sta_scoring_params.base_extending_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_sta_ppf_mat_4_2_loops_and_bgs_2[n - i][p - k] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let eps_of_terms_4_log_pf = [
        log_sta_ppf_mat_4_2_loops_and_bas[n - i][p - k],
        log_sta_ppf_mat_4_2_loops_and_base_gaps_1[n - i][p - k],
        log_sta_ppf_mat_4_2_loops_and_bgs_2[n - i][p - k],
      ];
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_ppf_mat_4_2_loops[n - i][p - k] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
    }
  }
  log_sta_ppf_mat_4_2_loops
}

#[inline]
fn get_begp(pos_pair: &PosPair, sta_scoring_params: &StaScoringParams) -> StaScore {
  (pos_pair.1 + 1 - pos_pair.0) as StaScore * sta_scoring_params.base_extending_gap_penalty
}

#[inline]
fn get_logp(pos_pair: &PosPair, lbpp_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams) -> StaScore {
  get_bp_score(pos_pair, lbpp_mat, sta_scoring_params) + 2. * sta_scoring_params.loop_opening_gap_penalty
}

#[inline]
fn get_legp(pos_pair: &PosPair, sta_scoring_params: &StaScoringParams) -> StaScore {
  (pos_pair.1 + 1 - pos_pair.0) as StaScore * sta_scoring_params.loop_extending_gap_penalty
}

#[inline]
fn get_ba_score(pos_pair: &PosPair, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams) -> StaScore {
  lbap_mat[pos_pair] - sta_scoring_params.log_nh_bap
}

#[inline]
fn get_bp_score(pos_pair: &PosPair, lbpp_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams) -> StaScore {
  lbpp_mat[pos_pair] - sta_scoring_params.log_null_hypothesis_bpp
}

#[inline]
fn get_bpa_score(pos_quadruple: &PosQuadruple, seq_len_pair: &(usize, usize), lbpp_mat_pair: &Ref2LogProbMatPair, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams, bpp_mat_pair: &SparseProbMatPair, not_bpp_mat_pair: &SparseProbMatPair, log_nbpp_mat_pair: &LogProbMatPair) -> StaScore {
  if *pos_quadruple == (0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1) {0.} else {
    let bp_pos_pair_1 = (pos_quadruple.0, pos_quadruple.1);
    let bp_pos_pair_2 = (pos_quadruple.2, pos_quadruple.3);
    let aligned_pos_pair_1 = (pos_quadruple.0, pos_quadruple.2);
    let aligned_pos_pair_2 = (pos_quadruple.1, pos_quadruple.3);
    let bin_prob_dist_1 = [bpp_mat_pair.0[&bp_pos_pair_1], not_bpp_mat_pair.0[&bp_pos_pair_1]];
    let bin_prob_dist_2 = [bpp_mat_pair.1[&bp_pos_pair_2], not_bpp_mat_pair.1[&bp_pos_pair_2]];
    let bin_log_prob_dist_1 = [lbpp_mat_pair.0[&bp_pos_pair_1], log_nbpp_mat_pair.0[&bp_pos_pair_1]];
    let bin_log_prob_dist_2 = [lbpp_mat_pair.1[&bp_pos_pair_2], log_nbpp_mat_pair.1[&bp_pos_pair_2]];
    get_bp_score(&bp_pos_pair_1, &lbpp_mat_pair.0, sta_scoring_params)
    + get_bp_score(&bp_pos_pair_2, &lbpp_mat_pair.1, sta_scoring_params)
    + get_ba_score(&aligned_pos_pair_1, lbap_mat, sta_scoring_params)
    + get_ba_score(&aligned_pos_pair_2, lbap_mat, sta_scoring_params)
    + sta_scoring_params.scale_param_4_bpa_score * (sta_scoring_params.offset_bpa_score - get_jensen_shannon_dist(&bin_prob_dist_1[..], &bin_prob_dist_2[..], &bin_log_prob_dist_1[..], &bin_log_prob_dist_2[..]))
  }
}

#[inline]
fn get_jensen_shannon_dist(prob_dist_1: ProbDistSlice, prob_dist_2: ProbDistSlice, log_prob_dist_1: LogProbDistSlice, log_prob_dist_2: LogProbDistSlice) -> ProbDistDist {
  let mid_prob_dist = prob_dist_1.iter().zip(prob_dist_2).map(|(&prob_1, &prob_2)| (prob_1 + prob_2) / 2.).collect::<ProbDist>();
  let log_mid_prob_dist = mid_prob_dist.iter().map(|&mid_prob| mid_prob.log2()).collect::<LogProbDist>();
  let jsd = ((get_kullback_leibler_div(prob_dist_1, log_prob_dist_1, &log_mid_prob_dist)
  + get_kullback_leibler_div(prob_dist_2, log_prob_dist_2, &log_mid_prob_dist)) / 2.).sqrt();
  if jsd.is_finite() {jsd} else {0.}
}

#[inline]
fn get_kullback_leibler_div(prob_dist: ProbDistSlice, log_prob_dist_1: LogProbDistSlice, log_prob_dist_2: LogProbDistSlice) -> ProbDistDiv {
  multizip((prob_dist.iter(), log_prob_dist_1.iter(), log_prob_dist_2.iter())).map(|(&prob, &log_prob_1, &log_prob_2)| prob * (log_prob_1 - log_prob_2)).sum()
}

#[inline]
fn get_lstapmq(log_sta_ppf_mat_4_bpas_1: &LogPpf4dMat, log_sta_ppf_mat_4_bpas_2: &LogPpf4dMat, seq_len_pair: &(usize, usize), lbpp_mat_pair: &Ref2LogProbMatPair, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams, bpp_mat_pair: &SparseProbMatPair, not_bpp_mat_pair: &SparseProbMatPair, log_nbpp_mat_pair: &LogProbMatPair, pos_seqs_with_poss_4_bas: &PosSeqsWithPoss, pos_pair_seqs_with_pos_span_pairs: &PosPairSeqsWithPosSpanPairs, pos_pair_seqs_with_pos_pairs_4_forward_bpas: &PosPairSeqsWithPosPairs, pos_pair_seqs_with_pos_pairs_4_backward_bpas: &PosPairSeqsWithPosPairs) -> LogStapmq {
  let mut lstapmq = LogStapmq::new();
  let pos_quadruple_4_2_external_loops = (0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1);
  let log_sta_pf_4_bpa = log_sta_ppf_mat_4_bpas_1[&pos_quadruple_4_2_external_loops];
  let log_sta_ppf_mat_4_2_els_1 = get_log_sta_ppf_mat_4_2_loops_1(&pos_quadruple_4_2_external_loops, lbap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_forward_bpas, log_sta_ppf_mat_4_bpas_1);
  let log_sta_ppf_mat_4_2_els_2 = get_log_sta_ppf_mat_4_2_loops_2(&pos_quadruple_4_2_external_loops, lbap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_backward_bpas, log_sta_ppf_mat_4_bpas_1);
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_quadruples = SeqsOfEpsOfTerms4LogProbsWithPosQuadruples::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_quadruples = EpsOfTerms4LogProbsWithPosQuadruples::default();
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_baps = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs_4_baps = EpsOfTerms4LogProbsWithPosPairs::default();
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1 = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1 = EpsOfTerms4LogProbsWithPosPairs::default();
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2 = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2 = EpsOfTerms4LogProbsWithPosPairs::default();
  for substr_len_1 in (2 .. seq_len_pair.0 + 1).rev() {
    for substr_len_2 in (2 .. seq_len_pair.1 + 1).rev() {
      match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let j = i + substr_len_1 - 1;
            let l = k + substr_len_2 - 1;
            let ep_of_term_4_log_prob = log_sta_ppf_mat_4_2_els_1[i - 1][k - 1] + log_sta_ppf_mat_4_bpas_1[&(i, j, k, l)] + log_sta_ppf_mat_4_2_els_2[j - 1 + 1][l - 1 + 1] - log_sta_pf_4_bpa;
            seqs_of_eps_of_terms_4_log_probs_with_pos_quadruples.insert((i, j, k, l), vec![ep_of_term_4_log_prob]);
            max_eps_of_terms_4_log_probs_with_pos_quadruples.insert((i, j, k, l), ep_of_term_4_log_prob);
          }
        }, None => {},
      }
    }
  }
  for i in 1 .. seq_len_pair.0 + 1 {
    match pos_seqs_with_poss_4_bas.get(&i) {
      Some(poss) => {
        for &k in poss {
          let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
          let mut max_ep_of_term_4_log_prob = log_sta_ppf_mat_4_2_els_1[i - 1][k - 1]
          + get_ba_score(&(i, k), lbap_mat, sta_scoring_params)
          + log_sta_ppf_mat_4_2_els_2[i - 1 + 1][k - 1 + 1]
          - log_sta_pf_4_bpa;
          if max_ep_of_term_4_log_prob.is_finite() {
            eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
          }
          match pos_pair_seqs_with_pos_pairs_4_forward_bpas.get(&(i, k)) {
            Some(pos_pairs) => {
              for &(j, l) in pos_pairs {
                let ep_of_term_4_log_prob = log_sta_ppf_mat_4_2_els_1[j - 1][l - 1]
                + log_sta_ppf_mat_4_bpas_1[&(j, i, l, k)]
                + log_sta_ppf_mat_4_2_els_2[i - 1 + 1][k - 1 + 1]
                - log_sta_pf_4_bpa;
                if ep_of_term_4_log_prob.is_finite() {
                  eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                  if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
                    max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
                }
              }
            }, None => {},
          }
          match pos_pair_seqs_with_pos_pairs_4_backward_bpas.get(&(i, k)) {
            Some(pos_pairs) => {
              for &(j, l) in pos_pairs {
                let ep_of_term_4_log_prob = log_sta_ppf_mat_4_2_els_1[i - 1][k - 1]
                + log_sta_ppf_mat_4_bpas_1[&(i, j, k, l)]
                + log_sta_ppf_mat_4_2_els_2[j - 1 + 1][l - 1 + 1]
                - log_sta_pf_4_bpa;
                if ep_of_term_4_log_prob.is_finite() {
                  eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                  if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
                    max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
                }
              }
            }, None => {},
          }
          seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_baps.insert((i, k), eps_of_terms_4_log_prob);
          max_eps_of_terms_4_log_probs_with_pos_pairs_4_baps.insert((i, k), max_ep_of_term_4_log_prob);
        }
      }, None => {},
    }
  }
  for substr_len_1 in (3 .. seq_len_pair.0 + 1).rev() {
    for substr_len_2 in (3 .. seq_len_pair.1 + 1).rev() {
      match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let j = i + substr_len_1 - 1;
            let l = k + substr_len_2 - 1;
            lstapmq.log_bpap_mat.insert((i, j, k, l), get_valid_log_prob(logsumexp(&seqs_of_eps_of_terms_4_log_probs_with_pos_quadruples[&(i, j, k, l)], max_eps_of_terms_4_log_probs_with_pos_quadruples[&(i, j, k, l)])));
            seqs_of_eps_of_terms_4_log_probs_with_pos_quadruples.remove(&(i, j, k, l));
            let lbpap = lstapmq.log_bpap_mat[&(i, j, k, l)];
            let log_sta_pf_4_bpa = log_sta_ppf_mat_4_bpas_1[&(i, j, k, l)];
            let sum_of_terms_4_ep_of_term_4_log_prob = lbpap
            + get_bpa_score(&(i, j, k, l), seq_len_pair, lbpp_mat_pair, lbap_mat, sta_scoring_params, bpp_mat_pair, not_bpp_mat_pair, log_nbpp_mat_pair)
            - log_sta_pf_4_bpa;
            let log_sta_ppf_mat_4_2_loops_1 = get_log_sta_ppf_mat_4_2_loops_1(&(i, j, k, l), lbap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_forward_bpas, log_sta_ppf_mat_4_bpas_1);
            let log_sta_ppf_mat_4_2_loops_2 = get_log_sta_ppf_mat_4_2_loops_2(&(i, j, k, l), lbap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_backward_bpas, log_sta_ppf_mat_4_bpas_1);
            for substr_len_3 in (2 .. substr_len_1 - 1).rev() {
              for substr_len_4 in (2 .. substr_len_2 - 1).rev() {
                match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_3, substr_len_4)) {
                  Some(pos_pairs) => {
                    for &(m, o) in pos_pairs {
                      let n = m + substr_len_3 - 1;
                      let p = o + substr_len_4 - 1;
                      if m <= i || n >= j || o <= k || p >= l {continue;}
                      let ep_of_term_4_log_prob = sum_of_terms_4_ep_of_term_4_log_prob
                      + log_sta_ppf_mat_4_2_loops_1[m - i - 1][o - k - 1]
                      + log_sta_ppf_mat_4_bpas_1[&(m, n, o, p)]
                      + log_sta_ppf_mat_4_2_loops_2[n - (i + 1) + 1][p - (k + 1) + 1];
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_log_probs_with_pos_quadruples.get_mut(&(m, n, o, p)).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_quadruples.get_mut(&(m, n, o, p)).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob > *max_ep_of_term_4_log_prob {
                          *max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                  }, None => {},
                }
              }
            }
            for m in i + 1 .. j {
              match pos_seqs_with_poss_4_bas.get(&m) {
                Some(poss) => {
                  for &o in poss {
                    if o <= k || l <= o {continue;}
                    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
                    let mut max_ep_of_term_4_log_prob_1 = sum_of_terms_4_ep_of_term_4_log_prob
                    + log_sta_ppf_mat_4_2_loops_1[m - i - 1][o - k - 1]
                    + get_ba_score(&(m, o), lbap_mat, sta_scoring_params)
                    + log_sta_ppf_mat_4_2_loops_2[m - (i + 1) + 1][o - (k + 1) + 1];
                    if max_ep_of_term_4_log_prob_1.is_finite() {
                      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob_1);
                    }
                    if m > i + 1 && o > k + 1 {
                      match pos_pair_seqs_with_pos_pairs_4_forward_bpas.get(&(m, o)) {
                        Some(pos_pairs) => {
                          for &(n, p) in pos_pairs {
                            if n <= i || p <= k {continue;}
                            let ep_of_term_4_log_prob = sum_of_terms_4_ep_of_term_4_log_prob
                            + log_sta_ppf_mat_4_2_loops_1[n - i - 1][p - k - 1]
                            + log_sta_ppf_mat_4_bpas_1[&(n, m, p, o)]
                            + log_sta_ppf_mat_4_2_loops_2[m - (i + 1) + 1][o - (k + 1) + 1];
                            if ep_of_term_4_log_prob.is_finite() {
                              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                              if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob_1 {
                                max_ep_of_term_4_log_prob_1 = ep_of_term_4_log_prob;
                              }
                            }
                          }
                        }, None => {},
                      }
                    }
                    if m < j - 1 && o < l - 1 {
                      match pos_pair_seqs_with_pos_pairs_4_backward_bpas.get(&(m, o)) {
                        Some(pos_pairs) => {
                          for &(n, p) in pos_pairs {
                            if n >= j || p >= l {continue;}
                            let ep_of_term_4_log_prob = sum_of_terms_4_ep_of_term_4_log_prob
                            + log_sta_ppf_mat_4_2_loops_1[m - i - 1][o - k - 1]
                            + log_sta_ppf_mat_4_bpas_1[&(m, n, o, p)]
                            + log_sta_ppf_mat_4_2_loops_2[n - (i + 1) + 1][p - (k + 1) + 1];
                            if ep_of_term_4_log_prob.is_finite() {
                              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                              if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob_1 {
                                max_ep_of_term_4_log_prob_1 = ep_of_term_4_log_prob;
                              }
                            }
                          }
                        }, None => {},
                      }
                    }
                    if eps_of_terms_4_log_prob.len() > 0 {
                      seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_baps.get_mut(&(m, o)).expect("Failed to get an element of a hash map.").extend_from_slice(&eps_of_terms_4_log_prob);
                      let max_ep_of_term_4_log_prob_2 = max_eps_of_terms_4_log_probs_with_pos_pairs_4_baps.get_mut(&(m, o)).expect("Failed to get an element of a hash map.");
                      if max_ep_of_term_4_log_prob_1 > *max_ep_of_term_4_log_prob_2 {
                        *max_ep_of_term_4_log_prob_2 = max_ep_of_term_4_log_prob_1;
                      }
                    }
                  }
                }, None => {},
              }
            }
            for substr_len_3 in (2 .. substr_len_1 - 1).rev() {
              match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_3, substr_len_2)) {
                Some(pos_pairs) => {
                  for &(m, o) in pos_pairs {
                    let n = m + substr_len_3 - 1;
                    let p = o + substr_len_2 - 1;
                    if m <= i || n >= j || o != k || p != l {continue;}
                    let max_ep_of_term_4_log_prob_exists = match max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1.get(&(m, n)) {
                      Some(_) => {true},
                      None => {false},
                    };
                    let ep_of_term_4_log_prob = lbpap
                    + get_legp(&(i + 1, m - 1), sta_scoring_params)
                    + get_logp(&(m, n), &lbpp_mat_pair.0, sta_scoring_params)
                    + log_sta_ppf_mat_4_bpas_2[&(m, n, k, l)]
                    + get_legp(&(n + 1, j - 1), sta_scoring_params)
                    - log_sta_pf_4_bpa;
                    if max_ep_of_term_4_log_prob_exists {
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1.get_mut(&(m, n)).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1.get_mut(&(m, n)).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob > *max_ep_of_term_4_log_prob {
                          *max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    } else {
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1.insert((m, n), vec![ep_of_term_4_log_prob]);
                        max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1.insert((m, n), ep_of_term_4_log_prob);
                      }
                    }
                  }
                }, None => {},
              }
            }
            for substr_len_3 in (2 .. substr_len_2 - 1).rev() {
              match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_1, substr_len_3)) {
                Some(pos_pairs) => {
                  for &(m, o) in pos_pairs {
                    let n = m + substr_len_1 - 1;
                    let p = o + substr_len_3 - 1;
                    if o <= k || p >= l || m != i || n != j {continue;}
                    let max_ep_of_term_4_log_prob_exists = match max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2.get(&(o, p)) {
                      Some(_) => {true},
                      None => {false},
                    };
                    let ep_of_term_4_log_prob = lbpap
                    + get_legp(&(k + 1, o - 1), sta_scoring_params)
                    + get_logp(&(o, p), &lbpp_mat_pair.1, sta_scoring_params)
                    + log_sta_ppf_mat_4_bpas_2[&(i, j, o, p)]
                    + get_legp(&(p + 1, l - 1), sta_scoring_params)
                    - log_sta_pf_4_bpa;
                    if max_ep_of_term_4_log_prob_exists {
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2.get_mut(&(o, p)).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2.get_mut(&(o, p)).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob > *max_ep_of_term_4_log_prob {
                          *max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    } else {
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2.insert((o, p), vec![ep_of_term_4_log_prob]);
                        max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2.insert((o, p), ep_of_term_4_log_prob);
                      }
                    }
                  }
                }, None => {},
              }
            }
          }
        }, None => {},
      }
    }
  }
  for i in 1 .. seq_len_pair.0 + 1 {
    match pos_seqs_with_poss_4_bas.get(&i) {
      Some(poss) => {
        for &k in poss {
          lstapmq.log_bap_mat.insert((i, k), get_valid_log_prob(logsumexp(&seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_baps[&(i, k)], max_eps_of_terms_4_log_probs_with_pos_pairs_4_baps[&(i, k)])));
          seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_baps.remove(&(i, k));
        }
      }, None => {},
    }
  }
  for (pos_pair, &max_ep_of_term_4_log_prob) in &max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1 {
    let eps_of_terms_4_log_prob = seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_1.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
    lstapmq.log_bpip_mat_1.insert(*pos_pair, get_valid_log_prob(logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob)));
     eps_of_terms_4_log_prob.drain(..);
  }
  for (pos_pair, &max_ep_of_term_4_log_prob) in &max_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2 {
    let eps_of_terms_4_log_prob = seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_4_bpips_2.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
    lstapmq.log_bpip_mat_2.insert(*pos_pair, get_valid_log_prob(logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob)));
     eps_of_terms_4_log_prob.drain(..);
  }
  lstapmq
}

#[inline]
fn get_log_sta_ppf_mat_4_2_loops_2(pos_quadruple: &PosQuadruple, lbap_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams, pos_pair_seqs_with_pos_pairs_4_backward_bpas: &PosPairSeqsWithPosPairs, log_sta_ppf_mat_4_bpas_1: &LogPpf4dMat) -> LogPpfMat {
  let &(i, j, k, l) = pos_quadruple;
  let substr_len_1 = j - i + 1;
  let substr_len_2 = l - k + 1;
  let mut log_sta_ppf_mat_4_2_loops = vec![vec![NEG_INFINITY; substr_len_2 - 1]; substr_len_1 - 1];
  let mut log_sta_ppf_mat_4_2_loops_and_bas = log_sta_ppf_mat_4_2_loops.clone();
  let mut log_sta_ppf_mat_4_2_loops_and_bgs_1 = log_sta_ppf_mat_4_2_loops.clone();
  let mut log_sta_ppf_mat_4_2_loops_and_bgs_2 = log_sta_ppf_mat_4_2_loops.clone();
  for m in i + 1 .. j {
    let bgp = sta_scoring_params.base_opening_gap_penalty + get_begp(&(m + 1, j - 1), sta_scoring_params);
    log_sta_ppf_mat_4_2_loops_and_bgs_1[m - (i + 1)][l - (k + 1)] = bgp;
    log_sta_ppf_mat_4_2_loops[m - (i + 1)][l - (k + 1)] = log_sta_ppf_mat_4_2_loops_and_bgs_1[m - (i + 1)][l - (k + 1)];
  }
  for o in k + 1 .. l {
    let bgp = sta_scoring_params.base_opening_gap_penalty + get_begp(&(o + 1, l - 1), sta_scoring_params);
    log_sta_ppf_mat_4_2_loops_and_bgs_2[j - (i + 1)][o - (k + 1)] = bgp;
    log_sta_ppf_mat_4_2_loops[j - (i + 1)][o - (k + 1)] = log_sta_ppf_mat_4_2_loops_and_bgs_2[j - (i + 1)][o - (k + 1)];
  }
  for m in (i + 1 .. j).rev() {
    for o in (k + 1 .. l).rev() {
      match lbap_mat.get(&(m, o)) {
        Some(_) => {
          let ba_score = get_ba_score(&(m, o), lbap_mat, sta_scoring_params);
          if m == j - 1 && o == l - 1 {
            log_sta_ppf_mat_4_2_loops_and_bas[m - (i + 1)][o - (k + 1)] = ba_score;
          } else {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops[m - (i + 1) + 1][o - (k + 1) + 1] + ba_score;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            match pos_pair_seqs_with_pos_pairs_4_backward_bpas.get(&(m, o)) {
              Some(pos_pairs) => {
                for &(n, p) in pos_pairs {
                  if n >= j || p >= l {continue;}
                  let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops[n - (i + 1) + 1][p - (k + 1) + 1] + log_sta_ppf_mat_4_bpas_1[&(m, n, o, p)];
                  if ep_of_term_4_log_pf.is_finite() {
                    eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                      max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    }
                  }
                }
              }, None => {},
            }
            if eps_of_terms_4_log_pf.len() > 0 {
              log_sta_ppf_mat_4_2_loops_and_bas[m - (i + 1)][o - (k + 1)] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
            }
          }
        }, None => {},
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bas[m - (i + 1)][o - (k + 1)] + sta_scoring_params.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bgs_2[m - (i + 1) + 1][o - (k + 1)] + sta_scoring_params.base_opening_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bgs_1[m - (i + 1) + 1][o - (k + 1)] + sta_scoring_params.base_extending_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_sta_ppf_mat_4_2_loops_and_bgs_1[m - (i + 1)][o - (k + 1)] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bas[m - (i + 1)][o - (k + 1) + 1] + sta_scoring_params.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bgs_1[m - (i + 1)][o - (k + 1) + 1] + sta_scoring_params.base_opening_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      let ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bgs_2[m - (i + 1)][o - (k + 1) + 1] + sta_scoring_params.base_extending_gap_penalty;
      if ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      if eps_of_terms_4_log_pf.len() > 0 {
        log_sta_ppf_mat_4_2_loops_and_bgs_2[m - (i + 1)][o - (k + 1)] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
      let eps_of_terms_4_log_pf = [
        log_sta_ppf_mat_4_2_loops_and_bas[m - (i + 1)][o - (k + 1)],
        log_sta_ppf_mat_4_2_loops_and_bgs_1[m - (i + 1)][o - (k + 1)],
        log_sta_ppf_mat_4_2_loops_and_bgs_2[m - (i + 1)][o - (k + 1)],
      ];
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
        if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
          max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_ppf_mat_4_2_loops[m - (i + 1)][o - (k + 1)] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
      }
    }
  }
  log_sta_ppf_mat_4_2_loops
}

#[inline]
pub fn prob_cons_transformation_of_lstapmq(lstapmqs_with_rna_id_pairs: &LstapmqsWithRnaIdPairs, rna_id_pair: &RnaIdPair, num_of_rnas: usize) -> LogStapmq {
  let mut lstapmq = lstapmqs_with_rna_id_pairs[rna_id_pair].clone();
  let log_coefficient = -fast_ln((num_of_rnas - 1) as Prob);
  for (pos_quadruple, lbpap) in lstapmq.log_bpap_mat.iter_mut() {
    let (i, j, k, l) = *pos_quadruple;
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = *lbpap;
    if max_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      for (&(m, n, o, p), &lbpap_1) in lstapmqs_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].log_bpap_mat.iter() {
        if (rna_id_pair.0 < rna_id && m != i && n != j) || (rna_id_pair.0 > rna_id && o != i && p != j) {continue;}
        for (&(m, n, o, p), &lbpap_2) in lstapmqs_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].log_bpap_mat.iter() {
          if (rna_id_pair.1 < rna_id && m != k && n != l) || (rna_id_pair.1 > rna_id && o != k && p != l) {continue;}
          let ep_of_term_4_log_prob = lbpap_1 + lbpap_2;
          if ep_of_term_4_log_prob.is_finite() {
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
              max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
            }
          }
        }
      }
    }
    if eps_of_terms_4_log_prob.len() > 0 {
      *lbpap = get_valid_log_prob(log_coefficient + logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob));
    }
  }
  for (pos_pair, lbap) in lstapmq.log_bap_mat.iter_mut() {
    let (i, j) = *pos_pair;
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = *lbap;
    if max_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      for (&(k, l), &lbap_1) in lstapmqs_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].log_bap_mat.iter() {
        if (rna_id_pair.0 < rna_id && k != i) || (rna_id_pair.0 > rna_id && l != i) {continue;}
        for (&(k, l), &lbap_2) in lstapmqs_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].log_bap_mat.iter() {
          if (rna_id_pair.1 < rna_id && k != j) || (rna_id_pair.1 > rna_id && l != j) {continue;}
          let ep_of_term_4_log_prob = lbap_1 + lbap_2;
          if ep_of_term_4_log_prob.is_finite() {
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
              max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
            }
          }
        }
      }
    }
    if eps_of_terms_4_log_prob.len() > 0 {
      *lbap = get_valid_log_prob(log_coefficient + logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob));
    }
  }
  lstapmq
}

#[inline]
pub fn get_stapmqs_with_rna_id_pairs(lstapmqs_with_rna_id_pairs: &LstapmqsWithRnaIdPairs) -> StapmqsWithRnaIdPairs {
  let mut stapmqs_with_rna_id_pairs = StapmqsWithRnaIdPairs::default();
  for (rna_id_pair, lstapmq) in lstapmqs_with_rna_id_pairs.iter() {
    stapmqs_with_rna_id_pairs.insert(*rna_id_pair, get_stapmq(lstapmq));
  }
  stapmqs_with_rna_id_pairs
}

#[inline]
pub fn pct_of_lbap_mat(lbap_mats_with_rna_id_pairs: &LogProbMatsWithRnaIdPairs, rna_id_pair: &RnaIdPair, num_of_rnas: usize) -> SparseLogProbMat {
  let mut lbap_mat = lbap_mats_with_rna_id_pairs[rna_id_pair].clone();
  let log_coefficient = -((num_of_rnas - 1) as Prob).log2();
  for (&(i, j), lbap) in lbap_mat.iter_mut() {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = *lbap;
    if max_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      for (&(k, l), &lbap_1) in lbap_mats_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].iter() {
        if (rna_id_pair.0 < rna_id && k != i) || (rna_id_pair.0 > rna_id && l != i) {continue;}
        for (&(k, l), &lbap_2) in lbap_mats_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].iter() {
          if (rna_id_pair.1 < rna_id && k != j) || (rna_id_pair.1 > rna_id && l != j) {continue;}
          let ep_of_term_4_log_prob = lbap_1 + lbap_2;
          if ep_of_term_4_log_prob.is_finite() {
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
              max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
            }
          }
        }
      }
    }
    *lbap = get_valid_log_prob(log_coefficient + logsumexp2(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob));
  }
  lbap_mat
}

#[inline]
pub fn remove_little_lbpps_from_sparse_lbpp_mat(lbpp_mat: &SparseLogProbMat, min_lbpp: LogProb) -> SparseLogProbMat {
  lbpp_mat.iter().filter(|&(_, &lbpp)| {lbpp >= min_lbpp && lbpp > NEG_INFINITY}).map(|(pos_pair, &lbpp)| {(*pos_pair, lbpp)}).collect::<SparseLogProbMat>()
}

#[inline]
pub fn get_sparse_lbap_mat(lbap_mat: &LogProbMat, min_lbap: LogProb) -> SparseLogProbMat {
  let mut sparse_lbap_mat = SparseLogProbMat::default();
  let lbap_mat_dims = (lbap_mat.len(), lbap_mat[0].len());
  for i in 0 .. lbap_mat_dims.0 {
    for j in 0 .. lbap_mat_dims.1 {
      let lbap = get_valid_log_prob(lbap_mat[i][j] / LN_2);
      if lbap >= min_lbap && lbap > NEG_INFINITY {
        sparse_lbap_mat.insert((i + 1, j + 1), lbap);
      }
    }
  }
  sparse_lbap_mat
}

#[inline]
pub fn remove_little_lbaps_from_sparse_lbap_mat(lbap_mat: &SparseLogProbMat, min_lbap: LogProb) -> SparseLogProbMat {
  lbap_mat.iter().filter(|&(_, &lbap)| {lbap >= min_lbap && lbap > NEG_INFINITY}).map(|(pos_pair, &lbap)| {(*pos_pair, lbap)}).collect::<SparseLogProbMat>()
}

#[inline]
pub fn convert_log_base_of_sparse_lbap_mat(lbap_mat: &SparseLogProbMat) -> SparseLogProbMat {
  lbap_mat.iter().map(|(pos_pair, &lbap)| {(*pos_pair, lbap / LN_2)}).collect::<SparseLogProbMat>()
}

#[inline]
pub fn is_rna_base(base: Base) -> bool {
  match base {
    A => true,
    U => true,
    G => true,
    C => true,
    _ => false,
  }
}

#[inline]
pub fn get_lbpp_mat_pair(lstapmqs_with_rna_id_pairs: &LstapmqsWithRnaIdPairs, rna_id_pair: &RnaIdPair) -> LogProbMatPair {
  let ref lstapmq = lstapmqs_with_rna_id_pairs[rna_id_pair];
  let mut lbpp_mat_pair = (SparseLogProbMat::default(), SparseLogProbMat::default());
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_1 = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_2 = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs_1 = EpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs_2 = EpsOfTerms4LogProbsWithPosPairs::default();
  for (pos_pair, &lbpip) in &lstapmq.log_bpip_mat_1 {
    seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_1.insert(*pos_pair, vec![lbpip]);
    max_eps_of_terms_4_log_probs_with_pos_pairs_1.insert(*pos_pair, lbpip);
  }
  for (pos_pair, &lbpip) in &lstapmq.log_bpip_mat_2 {
    seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_2.insert(*pos_pair, vec![lbpip]);
    max_eps_of_terms_4_log_probs_with_pos_pairs_2.insert(*pos_pair, lbpip);
  }
  for (&(i, j, k, l), &lbpap) in &lstapmq.log_bpap_mat {
    if lbpap.is_finite() {
      let max_ep_of_term_4_log_prob_exists = match max_eps_of_terms_4_log_probs_with_pos_pairs_1.get(&(i, j)) {
        Some(_) => {true},
        None => {false},
      };
      if max_ep_of_term_4_log_prob_exists {
        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_1.get_mut(&(i, j)).expect("Failed to get an element of a hash map.").push(lbpap);
        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs_1.get_mut(&(i, j)).expect("Failed to get an element of a hash map.");
        if lbpap > *max_ep_of_term_4_log_prob {
          *max_ep_of_term_4_log_prob = lbpap;
        }
      } else {
        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_1.insert((i, j), vec![lbpap]);
        max_eps_of_terms_4_log_probs_with_pos_pairs_1.insert((i, j), lbpap);
      }
      let max_ep_of_term_4_log_prob_exists = match max_eps_of_terms_4_log_probs_with_pos_pairs_2.get(&(k, l)) {
        Some(_) => {true},
        None => {false},
      };
      if max_ep_of_term_4_log_prob_exists {
        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_2.get_mut(&(k, l)).expect("Failed to get an element of a hash map.").push(lbpap);
        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs_2.get_mut(&(k, l)).expect("Failed to get an element of a hash map.");
        if lbpap > *max_ep_of_term_4_log_prob {
          *max_ep_of_term_4_log_prob = lbpap;
        }
      } else {
        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_2.insert((k, l), vec![lbpap]);
        max_eps_of_terms_4_log_probs_with_pos_pairs_2.insert((k, l), lbpap);
      }
    }
  }
  for (pos_pair, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_1.iter() {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs_1[pos_pair];
    let lbpp = get_valid_log_prob(logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob));
    lbpp_mat_pair.0.insert(*pos_pair, lbpp);
  }
  for (pos_pair, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_pos_pairs_2.iter() {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs_2[pos_pair];
    let lbpp = get_valid_log_prob(logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob));
    lbpp_mat_pair.1.insert(*pos_pair, lbpp);
  }
  lbpp_mat_pair
}

#[inline]
pub fn pct_of_lbpp_mat(lbpp_mat_pairs_with_rna_id_pairs: &LogProbMatPairsWithRnaIdPairs, rna_id_1: RnaId, num_of_rnas: usize) -> SparseLogProbMat {
  let mut lbpp_mat = SparseLogProbMat::default();
  let log_coefficient = - fast_ln((num_of_rnas - 1) as Prob);
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs = EpsOfTerms4LogProbsWithPosPairs::default();
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id_1 == rna_id_2 {continue;}
    let rna_id_pair = if rna_id_1 < rna_id_2 {(rna_id_1, rna_id_2)} else {(rna_id_2, rna_id_1)};
    let ref lbpp_mat_pair = lbpp_mat_pairs_with_rna_id_pairs[&rna_id_pair];
    for (pos_pair, &lbpp) in if rna_id_1 < rna_id_2 {lbpp_mat_pair.0.iter()} else {lbpp_mat_pair.1.iter()} {
      let max_ep_of_term_4_log_prob_exists = match max_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(&pos_pair) {
        Some(_) => {true},
        None => {false},
      };
      if max_ep_of_term_4_log_prob_exists {
        let eps_of_terms_4_log_prob = seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
        if lbpp.is_finite() {
          eps_of_terms_4_log_prob.push(lbpp);
          let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
          if lbpp > *max_ep_of_term_4_log_prob {
            *max_ep_of_term_4_log_prob = lbpp;
          }
        }
      } else {
        seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, vec![lbpp]);
        max_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, lbpp);
      }
    }
  }
  for (pos_pair, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.iter() {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs[pos_pair];
    let lbpp = get_valid_log_prob((log_coefficient + logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob)) / LN_2);
    lbpp_mat.insert(*pos_pair, lbpp);
  }
  lbpp_mat
}

#[inline]
pub fn logsumexp2(xs: SliceOfEpsOfTerms4LogPf, max: ExpPartOfTerm4LogPf) -> LogPf {
  if !max.is_finite() {
    xs.iter().fold(0., |acc, &x| acc + x.exp2()).log2()
  } else {
    xs.iter().fold(0., |acc, &x| acc + (x - max).exp2()).log2() + max
  }
}

#[inline]
pub fn get_valid_log_prob(log_prob: LogProb) -> LogProb {
  if log_prob <= 0. {log_prob} else {0.}
}

#[inline]
pub fn get_lbpp_mat_from_parasor(seq: SeqSlice, max_bp_span: usize, min_bpp: Prob) -> SparseProbMat {
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
      get_valid_log_prob(strings[3].parse::<Prob>().expect("Failed to parse a string.").log2())
    )
  }).collect::<SparseLogProbMat>()
}

#[inline]
pub fn run_command(command: &str, args: &[&str], expect: &str) -> Output {
  Command::new(command).args(args).output().expect(expect)
}
