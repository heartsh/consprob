extern crate bio_seq_algos;
extern crate rna_algos;
extern crate itertools;

pub use rna_algos::utils::*;
use itertools::multizip;

pub type PosQuadruple = (Pos, Pos, Pos, Pos);
pub type Prob4dMat = HashMap<PosQuadruple, Prob, Hasher>;
pub type SparseProbMat = HashMap<PosPair, Prob, Hasher>;
pub struct Stapmp {
  pub base_pair_align_prob_mat: Prob4dMat, 
  pub base_align_prob_mat: ProbMat,
}
pub type LogProb4dMat = HashMap<PosQuadruple, LogProb, Hasher>;
pub type SparseLogProbMat = HashMap<PosPair, LogProb, Hasher>;
pub struct LogStapmp {
  pub log_bpap_mat: LogProb4dMat,
  pub log_bap_mat: LogProbMat,
}
type LogPpf4dMat = HashMap<PosQuadruple, LogPf, Hasher>;
pub type ProbMatPair = (ProbMat, ProbMat);
pub type StaScore = LogProb;
pub struct StaScoringParams {
  pub log_null_hypothesis_bpp: LogProb,
  pub log_nh_bap: LogProb,
  pub scaling_param_4_bpa_score: StaScore,
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
type Poss = Vec<Pos>;
type PosSeqsWithPosSpans = HashMap<PosSpan, Poss, Hasher>;
type PosSpanPair = (PosSpan, PosSpan);
type PosPairs = Vec<PosPair>;
type PosPairSeqsWithPosSpanPairs = HashMap<PosSpanPair, PosPairs, Hasher>;
type SparseProbMatPair = (SparseProbMat, SparseProbMat);
type LogProbMatPair = (SparseLogProbMat, SparseLogProbMat);
type PosSeqsWithPoss = HashMap<Pos, Poss, Hasher>;
type PosSeqSetPairWithPoss = (PosSeqsWithPoss, PosSeqsWithPoss);
type PosPairSeqsWithPosPairs = HashMap<PosPair, PosPairs, Hasher>;

impl LogStapmp {
  fn new(seq_len_pair: &(usize, usize)) -> LogStapmp {
    LogStapmp {
      log_bpap_mat: LogProb4dMat::default(),
      log_bap_mat: vec![vec![NEG_INFINITY; seq_len_pair.1 + 2]; seq_len_pair.0 + 2], 
    }
  }
}

impl StaScoringParams {
  pub fn new(lnhbpp: LogProb, lnhbap: LogProb, sp4bpas: StaScore, obpas: StaScore, bogp: StaScore, begp: StaScore, logp: StaScore, legp: StaScore) -> StaScoringParams {
    StaScoringParams {
      log_null_hypothesis_bpp: lnhbpp,
      log_nh_bap: lnhbap,
      scaling_param_4_bpa_score: sp4bpas,
      offset_bpa_score: obpas,
      base_opening_gap_penalty: bogp,
      base_extending_gap_penalty: begp,
      loop_opening_gap_penalty: logp,
      loop_extending_gap_penalty: legp,
    }
  }
}

#[inline]
pub fn io_algo_4_rna_stapmp(seq_len_pair: &(usize, usize), bpp_mat_pair: &ProbMatPair, log_bap_mat: &LogProbMat, sta_scoring_params: &StaScoringParams, min_bpp: Prob) -> Stapmp {
  let mut pos_seq_set_pair_with_pos_spans = (PosSeqsWithPosSpans::default(), PosSeqsWithPosSpans::default());
  let mut sparse_bpp_mat_pair = (SparseProbMat::default(), SparseProbMat::default());
  let mut sparse_log_bpp_mat_pair = (SparseLogProbMat::default(), SparseLogProbMat::default());
  let mut sparse_not_bpp_mat_pair = (SparseProbMat::default(), SparseProbMat::default());
  let mut sparse_log_nbpp_mat_pair = (SparseLogProbMat::default(), SparseLogProbMat::default());
  for substr_len in 2 .. seq_len_pair.0 + 1 {
    let mut poss = Poss::new();
    for i in 0 .. seq_len_pair.0 - substr_len + 1 {
      let j = i + substr_len - 1;
      let bpp = bpp_mat_pair.0[i][j];
      if bpp >= min_bpp && bpp > 0. {
        let lbpp = bpp.log2();
        let nbpp = 1. - bpp;
        let lnbpp = nbpp.log2();
        poss.push(i + 1);
        let pos_pair = (i + 1, j + 1);
        sparse_bpp_mat_pair.0.insert(pos_pair, bpp);
        sparse_log_bpp_mat_pair.0.insert(pos_pair, lbpp);
        sparse_not_bpp_mat_pair.0.insert(pos_pair, nbpp);
        sparse_log_nbpp_mat_pair.0.insert(pos_pair, lnbpp);
      }
    }
    if poss.len() > 0 {
      pos_seq_set_pair_with_pos_spans.0.insert(substr_len, poss);
    }
  }
  for substr_len in 2 .. seq_len_pair.1 + 1 {
    let mut poss = Poss::new();
    for i in 0 .. seq_len_pair.1 - substr_len + 1 {
      let j = i + substr_len - 1;
      let bpp = bpp_mat_pair.1[i][j];
      if bpp >= min_bpp && bpp > 0. {
        let lbpp = bpp.log2();
        let nbpp = 1. - bpp;
        let lnbpp = nbpp.log2();
        poss.push(i + 1);
        let pos_pair = (i + 1, j + 1);
        sparse_bpp_mat_pair.1.insert(pos_pair, bpp);
        sparse_log_bpp_mat_pair.1.insert(pos_pair, lbpp);
        sparse_not_bpp_mat_pair.1.insert(pos_pair, nbpp);
        sparse_log_nbpp_mat_pair.1.insert(pos_pair, lnbpp);
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
                  pos_pairs.push((i, k));
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
  for i in 1 .. seq_len_pair.0 + 1 {
    let mut forward_poss = Poss::new();
    let mut backward_poss = Poss::new();
    let mut num_of_allowed_bps = 0.;
    let max_num_of_allowed_bps = (1. / min_bpp).floor();
    for j in 1 .. seq_len_pair.0 + 1 {
      if i == j {continue;}
      let pos_pair = (min(i, j), max(i, j));
      match sparse_bpp_mat_pair.0.get(&pos_pair) {
        Some(_) => {
          if j < i {
            forward_poss.push(j);
          } else {
            backward_poss.push(j);
          }
          num_of_allowed_bps += 1.;
          if num_of_allowed_bps == max_num_of_allowed_bps {
            break;
          }
        }, None => {},
      }
    }
    if forward_poss.len() > 0 {
      pos_seq_set_pair_with_poss_4_forward_bps.0.insert(i, forward_poss);
    }
    if backward_poss.len() > 0 {
      pos_seq_set_pair_with_poss_4_backward_bps.0.insert(i, backward_poss);
    }
  }
  for i in 1 .. seq_len_pair.1 + 1 {
    let mut forward_poss = Poss::new();
    let mut backward_poss = Poss::new();
    let mut num_of_allowed_bps = 0.;
    let max_num_of_allowed_bps = (1. / min_bpp).floor();
    for j in 1 .. seq_len_pair.1 + 1 {
      if i == j {continue;}
      let pos_pair = (min(i, j), max(i, j));
      match sparse_bpp_mat_pair.1.get(&pos_pair) {
        Some(_) => {
          if j < i {
            forward_poss.push(j);
          } else {
            backward_poss.push(j);
          }
          num_of_allowed_bps += 1.;
          if num_of_allowed_bps == max_num_of_allowed_bps {
            break;
          }
        }, None => {},
      }
    }
    if forward_poss.len() > 0 {
      pos_seq_set_pair_with_poss_4_forward_bps.1.insert(i, forward_poss);
    }
    if backward_poss.len() > 0 {
      pos_seq_set_pair_with_poss_4_backward_bps.1.insert(i, backward_poss);
    }
  }
  let mut pos_pair_seqs_with_pos_pairs_4_forward_bpas = PosPairSeqsWithPosPairs::default();
  let mut pos_pair_seqs_with_pos_pairs_4_backward_bpas = PosPairSeqsWithPosPairs::default();
  for i in 1 .. seq_len_pair.0 + 1 {
    for k in 1 .. seq_len_pair.1 + 1 {
      let mut forward_pos_pairs = PosPairs::new();
      match pos_seq_set_pair_with_poss_4_forward_bps.0.get(&i) {
        Some(poss_1) => {
          match pos_seq_set_pair_with_poss_4_forward_bps.1.get(&k) {
            Some(poss_2) => {
              for &j in poss_1 {
                for &l in poss_2 {
                  forward_pos_pairs.push((j, l));
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
                  backward_pos_pairs.push((j, l));
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
  }
  let log_sta_ppf_mat_4_bpas_1 = get_log_sta_ppf_mat_4_bpas_1(&seq_len_pair, &sparse_log_bpp_mat_pair, log_bap_mat, sta_scoring_params, &sparse_bpp_mat_pair, &sparse_not_bpp_mat_pair, &sparse_log_nbpp_mat_pair, &pos_pair_seqs_with_pos_span_pairs, &pos_seq_set_pair_with_poss_4_forward_bps, &pos_pair_seqs_with_pos_pairs_4_forward_bpas);
  let lstapmp = get_lstapmp(&log_sta_ppf_mat_4_bpas_1, &seq_len_pair, &sparse_log_bpp_mat_pair, &log_bap_mat, sta_scoring_params, &sparse_bpp_mat_pair, &sparse_not_bpp_mat_pair, &sparse_log_nbpp_mat_pair, &pos_pair_seqs_with_pos_span_pairs, &pos_pair_seqs_with_pos_pairs_4_forward_bpas, &pos_pair_seqs_with_pos_pairs_4_backward_bpas);
  get_stapmp(&lstapmp)
}

#[inline]
fn get_stapmp(lstapmp: &LogStapmp) -> Stapmp {
  Stapmp {
    base_pair_align_prob_mat: lstapmp.log_bpap_mat.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect(),
    base_align_prob_mat: lstapmp.log_bap_mat.iter().map(|lbaps| lbaps.iter().map(|&lbap| lbap.exp()).collect()).collect(),
  }
}

#[inline]
pub fn get_log_sta_ppf_mat_4_bpas_1(seq_len_pair: &(usize, usize), sparse_log_bpp_mat_pair: &LogProbMatPair, log_bap_mat: &LogProbMat, sta_scoring_params: &StaScoringParams, sparse_bpp_mat_pair: &SparseProbMatPair, sparse_not_bpp_mat_pair: &SparseProbMatPair, sparse_log_nbpp_mat_pair: &LogProbMatPair, pos_pair_seqs_with_pos_span_pairs: &PosPairSeqsWithPosSpanPairs, pos_seq_set_pair_with_poss_4_forward_bps: &PosSeqSetPairWithPoss, pos_pair_seqs_with_pos_pairs_4_forward_bpas: &PosPairSeqsWithPosPairs, ) -> LogPpf4dMat {
  let mut log_sta_ppf_mat_4_bpas_1 = LogPpf4dMat::default();
  let mut log_sta_ppf_mat_4_bpas_2 = log_sta_ppf_mat_4_bpas_1.clone();
  for substr_len_1 in 2 .. seq_len_pair.0 + 3 {
    for substr_len_2 in 2 .. seq_len_pair.1 + 3 {
      match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let j = i + substr_len_1 - 1;
            let l = k + substr_len_2 - 1;
            let log_sta_ppf_mat_4_2_loops = get_log_sta_ppf_mat_4_2_loops_1(&(i, j, k, l), log_bap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_forward_bpas, &log_sta_ppf_mat_4_bpas_1);
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
                        let ep_of_term_4_log_pf = get_legp(&(i + 1, m - 1), sta_scoring_params) + get_logp(&(m, n), &sparse_log_bpp_mat_pair.0, sta_scoring_params) + log_sta_pf_4_bpa;
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
                        let ep_of_term_4_log_pf = get_legp(&(k + 1, o - 1), sta_scoring_params) + get_logp(&(o, p), &sparse_log_bpp_mat_pair.1, sta_scoring_params) + log_sta_pf_4_bpa;
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
            log_sta_ppf_mat_4_bpas_1.insert((i, j, k, l), get_bpa_score(&(i, j, k, l), seq_len_pair, sparse_log_bpp_mat_pair, log_bap_mat, sta_scoring_params, sparse_bpp_mat_pair, sparse_not_bpp_mat_pair, sparse_log_nbpp_mat_pair) + log_sta_ppf_mat_4_bpas_2[&(i, j, k, l)]);
          }
        }, None => {},
      }
    }
  }
  log_sta_ppf_mat_4_bpas_1
}

#[inline]
fn get_log_sta_ppf_mat_4_2_loops_1(pos_quadruple: &PosQuadruple, log_bap_mat: &LogProbMat, sta_scoring_params: &StaScoringParams, pos_pair_seqs_with_pos_pairs_4_forward_bpas: &PosPairSeqsWithPosPairs, log_sta_ppf_mat_4_bpas_1: &LogPpf4dMat) -> LogPpfMat {
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
      let ba_score = get_ba_score(&(n, p), log_bap_mat, sta_scoring_params);
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
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bas[n - i - 1][p - k] + sta_scoring_params.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
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
fn get_logp(pos_pair: &PosPair, sparse_log_bpp_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams) -> StaScore {
  get_bp_score(pos_pair, sparse_log_bpp_mat, sta_scoring_params) + 2. * sta_scoring_params.loop_opening_gap_penalty
}

#[inline]
fn get_legp(pos_pair: &PosPair, sta_scoring_params: &StaScoringParams) -> StaScore {
  (pos_pair.1 + 1 - pos_pair.0) as StaScore * sta_scoring_params.loop_extending_gap_penalty
}

#[inline]
fn get_ba_score(pos_pair: &PosPair, log_bap_mat: &LogProbMat, sta_scoring_params: &StaScoringParams) -> StaScore {
  log_bap_mat[pos_pair.0 - 1][pos_pair.1 - 1] - sta_scoring_params.log_nh_bap
}

#[inline]
fn get_bp_score(pos_pair: &PosPair, sparse_log_bpp_mat: &SparseLogProbMat, sta_scoring_params: &StaScoringParams) -> StaScore {
  sparse_log_bpp_mat[pos_pair] - sta_scoring_params.log_null_hypothesis_bpp
}

#[inline]
fn get_bpa_score(pos_quadruple: &PosQuadruple, seq_len_pair: &(usize, usize), sparse_log_bpp_mat_pair: &LogProbMatPair, log_bap_mat: &LogProbMat, sta_scoring_params: &StaScoringParams, sparse_bpp_mat_pair: &SparseProbMatPair, sparse_not_bpp_mat_pair: &SparseProbMatPair, sparse_log_nbpp_mat_pair: &LogProbMatPair) -> StaScore {
  if *pos_quadruple == (0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1) {0.} else {
    let bp_pos_pair_1 = (pos_quadruple.0, pos_quadruple.1);
    let bp_pos_pair_2 = (pos_quadruple.2, pos_quadruple.3);
    let aligned_pos_pair_1 = (pos_quadruple.0, pos_quadruple.2);
    let aligned_pos_pair_2 = (pos_quadruple.1, pos_quadruple.3);
    let bin_prob_dist_1 = [sparse_bpp_mat_pair.0[&bp_pos_pair_1], sparse_not_bpp_mat_pair.0[&bp_pos_pair_1]];
    let bin_prob_dist_2 = [sparse_bpp_mat_pair.1[&bp_pos_pair_2], sparse_not_bpp_mat_pair.1[&bp_pos_pair_2]];
    let bin_log_prob_dist_1 = [sparse_log_bpp_mat_pair.0[&bp_pos_pair_1], sparse_log_nbpp_mat_pair.0[&bp_pos_pair_1]];
    let bin_log_prob_dist_2 = [sparse_log_bpp_mat_pair.1[&bp_pos_pair_2], sparse_log_nbpp_mat_pair.1[&bp_pos_pair_2]];
    get_bp_score(&bp_pos_pair_1, &sparse_log_bpp_mat_pair.0, sta_scoring_params)
    + get_bp_score(&bp_pos_pair_2, &sparse_log_bpp_mat_pair.1, sta_scoring_params)
    + get_ba_score(&aligned_pos_pair_1, log_bap_mat, sta_scoring_params)
    + get_ba_score(&aligned_pos_pair_2, log_bap_mat, sta_scoring_params)
    + sta_scoring_params.scaling_param_4_bpa_score * (sta_scoring_params.offset_bpa_score - get_jensen_shannon_dist(&bin_prob_dist_1[..], &bin_prob_dist_2[..], &bin_log_prob_dist_1[..], &bin_log_prob_dist_2[..]))
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
fn get_lstapmp(log_sta_ppf_mat_4_bpas_1: &LogPpf4dMat, seq_len_pair: &(usize, usize), sparse_log_bpp_mat_pair: &LogProbMatPair, log_bap_mat: &LogProbMat, sta_scoring_params: &StaScoringParams, sparse_bpp_mat_pair: &SparseProbMatPair, sparse_not_bpp_mat_pair: &SparseProbMatPair, sparse_log_nbpp_mat_pair: &LogProbMatPair, pos_pair_seqs_with_pos_span_pairs: &PosPairSeqsWithPosSpanPairs, pos_pair_seqs_with_pos_pairs_4_forward_bpas: &PosPairSeqsWithPosPairs, pos_pair_seqs_with_pos_pairs_4_backward_bpas: &PosPairSeqsWithPosPairs) -> LogStapmp {
  let mut lstapmp = LogStapmp::new(seq_len_pair);
  let pos_quadruple_4_2_external_loops = (0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1);
  let log_sta_pf_4_bpa = log_sta_ppf_mat_4_bpas_1[&pos_quadruple_4_2_external_loops];
  let log_sta_ppf_mat_4_2_els_1 = get_log_sta_ppf_mat_4_2_loops_1(&pos_quadruple_4_2_external_loops, log_bap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_forward_bpas, log_sta_ppf_mat_4_bpas_1);
  let log_sta_ppf_mat_4_2_els_2 = get_log_sta_ppf_mat_4_2_loops_2(&pos_quadruple_4_2_external_loops, log_bap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_backward_bpas, log_sta_ppf_mat_4_bpas_1);
  for substr_len_1 in (2 .. seq_len_pair.0 + 1).rev() {
    for substr_len_2 in (2 .. seq_len_pair.1 + 1).rev() {
      match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let j = i + substr_len_1 - 1;
            let l = k + substr_len_2 - 1;
            lstapmp.log_bpap_mat.insert((i, j, k, l), log_sta_ppf_mat_4_2_els_1[i - 1][k - 1] + log_sta_ppf_mat_4_bpas_1[&(i, j, k, l)] + log_sta_ppf_mat_4_2_els_2[j - 1 + 1][l - 1 + 1] - log_sta_pf_4_bpa);
          }
        }, None => {},
      }
    }
  }
  for i in 1 .. seq_len_pair.0 + 1 {
    for k in 1 .. seq_len_pair.1 + 1 {
      let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
      let mut max_ep_of_term_4_log_prob = log_sta_ppf_mat_4_2_els_1[i - 1][k - 1]
      + get_ba_score(&(i, k), log_bap_mat, sta_scoring_params)
      + log_sta_ppf_mat_4_2_els_2[i - 1 + 1][k - 1 + 1];
      if max_ep_of_term_4_log_prob.is_finite() {
        eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
      }
      match pos_pair_seqs_with_pos_pairs_4_forward_bpas.get(&(i, k)) {
        Some(pos_pairs) => {
          for &(j, l) in pos_pairs {
            let ep_of_term_4_log_prob = log_sta_ppf_mat_4_2_els_1[j - 1][l - 1]
            + log_sta_ppf_mat_4_bpas_1[&(j, i, l, k)]
            + log_sta_ppf_mat_4_2_els_2[i - 1 + 1][k - 1 + 1];
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
            + log_sta_ppf_mat_4_2_els_2[j - 1 + 1][l - 1 + 1];
            if ep_of_term_4_log_prob.is_finite() {
              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
              if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
                max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
              }
            }
          }
        }, None => {},
      }
      lstapmp.log_bap_mat[i][k] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob) - log_sta_pf_4_bpa;
    }
  }
  for substr_len_1 in (3 .. seq_len_pair.0 + 1).rev() {
    for substr_len_2 in (3 .. seq_len_pair.1 + 1).rev() {
      match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let j = i + substr_len_1 - 1;
            let l = k + substr_len_2 - 1;
            let lbpap = lstapmp.log_bpap_mat[&(i, j, k, l)];
            let log_sta_pf_4_bpa = log_sta_ppf_mat_4_bpas_1[&(i, j, k, l)];
            let sum_of_terms_4_ep_of_term_4_log_prob = lbpap
            + get_bpa_score(&(i, j, k, l), seq_len_pair, sparse_log_bpp_mat_pair, log_bap_mat, sta_scoring_params, sparse_bpp_mat_pair, sparse_not_bpp_mat_pair, sparse_log_nbpp_mat_pair)
            - log_sta_pf_4_bpa;
            let log_sta_ppf_mat_4_2_loops_1 = get_log_sta_ppf_mat_4_2_loops_1(&(i, j, k, l), log_bap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_forward_bpas, log_sta_ppf_mat_4_bpas_1);
            let log_sta_ppf_mat_4_2_loops_2 = get_log_sta_ppf_mat_4_2_loops_2(&(i, j, k, l), log_bap_mat, sta_scoring_params, pos_pair_seqs_with_pos_pairs_4_backward_bpas, log_sta_ppf_mat_4_bpas_1);
            for substr_len_3 in (2 .. substr_len_1 - 1).rev() {
              for substr_len_4 in (2 .. substr_len_2 - 1).rev() {
                match pos_pair_seqs_with_pos_span_pairs.get(&(substr_len_3, substr_len_4)) {
                  Some(pos_pairs) => {
                    for &(m, o) in pos_pairs {
                      let n = m + substr_len_3 - 1;
                      let p = o + substr_len_4 - 1;
                      if m <= i || n >= j || o <= k || p >= l {continue;}
                      let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
                      let mut max_ep_of_term_4_log_prob = lstapmp.log_bpap_mat[&(m, n, o, p)];
                      if max_ep_of_term_4_log_prob.is_finite() {
                        eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
                      }
                      let ep_of_term_4_log_prob = sum_of_terms_4_ep_of_term_4_log_prob
                      + log_sta_ppf_mat_4_2_loops_1[m - i - 1][o - k - 1]
                      + log_sta_ppf_mat_4_bpas_1[&(m, n, o, p)]
                      + log_sta_ppf_mat_4_2_loops_2[n - (i + 1) + 1][p - (k + 1) + 1];
                      if ep_of_term_4_log_prob.is_finite() {
                        eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                        if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
                          max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                        if eps_of_terms_4_log_prob.len() > 0 {
                          *lstapmp.log_bpap_mat.get_mut(&(m, n, o, p)).expect("Failed to index a hash map.") = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
                        }
                      }
                    }
                  }, None => {},
                }
              }
            }
            for m in i + 1 .. j {
              for o in k + 1 .. l {
                let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
                let mut max_ep_of_term_4_log_prob = lstapmp.log_bap_mat[m][o];
                if max_ep_of_term_4_log_prob.is_finite() {
                  eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
                }
                let ep_of_term_4_log_prob = sum_of_terms_4_ep_of_term_4_log_prob
                + log_sta_ppf_mat_4_2_loops_1[m - i - 1][o - k - 1]
                + get_ba_score(&(m, o), log_bap_mat, sta_scoring_params)
                + log_sta_ppf_mat_4_2_loops_2[m - (i + 1) + 1][o - (k + 1) + 1];
                if ep_of_term_4_log_prob.is_finite() {
                  eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                  if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
                    max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
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
                          if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
                            max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
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
                          if ep_of_term_4_log_prob > max_ep_of_term_4_log_prob {
                            max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                          }
                        }
                      }
                    }, None => {},
                  }
                }
                if eps_of_terms_4_log_prob.len() > 0 {
                  lstapmp.log_bap_mat[m][o] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
                }
              }
            }
          }
        }, None => {},
      }
    }
  }
  lstapmp.log_bpap_mat.insert((0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1), 0.);
  lstapmp.log_bap_mat[0][0] = 0.;
  lstapmp.log_bap_mat[seq_len_pair.0 + 1][seq_len_pair.1 + 1] = 0.;
  lstapmp
}

#[inline]
fn get_log_sta_ppf_mat_4_2_loops_2(pos_quadruple: &PosQuadruple, log_bap_mat: &LogProbMat, sta_scoring_params: &StaScoringParams, pos_pair_seqs_with_pos_pairs_4_backward_bpas: &PosPairSeqsWithPosPairs, log_sta_ppf_mat_4_bpas_1: &LogPpf4dMat) -> LogPpfMat {
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
      let ba_score = get_ba_score(&(m, o), log_bap_mat, sta_scoring_params);
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
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_4_2_loops_and_bas[m - (i + 1)][o - (k + 1)] + sta_scoring_params.base_opening_gap_penalty;
      if max_ep_of_term_4_log_pf.is_finite() {
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
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
