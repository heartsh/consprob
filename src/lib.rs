extern crate bio_seq_algos;
extern crate rna_algos;
extern crate itertools;

use bio_seq_algos::utils::fast_ln;
pub use rna_algos::utils::*;
use itertools::multizip;

pub type ThreeDPm = Vec<ProbMatrix>;
pub type FourDPm = Vec<ThreeDPm>;
pub struct Stapmp {
  pub bpap_matrix: FourDPm, 
  pub bap_matrix: ProbMatrix,
}
pub type Log3DPm = Vec<LogProbMatrix>;
pub type Log4DPm = Vec<Log3DPm>;
pub struct LogStapmp {
  pub log_bpap_matrix: Log4DPm, 
  pub log_bap_matrix: LogProbMatrix,
}
type LogPpf3DMatrix = Vec<LogPpfMatrix>;
type LogPpf4DMatrix = Vec<LogPpf3DMatrix>;
pub struct LogStaPpfMatrices {
  pub log_sta_ppf_matrix: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_rightmost_bpas: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_bpas: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_bas: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_base_gaps_1: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_bgs_2: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_2_loop_deletions_1: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_2lds_2: LogPpf4DMatrix,
}
pub type PmPair = (ProbMatrix, ProbMatrix);
pub type LpmPair = (LogProbMatrix, LogProbMatrix);
pub type StaScore = LogProb;
pub struct StaScoringParams {
  pub log_null_hypothesis_bpp: StaScore,
  pub log_nh_bap: StaScore,
  pub scaling_parameter_4_ss_distance: StaScore,
  pub base_opening_gap_penalty: StaScore,
  pub base_extending_gap_penalty: StaScore,
  pub loop_opening_gap_penalty: StaScore,
  pub loop_extending_gap_penalty: StaScore,
}
type ProbDist = Probs;
type LogProbDist = LogProbs;
type PdSlice<'a> = &'a[Prob];
type LpdSlice<'a> = &'a[LogProb];
type KullbackLeiblerDiv = LogProb;
type JensenShannonDist = KullbackLeiblerDiv;
type Index = usize;
type IndexPair = (Index, Index);
type Ips = Vec<IndexPair>;
type Ssl = usize;
type IpsWithSsls = HashMap<Ssl, Ips, Hasher>;
type PairOfIpsWithSsls = (IpsWithSsls, IpsWithSsls);
type ProbsWithBppt4StaSs = HashMap<Index, Prob, Hasher>;
type ProbMatWithBppt4StaSs = HashMap<Index, ProbsWithBppt4StaSs, Hasher>;
type PmPairWithBppt4StaSs = (ProbMatWithBppt4StaSs, ProbMatWithBppt4StaSs);
type LpsWithBppt4StaSs = HashMap<Index, LogProb, Hasher>;
type LpMatWithBppt4StaSs = HashMap<Index, LpsWithBppt4StaSs, Hasher>;
type LpmPairWithBppt4StaSs = (LpMatWithBppt4StaSs, LpMatWithBppt4StaSs);
type IndicesWithBppt4StaSs = HashMap<Index, bool, Hasher>;
type IndexMatWithBppt4StaSs = HashMap<Index, IndicesWithBppt4StaSs, Hasher>;
type ImPairWithBppt4StaSs = (IndexMatWithBppt4StaSs, IndexMatWithBppt4StaSs);

impl LogStaPpfMatrices {
  fn new(slp: &(usize, usize)) -> LogStaPpfMatrices {
    let ni_matrix = vec![vec![vec![vec![NEG_INFINITY; slp.1 + 2]; slp.1 + 2]; slp.0 + 2]; slp.0 + 2];
    LogStaPpfMatrices {
      log_sta_ppf_matrix: ni_matrix.clone(),
      log_sta_ppf_matrix_4_rightmost_bpas: ni_matrix.clone(),
      log_sta_ppf_matrix_4_bas: ni_matrix.clone(),
      log_sta_ppf_matrix_4_base_gaps_1: ni_matrix.clone(),
      log_sta_ppf_matrix_4_bgs_2: ni_matrix.clone(),
      log_sta_ppf_matrix_4_bpas: ni_matrix.clone(),
      log_sta_ppf_matrix_4_2_loop_deletions_1: ni_matrix.clone(),
      log_sta_ppf_matrix_4_2lds_2: ni_matrix,
    }
  }
}

impl LogStapmp {
  fn new(slp: &(usize, usize)) -> LogStapmp {
    LogStapmp {
      log_bpap_matrix: vec![vec![vec![vec![NEG_INFINITY; slp.1 + 2]; slp.1 + 2]; slp.0 + 2]; slp.0 + 2],
      log_bap_matrix: vec![vec![NEG_INFINITY; slp.1 + 2]; slp.0 + 2]
    }
  }
}

impl StaScoringParams {
  pub fn new(log_nh_bpp: LogProb, log_nh_bap: LogProb, sp_4_ssd: StaScore, bogp: StaScore, begp: StaScore, logp: StaScore, legp: StaScore) -> StaScoringParams {
    StaScoringParams {
      log_null_hypothesis_bpp: log_nh_bpp,
      log_nh_bap: log_nh_bap,
      scaling_parameter_4_ss_distance: sp_4_ssd,
      base_opening_gap_penalty: bogp,
      base_extending_gap_penalty: begp,
      loop_opening_gap_penalty: logp,
      loop_extending_gap_penalty: legp,
    }
  }
}

#[inline]
pub fn io_algo_4_rna_stapmp(seq_len_pair: &(usize, usize), bpp_mp: &PmPair, log_bap_matrix: &LogProbMatrix, sta_sps: &StaScoringParams, bpp_thres_4_sta_space_sparsification: Prob) -> Stapmp {
  let mut poips_with_ssls_with_bppt_4_sta_ss = (IpsWithSsls::default(), IpsWithSsls::default());
  for ssl in 2 .. seq_len_pair.0 + 1 {
    let mut ips = Ips::new();
    for i in 0 .. seq_len_pair.0 - ssl + 1 {
      let j = i + ssl - 1;
      let bpp = bpp_mp.0[i][j];
      if bpp >= bpp_thres_4_sta_space_sparsification {
        ips.push((i + 1, j + 1));
      }
    }
    if ips.len() > 0 {
      poips_with_ssls_with_bppt_4_sta_ss.0.insert(ssl, ips);
    }
  }
  for ssl in 2 .. seq_len_pair.1 + 1 {
    let mut ips = Ips::new();
    for i in 0 .. seq_len_pair.1 - ssl + 1 {
      let j = i + ssl - 1;
      let bpp = bpp_mp.1[i][j];
      if bpp >= bpp_thres_4_sta_space_sparsification {
        ips.push((i + 1, j + 1));
      }
    }
    if ips.len() > 0 {
      poips_with_ssls_with_bppt_4_sta_ss.1.insert(ssl, ips);
    }
  }
  let mut imp_with_bppt_4_sta_ss_1 = (IndexMatWithBppt4StaSs::default(), IndexMatWithBppt4StaSs::default());
  let mut imp_with_bppt_4_sta_ss_2 = (IndexMatWithBppt4StaSs::default(), IndexMatWithBppt4StaSs::default());
  let mut bppmp_with_bppt_4_sta_ss = (ProbMatWithBppt4StaSs::default(), ProbMatWithBppt4StaSs::default());
  let mut nbppmp_with_bppt_4_sta_ss = (ProbMatWithBppt4StaSs::default(), ProbMatWithBppt4StaSs::default());
  let mut log_bppmp_with_bppt_4_sta_ss = (LpMatWithBppt4StaSs::default(), LpMatWithBppt4StaSs::default());
  let mut lnbppmp_with_bppt_4_sta_ss = (LpMatWithBppt4StaSs::default(), LpMatWithBppt4StaSs::default());
  for i in 0 .. seq_len_pair.0 {
    let mut indices_with_bppt_4_sta_ss = IndicesWithBppt4StaSs::default();
    if i > 0 {
      for j in 0 .. i {
        let bpp = bpp_mp.0[j][i];
        if bpp >= bpp_thres_4_sta_space_sparsification {
          indices_with_bppt_4_sta_ss.insert(j + 1, true);
        }
      }
      if indices_with_bppt_4_sta_ss.len() > 0 {
        imp_with_bppt_4_sta_ss_1.0.insert(i + 1, indices_with_bppt_4_sta_ss.clone());
      }
    }
    let mut bpps_with_bppt_4_sta_ss = ProbsWithBppt4StaSs::default();
    let mut nbpps_with_bppt_4_sta_ss = ProbsWithBppt4StaSs::default();
    let mut lbpps_with_bppt_4_sta_ss = LpsWithBppt4StaSs::default();
    let mut lnbpps_with_bppt_4_sta_ss = LpsWithBppt4StaSs::default();
    for j in i + 1 .. seq_len_pair.0 {
      let bpp = bpp_mp.0[i][j];
      let nbpp = 1. - bpp;
      let lbpp = fast_ln(bpp);
      let lnbpp = fast_ln(nbpp);
      if bpp >= bpp_thres_4_sta_space_sparsification {
        bpps_with_bppt_4_sta_ss.insert(j + 1, bpp);
        nbpps_with_bppt_4_sta_ss.insert(j + 1, nbpp);
        lbpps_with_bppt_4_sta_ss.insert(j + 1, lbpp);
        lnbpps_with_bppt_4_sta_ss.insert(j + 1, lnbpp);
        indices_with_bppt_4_sta_ss.insert(j + 1, true);
      }
    }
    if bpps_with_bppt_4_sta_ss.len() > 0 {
      bppmp_with_bppt_4_sta_ss.0.insert(i + 1, bpps_with_bppt_4_sta_ss);
      nbppmp_with_bppt_4_sta_ss.0.insert(i + 1, nbpps_with_bppt_4_sta_ss);
      log_bppmp_with_bppt_4_sta_ss.0.insert(i + 1, lbpps_with_bppt_4_sta_ss);
      lnbppmp_with_bppt_4_sta_ss.0.insert(i + 1, lnbpps_with_bppt_4_sta_ss);
    }
    if indices_with_bppt_4_sta_ss.len() > 0 {
      imp_with_bppt_4_sta_ss_2.0.insert(i + 1, indices_with_bppt_4_sta_ss);
    }
  }
  for i in 0 .. seq_len_pair.1 {
    let mut indices_with_bppt_4_sta_ss = IndicesWithBppt4StaSs::default();
    if i > 0 {
      for j in 0 .. i {
        let bpp = bpp_mp.1[j][i];
        if bpp >= bpp_thres_4_sta_space_sparsification {
          indices_with_bppt_4_sta_ss.insert(j + 1, true);
        }
      }
      if indices_with_bppt_4_sta_ss.len() > 0 {
        imp_with_bppt_4_sta_ss_1.1.insert(i + 1, indices_with_bppt_4_sta_ss.clone());
      }
    }
    let mut bpps_with_bppt_4_sta_ss = ProbsWithBppt4StaSs::default();
    let mut nbpps_with_bppt_4_sta_ss = ProbsWithBppt4StaSs::default();
    let mut lbpps_with_bppt_4_sta_ss = LpsWithBppt4StaSs::default();
    let mut lnbpps_with_bppt_4_sta_ss = LpsWithBppt4StaSs::default();
    for j in i + 1 .. seq_len_pair.1 {
      let bpp = bpp_mp.1[i][j];
      let nbpp = 1. - bpp;
      let lbpp = fast_ln(bpp);
      let lnbpp = fast_ln(nbpp);
      if bpp >= bpp_thres_4_sta_space_sparsification {
        bpps_with_bppt_4_sta_ss.insert(j + 1, bpp);
        nbpps_with_bppt_4_sta_ss.insert(j + 1, nbpp);
        lbpps_with_bppt_4_sta_ss.insert(j + 1, lbpp);
        lnbpps_with_bppt_4_sta_ss.insert(j + 1, lnbpp);
        indices_with_bppt_4_sta_ss.insert(j + 1, true);
      }
    }
    if bpps_with_bppt_4_sta_ss.len() > 0 {
      bppmp_with_bppt_4_sta_ss.1.insert(i + 1, bpps_with_bppt_4_sta_ss);
      nbppmp_with_bppt_4_sta_ss.1.insert(i + 1, nbpps_with_bppt_4_sta_ss);
      log_bppmp_with_bppt_4_sta_ss.1.insert(i + 1, lbpps_with_bppt_4_sta_ss);
      lnbppmp_with_bppt_4_sta_ss.1.insert(i + 1, lnbpps_with_bppt_4_sta_ss);
    }
    if indices_with_bppt_4_sta_ss.len() > 0 {
      imp_with_bppt_4_sta_ss_2.1.insert(i + 1, indices_with_bppt_4_sta_ss);
    }
  }
  let log_sta_ppf_matrices = get_log_sta_ppf_matrices(&seq_len_pair, log_bap_matrix, sta_sps, &bppmp_with_bppt_4_sta_ss, &nbppmp_with_bppt_4_sta_ss, &log_bppmp_with_bppt_4_sta_ss, &lnbppmp_with_bppt_4_sta_ss, &imp_with_bppt_4_sta_ss_1);
  let log_stapmp = get_log_stapmp(&log_sta_ppf_matrices, &seq_len_pair, log_bap_matrix, sta_sps, &poips_with_ssls_with_bppt_4_sta_ss, &bppmp_with_bppt_4_sta_ss, &nbppmp_with_bppt_4_sta_ss, &log_bppmp_with_bppt_4_sta_ss, &lnbppmp_with_bppt_4_sta_ss, &imp_with_bppt_4_sta_ss_1, &imp_with_bppt_4_sta_ss_2);
  get_stapmp(&log_stapmp)
}

#[inline]
fn get_stapmp(log_stapmp: &LogStapmp) -> Stapmp {
  Stapmp {
    bpap_matrix: log_stapmp.log_bpap_matrix.iter().map(|three_d_matrix| three_d_matrix.iter().map(|two_d_matrix| two_d_matrix.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect()).collect()).collect(),
    bap_matrix: log_stapmp.log_bap_matrix.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect(),
  }
}

#[inline]
pub fn get_log_sta_ppf_matrices(slp: &(usize, usize), log_bap_matrix: &LogProbMatrix, sta_sps: &StaScoringParams, bppmp_with_bppt_4_sta_ss: &PmPairWithBppt4StaSs, nbppmp_with_bppt_4_sta_ss: &PmPairWithBppt4StaSs, lbppmp_with_bppt_4_sta_ss: &LpmPairWithBppt4StaSs, lnbppmp_with_bppt_4_sta_ss: &LpmPairWithBppt4StaSs, imp_with_bppt_4_sta_ss_1: &ImPairWithBppt4StaSs) -> LogStaPpfMatrices {
  let mut log_sta_ppf_matrices = LogStaPpfMatrices::new(slp);
  for i in 1 .. slp.0 + 1 {
    for j in i .. slp.0 + 1 {
      for k in 0 .. slp.1 + 1 {
        let bgp = sta_sps.base_opening_gap_penalty + get_begp(&(i + 1, j), sta_sps);
        log_sta_ppf_matrices.log_sta_ppf_matrix_4_base_gaps_1[i][j][k + 1][k] = bgp;
        log_sta_ppf_matrices.log_sta_ppf_matrix[i][j][k + 1][k] = bgp;
      }
    }
  }
  for k in 1 .. slp.1 + 1 {
    for l in k .. slp.1 + 1 {
      for i in 0 .. slp.0 + 1 {
        let bgp = sta_sps.base_opening_gap_penalty + get_begp(&(k + 1, l), sta_sps);
        log_sta_ppf_matrices.log_sta_ppf_matrix_4_bgs_2[i + 1][i][k][l] = bgp;
        log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][i][k][l] = bgp;
      }
    }
  }
  for i in 0 .. slp.0 + 1 {
    for j in 0 .. slp.1 + 1 {
      log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][i][j + 1][j] = 0.;
    }
  }
  for sub_seq_len_1 in 1 .. slp.0 + 1 {
    for i in 1 .. slp.0 - sub_seq_len_1 + 2 {
      let j = i + sub_seq_len_1 - 1;
      let pp_closing_loop_1 = (i, j);
      for ssl_2 in 1 .. slp.1 + 1 {
        for k in 1 .. slp.1 - ssl_2 + 2 {
          let l = k + ssl_2 - 1;
          let pp_closing_loop_2 = (k, l);
          match bppmp_with_bppt_4_sta_ss.0.get(&pp_closing_loop_1.0) {
            Some(bpps_with_bppt_4_sta_ss) => {
              match bpps_with_bppt_4_sta_ss.get(&pp_closing_loop_1.1) {
                Some(_) => {
                  match bppmp_with_bppt_4_sta_ss.1.get(&pp_closing_loop_2.0) {
                    Some(bpps_with_bppt_4_sta_ss) => {
                      match bpps_with_bppt_4_sta_ss.get(&pp_closing_loop_2.1) {
                        Some(_) => {
                          let bpa_score = get_bpa_score(&pp_closing_loop_1, &pp_closing_loop_2, log_bap_matrix, sta_sps, bppmp_with_bppt_4_sta_ss, nbppmp_with_bppt_4_sta_ss, lbppmp_with_bppt_4_sta_ss, lnbppmp_with_bppt_4_sta_ss);
                          if bpa_score.is_finite() {
                            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
                            let mut max_ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix[pp_closing_loop_1.0 + 1][pp_closing_loop_1.1 - 1][pp_closing_loop_2.0 + 1][pp_closing_loop_2.1 - 1];
                            if max_ep_of_term_4_log_pf.is_finite() {
                              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
                            }
                            for h in k + 1 .. l - 1 {
                              let ep_of_term_4_log_pf = get_legp(&(pp_closing_loop_2.0 + 1, h - 1), sta_sps) + log_sta_ppf_matrices.log_sta_ppf_matrix_4_2_loop_deletions_1[pp_closing_loop_1.0 + 1][pp_closing_loop_1.1 - 1][h][pp_closing_loop_2.1 - 1];
                              if ep_of_term_4_log_pf.is_finite() {
                                if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                              }
                            }
                            for h in i + 1 .. j - 1 {
                              let ep_of_term_4_log_pf = get_legp(&(pp_closing_loop_1.0 + 1, h - 1), sta_sps) + log_sta_ppf_matrices.log_sta_ppf_matrix_4_2lds_2[h][pp_closing_loop_1.1 - 1][pp_closing_loop_2.0 + 1][pp_closing_loop_2.1 - 1];
                              if ep_of_term_4_log_pf.is_finite() {
                                if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                              }
                            }
                            if eps_of_terms_4_log_pf.len() > 0 {
                              log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[pp_closing_loop_1.0][pp_closing_loop_1.1][pp_closing_loop_2.0][pp_closing_loop_2.1] = bpa_score + logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
                            }
                          }
                        }, None => {},
                      }
                    }, None => {},
                  }
                }, None => {},
              }
            }, None => {},
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          match imp_with_bppt_4_sta_ss_1.0.get(&j) {
            Some(indices_with_bppt_4_sta_ss) => {
              for &h_1 in indices_with_bppt_4_sta_ss.keys() {
                if h_1 >= i {
                  match imp_with_bppt_4_sta_ss_1.1.get(&l) {
                    Some(indices_with_bppt_4_sta_ss) => {
                      for &h_2 in indices_with_bppt_4_sta_ss.keys() {
                        if h_2 >= k {
                          let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[h_1][j][h_2][l];
                          if log_sta_ppf_4_bpa.is_finite() {
                            let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix[i][h_1 - 1][k][h_2 - 1] + log_sta_ppf_4_bpa;
                            if ep_of_term_4_log_pf.is_finite() {
                              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                            }
                          }
                        }
                      }
                    }, None => {},
                  }
                }
              }
            }, None => {},
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_rightmost_bpas[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          log_sta_ppf_matrices.log_sta_ppf_matrix_4_bas[i][j][k][l] = if i == j && k == l {0.} else {log_sta_ppf_matrices.log_sta_ppf_matrix[i][j - 1][k][l - 1]} + get_ba_score(&(j, l), log_bap_matrix, sta_sps);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_rightmost_bpas[i][j - 1][k][l] + sta_sps.base_opening_gap_penalty;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bas[i][j - 1][k][l] + sta_sps.base_opening_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
          let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_base_gaps_1[i][j - 1][k][l] + sta_sps.base_extending_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_base_gaps_1[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_rightmost_bpas[i][j][k][l - 1] + sta_sps.base_opening_gap_penalty;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bas[i][j][k][l - 1] + sta_sps.base_opening_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
          let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bgs_2[i][j][k][l - 1] + sta_sps.base_extending_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_bgs_2[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          let eps_of_terms_4_log_pf = [
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_rightmost_bpas[i][j][k][l],
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_bas[i][j][k][l],
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_base_gaps_1[i][j][k][l],
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_bgs_2[i][j][k][l],
          ];
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_ppf_matrices.log_sta_ppf_matrix[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          match bppmp_with_bppt_4_sta_ss.1.get(&k) {
            Some(bpps_with_bppt_4_sta_ss) => {
              for (&h, &bpp) in bpps_with_bppt_4_sta_ss {
                if h <= l {
                  let pp_cl = (k, h);
                  if bpp > 0. {
                    let ep_of_term_4_log_pf = get_logp(&pp_cl, sta_sps, &lbppmp_with_bppt_4_sta_ss.1) + log_sta_ppf_matrices.log_sta_ppf_matrix[i][j][k + 1][h - 1] + get_legp(&(h + 1, l), sta_sps);
                    if ep_of_term_4_log_pf.is_finite() {
                      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    }
                  }
                }
              }
            }, None => {},
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_2_loop_deletions_1[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          match bppmp_with_bppt_4_sta_ss.0.get(&i) {
            Some(bpps_with_bppt_4_sta_ss) => {
              for (&h, &bpp) in bpps_with_bppt_4_sta_ss {
                if h <= j {
                  let pp_cl = (i, h);
                  if bpp > 0. {
                    let ep_of_term_4_log_pf = get_logp(&pp_cl, sta_sps, &lbppmp_with_bppt_4_sta_ss.0) + log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][h - 1][k][l] + get_legp(&(h + 1, j), sta_sps);
                    if ep_of_term_4_log_pf.is_finite() {
                      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    }
                  }
                }
              }
            }, None => {},
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_2lds_2[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
        }
      }
    }
  }
  log_sta_ppf_matrices
}

#[inline]
fn get_begp(pp: &PosPair, sta_sps: &StaScoringParams) -> StaScore {
  (pp.1 + 1 - pp.0) as StaScore * sta_sps.base_extending_gap_penalty
}

#[inline]
fn get_logp(pp: &PosPair, sta_sps: &StaScoringParams, lbppm_with_bppt_4_sta_ss: &LpMatWithBppt4StaSs) -> StaScore {
  lbppm_with_bppt_4_sta_ss[&pp.0][&pp.1] - sta_sps.log_null_hypothesis_bpp + 2. * sta_sps.loop_opening_gap_penalty
}

#[inline]
fn get_legp(pp: &PosPair, sta_sps: &StaScoringParams) -> StaScore {
  (pp.1 + 1 - pp.0) as StaScore * sta_sps.loop_extending_gap_penalty
}

#[inline]
fn get_ba_score(pp: &PosPair, log_bap_matrix: &LogProbMatrix, sta_sps: &StaScoringParams) -> StaScore {
  log_bap_matrix[pp.0 - 1][pp.1 - 1] - sta_sps.log_nh_bap
}

#[inline]
fn get_bpa_score(pp_1: &PosPair, pp_2: &PosPair, log_bap_matrix: &LogProbMatrix, sta_sps: &StaScoringParams, bppmp_with_bppt_4_sta_ss: &PmPairWithBppt4StaSs, nbppmp_with_bppt_4_sta_ss: &PmPairWithBppt4StaSs, lbppmp_with_bppt_4_sta_ss: &LpmPairWithBppt4StaSs, lnbppmp_with_bppt_4_sta_ss: &LpmPairWithBppt4StaSs) -> StaScore {
  let bin_pd_1 = [bppmp_with_bppt_4_sta_ss.0[&pp_1.0][&pp_1.1], nbppmp_with_bppt_4_sta_ss.0[&pp_1.0][&pp_1.1]];
  let bin_pd_2 = [bppmp_with_bppt_4_sta_ss.1[&pp_2.0][&pp_2.1], nbppmp_with_bppt_4_sta_ss.1[&pp_2.0][&pp_2.1]];
  let bin_lpd_1 = [lbppmp_with_bppt_4_sta_ss.0[&pp_1.0][&pp_1.1], lnbppmp_with_bppt_4_sta_ss.0[&pp_1.0][&pp_1.1]];
  let bin_lpd_2 = [lbppmp_with_bppt_4_sta_ss.1[&pp_2.0][&pp_2.1], lnbppmp_with_bppt_4_sta_ss.1[&pp_2.0][&pp_2.1]];
  lbppmp_with_bppt_4_sta_ss.0[&pp_1.0][&pp_1.1] + lbppmp_with_bppt_4_sta_ss.1[&pp_2.0][&pp_2.1] - 2. * sta_sps.log_null_hypothesis_bpp
  + get_ba_score(&(pp_1.0, pp_2.0), log_bap_matrix, sta_sps) + get_ba_score(&(pp_1.1, pp_2.1), log_bap_matrix, sta_sps)
  - sta_sps.scaling_parameter_4_ss_distance * get_jsd(&bin_pd_1[..], &bin_pd_2[..], &bin_lpd_1[..], &bin_lpd_2[..])
}

#[inline]
fn get_jsd(pd_1: PdSlice, pd_2: PdSlice, lpd_1: LpdSlice, lpd_2: LpdSlice) -> JensenShannonDist {
  let mid_pd = pd_1.iter().zip(pd_2).map(|(&x, &y)| (x + y) / 2.).collect::<ProbDist>();
  let log_mid_prob_dist = mid_pd.iter().map(|&x| fast_ln(x)).collect::<LogProbDist>();
  let jsd = ((get_kld(pd_1, lpd_1, &log_mid_prob_dist)
  + get_kld(pd_2, lpd_2, &log_mid_prob_dist)) / 2.).sqrt();
  if jsd.is_finite() {jsd} else {0.}
}

#[inline]
fn get_kld(pd: PdSlice, lpd_1: LpdSlice, lpd_2: LpdSlice) -> JensenShannonDist {
  multizip((pd.iter(), lpd_1.iter(), lpd_2.iter())).map(|(&x, &y, &z)| x * (y - z)).sum()
}

#[inline]
fn get_log_stapmp(log_sta_ppf_matrices: &LogStaPpfMatrices, slp: &(usize, usize), log_bap_matrix: &LogProbMatrix, sta_sps: &StaScoringParams, poips_with_ssls_with_bppt_4_sta_ss: &PairOfIpsWithSsls, bppmp_with_bppt_4_sta_ss: &PmPairWithBppt4StaSs, nbppmp_with_bppt_4_sta_ss: &PmPairWithBppt4StaSs, lbppmp_with_bppt_4_sta_ss: &LpmPairWithBppt4StaSs, lnbppmp_with_bppt_4_sta_ss: &LpmPairWithBppt4StaSs, imp_with_bppt_4_sta_ss_1: &ImPairWithBppt4StaSs, imp_with_bppt_4_sta_ss_2: &ImPairWithBppt4StaSs) -> LogStapmp {
  let mut log_stapmp = LogStapmp::new(slp);
  let log_sta_ppf = log_sta_ppf_matrices.log_sta_ppf_matrix[1][slp.0][1][slp.1];
  let mut log_pseudo_prob_matrix = log_stapmp.log_bpap_matrix.clone();
  for ssl_1 in (2 .. slp.0 + 1).rev() {
    match poips_with_ssls_with_bppt_4_sta_ss.0.get(&ssl_1) {
      Some(ips) => {
        for &(i, j) in ips {
          let accessible_pp_1 = (i, j);
          for ssl_2 in (2 .. slp.1 + 1).rev() {
            match poips_with_ssls_with_bppt_4_sta_ss.1.get(&ssl_2) {
              Some(ips) => {
                for &(k, l) in ips {
                  let app_2 = (k, l);
                  let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[accessible_pp_1.0][accessible_pp_1.1][app_2.0][app_2.1];
                  if log_sta_ppf_4_bpa.is_finite() {
                    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogPf::new();
                    let mut max_ep_of_term_4_log_prob = log_sta_ppf_matrices.log_sta_ppf_matrix[1][accessible_pp_1.0 - 1][1][app_2.0 - 1]
                    + log_sta_ppf_matrices.log_sta_ppf_matrix[accessible_pp_1.1 + 1][slp.0][app_2.1 + 1][slp.1]
                    - log_sta_ppf;
                    if max_ep_of_term_4_log_prob.is_finite() {
                      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
                    }
                    for h_1 in 1 .. accessible_pp_1.0 {
                      for h_2 in 1 .. app_2.0 {
                        let ep_of_term_4_log_prob = log_pseudo_prob_matrix[h_1][accessible_pp_1.1][h_2][app_2.1] + log_sta_ppf_matrices.log_sta_ppf_matrix[h_1 + 1][accessible_pp_1.0 - 1][h_2 + 1][app_2.0 - 1];
                        if ep_of_term_4_log_prob.is_finite() {
                          if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                          eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                        }
                      }
                    }
                    if eps_of_terms_4_log_prob.len() > 0 {
                      log_stapmp.log_bpap_matrix[accessible_pp_1.0][accessible_pp_1.1][app_2.0][app_2.1] = log_sta_ppf_4_bpa + logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
                    }
                  }
                  let mut eps_of_terms_4_log_prob = EpsOfTerms4LogPf::new();
                  let mut max_ep_of_term_4_log_prob = NEG_INFINITY;
                  match bppmp_with_bppt_4_sta_ss.0.get(&i) {
                    Some(bpps_with_bppt_4_sta_ss) => {
                      for &h_1 in bpps_with_bppt_4_sta_ss.keys() {
                        if h_1 > j {
                          let app_1 = (i, h_1);
                          match bppmp_with_bppt_4_sta_ss.1.get(&k) {
                            Some(bpps_with_bppt_4_sta_ss) => {
                              for &h_2 in bpps_with_bppt_4_sta_ss.keys() {
                                if h_2 > l {
                                  let app_2 = (k, h_2);
                                  let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
                                  if log_sta_ppf_4_bpa.is_finite() {
                                    let log_bpap = log_stapmp.log_bpap_matrix[app_1.0][app_1.1][app_2.0][app_2.1];
                                    let bpa_score = get_bpa_score(&app_1, &app_2, log_bap_matrix, sta_sps, bppmp_with_bppt_4_sta_ss, nbppmp_with_bppt_4_sta_ss, lbppmp_with_bppt_4_sta_ss, lnbppmp_with_bppt_4_sta_ss);
                                    let ep_of_term_4_log_prob = log_bpap + bpa_score
                                    + log_sta_ppf_matrices.log_sta_ppf_matrix[j + 1][app_1.1 - 1][l + 1][app_2.1 - 1]
                                    - log_sta_ppf_4_bpa;
                                    if ep_of_term_4_log_prob.is_finite() {
                                      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                                      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                                    }
                                  }
                                }
                              }
                            }, None => {},
                          }
                        }
                      }
                    }, None => {},
                  }
                  if eps_of_terms_4_log_prob.len() > 0 {
                    log_pseudo_prob_matrix[i][j][k][l] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
                  }
                }
              }, None => {},
            }
          }
        }
      }, None => {}, 
    }
  }
  for i in 1 .. slp.0 + 1 {
    for j in 1 .. slp.1 + 1 {
      let mut eps_of_terms_4_log_prob = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_prob = log_sta_ppf_matrices.log_sta_ppf_matrix[1][i - 1][1][j - 1]
      + get_ba_score(&(i, j), log_bap_matrix, sta_sps)
      + log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][slp.0][j + 1][slp.1]
      - log_sta_ppf;
      if max_ep_of_term_4_log_prob.is_finite() {
        eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
      }
      let log_sta_ppf_2 = log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][slp.0][j + 1][slp.1];
      match imp_with_bppt_4_sta_ss_1.0.get(&i) {
        Some(indices_with_bppt_4_sta_ss) => {
          for &k in indices_with_bppt_4_sta_ss.keys() {
            let app_1 = (k, i);
            match imp_with_bppt_4_sta_ss_1.1.get(&j) {
              Some(indices_with_bppt_4_sta_ss) => {
                for &l in indices_with_bppt_4_sta_ss.keys() {
                  let app_2 = (l, j);
                  let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
                  if log_sta_ppf_4_bpa.is_finite() {
                    let ep_of_term_4_log_prob = log_sta_ppf_2 + log_sta_ppf_4_bpa + log_sta_ppf_matrices.log_sta_ppf_matrix[1][app_1.0 - 1][1][app_2.0 - 1]
                    - log_sta_ppf;
                    if ep_of_term_4_log_prob.is_finite() {
                      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                    }
                  }
                }
              }, None => {},
            }
            match bppmp_with_bppt_4_sta_ss.1.get(&j) {
              Some(bpps_with_bppt_4_sta_ss) => {
                for &l in bpps_with_bppt_4_sta_ss.keys() {
                  let app_2 = (j, l);
                  let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
                  if log_sta_ppf_4_bpa.is_finite() {
                    let ep_of_term_4_log_prob = log_sta_ppf_matrices.log_sta_ppf_matrix[1][app_1.0 - 1][1][app_2.0 - 1]
                    + log_sta_ppf_4_bpa
                    + log_sta_ppf_matrices.log_sta_ppf_matrix[app_1.1 + 1][slp.0][app_2.1 + 1][slp.1]
                    - log_sta_ppf;
                    if ep_of_term_4_log_prob.is_finite() {
                      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                    }
                  }
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      let log_sta_ppf_2 = log_sta_ppf_matrices.log_sta_ppf_matrix[1][i - 1][1][j - 1];
      match bppmp_with_bppt_4_sta_ss.0.get(&i) {
        Some(bpps_with_bppt_4_sta_ss) => {
          for &k in bpps_with_bppt_4_sta_ss.keys() {
            let app_1 = (i, k);
            match imp_with_bppt_4_sta_ss_1.1.get(&j) {
              Some(indices_with_bppt_4_sta_ss) => {
                for &l in indices_with_bppt_4_sta_ss.keys() {
                  let app_2 = (l, j);
                  let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
                  if log_sta_ppf_4_bpa.is_finite() {
                    let ep_of_term_4_log_prob = log_sta_ppf_matrices.log_sta_ppf_matrix[1][app_1.0 - 1][1][app_2.0 - 1]
                    + log_sta_ppf_4_bpa + log_sta_ppf_matrices.log_sta_ppf_matrix[app_1.1 + 1][slp.0][app_2.1 + 1][slp.1]
                    - log_sta_ppf;
                    if ep_of_term_4_log_prob.is_finite() {
                      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                    }
                  }
                }
              }, None => {},
            }
            match bppmp_with_bppt_4_sta_ss.1.get(&j) {
              Some(bpps_with_bppt_4_sta_ss) => {
                for &l in bpps_with_bppt_4_sta_ss.keys() {
                  let app_2 = (j, l);
                  let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
                  if log_sta_ppf_4_bpa.is_finite() {
                    let ep_of_term_4_log_prob = log_sta_ppf_2 + log_sta_ppf_4_bpa + log_sta_ppf_matrices.log_sta_ppf_matrix[app_1.1 + 1][slp.0][app_2.1 + 1][slp.1]
                    - log_sta_ppf;
                    if ep_of_term_4_log_prob.is_finite() {
                      if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                      eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                    }
                  }
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      for k in 1 .. i {
        for l in 1 .. j {
          let ep_of_term_4_log_prob = log_pseudo_prob_matrix[k][i][l][j] + get_ba_score(&(i, j), log_bap_matrix, sta_sps) + log_sta_ppf_matrices.log_sta_ppf_matrix[k + 1][i - 1][l + 1][j - 1];
          if ep_of_term_4_log_prob.is_finite() {
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
          match imp_with_bppt_4_sta_ss_2.0.get(&i) {
            Some(indices_with_bppt_4_sta_ss) => {
              for &s in indices_with_bppt_4_sta_ss.keys() {
                if k < s {
                  let app_1 = (min(s, i), max(s, i));
                  match imp_with_bppt_4_sta_ss_2.1.get(&j) {
                    Some(indices_with_bppt_4_sta_ss) => {
                      for &t in indices_with_bppt_4_sta_ss.keys() {
                        if l < t {
                          let app_2 = (min(t, j), max(t, j));
                          let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
                          if log_sta_ppf_4_bpa.is_finite() {
                            let ep_of_term_4_log_prob = log_pseudo_prob_matrix[k][app_1.1][l][app_2.1] + log_sta_ppf_matrices.log_sta_ppf_matrix[k + 1][app_1.0 - 1][l + 1][app_2.0 - 1] + log_sta_ppf_4_bpa;
                            if ep_of_term_4_log_prob.is_finite() {
                              if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                              eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                            }
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
      log_stapmp.log_bap_matrix[i][j] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
    }
  }
  log_pseudo_prob_matrix[0][slp.0 + 1][0][slp.1 + 1] = 0.;
  log_stapmp
}
