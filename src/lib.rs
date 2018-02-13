extern crate bio_seq_algos;
extern crate rna_algos;
extern crate itertools;

use bio_seq_algos::utils::fast_ln;
pub use rna_algos::utils::*;
use itertools::multizip;

pub type ThreeDPm = Vec<ProbMatrix>;
pub type FourDPm = Vec<ThreeDPm>;
pub struct StapMp {
  pub bpap_matrix: FourDPm, 
  pub bap_matrix: ProbMatrix,
}
pub type Log3DPm = Vec<LogProbMatrix>;
pub type Log4DPm = Vec<Log3DPm>;
pub struct LogStapMp {
  pub log_bpap_matrix: Log4DPm, 
  pub log_bap_matrix: LogProbMatrix,
}
type LogPpf3DMatrix = Vec<LogPpfMatrix>;
type LogPpf4DMatrix = Vec<LogPpf3DMatrix>;
pub struct LogStaPpfMatrices {
  pub log_sta_ppf_matrix: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_rightmost_bpas: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_bpas: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_2_loop_deletions_1: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_2lds_2: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_1_or_more_accessible_loops: LogPpf4DMatrix,
  pub log_sta_ppf_matrix_4_1_al: LogPpf4DMatrix,
  pub log_sa_ppf_matrices: LogSaPpfMatrices,
}
pub struct LogSaPpfMatrices {
  pub log_sa_ppf_matrix: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_ba_and_right_ba: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_base_gap_and_right_ba_1: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_bg_and_right_ba_2: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_ba_and_right_bg_1: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_ba_and_right_bg_2: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_bg_and_right_bg_1: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_bg_and_right_bg_2: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_bg_and_right_bg_3: LogPpf4DMatrix,
  pub log_sa_ppf_matrix_4_left_bg_and_right_bg_4: LogPpf4DMatrix,
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

impl LogStaPpfMatrices {
  fn new(slp: &(usize, usize)) -> LogStaPpfMatrices {
    let ni_matrix = vec![vec![vec![vec![NEG_INFINITY; slp.1]; slp.1]; slp.0]; slp.0];
    LogStaPpfMatrices {
      log_sta_ppf_matrix: vec![vec![vec![vec![0.; slp.1]; slp.1]; slp.0]; slp.0],
      log_sta_ppf_matrix_4_rightmost_bpas: ni_matrix.clone(),
      log_sta_ppf_matrix_4_bpas: ni_matrix.clone(),
      log_sta_ppf_matrix_4_2_loop_deletions_1: ni_matrix.clone(),
      log_sta_ppf_matrix_4_2lds_2: ni_matrix.clone(),
      log_sta_ppf_matrix_4_1_or_more_accessible_loops: ni_matrix.clone(),
      log_sta_ppf_matrix_4_1_al: ni_matrix,
      log_sa_ppf_matrices: LogSaPpfMatrices::new(slp),
    }
  }
}

impl LogSaPpfMatrices {
  fn new(slp: &(usize, usize)) -> LogSaPpfMatrices {
    let ni_matrix = vec![vec![vec![vec![NEG_INFINITY; slp.1]; slp.1]; slp.0]; slp.0];
    LogSaPpfMatrices {
      log_sa_ppf_matrix: vec![vec![vec![vec![0.; slp.1]; slp.1]; slp.0]; slp.0],
      log_sa_ppf_matrix_4_left_ba_and_right_ba: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_base_gap_and_right_ba_1: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_bg_and_right_ba_2: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_ba_and_right_bg_1: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_ba_and_right_bg_2: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_bg_and_right_bg_1: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_bg_and_right_bg_2: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_bg_and_right_bg_3: ni_matrix.clone(),
      log_sa_ppf_matrix_4_left_bg_and_right_bg_4: ni_matrix,
    }
  }
}

impl LogStapMp {
  fn new(slp: &(usize, usize)) -> LogStapMp {
    LogStapMp {
      log_bpap_matrix: vec![vec![vec![vec![NEG_INFINITY; slp.1]; slp.1]; slp.0]; slp.0],
      log_bap_matrix: vec![vec![NEG_INFINITY; slp.1]; slp.0]
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
pub fn io_algo_4_rna_stap_mp(sp: &SsPair, log_bpp_mp: &LpmPair, log_bap_matrix: &LogProbMatrix, bpp_mp: &PmPair, nbpp_mp: &PmPair, log_nbpp_mp: &LpmPair, sta_sps: &StaScoringParams) -> StapMp {
  let seq_len_pair = (sp.0.len(), sp.1.len());
  let log_sta_ppf_matrices = get_log_sta_ppf_matrices(sp, &seq_len_pair, log_bpp_mp, log_bap_matrix, bpp_mp, nbpp_mp, log_nbpp_mp, sta_sps);
  let log_stap_mp = get_log_stap_mp(sp, &log_sta_ppf_matrices, &seq_len_pair, log_bpp_mp, log_bap_matrix, bpp_mp, nbpp_mp, log_nbpp_mp, sta_sps);
  get_stap_mp(&log_stap_mp)
}

#[inline]
fn get_stap_mp(log_stap_mp: &LogStapMp) -> StapMp {
  StapMp {
    bpap_matrix: log_stap_mp.log_bpap_matrix.iter().map(|three_d_matrix| three_d_matrix.iter().map(|two_d_matrix| two_d_matrix.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect()).collect()).collect(),
    bap_matrix: log_stap_mp.log_bap_matrix.iter().map(|xs| xs.iter().map(|&x| x.exp()).collect()).collect(),
  }
}

#[inline]
pub fn get_log_sta_ppf_matrices(sp: &SsPair, slp: &(usize, usize), log_bpp_mp: &LpmPair, log_bap_matrix: &LogProbMatrix, bpp_mp: &PmPair, nbpp_mp: &PmPair, log_nbpp_mp: &LpmPair, sta_sps: &StaScoringParams) -> LogStaPpfMatrices {
  let mut log_sta_ppf_matrices = LogStaPpfMatrices::new(slp);
  for i in 0 .. slp.0 {
    for j in i .. slp.0 {
      for k in 0 .. slp.1 - 1 {
        log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][j][k + 1][k] = sta_sps.base_opening_gap_penalty + get_begp(&(i + 1, j), sta_sps);
        log_sta_ppf_matrices.log_sta_ppf_matrix[i][j][k + 1][k] = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][j][k + 1][k];
      }
    }
  }
  for k in 0 .. slp.1 {
    for l in k .. slp.1 {
      for i in 0 .. slp.0 - 1 {
        log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i + 1][i][k][l] = sta_sps.base_opening_gap_penalty + get_begp(&(k + 1, l), sta_sps);
        log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][i][k][l] = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i + 1][i][k][l];
      }
    }
  }
  for i in 0 .. slp.0 {
    for j in 0 .. slp.1 {
      log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][i][j][j] = get_ba_score(&(i, j), log_bap_matrix, sta_sps);
        log_sta_ppf_matrices.log_sta_ppf_matrix[i][i][j][j] = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][i][j][j];
    }
  }
  for i in 0 .. slp.0 - 1 {
    for j in 0 .. slp.1 {
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = get_ba_score(&(i, j), log_bap_matrix, sta_sps);
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = get_ba_score(&(i + 1, j), log_bap_matrix, sta_sps);
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][i + 1][j][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + sta_sps.base_opening_gap_penalty;
      log_sta_ppf_matrices.log_sta_ppf_matrix[i][i + 1][j][j] = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][i + 1][j][j];
    }
  }
  for i in 0 .. slp.0 {
    for j in 0 .. slp.1 - 1 {
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = get_ba_score(&(i, j), log_bap_matrix, sta_sps);
      eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
      let ep_of_term_4_log_pf = get_ba_score(&(i, j + 1), log_bap_matrix, sta_sps);
      if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
      log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][i][j][j + 1] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + sta_sps.base_opening_gap_penalty;
      log_sta_ppf_matrices.log_sta_ppf_matrix[i][i][j][j + 1] = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][i][j][j + 1];
    }
  }
  for i in 0 .. slp.0 - 2 {
    for j in i + 2 .. slp.0 {
      for k in 0 .. slp.1 {
        let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
        let mut max_ep_of_term_4_log_pf = get_ba_score(&(i, k), log_bap_matrix, sta_sps) + get_begp(&(i + 2, j), sta_sps);
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
        let ep_of_term_4_log_pf = get_ba_score(&(j, k), log_bap_matrix, sta_sps) + get_begp(&(i, j - 2), sta_sps);
        if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        for l in i + 1 .. j {
          let ep_of_term_4_log_pf = get_ba_score(&(l, k), log_bap_matrix, sta_sps) + sta_sps.base_opening_gap_penalty + get_begp(&(i + 1, l - 1), sta_sps) + get_begp(&(l + 2, j), sta_sps);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][j][k][k] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + sta_sps.base_opening_gap_penalty;
        log_sta_ppf_matrices.log_sta_ppf_matrix[i][j][k][k] = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][j][k][k];
      }
    }
  }
  for i in 0 .. slp.1 - 2 {
    for j in i + 2 .. slp.1 {
      for k in 0 .. slp.0 {
        let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
        let mut max_ep_of_term_4_log_pf = get_ba_score(&(k, i), log_bap_matrix, sta_sps) + get_begp(&(i + 2, j), sta_sps);
        eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
        let ep_of_term_4_log_pf = get_ba_score(&(k, j), log_bap_matrix, sta_sps) + get_begp(&(i, j - 2), sta_sps);
        if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
        eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        for l in i + 1 .. j {
          let ep_of_term_4_log_pf = get_ba_score(&(k, l), log_bap_matrix, sta_sps) + sta_sps.base_opening_gap_penalty + get_begp(&(i + 1, l - 1), sta_sps) + get_begp(&(l + 2, j), sta_sps);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
        }
        log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[k][k][i][j] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + sta_sps.base_opening_gap_penalty;
        log_sta_ppf_matrices.log_sta_ppf_matrix[k][k][i][j] = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[k][k][i][j];
      }
    }
  }
  for sub_seq_len_1 in 2 .. slp.0 + 1 {
    for i in 0 .. slp.0 - sub_seq_len_1 + 1 {
      let j = i + sub_seq_len_1 - 1;
      let pp_closing_loop_1 = (i, j);
      let bp_closing_loop_1 = (sp.0[pp_closing_loop_1.0], sp.0[pp_closing_loop_1.1]);
      for ssl_2 in 2 .. slp.1 + 1 {
        for k in 0 .. slp.1 - ssl_2 + 1 {
          let l = k + ssl_2 - 1;
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i][j][k][l] = get_ba_score(&(i, k), log_bap_matrix, sta_sps) + get_ba_score(&(j, l), log_bap_matrix, sta_sps) + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i + 1][j - 1][k + 1][l - 1];
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = sta_sps.base_opening_gap_penalty + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i + 1][j - 1][k][l - 1];
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_base_gap_and_right_ba_1[i + 1][j - 1][k][l - 1];
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_1[i + 1][j - 1][k][l - 1];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_3[i + 1][j - 1][k][l - 1];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = sta_sps.base_extending_gap_penalty + logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          let ba_score = get_ba_score(&(j, l), log_bap_matrix, sta_sps);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_base_gap_and_right_ba_1[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + ba_score;
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = sta_sps.base_opening_gap_penalty + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i][j - 1][k + 1][l - 1];
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_ba_2[i][j - 1][k + 1][l - 1];
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_2[i][j - 1][k + 1][l - 1];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_4[i][j - 1][k + 1][l - 1];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = sta_sps.base_extending_gap_penalty + logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_ba_2[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf) + ba_score;
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i + 1][j - 1][k + 1][l] + sta_sps.base_opening_gap_penalty;
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_1[i + 1][j - 1][k + 1][l];
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_1[i + 1][j - 1][k + 1][l];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_2[i + 1][j - 1][k + 1][l];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf) + sta_sps.base_extending_gap_penalty;
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          let ba_score = get_ba_score(&(i, k), log_bap_matrix, sta_sps);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_1[i][j][k][l] = ba_score + logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i + 1][j][k + 1][l - 1] + sta_sps.base_opening_gap_penalty;
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_2[i + 1][j][k + 1][l - 1];
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_3[i + 1][j][k + 1][l - 1];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_4[i + 1][j][k + 1][l - 1];
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf) + sta_sps.base_extending_gap_penalty;
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_1[i][j][k][l] = ba_score + logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = 2. * sta_sps.base_extending_gap_penalty + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_1[i + 1][j - 1][k][l];
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i + 1][j - 1][k][l] + sta_sps.base_opening_gap_penalty;
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_base_gap_and_right_ba_1[i + 1][j - 1][k][l] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_1[i + 1][j - 1][k][l] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = sta_sps.base_opening_gap_penalty + logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_1[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = 2. * sta_sps.base_extending_gap_penalty + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_2[i][j - 1][k + 1][l];
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i][j - 1][k + 1][l] + sta_sps.base_opening_gap_penalty;
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_ba_2[i][j - 1][k + 1][l] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_1[i][j - 1][k + 1][l] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = sta_sps.base_opening_gap_penalty + logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_2[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = 2. * sta_sps.base_extending_gap_penalty + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_3[i + 1][j][k][l - 1];
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i + 1][j][k][l - 1] + sta_sps.base_opening_gap_penalty;
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_base_gap_and_right_ba_1[i + 1][j][k][l - 1] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_2[i + 1][j][k][l - 1] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = sta_sps.base_opening_gap_penalty + logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_3[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = 2. * sta_sps.base_extending_gap_penalty + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_4[i][j][k + 1][l - 1];
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          let mut eps_of_sub_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i][j][k + 1][l - 1] + sta_sps.base_opening_gap_penalty;
          eps_of_sub_terms_4_log_pf.push(max_ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_ba_2[i][j][k + 1][l - 1] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_sub_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_2[i][j][k + 1][l - 1] + sta_sps.base_extending_gap_penalty;
          if max_ep_of_sub_term_4_log_pf < ep_of_sub_term_4_log_pf {max_ep_of_sub_term_4_log_pf = ep_of_sub_term_4_log_pf;}
          eps_of_sub_terms_4_log_pf.push(ep_of_sub_term_4_log_pf);
          let ep_of_term_4_log_pf = sta_sps.base_opening_gap_penalty + logsumexp(&eps_of_sub_terms_4_log_pf[..], max_ep_of_sub_term_4_log_pf);
          if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_4[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let eps_of_terms_4_log_pf = [
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_ba[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_base_gap_and_right_ba_1[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_ba_2[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_1[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_ba_and_right_bg_2[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_1[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_2[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_3[i][j][k][l],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix_4_left_bg_and_right_bg_4[i][j][k][l],
          ];
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
            if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
          }
          log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let pp_closing_loop_2 = (k, l);
          let bp_closing_loop_2 = (sp.1[pp_closing_loop_2.0], sp.1[pp_closing_loop_2.1]);
          if is_canonical(&bp_closing_loop_1) && is_canonical(&bp_closing_loop_2) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[pp_closing_loop_1.0 + 1][pp_closing_loop_1.1 - 1][pp_closing_loop_2.0 + 1][pp_closing_loop_2.1 - 1];
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            for h in k + 1 .. l - 1 {
              let ep_of_term_4_log_pf = get_legp(&(pp_closing_loop_2.0 + 1, h - 1), sta_sps) + log_sta_ppf_matrices.log_sta_ppf_matrix_4_2_loop_deletions_1[pp_closing_loop_1.0 + 1][pp_closing_loop_1.1 - 1][h][pp_closing_loop_2.1 - 1];
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
            for h in i + 1 .. j - 1 {
              let ep_of_term_4_log_pf = get_legp(&(pp_closing_loop_1.0 + 1, h - 1), sta_sps) + log_sta_ppf_matrices.log_sta_ppf_matrix_4_2_loop_deletions_1[h][pp_closing_loop_1.1 - 1][pp_closing_loop_2.0 + 1][pp_closing_loop_2.1 - 1];
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
            for h_1 in i + 1 .. j - 1 {
              for h_2 in k + 1 .. l - 1 {
                let log_sta_ppf_matrix_4_1_al = log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_al[h_1][pp_closing_loop_1.1 - 1][h_2][pp_closing_loop_2.1 - 1];
                let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_or_more_accessible_loops[pp_closing_loop_1.0 + 1][h_1 - 1][pp_closing_loop_2.0 + 1][h_2 - 1] + log_sta_ppf_matrix_4_1_al;
                if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[pp_closing_loop_1.0 + 1][h_1 - 1][pp_closing_loop_2.0 + 1][h_2 - 1] + log_sta_ppf_matrix_4_1_al;
                if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              }
            }
            if eps_of_terms_4_log_pf.len() > 0 {
              log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[pp_closing_loop_1.0][pp_closing_loop_1.1][pp_closing_loop_2.0][pp_closing_loop_2.1] = get_bpa_score(&pp_closing_loop_1, &pp_closing_loop_2, log_bpp_mp, log_bap_matrix, bpp_mp, nbpp_mp, log_nbpp_mp, sta_sps)
              + logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
            }
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          for h_1 in i + 1 .. j + 1 {
            for h_2 in k + 1 .. l + 1 {
              let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[i][h_1][k][h_2];
              if log_sta_ppf_4_bpa.is_finite() {
                let ep_of_term_4_log_pf = log_sta_ppf_4_bpa
                + if h_1 > slp.0 - 2 && h_2 > slp.1 - 2 {0.}
                else if h_1 > slp.0 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(h_2 + 2, l), sta_sps)}
                else if h_2 > slp.1 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(h_1 + 2, j), sta_sps)}
                else {log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[h_1 + 1][j][h_2 + 1][l]};
                if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              }
            }
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_rightmost_bpas[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][j][k][l];
          eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          for h_1 in i .. j {
            for h_2 in k .. l {
              let ep_of_term_4_log_pf = if h_1 < 1 && h_2 < 1 {0.}
              else if h_1 < 1 {
                if k > h_2 - 1 {0.}
                else {sta_sps.base_opening_gap_penalty + get_begp(&(k + 1, h_2 - 1), sta_sps)}
              }
              else if h_2 < 1 {
                if i > h_1 - 1 {0.}
                else {sta_sps.base_opening_gap_penalty + get_begp(&(i + 1, h_1 - 1), sta_sps)}
              }
              else {log_sta_ppf_matrices.log_sta_ppf_matrix[i][h_1 - 1][k][h_2 - 1]}
              + log_sta_ppf_matrices.log_sta_ppf_matrix_4_rightmost_bpas[h_1][j][h_2][l];
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
          }
          log_sta_ppf_matrices.log_sta_ppf_matrix[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          for h in k + 1 .. l + 1 {
            let pp_cl = (k, h);
            let bp_cl = (sp.1[pp_cl.0], sp.1[pp_cl.1]);
            if is_canonical(&bp_cl) {
              let ep_of_term_4_log_pf = get_logp(&pp_cl, &log_bpp_mp.1, sta_sps) + log_sta_ppf_matrices.log_sta_ppf_matrix[i][j][k + 1][h - 1] + get_legp(&(h + 1, l), sta_sps);
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_2_loop_deletions_1[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          for h in i + 1 .. j + 1 {
            let pp_cl = (i, h);
            let bp_cl = (sp.0[pp_cl.0], sp.0[pp_cl.1]);
            if is_canonical(&bp_cl) {
              let ep_of_term_4_log_pf = get_logp(&pp_cl, &log_bpp_mp.0, sta_sps) + log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][h - 1][k][l] + get_legp(&(h + 1, j), sta_sps);
              if max_ep_of_term_4_log_pf < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            }
          }
          if eps_of_terms_4_log_pf.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_2lds_2[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf);
          }
          let mut eps_of_terms_4_log_pf_1 = EpsOfTerms4LogPf::new();
          let mut eps_of_terms_4_log_pf_2 = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf_1 = NEG_INFINITY;
          let mut max_ep_of_term_4_log_pf_2 = NEG_INFINITY;
          for h_1 in i + 1 .. j + 1 {
            let pp_cl_1 = (i, h_1);
            for h_2 in k + 1 .. l + 1 {
              let log_sta_ppf_4_1_al = log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_al[h_1][j][h_2][l];
              let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[i][h_1 - 1][k][h_2 - 1] + log_sta_ppf_4_1_al;
              if max_ep_of_term_4_log_pf_1 < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf_1 = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf_1.push(ep_of_term_4_log_pf);
              let ep_of_term_4_log_pf = log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_or_more_accessible_loops[i][h_1 - 1][k][h_2 - 1] + log_sta_ppf_4_1_al;
              if max_ep_of_term_4_log_pf_1 < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf_1 = ep_of_term_4_log_pf;}
              eps_of_terms_4_log_pf_1.push(ep_of_term_4_log_pf);
              let pp_cl_2 = (k, h_2);
              let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[pp_cl_1.0][pp_cl_1.1][pp_cl_2.0][pp_cl_2.1];
              if log_sta_ppf_4_bpa.is_finite() {
                let ep_of_term_4_log_pf = log_sta_ppf_4_bpa + if pp_cl_1.1 > slp.0 - 2 && pp_cl_2.1 > slp.1 - 2 {0.}
                else if pp_cl_1.1 > slp.0 - 2 {
                  if pp_cl_2.1 + 1 > l {0.}
                  else {sta_sps.base_opening_gap_penalty + get_begp(&(pp_cl_2.1 + 2, l), sta_sps)}
                }
                else if pp_cl_2.1 > slp.1 - 2 {
                  if pp_cl_1.1 + 1 > j {0.}
                  else {sta_sps.base_opening_gap_penalty + get_begp(&(pp_cl_1.1 + 2, j), sta_sps)}
                }
                else {log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[pp_cl_1.1 + 1][j][pp_cl_2.1 + 1][l]};
                if max_ep_of_term_4_log_pf_2 < ep_of_term_4_log_pf {max_ep_of_term_4_log_pf_2 = ep_of_term_4_log_pf;}
                eps_of_terms_4_log_pf_2.push(ep_of_term_4_log_pf);
              }
            }
          }
          if eps_of_terms_4_log_pf_1.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_or_more_accessible_loops[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf_1[..], max_ep_of_term_4_log_pf_1);
          }
          if eps_of_terms_4_log_pf_2.len() > 0 {
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_al[i][j][k][l] = logsumexp(&eps_of_terms_4_log_pf_2[..], max_ep_of_term_4_log_pf_2);
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
fn get_logp(pp: &PosPair, log_bpp_matrix: &LogProbMatrix, sta_sps: &StaScoringParams) -> StaScore {
  log_bpp_matrix[pp.0][pp.1] - sta_sps.log_null_hypothesis_bpp + 2. * sta_sps.loop_opening_gap_penalty
}

#[inline]
fn get_legp(pp: &PosPair, sta_sps: &StaScoringParams) -> StaScore {
  (pp.1 + 1 - pp.0) as StaScore * sta_sps.loop_extending_gap_penalty
}

#[inline]
fn get_ba_score(pp: &PosPair, log_bap_matrix: &LogProbMatrix, sta_sps: &StaScoringParams) -> StaScore {
  log_bap_matrix[pp.0][pp.1] - sta_sps.log_nh_bap
}

#[inline]
fn get_bpa_score(pp_1: &PosPair, pp_2: &PosPair, log_bpp_mp: &LpmPair, log_bap_matrix: &LogProbMatrix, bpp_mp: &PmPair, nbpp_mp: &PmPair, log_nbpp_mp: &LpmPair, sta_sps: &StaScoringParams) -> StaScore {
  let bin_pd_1 = [bpp_mp.0[pp_1.0][pp_1.1], nbpp_mp.0[pp_1.0][pp_1.1]];
  let bin_pd_2 = [bpp_mp.1[pp_2.0][pp_2.1], nbpp_mp.1[pp_2.0][pp_2.1]];
  let bin_lpd_1 = [log_bpp_mp.0[pp_1.0][pp_1.1], log_nbpp_mp.0[pp_1.0][pp_1.1]];
  let bin_lpd_2 = [log_bpp_mp.1[pp_2.0][pp_2.1], log_nbpp_mp.1[pp_2.0][pp_2.1]];
  log_bpp_mp.0[pp_1.0][pp_1.1] + log_bpp_mp.1[pp_2.0][pp_2.1] - 2. * sta_sps.log_null_hypothesis_bpp
  + get_ba_score(&(pp_1.0, pp_2.0), log_bap_matrix, sta_sps) + get_ba_score(&(pp_1.1, pp_2.1), log_bap_matrix, sta_sps)
  + sta_sps.scaling_parameter_4_ss_distance * get_jsd(&bin_pd_1[..], &bin_pd_2[..], &bin_lpd_1[..], &bin_lpd_2[..])
}

#[inline]
fn get_jsd(pd_1: PdSlice, pd_2: PdSlice, lpd_1: LpdSlice, lpd_2: LpdSlice) -> JensenShannonDist {
  let mid_pd = pd_1.iter().zip(pd_2).map(|(&x, &y)| (x + y) / 2.).collect::<ProbDist>();
  let log_mid_prob_dist = mid_pd.iter().map(|&x| fast_ln(x)).collect::<LogProbDist>();
  ((get_kld(pd_1, lpd_1, &log_mid_prob_dist)
  + get_kld(pd_2, lpd_2, &log_mid_prob_dist)) / 2.).sqrt()
}

#[inline]
fn get_kld(pd: PdSlice, lpd_1: LpdSlice, lpd_2: LpdSlice) -> JensenShannonDist {
  multizip((pd.iter(), lpd_1.iter(), lpd_2.iter())).map(|(&x, &y, &z)| x * (y - z)).sum()
}

#[inline]
fn get_log_stap_mp(sp: &SsPair, log_sta_ppf_matrices: &LogStaPpfMatrices, slp: &(usize, usize), log_bpp_mp: &LpmPair, log_bap_matrix: &LogProbMatrix, bpp_mp: &PmPair, nbpp_mp: &PmPair, log_nbpp_mp: &LpmPair, sta_sps: &StaScoringParams) -> LogStapMp {
  let mut log_stap_mp = LogStapMp::new(slp);
  let log_sta_ppf = log_sta_ppf_matrices.log_sta_ppf_matrix[0][slp.0 - 1][0][slp.1 - 1];
  let mut log_prob_matrix_4_ils_1 = log_stap_mp.log_bpap_matrix.clone();
  let mut log_prob_matrix_4_ils_2 = log_stap_mp.log_bpap_matrix.clone();
  for ssl_1 in (2 .. slp.0 + 1).rev() {
    for i in 0 .. slp.0 - ssl_1 + 1 {
      let j = i + ssl_1 - 1;
      let accessible_pp_1 = (i, j);
      for ssl_2 in (2 .. slp.1 + 1).rev() {
        for k in 0 .. slp.1 - ssl_2 + 1 {
          let l = k + ssl_2 - 1;
          let app_2 = (k, l);
          let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[accessible_pp_1.0][accessible_pp_1.1][app_2.0][app_2.1];
          if log_sta_ppf_4_bpa.is_finite() {
            let mut eps_of_terms_4_log_prob = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_prob = if accessible_pp_1.0 < 1 && app_2.0 < 1 {0.}
            else if accessible_pp_1.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, app_2.0 - 1), sta_sps)}
            else if app_2.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, accessible_pp_1.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[0][accessible_pp_1.0 - 1][0][app_2.0 - 1]}
            + if accessible_pp_1.1 > slp.0 - 2 {0.}
            else if accessible_pp_1.1 < slp.0 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(app_2.1 + 2, slp.1 - 1), sta_sps)}
            else if app_2.1 > slp.1 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(accessible_pp_1.1 + 2, slp.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[accessible_pp_1.1 + 1][slp.0 - 1][app_2.1 + 1][slp.1 - 1]}
            - log_sta_ppf;
            eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
            for h_1 in 0 .. i {
              for h_2 in 0 .. k {
                let eps_of_sub_terms_4_log_prob_1 = [log_prob_matrix_4_ils_1[h_1][accessible_pp_1.1][h_2][app_2.1], log_prob_matrix_4_ils_2[h_1][accessible_pp_1.1][h_2][app_2.1]];
                let max_ep_of_sub_term_4_log_prob_1 = if eps_of_sub_terms_4_log_prob_1[0] > eps_of_sub_terms_4_log_prob_1[1] {eps_of_sub_terms_4_log_prob_1[0]} else {eps_of_sub_terms_4_log_prob_1[1]};
                let eps_of_sub_terms_4_log_prob_2 = [
                  log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_or_more_accessible_loops[h_1 + 1][accessible_pp_1.0 - 1][h_2 + 1][app_2.0 - 1],
                  log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[h_1 + 1][accessible_pp_1.0 - 1][h_2 + 1][app_2.0 - 1],
                ];
                let max_ep_of_sub_term_4_log_prob_2 = if eps_of_sub_terms_4_log_prob_2[0] > eps_of_sub_terms_4_log_prob_2[1] {eps_of_sub_terms_4_log_prob_2[0]} else {eps_of_sub_terms_4_log_prob_2[1]};
                let ep_of_term_4_log_prob = logsumexp(&eps_of_sub_terms_4_log_prob_1[..], max_ep_of_sub_term_4_log_prob_1) + logsumexp(&eps_of_sub_terms_4_log_prob_2[..], max_ep_of_sub_term_4_log_prob_2);
                if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
              }
            }
            log_stap_mp.log_bpap_matrix[accessible_pp_1.0][accessible_pp_1.1][app_2.0][app_2.1] = log_sta_ppf_4_bpa + logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
          }
          let mut eps_of_terms_4_log_prob_1 = EpsOfTerms4LogPf::new();
          let mut eps_of_terms_4_log_prob_2 = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_prob_1 = NEG_INFINITY;
          let mut max_ep_of_term_4_log_prob_2 = NEG_INFINITY;
          for h_1 in j + 1 .. slp.0 {
            let app_1 = (i, h_1);
            for h_2 in l + 1 .. slp.1 {
              let app_2 = (k, h_2);
              let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
              if log_sta_ppf_4_bpa.is_finite() {
                let log_bpap = log_stap_mp.log_bpap_matrix[app_1.0][app_1.1][app_2.0][app_2.1];
                let bpa_score = get_bpa_score(&app_1, &app_2, log_bpp_mp, log_bap_matrix, bpp_mp, nbpp_mp, log_nbpp_mp, sta_sps);
                let ep_of_term_4_log_prob = log_bpap + bpa_score
                + log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_or_more_accessible_loops[j + 1][app_1.1 - 1][l + 1][app_2.1 - 1]
                - log_sta_ppf_4_bpa;
                if max_ep_of_term_4_log_prob_1 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_1 = ep_of_term_4_log_prob;}
                eps_of_terms_4_log_prob_1.push(ep_of_term_4_log_prob);
                let ep_of_term_4_log_prob = log_bpap + bpa_score
                + log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[j + 1][app_1.1 - 1][l + 1][app_2.1 - 1]
                - log_sta_ppf_4_bpa;
                if max_ep_of_term_4_log_prob_2 < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob_2 = ep_of_term_4_log_prob;}
                eps_of_terms_4_log_prob_2.push(ep_of_term_4_log_prob);
              }
            }
          }
          if log_prob_matrix_4_ils_1.len() > 0 {
            log_prob_matrix_4_ils_1[i][j][k][l] = logsumexp(&eps_of_terms_4_log_prob_1[..], max_ep_of_term_4_log_prob_1);
            log_prob_matrix_4_ils_2[i][j][k][l] = logsumexp(&eps_of_terms_4_log_prob_2[..], max_ep_of_term_4_log_prob_2);
          }
        }
      }
    }
  }
  for i in 0 .. slp.0 {
    for j in 0 .. slp.1 {
      let mut eps_of_terms_4_log_prob = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_prob = if i < 1 && j < 1 {0.}
      else if i < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, j - 1), sta_sps)}
      else if j < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, i - 1), sta_sps)}
      else {log_sta_ppf_matrices.log_sta_ppf_matrix[0][i - 1][0][j - 1]}
      + get_ba_score(&(i, j), log_bap_matrix, sta_sps)
      + if i > slp.0 - 2 && j > slp.1 - 2 {0.}
      else if i > slp.0 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(j + 2, slp.1 - 1), sta_sps)}
      else if j > slp.1 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(i + 2, slp.0 - 1), sta_sps)}
      else {log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][slp.0 - 1][j + 1][slp.1 - 1]} - log_sta_ppf;
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
      let log_sta_ppf_2 = if i > slp.0 - 2 && j > slp.1 - 2 {0.}
      else if i > slp.0 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(j + 2, slp.1 - 1), sta_sps)}
      else if j > slp.1 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(i + 2, slp.0 - 1), sta_sps)}
      else {log_sta_ppf_matrices.log_sta_ppf_matrix[i + 1][slp.0 - 1][j + 1][slp.1 - 1]};
      for k in 0 .. i {
        let app_1 = (k, i);
        for l in 0 .. j {
          let app_2 = (l, j);
          let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
          if log_sta_ppf_4_bpa.is_finite() {
            let ep_of_term_4_log_prob = log_sta_ppf_2 + log_sta_ppf_4_bpa + if app_1.0 < 1 && app_2.0 < 1 {0.}
            else if app_1.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, app_2.0 - 1), sta_sps)}
            else if app_2.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, app_1.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[0][app_1.0 - 1][0][app_2.0 - 1]}
            - log_sta_ppf;
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
        for l in j + 1 .. slp.1 {
          let app_2 = (j, l);
          let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
          if log_sta_ppf_4_bpa.is_finite() {
            let ep_of_term_4_log_prob = if app_1.0 < 1 && app_2.0 < 1 {0.}
            else if app_1.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, app_2.0 - 1), sta_sps)}
            else if app_2.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, app_1.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[0][app_1.0 - 1][0][app_2.0 - 1]}
            + log_sta_ppf_4_bpa
            + if app_1.1 > slp.0 - 2 && app_2.1 > slp.1 - 2 {0.}
            else if app_1.1 > slp.0 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(app_2.1 + 2, slp.1 - 1), sta_sps)}
            else if app_2.1 > slp.1 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(app_1.1 + 2, slp.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[app_1.1 + 1][slp.0 - 1][app_2.1 + 1][slp.1 - 1]}
            - log_sta_ppf;
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
      }
      let log_sta_ppf_2 = if i < 1 && j < 1 {0.}
      else if i < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, j - 1), sta_sps)}
      else if j < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, i - 1), sta_sps)}
      else {log_sta_ppf_matrices.log_sta_ppf_matrix[0][i - 1][0][j - 1]};
      for k in i + 1 .. slp.0 {
        let app_1 = (i, k);
        for l in 0 .. j {
          let app_2 = (l, j);
          let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
          if log_sta_ppf_4_bpa.is_finite() {
            let ep_of_term_4_log_prob = if app_1.0 < 1 && app_2.0 < 1 {0.}
            else if app_1.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, app_2.0 - 1), sta_sps)}
            else if app_2.0 < 1 {sta_sps.base_opening_gap_penalty + get_begp(&(1, app_1.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[0][app_1.0 - 1][0][app_2.0 - 1]}
            + log_sta_ppf_4_bpa + if app_1.1 > slp.0 - 2 && app_2.1 > slp.1 - 2 {0.}
            else if app_1.1 > slp.0 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(app_2.1 + 2, slp.1 - 1), sta_sps)}
            else if app_2.1 > slp.1 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(app_1.1 + 2, slp.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[app_1.1 + 1][slp.0 - 1][app_2.1 + 1][slp.1 - 1]}
            - log_sta_ppf;
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
        for l in j + 1 .. slp.1 {
          let app_2 = (j, l);
          let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
          if log_sta_ppf_4_bpa.is_finite() {
            let ep_of_term_4_log_prob = log_sta_ppf_2 + log_sta_ppf_4_bpa + if app_1.1 > slp.0 - 2 && app_2.1 > slp.1 - 2 {0.}
            else if app_1.1 > slp.0 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(app_2.1 + 2, slp.1 - 1), sta_sps)}
            else if app_2.1 > slp.1 - 2 {sta_sps.base_opening_gap_penalty + get_begp(&(app_1.1 + 2, slp.0 - 1), sta_sps)}
            else {log_sta_ppf_matrices.log_sta_ppf_matrix[app_1.1 + 1][slp.0 - 1][app_2.1 + 1][slp.1 - 1]}
            - log_sta_ppf;
            if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          }
        }
      }
      for k in 0 .. i {
        for l in 0 .. j {
          let eps_of_sub_terms_4_log_prob_1 = [log_prob_matrix_4_ils_1[k][i][l][j], log_prob_matrix_4_ils_2[k][i][l][j]];
          let max_ep_of_sub_term_4_log_prob_1 = if eps_of_sub_terms_4_log_prob_1[0] > eps_of_sub_terms_4_log_prob_1[1] {eps_of_sub_terms_4_log_prob_1[0]} else {eps_of_sub_terms_4_log_prob_1[1]};
          let eps_of_sub_terms_4_log_prob_2 = [
            log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_or_more_accessible_loops[k + 1][i - 1][l + 1][j - 1],
            log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[k + 1][i - 1][l + 1][j - 1],
          ];
          let max_ep_of_sub_term_4_log_prob_2 = if eps_of_sub_terms_4_log_prob_2[0] > eps_of_sub_terms_4_log_prob_2[1] {eps_of_sub_terms_4_log_prob_2[0]} else {eps_of_sub_terms_4_log_prob_2[1]};
          let ep_of_term_4_log_prob = get_ba_score(&(i, j), log_bap_matrix, sta_sps) + logsumexp(&eps_of_sub_terms_4_log_prob_1[..], max_ep_of_sub_term_4_log_prob_1) + logsumexp(&eps_of_sub_terms_4_log_prob_2[..], max_ep_of_sub_term_4_log_prob_2);
          if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
          eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          for s in k + 1 .. slp.0 - 1 {
            if s == i {continue;}
            let app_1 = (min(s, i), max(s, i));
            for t in l + 1 .. slp.1 - 1 {
              if t == j {continue;}
              let app_2 = (min(t, j), max(t, j));
              let log_sta_ppf_4_bpa = log_sta_ppf_matrices.log_sta_ppf_matrix_4_bpas[app_1.0][app_1.1][app_2.0][app_2.1];
              if log_sta_ppf_4_bpa.is_finite() {
                let eps_of_sub_terms_4_log_prob_1 = [log_prob_matrix_4_ils_1[k][app_1.1][l][app_2.1], log_prob_matrix_4_ils_2[k][app_1.1][l][app_2.1]];
                let max_ep_of_sub_term_4_log_prob_1 = if eps_of_sub_terms_4_log_prob_1[0] > eps_of_sub_terms_4_log_prob_1[1] {eps_of_sub_terms_4_log_prob_1[0]} else {eps_of_sub_terms_4_log_prob_1[1]};
                let eps_of_sub_terms_4_log_prob_2 = [
                  log_sta_ppf_matrices.log_sta_ppf_matrix_4_1_or_more_accessible_loops[k + 1][app_1.0 - 1][l + 1][app_2.0 - 1],
                  log_sta_ppf_matrices.log_sa_ppf_matrices.log_sa_ppf_matrix[k + 1][app_1.0 - 1][l + 1][app_2.0 - 1],
                ];
                let max_ep_of_sub_term_4_log_prob_2 = if eps_of_sub_terms_4_log_prob_2[0] > eps_of_sub_terms_4_log_prob_2[1] {eps_of_sub_terms_4_log_prob_2[0]} else {eps_of_sub_terms_4_log_prob_2[1]};
                let ep_of_term_4_log_prob = log_sta_ppf_4_bpa + logsumexp(&eps_of_sub_terms_4_log_prob_1[..], max_ep_of_sub_term_4_log_prob_1) + logsumexp(&eps_of_sub_terms_4_log_prob_2[..], max_ep_of_sub_term_4_log_prob_2);
                if max_ep_of_term_4_log_prob < ep_of_term_4_log_prob {max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;}
                eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
              }
            }
          }
        }
      }
      log_stap_mp.log_bap_matrix[i][j] = logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
    }
  }
  log_stap_mp
}
