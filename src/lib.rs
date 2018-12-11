extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate getopts;

pub mod utils;
pub mod stem_params;

pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use utils::*;
pub use stem_params::*;

pub type PosTriple = (Pos, Pos, Pos);
pub type PosQuadruple = (Pos, Pos, Pos, Pos);
pub type Prob1dMat = HashMap<Pos, Prob, Hasher>;
pub type SparseProbMat = HashMap<PosPair, Prob, Hasher>;
pub type Prob3dMat = HashMap<PosTriple, Prob, Hasher>;
pub type Prob4dMat = HashMap<PosQuadruple, Prob, Hasher>;
pub struct Stapmt {
  pub base_align_prob_mat: SparseProbMat,
  pub opening_gap_prob_mat_1: Probs,
  pub ogp_mat_2: Probs,
  pub extending_gap_prob_mat_1: Probs,
  pub egp_mat_2: Probs,
  pub base_pair_align_prob_mat: Prob4dMat,
}
pub type SparseLogProbMat = HashMap<PosPair, LogProb, Hasher>;
pub type LogProb3dMat = HashMap<PosTriple, LogProb, Hasher>;
pub type LogProb4dMat = HashMap<PosQuadruple, LogProb, Hasher>;
#[derive(Clone)]
pub struct LogStapmt {
  pub lbap_mat: SparseLogProbMat,
  pub logp_mat_1: LogProbs,
  pub logp_mat_2: LogProbs,
  pub legp_mat_1: LogProbs,
  pub legp_mat_2: LogProbs,
  pub lbpap_mat: LogProb4dMat,
}
type SparseLogPpfMat = HashMap<PosPair, LogPf, Hasher>;
type LogPpf4dMat = HashMap<PosQuadruple, LogPf, Hasher>;
#[derive(Clone)]
pub struct LogStaPpfMats {
  pub log_ppf_mat: SparseLogPpfMat,
  pub log_ppf_mat_4_bas: SparseLogPpfMat,
  pub log_ppf_mat_4_ogs_1: SparseLogPpfMat,
  pub log_ppf_mat_4_ogs_2: SparseLogPpfMat,
  pub log_ppf_mat_4_egs_1: SparseLogPpfMat,
  pub log_ppf_mat_4_egs_2: SparseLogPpfMat,
  pub log_ppf_mat_4_gaps_1: SparseLogPpfMat,
  pub log_ppf_mat_4_gaps_2: SparseLogPpfMat,
}
pub type StaFreeEnergy = LogProb;
pub struct StaFeParams {
  pub lstapmt: LogStapmt,
  pub lstapmt_on_random_assump: LogStapmt,
}
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type StaFeParamSetsWithRnaIdPairs = HashMap<RnaIdPair, StaFeParams, Hasher>;
type StapmtsWithRnaIdPairs = HashMap<RnaIdPair, Stapmt, Hasher>;
pub type LstapmtsWithRnaIdPairs = HashMap<RnaIdPair, LogStapmt, Hasher>;
pub type LogProbMats = Vec<SparseLogProbMat>;
pub type LogProbSeqs = Vec<LogProbs>;
type SeqsOfEpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, ExpPartOfTerm4LogProb, Hasher>;
type SeqsOfEpsOfTerms4LogProbsWithPosQuadruples = HashMap<PosQuadruple, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosQuadruples = HashMap<PosQuadruple, ExpPartOfTerm4LogProb, Hasher>;
type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
#[derive(Clone)]
pub struct FastaRecord {
  pub fasta_id: FastaId,
  pub seq: Seq,
}
pub type FastaRecords = Vec<FastaRecord>;
pub type SeqPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub struct LogStaPpfMatSets {
  pub log_ppf_mats_on_sas: LogStaPpfMats,
  pub log_ppf_mats_on_2ls: LogStaPpfMats,
  pub log_ppf_mats_on_mls: LogStaPpfMats,
  pub log_ppf_mats_4_first_bpas_on_mls: LogStaPpfMats,
}

pub const INVERSE_TEMPERATURE: FreeEnergy = 1. / (2. * GAS_CONST * TEMPERATURE);

impl LogStapmt {
  pub fn new(seq_len_pair: &(usize, usize)) -> LogStapmt {
    let log_prob_seq_pair = (vec![NEG_INFINITY; seq_len_pair.0], vec![NEG_INFINITY; seq_len_pair.1]);
    let log_prob_mat = SparseLogProbMat::default();
    LogStapmt {
      lbap_mat: log_prob_mat.clone(),
      logp_mat_1: log_prob_seq_pair.0.clone(),
      logp_mat_2: log_prob_seq_pair.1.clone(),
      legp_mat_1: log_prob_seq_pair.0,
      legp_mat_2: log_prob_seq_pair.1,
      lbpap_mat: LogProb4dMat::default(),
    }
  }
  pub fn origin() -> LogStapmt {
    let log_probs = LogProbs::new();
    let log_prob_mat = SparseLogProbMat::default();
    LogStapmt {
      lbap_mat: log_prob_mat.clone(),
      logp_mat_1: log_probs.clone(),
      logp_mat_2: log_probs.clone(),
      legp_mat_1: log_probs.clone(),
      legp_mat_2: log_probs,
      lbpap_mat: LogProb4dMat::default(),
    }
  }
}

impl StaFeParams {
  pub fn origin() -> StaFeParams {
    let lstapmt = LogStapmt::origin();
    StaFeParams {
      lstapmt: lstapmt.clone(),
      lstapmt_on_random_assump: lstapmt,
    }
  }
  pub fn new(rna_id_pair: &RnaIdPair, fasta_records: &FastaRecords, max_gap_num: usize, lbpp_mats: &LogProbMats) -> StaFeParams {
    let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
    let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
    let lstapmt = LogStapmt::new(&seq_len_pair);
    let mut sta_fe_params = StaFeParams {
      lstapmt: lstapmt.clone(),
      lstapmt_on_random_assump: lstapmt,
    };
    let lbpp_mat_pair = (&lbpp_mats[rna_id_pair.0], &lbpp_mats[rna_id_pair.1]);
    let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
    sta_fe_params.lstapmt.lbpap_mat.insert(pseudo_pos_quadruple, 0.);
    sta_fe_params.lstapmt_on_random_assump.lbpap_mat.insert(pseudo_pos_quadruple, 0.);
    for i in 1 .. seq_len_pair.0 - 1 {
      let base = seq_pair.0[i];
      let lnbpp = STEM_PARAMS.lnbpps_with_bases[&base];
      sta_fe_params.lstapmt.logp_mat_1[i] = STEM_PARAMS.logps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.logp_mat_1[i] = lnbpp;
      sta_fe_params.lstapmt.legp_mat_1[i] = STEM_PARAMS.legps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.legp_mat_1[i] = lnbpp;
      for j in 1 .. seq_len_pair.1 - 1 {
        let pos_pair = (i, j);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
        let base_pair = (base, seq_pair.1[j]);
        sta_fe_params.lstapmt.lbap_mat.insert(pos_pair, STEM_PARAMS.lbaps_with_base_pairs[&base_pair]);
        sta_fe_params.lstapmt_on_random_assump.lbap_mat.insert(pos_pair, lnbpp + STEM_PARAMS.lnbpps_with_bases[&base_pair.1]);
      }
      for j in i + 1 .. seq_len_pair.0 - 1 {
        let pos_pair = (i, j);
        let base_pair = (base, seq_pair.0[j]);
        if !lbpp_mat_pair.0.contains_key(&pos_pair) {continue;}
        let lbpp = STEM_PARAMS.lbpps_with_base_pairs[&base_pair];
        for k in 1 .. seq_len_pair.1 - 1 {
          if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num) {continue;}
          for l in k + 1 .. seq_len_pair.1 - 1 {
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num) {continue;}
            let pos_pair_2 = (k, l);
            let pos_quadruple = (pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1);
            if is_min_gap_ok_2(&pos_quadruple, max_gap_num) && lbpp_mat_pair.0.contains_key(&pos_pair) && lbpp_mat_pair.1.contains_key(&pos_pair_2) {
              let base_quadruple = (base_pair.0, base_pair.1, seq_pair.1[k], seq_pair.1[l]);
              let base_pair_2 = (base_quadruple.2, base_quadruple.3);
              sta_fe_params.lstapmt.lbpap_mat.insert(pos_quadruple, STEM_PARAMS.lbpaps_with_base_quadruples[&base_quadruple]);
              sta_fe_params.lstapmt_on_random_assump.lbpap_mat.insert(pos_quadruple, lbpp + STEM_PARAMS.lbpps_with_base_pairs[&base_pair_2]);
            }
          }
        }
      }
    }
    for i in 1 .. seq_len_pair.1 - 1 {
      let base = seq_pair.1[i];
      sta_fe_params.lstapmt.logp_mat_2[i] = STEM_PARAMS.logps_with_bases[&base];
      let lnbpp = STEM_PARAMS.lnbpps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.logp_mat_2[i] = lnbpp;
      sta_fe_params.lstapmt.legp_mat_2[i] = STEM_PARAMS.legps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.legp_mat_2[i] = lnbpp;
    }
    sta_fe_params
  }
}

impl LogStaPpfMats {
  pub fn new() -> LogStaPpfMats {
    let log_ppf_mat = SparseLogPpfMat::default();
    LogStaPpfMats {
      log_ppf_mat: log_ppf_mat.clone(),
      log_ppf_mat_4_bas: log_ppf_mat.clone(),
      log_ppf_mat_4_ogs_1: log_ppf_mat.clone(),
      log_ppf_mat_4_ogs_2: log_ppf_mat.clone(),
      log_ppf_mat_4_egs_1: log_ppf_mat.clone(),
      log_ppf_mat_4_egs_2: log_ppf_mat.clone(),
      log_ppf_mat_4_gaps_1: log_ppf_mat.clone(),
      log_ppf_mat_4_gaps_2: log_ppf_mat,
    }
  }
}

impl FastaRecord {
  pub fn new(fasta_id: FastaId, seq: Seq) -> FastaRecord {
    FastaRecord {
      fasta_id: fasta_id,
      seq: seq,
    }
  }
}

impl LogStaPpfMatSets {
  pub fn new() -> LogStaPpfMatSets {
    let log_ppf_mats = LogStaPpfMats::new();
    LogStaPpfMatSets {
      log_ppf_mats_on_sas: log_ppf_mats.clone(),
      log_ppf_mats_on_2ls: log_ppf_mats.clone(),
      log_ppf_mats_on_mls: log_ppf_mats.clone(),
      log_ppf_mats_4_first_bpas_on_mls: log_ppf_mats,
    }
  }
}

#[inline]
pub fn io_algo_4_rna_stapmt(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb, max_gap_num: usize) -> Stapmt {
  let lstapmt = io_algo_4_rna_lstapmt(seq_pair, seq_len_pair, sta_fe_params, sta_fe_scale_param, max_gap_num);
  get_stapmt(&lstapmt)
}

#[inline]
pub fn io_algo_4_rna_lstapmt(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb, max_gap_num: usize) -> LogStapmt {
  let log_sta_ppf_mats = get_log_sta_ppf_4d_mats(seq_pair, seq_len_pair, sta_fe_params, sta_fe_scale_param, max_gap_num);
  get_lstapmt(seq_pair, seq_len_pair, sta_fe_params, &log_sta_ppf_mats, sta_fe_scale_param, max_gap_num)
}

#[inline]
fn get_stapmt(lstapmt: &LogStapmt) -> Stapmt {
  Stapmt {
    base_align_prob_mat: lstapmt.lbap_mat.iter().map(|(pos_pair, &lbap)| (*pos_pair, lbap.exp())).collect(),
    opening_gap_prob_mat_1: lstapmt.logp_mat_1.iter().map(|&logp| logp.exp()).collect(),
    ogp_mat_2: lstapmt.logp_mat_2.iter().map(|&logp| logp.exp()).collect(),
    extending_gap_prob_mat_1: lstapmt.legp_mat_1.iter().map(|&legp| legp.exp()).collect(),
    egp_mat_2: lstapmt.legp_mat_2.iter().map(|&legp| legp.exp()).collect(),
    base_pair_align_prob_mat: lstapmt.lbpap_mat.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect(),
  }
}

#[inline]
pub fn get_log_sta_ppf_4d_mats(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb, max_gap_num: usize) -> LogPpf4dMat {
  let mut log_sta_ppf_mat = LogPpf4dMat::default();
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  for substr_len_1 in 2 .. seq_len_pair.0 + 1 {
    for substr_len_2 in 2 .. seq_len_pair.1 + 1 {
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_gap_num {continue;}
      for i in 0 .. seq_len_pair.0 + 1 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 0 .. seq_len_pair.1 + 1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple) {continue;}
          let log_sta_ppf_mat_sets = get_log_sta_forward_ppf_mat_sets(&pos_quadruple, &pseudo_pos_quadruple, seq_pair, seq_len_pair, sta_fe_params, &log_sta_ppf_mat, sta_fe_scale_param, max_gap_num);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          let bpap_lor = get_bpa_lor(&pos_quadruple, sta_fe_params, sta_fe_scale_param);
          let ep_of_term_4_log_pf = bpap_lor - INVERSE_TEMPERATURE * (if pos_quadruple == pseudo_pos_quadruple {0.} else {get_hl_fe(seq_pair.0, &(i, j)) + get_hl_fe(seq_pair.1, &(k, l))}) + log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(j - 1, l - 1)];
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
          let ep_of_term_4_log_pf = bpap_lor + log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat[&(j - 1, l - 1)];
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
          let bp_closing_loop_pair = (
            (seq_pair.0[i], seq_pair.0[j]),
            (seq_pair.1[k], seq_pair.1[l]),
          );
          let ml_tm_delta_fe_pair = if pos_quadruple == pseudo_pos_quadruple {
            (0., 0.)
          } else {
            (
              - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[j - 1])))],
              - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[l - 1])))],
            )
          };
          let au_or_gu_end_penalty_delta_fe_pair = (
            - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
            - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
          );
          let ep_of_term_4_log_pf = bpap_lor + if pos_quadruple == pseudo_pos_quadruple {0.} else {- 2. * INVERSE_TEMPERATURE * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1} + log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(j - 1, l - 1)];
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
          log_sta_ppf_mat.insert(pos_quadruple, if max_ep_of_term_4_log_pf.is_finite() {logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)} else {NEG_INFINITY});
        }
      }
    }
  }
  log_sta_ppf_mat
}

#[inline]
fn get_log_sta_forward_ppf_mat_sets(pos_quadruple: &PosQuadruple, pseudo_pos_quadruple: &PosQuadruple, seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, log_sta_ppf_mat_4_bpas: &LogPpf4dMat, sta_fe_scale_param: LogProb, max_gap_num: usize) -> LogStaPpfMatSets {
  let &(i, j, k, l) = pos_quadruple;
  let mut log_sta_ppf_mat_sets = LogStaPpfMatSets::new();
  for n in i .. j {
    for p in k .. l {
      let pos_pair_1 = (n, p);
      if !(is_min_gap_ok_1(&pos_pair_1, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_1, pos_quadruple, max_gap_num)) {continue;}
      if n == i && p == k {
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
      } else {
        if n == i || p == k {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (n - 1, p - 1);
          let ba_lor = get_ba_log_odds_ratio(&pos_pair_1, sta_fe_params, sta_fe_scale_param);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair_2] + ba_lor);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat[&pos_pair_2] + ba_lor;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (m - 1, o - 1);
              if pos_quadruple != pseudo_pos_quadruple && (j - n - 1 + m - i - 1 > MAX_2_LOOP_LEN || l - p - 1 + o - k - 1 > MAX_2_LOOP_LEN) {continue;}
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple_2)) {continue;}
              let two_loop_fe_pair = if pos_quadruple == pseudo_pos_quadruple {
                (0., 0.)
              } else {
                (
                  get_2_loop_fe(seq_pair.0, &(i, j), &(m, n)),
                  get_2_loop_fe(seq_pair.1, &(k, l), &(o, p)),
                )
              };
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (two_loop_fe_pair.0 + two_loop_fe_pair.1) + log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)} else {NEG_INFINITY});
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair_2] + ba_lor;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (m - 1, o - 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple_2)) {continue;}
              let log_sta_pf_4_bpa = log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
              let accessible_bp_pair = (
                (seq_pair.0[m], seq_pair.0[n]),
                (seq_pair.1[o], seq_pair.1[p]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[m - 1], seq_pair.0[n + 1]),
                (seq_pair.1[o - 1], seq_pair.1[p + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                if m > 1 && n < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if m > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if n < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                if o > 1 && p < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if o > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if p < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1) + log_sta_pf_4_bpa;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1) + log_sta_pf_4_bpa;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)} else {NEG_INFINITY});
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair_2] + ba_lor;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (m - 1, o - 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple_2)) {continue;}
              let accessible_bp_pair = (
                (seq_pair.0[m], seq_pair.0[n]),
                (seq_pair.1[o], seq_pair.1[p]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[m - 1], seq_pair.0[n + 1]),
                (seq_pair.1[o - 1], seq_pair.1[p + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                if m > 1 && n < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if m > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if n < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                if o > 1 && p < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if o > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if p < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1) + log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)} else {NEG_INFINITY});
        }
        if n == i {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (n - 1, p);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_1(n, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          }
          let eg_lor = get_eg_lor_1(n, sta_fe_params, sta_fe_scale_param);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_1(n, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1]);
          }
        }
        if p == k {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (n, p - 1);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_2(p, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          }
          let eg_lor = get_eg_lor_2(p, sta_fe_params, sta_fe_scale_param);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor 
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_2(p, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1]);
          }
        }
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
      }
    }
  }
  log_sta_ppf_mat_sets
}

#[inline]
fn get_ba_log_odds_ratio(pos_pair: &PosPair, sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb) -> StaFreeEnergy {
  INVERSE_TEMPERATURE * sta_fe_scale_param * (sta_fe_params.lstapmt.lbap_mat[pos_pair] - sta_fe_params.lstapmt_on_random_assump.lbap_mat[pos_pair])
}

#[inline]
fn get_og_lor_1(pos: Pos, sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb) -> StaFreeEnergy {
  INVERSE_TEMPERATURE * sta_fe_scale_param * (sta_fe_params.lstapmt.logp_mat_1[pos] - sta_fe_params.lstapmt_on_random_assump.logp_mat_1[pos])
}

#[inline]
fn get_og_lor_2(pos: Pos, sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb) -> StaFreeEnergy {
  INVERSE_TEMPERATURE * sta_fe_scale_param * (sta_fe_params.lstapmt.logp_mat_2[pos] - sta_fe_params.lstapmt_on_random_assump.logp_mat_2[pos])
}

#[inline]
fn get_eg_lor_1(pos: Pos, sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb) -> StaFreeEnergy {
  INVERSE_TEMPERATURE * sta_fe_scale_param * (sta_fe_params.lstapmt.legp_mat_1[pos] - sta_fe_params.lstapmt_on_random_assump.legp_mat_1[pos])
}

#[inline]
fn get_eg_lor_2(pos: Pos, sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb) -> StaFreeEnergy {
  INVERSE_TEMPERATURE * sta_fe_scale_param * (sta_fe_params.lstapmt.legp_mat_2[pos] - sta_fe_params.lstapmt_on_random_assump.legp_mat_2[pos])
}

#[inline]
fn get_bpa_lor(pos_quadruple: &PosQuadruple, sta_fe_params: &StaFeParams, sta_fe_scale_param: LogProb) -> StaFreeEnergy {
  INVERSE_TEMPERATURE * sta_fe_scale_param * (sta_fe_params.lstapmt.lbpap_mat[pos_quadruple]
  - sta_fe_params.lstapmt_on_random_assump.lbpap_mat[pos_quadruple])
}

#[inline]
fn get_lstapmt(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, log_sta_ppf_mat_4_bpas: &LogPpf4dMat, sta_fe_scale_param: LogProb, max_gap_num: usize) -> LogStapmt {
  let mut lstapmt = LogStapmt::new(seq_len_pair);
  let mut seqs_of_eps_of_terms_4_lbaps_with_pos_pairs = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_lbaps_with_pos_pairs = EpsOfTerms4LogProbsWithPosPairs::default();
  let mut seqs_of_eps_of_terms_4_logps_1 = vec![EpsOfTerms4LogProb::new(); seq_len_pair.0];
  let mut max_eps_of_terms_4_logps_1 = vec![NEG_INFINITY; seq_len_pair.0];
  let mut seqs_of_eps_of_terms_4_logps_2 = vec![EpsOfTerms4LogProb::new(); seq_len_pair.1];
  let mut max_eps_of_terms_4_logps_2 = vec![NEG_INFINITY; seq_len_pair.1];
  let mut seqs_of_eps_of_terms_4_legps_1 = seqs_of_eps_of_terms_4_logps_1.clone();
  let mut max_eps_of_terms_4_legps_1 = max_eps_of_terms_4_logps_1.clone();
  let mut seqs_of_eps_of_terms_4_legps_2 = seqs_of_eps_of_terms_4_logps_2.clone();
  let mut max_eps_of_terms_4_legps_2 = max_eps_of_terms_4_logps_2.clone();
  let mut seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples = SeqsOfEpsOfTerms4LogProbsWithPosQuadruples::default();
  let mut max_eps_of_terms_4_lbpaps_with_pos_quadruples = EpsOfTerms4LogProbsWithPosQuadruples::default();
  for pos_pair in sta_fe_params.lstapmt.lbap_mat.keys() {
    seqs_of_eps_of_terms_4_lbaps_with_pos_pairs.insert(*pos_pair, EpsOfTerms4LogProb::new());
    max_eps_of_terms_4_lbaps_with_pos_pairs.insert(*pos_pair, NEG_INFINITY);
  }
  for pos_quadruple in sta_fe_params.lstapmt.lbpap_mat.keys() {
    seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples.insert(*pos_quadruple, EpsOfTerms4LogProb::new());
    max_eps_of_terms_4_lbpaps_with_pos_quadruples.insert(*pos_quadruple, NEG_INFINITY);
  }
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  for substr_len_1 in (3 .. seq_len_pair.0 + 1).rev() {
    for substr_len_2 in (3 .. seq_len_pair.1 + 1).rev() {
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_gap_num {continue;}
      for i in 0 .. seq_len_pair.0 + 1 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 0 .. seq_len_pair.1 + 1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple) {continue;}
          if pos_quadruple == pseudo_pos_quadruple {
            lstapmt.lbpap_mat.insert(pos_quadruple, 0.);
          } else {
            let max_ep_of_term_4_log_prob = max_eps_of_terms_4_lbpaps_with_pos_quadruples[&pos_quadruple];
            if max_ep_of_term_4_log_prob.is_finite() {
              let lbpap = logsumexp(&seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples[&pos_quadruple][..], max_ep_of_term_4_log_prob);
              lstapmt.lbpap_mat.insert(pos_quadruple, lbpap);
              if lbpap.is_finite() {
                let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_1[i];
                if lbpap > max_ep_of_term_4_log_prob {
                  max_eps_of_terms_4_logps_1[i] = lbpap;
                }
                let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_1[j];
                if lbpap > max_ep_of_term_4_log_prob {
                  max_eps_of_terms_4_logps_1[j] = lbpap;
                }
                let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_2[k];
                if lbpap > max_ep_of_term_4_log_prob {
                  max_eps_of_terms_4_logps_2[k] = lbpap;
                }
                let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_2[l];
                if lbpap > max_ep_of_term_4_log_prob {
                  max_eps_of_terms_4_logps_2[l] = lbpap;
                }
              }
            } else {
              lstapmt.lbpap_mat.insert(pos_quadruple, NEG_INFINITY);
            }
          }
          let log_coefficient = lstapmt.lbpap_mat[&pos_quadruple] + get_bpa_lor(&pos_quadruple, sta_fe_params, sta_fe_scale_param) - log_sta_ppf_mat_4_bpas[&pos_quadruple];
          if !log_coefficient.is_finite() {continue;}
          let bp_closing_loop_pair = (
            (seq_pair.0[i], seq_pair.0[j]),
            (seq_pair.1[k], seq_pair.1[l]),
          );
          let ml_tm_delta_fe_pair = if pos_quadruple == pseudo_pos_quadruple {
            (0., 0.)
          } else {
            (
              - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[j - 1])))],
              - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[l - 1])))],
            )
          };
          let au_or_gu_end_penalty_delta_fe_pair = (
            - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
            - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
          );
          let log_sta_forward_ppf_mat_sets = get_log_sta_forward_ppf_mat_sets(&pos_quadruple, &pseudo_pos_quadruple, seq_pair, seq_len_pair, sta_fe_params, &log_sta_ppf_mat_4_bpas, sta_fe_scale_param, max_gap_num);
          let log_sta_backward_ppf_mat_sets = get_log_sta_backward_ppf_mat_sets(&pos_quadruple, &pseudo_pos_quadruple, seq_pair, seq_len_pair, sta_fe_params, &log_sta_ppf_mat_4_bpas, sta_fe_scale_param, max_gap_num);
          let hl_fe_pair = if pos_quadruple == pseudo_pos_quadruple {
            (0., 0.)
          } else {
            (
              - INVERSE_TEMPERATURE * get_hl_fe(seq_pair.0, &(i, j)),
              - INVERSE_TEMPERATURE * get_hl_fe(seq_pair.1, &(k, l)),
            )
          };
          for n in i + 1 .. j {
            for p in k + 1 .. l {
              let pos_pair = (n, p);
              if !(is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num)) {continue;}
              let ba_lor = get_ba_log_odds_ratio(&pos_pair, sta_fe_params, sta_fe_scale_param);
              let log_coefficient_1 = ba_lor + log_coefficient;
              let log_coefficient_2 = log_coefficient_1 - 2. * INVERSE_TEMPERATURE * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1;
              let eps_of_terms_4_log_prob = [
                log_coefficient_1 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n - 1, p - 1)] + hl_fe_pair.0 + hl_fe_pair.1 + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_1 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_1 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(n - 1, p - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(n + 1, p + 1)],
              ];
              for &ep_of_term_4_log_prob in eps_of_terms_4_log_prob.iter() {
                if ep_of_term_4_log_prob.is_finite() {
                  seqs_of_eps_of_terms_4_lbaps_with_pos_pairs.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                  let max_ep_of_term_4_log_prob = max_eps_of_terms_4_lbaps_with_pos_pairs.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                  if ep_of_term_4_log_prob > *max_ep_of_term_4_log_prob {
                    *max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
                }
              }
              let log_coefficient_2 = log_coefficient - 2. * INVERSE_TEMPERATURE * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1;
              let eps_of_terms_4_log_prob = [
                log_coefficient + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1[&pos_pair] + hl_fe_pair.0 + hl_fe_pair.1 + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)],
              ];
              for &ep_of_term_4_log_prob in eps_of_terms_4_log_prob.iter() {
                if ep_of_term_4_log_prob.is_finite() {
                  seqs_of_eps_of_terms_4_legps_1[n].push(ep_of_term_4_log_prob);
                  let ref mut max_ep_of_term_4_log_prob = max_eps_of_terms_4_legps_1[n];
                  if ep_of_term_4_log_prob > *max_ep_of_term_4_log_prob {
                    *max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
                }
              }
              let eps_of_terms_4_log_prob = [
                log_coefficient + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2[&pos_pair] + hl_fe_pair.0 + hl_fe_pair.1 + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
                log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2[&pos_pair] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)],
              ];
              for &ep_of_term_4_log_prob in eps_of_terms_4_log_prob.iter() {
                if ep_of_term_4_log_prob.is_finite() {
                  seqs_of_eps_of_terms_4_legps_2[p].push(ep_of_term_4_log_prob);
                  let ref mut max_ep_of_term_4_log_prob = max_eps_of_terms_4_legps_2[p];
                  if ep_of_term_4_log_prob > *max_ep_of_term_4_log_prob {
                    *max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
                }
              }
            }
          }
          for m in i + 1 .. j {
            for o in k + 1 .. l {
              let pos_pair = (m, o);
              if !(is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num)) {continue;}
              for n in m + 1 .. j {
                for p in o + 1 .. l {
                  let pos_quadruple_2 = (m, n, o, p);
                  let pos_pair = (n, p);
                  if !(is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple_2)) {continue;}
                  let log_sta_pf_4_bpa = log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
                  let accessible_bp_pair = (
                    (seq_pair.0[m], seq_pair.0[n]),
                    (seq_pair.1[o], seq_pair.1[p]),
                  );
                  let stacking_bp_pair = (
                    (seq_pair.0[m - 1], seq_pair.0[n + 1]),
                    (seq_pair.1[o - 1], seq_pair.1[p + 1]),
                  );
                  let ml_tm_or_de_delta_fe_pair = (
                    - INVERSE_TEMPERATURE * if m > 1 && n < seq_len_pair.0 - 2 {
                      ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                    } else if m > 1 {
                      FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                    } else if n < seq_len_pair.0 - 2 {
                      THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                    } else {
                      0.
                    },
                    - INVERSE_TEMPERATURE * if o > 1 && p < seq_len_pair.1 - 2 {
                      ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                    } else if o > 1 {
                      FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                    } else if p < seq_len_pair.1 - 2 {
                      THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                    } else {
                      0.
                    },
                  );
                  let au_or_gu_end_penalty_delta_fe_pair_2 = (
                    - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                    - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                  );
                  let log_coefficient_1 = log_sta_pf_4_bpa + log_coefficient;
                  let log_coefficient_2 = log_coefficient_1 - 2. * INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair_2.0 + au_or_gu_end_penalty_delta_fe_pair_2.1;
                  let two_loop_fe_pair = if pos_quadruple != pseudo_pos_quadruple && (j - n - 1 + m - i - 1 > MAX_2_LOOP_LEN || l - p - 1 + o - k - 1 > MAX_2_LOOP_LEN) {
                    (NEG_INFINITY, NEG_INFINITY)
                  } else {
                    if pos_quadruple == pseudo_pos_quadruple {
                      (0., 0.)
                    } else {
                      (
                        - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.0, &(i, j), &(m, n)),
                        - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.1, &(k, l), &(o, p)),
                      )
                    }
                  };
                  let eps_of_terms_4_log_prob = [
                    log_coefficient_1 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(m - 1, o - 1)] + two_loop_fe_pair.0 + two_loop_fe_pair.1 + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                    log_coefficient_2 + log_sta_forward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(m - 1, o - 1)] + log_sta_backward_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(n + 1, p + 1)],
                  ];
                  for &ep_of_term_4_log_prob in eps_of_terms_4_log_prob.iter() {
                    if ep_of_term_4_log_prob.is_finite() {
                      seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples.get_mut(&pos_quadruple_2).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                      let max_ep_of_term_4_log_prob = max_eps_of_terms_4_lbpaps_with_pos_quadruples.get_mut(&pos_quadruple_2).expect("Failed to get an element of a hash map.");
                      if ep_of_term_4_log_prob > *max_ep_of_term_4_log_prob {
                        *max_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  for pos_pair in sta_fe_params.lstapmt.lbap_mat.keys() {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_lbaps_with_pos_pairs[pos_pair];
    if max_ep_of_term_4_log_prob.is_finite() {
      let lbap = logsumexp(&seqs_of_eps_of_terms_4_lbaps_with_pos_pairs[pos_pair][..], max_ep_of_term_4_log_prob);
      if lbap.is_finite() {
        lstapmt.lbap_mat.insert(*pos_pair, lbap);
        let &(i, j) = pos_pair;
        seqs_of_eps_of_terms_4_logps_1[i].push(lbap);
        let ref mut max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_1[i];
        if lbap > *max_ep_of_term_4_log_prob {
          *max_ep_of_term_4_log_prob = lbap;
        }
        seqs_of_eps_of_terms_4_logps_2[j].push(lbap);
        let ref mut max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_2[j];
        if lbap > *max_ep_of_term_4_log_prob {
          *max_ep_of_term_4_log_prob = lbap;
        }
      }
    }
  }
  for i in 0 .. seq_len_pair.0 {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_legps_1[i];
    let legp = if max_ep_of_term_4_log_prob.is_finite() {
      logsumexp(&seqs_of_eps_of_terms_4_legps_1[i][..], max_ep_of_term_4_log_prob)
    } else {
      NEG_INFINITY
    };
    if legp.is_finite() {
      lstapmt.legp_mat_1[i] = legp;
      seqs_of_eps_of_terms_4_logps_1[i].push(legp);
      let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_1[i];
      if legp > max_ep_of_term_4_log_prob {
        max_eps_of_terms_4_logps_1[i]= legp;
      }
    }
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_1[i];
    lstapmt.logp_mat_1[i] = if max_ep_of_term_4_log_prob.is_finite() {(1. - logsumexp(&seqs_of_eps_of_terms_4_logps_1[i][..], max_ep_of_term_4_log_prob).exp()).ln()} else {NEG_INFINITY};
  }
  for i in 0 .. seq_len_pair.1 {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_legps_2[i];
    let legp = if max_ep_of_term_4_log_prob.is_finite() {
      logsumexp(&seqs_of_eps_of_terms_4_legps_2[i][..], max_ep_of_term_4_log_prob)
    } else {
      NEG_INFINITY
    };
    if legp.is_finite() {
      lstapmt.legp_mat_2[i] = legp;
      seqs_of_eps_of_terms_4_logps_2[i].push(legp);
      let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_2[i];
      if legp > max_ep_of_term_4_log_prob {
        max_eps_of_terms_4_logps_2[i] = legp;
      }
    }
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_logps_2[i];
    lstapmt.logp_mat_2[i] = if max_ep_of_term_4_log_prob.is_finite() {(1. - logsumexp(&seqs_of_eps_of_terms_4_logps_2[i][..], max_ep_of_term_4_log_prob).exp()).ln()} else {NEG_INFINITY};
  }
  lstapmt
}

#[inline]
fn get_log_sta_backward_ppf_mat_sets(pos_quadruple: &PosQuadruple, pseudo_pos_quadruple: &PosQuadruple, seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, log_sta_ppf_mat_4_bpas: &LogPpf4dMat, sta_fe_scale_param: LogProb, max_gap_num: usize) -> LogStaPpfMatSets {
  let &(i, j, k, l) = pos_quadruple;
  let mut log_sta_ppf_mat_sets = LogStaPpfMatSets::new();
  for m in (i + 1 .. j + 1).rev() {
    for o in (k + 1 .. l + 1).rev() {
      let pos_pair_1 = (m, o);
      if !(is_min_gap_ok_1(&pos_pair_1, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_1, pos_quadruple, max_gap_num)) {continue;}
      if m == j && o == l {
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
      } else {
        if m == j || o == l {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (m + 1, o + 1);
          let ba_lor = get_ba_log_odds_ratio(&pos_pair_1, sta_fe_params, sta_fe_scale_param);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair_2] + ba_lor);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat[&pos_pair_2] + ba_lor;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (n + 1, p + 1);
              if pos_quadruple != pseudo_pos_quadruple && (j - n - 1 + m - i - 1 > MAX_2_LOOP_LEN || l - p - 1 + o - k - 1 > MAX_2_LOOP_LEN) {continue;}
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple_2)) {continue;}
              let two_loop_fe_pair = if pos_quadruple == pseudo_pos_quadruple {
                (0., 0.)
              } else {
                (
                  get_2_loop_fe(seq_pair.0, &(i, j), &(m, n)),
                  get_2_loop_fe(seq_pair.1, &(k, l), &(o, p)),
                )
              };
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (two_loop_fe_pair.0 + two_loop_fe_pair.1) + log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)} else {NEG_INFINITY});
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair_2] + ba_lor;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (n + 1, p + 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple_2)) {continue;}
              let log_sta_pf_4_bpa = log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
              let accessible_bp_pair = (
                (seq_pair.0[m], seq_pair.0[n]),
                (seq_pair.1[o], seq_pair.1[p]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[m - 1], seq_pair.0[n + 1]),
                (seq_pair.1[o - 1], seq_pair.1[p + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                if m > 1 && n < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if m > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if n < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                if o > 1 && p < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if o > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if p < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1) + log_sta_pf_4_bpa;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1) + log_sta_pf_4_bpa;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)} else {NEG_INFINITY});
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair_2] + ba_lor;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (n + 1, p + 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.lstapmt.lbpap_mat.contains_key(&pos_quadruple_2)) {continue;}
              let accessible_bp_pair = (
                (seq_pair.0[m], seq_pair.0[n]),
                (seq_pair.1[o], seq_pair.1[p]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[m - 1], seq_pair.0[n + 1]),
                (seq_pair.1[o - 1], seq_pair.1[p + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                if m > 1 && n < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if m > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if n < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                if o > 1 && p < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if o > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if p < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1) + log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)} else {NEG_INFINITY});
        }
        if m == j {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (m + 1, o);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_1(m, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
          }
          let eg_lor = get_eg_lor_1(m, sta_fe_params, sta_fe_scale_param);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_1(m, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_1[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_1[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_1[&pos_pair_1]);
          }
        }
        if o == l {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (m, o + 1);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_2(o, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
          }
          let eg_lor = get_eg_lor_2(o, sta_fe_params, sta_fe_scale_param);
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor 
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2.insert(pos_pair_1, if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + eg_lor
          } else {
            NEG_INFINITY
          });
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let og_lor = get_og_lor_2(o, sta_fe_params, sta_fe_scale_param);
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
              logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
            } else {
              NEG_INFINITY
            });
          } else {
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_egs_2[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_egs_2[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1]);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_egs_2[&pos_pair_1]);
          }
        }
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
      }
    }
  }
  log_sta_ppf_mat_sets
}
 
#[inline]
pub fn prob_cons_transformation_of_lstapmt(lstapmts_with_rna_id_pairs: &LstapmtsWithRnaIdPairs, rna_id_pair: &RnaIdPair, num_of_rnas: usize) -> LogStapmt {
  let mut lstapmt = lstapmts_with_rna_id_pairs[rna_id_pair].clone();
  let log_coefficient = -((num_of_rnas - 1) as Prob).ln();
  for (pos_quadruple, lbpap) in lstapmt.lbpap_mat.iter_mut() {
    let (i, j, k, l) = *pos_quadruple;
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = *lbpap;
    if max_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      for (&(m, n, o, p), &lbpap_1) in lstapmts_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].lbpap_mat.iter() {
        if (rna_id_pair.0 < rna_id && m != i && n != j) || (rna_id_pair.0 > rna_id && o != i && p != j) {continue;}
        for (&(m, n, o, p), &lbpap_2) in lstapmts_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].lbpap_mat.iter() {
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
      *lbpap = log_coefficient + logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
    }
  }
  for (&(i, j), lbap) in lstapmt.lbap_mat.iter_mut() {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = *lbap;
    if max_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      for (&(k, l), lbap_1) in lstapmts_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].lbap_mat.iter() {
        if (rna_id_pair.0 < rna_id && k != i) || (rna_id_pair.0 > rna_id && l != i) {continue;}
        for (&(k, l), lbap_2) in lstapmts_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].lbap_mat.iter() {
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
      *lbap = log_coefficient + logsumexp(&eps_of_terms_4_log_prob[..], max_ep_of_term_4_log_prob);
    }
  }
  lstapmt
}

#[inline]
pub fn pct_of_lbpp_mat(lstapmts_with_rna_id_pairs: &LstapmtsWithRnaIdPairs, rna_id: RnaId, num_of_rnas: usize, lbpp_mat: &SparseLogProbMat) -> SparseLogProbMat {
  let log_coefficient = -((num_of_rnas - 1) as Prob).ln();
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs = EpsOfTerms4LogProbsWithPosPairs::default();
  for pos_pair in lbpp_mat.keys() {
    if pos_pair.0 == 0 {continue;}
    seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, EpsOfTerms4LogProb::new());
    max_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, NEG_INFINITY);
  }
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id == rna_id_2 {continue;}
    let rna_id_pair = if rna_id < rna_id_2 {(rna_id, rna_id_2)} else {(rna_id_2, rna_id)};
    let ref ref_2_lstapmt = lstapmts_with_rna_id_pairs[&rna_id_pair];
    for (pos_quadruple, &lbpap) in ref_2_lstapmt.lbpap_mat.iter() {
      if pos_quadruple.0 == 0 {continue;}
      if lbpap.is_finite() {
        let pos_pair = if rna_id < rna_id_2 {(pos_quadruple.0, pos_quadruple.1)} else {(pos_quadruple.2, pos_quadruple.3)};
        let eps_of_terms_4_log_prob = seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
        eps_of_terms_4_log_prob.push(lbpap);
        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
        if lbpap > *max_ep_of_term_4_log_prob {
          *max_ep_of_term_4_log_prob = lbpap;
        }
      }
    }
  }
  let mut lbpp_mat = lbpp_mat.clone();
  for (pos_pair, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.iter() {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs[pos_pair];
    if max_ep_of_term_4_log_prob.is_finite() {
      let new_lbpp = log_coefficient + logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob);
      let lbpp = lbpp_mat.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
      if new_lbpp > *lbpp {
        *lbpp = new_lbpp;
      }
    }
  }
  lbpp_mat
}

#[inline]
pub fn get_stapmts_with_rna_id_pairs(lstapmts_with_rna_id_pairs: &LstapmtsWithRnaIdPairs) -> StapmtsWithRnaIdPairs {
  let mut stapmts_with_rna_id_pairs = StapmtsWithRnaIdPairs::default();
  for (rna_id_pair, lstapmt) in lstapmts_with_rna_id_pairs.iter() {
    stapmts_with_rna_id_pairs.insert(*rna_id_pair, get_stapmt(lstapmt));
  }
  stapmts_with_rna_id_pairs
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
pub fn remove_little_lbpps_from_lbpp_mat(lbpp_mat: &LogProbMat, min_lbpp: LogProb) -> SparseLogProbMat {
  let mut sparse_lbpp_mat = SparseLogProbMat::default();
  for (i, lbpps) in lbpp_mat.iter().enumerate() {
    for (j, &lbpp) in lbpps.iter().enumerate() {
      if lbpp >= min_lbpp && lbpp > NEG_INFINITY {
        sparse_lbpp_mat.insert((i + 1, j + 1), lbpp);
      }
    }
  }
  sparse_lbpp_mat
}

#[inline]
pub fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
}

#[inline]
pub fn get_seq_len_diff(seq_len_pair: &(usize, usize)) -> usize {
  max(seq_len_pair.0, seq_len_pair.1) - min(seq_len_pair.0, seq_len_pair.1)
}

#[inline]
fn is_min_gap_ok_1(pos_pair: &PosPair, pos_quadruple: &PosQuadruple, max_gap_num: usize) -> bool {
  let min_gap_num_1 = get_min_gap_num(&(pos_quadruple.0, pos_pair.0, pos_quadruple.2, pos_pair.1));
  let min_gap_num_2 = get_min_gap_num(&(pos_pair.0, pos_quadruple.1, pos_pair.1, pos_quadruple.3));
  if min_gap_num_1 <= max_gap_num && min_gap_num_2 <= max_gap_num {
    true
  } else {
    false
  }
}

#[inline]
fn is_min_gap_ok_2(pos_quadruple: &PosQuadruple, max_gap_num: usize) -> bool {
  let min_gap_num = get_min_gap_num(&pos_quadruple);
  if min_gap_num <= max_gap_num {
    true
  } else {
    false
  }
}

#[inline]
fn get_min_gap_num(pos_quadruple: &PosQuadruple) -> usize {
  let substr_len_pair = (pos_quadruple.1 - pos_quadruple.0, pos_quadruple.3 - pos_quadruple.2);
  get_seq_len_diff(&substr_len_pair)
}
