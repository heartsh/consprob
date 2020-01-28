extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate getopts;

pub mod utils;

pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use utils::*;

pub type PosQuadruple = (Pos, Pos, Pos, Pos);
pub type SparseProbMat = HashMap<PosPair, Prob, Hasher>;
pub type Prob4dMat = HashMap<PosQuadruple, Prob, Hasher>;
pub type SparseLogProbMat = HashMap<PosPair, LogProb, Hasher>;
pub type LogProb4dMat = HashMap<PosQuadruple, LogProb, Hasher>;
type SparseLogPpfMat = HashMap<PosPair, LogPf, Hasher>;
type LogPpf4dMat = HashMap<PosQuadruple, LogPf, Hasher>;
#[derive(Clone)]
pub struct LogStaInsidePpf4dMats {
  pub log_ppf_mat: LogPpf4dMat,
  pub log_forward_ppf_mat_4_bas: LogPpf4dMat,
  pub log_forward_ppf_mat_4_gaps_1: LogPpf4dMat,
  pub log_forward_ppf_mat_4_gaps_2: LogPpf4dMat,
  pub log_backward_ppf_mat_4_bas: LogPpf4dMat,
  pub log_backward_ppf_mat_4_gaps_1: LogPpf4dMat,
  pub log_backward_ppf_mat_4_gaps_2: LogPpf4dMat,
}
#[derive(Clone)]
pub struct LogStaInsidePpfMats {
  pub log_forward_ppf_mat: SparseLogPpfMat,
  pub log_forward_ppf_mat_4_bas: SparseLogPpfMat,
  pub log_forward_ppf_mat_4_gaps_1: SparseLogPpfMat,
  pub log_forward_ppf_mat_4_gaps_2: SparseLogPpfMat,
  pub log_backward_ppf_mat: SparseLogPpfMat,
  pub log_backward_ppf_mat_4_bas: SparseLogPpfMat,
  pub log_backward_ppf_mat_4_gaps_1: SparseLogPpfMat,
  pub log_backward_ppf_mat_4_gaps_2: SparseLogPpfMat,
}
#[derive(Clone)]
pub struct LogStaInsidePpfMatSets {
  pub log_ppf_mat_4_bpas: LogPpf4dMat,
  pub log_ppf_mats_on_sa: LogStaInsidePpf4dMats,
  pub log_ppf_mats_4_internal_multiloop: LogStaInsidePpf4dMats,
  pub log_ppf_mats_4_first_bpas_on_mls: LogStaInsidePpf4dMats,
  pub log_ppf_mats_4_external_loop: LogStaInsidePpfMats,
}
#[derive(Clone)]
pub struct LogStaOutsidePpf4dMats {
  pub log_ppf_mat_4_bpas: LogPpf4dMat,
  pub log_ppf_mat_4_bpas_on_el: LogPpf4dMat,
  pub log_ppf_mat_4_bpas_on_internal_2loops: LogPpf4dMat,
  pub log_ppf_mat_4_bpas_on_internal_mls: LogPpf4dMat,
  pub log_ppf_mat_4_right: LogPpf4dMat,
  pub log_ppf_mat_4_right_2: LogPpf4dMat,
}
pub struct StaFeParams {
  pub ba_score_mat: SparseLogProbMat,
  pub bpa_score_mat: LogProb4dMat,
  pub opening_gap_penalty: StaFreeEnergy,
  pub extending_gap_penalty: StaFreeEnergy,
}
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type StaFeParamSetsWithRnaIdPairs = HashMap<RnaIdPair, StaFeParams, Hasher>;
type Prob4dMatsWithRnaIdPairs = HashMap<RnaIdPair, Prob4dMat, Hasher>;
pub type LogProb4dMatsWithRnaIdPairs = HashMap<RnaIdPair, LogProb4dMat, Hasher>;
pub type LogProbMats = Vec<SparseLogProbMat>;
pub type ProbMats = Vec<SparseProbMat>;
type SeqsOfEpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, ExpPartOfTerm4LogProb, Hasher>;
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

pub const INVERSE_TEMPERATURE: FreeEnergy = 1. / (GAS_CONST * TEMPERATURE);

impl StaFeParams {
  pub fn origin() -> StaFeParams {
    StaFeParams {
      ba_score_mat: SparseLogProbMat::default(),
      bpa_score_mat: LogProb4dMat::default(),
      opening_gap_penalty: 0.,
      extending_gap_penalty: 0.,
    }
  }
  pub fn new(rna_id_pair: &RnaIdPair, fasta_records: &FastaRecords, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_pos_dist_4_el: usize, max_substr_dist: usize, bpp_mats: &ProbMats, sta_fe_scale_param: StaFreeEnergy, opening_gap_penalty: StaFreeEnergy, extending_gap_penalty: StaFreeEnergy) -> StaFeParams {
    let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
    let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
    let mut sta_fe_params = StaFeParams::origin();
    sta_fe_params.opening_gap_penalty = sta_fe_scale_param * opening_gap_penalty;
    sta_fe_params.extending_gap_penalty = sta_fe_scale_param * extending_gap_penalty;
    let bpp_mat_pair = (&bpp_mats[rna_id_pair.0], &bpp_mats[rna_id_pair.1]);
    for i in 1 .. seq_len_pair.0 - 1 {
      let base = seq_pair.0[i];
      for j in 1 .. seq_len_pair.1 - 1 {
        let pos_pair = (i, j);
        if !is_pos_dist_ok(&pos_pair, max_pos_dist_4_el) {continue;}
        let base_pair = (base, seq_pair.1[j]);
        sta_fe_params.ba_score_mat.insert(pos_pair, INVERSE_TEMPERATURE * sta_fe_scale_param * RIBOSUM_85_60_BA_SCORE_MAT[&base_pair]);
      }
      let upper_j = if i + max_bp_span_pair.0 >= seq_len_pair.0 - 1 {seq_len_pair.0 - 1} else {i + max_bp_span_pair.0};
      for j in i + 1 .. upper_j {
        let pos_pair = (i, j);
        let base_pair = (base, seq_pair.0[j]);
        if !bpp_mat_pair.0.contains_key(&pos_pair) {continue;}
        for k in 1 .. seq_len_pair.1 - 1 {
          if !is_pos_dist_ok(&(i, k), max_pos_dist_4_il) {continue;}
          let upper_l = if k + max_bp_span_pair.1 >= seq_len_pair.1 - 1 {seq_len_pair.1 - 1} else {k + max_bp_span_pair.1};
          for l in k + 1 .. upper_l {
            if !is_pos_dist_ok(&(j, l), max_pos_dist_4_il) {continue;}
            let pos_pair_2 = (k, l);
            let pos_quadruple = (pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1);
            if !(is_substr_dist_ok(&pos_quadruple, max_substr_dist, max_pos_dist_4_il) && bpp_mat_pair.1.contains_key(&pos_pair_2)) {continue;}
            let base_quadruple = (base_pair, (seq_pair.1[k], seq_pair.1[l]));
            sta_fe_params.bpa_score_mat.insert(pos_quadruple, INVERSE_TEMPERATURE * sta_fe_scale_param * RIBOSUM_85_60_BPA_SCORE_MAT[&base_quadruple]);
          }
        }
      }
    }
    sta_fe_params
  }
}

impl LogStaInsidePpf4dMats {
  pub fn new() -> LogStaInsidePpf4dMats {
    let log_ppf_4d_mat = LogPpf4dMat::default();
    LogStaInsidePpf4dMats {
      log_ppf_mat: log_ppf_4d_mat.clone(),
      log_forward_ppf_mat_4_bas: log_ppf_4d_mat.clone(),
      log_forward_ppf_mat_4_gaps_1: log_ppf_4d_mat.clone(),
      log_forward_ppf_mat_4_gaps_2: log_ppf_4d_mat.clone(),
      log_backward_ppf_mat_4_bas: log_ppf_4d_mat.clone(),
      log_backward_ppf_mat_4_gaps_1: log_ppf_4d_mat.clone(),
      log_backward_ppf_mat_4_gaps_2: log_ppf_4d_mat,
    }
  }
}

impl LogStaInsidePpfMats {
  pub fn new() -> LogStaInsidePpfMats {
    let log_ppf_mat = SparseLogPpfMat::default();
    LogStaInsidePpfMats {
      log_forward_ppf_mat: log_ppf_mat.clone(),
      log_forward_ppf_mat_4_bas: log_ppf_mat.clone(),
      log_forward_ppf_mat_4_gaps_1: log_ppf_mat.clone(),
      log_forward_ppf_mat_4_gaps_2: log_ppf_mat.clone(),
      log_backward_ppf_mat: log_ppf_mat.clone(),
      log_backward_ppf_mat_4_bas: log_ppf_mat.clone(),
      log_backward_ppf_mat_4_gaps_1: log_ppf_mat.clone(),
      log_backward_ppf_mat_4_gaps_2: log_ppf_mat,
    }
  }
}

impl LogStaInsidePpfMatSets {
  pub fn new() -> LogStaInsidePpfMatSets {
    let log_ppf_4d_mat = LogPpf4dMat::default();
    let log_ppf_4d_mats = LogStaInsidePpf4dMats::new();
    LogStaInsidePpfMatSets {
      log_ppf_mat_4_bpas: log_ppf_4d_mat.clone(),
      log_ppf_mats_on_sa: log_ppf_4d_mats.clone(),
      log_ppf_mats_4_internal_multiloop: log_ppf_4d_mats.clone(),
      log_ppf_mats_4_first_bpas_on_mls: log_ppf_4d_mats,
      log_ppf_mats_4_external_loop: LogStaInsidePpfMats::new(),
    }
  }
}

impl LogStaOutsidePpf4dMats {
  pub fn new() -> LogStaOutsidePpf4dMats {
    let log_ppf_4d_mat = LogPpf4dMat::default();
    LogStaOutsidePpf4dMats {
      log_ppf_mat_4_bpas: log_ppf_4d_mat.clone(),
      log_ppf_mat_4_bpas_on_el: log_ppf_4d_mat.clone(),
      log_ppf_mat_4_bpas_on_internal_2loops: log_ppf_4d_mat.clone(),
      log_ppf_mat_4_bpas_on_internal_mls: log_ppf_4d_mat.clone(),
      log_ppf_mat_4_right: log_ppf_4d_mat.clone(),
      log_ppf_mat_4_right_2: log_ppf_4d_mat,
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

#[inline]
pub fn io_algo_4_base_pair_align_prob_mat(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_pos_dist_4_el: usize, max_substr_dist: usize) -> Prob4dMat {
  let lbpap_mat = io_algo_4_lbpap_mat(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_pos_dist_4_il, max_pos_dist_4_el, max_substr_dist);
  get_bpap_mat(&lbpap_mat)
}

#[inline]
pub fn io_algo_4_lbpap_mat(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_pos_dist_4_el: usize, max_substr_dist: usize) -> LogProb4dMat {
  let log_sta_inside_ppf_mat_sets = get_log_sta_inside_ppf_mat_sets(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_pos_dist_4_il, max_pos_dist_4_el, max_substr_dist);
  let log_sta_outside_ppf_4d_mats = get_log_sta_outside_ppf_4d_mats(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_pos_dist_4_il, max_substr_dist, &log_sta_inside_ppf_mat_sets);
  get_lbpap_mat(seq_len_pair, &log_sta_inside_ppf_mat_sets, &log_sta_outside_ppf_4d_mats)
}

#[inline]
fn get_bpap_mat(lbpap_mat: &LogProb4dMat) -> Prob4dMat {
  lbpap_mat.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect()
}

#[inline]
pub fn get_log_sta_inside_ppf_mat_sets(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_pos_dist_4_el: usize, max_substr_dist: usize) -> LogStaInsidePpfMatSets {
  let mut log_sta_inside_ppf_mat_sets = LogStaInsidePpfMatSets::new();
  for i in 1 .. seq_len_pair.0 {
    for j in 1 .. seq_len_pair.1 {
      let pos_pair = (i, j);
      let pos_quadruple = (i, i - 1, j, j - 1);
      if !is_pos_dist_ok(&pos_pair, max_pos_dist_4_il) {continue;}
      log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.insert(pos_quadruple, 0.);
      log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas.insert(pos_quadruple, 0.);
      log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas.insert(pos_quadruple, 0.);
    }
  }
  let leftmost_pos_pair = (0, 0);
  let rightmost_pos_pair = (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
  log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat.insert(leftmost_pos_pair, 0.);
  log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas.insert(leftmost_pos_pair, 0.);
  log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat.insert(rightmost_pos_pair, 0.);
  log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas.insert(rightmost_pos_pair, 0.);
  for substr_len_1 in 0 .. max_bp_span_pair.0 + 1 {
    for substr_len_2 in 0 .. max_bp_span_pair.1 + 1 {
      if substr_len_1 == substr_len_2 && substr_len_1 == 0 {
        continue;
      }
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_substr_dist {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !is_substr_dist_ok(&pos_quadruple, max_substr_dist, max_pos_dist_4_il) {continue;}
          if sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
            let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(i + 1, j - 1, k + 1, l - 1)) {
              let ep_of_term_4_log_pf = bpa_score - INVERSE_TEMPERATURE * (get_hl_fe(seq_pair.0, &(i, j)) + get_hl_fe(seq_pair.1, &(k, l))) + log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(i + 1, j - 1, k + 1, l - 1)];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
            for m in i + 1 .. j - 1 {
              for n in m + 1 .. j {
                if j - n - 1 + m - i - 1 > MAX_2_LOOP_LEN {continue;}
                for o in k + 1 .. l - 1 {
                  for p in o + 1 .. l {
                    if l - p - 1 + o - k - 1 > MAX_2_LOOP_LEN {continue;}
                    let pos_quadruple_2 = (m, n, o, p);
                    if !log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
                    if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(i + 1, m - 1, k + 1, o - 1)) && log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(n + 1, j - 1, p + 1, l - 1))) {continue;}
                    let two_loop_fe_pair = (
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.0, &(i, j), &(m , n)) as FreeEnergy,
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.1, &(k, l), &(o , p)) as FreeEnergy,
                    );
                    let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(i + 1, m - 1, k + 1, o - 1)] + bpa_score + two_loop_fe_pair.0 + two_loop_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple_2] + log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(n + 1, j - 1, p + 1, l - 1)];
                    if ep_of_term_4_log_pf.is_finite() {
                      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                      if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                        max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                      }
                    }
                  }
                }
              }
            }
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.contains_key(&(i + 1, j - 1, k + 1, l - 1)) {
              let bp_closing_loop_pair = (
                (seq_pair.0[i], seq_pair.0[j]),
                (seq_pair.1[k], seq_pair.1[l]),
              );
              let ml_tm_delta_fe_pair = (
                - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[j - 1])))],
                - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[l - 1])))],
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              let ep_of_term_4_log_pf = bpa_score - 2. * INVERSE_TEMPERATURE * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat[&(i + 1, j - 1, k + 1, l - 1)];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
            if max_ep_of_term_4_log_pf.is_finite() {
              log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
            }
          }
          // Compute "log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa"
          if sta_fe_params.ba_score_mat.contains_key(&(j, l)) {
            let ba_score = sta_fe_params.ba_score_mat[&(j, l)];
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(i, j - 1, k, l - 1)) {
              log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas.insert(pos_quadruple,
                log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(i, j - 1, k, l - 1)] + ba_score
              );
            }
          }
          let pos_quadruple_2 = (i, j - 1, k, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
              let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_1[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_2[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          if sta_fe_params.ba_score_mat.contains_key(&(i, k)) {
            let ba_score = sta_fe_params.ba_score_mat[&(i, k)];
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(i + 1, j, k + 1, l)) {
              log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas.insert(pos_quadruple,
                log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(i + 1, j, k + 1, l)] + ba_score
              );
            }
          }
          let pos_quadruple_2 = (i + 1, j, k, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_1[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_2[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_1[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_2[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_1[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_2[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          // Compute "log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop"
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          for m in i .. j {
            for n in k .. l {
              let pos_quadruple_2 = (i, m - 1, k, n - 1);
              let pos_quadruple_3 = (m, j, n, l);
              if !log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_3) {
                continue;
              }
              let accessible_bp_pair = (
                (seq_pair.0[m], seq_pair.0[j]),
                (seq_pair.1[n], seq_pair.1[l]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[m - 1], seq_pair.0[j + 1]),
                (seq_pair.1[n - 1], seq_pair.1[l + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if m > 1 && j < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if m > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if j < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                - INVERSE_TEMPERATURE * if n > 1 && l < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if n > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if l < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.contains_key(&pos_quadruple_2) {
                let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_quadruple_2] + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 - 2. * INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple_3];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.contains_key(&pos_quadruple_2) {
                let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat[&pos_quadruple_2] + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 - 2. * INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple_3];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          if sta_fe_params.ba_score_mat.contains_key(&(j, l)) {
            let ba_score = sta_fe_params.ba_score_mat[&(j, l)];
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.contains_key(&(i, j - 1, k, l - 1)) {
              let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat[&(i, j - 1, k, l - 1)] + ba_score;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_bas.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j - 1, k, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_1[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_2[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          let mut eps_of_terms_4_log_pf_2 = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf_2 = NEG_INFINITY;
          for m in i + 1 .. j + 1 {
            for n in k + 1 .. l + 1 {
              let pos_quadruple_2 = (i, m, k, n);
              let pos_quadruple_3 = (m + 1, j, n + 1, l);
              if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.contains_key(&pos_quadruple_3) && log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_2)) {
                continue;
              }
              let accessible_bp_pair = (
                (seq_pair.0[i], seq_pair.0[m]),
                (seq_pair.1[k], seq_pair.1[n]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[i - 1], seq_pair.0[m + 1]),
                (seq_pair.1[k - 1], seq_pair.1[n + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if i > 1 && m < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if i > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if m < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                - INVERSE_TEMPERATURE * if k > 1 && n < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if k > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if n < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.contains_key(&pos_quadruple_3) {
                let ep_of_term_4_log_pf = ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 - 2. * INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple_2] + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_quadruple_3];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if m != j && n != l {
                    eps_of_terms_4_log_pf_2.push(ep_of_term_4_log_pf);
                  }
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    if m != j && n != l {
                      max_ep_of_term_4_log_pf_2 = ep_of_term_4_log_pf;
                    }
                  }
                }
              }
              if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.contains_key(&pos_quadruple_3) {
                let ep_of_term_4_log_pf = ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 - 2. * INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple_2] + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat[&pos_quadruple_3];
                if ep_of_term_4_log_pf.is_finite() {
                  if m != j && n != l {
                    eps_of_terms_4_log_pf_2.push(ep_of_term_4_log_pf);
                  }
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    if m != j && n != l {
                      max_ep_of_term_4_log_pf_2 = ep_of_term_4_log_pf;
                    }
                  }
                }
              }
            }
          }
          if sta_fe_params.ba_score_mat.contains_key(&(i, k)) {
            let ba_score = sta_fe_params.ba_score_mat[&(i, k)];
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.contains_key(&(i + 1, j, k + 1, l)) {
              let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat[&(i + 1, j, k + 1, l)] + ba_score;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                eps_of_terms_4_log_pf_2.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  max_ep_of_term_4_log_pf_2 = ep_of_term_4_log_pf;
                }
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_bas.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let tmp_ep_of_term_4_log_pf = if max_ep_of_term_4_log_pf_2.is_finite() {logsumexp(&eps_of_terms_4_log_pf_2[..], max_ep_of_term_4_log_pf_2)} else {NEG_INFINITY};
          let pos_quadruple_2 = (i + 1, j, k, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_1[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_2[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_bas[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_1[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_2[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if tmp_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(tmp_ep_of_term_4_log_pf);
            if tmp_ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = tmp_ep_of_term_4_log_pf;
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_1[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_2[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          // Compute "log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls"
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          for m in i .. j {
            for n in k .. l {
              let pos_quadruple_2 = (i, m - 1, k, n - 1);
              let pos_quadruple_3 = (m, j, n, l);
              if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&pos_quadruple_2) && log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_3)) {
                continue;
              }
              let accessible_bp_pair = (
                (seq_pair.0[m], seq_pair.0[j]),
                (seq_pair.1[n], seq_pair.1[l]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[m - 1], seq_pair.0[j + 1]),
                (seq_pair.1[n - 1], seq_pair.1[l + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if m > 1 && j < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if m > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if j < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                - INVERSE_TEMPERATURE * if n > 1 && l < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if n > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if l < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&pos_quadruple_2] + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 - 2. * INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple_3];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          if sta_fe_params.ba_score_mat.contains_key(&(j, l)) {
            let ba_score = sta_fe_params.ba_score_mat[&(j, l)];
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.contains_key(&(i, j - 1, k, l - 1)) {
              let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(i, j - 1, k, l - 1)] + ba_score;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_bas.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j - 1, k, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_1[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_2[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          let mut eps_of_terms_4_log_pf_2 = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf_2 = NEG_INFINITY;
          for m in i + 1 .. j + 1 {
            for n in k + 1 .. l + 1 {
              let pos_quadruple_2 = (i, m, k, n);
              let pos_quadruple_3 = (m + 1, j, n + 1, l);
              if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&pos_quadruple_3) && log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_2)) {
                continue;
              }
              let accessible_bp_pair = (
                (seq_pair.0[i], seq_pair.0[m]),
                (seq_pair.1[k], seq_pair.1[n]),
              );
              let stacking_bp_pair = (
                (seq_pair.0[i - 1], seq_pair.0[m + 1]),
                (seq_pair.1[k - 1], seq_pair.1[n + 1]),
              );
              let ml_tm_or_de_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if i > 1 && m < seq_len_pair.0 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
                } else if i > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
                } else if m < seq_len_pair.0 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
                } else {
                  0.
                },
                - INVERSE_TEMPERATURE * if k > 1 && n < seq_len_pair.1 - 2 {
                  ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
                } else if k > 1 {
                  FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
                } else if n < seq_len_pair.1 - 2 {
                  THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
                } else {
                  0.
                },
              );
              let au_or_gu_end_penalty_delta_fe_pair = (
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              );
              let ep_of_term_4_log_pf = ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 - 2. * INVERSE_TEMPERATURE * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple_2] + log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&pos_quadruple_3];
              if ep_of_term_4_log_pf.is_finite() {
                if m != j && n != l {
                  eps_of_terms_4_log_pf_2.push(ep_of_term_4_log_pf);
                }
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  if m != j && n != l {
                    max_ep_of_term_4_log_pf_2 = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          if sta_fe_params.ba_score_mat.contains_key(&(i, k)) {
            let ba_score = sta_fe_params.ba_score_mat[&(i, k)];
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.contains_key(&(i + 1, j, k + 1, l)) {
              let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(i + 1, j, k + 1, l)] + ba_score;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                eps_of_terms_4_log_pf_2.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  max_ep_of_term_4_log_pf_2 = ep_of_term_4_log_pf;
                }
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_bas.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let tmp_ep_of_term_4_log_pf = if max_ep_of_term_4_log_pf_2.is_finite() {logsumexp(&eps_of_terms_4_log_pf_2[..], max_ep_of_term_4_log_pf_2)} else {NEG_INFINITY};
          let pos_quadruple_2 = (i + 1, j, k, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_1[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_2[&pos_quadruple_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_bas[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_1[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_2[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if tmp_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(tmp_ep_of_term_4_log_pf);
            if tmp_ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = tmp_ep_of_term_4_log_pf;
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_1[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_2[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
        }
      }
    }
  }
  for i in 0 .. seq_len_pair.0 - 1 {
    for j in 0 .. seq_len_pair.1 - 1 {
      let pos_pair = (i, j);
      if pos_pair == (0, 0) {continue;}
      if !is_pos_dist_ok(&pos_pair, max_pos_dist_4_el) {continue;}
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in 1 .. i {
        for l in 1 .. j {
          let pos_pair_2 = (k - 1, l - 1);
          let pos_quadruple = (k, i, l, j);
          if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat.contains_key(&pos_pair_2) && log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple)) {
            continue;
          }
          let accessible_bp_pair = (
            (seq_pair.0[k], seq_pair.0[i]),
            (seq_pair.1[l], seq_pair.1[j]),
          );
          let stacking_bp_pair = (
            (seq_pair.0[k - 1], seq_pair.0[i + 1]),
            (seq_pair.1[l - 1], seq_pair.1[j + 1]),
          );
          let ml_tm_or_de_delta_fe_pair = (
            - INVERSE_TEMPERATURE * if k > 1 && i < seq_len_pair.0 - 2 {
              ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
            } else if k > 1 {
              FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
            } else if i < seq_len_pair.0 - 2 {
              THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
            } else {
              0.
            },
            - INVERSE_TEMPERATURE * if l > 1 && j < seq_len_pair.1 - 2 {
              ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
            } else if l > 1 {
              FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
            } else if j < seq_len_pair.1 - 2 {
              THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
            } else {
              0.
            },
          );
          let au_or_gu_end_penalty_delta_fe_pair = (
            - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
            - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
          );
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat[&pos_pair_2] + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple];
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
      }
      if i > 0 && j > 0 && log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat.contains_key(&(i - 1, j - 1)) {
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat[&(i - 1, j - 1)] + ba_score;
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if i > 0 {
        let pos_pair_2 = (i - 1, j);
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_1.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if j > 0 {
        let pos_pair_2 = (i, j - 1);
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_2.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas.contains_key(&pos_pair) {
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_bas[&pos_pair];
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_1.contains_key(&pos_pair) {
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_1[&pos_pair];
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_2.contains_key(&pos_pair) {
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat_4_gaps_2[&pos_pair];
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
    }
  }
  for i in (1 .. seq_len_pair.0).rev() {
    for j in (1 .. seq_len_pair.1).rev() {
      let pos_pair = (i, j);
      if pos_pair == (seq_len_pair.0 - 1, seq_len_pair.1 - 1) {continue;}
      if !is_pos_dist_ok(&pos_pair, max_pos_dist_4_el) {continue;}
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in i + 1 .. seq_len_pair.0 - 1 {
        for l in j + 1 .. seq_len_pair.1 - 1 {
          let pos_pair_2 = (k + 1, l + 1);
          let pos_quadruple = (i, k, j, l);
          if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat.contains_key(&pos_pair_2) && log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple)) {
            continue;
          }
          let accessible_bp_pair = (
            (seq_pair.0[i], seq_pair.0[k]),
            (seq_pair.1[j], seq_pair.1[l]),
          );
          let stacking_bp_pair = (
            (seq_pair.0[i - 1], seq_pair.0[k + 1]),
            (seq_pair.1[j - 1], seq_pair.1[l + 1]),
          );
          let ml_tm_or_de_delta_fe_pair = (
            - INVERSE_TEMPERATURE * if i > 1 && k < seq_len_pair.0 - 2 {
              ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
            } else if i > 1 {
              FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
            } else if k < seq_len_pair.0 - 2 {
              THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
            } else {
              0.
            },
            - INVERSE_TEMPERATURE * if j > 1 && l < seq_len_pair.1 - 2 {
              ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
            } else if j > 1 {
              FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
            } else if l < seq_len_pair.1 - 2 {
              THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
            } else {
              0.
            },
          );
          let au_or_gu_end_penalty_delta_fe_pair = (
            - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
            - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
          );
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat[&pos_pair_2] + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[&pos_quadruple];
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
      }
      if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 && log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat.contains_key(&(i + 1, j + 1)) {
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat[&(i + 1, j + 1)] + ba_score;
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if i < seq_len_pair.0 - 1 {
        let pos_pair_2 = (i + 1, j);
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_1.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if j < seq_len_pair.1 - 1 {
        let pos_pair_2 = (i, j + 1);
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
        if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_pair_2) {
          let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
            if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
              max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
            }
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_2.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas.contains_key(&pos_pair) {
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_bas[&pos_pair];
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_1.contains_key(&pos_pair) {
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_1[&pos_pair];
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_2.contains_key(&pos_pair) {
        let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat_4_gaps_2[&pos_pair];
        if ep_of_term_4_log_pf.is_finite() {
          eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
      }
      if max_ep_of_term_4_log_pf.is_finite() {
        log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat.insert(pos_pair, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
      }
    }
  }
  log_sta_inside_ppf_mat_sets
}

#[inline]
pub fn get_log_sta_outside_ppf_4d_mats(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_substr_dist: usize, log_sta_inside_ppf_mat_sets: &LogStaInsidePpfMatSets) -> LogStaOutsidePpf4dMats {
  let mut log_sta_outside_ppf_4d_mats = LogStaOutsidePpf4dMats::new();
  for substr_len_1 in (2 .. max_bp_span_pair.0 + 1).rev() {
    for substr_len_2 in (2 .. max_bp_span_pair.1 + 1).rev() {
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_substr_dist {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !is_substr_dist_ok(&pos_quadruple, max_substr_dist, max_pos_dist_4_il) {continue;}
          // Compute "log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas"
          if log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple) {
            let accessible_bp_pair = (
              (seq_pair.0[i], seq_pair.0[j]),
              (seq_pair.1[k], seq_pair.1[l]),
            );
            let stacking_bp_pair = (
              (seq_pair.0[i - 1], seq_pair.0[j + 1]),
              (seq_pair.1[k - 1], seq_pair.1[l + 1]),
            );
            let ml_tm_or_de_delta_fe_pair = (
              - INVERSE_TEMPERATURE * if i > 1 && j < seq_len_pair.0 - 2 {
                ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
              } else if i > 1 {
                FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
              } else if j < seq_len_pair.0 - 2 {
                THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
              } else {
                0.
              },
              - INVERSE_TEMPERATURE * if k > 1 && l < seq_len_pair.1 - 2 {
                ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
              } else if k > 1 {
                FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
              } else if l < seq_len_pair.1 - 2 {
                THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
              } else {
                0.
              },
            );
            let au_or_gu_end_penalty_delta_fe_pair = (
              - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              - INVERSE_TEMPERATURE * if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
            );
            if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat.contains_key(&(i - 1, k - 1)) && log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat.contains_key(&(j + 1, l + 1)) {
              log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_el.insert(pos_quadruple, log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat[&(i - 1, k - 1)] + ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat[&(j + 1, l + 1)]);
            }
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
            for m in 1 .. i {
              for n in j + 1 .. seq_len_pair.0 - 1 {
                if n - j - 1 + i - m - 1 > MAX_2_LOOP_LEN {continue;}
                for o in 1 .. k {
                  for p in l + 1 .. seq_len_pair.1 - 1 {
                    if p - l - 1 + k - o - 1 > MAX_2_LOOP_LEN {continue;}
                    let pos_quadruple_2 = (m, n, o, p);
                    if !log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
                    if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(m + 1, i - 1, o + 1, k - 1)) && log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(j + 1, n - 1, l + 1, p - 1))) {continue;}
                    let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                    let two_loop_fe_pair = (
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.0, &(m, n), &(i , j)) as FreeEnergy,
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.1, &(o, p), &(k , l)) as FreeEnergy,
                    );
                    let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(m + 1, i - 1, o + 1, k - 1)] + bpa_score + two_loop_fe_pair.0 + two_loop_fe_pair.1 + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas[&pos_quadruple_2] + log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(j + 1, n - 1, l + 1, p - 1)];
                    if ep_of_term_4_log_pf.is_finite() {
                      eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                      if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                        max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                      }
                    }
                  }
                }
              }
            }
            if max_ep_of_term_4_log_pf.is_finite() {
              log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_internal_2loops.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
            }
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
            for m in j + 1 .. seq_len_pair.0 - 1 {
              for n in l + 1 .. seq_len_pair.1 - 1 {
                let pos_quadruple_2 = (i, m, k, n);
                if !log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
                let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                let log_coefficient = bpa_score - 2. * INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
                let bp_closing_loop_pair = (
                  (seq_pair.0[i], seq_pair.0[m]),
                  (seq_pair.1[k], seq_pair.1[n]),
                );
                let ml_tm_delta_fe_pair = (
                  - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[m - 1])))],
                  - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[n - 1])))],
                );
                let au_or_gu_end_penalty_delta_fe_pair = (
                  - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                  - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                );
                if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.contains_key(&(j + 1, m - 1, l + 1, n - 1)) {
                  let ep_of_term_4_log_pf = log_coefficient + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas[&pos_quadruple_2] + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat[&(j + 1, m - 1, l + 1, n - 1)];
                  if ep_of_term_4_log_pf.is_finite() {
                    eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                      max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    }
                  }
                }
                if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.contains_key(&(j + 1, m - 1, l + 1, n - 1)) {
                  let ep_of_term_4_log_pf = log_coefficient + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas[&pos_quadruple_2] + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(j + 1, m - 1, l + 1, n - 1)];
                  if ep_of_term_4_log_pf.is_finite() {
                    eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                      max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    }
                  }
                }
              }
            }
            if max_ep_of_term_4_log_pf.is_finite() {
              log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
            }
            for m in j + 1 .. seq_len_pair.0 - 1 {
              for n in l + 1 .. seq_len_pair.1 - 1 {
                let pos_quadruple_2 = (i, m, k, n);
                if !(log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas.contains_key(&pos_quadruple_2) && log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(j + 1, m - 1, l + 1, n - 1))) {continue;}
                let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                let log_coefficient = bpa_score - 2. * INVERSE_TEMPERATURE * (CONST_4_INIT_ML_DELTA_FE + COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
                let bp_closing_loop_pair = (
                  (seq_pair.0[i], seq_pair.0[m]),
                  (seq_pair.1[k], seq_pair.1[n]),
                );
                let ml_tm_delta_fe_pair = (
                  - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[m - 1])))],
                  - INVERSE_TEMPERATURE * ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[n - 1])))],
                );
                let au_or_gu_end_penalty_delta_fe_pair = (
                  - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                  - INVERSE_TEMPERATURE * if is_au_or_gu(&bp_closing_loop_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
                );
                let ep_of_term_4_log_pf = log_coefficient + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas[&pos_quadruple_2] + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(j + 1, m - 1, l + 1, n - 1)];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
            if max_ep_of_term_4_log_pf.is_finite() {
              log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
            }
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
            for m in 1 .. i {
              for n in 1 .. k {
                let pos_quadruple_2 = (m, j, n, l);
                if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat.contains_key(&(m + 1, i - 1, n + 1, k - 1)) && log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right.contains_key(&pos_quadruple_2) {
                  let ep_of_term_4_log_pf = ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_ppf_mat[&(m + 1, i - 1, n + 1, k - 1)] + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right[&pos_quadruple_2];
                  if ep_of_term_4_log_pf.is_finite() {
                    eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                      max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    }
                  }
                }
                if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.contains_key(&(m + 1, i - 1, n + 1, k - 1)) && log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right.contains_key(&pos_quadruple_2) {
                  let ep_of_term_4_log_pf = ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&(m + 1, i - 1, n + 1, k - 1)] + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right[&pos_quadruple_2];
                  if ep_of_term_4_log_pf.is_finite() {
                    eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                      max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    }
                  }
                }
                if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(m + 1, i - 1, n + 1, k - 1)) && log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right_2.contains_key(&pos_quadruple_2) {
                  let ep_of_term_4_log_pf = ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(m + 1, i - 1, n + 1, k - 1)] + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_right_2[&pos_quadruple_2];
                  if ep_of_term_4_log_pf.is_finite() {
                    eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                    if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                      max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                    }
                  }
                }
              }
            }
            if max_ep_of_term_4_log_pf.is_finite() {
              log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_internal_mls.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
            }
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
            if log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_el.contains_key(&pos_quadruple) {
              let ep_of_term_4_log_pf = log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_el[&pos_quadruple];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
            if log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_internal_2loops.contains_key(&pos_quadruple) {
              let ep_of_term_4_log_pf = log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_internal_2loops[&pos_quadruple];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
            if log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_internal_mls.contains_key(&pos_quadruple) {
              let ep_of_term_4_log_pf = log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas_on_internal_mls[&pos_quadruple];
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
            }
            if max_ep_of_term_4_log_pf.is_finite() {
              log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
            }
          }
        }
      }
    }
  }
  log_sta_outside_ppf_4d_mats
}

#[inline]
fn get_lbpap_mat(seq_len_pair: &(usize, usize), log_sta_inside_ppf_mat_sets: &LogStaInsidePpfMatSets, log_sta_outside_ppf_4d_mats: &LogStaOutsidePpf4dMats) -> LogProb4dMat {
  let mut lbpap_mat = LogProb4dMat::default();
  let log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat[&(seq_len_pair.0 - 2, seq_len_pair.1 - 2)];
  for pos_quadruple in log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.keys() {
    if !log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas.contains_key(pos_quadruple) {continue;}
    let lbpap = log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[pos_quadruple] + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas[pos_quadruple] - log_pf;
    if lbpap.is_finite() {
      debug_assert!(NEG_INFINITY < lbpap && lbpap <= 0.);
      lbpap_mat.insert(*pos_quadruple, lbpap);
    }
  }
  lbpap_mat
}
 
#[inline]
pub fn pct_of_bpp_and_upp_mat(lbpap_mats_with_rna_id_pairs: &LogProb4dMatsWithRnaIdPairs, rna_id: RnaId, num_of_rnas: usize, bpp_mat: &SparseProbMat, upp_mat: &Probs) -> (SparseProbMat, Probs) {
  let log_coefficient = -((num_of_rnas - 1) as Prob).ln();
  let upp_mat_len = upp_mat.len();
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut max_eps_of_terms_4_log_probs_with_pos_pairs = EpsOfTerms4LogProbsWithPosPairs::default();
  let mut seqs_of_eps_of_terms_4_log_probs_with_poss = vec![EpsOfTerms4LogProb::new(); upp_mat_len];
  let mut max_eps_of_terms_4_log_probs_with_poss = vec![NEG_INFINITY; upp_mat_len];
  for pos_pair in bpp_mat.keys() {
    if pos_pair.0 == 0 {continue;}
    seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, EpsOfTerms4LogProb::new());
    max_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, NEG_INFINITY);
  }
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id == rna_id_2 {continue;}
    let rna_id_pair = if rna_id < rna_id_2 {(rna_id, rna_id_2)} else {(rna_id_2, rna_id)};
    let ref ref_2_lbpap_mat = lbpap_mats_with_rna_id_pairs[&rna_id_pair];
    for (pos_quadruple, &lbpap) in ref_2_lbpap_mat.iter() {
      if pos_quadruple.0 == 0 {continue;}
      if lbpap.is_finite() {
        let pos_pair = if rna_id < rna_id_2 {(pos_quadruple.0, pos_quadruple.1)} else {(pos_quadruple.2, pos_quadruple.3)};
        let eps_of_terms_4_log_prob = seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
        eps_of_terms_4_log_prob.push(lbpap);
        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
        if lbpap > *max_ep_of_term_4_log_prob {
          *max_ep_of_term_4_log_prob = lbpap;
        }
        seqs_of_eps_of_terms_4_log_probs_with_poss[pos_pair.0].push(lbpap);
        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_poss[pos_pair.0];
        if lbpap > max_ep_of_term_4_log_prob {
          max_eps_of_terms_4_log_probs_with_poss[pos_pair.0] = lbpap;
        }
        seqs_of_eps_of_terms_4_log_probs_with_poss[pos_pair.1].push(lbpap);
        let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_poss[pos_pair.1];
        if lbpap > max_ep_of_term_4_log_prob {
          max_eps_of_terms_4_log_probs_with_poss[pos_pair.1] = lbpap;
        }
      }
    }
  }
  let mut bpp_mat = bpp_mat.clone();
  for (pos_pair, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.iter() {
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_pos_pairs[pos_pair];
    if max_ep_of_term_4_log_prob.is_finite() {
      let new_bpp = (log_coefficient + logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob)).exp();
      *bpp_mat.get_mut(pos_pair).expect("Failed to get an element of a hash map.") = new_bpp;
    }
  }
  let mut upp_mat = vec![1.; upp_mat_len];
  for (i, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_poss.iter().enumerate() {
    if i == 0 || i == upp_mat_len - 1 {continue;}
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_poss[i];
    if max_ep_of_term_4_log_prob.is_finite() {
      let new_upp = 1. - (log_coefficient + logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob)).exp();
      upp_mat[i] = new_upp;
    }
  }
  (bpp_mat, upp_mat)
}

#[inline]
pub fn get_bpap_mats_with_rna_id_pairs(lbpap_mats_with_rna_id_pairs: &LogProb4dMatsWithRnaIdPairs) -> Prob4dMatsWithRnaIdPairs {
  let mut bpap_mats_with_rna_id_pairs = Prob4dMatsWithRnaIdPairs::default();
  for (rna_id_pair, lbpap_mat) in lbpap_mats_with_rna_id_pairs.iter() {
    bpap_mats_with_rna_id_pairs.insert(*rna_id_pair, get_bpap_mat(lbpap_mat));
  }
  bpap_mats_with_rna_id_pairs
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
pub fn remove_zeros_from_bpp_mat(bpp_mat: &ProbMat) -> SparseLogProbMat {
  let mut sparse_bpp_mat = SparseProbMat::default();
  for (i, bpps) in bpp_mat.iter().enumerate() {
    for (j, &bpp) in bpps.iter().enumerate() {
      if bpp > 0. {
        sparse_bpp_mat.insert((i + 1, j + 1), bpp);
      }
    }
  }
  sparse_bpp_mat
}

#[inline]
pub fn remove_small_bpps_from_bpp_mat(sparse_bpp_mat: &SparseProbMat, min_bpp: Prob) -> SparseProbMat {
  sparse_bpp_mat.iter().filter(|(_, &bpp)| {bpp >= min_bpp}).map(|(pos_pair, &bpp)| {(*pos_pair, bpp)}).collect()
}

#[inline]
pub fn get_max_bp_span(sparse_bpp_mat: &SparseProbMat) -> usize {
  sparse_bpp_mat.iter().map(|(pos_pair, _)| {pos_pair.1 - pos_pair.0 + 1}).max().unwrap()
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
fn is_pos_dist_ok(pos_pair: &PosPair, max_pos_dist: usize) -> bool {
  if max(pos_pair.0, pos_pair.1) - min(pos_pair.0, pos_pair.1) <= max_pos_dist {true} else {false}
}

#[inline]
fn is_substr_dist_ok(pos_quadruple: &PosQuadruple, max_substr_dist: usize, max_pos_dist_4_il: usize) -> bool {
  let pos_pair_1 = (pos_quadruple.0, pos_quadruple.2);
  let pos_pair_2 = (pos_quadruple.1, pos_quadruple.3);
  if pos_quadruple.1 - pos_quadruple.0 - (pos_quadruple.3 - pos_quadruple.2) <= max_substr_dist  && is_pos_dist_ok(&pos_pair_1, max_pos_dist_4_il) && is_pos_dist_ok(&pos_pair_2, max_pos_dist_4_il) {true} else {false}
}
