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
  pub fn new(rna_id_pair: &RnaIdPair, fasta_records: &FastaRecords, max_bp_span_pair: &(usize, usize), max_gap_num: usize, bpp_mats: &ProbMats, sta_fe_scale_param: StaFreeEnergy, opening_gap_penalty: StaFreeEnergy, extending_gap_penalty: StaFreeEnergy) -> StaFeParams {
    let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
    let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
    let max_gap_num_2 = get_seq_len_diff(&seq_len_pair) + max_gap_num;
    let mut sta_fe_params = StaFeParams::origin();
    sta_fe_params.opening_gap_penalty = sta_fe_scale_param * opening_gap_penalty;
    sta_fe_params.extending_gap_penalty = sta_fe_scale_param * extending_gap_penalty;
    let bpp_mat_pair = (&bpp_mats[rna_id_pair.0], &bpp_mats[rna_id_pair.1]);
    let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
    for i in 1 .. seq_len_pair.0 - 1 {
      let base = seq_pair.0[i];
      for j in 1 .. seq_len_pair.1 - 1 {
        let pos_pair = (i, j);
        if !is_min_gap_ok_1(&pos_pair, /* &pseudo_pos_quadruple, */ max_gap_num_2) {continue;}
        let base_pair = (base, seq_pair.1[j]);
        sta_fe_params.ba_score_mat.insert(pos_pair, INVERSE_TEMPERATURE * sta_fe_scale_param * RIBOSUM_85_60_BA_SCORE_MAT[&base_pair]);
      }
      let upper_j = if i + max_bp_span_pair.0 >= seq_len_pair.0 - 1 {seq_len_pair.0 - 1} else {i + max_bp_span_pair.0};
      for j in i + 1 .. upper_j {
        let pos_pair = (i, j);
        let base_pair = (base, seq_pair.0[j]);
        if !bpp_mat_pair.0.contains_key(&pos_pair) {continue;}
        for k in 1 .. seq_len_pair.1 - 1 {
          if !is_min_gap_ok_1(&(i, k), /* &pseudo_pos_quadruple, */ max_gap_num) {continue;}
          let upper_l = if k + max_bp_span_pair.1 >= seq_len_pair.1 - 1 {seq_len_pair.1 - 1} else {k + max_bp_span_pair.1};
          for l in k + 1 .. upper_l {
            if !is_min_gap_ok_1(&(j, l), /* &pseudo_pos_quadruple, */ max_gap_num) {continue;}
            let pos_pair_2 = (k, l);
            let pos_quadruple = (pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1);
            if !(is_min_gap_ok_2(&pos_quadruple, /* &pseudo_pos_quadruple, */ max_gap_num) && bpp_mat_pair.1.contains_key(&pos_pair_2)) {continue;}
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

/* impl LogStaPpfMatSets {
  pub fn new() -> LogStaPpfMatSets {
    let log_ppf_mats = LogStaPpfMats::new();
    LogStaPpfMatSets {
      log_ppf_mats_on_sas: log_ppf_mats.clone(),
      log_ppf_mats_on_2ls: log_ppf_mats.clone(),
      log_ppf_mats_on_mls: log_ppf_mats.clone(),
      log_ppf_mats_4_first_bpas_on_mls: log_ppf_mats,
    }
  }
} */

#[inline]
pub fn io_algo_4_base_pair_align_prob_mat(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize) -> Prob4dMat {
  let lbpap_mat = io_algo_4_lbpap_mat(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num);
  get_bpap_mat(&lbpap_mat)
}

#[inline]
pub fn io_algo_4_lbpap_mat(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize) -> LogProb4dMat {
  let log_sta_inside_ppf_mat_sets = get_log_sta_inside_ppf_mat_sets(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num);
  println!("An inside step finished.");
  let log_sta_outside_ppf_4d_mats = get_log_sta_outside_ppf_4d_mats(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, &log_sta_inside_ppf_mat_sets);
  println!("An outside step finished.");
  get_lbpap_mat(seq_len_pair, &log_sta_inside_ppf_mat_sets, &log_sta_outside_ppf_4d_mats)
}

#[inline]
fn get_bpap_mat(lbpap_mat: &LogProb4dMat) -> Prob4dMat {
  lbpap_mat.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect()
}

#[inline]
pub fn get_log_sta_inside_ppf_mat_sets(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize) -> LogStaInsidePpfMatSets {
  let mut log_sta_inside_ppf_mat_sets = LogStaInsidePpfMatSets::new();
  let max_gap_num_2 = get_seq_len_diff(&seq_len_pair) + max_gap_num;
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  for i in 1 .. seq_len_pair.0 {
    for j in 1 .. seq_len_pair.1 {
      let pos_pair = (i, j);
      let pos_quadruple = (i, i - 1, j, j - 1);
      if !is_min_gap_ok_1(&pos_pair, /* &pseudo_pos_quadruple, */ max_gap_num_2) {continue;}
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
  println!("seq_len_pair: {:?}.", &seq_len_pair);
  println!("max_bp_span_pair: {:?}.", &max_bp_span_pair);
  for substr_len_1 in 0 .. max_bp_span_pair.0 + 1 {
    for substr_len_2 in 0 .. max_bp_span_pair.1 + 1 {
      if substr_len_1 == substr_len_2 && substr_len_1 == 0 {
        continue;
      }
      let remain_substr_len_pair = (
        seq_len_pair.0 - substr_len_1,
        seq_len_pair.1 - substr_len_2,
      );
      // if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) + max(remain_substr_len_pair.0, remain_substr_len_pair.1) - min(remain_substr_len_pair.0, remain_substr_len_pair.1) > max_gap_num {continue;}
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_gap_num {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          // if !(sta_fe_params.ba_score_mat.contains_key(&(i, k)) && sta_fe_params.ba_score_mat.contains_key(&(j, l))) {continue;}
          let pos_quadruple = (i, j, k, l);
          /* let min_gap_num_1 = get_min_gap_num(&(0, i, 0, k));
          let min_gap_num_2 = get_min_gap_num(&pos_quadruple);
          let min_gap_num_3 = get_min_gap_num(&(j, seq_len_pair.0, l, seq_len_pair.1));
          let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3;
          if total_min_gap_num > max_gap_num {continue;} */
          if !is_min_gap_ok_2(&pos_quadruple, max_gap_num) {continue;}
          // println!("total_min_gap_num: {}, max_gap_num: {}.", total_min_gap_num, max_gap_num);
          if sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
            let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
            // let ep_of_term_4_log_pf = bpa_score - INVERSE_TEMPERATURE * (if pos_quadruple == pseudo_pos_quadruple {0.} else {get_hl_fe(seq_pair.0, &(i, j)) + get_hl_fe(seq_pair.1, &(k, l))}) + log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&(j - 1, l - 1)];
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
                    /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple_2);
                    let min_gap_num_2 = get_min_gap_num(&(i + 1, m - 1, k + 1, o - 1));
                    let min_gap_num_3 = get_min_gap_num(&(n + 1, j - 1, p + 1, l - 1));
                    let min_gap_num_4 = get_min_gap_num(&(0, i, 0, k));
                    let min_gap_num_5 = get_min_gap_num(&(j, seq_len_pair.0, l, seq_len_pair.1));
                    let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4 + min_gap_num_5;
                    if total_min_gap_num > max_gap_num {continue;} */
                    let two_loop_fe_pair = (
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.0, &(i, j), &(m , n)) as Energy,
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.1, &(k, l), &(o , p)) as Energy,
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
          // println!("i, j, k, l: {:?}.", &pos_quadruple);
          // println!("OK.");
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
          /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple_2);
          let min_gap_num_2 = get_min_gap_num(&(0, i, 0, k));
          let min_gap_num_3 = get_min_gap_num(&(j - 1, seq_len_pair.0, l, seq_len_pair.1));
          let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + 1;
          if total_min_gap_num <= max_gap_num { */
          if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas.contains_key(&pos_quadruple_2) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas[&pos_quadruple_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            // }
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
          // if total_min_gap_num <= max_gap_num {
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
          // }
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
          /* log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas.insert(pos_quadruple, if log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat.contains_key(&(i + 1, j, k + 1, l)) {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_ppf_mat[&(i + 1, j, k + 1, l)] + ba_score_pair.0
          } else {
            NEG_INFINITY
          }); */
          let pos_quadruple_2 = (i + 1, j, k, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          // if total_min_gap_num <= max_gap_num {
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
          // }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          // if total_min_gap_num <= max_gap_num {
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
          // }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_2.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          /* let eps_of_terms_4_log_pf = [
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_bas[&pos_quadruple],
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_1[&pos_quadruple],
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_forward_ppf_mat_4_gaps_2[&pos_quadruple],
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_bas[&pos_quadruple],
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_1[&pos_quadruple],
            log_sta_inside_ppf_mat_sets.log_ppf_mats_on_sa.log_backward_ppf_mat_4_gaps_2[&pos_quadruple],
          ]; */
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
              /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple_2);
              let min_gap_num_2 = get_min_gap_num(&pos_quadruple_3);
              let min_gap_num_3 = get_min_gap_num(&(0, i, 0, k));
              let min_gap_num_4 = get_min_gap_num(&(j, seq_len_pair.0, l, seq_len_pair.1));
              let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4;
              if total_min_gap_num > max_gap_num {continue;} */
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
          // if total_min_gap_num <= max_gap_num {
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
          // }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_forward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          // if total_min_gap_num <= max_gap_num {
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
          // }
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
              /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple_2);
              let min_gap_num_2 = get_min_gap_num(&pos_quadruple_3);
              let min_gap_num_3 = get_min_gap_num(&(0, i, 0, k));
              let min_gap_num_4 = get_min_gap_num(&(j, seq_len_pair.0, l, seq_len_pair.1));
              let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4;
              if total_min_gap_num > max_gap_num {continue;} */
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
          // if total_min_gap_num <= max_gap_num {
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
          // }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_internal_multiloop.log_backward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          // if total_min_gap_num <= max_gap_num {
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
          // }
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
              /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple_2);
              let min_gap_num_2 = get_min_gap_num(&pos_quadruple_3);
              let min_gap_num_3 = get_min_gap_num(&(0, i, 0, k));
              let min_gap_num_4 = get_min_gap_num(&(j, seq_len_pair.0, l, seq_len_pair.1));
              let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4;
              if total_min_gap_num > max_gap_num {continue;} */
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
          // if total_min_gap_num <= max_gap_num {
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
          // }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_forward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          // if total_min_gap_num <= max_gap_num {
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
          // }
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
              /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple_2);
              let min_gap_num_2 = get_min_gap_num(&pos_quadruple_3);
              let min_gap_num_3 = get_min_gap_num(&(0, i, 0, k));
              let min_gap_num_4 = get_min_gap_num(&(j, seq_len_pair.0, l, seq_len_pair.1));
              let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4;
              if total_min_gap_num > max_gap_num {continue;} */
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
          // if total_min_gap_num <= max_gap_num {
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
          // }
          if max_ep_of_term_4_log_pf.is_finite() {
            log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_gaps_1.insert(pos_quadruple, logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf));
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
          // if total_min_gap_num <= max_gap_num {
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
          // }
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
          /* if log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_bas.contains_key(&pos_quadruple) {
            let ep_of_term_4_log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_backward_ppf_mat_4_bas[&pos_quadruple];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          } */
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
      if !is_min_gap_ok_1(&pos_pair, /* &pseudo_pos_quadruple, */ max_gap_num_2) {continue;}
      // println!("pos_pair: {:?}.", &pos_pair);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in 1 .. i {
        for l in 1 .. j {
          let pos_pair_2 = (k - 1, l - 1);
          let pos_quadruple = (k, i, l, j);
          if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat.contains_key(&pos_pair_2) && log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple)) {
            continue;
          }
          /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple);
          let min_gap_num_2 = get_min_gap_num(&(0, k, 0, l));
          let min_gap_num_3 = get_min_gap_num(&(i, seq_len_pair.0, j, seq_len_pair.1));
          let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3;
          if total_min_gap_num > max_gap_num {continue;} */
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
      /* let min_gap_num_1 = get_min_gap_num(&(0, i - 1, 0, j));
      let min_gap_num_2 = get_min_gap_num(&(i - 1, seq_len_pair.0, j, seq_len_pair.1));
      let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + 1; */
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
      if !is_min_gap_ok_1(&pos_pair, /* &pseudo_pos_quadruple, */ max_gap_num_2) {continue;}
      // println!("pos_pair: {:?}.", &pos_pair);
      let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
      let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
      for k in i + 1 .. seq_len_pair.0 - 1 {
        for l in j + 1 .. seq_len_pair.1 - 1 {
          let pos_pair_2 = (k + 1, l + 1);
          let pos_quadruple = (i, k, j, l);
          if !(log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_backward_ppf_mat.contains_key(&pos_pair_2) && log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.contains_key(&pos_quadruple)) {
            continue;
          }
          /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple);
          let min_gap_num_2 = get_min_gap_num(&(0, i, 0, j));
          let min_gap_num_3 = get_min_gap_num(&(k, seq_len_pair.0, l, seq_len_pair.1));
          let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3;
          if total_min_gap_num > max_gap_num {continue;} */
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
      /* let min_gap_num_1 = get_min_gap_num(&(0, i + 1, 0, j));
      let min_gap_num_2 = get_min_gap_num(&(i + 1, seq_len_pair.0, j, seq_len_pair.1));
      let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + 1; */
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
pub fn get_log_sta_outside_ppf_4d_mats(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize, log_sta_inside_ppf_mat_sets: &LogStaInsidePpfMatSets) -> LogStaOutsidePpf4dMats {
  let mut log_sta_outside_ppf_4d_mats = LogStaOutsidePpf4dMats::new();
  for substr_len_1 in (2 .. max_bp_span_pair.0 + 1).rev() {
    for substr_len_2 in (2 .. max_bp_span_pair.1 + 1).rev() {
      /* let remain_substr_len_pair = (
        seq_len_pair.0 - substr_len_1,
        seq_len_pair.1 - substr_len_2,
      ); */
      // if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) + max(remain_substr_len_pair.0, remain_substr_len_pair.1) - min(remain_substr_len_pair.0, remain_substr_len_pair.1) > max_gap_num {continue;}
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_gap_num {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !is_min_gap_ok_2(&pos_quadruple, max_gap_num) {continue;}
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
                    /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple);
                    let min_gap_num_2 = get_min_gap_num(&(m + 1, i - 1, o + 1, k - 1));
                    let min_gap_num_3 = get_min_gap_num(&(j + 1, n - 1, l + 1, p - 1));
                    let min_gap_num_4 = get_min_gap_num(&(0, m, 0, o));
                    let min_gap_num_5 = get_min_gap_num(&(n, seq_len_pair.0, p, seq_len_pair.1));
                    let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4 + min_gap_num_5;
                    if total_min_gap_num > max_gap_num {continue;} */
                    let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                    let two_loop_fe_pair = (
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.0, &(m, n), &(i , j)) as Energy,
                      - INVERSE_TEMPERATURE * get_2_loop_fe(seq_pair.1, &(o, p), &(k , l)) as Energy,
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
                /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple);
                let min_gap_num_2 = get_min_gap_num(&(j + 1, m - 1, l + 1, n - 1));
                let min_gap_num_3 = get_min_gap_num(&(0, i, 0, k));
                let min_gap_num_4 = get_min_gap_num(&(m, seq_len_pair.0, n, seq_len_pair.1));
                let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4;
                if total_min_gap_num > max_gap_num {continue;} */
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
                /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple);
                let min_gap_num_2 = get_min_gap_num(&(j + 1, m - 1, l + 1, n - 1));
                let min_gap_num_3 = get_min_gap_num(&(0, i, 0, k));
                let min_gap_num_4 = get_min_gap_num(&(m, seq_len_pair.0, n, seq_len_pair.1));
                let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4;
                if total_min_gap_num > max_gap_num {continue;} */
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
                /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple);
                let min_gap_num_2 = get_min_gap_num(&(m + 1, i - 1, n + 1, k - 1));
                let min_gap_num_3 = get_min_gap_num(&(0, m, 0, n));
                let min_gap_num_4 = get_min_gap_num(&(j, seq_len_pair.0, l, seq_len_pair.1));
                let total_min_gap_num = min_gap_num_1 + min_gap_num_2 + min_gap_num_3 + min_gap_num_4;
                if total_min_gap_num > max_gap_num {continue;} */
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

/* #[inline]
fn get_log_sta_forward_ppf_mat_sets(pos_quadruple: &PosQuadruple, pseudo_pos_quadruple: &PosQuadruple, seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, log_sta_ppf_mat_4_bpas: &LogPpf4dMat, max_gap_num: usize) -> LogStaPpfMatSets {
  let &(i, j, k, l) = pos_quadruple;
  let mut log_sta_ppf_mat_sets = LogStaPpfMatSets::new();
  for n in i .. j {
    for p in k .. l {
      let pos_pair_1 = (n, p);
      if !(is_min_gap_ok_1(&pos_pair_1, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_1, pos_quadruple, max_gap_num)) {continue;}
      if n == i && p == k {
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
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
          let ba_score = sta_fe_params.ba_score_mat[&pos_pair_1];
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair_2] + ba_score);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat[&pos_pair_2] + ba_score;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (m - 1, o - 1);
              if pos_quadruple != pseudo_pos_quadruple && (j - n - 1 + m - i - 1 > MAX_2_LOOP_LEN || l - p - 1 + o - k - 1 > MAX_2_LOOP_LEN) {continue;}
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple_2)) {continue;}
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
              let two_loop_fe_pair = if pos_quadruple == pseudo_pos_quadruple {
                (ml_tm_or_de_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.0, ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.1)
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
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair_2] + ba_score;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (m - 1, o - 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple_2)) {continue;}
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
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + if pos_quadruple == pseudo_pos_quadruple {0.} else {2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE}) + log_sta_pf_4_bpa;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + if pos_quadruple == pseudo_pos_quadruple {0.} else {2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE}) + log_sta_pf_4_bpa;
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
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair_2] + ba_score;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (m - 1, o - 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple_2)) {continue;}
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
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + if pos_quadruple == pseudo_pos_quadruple {0.} else {2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE}) + log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
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
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (n - 1, p);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          }
        }
        if p == k {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (n, p - 1);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
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
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          }
        }
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
} */

#[inline]
fn get_lbpap_mat(/* seq_pair: &SeqPair, */ seq_len_pair: &(usize, usize), /* sta_fe_params: &StaFeParams, */ log_sta_inside_ppf_mat_sets: &LogStaInsidePpfMatSets, log_sta_outside_ppf_4d_mats: &LogStaOutsidePpf4dMats, /* max_bp_span: usize, max_gap_num: usize */) -> LogProb4dMat {
  let mut lbpap_mat = LogProb4dMat::default();
  let log_pf = log_sta_inside_ppf_mat_sets.log_ppf_mats_4_external_loop.log_forward_ppf_mat[&(seq_len_pair.0 - 2, seq_len_pair.1 - 2)];
  for pos_quadruple in log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas.keys() {
    if !log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas.contains_key(pos_quadruple) {continue;}
    let lbpap = log_sta_inside_ppf_mat_sets.log_ppf_mat_4_bpas[pos_quadruple] + log_sta_outside_ppf_4d_mats.log_ppf_mat_4_bpas[pos_quadruple] - log_pf;
    if lbpap.is_finite() {
      // assert!(NEG_INFINITY < lbpap && lbpap <= 0.);
      lbpap_mat.insert(*pos_quadruple, lbpap);
    }
  }
  lbpap_mat
}

/* #[inline]
fn get_log_sta_backward_ppf_mat_sets(pos_quadruple: &PosQuadruple, pseudo_pos_quadruple: &PosQuadruple, seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, log_sta_ppf_mat_4_bpas: &LogPpf4dMat, max_gap_num: usize) -> LogStaPpfMatSets {
  let &(i, j, k, l) = pos_quadruple;
  let mut log_sta_ppf_mat_sets = LogStaPpfMatSets::new();
  for m in (i + 1 .. j + 1).rev() {
    for o in (k + 1 .. l + 1).rev() {
      let pos_pair_1 = (m, o);
      if !(is_min_gap_ok_1(&pos_pair_1, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_1, pos_quadruple, max_gap_num)) {continue;}
      if m == j && o == l {
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat.insert(pos_pair_1, 0.);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat.insert(pos_pair_1, NEG_INFINITY);
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
          let ba_score = sta_fe_params.ba_score_mat[&pos_pair_1];
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas.insert(pos_pair_1, log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair_2] + ba_score);
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat[&pos_pair_2] + ba_score;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (n + 1, p + 1);
              if pos_quadruple != pseudo_pos_quadruple && (j - n - 1 + m - i - 1 > MAX_2_LOOP_LEN || l - p - 1 + o - k - 1 > MAX_2_LOOP_LEN) {continue;}
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple_2)) {continue;}
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
              let two_loop_fe_pair = if pos_quadruple == pseudo_pos_quadruple {
                (ml_tm_or_de_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.0, ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.1)
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
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair_2] + ba_score;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (n + 1, p + 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple_2)) {continue;}
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
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + if pos_quadruple == pseudo_pos_quadruple {0.} else {2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE}) + log_sta_pf_4_bpa;
              if ep_of_term_4_log_pf.is_finite() {
                eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                  max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                }
              }
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + if pos_quadruple == pseudo_pos_quadruple {0.} else {2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE}) + log_sta_pf_4_bpa;
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
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat[&pos_pair_2] + ba_score;
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple_2 = (m, n, o, p);
              let pos_pair = (n + 1, p + 1);
              if !(is_min_gap_ok_1(&pos_pair, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num) && sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple_2)) {continue;}
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
              let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat[&pos_pair] - INVERSE_TEMPERATURE * (ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + if pos_quadruple == pseudo_pos_quadruple {0.} else {2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE}) + log_sta_ppf_mat_4_bpas[&pos_quadruple_2];
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
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (m + 1, o);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
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
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
          }
        }
        if o == l {
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let pos_pair_2 = (m, o + 1);
          if is_min_gap_ok_1(&pos_pair_2, pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, pos_quadruple, max_gap_num) {
            let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
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
            let mut max_ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
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
            log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
            log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
          }
        }
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_sas.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_2ls.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_1[&pos_pair_1],
          log_sta_ppf_mat_sets.log_ppf_mats_4_first_bpas_on_mls.log_ppf_mat_4_gaps_2[&pos_pair_1],
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
} */
 
/* #[inline]
// pub fn prob_cons_transformation_of_log_prob_mat_pair_on_sta(log_prob_mat_pairs_on_sta_with_rna_id_pairs: &LogProbMatPairsOnStaWithRnaIdPairs, rna_id_pair: &RnaIdPair, num_of_rnas: usize) -> LogProbMatPairOnSta {
pub fn prob_cons_transformation_of_lbpap_mat(lbpap_mats_with_rna_id_pairs: &LogProb4dMatsWithRnaIdPairs, rna_id_pair: &RnaIdPair, num_of_rnas: usize) -> LogProb4dMat {
  let mut lbpap_mat_pair_on_sta = log_prob_mat_pairs_on_sta_with_rna_id_pairs[rna_id_pair].clone();
  let log_coefficient = -((num_of_rnas - 1) as Prob).ln();
  // for (pos_quadruple, lbpap) in log_prob_mat_pair_on_sta.lbpap_mat.iter_mut() {
  for (pos_quadruple, lbpap) in lbpap_mats_with_rna_id_pairs.iter_mut() {
    let (i, j, k, l) = *pos_quadruple;
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = *lbpap;
    if max_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      // for (&(m, n, o, p), &lbpap_1) in log_prob_mat_pairs_on_sta_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].lbpap_mat.iter() {
      for (&(m, n, o, p), &lbpap_1) in lbpap_mats_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].iter() {
        if (rna_id_pair.0 < rna_id && m != i && n != j) || (rna_id_pair.0 > rna_id && o != i && p != j) {continue;}
        // for (&(m, n, o, p), &lbpap_2) in log_prob_mat_pairs_on_sta_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].lbpap_mat.iter() {
        for (&(m, n, o, p), &lbpap_2) in lbpap_mats_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].iter() {
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
  for (&(i, j), lbap) in log_prob_mat_pair_on_sta.lbap_mat.iter_mut() {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut max_ep_of_term_4_log_prob = *lbap;
    if max_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(max_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      for (&(k, l), lbap_1) in log_prob_mat_pairs_on_sta_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].lbap_mat.iter() {
        if (rna_id_pair.0 < rna_id && k != i) || (rna_id_pair.0 > rna_id && l != i) {continue;}
        for (&(k, l), lbap_2) in log_prob_mat_pairs_on_sta_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].lbap_mat.iter() {
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
  log_prob_mat_pair_on_sta
} */

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
    // let ref ref_2_log_prob_mat_pair_on_sta = log_prob_mat_pairs_on_sta_with_rna_id_pairs[&rna_id_pair];
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
      // let bpp = bpp_mat.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
      // if new_bpp > *bpp {
        // *bpp = new_bpp;
      // }
    }
  }
  // let mut upp_mat = upp_mat.clone();
  let mut upp_mat = vec![1.; upp_mat_len];
  for (i, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_poss.iter().enumerate() {
    if i == 0 || i == upp_mat_len - 1 {continue;}
    let max_ep_of_term_4_log_prob = max_eps_of_terms_4_log_probs_with_poss[i];
    if max_ep_of_term_4_log_prob.is_finite() {
      let new_upp = 1. - (log_coefficient + logsumexp(eps_of_terms_4_log_prob, max_ep_of_term_4_log_prob)).exp();
      upp_mat[i] = new_upp;
      // let ref mut up = upp_mat[i];
      // if new_up > *up {
        // *up = new_up;
      // }
    }
  }
  (bpp_mat, upp_mat)
}

#[inline]
pub fn get_bpap_mats_with_rna_id_pairs(lbpap_mats_with_rna_id_pairs: &LogProb4dMatsWithRnaIdPairs) -> Prob4dMatsWithRnaIdPairs {
  /* let mut prob_mat_pairs_on_sta_with_rna_id_pairs = ProbMatPairsOnStaWithRnaIdPairs::default();
  for (rna_id_pair, log_prob_mat_pair_on_sta) in log_prob_mat_pairs_on_sta_with_rna_id_pairs.iter() {
    prob_mat_pairs_on_sta_with_rna_id_pairs.insert(*rna_id_pair, get_prob_mat_pair_on_sta(log_prob_mat_pair_on_sta));
  }
  prob_mat_pairs_on_sta_with_rna_id_pairs */
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
pub fn remove_zeros_from_bpp_mat(bpp_mat: &ProbMat, seq_len: usize) -> SparseLogProbMat {
  let mut sparse_bpp_mat = SparseProbMat::default();
  for (i, bpps) in bpp_mat.iter().enumerate() {
    for (j, &bpp) in bpps.iter().enumerate() {
      if bpp > 0. {
        sparse_bpp_mat.insert((i + 1, j + 1), bpp);
      }
    }
  }
  // sparse_bpp_mat.insert((0, seq_len - 1), 1.);
  sparse_bpp_mat
}

#[inline]
pub fn remove_small_bpps_from_bpp_mat(sparse_bpp_mat: &SparseProbMat, min_bpp: Prob) -> SparseProbMat {
  sparse_bpp_mat.iter().filter(|(_, &bpp)| {bpp >= min_bpp}).map(|(pos_pair, &bpp)| {(*pos_pair, bpp)}).collect()
}

#[inline]
pub fn get_max_bp_span(sparse_bpp_mat: &SparseProbMat) -> usize {
  sparse_bpp_mat.iter().map(|(pos_pair, &bpp)| {pos_pair.1 - pos_pair.0 + 1}).max().unwrap()
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
fn is_min_gap_ok_1(pos_pair: &PosPair, /* pos_quadruple: &PosQuadruple, */ max_gap_num: usize) -> bool {
  /* let min_gap_num_1 = get_min_gap_num(&(pos_quadruple.0, pos_pair.0, pos_quadruple.2, pos_pair.1));
  let min_gap_num_2 = get_min_gap_num(&(pos_pair.0, pos_quadruple.1, pos_pair.1, pos_quadruple.3));
  if min_gap_num_1 + min_gap_num_2 <= max_gap_num {
    true
  } else {
    false
  } */
  if max(pos_pair.0, pos_pair.1) - min(pos_pair.0, pos_pair.1) <= max_gap_num {true} else {false}
}

#[inline]
fn is_min_gap_ok_2(pos_quadruple: &PosQuadruple, /* pseudo_pos_quadruple: &PosQuadruple, */ max_gap_num: usize) -> bool {
  /* let min_gap_num_1 = get_min_gap_num(&pos_quadruple);
  let min_gap_num_2 = get_min_gap_num(&(pseudo_pos_quadruple.0, pos_quadruple.0, pseudo_pos_quadruple.2, pos_quadruple.2));
  let min_gap_num_3 = get_min_gap_num(&(pos_quadruple.1, pseudo_pos_quadruple.1, pos_quadruple.3, pseudo_pos_quadruple.3));
  if min_gap_num_1 + min_gap_num_2 + min_gap_num_3 <= max_gap_num {
    true
  } else {
    false
  } */
  let pos_pair_1 = (pos_quadruple.0, pos_quadruple.2);
  let pos_pair_2 = (pos_quadruple.1, pos_quadruple.3);
  if pos_quadruple.1 - pos_quadruple.0 - (pos_quadruple.3 - pos_quadruple.2) <= max_gap_num  && is_min_gap_ok_1(&pos_pair_1, max_gap_num) && is_min_gap_ok_1(&pos_pair_2, max_gap_num) {true} else {false}
}

#[inline]
fn get_min_gap_num(pos_quadruple: &PosQuadruple) -> usize {
  let substr_len_pair = (pos_quadruple.1 - pos_quadruple.0, pos_quadruple.3 - pos_quadruple.2);
  get_seq_len_diff(&substr_len_pair)
}
