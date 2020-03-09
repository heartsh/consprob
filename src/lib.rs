extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate getopts;
extern crate scoped_threadpool;
extern crate itertools;

pub mod utils;

pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use self::scoped_threadpool::Pool;
pub use itertools::multizip;
pub use utils::*;

pub type Prob4dMat = FxHashMap<PosQuadruple, Prob>;
type PartFunc4dMat = FxHashMap<PosQuadruple, PartFunc>;
#[derive(Clone)]
pub struct StaInsidePartFuncMats {
  pub part_func_mat: SparsePartFuncMat,
  pub part_func_mat_4_bas: SparsePartFuncMat,
  pub part_func_mat_4_gaps_1: SparsePartFuncMat,
  pub part_func_mat_4_gaps_2: SparsePartFuncMat,
}
#[derive(Clone)]
pub struct TmpStaInsidePartFuncMatSets {
  pub part_func_mats_on_sa: StaInsidePartFuncMats,
  pub part_func_mats_4_internal_multiloop: StaInsidePartFuncMats,
  pub part_func_mats_4_first_bpas_on_mls: StaInsidePartFuncMats,
}
#[derive(Clone)]
pub struct StaInsidePartFuncMatSets {
  pub part_func_4d_mat_4_bpas: PartFunc4dMat,
  pub part_func_4d_mat_4_bpas_accessible_on_els: PartFunc4dMat,
  pub part_func_4d_mat_4_bpas_accessible_on_mls: PartFunc4dMat,
  pub forward_part_func_mats_4_external_loop: StaInsidePartFuncMats,
  pub backward_part_func_mats_4_external_loop: StaInsidePartFuncMats,
  pub forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples: TmpStaInsidePartFuncMatSetsWithPosQuadruples,
  pub backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples: TmpStaInsidePartFuncMatSetsWithPosQuadruples,
}
pub type TmpStaInsidePartFuncMatSetsWithPosQuadruples = FxHashMap<PosQuadruple, TmpStaInsidePartFuncMatSets>;
pub struct StaFeParams {
  pub ba_score_mat: SparseFreeEnergyMat,
  pub bpa_score_mat: FreeEnergy4dMat,
  pub exp_ba_score_mat: SparseFreeEnergyMat,
  pub exp_bpa_score_mat: FreeEnergy4dMat,
  pub opening_gap_penalty: FreeEnergy,
  pub extending_gap_penalty: FreeEnergy,
  pub exp_opening_gap_penalty: FreeEnergy,
  pub exp_extending_gap_penalty: FreeEnergy,
  scaler: FreeEnergy,
  pub is_logsumexp_applied: bool,
}
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type StaFeParamSetsWithRnaIdPairs = FxHashMap<RnaIdPair, StaFeParams>;
pub type Prob4dMatsWithRnaIdPairs = FxHashMap<RnaIdPair, Prob4dMat>;
pub type ProbMats = Vec<SparseProbMat>;
pub type Prob1dMats = Vec<Probs>;
pub type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
#[derive(Clone)]
pub struct FastaRecord {
  pub fasta_id: FastaId,
  pub seq: Seq,
}
pub type FastaRecords = Vec<FastaRecord>;
pub type SeqPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type FreeEnergyPair = (FreeEnergy, FreeEnergy);
pub type SparseFreeEnergyMat = FxHashMap<PosPair, FreeEnergy>;
pub type PosPairsWithPosPairs = FxHashMap<PosPair, PosPair>;
pub type BoolsWithPosPairs = FxHashMap<PosPair, bool>;
pub type ProbMatPair<'a> =  (&'a SparseProbMat, &'a SparseProbMat);
pub type SsFreeEnergyMatSetPair<'a> =  (&'a SsFreeEnergyMats, &'a SsFreeEnergyMats);
pub type NumOfThreads = u32;

impl StaFeParams {
  pub fn origin() -> StaFeParams {
    let sparse_free_energy_mat = SparseFreeEnergyMat::default();
    let free_energy_4d_mat = FreeEnergy4dMat::default();
    StaFeParams {
      ba_score_mat: sparse_free_energy_mat.clone(),
      bpa_score_mat: free_energy_4d_mat.clone(),
      exp_ba_score_mat: sparse_free_energy_mat,
      exp_bpa_score_mat: free_energy_4d_mat,
      opening_gap_penalty: 0.,
      extending_gap_penalty: 0.,
      exp_opening_gap_penalty: 0.,
      exp_extending_gap_penalty: 0.,
      scaler: NEG_INFINITY,
      is_logsumexp_applied: true,
    }
  }
  pub fn new(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), max_bp_span_pair: &PosPair, max_gap_num: Pos, max_gap_num_4_il: Pos, bpp_mat_pair: &ProbMatPair, max_free_energy: FreeEnergy, opening_gap_penalty: FreeEnergy, extending_gap_penalty: FreeEnergy, exp_opening_gap_penalty: FreeEnergy, exp_extending_gap_penalty: FreeEnergy) -> StaFeParams {
    let max = max(max_gap_num, max_gap_num_4_il);
    let mut sta_fe_params = StaFeParams::origin();
    sta_fe_params.is_logsumexp_applied = !(seq_len_pair.0 - 2 < MIN_SEQ_LENGTH_4_LOGSUMEXP_METHOD && seq_len_pair.1 - 2 < MIN_SEQ_LENGTH_4_LOGSUMEXP_METHOD);
    sta_fe_params.opening_gap_penalty = opening_gap_penalty;
    sta_fe_params.extending_gap_penalty = extending_gap_penalty;
    sta_fe_params.exp_opening_gap_penalty = exp_opening_gap_penalty;
    sta_fe_params.exp_extending_gap_penalty = exp_extending_gap_penalty;
    let seq_len_pair = (seq_len_pair.0 as Pos, seq_len_pair.1 as Pos);
    let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
    for i in 1 .. seq_len_pair.0 - 1 {
      let long_i = i as usize;
      let base = seq_pair.0[long_i];
      for j in 1 .. seq_len_pair.1 - 1 {
        let long_j = j as usize;
        let pos_pair = (i, j);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max) {continue;}
        let base_pair = (base, seq_pair.1[long_j]);
        sta_fe_params.ba_score_mat.insert(pos_pair, RIBOSUM_85_60_BA_SCORE_MAT[&base_pair]);
        sta_fe_params.exp_ba_score_mat.insert(pos_pair, EXP_RIBOSUM_85_60_BA_SCORE_MAT[&base_pair]);
      }
      let upper_j = if i + max_bp_span_pair.0 >= seq_len_pair.0 - 1 {seq_len_pair.0 - 1} else {i + max_bp_span_pair.0};
      for j in i + 1 .. upper_j {
        let long_j = j as usize;
        let pos_pair = (i, j);
        let base_pair = (base, seq_pair.0[long_j]);
        if !is_canonical(&base_pair) {continue;}
        if !bpp_mat_pair.0.contains_key(&pos_pair) {continue;}
        for k in 1 .. seq_len_pair.1 - 1 {
          if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
          let long_k = k as usize;
          let upper_l = if k + max_bp_span_pair.1 >= seq_len_pair.1 - 1 {seq_len_pair.1 - 1} else {k + max_bp_span_pair.1};
          for l in k + 1 .. upper_l {
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_pair_2 = (k, l);
            let pos_quadruple = (pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1);
            let long_l = l as usize;
            let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
            if !is_canonical(&base_pair_2) {continue;}
            if !(is_min_gap_ok_2(&pos_quadruple, max_gap_num_4_il) && bpp_mat_pair.1.contains_key(&pos_pair_2)) {continue;}
            let base_quadruple = (base_pair, base_pair_2);
            sta_fe_params.bpa_score_mat.insert(pos_quadruple, RIBOSUM_85_60_BPA_SCORE_MAT[&base_quadruple]);
            sta_fe_params.exp_bpa_score_mat.insert(pos_quadruple, EXP_RIBOSUM_85_60_BPA_SCORE_MAT[&base_quadruple]);
          }
        }
      }
    }
    sta_fe_params.scaler = if sta_fe_params.is_logsumexp_applied {
      - max_free_energy
    } else {
      1. / max_free_energy.exp()
    };
    sta_fe_params
  }
}

impl StaInsidePartFuncMats {
  pub fn new() -> StaInsidePartFuncMats {
    let part_func_mat = SparsePartFuncMat::default();
    StaInsidePartFuncMats {
      part_func_mat: part_func_mat.clone(),
      part_func_mat_4_bas: part_func_mat.clone(),
      part_func_mat_4_gaps_1: part_func_mat.clone(),
      part_func_mat_4_gaps_2: part_func_mat,
    }
  }
}

impl TmpStaInsidePartFuncMatSets {
  pub fn new() -> TmpStaInsidePartFuncMatSets {
    let part_func_mats = StaInsidePartFuncMats::new();
    TmpStaInsidePartFuncMatSets {
      part_func_mats_on_sa: part_func_mats.clone(),
      part_func_mats_4_internal_multiloop: part_func_mats.clone(),
      part_func_mats_4_first_bpas_on_mls: part_func_mats,
    }
  }
}

impl StaInsidePartFuncMatSets {
  pub fn new() -> StaInsidePartFuncMatSets {
    let part_func_4d_mat = PartFunc4dMat::default();
    let sta_inside_part_func_mats = StaInsidePartFuncMats::new();
    let tmp_sta_inside_part_func_mat_sets_with_pos_quadruples = TmpStaInsidePartFuncMatSetsWithPosQuadruples::default();
    StaInsidePartFuncMatSets {
      part_func_4d_mat_4_bpas: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_els: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_mls: part_func_4d_mat,
      forward_part_func_mats_4_external_loop: sta_inside_part_func_mats.clone(),
      backward_part_func_mats_4_external_loop: sta_inside_part_func_mats,
      forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples: tmp_sta_inside_part_func_mat_sets_with_pos_quadruples.clone(),
      backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples: tmp_sta_inside_part_func_mat_sets_with_pos_quadruples,
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

pub const MIN_SEQ_LENGTH_4_LOGSUMEXP_METHOD: usize = 200;
pub const MAX_GAP_NUM_4_IL: Pos = 100;
pub const MIN_GAP_NUM_4_IL: Pos = 5;
pub const DEFAULT_OPENING_GAP_PENALTY: FreeEnergy = 0.;
pub const DEFAULT_EXTENDING_GAP_PENALTY: FreeEnergy = 0.;
pub const DEFAULT_MIN_BPP: Prob = 0.005;
pub const DEFAULT_OFFSET_4_MAX_GAP_NUM: Pos = 0;

pub fn io_algo_4_bpap_mat(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &PosPair, max_gap_num: Pos, max_gap_num_4_il: Pos, bpp_mat_pair: &ProbMatPair, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair) -> Prob4dMat {
  let sta_inside_part_func_mat_sets = get_sta_inside_part_func_mat_sets(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, max_gap_num_4_il, bpp_mat_pair, ss_free_energy_mat_set_pair);
  let sta_outside_part_func_4d_mat_4_bpas = get_sta_outside_part_func_4d_mat_4_bpas(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num_4_il, &sta_inside_part_func_mat_sets, bpp_mat_pair, ss_free_energy_mat_set_pair);
  get_bpap_mat(seq_len_pair, &sta_inside_part_func_mat_sets, &sta_outside_part_func_4d_mat_4_bpas, sta_fe_params)
}

pub fn get_sta_inside_part_func_mat_sets(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &PosPair, max_gap_num: Pos, max_gap_num_4_il: Pos, bpp_mat_pair: &ProbMatPair, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair) -> StaInsidePartFuncMatSets {
  let seq_len_pair = (seq_len_pair.0 as Pos, seq_len_pair.1 as Pos);
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut sta_inside_part_func_mat_sets = StaInsidePartFuncMatSets::new();
  let scaler = sta_fe_params.scaler;
  if !sta_fe_params.is_logsumexp_applied {
    for substr_len_1 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.0 + 1 {
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        let (long_i, long_j) = (i as usize, j as usize);
        let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
        if !is_canonical(&base_pair) {continue;}
        if !bpp_mat_pair.0.contains_key(&(i, j)) {continue;}
        let invert_base_pair = invert_bp(&base_pair);
        let invert_stacking_base_pair = invert_bp(&(seq_pair.0[long_i + 1], seq_pair.0[long_j - 1]));
        let exp_ml_tm_delta_fe = EXP_ML_TM_DELTA_FES[invert_base_pair.0][invert_base_pair.1][invert_stacking_base_pair.0][invert_stacking_base_pair.1];
        let stacking_bp = (seq_pair.0[long_i - 1], seq_pair.0[long_j + 1]);
        let exp_ml_tm_or_de_delta_fe = if i > 1 && j < seq_len_pair.0 - 2 {
          EXP_ML_TM_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.0][stacking_bp.1]
        } else if i > 1 {
          EXP_FIVE_PRIME_DE_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.0]
        } else if j < seq_len_pair.0 - 2 {
          EXP_THREE_PRIME_DE_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.1]
        } else {
          1.
        };
        let exp_au_or_gu_end_penalty_delta_fe = if is_au_or_gu(&base_pair) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
        let exp_hl_fe = ss_free_energy_mat_set_pair.0.exp_hl_fe_mat[&(i, j)];
        for substr_len_2 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.1 + 1 {
          let min_gap_num_4_il = max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2);
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          for k in 1 .. seq_len_pair.1 - substr_len_2 {
            if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let l = k + substr_len_2 - 1;
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let (long_k, long_l) = (k as usize, l as usize);
            let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
            if !is_canonical(&base_pair_2) {continue;}
            let pos_quadruple = (i, j, k, l);
            match sta_fe_params.exp_bpa_score_mat.get(&pos_quadruple) {
              Some(exp_bpa_score) => {
                let forward_tmp_sta_inside_part_func_mat_sets = get_forward_tmp_sta_inside_part_func_mat_sets(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_inside_part_func_mat_sets, bpp_mat_pair);
                let backward_tmp_sta_inside_part_func_mat_sets = get_backward_tmp_sta_inside_part_func_mat_sets(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_inside_part_func_mat_sets, bpp_mat_pair);
                let mut sum = 0.;
                let pos_pair = (j - 1, l - 1);
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair) {
                  Some(&part_func) => {
                    sum += exp_bpa_score * exp_hl_fe * ss_free_energy_mat_set_pair.1.exp_hl_fe_mat[&(k, l)] * part_func;
                  }, None => {},
                }
                for m in i + 1 .. j {
                  let long_m = m as usize;
                  for n in m + 1 .. j {
                    let long_n = n as usize;
                    let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    if !is_canonical(&base_pair_3) {continue;}
                    if !bpp_mat_pair.0.contains_key(&(m, n)) {continue;}
                    if long_m - long_i - 1 + long_j - long_n - 1 > MAX_2_LOOP_LEN {continue;}
                    let exp_2loop_fe = ss_free_energy_mat_set_pair.0.exp_2loop_fe_4d_mat[&(i, j, m, n)];
                    for o in k + 1 .. l {
                      if !is_min_gap_ok_1(&(m, o), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                      let long_o = o as usize;
                      for p in o + 1 .. l {
                        if !is_min_gap_ok_1(&(n, p), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                        let long_p = p as usize;
                        let base_pair_4 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                        if !is_canonical(&base_pair_4) {continue;}
                        if long_o - long_k - 1 + long_l - long_p - 1 > MAX_2_LOOP_LEN {continue;}
                        let pos_quadruple_2 = (m, n, o, p);
                        match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                          Some(&part_func) => {
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(m - 1, o - 1)) {
                              Some(&part_func_2) => {
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(n + 1, p + 1)) {
                                  Some(&part_func_3) => {
                                    let exp_2loop_fe_2 = ss_free_energy_mat_set_pair.1.exp_2loop_fe_4d_mat[&(k, l, o, p)];
                                    sum += exp_bpa_score * part_func_2 * exp_2loop_fe * exp_2loop_fe_2 * part_func / scaler * part_func_3 / scaler;
                                  }, None => {},
                                }
                              }, None => {},
                            }
                          }, None => {},
                        }
                      }
                    }
                  }
                }
                let exp_au_or_gu_end_penalty_delta_fe_2 = if is_au_or_gu(&base_pair_2) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair) {
                  Some(&part_func) => {
                    let invert_base_pair_2 = invert_bp(&base_pair_2);
                    let invert_stacking_base_pair_2 = invert_bp(&(seq_pair.1[long_k + 1], seq_pair.1[long_l - 1]));
                    let exp_ml_tm_delta_fe_2 = EXP_ML_TM_DELTA_FES[invert_base_pair_2.0][invert_base_pair_2.1][invert_stacking_base_pair_2.0][invert_stacking_base_pair_2.1];
                    sum += exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe * exp_ml_tm_delta_fe_2 * exp_au_or_gu_end_penalty_delta_fe * exp_au_or_gu_end_penalty_delta_fe_2 * part_func;
                  }, None => {},
                }
                if sum > 0. {
                  sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
                  let stacking_bp_2 = (seq_pair.1[long_k - 1], seq_pair.1[long_l + 1]);
                  let exp_ml_tm_or_de_delta_fe_2 = if k > 1 && l < seq_len_pair.1 - 2 {
                    EXP_ML_TM_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.0][stacking_bp_2.1]
                  } else if k > 1 {
                    EXP_FIVE_PRIME_DE_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.0]
                  } else if l < seq_len_pair.1 - 2 {
                    EXP_THREE_PRIME_DE_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.1]
                  } else {
                    1.
                  };
                  sum *= exp_ml_tm_or_de_delta_fe * exp_ml_tm_or_de_delta_fe_2 * exp_au_or_gu_end_penalty_delta_fe * exp_au_or_gu_end_penalty_delta_fe_2;
                  sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.insert(pos_quadruple, sum);
                  sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, sum * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE));
                }
                sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples.insert(pos_quadruple, forward_tmp_sta_inside_part_func_mat_sets);
                sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples.insert(pos_quadruple, backward_tmp_sta_inside_part_func_mat_sets);
              }, None => {},
            }
          }
        }
      }
    }
    let leftmost_pos_pair = (0, 0);
    let rightmost_pos_pair = (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
    sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.insert(leftmost_pos_pair, scaler);
    sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(leftmost_pos_pair, scaler);
    sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.insert(rightmost_pos_pair, scaler);
    sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(rightmost_pos_pair, scaler);
    for i in 0 .. seq_len_pair.0 - 1 {
      for j in 0 .. seq_len_pair.1 - 1 {
        let pos_pair = (i, j);
        if pos_pair == (0, 0) {continue;}
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
        let mut sum = 0.;
        for k in 1 .. i {
          if !bpp_mat_pair.0.contains_key(&(k, i)) {continue;}
          for l in 1 .. j {
            let pos_pair_2 = (k - 1, l - 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num) {continue;}
            let pos_quadruple = (k, i, l, j);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        if i > 0 && j > 0 {
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.get(&(i - 1, j - 1)) {
            Some(&part_func) => {
              let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
              sum += part_func * exp_ba_score;
            }, None => {},
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        sum = 0.;
        if i > 0 {
          let pos_pair_2 = (i - 1, j);
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_extending_gap_penalty;
            }, None => {},
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        sum = 0.;
        if j > 0 {
          let pos_pair_2 = (i, j - 1);
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_extending_gap_penalty;
            }, None => {},
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
    for i in (1 .. seq_len_pair.0).rev() {
      for j in (1 .. seq_len_pair.1).rev() {
        let pos_pair = (i, j);
        if pos_pair == (seq_len_pair.0 - 1, seq_len_pair.1 - 1) {continue;}
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
        let mut sum = 0.;
        for k in i + 1 .. seq_len_pair.0 - 1 {
          if !bpp_mat_pair.0.contains_key(&(i, k)) {continue;}
          for l in j + 1 .. seq_len_pair.1 - 1 {
            let pos_pair_2 = (k + 1, l + 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num) {continue;}
            let pos_quadruple = (i, k, j, l);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 {
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.get(&(i + 1, j + 1)) {
            Some(&part_func) => {
              let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
              sum += part_func * exp_ba_score;
            }, None => {},
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        sum = 0.;
        if i < seq_len_pair.0 - 1 {
          let pos_pair_2 = (i + 1, j);
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_extending_gap_penalty;
            }, None => {},
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        sum = 0.;
        if j < seq_len_pair.1 - 1 {
          let pos_pair_2 = (i, j + 1);
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_opening_gap_penalty;
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              sum += part_func * sta_fe_params.exp_extending_gap_penalty;
            }, None => {},
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None  => {},
        }
        match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
  } else {
    for substr_len_1 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.0 + 1 {
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        let (long_i, long_j) = (i as usize, j as usize);
        let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
        if !is_canonical(&base_pair) {continue;}
        if !bpp_mat_pair.0.contains_key(&(i, j)) {continue;}
        let invert_base_pair = invert_bp(&base_pair);
        let invert_stacking_base_pair = invert_bp(&(seq_pair.0[long_i + 1], seq_pair.0[long_j - 1]));
        let ml_tm_delta_fe = ML_TM_DELTA_FES[invert_base_pair.0][invert_base_pair.1][invert_stacking_base_pair.0][invert_stacking_base_pair.1];
        let stacking_bp = (seq_pair.0[long_i - 1], seq_pair.0[long_j + 1]);
        let ml_tm_or_de_delta_fe = if i > 1 && j < seq_len_pair.0 - 2 {
          ML_TM_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.0][stacking_bp.1]
        } else if i > 1 {
          FIVE_PRIME_DE_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.0]
        } else if j < seq_len_pair.0 - 2 {
          THREE_PRIME_DE_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.1]
        } else {
          0.
        };
        let au_or_gu_end_penalty_delta_fe = if is_au_or_gu(&base_pair) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
        let hl_fe = ss_free_energy_mat_set_pair.0.hl_fe_mat[&(i, j)];
        for substr_len_2 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.1 + 1 {
          let min_gap_num_4_il = max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2);
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          for k in 1 .. seq_len_pair.1 - substr_len_2 {
            if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let l = k + substr_len_2 - 1;
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let (long_k, long_l) = (k as usize, l as usize);
            let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
            if !is_canonical(&base_pair_2) {continue;}
            let pos_quadruple = (i, j, k, l);
            match sta_fe_params.bpa_score_mat.get(&pos_quadruple) {
              Some(&bpa_score) => {
                let forward_tmp_sta_inside_part_func_mat_sets = get_forward_tmp_sta_inside_part_func_mat_sets(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_inside_part_func_mat_sets, bpp_mat_pair);
                let backward_tmp_sta_inside_part_func_mat_sets = get_backward_tmp_sta_inside_part_func_mat_sets(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_inside_part_func_mat_sets, bpp_mat_pair);
                let mut sum = NEG_INFINITY;
                let pos_pair = (j - 1, l - 1);
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair) {
                  Some(&part_func) => {
                    logsumexp(&mut sum, bpa_score + hl_fe + ss_free_energy_mat_set_pair.1.hl_fe_mat[&(k, l)] + part_func);
                  }, None => {},
                }
                for m in i + 1 .. j {
                  let long_m = m as usize;
                  for n in m + 1 .. j {
                    let long_n = n as usize;
                    let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    if !is_canonical(&base_pair_3) {continue;}
                    if !bpp_mat_pair.0.contains_key(&(m, n)) {continue;}
                    if long_m - long_i - 1 + long_j - long_n - 1 > MAX_2_LOOP_LEN {continue;}
                    let twoloop_fe = ss_free_energy_mat_set_pair.0.twoloop_fe_4d_mat[&(i, j, m, n)];
                    for o in k + 1 .. l {
                      if !is_min_gap_ok_1(&(m, o), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                      let long_o = o as usize;
                      for p in o + 1 .. l {
                         if !is_min_gap_ok_1(&(n, p), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                        let long_p = p as usize;
                        let base_pair_4 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                        if !is_canonical(&base_pair_4) {continue;}
                        if long_o - long_k - 1 + long_l - long_p - 1 > MAX_2_LOOP_LEN {continue;}
                        let pos_quadruple_2 = (m, n, o, p);
                        match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                          Some(&part_func) => {
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(m - 1, o - 1)) {
                              Some(&part_func_2) => {
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(n + 1, p + 1)) {
                                  Some(&part_func_3) => {
                                    let twoloop_fe_2 = ss_free_energy_mat_set_pair.1.twoloop_fe_4d_mat[&(k, l, o, p)];
                                    logsumexp(&mut sum, bpa_score + part_func_2 + twoloop_fe + twoloop_fe_2 + part_func - scaler + part_func_3 - scaler);
                                  }, None => {},
                                }
                              }, None => {},
                            }
                          }, None => {},
                        }
                      }
                    }
                  }
                }
                let au_or_gu_end_penalty_delta_fe_2 = if is_au_or_gu(&base_pair_2) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair) {
                  Some(&part_func) => {
                    let invert_base_pair_2 = invert_bp(&base_pair_2);
                    let invert_stacking_base_pair_2 = invert_bp(&(seq_pair.1[long_k + 1], seq_pair.1[long_l - 1]));
                    let ml_tm_delta_fe_2 = ML_TM_DELTA_FES[invert_base_pair_2.0][invert_base_pair_2.1][invert_stacking_base_pair_2.0][invert_stacking_base_pair_2.1];
                    logsumexp(&mut sum, bpa_score + 2. * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe + ml_tm_delta_fe_2 + au_or_gu_end_penalty_delta_fe + au_or_gu_end_penalty_delta_fe_2 + part_func);
                  }, None => {},
                }
                if sum > NEG_INFINITY {
                  sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
                  let stacking_bp_2 = (seq_pair.1[long_k - 1], seq_pair.1[long_l + 1]);
                  let ml_tm_or_de_delta_fe_2 = if k > 1 && l < seq_len_pair.1 - 2 {
                    ML_TM_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.0][stacking_bp_2.1]
                  } else if k > 1 {
                    FIVE_PRIME_DE_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.0]
                  } else if l < seq_len_pair.1 - 2 {
                    THREE_PRIME_DE_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.1]
                  } else {
                    0.
                  };
                  sum += ml_tm_or_de_delta_fe + ml_tm_or_de_delta_fe_2 + au_or_gu_end_penalty_delta_fe + au_or_gu_end_penalty_delta_fe_2;
                  sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.insert(pos_quadruple, sum);
                  sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, sum + 2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
                }
                sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples.insert(pos_quadruple, forward_tmp_sta_inside_part_func_mat_sets);
                sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples.insert(pos_quadruple, backward_tmp_sta_inside_part_func_mat_sets);
              }, None => {},
            }
          }
        }
      }
    }
    let leftmost_pos_pair = (0, 0);
    let rightmost_pos_pair = (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
    sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.insert(leftmost_pos_pair, scaler);
    sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(leftmost_pos_pair, scaler);
    sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.insert(rightmost_pos_pair, scaler);
    sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(rightmost_pos_pair, scaler);
    for i in 0 .. seq_len_pair.0 - 1 {
      for j in 0 .. seq_len_pair.1 - 1 {
        let pos_pair = (i, j);
        if pos_pair == (0, 0) {continue;}
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
        let mut sum = NEG_INFINITY;
        for k in 1 .. i {
          if !bpp_mat_pair.0.contains_key(&(k, i)) {continue;}
          for l in 1 .. j {
            let pos_pair_2 = (k - 1, l - 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num) {continue;}
            let pos_quadruple = (k, i, l, j);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        if i > 0 && j > 0 {
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.get(&(i - 1, j - 1)) {
            Some(&part_func) => {
              let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
              logsumexp(&mut sum, part_func + ba_score);
            }, None => {},
          }
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        if i > 0 {
          let pos_pair_2 = (i - 1, j);
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
            }, None => {},
          }
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        if j > 0 {
          let pos_pair_2 = (i, j - 1);
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
            }, None => {},
          }
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
    for i in (1 .. seq_len_pair.0).rev() {
      for j in (1 .. seq_len_pair.1).rev() {
        let pos_pair = (i, j);
        if pos_pair == (seq_len_pair.0 - 1, seq_len_pair.1 - 1) {continue;}
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
        let mut sum = NEG_INFINITY;
        for k in i + 1 .. seq_len_pair.0 - 1 {
          if !bpp_mat_pair.0.contains_key(&(i, k)) {continue;}
          for l in j + 1 .. seq_len_pair.1 - 1 {
            let pos_pair_2 = (k + 1, l + 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num) {continue;}
            let pos_quadruple = (i, k, j, l);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 {
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.get(&(i + 1, j + 1)) {
            Some(&part_func) => {
              let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
              logsumexp(&mut sum, part_func + ba_score);
            }, None => {},
          }
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        if i < seq_len_pair.0 - 1 {
          let pos_pair_2 = (i + 1, j);
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
            }, None => {},
          }
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        if j < seq_len_pair.1 - 1 {
          let pos_pair_2 = (i, j + 1);
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
            }, None => {},
          }
          match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
            Some(&part_func) => {
              logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
            }, None => {},
          }
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
  }
  sta_inside_part_func_mat_sets
}

pub fn get_forward_tmp_sta_inside_part_func_mat_sets(seq_len_pair: &PosPair, sta_fe_params: &StaFeParams, max_gap_num_4_il: Pos, pos_quadruple: &PosQuadruple, sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets, bpp_mat_pair: &ProbMatPair) -> TmpStaInsidePartFuncMatSets {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut forward_tmp_sta_inside_part_func_mat_sets = TmpStaInsidePartFuncMatSets::new();
  let &(i, j, k, l) = pos_quadruple;
  let scaler = sta_fe_params.scaler;
  if !sta_fe_params.is_logsumexp_applied {
    for u in i .. j {
      for v in k .. l {
        let pos_pair = (u, v);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
        if u == i && v == k {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair, scaler);
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, scaler);
          continue;
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
        let pos_pair_2 = (u - 1, v - 1);
        let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair,
              part_func * exp_ba_score
            );
          }, None => {},
        }
        let mut sum = 0.;
        let pos_pair_2 = (u - 1, v);
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
        sum = 0.;
        for m in i + 1 .. u {
          if !bpp_mat_pair.0.contains_key(&(m, u)) {continue;}
          for n in k + 1 .. v {
            let pos_pair_2 = (m - 1, n - 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (m, u, n, v);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u - 1, v - 1);
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * exp_ba_score;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u - 1, v);
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
        sum = 0.;
        for m in i + 1 .. u {
          if !bpp_mat_pair.0.contains_key(&(m, u)) {continue;}
          for n in k + 1 .. v {
            let pos_pair_2 = (m - 1, n - 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (m, u, n, v);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u - 1, v - 1);
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * exp_ba_score;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u - 1, v);
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
  } else {
    for u in i .. j {
      for v in k .. l {
        let pos_pair = (u, v);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
        if u == i && v == k {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair, scaler);
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, scaler);
          continue;
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
        let pos_pair_2 = (u - 1, v - 1);
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair,
              part_func + ba_score
            );
          }, None => {},
        }
        let mut sum = NEG_INFINITY;
        let pos_pair_2 = (u - 1, v);
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
        sum = NEG_INFINITY;
        for m in i + 1 .. u {
          if !bpp_mat_pair.0.contains_key(&(m, u)) {continue;}
          for n in k + 1 .. v {
            let pos_pair_2 = (m - 1, n - 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (m, u, n, v);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u - 1, v - 1);
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + ba_score);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u - 1, v);
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
        sum = NEG_INFINITY;
        for m in i + 1 .. u {
          if !bpp_mat_pair.0.contains_key(&(m, u)) {continue;}
          for n in k + 1 .. v {
            let pos_pair_2 = (m - 1, n - 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (m, u, n, v);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u - 1, v - 1);
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + ba_score);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u - 1, v);
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
  }
  forward_tmp_sta_inside_part_func_mat_sets
}

pub fn get_backward_tmp_sta_inside_part_func_mat_sets(seq_len_pair: &PosPair, sta_fe_params: &StaFeParams, max_gap_num_4_il: Pos, pos_quadruple: &PosQuadruple, sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets, bpp_mat_pair: &ProbMatPair) -> TmpStaInsidePartFuncMatSets {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut backward_tmp_sta_inside_part_func_mat_sets = TmpStaInsidePartFuncMatSets::new();
  let &(i, j, k, l) = pos_quadruple;
  let scaler = sta_fe_params.scaler;
  if !sta_fe_params.is_logsumexp_applied {
    for u in (i + 1 .. j + 1).rev() {
      for v in (k + 1 .. l + 1).rev() {
        let pos_pair = (u, v);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
        if u == j && v == l {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair, scaler);
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, scaler);
          continue;
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
        let pos_pair_2 = (u + 1, v + 1);
        let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair,
              part_func * exp_ba_score
            );
          }, None => {},
        }
        let mut sum = 0.;
        let pos_pair_2 = (u + 1, v);
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
        sum = 0.;
        for m in u + 1 .. j {
          if !bpp_mat_pair.0.contains_key(&(u, m)) {continue;}
          for n in v + 1 .. l {
            let pos_pair_2 = (m + 1, n + 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (u, m, v, n);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u + 1, v + 1);
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * exp_ba_score;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u + 1, v);
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
        sum = 0.;
        for m in u + 1 .. j {
          if !bpp_mat_pair.0.contains_key(&(u, m)) {continue;}
          for n in v + 1 .. l {
            let pos_pair_2 = (m + 1, n + 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (u, m, v, n);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    sum += part_func_2 * part_func / scaler;
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u + 1, v + 1);
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * exp_ba_score;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u + 1, v);
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_opening_gap_penalty;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            sum += part_func * sta_fe_params.exp_extending_gap_penalty;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = 0.;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            sum += part_func;
          }, None => {},
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
  } else {
    for u in (i + 1 .. j + 1).rev() {
      for v in (k + 1 .. l + 1).rev() {
        let pos_pair = (u, v);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
        if u == j && v == l {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair, scaler);
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, scaler);
          continue;
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
        let pos_pair_2 = (u + 1, v + 1);
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair,
              part_func + ba_score
            );
          }, None => {},
        }
        let mut sum = NEG_INFINITY;
        let pos_pair_2 = (u + 1, v);
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
        sum = NEG_INFINITY;
        for m in u + 1 .. j {
          if !bpp_mat_pair.0.contains_key(&(u, m)) {continue;}
          for n in v + 1 .. l {
            let pos_pair_2 = (m + 1, n + 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (u, m, v, n);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u + 1, v + 1);
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + ba_score);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u + 1, v);
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.insert(pos_pair, sum);
        }
        // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
        sum = NEG_INFINITY;
        for m in u + 1 .. j {
          if !bpp_mat_pair.0.contains_key(&(u, m)) {continue;}
          for n in v + 1 .. l {
            let pos_pair_2 = (m + 1, n + 1);
            if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_quadruple_2 = (u, m, v, n);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    logsumexp(&mut sum, part_func_2 + part_func - scaler);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
        let pos_pair_2 = (u + 1, v + 1);
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + ba_score);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u + 1, v);
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.opening_gap_penalty);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair_2) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func + sta_fe_params.extending_gap_penalty);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
        sum = NEG_INFINITY;
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.get(&pos_pair) {
          Some(&part_func) => {
            logsumexp(&mut sum, part_func);
          }, None => {},
        }
        if sum > NEG_INFINITY {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.insert(pos_pair, sum);
        }
      }
    }
  }
  backward_tmp_sta_inside_part_func_mat_sets
}

pub fn get_sta_outside_part_func_4d_mat_4_bpas(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &PosPair, max_gap_num_4_il: Pos, sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets, bpp_mat_pair: &ProbMatPair, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair) -> PartFunc4dMat {
  let seq_len_pair = (seq_len_pair.0 as Pos, seq_len_pair.1 as Pos);
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut sta_outside_part_func_4d_mat_4_bpas = PartFunc4dMat::default();
  let scaler = sta_fe_params.scaler;
  if !sta_fe_params.is_logsumexp_applied {
    for substr_len_1 in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.0 + 1).rev() {
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        let (long_i, long_j) = (i as usize, j as usize);
        let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
        if !is_canonical(&base_pair) {continue;}
        if !bpp_mat_pair.0.contains_key(&(i, j)) {continue;}
        for substr_len_2 in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.1 + 1).rev() {
          let min_gap_num_4_il = max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2);
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          for k in 1 .. seq_len_pair.1 - substr_len_2 {
            let l = k + substr_len_2 - 1;
            if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let (long_k, long_l) = (k as usize, l as usize);
            let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
            if !is_canonical(&base_pair_2) {continue;}
            let pos_quadruple = (i, j, k, l);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.get(&pos_quadruple) {
              Some(&part_func_4_bpa) => {
                let mut sum = 0.;
                match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.get(&(i - 1, k - 1)) {
                  Some(&part_func) => {
                    match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.get(&(j + 1, l + 1)) {
                      Some(&part_func_2) => {
                        sum += part_func * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els[&pos_quadruple] / part_func_4_bpa * part_func_2 / scaler;
                      }, None => {},
                    }
                  }, None => {},
                }
                for m in 1 .. i {
                  let long_m = m as usize;
                  for n in j + 1 .. seq_len_pair.0 - 1 {
                    let long_n = n as usize;
                    let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    if !is_canonical(&base_pair_3) {continue;}
                    if !bpp_mat_pair.0.contains_key(&(m, n)) {continue;}
                    if long_n - long_j - 1 + long_i - long_m - 1 > MAX_2_LOOP_LEN {continue;}
                    // let exp_2loop_fe = get_exp_2_loop_fe(seq_pair.0, &(long_m, long_n), &(long_i, long_j));
                    let exp_2loop_fe = ss_free_energy_mat_set_pair.0.exp_2loop_fe_4d_mat[&(m, n, i, j)];
                    for o in 1 .. k {
                      if !is_min_gap_ok_1(&(m, o), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                      let long_o = o as usize;
                      for p in l + 1 .. seq_len_pair.1 - 1 {
                        if !is_min_gap_ok_1(&(n, p), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                        let long_p = p as usize;
                        let base_pair_4 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                        if !is_canonical(&base_pair_4) {continue;}
                        if long_p - long_l - 1 + long_k - long_o - 1 > MAX_2_LOOP_LEN {continue;}
                        let pos_quadruple_2 = (m, n, o, p);
                        match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                          Some(&part_func) => {
                            let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            let ref backward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func_2) => {
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_3) => {
                                    let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple_2];
                                    // let exp_2loop_fe_2 = get_exp_2_loop_fe(seq_pair.1, &(long_o, long_p), &(long_k, long_l));
                                    let exp_2loop_fe_2 = ss_free_energy_mat_set_pair.1.exp_2loop_fe_4d_mat[&(o, p, k, l)];
                                    sum += part_func_2 * exp_bpa_score * exp_2loop_fe * exp_2loop_fe_2 * part_func / scaler * part_func_3 / scaler;
                                  }, None => {},
                                }
                              }, None => {},
                            }
                          }, None => {},
                        }
                      }
                    }
                  }
                }
                let part_func_ratio = sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple] / part_func_4_bpa;
                for m in 1 .. i {
                  let long_m = m as usize;
                  for n in j + 1 .. seq_len_pair.0 - 1 {
                    let long_n = n as usize;
                    let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    if !is_canonical(&base_pair_3) {continue;}
                    if !bpp_mat_pair.0.contains_key(&(m, n)) {continue;}
                    let base_pair = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    let invert_base_pair = invert_bp(&base_pair);
                    let invert_stacking_bp = invert_bp(&(seq_pair.0[long_m + 1], seq_pair.0[long_n - 1]));
                    let exp_ml_tm_delta_fe = EXP_ML_TM_DELTA_FES[invert_base_pair.0][invert_base_pair.1][invert_stacking_bp.0][invert_stacking_bp.1];
                    let exp_au_or_gu_end_penalty_delta_fe = if is_au_or_gu(&base_pair) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
                    for o in 1 .. k {
                      if !is_min_gap_ok_1(&(m, o), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                      let long_o = o as usize;
                      for p in l + 1 .. seq_len_pair.1 - 1 {
                        if !is_min_gap_ok_1(&(n, p), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                        let long_p = p as usize;
                        let base_pair_4 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                        if !is_canonical(&base_pair_4) {continue;}
                        let pos_quadruple_2 = (m, n, o, p);
                        match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                          Some(&part_func_4_bpa_2) => {
                            let base_pair_2 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                            let invert_base_pair_2 = invert_bp(&base_pair_2);
                            let invert_stacking_bp_2 = invert_bp(&(seq_pair.1[long_o + 1], seq_pair.1[long_p - 1]));
                            let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            let ref backward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple_2];
                            let exp_ml_tm_delta_fe_2 = EXP_ML_TM_DELTA_FES[invert_base_pair_2.0][invert_base_pair_2.1][invert_stacking_bp_2.0][invert_stacking_bp_2.1];
                            let exp_au_or_gu_end_penalty_delta_fe_2 = if is_au_or_gu(&base_pair_2) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.};
                            let coefficient = part_func_ratio * exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe * exp_ml_tm_delta_fe_2 * exp_au_or_gu_end_penalty_delta_fe * exp_au_or_gu_end_penalty_delta_fe_2 * part_func_4_bpa_2;
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func) => {
                                let coefficient = coefficient * part_func / scaler;
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                              }, None => {},
                            }
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func) => {
                                let coefficient = coefficient * part_func / scaler;
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                              }, None => {},
                            }
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func) => {
                                let coefficient = coefficient * part_func / scaler;
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    sum += coefficient * part_func_2 / scaler;
                                  }, None => {},
                                }
                              }, None => {},
                            }
                          }, None => {},
                        }
                      }
                    }
                  }
                }
                if sum > 0. {
                  sta_outside_part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
                }
              }, None => {},
            }
          }
        }
      }
    }
  } else {
    for substr_len_1 in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.0 + 1).rev() {
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        let (long_i, long_j) = (i as usize, j as usize);
        let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
        if !is_canonical(&base_pair) {continue;}
        if !bpp_mat_pair.0.contains_key(&(i, j)) {continue;}
        for substr_len_2 in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.1 + 1).rev() {
          let min_gap_num_4_il = max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2);
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          for k in 1 .. seq_len_pair.1 - substr_len_2 {
            let l = k + substr_len_2 - 1;
            if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let (long_k, long_l) = (k as usize, l as usize);
            let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
            if !is_canonical(&base_pair_2) {continue;}
            let pos_quadruple = (i, j, k, l);
            match sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.get(&pos_quadruple) {
              Some(&part_func_4_bpa) => {
                let mut sum = NEG_INFINITY;
                match sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.get(&(i - 1, k - 1)) {
                  Some(&part_func) => {
                    match sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.get(&(j + 1, l + 1)) {
                      Some(&part_func_2) => {
                        logsumexp(&mut sum, part_func + sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els[&pos_quadruple] - part_func_4_bpa + part_func_2 - scaler);
                      }, None => {},
                    }
                  }, None => {},
                }
                for m in 1 .. i {
                  let long_m = m as usize;
                  for n in j + 1 .. seq_len_pair.0 - 1 {
                    let long_n = n as usize;
                    let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    if !is_canonical(&base_pair_3) {continue;}
                    if !bpp_mat_pair.0.contains_key(&(m, n)) {continue;}
                    if long_n - long_j - 1 + long_i - long_m - 1 > MAX_2_LOOP_LEN {continue;}
                    // let twoloop_fe = get_2_loop_fe(seq_pair.0, &(long_m, long_n), &(long_i, long_j));
                    let twoloop_fe = ss_free_energy_mat_set_pair.0.twoloop_fe_4d_mat[&(m, n, i, j)];
                    for o in 1 .. k {
                      if !is_min_gap_ok_1(&(m, o), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                      let long_o = o as usize;
                      for p in l + 1 .. seq_len_pair.1 - 1 {
                        if !is_min_gap_ok_1(&(n, p), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                        let long_p = p as usize;
                        let base_pair_4 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                        if !is_canonical(&base_pair_4) {continue;}
                        if long_p - long_l - 1 + long_k - long_o - 1 > MAX_2_LOOP_LEN {continue;}
                        let pos_quadruple_2 = (m, n, o, p);
                        match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                          Some(&part_func) => {
                            let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            let ref backward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func_2) => {
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_3) => {
                                    let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                                    // let twoloop_fe_2 = get_2_loop_fe(seq_pair.1, &(long_o, long_p), &(long_k, long_l));
                                    let twoloop_fe_2 = ss_free_energy_mat_set_pair.1.twoloop_fe_4d_mat[&(o, p, k, l)];
                                    logsumexp(&mut sum, part_func_2 + bpa_score + twoloop_fe + twoloop_fe_2 + part_func - scaler + part_func_3 - scaler);
                                  }, None => {},
                                }
                              }, None => {},
                            }
                          }, None => {},
                        }
                      }
                    }
                  }
                }
                let part_func_ratio = sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple] - part_func_4_bpa;
                for m in 1 .. i {
                  let long_m = m as usize;
                  for n in j + 1 .. seq_len_pair.0 - 1 {
                    let long_n = n as usize;
                    let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    if !is_canonical(&base_pair_3) {continue;}
                    if !bpp_mat_pair.0.contains_key(&(m, n)) {continue;}
                    let base_pair = (seq_pair.0[long_m], seq_pair.0[long_n]);
                    let invert_base_pair = invert_bp(&base_pair);
                    let invert_stacking_bp = invert_bp(&(seq_pair.0[long_m + 1], seq_pair.0[long_n - 1]));
                    let ml_tm_delta_fe = ML_TM_DELTA_FES[invert_base_pair.0][invert_base_pair.1][invert_stacking_bp.0][invert_stacking_bp.1];
                    let au_or_gu_end_penalty_delta_fe = if is_au_or_gu(&base_pair) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
                    for o in 1 .. k {
                      if !is_min_gap_ok_1(&(m, o), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                      let long_o = o as usize;
                      for p in l + 1 .. seq_len_pair.1 - 1 {
                        if !is_min_gap_ok_1(&(n, p), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                        let long_p = p as usize;
                        let base_pair_4 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                        if !is_canonical(&base_pair_4) {continue;}
                        let pos_quadruple_2 = (m, n, o, p);
                        match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                          Some(&part_func_4_bpa_2) => {
                            let base_pair_2 = (seq_pair.1[long_o], seq_pair.1[long_p]);
                            let invert_base_pair_2 = invert_bp(&base_pair_2);
                            let invert_stacking_bp_2 = invert_bp(&(seq_pair.1[long_o + 1], seq_pair.1[long_p - 1]));
                            let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            let ref backward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                            let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                            let ml_tm_delta_fe_2 = ML_TM_DELTA_FES[invert_base_pair_2.0][invert_base_pair_2.1][invert_stacking_bp_2.0][invert_stacking_bp_2.1];
                            let au_or_gu_end_penalty_delta_fe_2 = if is_au_or_gu(&base_pair_2) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.};
                            let coefficient = part_func_ratio + bpa_score + 2. * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe + ml_tm_delta_fe_2 + au_or_gu_end_penalty_delta_fe + au_or_gu_end_penalty_delta_fe_2 + part_func_4_bpa_2;
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func) => {
                                let coefficient = coefficient + part_func - scaler;
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                              }, None => {},
                            }
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func) => {
                                let coefficient = coefficient + part_func - scaler;
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                              }, None => {},
                            }
                            match forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.get(&(i - 1, k - 1)) {
                              Some(&part_func) => {
                                let coefficient = coefficient + part_func - scaler;
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                                match backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.get(&(j + 1, l + 1)) {
                                  Some(&part_func_2) => {
                                    logsumexp(&mut sum, coefficient + part_func_2 - scaler);
                                  }, None => {},
                                }
                              }, None => {},
                            }
                          }, None => {},
                        }
                      }
                    }
                  }
                }
                if sum > NEG_INFINITY {
                  sta_outside_part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
                }
              }, None => {},
            }
          }
        }
      }
    }
  }
  sta_outside_part_func_4d_mat_4_bpas
}

fn get_bpap_mat(seq_len_pair: &(usize, usize), sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets, sta_outside_part_func_4d_mat_4_bpas: &PartFunc4dMat, sta_fe_params: &StaFeParams) -> Prob4dMat {
  let mut bpap_mat = Prob4dMat::default();
  let scaler = sta_fe_params.scaler;
  let part_func = sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat[&(seq_len_pair.0 as Pos - 2, seq_len_pair.1 as Pos - 2)];
  assert!(part_func.is_finite());
  for pos_quadruple in sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.keys() {
    match sta_outside_part_func_4d_mat_4_bpas.get(pos_quadruple) {
      Some(&outside_part_func) => {
        let bpap = if !sta_fe_params.is_logsumexp_applied {
          sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas[pos_quadruple] / part_func * outside_part_func / scaler
        } else {
          (sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas[pos_quadruple] - part_func + outside_part_func - scaler).exp()
        };
        debug_assert!(0. <= bpap && bpap <= 1.);
        bpap_mat.insert(*pos_quadruple, bpap);
      }, None => {},
    }
  }
  bpap_mat
}
 
pub fn pct_of_bpp_and_upp_mat(bpap_mats_with_rna_id_pairs: &Prob4dMatsWithRnaIdPairs, rna_id: RnaId, num_of_rnas: usize, bpp_mat: &SparseProbMat, upp_mat_len: usize) -> (SparseProbMat, Probs) {
  let weight = 1. / (num_of_rnas - 1) as Prob;
  let mut new_bpp_mat = SparseProbMat::default();
  let mut new_upp_mat = vec![1.; upp_mat_len];
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id == rna_id_2 {continue;}
    let rna_id_pair = if rna_id < rna_id_2 {(rna_id, rna_id_2)} else {(rna_id_2, rna_id)};
    let ref ref_2_bpap_mat = bpap_mats_with_rna_id_pairs[&rna_id_pair];
    for (pos_quadruple, &bpap) in ref_2_bpap_mat.iter() {
      let pos_pair = if rna_id < rna_id_2 {(pos_quadruple.0, pos_quadruple.1)} else {(pos_quadruple.2, pos_quadruple.3)};
      if !new_bpp_mat.contains_key(&pos_pair) {
        new_bpp_mat.insert(pos_pair, 0.);
      }
      let weighted_bpap = weight * bpap;
      *new_bpp_mat.get_mut(&pos_pair).unwrap() += weighted_bpap;
      new_upp_mat[pos_pair.0 as usize] -= weighted_bpap;
      new_upp_mat[pos_pair.1 as usize] -= weighted_bpap;
    }
  }
  for (&(i, j), &bpp) in bpp_mat {
    let pos_pair = (i + 1, j + 1);
    if !new_bpp_mat.contains_key(&pos_pair) {
      new_bpp_mat.insert(pos_pair, bpp);
    }
  }
  (new_bpp_mat, new_upp_mat)
}

pub fn remove_small_bpps_from_bpp_mat(sparse_bpp_mat: &SparseProbMat, min_bpp: Prob) -> SparseProbMat {
  sparse_bpp_mat.iter().filter(|(_, &bpp)| {bpp >= min_bpp}).map(|(&(i, j), &bpp)| {((i + 1, j + 1), bpp)}).collect()
}

pub fn get_max_bp_span(sparse_bpp_mat: &SparseProbMat) -> Pos {
  sparse_bpp_mat.iter().map(|(pos_pair, _)| {pos_pair.1 - pos_pair.0 + 1}).max().unwrap() as Pos
}

pub fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
}

pub fn get_seq_len_diff(pos_quadruple: &PosQuadruple) -> Pos {
  let seq_len_pair = (pos_quadruple.1 + 1 - pos_quadruple.0, pos_quadruple.3 + 1 - pos_quadruple.2);
  max(seq_len_pair.0, seq_len_pair.1) - min(seq_len_pair.0, seq_len_pair.1)
}

fn is_min_gap_ok_1(pos_pair: &PosPair, pos_quadruple: &PosQuadruple, max_gap_num: Pos) -> bool {
  let min_gap_num_1 = get_seq_len_diff(&(pos_quadruple.0, pos_pair.0, pos_quadruple.2, pos_pair.1));
  let min_gap_num_2 = get_seq_len_diff(&(pos_pair.0, pos_quadruple.1, pos_pair.1, pos_quadruple.3));
  if min_gap_num_1 <= max_gap_num && min_gap_num_2 <= max_gap_num {
    true
  } else {
    false
  }
}

fn is_min_gap_ok_2(pos_quadruple: &PosQuadruple, max_gap_num: Pos) -> bool {
  let min_gap_num = get_seq_len_diff(&pos_quadruple);
  if min_gap_num <= max_gap_num {
    true
  } else {
    false
  }
}

pub fn phyloprob(thread_pool: &mut Pool, fasta_records: &FastaRecords, opening_gap_penalty: FreeEnergy, extending_gap_penalty: FreeEnergy, min_bpp: Prob, offset_4_max_gap_num: Pos) -> (ProbMats, Prob1dMats) {
  let exp_opening_gap_penalty = opening_gap_penalty.exp();
  let exp_extending_gap_penalty = extending_gap_penalty.exp();
  let num_of_fasta_records = fasta_records.len();
  let mut bpp_mats = vec![SparseProbMat::default(); num_of_fasta_records];
  let mut sparse_bpp_mats = bpp_mats.clone();
  let mut upp_mats = vec![Probs::new(); num_of_fasta_records];
  let mut max_bp_spans = vec![0; num_of_fasta_records];
  let mut max_free_energies = vec![NEG_INFINITY; num_of_fasta_records];
  let mut ss_free_energy_mat_sets = vec![SsFreeEnergyMats::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (bpp_mat, sparse_bpp_mat, upp_mat, max_bp_span, fasta_record, max_free_energy, ss_free_energy_mats) in multizip((bpp_mats.iter_mut(), sparse_bpp_mats.iter_mut(), upp_mats.iter_mut(), max_bp_spans.iter_mut(), fasta_records.iter(), max_free_energies.iter_mut(), ss_free_energy_mat_sets.iter_mut())) {
      let seq_len = fasta_record.seq.len();
      scope.execute(move || {
        let (obtained_bpp_mat, obtained_upp_mat, obtained_max_free_energy, obtained_ss_free_energy_mats) = get_bpp_and_unpair_prob_mats(&fasta_record.seq[1 .. seq_len - 1]);
        *bpp_mat = obtained_bpp_mat;
        *ss_free_energy_mats = obtained_ss_free_energy_mats;
        ss_free_energy_mats.sparsify(bpp_mat, min_bpp);
        *sparse_bpp_mat = remove_small_bpps_from_bpp_mat(&bpp_mat, min_bpp);
        *upp_mat = obtained_upp_mat;
        *max_free_energy = obtained_max_free_energy;
        *max_bp_span = get_max_bp_span(sparse_bpp_mat);
        upp_mat.insert(0, 1.);
        upp_mat.push(1.);
      });
    }
  });
  let mut sta_fe_param_sets_with_rna_id_pairs = StaFeParamSetsWithRnaIdPairs::default();
  let mut bpap_mats_with_rna_id_pairs = Prob4dMatsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      sta_fe_param_sets_with_rna_id_pairs.insert(rna_id_pair, StaFeParams::origin());
      bpap_mats_with_rna_id_pairs.insert(rna_id_pair, Prob4dMat::default());
    }
  }
  thread_pool.scoped(|scope| {
    for (rna_id_pair, bpap_mat) in bpap_mats_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
      let max_gap_num = offset_4_max_gap_num + (max(seq_len_pair.0, seq_len_pair.1) - min(seq_len_pair.0, seq_len_pair.1)) as Pos;
      let max_gap_num_4_il = max(min(max_gap_num, MAX_GAP_NUM_4_IL), MIN_GAP_NUM_4_IL);
      let max_bp_span_pair = (max_bp_spans[rna_id_pair.0], max_bp_spans[rna_id_pair.1]);
      let bpp_mat_pair = (&sparse_bpp_mats[rna_id_pair.0], &sparse_bpp_mats[rna_id_pair.1]);
      let max_free_energy = max_free_energies[rna_id_pair.0] + max_free_energies[rna_id_pair.1];
      let ss_free_energy_mat_set_pair = (&ss_free_energy_mat_sets[rna_id_pair.0], &ss_free_energy_mat_sets[rna_id_pair.1]);
      scope.execute(move || {
        let sta_fe_params = StaFeParams::new(&seq_pair, &seq_len_pair, &max_bp_span_pair, max_gap_num, max_gap_num_4_il, &bpp_mat_pair, max_free_energy, opening_gap_penalty, extending_gap_penalty, exp_opening_gap_penalty, exp_extending_gap_penalty);
        *bpap_mat = io_algo_4_bpap_mat(&seq_pair, &seq_len_pair, &sta_fe_params, &max_bp_span_pair, max_gap_num, max_gap_num_4_il, &bpp_mat_pair, &ss_free_energy_mat_set_pair);
      });
    }
  });
  thread_pool.scoped(|scope| {
    for (rna_id, bpp_mat, upp_mat) in multizip((0 .. num_of_fasta_records, bpp_mats.iter_mut(), upp_mats.iter_mut())) {
      let ref ref_2_bpap_mats_with_rna_id_pairs = bpap_mats_with_rna_id_pairs;
      scope.execute(move || {
        let prob_mat_pair = pct_of_bpp_and_upp_mat(ref_2_bpap_mats_with_rna_id_pairs, rna_id, num_of_fasta_records, bpp_mat, upp_mat.len());
        *bpp_mat = prob_mat_pair.0;
        *upp_mat = prob_mat_pair.1;
      });
    }
  });
  (bpp_mats, upp_mats)
}
