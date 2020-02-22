extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate getopts;

pub mod utils;

pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use utils::*;

pub type PosQuadruple = (Pos, Pos, Pos, Pos);
pub type Prob4dMat = HashMap<PosQuadruple, Prob, Hasher>;
type PartFunc4dMat = HashMap<PosQuadruple, PartFunc, Hasher>;
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
  pub part_func_mats_4_internal_2loop: StaInsidePartFuncMats,
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
pub type TmpStaInsidePartFuncMatSetsWithPosQuadruples = HashMap<PosQuadruple, TmpStaInsidePartFuncMatSets, Hasher>;
pub struct StaFeParams {
  pub ba_score_mat: SparseFreeEnergyMat,
  pub bpa_score_mat: FreeEnergy4dMat,
  pub exp_ba_score_mat: SparseFreeEnergyMat,
  pub exp_bpa_score_mat: FreeEnergy4dMat,
  pub opening_gap_penalty: FreeEnergy,
  pub extending_gap_penalty: FreeEnergy,
  pub exp_opening_gap_penalty: FreeEnergy,
  pub exp_extending_gap_penalty: FreeEnergy,
  pub invert_exp_max_free_energy: FreeEnergy,
}
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type StaFeParamSetsWithRnaIdPairs = HashMap<RnaIdPair, StaFeParams, Hasher>;
pub type Prob4dMatsWithRnaIdPairs = HashMap<RnaIdPair, Prob4dMat, Hasher>;
pub type ProbMats = Vec<SparseProbMat>;
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
pub type FreeEnergyPair = (FreeEnergy, FreeEnergy);
pub type SparseFreeEnergyMat = HashMap<PosPair, FreeEnergy, Hasher>;
pub type FreeEnergy4dMat = HashMap<PosQuadruple, FreeEnergy, Hasher>;
pub type PosPairsWithPosPairs = HashMap<PosPair, PosPair, Hasher>;
pub type BoolsWithPosPairs = HashMap<PosPair, bool, Hasher>;
pub struct MaxStaFreeEnergies {
  pub free_energy_on_sa: FreeEnergy,
  pub free_energy_4_internal_2loop: FreeEnergy,
  pub free_energy_4_internal_multiloop: FreeEnergy,
}
#[derive(Clone)]
pub struct MaxStaFreeEnergyMatSets {
  pub free_energy_4d_mat_4_bpas: FreeEnergy4dMat,
  pub free_energy_4d_mat_4_bpas_accessible_on_els: FreeEnergy4dMat,
  pub free_energy_4d_mat_4_bpas_accessible_on_mls: FreeEnergy4dMat,
  pub free_energy_mats_4_external_loop: MaxStaFreeEnergyMats,
}
#[derive(Clone)]
pub struct MaxStaFreeEnergyMats {
  pub free_energy_mat: SparseFreeEnergyMat,
  pub free_energy_mat_4_bas: SparseFreeEnergyMat,
  pub free_energy_mat_4_gaps_1: SparseFreeEnergyMat,
  pub free_energy_mat_4_gaps_2: SparseFreeEnergyMat,
}
#[derive(Clone)]
pub struct TmpMaxStaFreeEnergyMatSets {
  pub free_energy_mats_on_sa: MaxStaFreeEnergyMats,
  pub free_energy_mats_4_internal_2loop: MaxStaFreeEnergyMats,
  pub free_energy_mats_4_internal_multiloop: MaxStaFreeEnergyMats,
  pub free_energy_mats_4_first_bpas_on_mls: MaxStaFreeEnergyMats,
}

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
      invert_exp_max_free_energy: 0.,
    }
  }
  pub fn new(rna_id_pair: &RnaIdPair, fasta_records: &FastaRecords, max_bp_span_pair: &(usize, usize), max_gap_num: usize, max_gap_num_4_il: usize, bpp_mats: &ProbMats, opening_gap_penalty: FreeEnergy, extending_gap_penalty: FreeEnergy, exp_opening_gap_penalty: FreeEnergy, exp_extending_gap_penalty: FreeEnergy) -> StaFeParams {
    let max = max(max(max_gap_num, max_gap_num_4_il), MAX_GAP_NUM_4_IL);
    let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
    let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
    let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
    let mut sta_fe_params = StaFeParams::origin();
    sta_fe_params.opening_gap_penalty = opening_gap_penalty;
    sta_fe_params.extending_gap_penalty = extending_gap_penalty;
    sta_fe_params.exp_opening_gap_penalty = exp_opening_gap_penalty;
    sta_fe_params.exp_extending_gap_penalty = exp_extending_gap_penalty;
    let bpp_mat_pair = (&bpp_mats[rna_id_pair.0], &bpp_mats[rna_id_pair.1]);
    for i in 1 .. seq_len_pair.0 - 1 {
      let base = seq_pair.0[i];
      for j in 1 .. seq_len_pair.1 - 1 {
        let pos_pair = (i, j);
        /* let min_gap_num = get_seq_len_diff(&(0, i, 0, j)) + get_seq_len_diff(&(i, seq_len_pair.0, j, seq_len_pair.1));
        if min_gap_num > max_gap_num {continue;} */
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max) {continue;}
        let base_pair = (base, seq_pair.1[j]);
        sta_fe_params.ba_score_mat.insert(pos_pair, RIBOSUM_85_60_BA_SCORE_MAT[&base_pair]);
        sta_fe_params.exp_ba_score_mat.insert(pos_pair, EXP_RIBOSUM_85_60_BA_SCORE_MAT[&base_pair]);
      }
      let upper_j = if i + max_bp_span_pair.0 >= seq_len_pair.0 - 1 {seq_len_pair.0 - 1} else {i + max_bp_span_pair.0};
      for j in i + 1 .. upper_j {
        let pos_pair = (i, j);
        let base_pair = (base, seq_pair.0[j]);
        if !bpp_mat_pair.0.contains_key(&pos_pair) {continue;}
        for k in 1 .. seq_len_pair.1 - 1 {
          if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
          let upper_l = if k + max_bp_span_pair.1 >= seq_len_pair.1 - 1 {seq_len_pair.1 - 1} else {k + max_bp_span_pair.1};
          for l in k + 1 .. upper_l {
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_pair_2 = (k, l);
            let pos_quadruple = (pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1);
            // if !(is_min_gap_ok_2(&pos_quadruple, max_gap_num_4_il) && bpp_mat_pair.1.contains_key(&pos_pair_2)) {continue;}
            if !(is_min_gap_ok_2(&pos_quadruple, max_gap_num_4_il) && bpp_mat_pair.1.contains_key(&pos_pair_2)) {continue;}
            /* let min_gap_num_4_il = get_seq_len_diff(&pos_quadruple);
            if min_gap_num_4_il > max_gap_num_4_il {continue;} */
            /* let min_gap_num = get_seq_len_diff(&(0, j, 0, l)) + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
            if min_gap_num > max_gap_num {continue;}
            if !bpp_mat_pair.1.contains_key(&pos_pair_2) {continue;} */
            let base_quadruple = (base_pair, (seq_pair.1[k], seq_pair.1[l]));
            sta_fe_params.bpa_score_mat.insert(pos_quadruple, RIBOSUM_85_60_BPA_SCORE_MAT[&base_quadruple]);
            sta_fe_params.exp_bpa_score_mat.insert(pos_quadruple, EXP_RIBOSUM_85_60_BPA_SCORE_MAT[&base_quadruple]);
          }
        }
      }
    }
    let max_sta_free_energy = get_max_sta_free_energy(&seq_pair, &seq_len_pair, &sta_fe_params, &max_bp_span_pair, max_gap_num, max_gap_num_4_il);
    // println!("Max energy computed.");
    sta_fe_params.invert_exp_max_free_energy = 1. / max_sta_free_energy.exp();
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
      part_func_mats_4_internal_2loop: part_func_mats.clone(),
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

impl TmpMaxStaFreeEnergyMatSets {
  pub fn new() -> TmpMaxStaFreeEnergyMatSets {
    let free_energy_mats = MaxStaFreeEnergyMats::new();
    TmpMaxStaFreeEnergyMatSets {
      free_energy_mats_on_sa: free_energy_mats.clone(),
      free_energy_mats_4_internal_2loop: free_energy_mats.clone(),
      free_energy_mats_4_internal_multiloop: free_energy_mats.clone(),
      free_energy_mats_4_first_bpas_on_mls: free_energy_mats.clone(),
    }
  }
}

impl MaxStaFreeEnergyMatSets {
  pub fn new() -> MaxStaFreeEnergyMatSets {
    let free_energy_4d_mat = FreeEnergy4dMat::default();
    MaxStaFreeEnergyMatSets {
      free_energy_4d_mat_4_bpas: free_energy_4d_mat.clone(),
      free_energy_4d_mat_4_bpas_accessible_on_els: free_energy_4d_mat.clone(),
      free_energy_4d_mat_4_bpas_accessible_on_mls: free_energy_4d_mat,
      free_energy_mats_4_external_loop: MaxStaFreeEnergyMats::new(),
    }
  }
}

impl MaxStaFreeEnergyMats {
  pub fn new() -> MaxStaFreeEnergyMats {
    let free_energy_mat = SparseFreeEnergyMat::default();
    MaxStaFreeEnergyMats {
      free_energy_mat: free_energy_mat.clone(),
      free_energy_mat_4_bas: free_energy_mat.clone(),
      free_energy_mat_4_gaps_1: free_energy_mat.clone(),
      free_energy_mat_4_gaps_2: free_energy_mat,
    }
  }
}

pub const MAX_GAP_NUM_4_IL: usize = 5;

#[inline]
pub fn io_algo_4_bpap_mat(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize, max_gap_num_4_il: usize) -> Prob4dMat {
  let (sta_inside_part_func_mat_sets, farthest_pos_pairs_with_pos_pairs, bools_with_pair_aligned_right_pos_pairs) = get_sta_inside_part_func_mat_sets(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, max_gap_num_4_il);
  // println!("Inside part funcs computed.");
  let sta_outside_part_func_4d_mat_4_bpas = get_sta_outside_part_func_4d_mat_4_bpas(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, max_gap_num_4_il, &sta_inside_part_func_mat_sets, &farthest_pos_pairs_with_pos_pairs, &bools_with_pair_aligned_right_pos_pairs);
  get_bpap_mat(seq_len_pair, &sta_inside_part_func_mat_sets, &sta_outside_part_func_4d_mat_4_bpas, sta_fe_params)
}

#[inline]
pub fn get_sta_inside_part_func_mat_sets(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize, max_gap_num_4_il: usize) -> (StaInsidePartFuncMatSets, PosPairsWithPosPairs, BoolsWithPosPairs) {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut sta_inside_part_func_mat_sets = StaInsidePartFuncMatSets::new();
  let invert_exp_max_free_energy = sta_fe_params.invert_exp_max_free_energy;
  let mut farthest_pos_pairs_with_pos_pairs = PosPairsWithPosPairs::default();
  let mut bools_with_pair_aligned_right_pos_pairs = BoolsWithPosPairs::default();
  for substr_len_1 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. max_bp_span_pair.0 + 1 {
    for substr_len_2 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. max_bp_span_pair.1 + 1 {
      let min_gap_num_4_il = max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2);
      // if min_gap_num_4_il > max_gap_num_4_il {continue;}
      if min_gap_num_4_il > max_gap_num_4_il {continue;}
      // if min_gap_num_4_il > max_gap_num {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !sta_fe_params.exp_bpa_score_mat.contains_key(&pos_quadruple) {continue;}
          // println!("Forward computing for {:?}.", &pos_quadruple);
          let forward_tmp_sta_inside_part_func_mat_sets = get_forward_tmp_sta_inside_part_func_mat_sets(seq_pair, seq_len_pair, sta_fe_params, max_gap_num, max_gap_num_4_il, &pos_quadruple, invert_exp_max_free_energy, &sta_inside_part_func_mat_sets);
          // println!("Forward computed for {:?}.", &pos_quadruple);
          let mut sum = 0.;
          let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple];
          let pos_pair = (j - 1, l - 1);
          if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_pair) {
            sum += exp_bpa_score * get_exp_hl_fe(seq_pair.0, &(i, j)) * get_exp_hl_fe(seq_pair.1, &(k, l)) * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_pair];
          }
          if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat.contains_key(&pos_pair) {
            sum += exp_bpa_score * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat[&pos_pair];
          }
          if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&pos_pair) {
            let bp_closing_loop_pair = (
              (seq_pair.0[i], seq_pair.0[j]),
              (seq_pair.1[k], seq_pair.1[l]),
            );
            let exp_ml_tm_delta_fe_pair = (
              EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[j - 1])))],
              EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[l - 1])))],
            );
            let exp_au_or_gu_end_penalty_delta_fe_pair = (
              if is_au_or_gu(&bp_closing_loop_pair.0) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
              if is_au_or_gu(&bp_closing_loop_pair.1) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
            );
            sum += exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe_pair.0 * exp_ml_tm_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1 * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&pos_pair];
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
            if !bools_with_pair_aligned_right_pos_pairs.contains_key(&(j, l)) {
              bools_with_pair_aligned_right_pos_pairs.insert((j, l), true);
            }
            let accessible_bp_pair = (
              (seq_pair.0[i], seq_pair.0[j]),
              (seq_pair.1[k], seq_pair.1[l]),
            );
            let stacking_bp_pair = (
              (seq_pair.0[i - 1], seq_pair.0[j + 1]),
              (seq_pair.1[k - 1], seq_pair.1[l + 1]),
            );
            let exp_ml_tm_or_de_delta_fe_pair = (
              if i > 1 && j < seq_len_pair.0 - 2 {
                EXP_ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
              } else if i > 1 {
                EXP_FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
              } else if j < seq_len_pair.0 - 2 {
                EXP_THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
              } else {
                1.
              },
              if k > 1 && l < seq_len_pair.1 - 2 {
                EXP_ML_TM_DELTA_FES[&(accessible_bp_pair.1, stacking_bp_pair.1)]
              } else if k > 1 {
                EXP_FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).0)]
              } else if l < seq_len_pair.1 - 2 {
                EXP_THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.1, (stacking_bp_pair.1).1)]
              } else {
                1.
              },
            );
            let exp_au_or_gu_end_penalty_delta_fe_pair = (
              if is_au_or_gu(&accessible_bp_pair.0) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
              if is_au_or_gu(&accessible_bp_pair.1) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
            );
            sum *= exp_ml_tm_or_de_delta_fe_pair.0 * exp_ml_tm_or_de_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1;
            sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.insert(pos_quadruple, sum);
            sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, sum * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE));
          }
          let backward_tmp_sta_inside_part_func_mat_sets = get_backward_tmp_sta_inside_part_func_mat_sets(seq_len_pair, sta_fe_params, max_gap_num, max_gap_num_4_il, &pos_quadruple, invert_exp_max_free_energy, &sta_inside_part_func_mat_sets);
          sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples.insert(pos_quadruple, forward_tmp_sta_inside_part_func_mat_sets);
          sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples.insert(pos_quadruple, backward_tmp_sta_inside_part_func_mat_sets);
          if !farthest_pos_pairs_with_pos_pairs.contains_key(&(i, k)) {
            farthest_pos_pairs_with_pos_pairs.insert((i, k), (j, l));
          } else {
            *farthest_pos_pairs_with_pos_pairs.get_mut(&(i, k)).expect("Failed to get an element from a hash map.") = (j, l);
          }
        }
      }
    }
  }
  let leftmost_pos_pair = (0, 0);
  let rightmost_pos_pair = (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
  sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.insert(leftmost_pos_pair, invert_exp_max_free_energy);
  sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(leftmost_pos_pair, invert_exp_max_free_energy);
  sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.insert(rightmost_pos_pair, invert_exp_max_free_energy);
  sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(rightmost_pos_pair, invert_exp_max_free_energy);
  for i in 0 .. seq_len_pair.0 - 1 {
    for j in 0 .. seq_len_pair.1 - 1 {
      let pos_pair = (i, j);
      if pos_pair == (0, 0) {continue;}
      /* let min_gap_num = get_seq_len_diff(&(0, i, 0, j)) + get_seq_len_diff(&(i, seq_len_pair.0, j, seq_len_pair.1));
      if min_gap_num > max_gap_num {continue;} */
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
      let mut sum = 0.;
      for k in 1 .. i {
        for l in 1 .. j {
          let pos_pair_2 = (k - 1, l - 1);
          let pos_quadruple = (k, i, l, j);
          /* let min_gap_num = get_seq_len_diff(&(0, k, 0, l)) + get_seq_len_diff(&pos_quadruple) + get_seq_len_diff(&(i, seq_len_pair.0, j, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !(sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.contains_key(&pos_pair_2) && sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.contains_key(&pos_quadruple)) {
            continue;
          }
          sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els[&pos_quadruple] / invert_exp_max_free_energy;
        }
      }
      if i > 0 && j > 0 && sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.contains_key(&(i - 1, j - 1)) {
        let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
        sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat[&(i - 1, j - 1)] * exp_ba_score;
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(pos_pair, sum);
      }
      /* let min_gap_num_2 = min_gap_num /* + 2 */;
      if min_gap_num_2 <= max_gap_num { */
        sum = 0.;
        if i > 0 {
          let pos_pair_2 = (i - 1, j);
          if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        sum = 0.;
        if j > 0 {
          let pos_pair_2 = (i, j - 1);
          if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_bas[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2[&pos_pair];
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
      /* let min_gap_num = get_seq_len_diff(&(0, i, 0, j)) + get_seq_len_diff(&(i, seq_len_pair.0, j, seq_len_pair.1));
      if min_gap_num > max_gap_num {continue;} */
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
      let mut sum = 0.;
      for k in i + 1 .. seq_len_pair.0 - 1 {
        for l in j + 1 .. seq_len_pair.1 - 1 {
          let pos_pair_2 = (k + 1, l + 1);
          let pos_quadruple = (i, k, j, l);
          /* let min_gap_num = get_seq_len_diff(&(0, i, 0, j)) + get_seq_len_diff(&pos_quadruple) + get_seq_len_diff(&(k, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !(sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.contains_key(&pos_pair_2) && sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els.contains_key(&pos_quadruple)) {
            continue;
          }
          sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els[&pos_quadruple] / invert_exp_max_free_energy;
        }
      }
      if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 && sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.contains_key(&(i + 1, j + 1)) {
        let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
        sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat[&(i + 1, j + 1)] * exp_ba_score;
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.insert(pos_pair, sum);
      }
      /* let min_gap_num_2 = min_gap_num /* + 2 */;
      if min_gap_num_2 <= max_gap_num { */
        sum = 0.;
        if i < seq_len_pair.0 - 1 {
          let pos_pair_2 = (i + 1, j);
          if sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        sum = 0.;
        if j < seq_len_pair.1 - 1 {
          let pos_pair_2 = (i, j + 1);
          if sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
            sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
          }
        }
        if sum > 0. {
          sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_bas[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_1[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.insert(pos_pair, sum);
      }
    }
  }
  (sta_inside_part_func_mat_sets, farthest_pos_pairs_with_pos_pairs, bools_with_pair_aligned_right_pos_pairs)
}

#[inline]
pub fn get_max_sta_free_energy(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize, max_gap_num_4_il: usize) -> FreeEnergy {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut max_sta_free_energy_mat_sets = MaxStaFreeEnergyMatSets::new();
  for substr_len_1 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. max_bp_span_pair.0 + 1 {
    for substr_len_2 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. max_bp_span_pair.1 + 1 {
      let min_gap_num_4_il = max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2);
      // if min_gap_num_4_il > max_gap_num_4_il {continue;}
      if min_gap_num_4_il > MAX_GAP_NUM_4_IL {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !sta_fe_params.bpa_score_mat.contains_key(&pos_quadruple) {continue;}
          // println!("i, j, k, l: {:?}.", &pos_quadruple);
          let max_sta_free_energies = get_max_sta_free_energies(seq_pair, seq_len_pair, sta_fe_params, max_gap_num, max_gap_num_4_il, &pos_quadruple, &max_sta_free_energy_mat_sets);
          let mut max = NEG_INFINITY;
          let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
          if max_sta_free_energies.free_energy_on_sa > NEG_INFINITY {
            let hl_free_energy = bpa_score + get_hl_fe(seq_pair.0, &(i, j)) + get_hl_fe(seq_pair.1, &(k, l)) + max_sta_free_energies.free_energy_on_sa;
            if hl_free_energy > max {max = hl_free_energy;}
          }
          if max_sta_free_energies.free_energy_4_internal_2loop > NEG_INFINITY {
            let twoloop_free_energy = bpa_score + max_sta_free_energies.free_energy_4_internal_2loop;
            if twoloop_free_energy > max {max = twoloop_free_energy;}
          }
          if max_sta_free_energies.free_energy_4_internal_multiloop > NEG_INFINITY {
            let bp_closing_loop_pair = (
              (seq_pair.0[i], seq_pair.0[j]),
              (seq_pair.1[k], seq_pair.1[l]),
            );
            let ml_tm_delta_fe_pair = (
              ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[j - 1])))],
              ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[l - 1])))],
            );
            let au_or_gu_end_penalty_delta_fe_pair = (
              if is_au_or_gu(&bp_closing_loop_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              if is_au_or_gu(&bp_closing_loop_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
            );
            let ml_free_energy = bpa_score + 2. * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe_pair.0 + ml_tm_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1 + max_sta_free_energies.free_energy_4_internal_multiloop;
            if ml_free_energy > max {max = ml_free_energy;}
          }
          if max > NEG_INFINITY {
            max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas.insert(pos_quadruple, max);
            let accessible_bp_pair = (
              (seq_pair.0[i], seq_pair.0[j]),
              (seq_pair.1[k], seq_pair.1[l]),
            );
            let stacking_bp_pair = (
              (seq_pair.0[i - 1], seq_pair.0[j + 1]),
              (seq_pair.1[k - 1], seq_pair.1[l + 1]),
            );
            let ml_tm_or_de_delta_fe_pair = (
              if i > 1 && j < seq_len_pair.0 - 2 {
                ML_TM_DELTA_FES[&(accessible_bp_pair.0, stacking_bp_pair.0)]
              } else if i > 1 {
                FIVE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).0)]
              } else if j < seq_len_pair.0 - 2 {
                THREE_PRIME_DE_DELTA_FES[&(accessible_bp_pair.0, (stacking_bp_pair.0).1)]
              } else {
                0.
              },
              if k > 1 && l < seq_len_pair.1 - 2 {
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
              if is_au_or_gu(&accessible_bp_pair.0) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
              if is_au_or_gu(&accessible_bp_pair.1) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.},
            );
            let el_free_energy = ml_tm_or_de_delta_fe_pair.0 + ml_tm_or_de_delta_fe_pair.1 + au_or_gu_end_penalty_delta_fe_pair.0 + au_or_gu_end_penalty_delta_fe_pair.1;
            max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_els.insert(pos_quadruple, el_free_energy);
            max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, el_free_energy + 2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
          }
        }
      }
    }
  }
  let leftmost_pos_pair = (0, 0);
  max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat.insert(leftmost_pos_pair, 0.);
  max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas.insert(leftmost_pos_pair, 0.);
  for i in 0 .. seq_len_pair.0 - 1 {
    for j in 0 .. seq_len_pair.1 - 1 {
      let pos_pair = (i, j);
      if pos_pair == (0, 0) {continue;}
      /* let min_gap_num = get_seq_len_diff(&(0, i, 0, j)) + get_seq_len_diff(&(i, seq_len_pair.0, j, seq_len_pair.1));
      if min_gap_num > max_gap_num {continue;} */
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
      let mut max = NEG_INFINITY;
      for k in 1 .. i {
        for l in 1 .. j {
          let pos_pair_2 = (k - 1, l - 1);
          let pos_quadruple = (k, i, l, j);
          /* let min_gap_num = get_seq_len_diff(&(0, k, 0, l)) + get_seq_len_diff(&pos_quadruple) + get_seq_len_diff(&(i, seq_len_pair.0, j, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !(max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat.contains_key(&pos_pair_2) && max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_els.contains_key(&pos_quadruple)) {
            continue;
          }
          let free_energy_4_bpa = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat[&pos_pair_2] + max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_els[&pos_quadruple];
          if free_energy_4_bpa > max {max = free_energy_4_bpa;}
        }
      }
      if i > 0 && j > 0 && max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat.contains_key(&(i - 1, j - 1)) {
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        let free_energy_4_ba = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat[&(i - 1, j - 1)] + ba_score;
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if max > NEG_INFINITY {
        max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas.insert(pos_pair, max);
      }
      /* let min_gap_num_2 = min_gap_num /* + 2 */;
      if min_gap_num_2 <= max_gap_num { */
        max = NEG_INFINITY;
        if i > 0 {
          let pos_pair_2 = (i - 1, j);
          if max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
            let free_energy_4_opening_gap = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
          }
          if max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_1.contains_key(&pos_pair_2) {
            let free_energy_4_extending_gap = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
          }
        }
        if max > NEG_INFINITY {
          max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_1.insert(pos_pair, max);
        }
        max = NEG_INFINITY;
        if j > 0 {
          let pos_pair_2 = (i, j - 1);
          if max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
            let free_energy_4_opening_gap = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
            if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
          }
          if max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_2.contains_key(&pos_pair_2) {
            let free_energy_4_extending_gap = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
            if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
          }
        }
        if max > NEG_INFINITY {
          max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_2.insert(pos_pair, max);
        }
      // }
      max = NEG_INFINITY;
      if max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas.contains_key(&pos_pair) {
        let free_energy_4_ba = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_bas[&pos_pair];
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_1.contains_key(&pos_pair) {
        let free_energy_4_gap = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_1[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_2.contains_key(&pos_pair) {
        let free_energy_4_gap = max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat_4_gaps_2[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if max > NEG_INFINITY {
        max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat.insert(pos_pair, max);
      }
    }
  }
  max_sta_free_energy_mat_sets.free_energy_mats_4_external_loop.free_energy_mat[&(seq_len_pair.0 - 2, seq_len_pair.1 - 2)]
}

#[inline]
pub fn get_max_sta_free_energies(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_gap_num: usize, max_gap_num_4_il: usize, pos_quadruple: &PosQuadruple, max_sta_free_energy_mat_sets: &MaxStaFreeEnergyMatSets) -> MaxStaFreeEnergies {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut tmp_max_sta_free_energy_mat_sets = TmpMaxStaFreeEnergyMatSets::new();
  let &(i, j, k, l) = pos_quadruple;
  for u in i .. j {
    for v in k .. l {
      let pos_pair = (u, v);
      /* let min_gap_num_4_il = get_seq_len_diff(&(i, u, k, v)) + get_seq_len_diff(&(u, j, v, l));
      if min_gap_num_4_il > max_gap_num_4_il {continue;}
      let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
      if min_gap_num > max_gap_num {continue;} */
      // if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
      // if !(is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num_4_il)) {continue;}
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, MAX_GAP_NUM_4_IL) {continue;}
      // println!("u, v: {:?}.", &pos_pair);
      if u == i && v == k {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas.insert(pos_pair, 0.);
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat.insert(pos_pair, 0.);
        continue;
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
      let pos_pair_2 = (u - 1, v - 1);
      let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat.contains_key(&pos_pair_2) {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas.insert(pos_pair,
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat[&pos_pair_2] + ba_score
        );
      }
      // println!("OK.");
      let mut max = NEG_INFINITY;
      /* let min_gap_num_4_il_2 = min_gap_num_4_il /* + 2 */;
      let min_gap_num_2 = min_gap_num /* + 2 */;
      if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num { */
        let pos_pair_2 = (u - 1, v);
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_1.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_1.insert(pos_pair, max);
        }
        let pos_pair_2 = (u, v - 1);
        max = NEG_INFINITY;
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_2.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_2.insert(pos_pair, max);
        }
      // }
      max = NEG_INFINITY;
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas.contains_key(&pos_pair) {
        let free_energy_4_ba = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_bas[&pos_pair];
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_1.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_1[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_2.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat_4_gaps_2[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if max > NEG_INFINITY {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat.insert(pos_pair, max);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop"
      max = NEG_INFINITY;
      for m in i + 1 .. u {
        if m - i - 1 + j - u - 1 > MAX_2_LOOP_LEN {continue;}
        for n in k + 1 .. v {
          if n - k - 1 + l - v - 1 > MAX_2_LOOP_LEN {continue;}
          let pos_pair_2 = (m - 1, n - 1);
          let pos_quadruple_2 = (m, u, n, v);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, m, k, n)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(u, j, v, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
          let twoloop_fe_pair = (
            get_2_loop_fe(seq_pair.0, &(i, j), &(m , u)),
            get_2_loop_fe(seq_pair.1, &(k, l), &(n , v)),
          );
          if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat.contains_key(&pos_pair_2) {
            let free_energy_4_bpa = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat[&pos_pair_2] + twoloop_fe_pair.0 + twoloop_fe_pair.1 + max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas[&pos_quadruple_2];
            if free_energy_4_bpa > max {max = free_energy_4_bpa;}
          }
        }
      }
      let pos_pair_2 = (u - 1, v - 1);
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat.contains_key(&pos_pair_2) {
        let free_energy_4_ba = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat[&pos_pair_2] + ba_score;
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if max > NEG_INFINITY {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_bas.insert(pos_pair, max);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u - 1, v);
        max = NEG_INFINITY;
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_1.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_1.insert(pos_pair, max);
        }
        let pos_pair_2 = (u, v - 1);
        max = NEG_INFINITY;
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_2.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_2.insert(pos_pair, max);
        }
      // }
      max = NEG_INFINITY;
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_bas.contains_key(&pos_pair) {
        let free_energy_4_ba = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_bas[&pos_pair];
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_1.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_1[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_2.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat_4_gaps_2[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if max > NEG_INFINITY {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat.insert(pos_pair, max);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
      max = NEG_INFINITY;
      for m in i + 1 .. u {
        for n in k + 1 .. v {
          let pos_pair_2 = (m - 1, n - 1);
          let pos_quadruple_2 = (m, u, n, v);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, m, k, n)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(u, j, v, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2) {
            continue;
          }
          if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat.contains_key(&pos_pair_2) {
            let free_energy_4_bpa = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat[&pos_pair_2] + max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2];
            if free_energy_4_bpa > max {max = free_energy_4_bpa;}
          }
          if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat.contains_key(&pos_pair_2) {
            let free_energy_4_bpa = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat[&pos_pair_2] + max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2];
            if free_energy_4_bpa > max {max = free_energy_4_bpa;}
          }
        }
      }
      let pos_pair_2 = (u - 1, v - 1);
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat.contains_key(&pos_pair_2) {
        let free_energy_4_ba = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat[&pos_pair_2] + ba_score;
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if max > NEG_INFINITY {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_bas.insert(pos_pair, max);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u - 1, v);
        max = NEG_INFINITY;
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_1.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_1.insert(pos_pair, max);
        }
        let pos_pair_2 = (u, v - 1);
        max = NEG_INFINITY;
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_2.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_2.insert(pos_pair, max);
        }
      // }
      max = NEG_INFINITY;
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_bas.contains_key(&pos_pair) {
        let free_energy_4_ba = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_bas[&pos_pair];
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_1.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_1[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_2.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat_4_gaps_2[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if max > NEG_INFINITY {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat.insert(pos_pair, max);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
      max = NEG_INFINITY;
      for m in i + 1 .. u {
        for n in k + 1 .. v {
          let pos_pair_2 = (m - 1, n - 1);
          let pos_quadruple_2 = (m, u, n, v);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, m, k, n)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(u, j, v, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !(tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat.contains_key(&pos_pair_2) && max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2)) {
            continue;
          }
          let free_energy_4_bpa = tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat[&pos_pair_2] + max_sta_free_energy_mat_sets.free_energy_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2];
          if free_energy_4_bpa > max {max = free_energy_4_bpa;}
        }
      }
      let pos_pair_2 = (u - 1, v - 1);
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat.contains_key(&pos_pair_2) {
        let free_energy_4_ba = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat[&pos_pair_2] + ba_score;
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if max > NEG_INFINITY {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_bas.insert(pos_pair, max);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u - 1, v);
        max = NEG_INFINITY;
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_1.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_1[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_1.insert(pos_pair, max);
        }
        let pos_pair_2 = (u, v - 1);
        max = NEG_INFINITY;
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_bas.contains_key(&pos_pair_2) {
          let free_energy_4_opening_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_bas[&pos_pair_2] + sta_fe_params.opening_gap_penalty;
          if free_energy_4_opening_gap > max {max = free_energy_4_opening_gap;}
        }
        if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_2.contains_key(&pos_pair_2) {
          let free_energy_4_extending_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_2[&pos_pair_2] + sta_fe_params.extending_gap_penalty;
          if free_energy_4_extending_gap > max {max = free_energy_4_extending_gap;}
        }
        if max > NEG_INFINITY {
          tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_2.insert(pos_pair, max);
        }
      // }
      max = NEG_INFINITY;
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_bas.contains_key(&pos_pair) {
        let free_energy_4_ba = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_bas[&pos_pair];
        if free_energy_4_ba > max {max = free_energy_4_ba;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_1.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_1[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_2.contains_key(&pos_pair) {
        let free_energy_4_gap = tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat_4_gaps_2[&pos_pair];
        if free_energy_4_gap > max {max = free_energy_4_gap;}
      }
      if max > NEG_INFINITY {
        tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_first_bpas_on_mls.free_energy_mat.insert(pos_pair, max);
      }
    }
  }
  let rightmost_pos_pair = (j - 1 , l - 1);
  MaxStaFreeEnergies {
    free_energy_on_sa: if tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat.contains_key(&rightmost_pos_pair) {tmp_max_sta_free_energy_mat_sets.free_energy_mats_on_sa.free_energy_mat[&rightmost_pos_pair]} else {NEG_INFINITY},
    free_energy_4_internal_2loop: if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat.contains_key(&rightmost_pos_pair) {tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_2loop.free_energy_mat[&rightmost_pos_pair]} else {NEG_INFINITY},
    free_energy_4_internal_multiloop: if tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat.contains_key(&rightmost_pos_pair) {tmp_max_sta_free_energy_mat_sets.free_energy_mats_4_internal_multiloop.free_energy_mat[&rightmost_pos_pair]} else {NEG_INFINITY},
  }
}

#[inline]
pub fn get_forward_tmp_sta_inside_part_func_mat_sets(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_gap_num: usize, max_gap_num_4_il: usize, pos_quadruple: &PosQuadruple, invert_exp_max_free_energy: FreeEnergy, sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets) -> TmpStaInsidePartFuncMatSets {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut forward_tmp_sta_inside_part_func_mat_sets = TmpStaInsidePartFuncMatSets::new();
  let &(i, j, k, l) = pos_quadruple;
  for u in i .. j {
    for v in k .. l {
      let pos_pair = (u, v);
      /* let min_gap_num_4_il = get_seq_len_diff(&(i, u, k, v)) + get_seq_len_diff(&(u, j, v, l));
      if min_gap_num_4_il > max_gap_num_4_il {continue;} */
      // if min_gap_num_4_il > max_gap_num {continue;}
      /* let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
      if min_gap_num > max_gap_num {continue;} */
      // if !(is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num_4_il)) {continue;}
      // if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
      // if !(is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num_4_il)) {continue;}
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
      if u == i && v == k {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair, invert_exp_max_free_energy);
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, invert_exp_max_free_energy);
        continue;
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
      let pos_pair_2 = (u - 1, v - 1);
      let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_pair_2) {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair,
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_pair_2] * exp_ba_score
        );
      }
      let mut sum = 0.;
      /* let min_gap_num_4_il_2 = min_gap_num_4_il /* + 2 */;
      let min_gap_num_2 = min_gap_num /* + 2 */;
      if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num { */
        let pos_pair_2 = (u - 1, v);
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = 0.;
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, sum);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop"
      sum = 0.;
      for m in i + 1 .. u {
        if m - i - 1 + j - u - 1 > MAX_2_LOOP_LEN {continue;}
        for n in k + 1 .. v {
          if n - k - 1 + l - v - 1 > MAX_2_LOOP_LEN {continue;}
          let pos_pair_2 = (m - 1, n - 1);
          let pos_quadruple_2 = (m, u, n, v);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, m, k, n)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(u, j, v, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;} */
          /* let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
          let exp_2loop_fe_pair = (
            get_exp_2_loop_fe(seq_pair.0, &(i, j), &(m , u)),
            get_exp_2_loop_fe(seq_pair.1, &(k, l), &(n , v)),
          );
          if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_pair_2) {
            sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_pair_2] * exp_2loop_fe_pair.0 * exp_2loop_fe_pair.1 * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas[&pos_quadruple_2] / invert_exp_max_free_energy;
          }
        }
      }
      let pos_pair_2 = (u - 1, v - 1);
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat.contains_key(&pos_pair_2) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat[&pos_pair_2] * exp_ba_score;
      }
      if sum > 0. {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_bas.insert(pos_pair, sum);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u - 1, v);
        sum = 0.;
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = 0.;
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_bas[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_1[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_2loop.part_func_mat.insert(pos_pair, sum);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
      sum = 0.;
      for m in i + 1 .. u {
        for n in k + 1 .. v {
          let pos_pair_2 = (m - 1, n - 1);
          let pos_quadruple_2 = (m, u, n, v);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, m, k, n)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(u, j, v, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2) {
            continue;
          }
          if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&pos_pair_2) {
            sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] / invert_exp_max_free_energy;
          }
          if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&pos_pair_2) {
            sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] / invert_exp_max_free_energy;
          }
        }
      }
      let pos_pair_2 = (u - 1, v - 1);
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&pos_pair_2) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&pos_pair_2] * exp_ba_score;
      }
      if sum > 0. {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.insert(pos_pair, sum);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u - 1, v);
        sum = 0.;
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = 0.;
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.insert(pos_pair, sum);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
      sum = 0.;
      for m in i + 1 .. u {
        for n in k + 1 .. v {
          let pos_pair_2 = (m - 1, n - 1);
          let pos_quadruple_2 = (m, u, n, v);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, m, k, n)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(u, j, v, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !(forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_pair_2) && sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2)) {
            continue;
          }
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] / invert_exp_max_free_energy;
        }
      }
      let pos_pair_2 = (u - 1, v - 1);
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&pos_pair_2) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&pos_pair_2] * exp_ba_score;
      }
      if sum > 0. {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.insert(pos_pair, sum);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u - 1, v);
        sum = 0.;
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v - 1);
        sum = 0.;
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1[&pos_pair];
      }
      if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.insert(pos_pair, sum);
      }
    }
  }
  forward_tmp_sta_inside_part_func_mat_sets
}

#[inline]
pub fn get_backward_tmp_sta_inside_part_func_mat_sets(seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_gap_num: usize, max_gap_num_4_il: usize, pos_quadruple: &PosQuadruple, invert_exp_max_free_energy: FreeEnergy, sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets) -> TmpStaInsidePartFuncMatSets {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut backward_tmp_sta_inside_part_func_mat_sets = TmpStaInsidePartFuncMatSets::new();
  let &(i, j, k, l) = pos_quadruple;
  for u in (i + 1 .. j + 1).rev() {
    for v in (k + 1 .. l + 1).rev() {
      let pos_pair = (u, v);
      /* let min_gap_num_4_il = get_seq_len_diff(&(i, u, k, v)) + get_seq_len_diff(&(u, j, v, l));
      if min_gap_num_4_il > max_gap_num_4_il {continue;}
      let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
      if min_gap_num > max_gap_num {continue;} */
      // if !(is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num_4_il)) {continue;}
      // if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
      // if !(is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) && is_min_gap_ok_1(&pos_pair, &pos_quadruple, max_gap_num_4_il)) {continue;}
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
      if u == j && v == l {
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair, invert_exp_max_free_energy);
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, invert_exp_max_free_energy);
        continue;
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
      let pos_pair_2 = (u + 1, v + 1);
      let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_pair_2) {
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.insert(pos_pair,
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_pair_2] * exp_ba_score
        );
      }
      let mut sum = 0.;
      /* let min_gap_num_4_il_2 = min_gap_num_4_il /* + 2 */;
      let min_gap_num_2 = min_gap_num /* + 2 */;
      if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num { */
        let pos_pair_2 = (u + 1, v);
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = 0.;
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_bas[&pos_pair];
      }
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_1[&pos_pair];
      }
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_pair, sum);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
      sum = 0.;
      for m in u + 1 .. j {
        for n in v + 1 .. l {
          let pos_pair_2 = (m + 1, n + 1);
          let pos_quadruple_2 = (u, m, v, n);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, u, k, v)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(m, j, n, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + get_seq_len_diff(&(i, u, k, v)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(m, j, n, l)) + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2) {
            continue;
          }
          if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&pos_pair_2) {
            sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] / invert_exp_max_free_energy;
          }
          if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&pos_pair_2) {
            sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] / invert_exp_max_free_energy;
          }
        }
      }
      let pos_pair_2 = (u + 1, v + 1);
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&pos_pair_2) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&pos_pair_2] * exp_ba_score;
      }
      if sum > 0. {
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.insert(pos_pair, sum);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u + 1, v);
        sum = 0.;
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = 0.;
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_bas[&pos_pair];
      }
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_1[&pos_pair];
      }
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.insert(pos_pair, sum);
      }
      // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
      sum = 0.;
      for m in u + 1 .. j {
        for n in v + 1 .. l {
          let pos_pair_2 = (m + 1, n + 1);
          let pos_quadruple_2 = (u, m, v, n);
          /* let min_gap_num_4_il = get_seq_len_diff(&(i, u, k, v)) + get_seq_len_diff(&pos_quadruple_2) + get_seq_len_diff(&(m, j, n, l));
          if min_gap_num_4_il > max_gap_num_4_il {continue;}
          let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
          if min_gap_num > max_gap_num {continue;} */
          if !(backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_pair_2) && sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2)) {
            continue;
          }
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] / invert_exp_max_free_energy;
        }
      }
      let pos_pair_2 = (u + 1, v + 1);
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&pos_pair_2) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&pos_pair_2] * exp_ba_score;
      }
      if sum > 0. {
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.insert(pos_pair, sum);
      }
      // if min_gap_num_4_il_2 <= max_gap_num_4_il && min_gap_num_2 <= max_gap_num {
        let pos_pair_2 = (u + 1, v);
        sum = 0.;
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.insert(pos_pair, sum);
        }
        let pos_pair_2 = (u, v + 1);
        sum = 0.;
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
        if sum > 0. {
          backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.insert(pos_pair, sum);
        }
      // }
      sum = 0.;
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_bas[&pos_pair];
      }
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_1[&pos_pair];
      }
      if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.insert(pos_pair, sum);
      }
    }
  }
  backward_tmp_sta_inside_part_func_mat_sets
}

#[inline]
pub fn get_sta_outside_part_func_4d_mat_4_bpas(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_gap_num: usize, max_gap_num_4_il: usize, sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets, farthest_pos_pairs_with_pos_pairs: &PosPairsWithPosPairs, bools_with_pair_aligned_right_pos_pairs: &BoolsWithPosPairs) -> PartFunc4dMat {
  let mut sta_outside_part_func_4d_mat_4_bpas = PartFunc4dMat::default();
  let invert_exp_max_free_energy = sta_fe_params.invert_exp_max_free_energy;
  let mut sta_outside_part_func_4d_mat_4_right = sta_outside_part_func_4d_mat_4_bpas.clone();
  let mut sta_outside_part_func_4d_mat_4_right_2 = sta_outside_part_func_4d_mat_4_bpas.clone();
  for substr_len_1 in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. max_bp_span_pair.0 + 1).rev() {
    for substr_len_2 in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL .. max_bp_span_pair.1 + 1).rev() {
      let min_gap_num_4_il = max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2);
      // if min_gap_num_4_il > max_gap_num_4_il {continue;}
      // if min_gap_num_4_il > max_gap_num_4_il {continue;}
      if min_gap_num_4_il > max_gap_num_4_il {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !bools_with_pair_aligned_right_pos_pairs.contains_key(&(j, l)) {continue;}
          let mut sum = 0.;
          let mut sum_2 = sum;
          for m in j + 1 .. seq_len_pair.0 - 1 {
            for n in l + 1 .. seq_len_pair.1 - 1 {
              let pos_quadruple_2 = (i, m, k, n);
              if !sta_outside_part_func_4d_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
              // let min_gap_num = get_seq_len_diff(&(0, m, 0, o)) + get_seq_len_diff(&(m, i, o, k)) + get_seq_len_diff(&pos_quadruple) + get_seq_len_diff(&(j, n, l, p)) + get_seq_len_diff(&(n, seq_len_pair.0, p, seq_len_pair.1));
              /* let min_gap_num_4_il = min_gap_num_4_il + get_seq_len_diff(&(j, m, l, n));
              if min_gap_num_4_il > max_gap_num_4_il {continue;}
              let min_gap_num = get_seq_len_diff(&(0, i, 0, k)) + min_gap_num_4_il + get_seq_len_diff(&(m, seq_len_pair.0, n, seq_len_pair.1));
              if min_gap_num > max_gap_num {continue;} */
              // let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
              let ref backward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
              let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple_2];
              let bp_closing_loop_pair = (
                (seq_pair.0[i], seq_pair.0[m]),
                (seq_pair.1[k], seq_pair.1[n]),
              );
              let exp_ml_tm_delta_fe_pair = (
                EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[i + 1], seq_pair.0[m - 1])))],
                EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[k + 1], seq_pair.1[n - 1])))],
              );
              let exp_au_or_gu_end_penalty_delta_fe_pair = (
                if is_au_or_gu(&bp_closing_loop_pair.0) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
                if is_au_or_gu(&bp_closing_loop_pair.1) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
              );
              // let coefficient = exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe_pair.0 * exp_ml_tm_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1 * sta_outside_part_func_4d_mat_4_bpas[&pos_quadruple_2] / invert_exp_max_free_energy;
              let coefficient = exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe_pair.0 * exp_ml_tm_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1 * sta_outside_part_func_4d_mat_4_bpas[&pos_quadruple_2] / invert_exp_max_free_energy;
              if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(j + 1, l + 1)) {
                sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(j + 1, l + 1)];
              }
              if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(j + 1, l + 1)) {
                sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(j + 1, l + 1)];
              }
              if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(j + 1, l + 1)) {
                sum_2 += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(j + 1, l + 1)];
              }
            }
          }
          if sum > 0. {
            sta_outside_part_func_4d_mat_4_right.insert(pos_quadruple, sum);
          }
          if sum_2 > 0. {
            sta_outside_part_func_4d_mat_4_right_2.insert(pos_quadruple, sum_2);
          }
          if !sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.contains_key(&pos_quadruple) {continue;}
          let mut sum = 0.;
          if sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat.contains_key(&(i - 1, k - 1)) && sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat.contains_key(&(j + 1, l + 1)) {
            sum += sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat[&(i - 1, k - 1)] * sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_els[&pos_quadruple] / sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas[&pos_quadruple] * sta_inside_part_func_mat_sets.backward_part_func_mats_4_external_loop.part_func_mat[&(j + 1, l + 1)] / invert_exp_max_free_energy;
          }
          for m in 1 .. i {
            for n in j + 1 .. seq_len_pair.0 - 1 {
              if n - j - 1 + i - m - 1 > MAX_2_LOOP_LEN {continue;}
              for o in 1 .. k {
                for p in l + 1 .. seq_len_pair.1 - 1 {
                  if p - l - 1 + k - o - 1 > MAX_2_LOOP_LEN {continue;}
                  let pos_quadruple_2 = (m, n, o, p);
                  if !sta_outside_part_func_4d_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
                  /* let min_gap_num_4_il = get_seq_len_diff(&(m, i, o, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, n, l, p));
                  if min_gap_num_4_il > max_gap_num_4_il {continue;}
                  let min_gap_num = get_seq_len_diff(&(0, m, 0, o)) + min_gap_num_4_il + get_seq_len_diff(&(n, seq_len_pair.0, p, seq_len_pair.1));
                  if min_gap_num > max_gap_num {continue;} */
                  let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                  let ref backward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                  if !(forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(i - 1, k - 1)) && backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(j + 1, l + 1))) {continue;}
                  let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple_2];
                  let exp_2loop_fe_pair = (
                    get_exp_2_loop_fe(seq_pair.0, &(m, n), &(i , j)),
                    get_exp_2_loop_fe(seq_pair.1, &(o, p), &(k , l)),
                  );
                  sum += forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(i - 1, k - 1)] / invert_exp_max_free_energy * exp_bpa_score * exp_2loop_fe_pair.0 * exp_2loop_fe_pair.1 * sta_outside_part_func_4d_mat_4_bpas[&pos_quadruple_2] * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(j + 1, l + 1)] / invert_exp_max_free_energy;
                }
              }
            }
          }
          let part_func_ratio = sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple] / sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas[&pos_quadruple];
          /* for m in 1 .. i {
            for n in j + 1 .. seq_len_pair.0 - 1 {
              for o in 1 .. k {
                for p in l + 1 .. seq_len_pair.1 - 1 {
                  let pos_quadruple_2 = (m, n, o, p);
                  /* let min_gap_num_4_il = get_seq_len_diff(&(m, i, o, k)) + min_gap_num_4_il + get_seq_len_diff(&(j, n, l, p));
                  if min_gap_num_4_il > max_gap_num_4_il {continue;}
                  let min_gap_num = get_seq_len_diff(&(0, m, 0, n)) + min_gap_num_4_il + get_seq_len_diff(&(j, seq_len_pair.0, l, seq_len_pair.1));
                  if min_gap_num > max_gap_num {continue;} */
                  if !sta_outside_part_func_4d_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
                  let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple_2];
                  let bp_closing_loop_pair = (
                    (seq_pair.0[m], seq_pair.0[n]),
                    (seq_pair.1[o], seq_pair.1[p]),
                  );
                  let exp_ml_tm_delta_fe_pair = (
                    EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.0), invert_bp(&(seq_pair.0[m + 1], seq_pair.0[n - 1])))],
                    EXP_ML_TM_DELTA_FES[&(invert_bp(&bp_closing_loop_pair.1), invert_bp(&(seq_pair.1[o + 1], seq_pair.1[p - 1])))],
                  );
                  let exp_au_or_gu_end_penalty_delta_fe_pair = (
                    if is_au_or_gu(&bp_closing_loop_pair.0) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
                    if is_au_or_gu(&bp_closing_loop_pair.1) {*EXP_HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {1.},
                  );
                  let coefficient = part_func_ratio * exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe_pair.0 * exp_ml_tm_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1 * sta_outside_part_func_4d_mat_4_bpas[&pos_quadruple_2] / invert_exp_max_free_energy;
                  /* if !farthest_pos_pairs_with_pos_pairs.contains_key(&(m, n)) {continue;}
                  let (o, p) = farthest_pos_pairs_with_pos_pairs[&(m, n)];
                  let pos_quadruple_2 = (m, o, n, p); */
                  let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                  let ref backward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.backward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
                  if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(i - 1, k - 1)) {
                    let coefficient = coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(i - 1, k - 1)] / invert_exp_max_free_energy;
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(j + 1, l + 1)];
                    }
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(j + 1, l + 1)];
                    }
                  }
                  if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(i - 1, k - 1)) {
                    let coefficient = coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(i - 1, k - 1)] / invert_exp_max_free_energy;
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(j + 1, l + 1)];
                    }
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(j + 1, l + 1)];
                    }
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(j + 1, l + 1)];
                    }
                  }
                  if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(i - 1, k - 1)) {
                    let coefficient = coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(i - 1, k - 1)] / invert_exp_max_free_energy;
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(j + 1, l + 1)];
                    }
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(j + 1, l + 1)];
                    }
                    if backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(j + 1, l + 1)) {
                      sum += coefficient * backward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(j + 1, l + 1)];
                    }
                  }
                }
              }
            }
          } */
          for m in 1 .. i {
            for n in 1 .. k {
              if !farthest_pos_pairs_with_pos_pairs.contains_key(&(m, n)) {continue;}
              let (o, p) = farthest_pos_pairs_with_pos_pairs[&(m, n)];
              let pos_quadruple_2 = (m, o, n, p);
              let ref forward_tmp_sta_inside_part_func_mat_sets = sta_inside_part_func_mat_sets.forward_tmp_sta_inside_part_func_mat_sets_with_pos_quadruples[&pos_quadruple_2];
              let pos_quadruple_3 = (m, j, n, l);
              if sta_outside_part_func_4d_mat_4_right.contains_key(&pos_quadruple_3) {
                let sta_outside_part_func_4_right = sta_outside_part_func_4d_mat_4_right[&pos_quadruple_3];
                let coefficient = sta_outside_part_func_4_right * part_func_ratio / invert_exp_max_free_energy;
                if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(i - 1, k - 1)) {
                  sum += coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(i - 1, k - 1)];
                }
                if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(i - 1, k - 1)) {
                  sum += coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(i - 1, k - 1)];
                }
                if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(i - 1, k - 1)) {
                  sum += coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(i - 1, k - 1)];
                }
              }
              if sta_outside_part_func_4d_mat_4_right_2.contains_key(&pos_quadruple_3) {
                let sta_outside_part_func_4_right_2 = sta_outside_part_func_4d_mat_4_right_2[&pos_quadruple_3];
                let coefficient = sta_outside_part_func_4_right_2 * part_func_ratio / invert_exp_max_free_energy;
                if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(i - 1, k - 1)) {
                  sum += coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(i - 1, k - 1)];
                }
                if forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(i - 1, k - 1)) {
                  sum += coefficient * forward_tmp_sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(i - 1, k - 1)];
                }
              }
            }
          }
          if sum > 0. {
            sta_outside_part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
          }
        }
      }
    }
  }
  sta_outside_part_func_4d_mat_4_bpas
}

#[inline]
fn get_bpap_mat(seq_len_pair: &(usize, usize), sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets, sta_outside_part_func_4d_mat_4_bpas: &PartFunc4dMat, sta_fe_params: &StaFeParams) -> Prob4dMat {
  let mut bpap_mat = Prob4dMat::default();
  let invert_exp_max_free_energy = sta_fe_params.invert_exp_max_free_energy;
  let part_func = sta_inside_part_func_mat_sets.forward_part_func_mats_4_external_loop.part_func_mat[&(seq_len_pair.0 - 2, seq_len_pair.1 - 2)];
  // println!("part_func: {}", part_func);
  for pos_quadruple in sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas.keys() {
    if !sta_outside_part_func_4d_mat_4_bpas.contains_key(pos_quadruple) {continue;}
    let bpap = sta_inside_part_func_mat_sets.part_func_4d_mat_4_bpas[pos_quadruple] * sta_outside_part_func_4d_mat_4_bpas[pos_quadruple] / (part_func * invert_exp_max_free_energy);
    debug_assert!(0. <= bpap && bpap <= 1.);
    bpap_mat.insert(*pos_quadruple, bpap);
  }
  bpap_mat
}
 
#[inline]
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
      *new_bpp_mat.get_mut(&pos_pair).expect("Failed to get an element from a hash map.") += weighted_bpap;
      new_upp_mat[pos_pair.0] -= weighted_bpap;
      new_upp_mat[pos_pair.1] -= weighted_bpap;
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
pub fn remove_small_bpps_from_bpp_mat(sparse_bpp_mat: &SparseProbMat, min_bpp: Prob) -> SparseProbMat {
  sparse_bpp_mat.iter().filter(|(_, &bpp)| {bpp >= min_bpp}).map(|(&(i, j), &bpp)| {((i + 1, j + 1), bpp)}).collect()
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
pub fn get_seq_len_diff(pos_quadruple: &PosQuadruple) -> usize {
  let seq_len_pair = (pos_quadruple.1 - pos_quadruple.0 + 1, pos_quadruple.3 - pos_quadruple.2 + 1);
  max(seq_len_pair.0, seq_len_pair.1) - min(seq_len_pair.0, seq_len_pair.1)
}
#[inline]
fn is_min_gap_ok_1(pos_pair: &PosPair, pos_quadruple: &PosQuadruple, max_gap_num: usize) -> bool {
  let min_gap_num_1 = get_seq_len_diff(&(pos_quadruple.0, pos_pair.0, pos_quadruple.2, pos_pair.1));
  let min_gap_num_2 = get_seq_len_diff(&(pos_pair.0, pos_quadruple.1, pos_pair.1, pos_quadruple.3));
  if min_gap_num_1 + min_gap_num_2 <= max_gap_num {
    true
  } else {
    false
  }
}

#[inline]
fn is_min_gap_ok_2(pos_quadruple: &PosQuadruple, max_gap_num: usize) -> bool {
  let min_gap_num = get_seq_len_diff(&pos_quadruple);
  if min_gap_num <= max_gap_num {
    true
  } else {
    false
  }
}
