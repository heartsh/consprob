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
pub struct StaInsidePartFunc4dMats {
  pub part_func_mat: PartFunc4dMat,
  pub forward_part_func_mat_4_bas: PartFunc4dMat,
  pub forward_part_func_mat_4_gaps_1: PartFunc4dMat,
  pub forward_part_func_mat_4_gaps_2: PartFunc4dMat,
  pub backward_part_func_mat_4_bas: PartFunc4dMat,
  pub backward_part_func_mat_4_gaps_1: PartFunc4dMat,
  pub backward_part_func_mat_4_gaps_2: PartFunc4dMat,
}
#[derive(Clone)]
pub struct StaInsidePartFuncMats {
  pub forward_part_func_mat: SparsePartFuncMat,
  pub forward_part_func_mat_4_bas: SparsePartFuncMat,
  pub forward_part_func_mat_4_gaps_1: SparsePartFuncMat,
  pub forward_part_func_mat_4_gaps_2: SparsePartFuncMat,
  pub backward_part_func_mat: SparsePartFuncMat,
  pub backward_part_func_mat_4_bas: SparsePartFuncMat,
  pub backward_part_func_mat_4_gaps_1: SparsePartFuncMat,
  pub backward_part_func_mat_4_gaps_2: SparsePartFuncMat,
}
#[derive(Clone)]
pub struct StaInsidePartFuncMatSets {
  pub part_func_mat_4_bpas: PartFunc4dMat,
  pub part_func_mat_4_bpas_accessible_on_els: PartFunc4dMat,
  pub part_func_mat_4_bpas_accessible_on_mls: PartFunc4dMat,
  pub part_func_mats_on_sa: StaInsidePartFunc4dMats,
  pub part_func_mats_4_internal_multiloop: StaInsidePartFunc4dMats,
  pub part_func_mats_4_first_bpas_on_mls: StaInsidePartFunc4dMats,
  pub part_func_mats_4_external_loop: StaInsidePartFuncMats,
}
#[derive(Clone)]
pub struct StaOutsidePartFunc4dMats {
  pub part_func_mat_4_bpas: PartFunc4dMat,
  pub part_func_mat_4_bpas_on_el: PartFunc4dMat,
  pub part_func_mat_4_bpas_on_internal_2loops: PartFunc4dMat,
  pub part_func_mat_4_bpas_on_internal_mls: PartFunc4dMat,
  pub part_func_mat_4_right: PartFunc4dMat,
  pub part_func_mat_4_right_2: PartFunc4dMat,
}
pub struct StaFeParams {
  pub exp_ba_score_mat: SparseFreeEnergyMat,
  pub exp_bpa_score_mat: FreeEnergy4dMat,
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

impl StaFeParams {
  pub fn origin() -> StaFeParams {
    StaFeParams {
      exp_ba_score_mat: SparseFreeEnergyMat::default(),
      exp_bpa_score_mat: FreeEnergy4dMat::default(),
      exp_opening_gap_penalty: 0.,
      exp_extending_gap_penalty: 0.,
      invert_exp_max_free_energy: 0.,
    }
  }
  pub fn new(rna_id_pair: &RnaIdPair, fasta_records: &FastaRecords, invert_exp_max_free_energies: &FreeEnergies, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_pos_dist_4_el: usize, max_substr_dist: usize, bpp_mats: &ProbMats, exp_opening_gap_penalty: FreeEnergy, exp_extending_gap_penalty: FreeEnergy) -> StaFeParams {
    let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
    let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
    let mut sta_fe_params = StaFeParams::origin();
    sta_fe_params.invert_exp_max_free_energy = invert_exp_max_free_energies[rna_id_pair.0] * invert_exp_max_free_energies[rna_id_pair.1];
    sta_fe_params.exp_opening_gap_penalty = exp_opening_gap_penalty;
    sta_fe_params.exp_extending_gap_penalty = exp_extending_gap_penalty;
    let bpp_mat_pair = (&bpp_mats[rna_id_pair.0], &bpp_mats[rna_id_pair.1]);
    for i in 1 .. seq_len_pair.0 - 1 {
      let base = seq_pair.0[i];
      for j in 1 .. seq_len_pair.1 - 1 {
        let pos_pair = (i, j);
        if !is_pos_dist_ok(&pos_pair, max_pos_dist_4_el) {continue;}
        let base_pair = (base, seq_pair.1[j]);
        sta_fe_params.exp_ba_score_mat.insert(pos_pair, EXP_RIBOSUM_85_60_BA_SCORE_MAT[&base_pair]);
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
            sta_fe_params.exp_bpa_score_mat.insert(pos_quadruple, EXP_RIBOSUM_85_60_BPA_SCORE_MAT[&base_quadruple]);
          }
        }
      }
    }
    sta_fe_params
  }
}

impl StaInsidePartFunc4dMats {
  pub fn new() -> StaInsidePartFunc4dMats {
    let part_func_4d_mat = PartFunc4dMat::default();
    StaInsidePartFunc4dMats {
      part_func_mat: part_func_4d_mat.clone(),
      forward_part_func_mat_4_bas: part_func_4d_mat.clone(),
      forward_part_func_mat_4_gaps_1: part_func_4d_mat.clone(),
      forward_part_func_mat_4_gaps_2: part_func_4d_mat.clone(),
      backward_part_func_mat_4_bas: part_func_4d_mat.clone(),
      backward_part_func_mat_4_gaps_1: part_func_4d_mat.clone(),
      backward_part_func_mat_4_gaps_2: part_func_4d_mat,
    }
  }
}

impl StaInsidePartFuncMats {
  pub fn new() -> StaInsidePartFuncMats {
    let part_func_mat = SparsePartFuncMat::default();
    StaInsidePartFuncMats {
      forward_part_func_mat: part_func_mat.clone(),
      forward_part_func_mat_4_bas: part_func_mat.clone(),
      forward_part_func_mat_4_gaps_1: part_func_mat.clone(),
      forward_part_func_mat_4_gaps_2: part_func_mat.clone(),
      backward_part_func_mat: part_func_mat.clone(),
      backward_part_func_mat_4_bas: part_func_mat.clone(),
      backward_part_func_mat_4_gaps_1: part_func_mat.clone(),
      backward_part_func_mat_4_gaps_2: part_func_mat,
    }
  }
}

impl StaInsidePartFuncMatSets {
  pub fn new() -> StaInsidePartFuncMatSets {
    let part_func_4d_mat = PartFunc4dMat::default();
    let part_func_4d_mats = StaInsidePartFunc4dMats::new();
    StaInsidePartFuncMatSets {
      part_func_mat_4_bpas: part_func_4d_mat.clone(),
      part_func_mat_4_bpas_accessible_on_els: part_func_4d_mat.clone(),
      part_func_mat_4_bpas_accessible_on_mls: part_func_4d_mat.clone(),
      part_func_mats_on_sa: part_func_4d_mats.clone(),
      part_func_mats_4_internal_multiloop: part_func_4d_mats.clone(),
      part_func_mats_4_first_bpas_on_mls: part_func_4d_mats,
      part_func_mats_4_external_loop: StaInsidePartFuncMats::new(),
    }
  }
}

impl StaOutsidePartFunc4dMats {
  pub fn new() -> StaOutsidePartFunc4dMats {
    let part_func_4d_mat = PartFunc4dMat::default();
    StaOutsidePartFunc4dMats {
      part_func_mat_4_bpas: part_func_4d_mat.clone(),
      part_func_mat_4_bpas_on_el: part_func_4d_mat.clone(),
      part_func_mat_4_bpas_on_internal_2loops: part_func_4d_mat.clone(),
      part_func_mat_4_bpas_on_internal_mls: part_func_4d_mat.clone(),
      part_func_mat_4_right: part_func_4d_mat.clone(),
      part_func_mat_4_right_2: part_func_4d_mat,
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
pub fn io_algo_4_bpap_mat(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_pos_dist_4_el: usize, max_substr_dist: usize) -> Prob4dMat {
  let sta_inside_part_func_mat_sets = get_sta_inside_part_func_mat_sets(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_pos_dist_4_il, max_pos_dist_4_el, max_substr_dist);
  let sta_outside_part_func_4d_mats = get_sta_outside_part_func_4d_mats(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_pos_dist_4_il, max_substr_dist, &sta_inside_part_func_mat_sets);
  get_bpap_mat(seq_len_pair, &sta_inside_part_func_mat_sets, &sta_outside_part_func_4d_mats, sta_fe_params)
}

#[inline]
pub fn get_sta_inside_part_func_mat_sets(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_pos_dist_4_el: usize, max_substr_dist: usize) -> StaInsidePartFuncMatSets {
  let mut sta_inside_part_func_mat_sets = StaInsidePartFuncMatSets::new();
  let invert_exp_max_free_energy = sta_fe_params.invert_exp_max_free_energy;
  for i in 1 .. seq_len_pair.0 {
    for j in 1 .. seq_len_pair.1 {
      let pos_pair = (i, j);
      let pos_quadruple = (i, i - 1, j, j - 1);
      if !is_pos_dist_ok(&pos_pair, max_pos_dist_4_il) {continue;}
      sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_quadruple, invert_exp_max_free_energy);
      sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas.insert(pos_quadruple, invert_exp_max_free_energy);
      sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas.insert(pos_quadruple, invert_exp_max_free_energy);
    }
  }
  let leftmost_pos_pair = (0, 0);
  let rightmost_pos_pair = (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
  sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat.insert(leftmost_pos_pair, invert_exp_max_free_energy);
  sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas.insert(leftmost_pos_pair, invert_exp_max_free_energy);
  sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat.insert(rightmost_pos_pair, invert_exp_max_free_energy);
  sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas.insert(rightmost_pos_pair, invert_exp_max_free_energy);
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
          let mut sum = 0.;
          if sta_fe_params.exp_bpa_score_mat.contains_key(&pos_quadruple) {
            let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple];
            if sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(i + 1, j - 1, k + 1, l - 1)) {
              sum += exp_bpa_score * get_exp_hl_fe(seq_pair.0, &(i, j)) * get_exp_hl_fe(seq_pair.1, &(k, l)) * sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(i + 1, j - 1, k + 1, l - 1)];
            }
            for m in i + 1 .. j - 1 {
              for n in m + 1 .. j {
                if j - n - 1 + m - i - 1 > MAX_2_LOOP_LEN {continue;}
                for o in k + 1 .. l - 1 {
                  for p in o + 1 .. l {
                    if l - p - 1 + o - k - 1 > MAX_2_LOOP_LEN {continue;}
                    let pos_quadruple_2 = (m, n, o, p);
                    if !sta_inside_part_func_mat_sets.part_func_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
                    if !(sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(i + 1, m - 1, k + 1, o - 1)) && sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(n + 1, j - 1, p + 1, l - 1))) {continue;}
                    let exp_2loop_fe_pair = (
                      get_exp_2_loop_fe(seq_pair.0, &(i, j), &(m , n)),
                      get_exp_2_loop_fe(seq_pair.1, &(k, l), &(o , p)),
                    );
                    sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(i + 1, m - 1, k + 1, o - 1)] / invert_exp_max_free_energy * exp_bpa_score * exp_2loop_fe_pair.0 * exp_2loop_fe_pair.1 * sta_inside_part_func_mat_sets.part_func_mat_4_bpas[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(n + 1, j - 1, p + 1, l - 1)] / invert_exp_max_free_energy;
                  }
                }
              }
            }
            if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(i + 1, j - 1, k + 1, l - 1)) {
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
              sum += exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe_pair.0 * exp_ml_tm_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1 * sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(i + 1, j - 1, k + 1, l - 1)];
            }
            if sum > 0. {
              sta_inside_part_func_mat_sets.part_func_mat_4_bpas.insert(pos_quadruple, sum);
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
              sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_els.insert(pos_quadruple, sum);
              sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, sum * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE));
            }
          }
          // Compute "sta_inside_part_func_mat_sets.part_func_mats_on_sa"
          if sta_fe_params.exp_ba_score_mat.contains_key(&(j, l)) {
            let exp_ba_score = sta_fe_params.exp_ba_score_mat[&(j, l)];
            if sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(i, j - 1, k, l - 1)) {
              sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas.insert(pos_quadruple,
                sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(i, j - 1, k, l - 1)] * exp_ba_score
              );
            }
          }
          let pos_quadruple_2 = (i, j - 1, k, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
            if sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
              sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_1[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
            }
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_1.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_2[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_2.insert(pos_quadruple, sum);
          }
          if sta_fe_params.exp_ba_score_mat.contains_key(&(i, k)) {
            let exp_ba_score = sta_fe_params.exp_ba_score_mat[&(i, k)];
            if sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(i + 1, j, k + 1, l)) {
              sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas.insert(pos_quadruple,
                sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(i + 1, j, k + 1, l)] * exp_ba_score
              );
            }
          }
          let pos_quadruple_2 = (i + 1, j, k, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_1[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_1.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_2[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_2.insert(pos_quadruple, sum);
          }
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_bas[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_1[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.forward_part_func_mat_4_gaps_2[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_bas[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_1[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.backward_part_func_mat_4_gaps_2[&pos_quadruple];
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.insert(pos_quadruple, sum);
          }
          // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop"
          sum = 0.;
          for m in i .. j {
            for n in k .. l {
              let pos_quadruple_2 = (i, m - 1, k, n - 1);
              let pos_quadruple_3 = (m, j, n, l);
              if !sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_3) {
                continue;
              }
              if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&pos_quadruple_2) {
                sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple_3] / invert_exp_max_free_energy;
              }
              if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&pos_quadruple_2) {
                sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple_3] / invert_exp_max_free_energy;
              }
            }
          }
          if sta_fe_params.exp_ba_score_mat.contains_key(&(j, l)) {
            let exp_ba_score = sta_fe_params.exp_ba_score_mat[&(j, l)];
            if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(i, j - 1, k, l - 1)) {
              sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(i, j - 1, k, l - 1)] * exp_ba_score;
            }
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_bas.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j - 1, k, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_1[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_1.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_2[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_2.insert(pos_quadruple, sum);
          }
          sum = 0.;
          let mut sum_2 = sum;
          for m in i + 1 .. j + 1 {
            for n in k + 1 .. l + 1 {
              let pos_quadruple_2 = (i, m, k, n);
              let pos_quadruple_3 = (m + 1, j, n + 1, l);
              if !sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2) {
                continue;
              }
              if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&pos_quadruple_3) {
                let new_term = sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&pos_quadruple_3] / invert_exp_max_free_energy;
                sum += new_term;
                if m != j && n != l {
                  sum_2 += new_term;
                }
              }
              if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&pos_quadruple_3) {
                let new_term = sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&pos_quadruple_3] / invert_exp_max_free_energy;
                sum += new_term;
                if m != j && n != l {
                  sum_2 += sum;
                }
              }
            }
          }
          if sta_fe_params.exp_ba_score_mat.contains_key(&(i, k)) {
            let exp_ba_score = sta_fe_params.exp_ba_score_mat[&(i, k)];
            if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(i + 1, j, k + 1, l)) {
              sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(i + 1, j, k + 1, l)] * exp_ba_score;
            }
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_bas.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i + 1, j, k, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_1[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_1.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_2[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_2.insert(pos_quadruple, sum);
          }
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_bas.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_bas[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_1[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.forward_part_func_mat_4_gaps_2[&pos_quadruple];
          }
          if sum_2 > 0. {
            sum += sum_2;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_1[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.backward_part_func_mat_4_gaps_2[&pos_quadruple];
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.insert(pos_quadruple, sum);
          }
          // Compute "sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls"
          sum = 0.;
          for m in i .. j {
            for n in k .. l {
              let pos_quadruple_2 = (i, m - 1, k, n - 1);
              let pos_quadruple_3 = (m, j, n, l);
              if !(sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_quadruple_2) && sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_3)) {
                continue;
              }
              sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple_3] / invert_exp_max_free_energy;
            }
          }
          if sta_fe_params.exp_ba_score_mat.contains_key(&(j, l)) {
            let exp_ba_score = sta_fe_params.exp_ba_score_mat[&(j, l)];
            if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(i, j - 1, k, l - 1)) {
              sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(i, j - 1, k, l - 1)] * exp_ba_score;
            }
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_bas.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j - 1, k, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_1[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_1.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j, k, l - 1);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_2[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_2.insert(pos_quadruple, sum);
          }
          sum = 0.;
          sum_2 = sum;
          for m in i + 1 .. j + 1 {
            for n in k + 1 .. l + 1 {
              let pos_quadruple_2 = (i, m, k, n);
              let pos_quadruple_3 = (m + 1, j, n + 1, l);
              if !(sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&pos_quadruple_3) && sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls.contains_key(&pos_quadruple_2)) {
                continue;
              }
              let new_term = sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&pos_quadruple_3] / invert_exp_max_free_energy;
              sum += new_term;
              if m != j && n != l {
                sum_2 += new_term;
              }
            }
          }
          if sta_fe_params.exp_ba_score_mat.contains_key(&(i, k)) {
            let exp_ba_score = sta_fe_params.exp_ba_score_mat[&(i, k)];
            if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(i + 1, j, k + 1, l)) {
              sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(i + 1, j, k + 1, l)] * exp_ba_score;
            }
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_bas.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i + 1, j, k, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_1[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_1.insert(pos_quadruple, sum);
          }
          let pos_quadruple_2 = (i, j, k + 1, l);
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_bas.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_bas[&pos_quadruple_2] * sta_fe_params.exp_opening_gap_penalty;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple_2) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_2[&pos_quadruple_2] * sta_fe_params.exp_extending_gap_penalty;
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_2.insert(pos_quadruple, sum);
          }
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_bas.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_bas[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_1[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.forward_part_func_mat_4_gaps_2[&pos_quadruple];
          }
          if sum_2 > 0. {
            sum += sum_2;
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_1.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_1[&pos_quadruple];
          }
          if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_2.contains_key(&pos_quadruple) {
            sum += sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.backward_part_func_mat_4_gaps_2[&pos_quadruple];
          }
          if sum > 0. {
            sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.insert(pos_quadruple, sum);
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
      let mut sum = 0.;
      for k in 1 .. i {
        for l in 1 .. j {
          let pos_pair_2 = (k - 1, l - 1);
          let pos_quadruple = (k, i, l, j);
          if !(sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat.contains_key(&pos_pair_2) && sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_els.contains_key(&pos_quadruple)) {
            continue;
          }
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_els[&pos_quadruple] / invert_exp_max_free_energy;
        }
      }
      if i > 0 && j > 0 && sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat.contains_key(&(i - 1, j - 1)) {
        let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat[&(i - 1, j - 1)] * exp_ba_score;
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas.insert(pos_pair, sum);
      }
      sum = 0.;
      if i > 0 {
        let pos_pair_2 = (i - 1, j);
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_1.insert(pos_pair, sum);
      }
      sum = 0.;
      if j > 0 {
        let pos_pair_2 = (i, j - 1);
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_2.insert(pos_pair, sum);
      }
      sum = 0.;
      if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_bas[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_1[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat.insert(pos_pair, sum);
      }
    }
  }
  for i in (1 .. seq_len_pair.0).rev() {
    for j in (1 .. seq_len_pair.1).rev() {
      let pos_pair = (i, j);
      if pos_pair == (seq_len_pair.0 - 1, seq_len_pair.1 - 1) {continue;}
      if !is_pos_dist_ok(&pos_pair, max_pos_dist_4_el) {continue;}
      let mut sum = 0.;
      for k in i + 1 .. seq_len_pair.0 - 1 {
        for l in j + 1 .. seq_len_pair.1 - 1 {
          let pos_pair_2 = (k + 1, l + 1);
          let pos_quadruple = (i, k, j, l);
          if !(sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat.contains_key(&pos_pair_2) && sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_els.contains_key(&pos_quadruple)) {
            continue;
          }
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat[&pos_pair_2] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_els[&pos_quadruple] / invert_exp_max_free_energy;
        }
      }
      if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 && sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat.contains_key(&(i + 1, j + 1)) {
        let exp_ba_score = sta_fe_params.exp_ba_score_mat[&pos_pair];
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat[&(i + 1, j + 1)] * exp_ba_score;
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas.insert(pos_pair, sum);
      }
      sum = 0.;
      if i < seq_len_pair.0 - 1 {
        let pos_pair_2 = (i + 1, j);
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_1.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_1[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_1.insert(pos_pair, sum);
      }
      sum = 0.;
      if j < seq_len_pair.1 - 1 {
        let pos_pair_2 = (i, j + 1);
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas[&pos_pair_2] * sta_fe_params.exp_opening_gap_penalty;
        }
        if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_2.contains_key(&pos_pair_2) {
          sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_2[&pos_pair_2] * sta_fe_params.exp_extending_gap_penalty;
        }
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_2.insert(pos_pair, sum);
      }
      sum = 0.;
      if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_bas[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_1.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_1[&pos_pair];
      }
      if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_2.contains_key(&pos_pair) {
        sum += sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat_4_gaps_2[&pos_pair];
      }
      if sum > 0. {
        sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat.insert(pos_pair, sum);
      }
    }
  }
  sta_inside_part_func_mat_sets
}

#[inline]
pub fn get_sta_outside_part_func_4d_mats(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &(usize, usize), max_pos_dist_4_il: usize, max_substr_dist: usize, sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets) -> StaOutsidePartFunc4dMats {
  let mut sta_outside_part_func_4d_mats = StaOutsidePartFunc4dMats::new();
  let invert_exp_max_free_energy = sta_fe_params.invert_exp_max_free_energy;
  for substr_len_1 in (2 .. max_bp_span_pair.0 + 1).rev() {
    for substr_len_2 in (2 .. max_bp_span_pair.1 + 1).rev() {
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_substr_dist {continue;}
      for i in 1 .. seq_len_pair.0 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 1 .. seq_len_pair.1 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !is_substr_dist_ok(&pos_quadruple, max_substr_dist, max_pos_dist_4_il) {continue;}
          // Compute "sta_outside_part_func_mat_sets.part_func_mat_4_bpas"
          let mut sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mat_4_bpas.contains_key(&pos_quadruple) {
            if sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat.contains_key(&(i - 1, k - 1)) && sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat.contains_key(&(j + 1, l + 1)) {
              sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_el.insert(pos_quadruple, sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat[&(i - 1, k - 1)] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_els[&pos_quadruple] / sta_inside_part_func_mat_sets.part_func_mat_4_bpas[&pos_quadruple] * sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.backward_part_func_mat[&(j + 1, l + 1)] / invert_exp_max_free_energy);
            }
            for m in 1 .. i {
              for n in j + 1 .. seq_len_pair.0 - 1 {
                if n - j - 1 + i - m - 1 > MAX_2_LOOP_LEN {continue;}
                for o in 1 .. k {
                  for p in l + 1 .. seq_len_pair.1 - 1 {
                    if p - l - 1 + k - o - 1 > MAX_2_LOOP_LEN {continue;}
                    let pos_quadruple_2 = (m, n, o, p);
                    if !sta_outside_part_func_4d_mats.part_func_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
                    if !(sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(m + 1, i - 1, o + 1, k - 1)) && sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(j + 1, n - 1, l + 1, p - 1))) {continue;}
                    let exp_bpa_score = sta_fe_params.exp_bpa_score_mat[&pos_quadruple_2];
                    let exp_2loop_fe_pair = (
                      get_exp_2_loop_fe(seq_pair.0, &(m, n), &(i , j)) as FreeEnergy,
                      get_exp_2_loop_fe(seq_pair.1, &(o, p), &(k , l)) as FreeEnergy,
                    );
                    sum += sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(m + 1, i - 1, o + 1, k - 1)] * exp_bpa_score * exp_2loop_fe_pair.0 * exp_2loop_fe_pair.1 * sta_outside_part_func_4d_mats.part_func_mat_4_bpas[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(j + 1, n - 1, l + 1, p - 1)] / invert_exp_max_free_energy;
                  }
                }
              }
            }
            if sum > 0. {
              sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_internal_2loops.insert(pos_quadruple, sum);
            }
            sum = 0.;
            for m in j + 1 .. seq_len_pair.0 - 1 {
              for n in l + 1 .. seq_len_pair.1 - 1 {
                let pos_quadruple_2 = (i, m, k, n);
                if !sta_outside_part_func_4d_mats.part_func_mat_4_bpas.contains_key(&pos_quadruple_2) {continue;}
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
                let coefficient = exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe_pair.0 * exp_ml_tm_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1;
                if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(j + 1, m - 1, l + 1, n - 1)) {
                  sum += coefficient * sta_outside_part_func_4d_mats.part_func_mat_4_bpas[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(j + 1, m - 1, l + 1, n - 1)] / invert_exp_max_free_energy;
                }
                if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(j + 1, m - 1, l + 1, n - 1)) {
                  sum += coefficient * sta_outside_part_func_4d_mats.part_func_mat_4_bpas[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(j + 1, m - 1, l + 1, n - 1)] / invert_exp_max_free_energy;
                }
              }
            }
            if sum > 0. {
              sta_outside_part_func_4d_mats.part_func_mat_4_right_2.insert(pos_quadruple, sum);
            }
          }
          for m in j + 1 .. seq_len_pair.0 - 1 {
            for n in l + 1 .. seq_len_pair.1 - 1 {
              let pos_quadruple_2 = (i, m, k, n);
              if !(sta_outside_part_func_4d_mats.part_func_mat_4_bpas.contains_key(&pos_quadruple_2) && sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(j + 1, m - 1, l + 1, n - 1))) {continue;}
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
              let coefficient = exp_bpa_score * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * (*EXP_CONST_4_INIT_ML_DELTA_FE) * (*EXP_COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE) * exp_ml_tm_delta_fe_pair.0 * exp_ml_tm_delta_fe_pair.1 * exp_au_or_gu_end_penalty_delta_fe_pair.0 * exp_au_or_gu_end_penalty_delta_fe_pair.1;
              sum += coefficient * sta_outside_part_func_4d_mats.part_func_mat_4_bpas[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(j + 1, m - 1, l + 1, n - 1)] / invert_exp_max_free_energy;
            }
          }
          if sum > 0. {
            sta_outside_part_func_4d_mats.part_func_mat_4_right.insert(pos_quadruple, sum);
          }
          sum = 0.;
          if sta_inside_part_func_mat_sets.part_func_mat_4_bpas.contains_key(&pos_quadruple) {
            for m in 1 .. i {
              for n in 1 .. k {
                let pos_quadruple_2 = (m, j, n, l);
                if sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat.contains_key(&(m + 1, i - 1, n + 1, k - 1)) && sta_outside_part_func_4d_mats.part_func_mat_4_right.contains_key(&pos_quadruple_2) {
                  sum += sta_outside_part_func_4d_mats.part_func_mat_4_right[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_4_internal_multiloop.part_func_mat[&(m + 1, i - 1, n + 1, k - 1)] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple] / sta_inside_part_func_mat_sets.part_func_mat_4_bpas[&pos_quadruple] / invert_exp_max_free_energy;
                }
                if sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat.contains_key(&(m + 1, i - 1, n + 1, k - 1)) && sta_outside_part_func_4d_mats.part_func_mat_4_right.contains_key(&pos_quadruple_2) {
                  sum += sta_outside_part_func_4d_mats.part_func_mat_4_right[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_4_first_bpas_on_mls.part_func_mat[&(m + 1, i - 1, n + 1, k - 1)] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple] / sta_inside_part_func_mat_sets.part_func_mat_4_bpas[&pos_quadruple] / invert_exp_max_free_energy;
                }
                if sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat.contains_key(&(m + 1, i - 1, n + 1, k - 1)) && sta_outside_part_func_4d_mats.part_func_mat_4_right_2.contains_key(&pos_quadruple_2) {
                  sum += sta_outside_part_func_4d_mats.part_func_mat_4_right_2[&pos_quadruple_2] * sta_inside_part_func_mat_sets.part_func_mats_on_sa.part_func_mat[&(m + 1, i - 1, n + 1, k - 1)] * sta_inside_part_func_mat_sets.part_func_mat_4_bpas_accessible_on_mls[&pos_quadruple] / sta_inside_part_func_mat_sets.part_func_mat_4_bpas[&pos_quadruple] / invert_exp_max_free_energy;
                }
              }
            }
            if sum > 0. {
              sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_internal_mls.insert(pos_quadruple, sum);
            }
            sum = 0.;
            if sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_el.contains_key(&pos_quadruple) {
              sum += sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_el[&pos_quadruple];
            }
            if sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_internal_2loops.contains_key(&pos_quadruple) {
              sum += sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_internal_2loops[&pos_quadruple];
            }
            if sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_internal_mls.contains_key(&pos_quadruple) {
              sum += sta_outside_part_func_4d_mats.part_func_mat_4_bpas_on_internal_mls[&pos_quadruple];
            }
            if sum > 0. {
              sta_outside_part_func_4d_mats.part_func_mat_4_bpas.insert(pos_quadruple, sum);
            }
          }
        }
      }
    }
  }
  sta_outside_part_func_4d_mats
}

#[inline]
fn get_bpap_mat(seq_len_pair: &(usize, usize), sta_inside_part_func_mat_sets: &StaInsidePartFuncMatSets, sta_outside_part_func_4d_mats: &StaOutsidePartFunc4dMats, sta_fe_params: &StaFeParams) -> Prob4dMat {
  let mut bpap_mat = Prob4dMat::default();
  let invert_exp_max_free_energy = sta_fe_params.invert_exp_max_free_energy;
  let part_func = sta_inside_part_func_mat_sets.part_func_mats_4_external_loop.forward_part_func_mat[&(seq_len_pair.0 - 2, seq_len_pair.1 - 2)];
  for pos_quadruple in sta_inside_part_func_mat_sets.part_func_mat_4_bpas.keys() {
    if !sta_outside_part_func_4d_mats.part_func_mat_4_bpas.contains_key(pos_quadruple) {continue;}
    let bpap = sta_inside_part_func_mat_sets.part_func_mat_4_bpas[pos_quadruple] * sta_outside_part_func_4d_mats.part_func_mat_4_bpas[pos_quadruple] / (part_func * invert_exp_max_free_energy);
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
      if pos_quadruple.0 == 0 {continue;}
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
  for (pos_pair, &bpp) in bpp_mat {
    if !new_bpp_mat.contains_key(&pos_pair) {
      new_bpp_mat.insert(*pos_pair, bpp);
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
