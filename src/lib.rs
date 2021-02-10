extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate getopts;
extern crate scoped_threadpool;
extern crate itertools;
extern crate bio;
extern crate num_cpus;
extern crate hashbrown;

pub mod utils;

pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use self::scoped_threadpool::Pool;
pub use itertools::multizip;
pub use utils::*;
pub use std::path::Path;
pub use bio::io::fasta::Reader;
pub use std::io::prelude::*;
pub use std::io::BufWriter;
pub use std::fs::File;
pub use std::fs::create_dir;
pub use std::cmp::Ord;
pub use std::marker::{Sync, Send};
pub use hashbrown::HashSet;

pub type PosQuadrupleMat<T> = HashSet<PosQuadruple<T>>;
pub type PosPairMatSet<T> = HashMap<PosPair<T>, PosPairMat<T>>;
pub type PosPairMat<T> = HashSet<PosPair<T>>;
pub type Prob4dMat<T> = HashMap<PosQuadruple<T>, Prob>;
pub type PartFunc4dMat<T> = HashMap<PosQuadruple<T>, PartFunc>;
pub type TmpPartFuncSetMat<T> = HashMap<PosPair<T>, TmpPartFuncSets>;
#[derive(Clone)]
pub struct TmpPartFuncSets {
  pub part_funcs_on_sa: TmpPartFuncs,
  pub part_funcs_4_ml: TmpPartFuncs,
  pub part_funcs_4_first_bpas_on_mls: TmpPartFuncs,
  pub part_funcs_4_bpas_on_mls: TmpPartFuncs,
  pub part_funcs_on_mls: TmpPartFuncs,
}
#[derive(Clone)]
pub struct TmpPartFuncs {
  pub part_func_4_align: PartFunc,
  pub part_func_4_insert: PartFunc,
  pub part_func_4_insert_2: PartFunc,
}
pub type TmpPartFuncSetMat4El<T> = HashMap<PosPair<T>, TmpPartFuncs>;
#[derive(Clone)]
pub struct StaPartFuncMats<T> {
  pub part_func_4d_mat_4_bpas: PartFunc4dMat<T>,
  pub part_func_4d_mat_4_bpas_accessible_on_els: PartFunc4dMat<T>,
  pub part_func_4d_mat_4_bpas_accessible_on_mls: PartFunc4dMat<T>,
  pub forward_part_func_set_mat_4_external_loop: TmpPartFuncSetMat4El<T>,
  pub backward_part_func_set_mat_4_external_loop: TmpPartFuncSetMat4El<T>,
  pub forward_tmp_part_func_set_mats_with_pos_pairs: TmpPartFuncSetMatsWithPosPairs<T>,
  pub backward_tmp_part_func_set_mats_with_pos_pairs: TmpPartFuncSetMatsWithPosPairs<T>,
}
pub struct StaFeParams<T> {
  pub ba_score_mat: SparseFreeEnergyMat<T>,
  pub bpa_score_mat: FreeEnergy4dMat<T>,
  pub insert_scores: FreeEnergies,
  pub insert_scores_2: FreeEnergies,
}
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type StaFeParamSetsWithRnaIdPairs<T> = HashMap<RnaIdPair, StaFeParams<T>>;
pub type Prob4dMatsWithRnaIdPairs<T> = HashMap<RnaIdPair, Prob4dMat<T>>;
pub type ProbMats<T> = Vec<SparseProbMat<T>>;
pub type Prob1dMats = Vec<Probs>;
pub type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
pub type SeqPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type FreeEnergyPair = (FreeEnergy, FreeEnergy);
pub type SparseFreeEnergyMat<T> = HashMap<PosPair<T>, FreeEnergy>;
pub type PosPairsWithPosPairs<T> = HashMap<PosPair<T>, PosPair<T>>;
pub type BoolsWithPosPairs<T> = HashMap<PosPair<T>, bool>;
pub type ProbMatPair<'a, T> =  (&'a SparseProbMat<T>, &'a SparseProbMat<T>);
pub type SsFreeEnergyMatSetPair<'a, T> =  (&'a SsFreeEnergyMats<T>, &'a SsFreeEnergyMats<T>);
pub type NumOfThreads = u32;
pub type FreeEnergySetPair<'a> = (&'a FreeEnergies, &'a FreeEnergies);
pub struct StaProbMats<T> {
  pub bpp_mat_pair: SparseProbMatPair<T>,
  pub access_bpp_mat_pair_4_2l: SparseProbMatPair<T>,
  pub access_bpp_mat_pair_4_ml: SparseProbMatPair<T>,
  pub bpp_mat_pair_4_el: SparseProbMatPair<T>,
  pub upp_mat_pair: ProbSetPair,
  pub upp_mat_pair_4_hl: ProbSetPair,
  pub upp_mat_pair_4_2l: ProbSetPair,
  pub upp_mat_pair_4_ml: ProbSetPair,
  pub upp_mat_pair_4_el: ProbSetPair,
}
#[derive(Clone)]
pub struct PctStaProbMats<T> {
  pub max_bpp_mat: SparseProbMat<T>,
  pub bpp_mat: SparseProbMat<T>,
  pub access_bpp_mat_4_2l: SparseProbMat<T>,
  pub access_bpp_mat_4_ml: SparseProbMat<T>,
  pub bpp_mat_4_el: SparseProbMat<T>,
  pub bpp_mat_on_ss: SparseProbMat<T>,
  pub max_upp_mat: Probs,
  pub upp_mat: Probs,
  pub upp_mat_4_hl: Probs,
  pub upp_mat_4_2l: Probs,
  pub upp_mat_4_ml: Probs,
  pub upp_mat_4_el: Probs,
}
pub type SparseProbMatPair<T> = (SparseProbMat<T>, SparseProbMat<T>);
pub type ProbSetPair = (Probs, Probs);
pub type ProbMatSets<T> = Vec<PctStaProbMats<T>>;
pub type StaProbMatsWithRnaIdPairs<T> = HashMap<RnaIdPair, StaProbMats<T>>;
pub type ProbSeqPair<'a> = (&'a Probs, &'a Probs);
pub type Poss<T> = Vec<T>;
pub type TmpPartFuncSetMatsWithPosPairs<T> = HashMap<PosPair<T>, TmpPartFuncSetMat<T>>;

impl<T: Hash + ToPrimitive + Clone> StaProbMats<T> {
  pub fn origin() -> StaProbMats<T> {
    let prob_mat_pair = (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    let prob_set_pair = (Vec::new(), Vec::new());
    StaProbMats {
      bpp_mat_pair: prob_mat_pair.clone(),
      access_bpp_mat_pair_4_2l: prob_mat_pair.clone(),
      access_bpp_mat_pair_4_ml: prob_mat_pair.clone(),
      bpp_mat_pair_4_el: prob_mat_pair,
      upp_mat_pair: prob_set_pair.clone(),
      upp_mat_pair_4_hl: prob_set_pair.clone(),
      upp_mat_pair_4_2l: prob_set_pair.clone(),
      upp_mat_pair_4_ml: prob_set_pair.clone(),
      upp_mat_pair_4_el: prob_set_pair,
    }
  }
  pub fn new(seq_len_pair: &PosPair<T>) -> StaProbMats<T> {
    let prob_mat_pair = (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    let prob_set_pair = (vec![NEG_INFINITY; seq_len_pair.0.to_usize().unwrap()], vec![NEG_INFINITY; seq_len_pair.1.to_usize().unwrap()]);
    StaProbMats {
      bpp_mat_pair: prob_mat_pair.clone(),
      access_bpp_mat_pair_4_2l: prob_mat_pair.clone(),
      access_bpp_mat_pair_4_ml: prob_mat_pair.clone(),
      bpp_mat_pair_4_el: prob_mat_pair,
      upp_mat_pair: (vec![1.; seq_len_pair.0.to_usize().unwrap()], vec![1.; seq_len_pair.1.to_usize().unwrap()]),
      upp_mat_pair_4_hl: prob_set_pair.clone(),
      upp_mat_pair_4_2l: prob_set_pair.clone(),
      upp_mat_pair_4_ml: prob_set_pair.clone(),
      upp_mat_pair_4_el: prob_set_pair,
    }
  }
}

impl<T: Hash + Clone> PctStaProbMats<T> {
  pub fn origin() -> PctStaProbMats<T> {
    let prob_mat = SparseProbMat::<T>::default();
    let probs = Vec::new();
    PctStaProbMats {
      max_bpp_mat: prob_mat.clone(),
      bpp_mat: prob_mat.clone(),
      access_bpp_mat_4_2l: prob_mat.clone(),
      access_bpp_mat_4_ml: prob_mat.clone(),
      bpp_mat_4_el: prob_mat.clone(),
      bpp_mat_on_ss: prob_mat,
      max_upp_mat: probs.clone(),
      upp_mat: probs.clone(),
      upp_mat_4_hl: probs.clone(),
      upp_mat_4_2l: probs.clone(),
      upp_mat_4_ml: probs.clone(),
      upp_mat_4_el: probs,
    }
  }
  pub fn new(seq_len: usize) -> PctStaProbMats<T> {
    let prob_mat = SparseProbMat::<T>::default();
    let probs = vec![0.; seq_len as usize];
    PctStaProbMats {
      max_bpp_mat: prob_mat.clone(),
      bpp_mat: prob_mat.clone(),
      access_bpp_mat_4_2l: prob_mat.clone(),
      access_bpp_mat_4_ml: prob_mat.clone(),
      bpp_mat_4_el: prob_mat.clone(),
      bpp_mat_on_ss: prob_mat,
      max_upp_mat: probs.clone(),
      upp_mat: probs.clone(),
      upp_mat_4_hl: probs.clone(),
      upp_mat_4_2l: probs.clone(),
      upp_mat_4_ml: probs.clone(),
      upp_mat_4_el: probs,
    }
  }
}

impl<T: Hash + Eq + Integer + FromPrimitive + PrimInt + Unsigned> StaFeParams<T> {
  pub fn origin() -> StaFeParams<T> {
    let sparse_free_energy_mat = SparseFreeEnergyMat::<T>::default();
    let free_energy_4d_mat = FreeEnergy4dMat::<T>::default();
    let insert_scores = Vec::new();
    StaFeParams {
      ba_score_mat: sparse_free_energy_mat.clone(),
      bpa_score_mat: free_energy_4d_mat.clone(),
      insert_scores: insert_scores.clone(),
      insert_scores_2: insert_scores,
    }
  }

  pub fn new(seq_pair: &SeqPair, seq_len_pair: &PosPair<T>, max: T, pos_quadruple_mat: &PosQuadrupleMat<T>,) -> StaFeParams<T> {
    let mut sta_fe_params = StaFeParams::<T>::origin();
    sta_fe_params.insert_scores = vec![0.; seq_len_pair.0.to_usize().unwrap()];
    sta_fe_params.insert_scores_2 = vec![0.; seq_len_pair.1.to_usize().unwrap()];
    for j in range(T::one(), seq_len_pair.1 - T::one()) {
      let long_j = j.to_usize().unwrap();
      let base = seq_pair.1[long_j];
      sta_fe_params.insert_scores_2[long_j] = INSERT_SCORES[base];
    }
    let pseudo_pos_quadruple = (T::zero(), seq_len_pair.0 - T::one(), T::zero(), seq_len_pair.1 - T::one());
    for i in range(T::one(), seq_len_pair.0 - T::one()) {
      let long_i = i.to_usize().unwrap();
      let base = seq_pair.0[long_i];
      sta_fe_params.insert_scores[long_i] = INSERT_SCORES[base];
      for j in range(T::one(), seq_len_pair.1 - T::one()) {
        let long_j = j.to_usize().unwrap();
        let pos_pair = (i, j);
        if !is_min_gap_ok(&pos_pair, &pseudo_pos_quadruple, max) {continue;}
        let base_pair = (base, seq_pair.1[long_j]);
        sta_fe_params.ba_score_mat.insert(pos_pair, BA_SCORE_MAT[&base_pair] + RIBOSUM_BA_SCORE_MAT[&base_pair]);
      }
    }
    for &(i, j, k, l) in pos_quadruple_mat {
      let (long_i, long_j, long_k, long_l) = (i.to_usize().unwrap(), j.to_usize().unwrap(), k.to_usize().unwrap(), l.to_usize().unwrap());
      let base_pair = (seq_pair.0[long_i], seq_pair.1[long_k]);
      let base_pair_2 = (seq_pair.0[long_j], seq_pair.1[long_l]);
      let base_pair_3 = (seq_pair.0[long_i], seq_pair.0[long_j]);
      let base_pair_4 = (seq_pair.1[long_k], seq_pair.1[long_l]);
      let pos_quadruple = (i, j, k, l);
      let align_score = BA_SCORE_MAT[&base_pair];
      sta_fe_params.bpa_score_mat.insert(pos_quadruple, align_score + BA_SCORE_MAT[&base_pair_2] + RIBOSUM_BPA_SCORE_MAT[&(base_pair_3, base_pair_4)]);
    }
    sta_fe_params
  }
}

impl TmpPartFuncSets {
  pub fn new() -> TmpPartFuncSets {
    let part_funcs = TmpPartFuncs::new();
    TmpPartFuncSets {
      part_funcs_on_sa: part_funcs.clone(),
      part_funcs_4_ml: part_funcs.clone(),
      part_funcs_4_first_bpas_on_mls: part_funcs.clone(),
      part_funcs_4_bpas_on_mls: part_funcs.clone(),
      part_funcs_on_mls: part_funcs,
    }
  }
}

impl TmpPartFuncs {
  pub fn new() -> TmpPartFuncs {
    TmpPartFuncs {
      part_func_4_align: NEG_INFINITY,
      part_func_4_insert: NEG_INFINITY,
      part_func_4_insert_2: NEG_INFINITY,
    }
  }
}

impl<T: Hash + Clone> StaPartFuncMats<T> {
  pub fn new() -> StaPartFuncMats<T> {
    let part_func_4d_mat = PartFunc4dMat::<T>::default();
    let part_func_set_mat = TmpPartFuncSetMat4El::<T>::new();
    let tmp_part_func_set_mats_with_pos_pairs = TmpPartFuncSetMatsWithPosPairs::<T>::default();
    StaPartFuncMats {
      part_func_4d_mat_4_bpas: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_els: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_mls: part_func_4d_mat,
      forward_part_func_set_mat_4_external_loop: part_func_set_mat.clone(),
      backward_part_func_set_mat_4_external_loop: part_func_set_mat,
      forward_tmp_part_func_set_mats_with_pos_pairs: tmp_part_func_set_mats_with_pos_pairs.clone(),
      backward_tmp_part_func_set_mats_with_pos_pairs: tmp_part_func_set_mats_with_pos_pairs,
    }
  }
}

pub const MAX_GAP_NUM_4_IL: usize = 15;
pub const MIN_GAP_NUM_4_IL: usize = 2;
pub const DEFAULT_MIN_BPP: Prob = 0.04;
pub const DEFAULT_OFFSET_4_MAX_GAP_NUM: usize = 1;
pub const BPP_MAT_FILE_NAME: &'static str = "bpp_mats.dat";
pub const MAX_BPP_MAT_FILE_NAME: &'static str = "max_bpp_mats.dat";
pub const ACCESS_BPP_MAT_ON_2L_FILE_NAME: &'static str = "access_bpp_mats_on_2l.dat";
pub const ACCESS_BPP_MAT_ON_ML_FILE_NAME: &'static str = "access_bpp_mats_on_ml.dat";
pub const BPP_MAT_ON_EL_FILE_NAME: &'static str = "bpp_mats_on_el.dat";
pub const BPP_MAT_ON_SS_FILE_NAME: &'static str = "bpp_mats_on_ss.dat";
pub const UPP_MAT_FILE_NAME: &'static str = "upp_mats.dat";
pub const MAX_UPP_MAT_FILE_NAME: &'static str = "max_upp_mats.dat";
pub const UPP_MAT_ON_HL_FILE_NAME: &'static str = "upp_mats_on_hl.dat";
pub const UPP_MAT_ON_2L_FILE_NAME: &'static str = "upp_mats_on_2l.dat";
pub const UPP_MAT_ON_ML_FILE_NAME: &'static str = "upp_mats_on_ml.dat";
pub const UPP_MAT_ON_EL_FILE_NAME: &'static str = "upp_mats_on_el.dat";

pub fn io_algo_4_prob_mats<T>(seq_pair: &SeqPair, seq_len_pair: &PosPair<T>, sta_fe_params: &StaFeParams<T>, max_bp_span_pair: &PosPair<T>, max_gap_num: T, max_gap_num_4_il: T, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair<T>, produces_access_probs: bool, forward_pos_pair_mat_set: &PosPairMatSet<T>, backward_pos_pair_mat_set: &PosPairMatSet<T>, pos_quadruple_mat: &PosQuadrupleMat<T>,) -> StaProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let (sta_part_func_mats, global_part_func) = get_sta_inside_part_func_mats::<T>(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, max_gap_num_4_il, ss_free_energy_mat_set_pair, forward_pos_pair_mat_set, backward_pos_pair_mat_set, pos_quadruple_mat,);
  get_sta_prob_mats::<T>(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, max_gap_num_4_il, &sta_part_func_mats, ss_free_energy_mat_set_pair, produces_access_probs, global_part_func, pos_quadruple_mat,)
}

pub fn get_sta_inside_part_func_mats<T>(seq_pair: &SeqPair, seq_len_pair: &PosPair<T>, sta_fe_params: &StaFeParams<T>, max_bp_span_pair: &PosPair<T>, max_gap_num: T, max_gap_num_4_il: T, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair<T>, forward_pos_pair_mat_set: &PosPairMatSet<T>, backward_pos_pair_mat_set: &PosPairMatSet<T>, pos_quadruple_mat: &PosQuadrupleMat<T>,) -> (StaPartFuncMats<T>, PartFunc)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let pseudo_pos_quadruple = (T::zero(), seq_len_pair.0 - T::one(), T::zero(), seq_len_pair.1 - T::one());
  let mut sta_part_func_mats = StaPartFuncMats::<T>::new();
  for substr_len_1 in range_inclusive(T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(), max_bp_span_pair.0) {
    for substr_len_2 in range_inclusive(T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(), max_bp_span_pair.1) {
      for &(i, j, k, l) in pos_quadruple_mat {
        if j - i + T::one() != substr_len_1 || l - k + T::one() != substr_len_2 {continue;}
        let (long_i, long_j, long_k, long_l) = (i.to_usize().unwrap(), j.to_usize().unwrap(), k.to_usize().unwrap(), l.to_usize().unwrap());
        let pos_quadruple = (i, j, k, l);
        match sta_fe_params.bpa_score_mat.get(&pos_quadruple) {
          Some(&bpa_score) => {
            let (forward_tmp_part_func_set_mat, part_func_on_sa, part_func_4_ml) = get_tmp_part_func_set_mat::<T>(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_part_func_mats, true, forward_pos_pair_mat_set, backward_pos_pair_mat_set);
            let (backward_tmp_part_func_set_mat, _, _) = get_tmp_part_func_set_mat::<T>(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_part_func_mats, false, forward_pos_pair_mat_set, backward_pos_pair_mat_set);
            let mut sum = NEG_INFINITY;
            let long_pos_pair = (long_i, long_j);
            let long_pos_pair_2 = (long_k, long_l);
            let score = bpa_score + ss_free_energy_mat_set_pair.0.hl_fe_mat[&(i, j)] + ss_free_energy_mat_set_pair.1.hl_fe_mat[&(k, l)] + part_func_on_sa;
            logsumexp(&mut sum, score);
            for &(m, n, o, p) in pos_quadruple_mat {
              if !(i < m && n < j) || !(k < o && p < l) {
                continue;
              }
              let (long_m, long_n, long_o, long_p) = (m.to_usize().unwrap(), n.to_usize().unwrap(), o.to_usize().unwrap(), p.to_usize().unwrap());
              if long_m - long_i - 1 + long_j - long_n - 1 > MAX_2_LOOP_LEN {continue;}
              if long_o - long_k - 1 + long_l - long_p - 1 > MAX_2_LOOP_LEN {continue;}
              let pos_quadruple_2 = (m, n, o, p);
              match sta_part_func_mats.part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                Some(&part_func) => {
                  let mut forward_term = NEG_INFINITY;
                  let mut backward_term = forward_term;
                  match forward_tmp_part_func_set_mat.get(&(m - T::one(), o - T::one())) {
                    Some(part_func_sets) => {
                      let ref part_funcs = part_func_sets.part_funcs_on_sa;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                      let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                      let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                    }, None => {},
                  }
                  match backward_tmp_part_func_set_mat.get(&(n + T::one(), p + T::one())) {
                    Some(part_func_sets) => {
                      let ref part_funcs = part_func_sets.part_funcs_on_sa;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut backward_term, term);
                      let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term, term);
                      let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term, term);
                    }, None => {},
                  }
                  let part_func_4_2l = forward_term + backward_term;
                  let twoloop_fe = ss_free_energy_mat_set_pair.0.twoloop_fe_4d_mat[&(i, j, m, n)];
                  let twoloop_fe_2 = ss_free_energy_mat_set_pair.1.twoloop_fe_4d_mat[&(k, l, o, p)];
                  let coefficient = bpa_score + twoloop_fe + twoloop_fe_2 + part_func;
                  logsumexp(&mut sum, coefficient + part_func_4_2l);
                }, None => {},
              }
            }
            let multi_loop_closing_basepairing_fe = get_ml_closing_basepairing_fe(
              seq_pair.0,
              &long_pos_pair,
            );
            let multi_loop_closing_basepairing_fe_2 =
              get_ml_closing_basepairing_fe(
                seq_pair.1,
                &long_pos_pair_2,
              );
            let score = bpa_score + multi_loop_closing_basepairing_fe + multi_loop_closing_basepairing_fe_2 + part_func_4_ml;
            logsumexp(&mut sum, score);
            if sum > NEG_INFINITY {
              sta_part_func_mats.part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
              let multi_loop_closing_basepairing_fe = get_ml_or_el_accessible_basepairing_fe(
                seq_pair.0,
                &long_pos_pair,
                true,
              );
              let multi_loop_closing_basepairing_fe_2 =
                get_ml_or_el_accessible_basepairing_fe(
                  seq_pair.1,
                  &long_pos_pair_2,
                  true,
                );
              sum += multi_loop_closing_basepairing_fe + multi_loop_closing_basepairing_fe_2;
              sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.insert(pos_quadruple, sum);
              sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, sum + 2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
            }
            match sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs.get_mut(&(i, k)) {
              Some(part_func_set_mat) => {
                *part_func_set_mat = forward_tmp_part_func_set_mat;
              },
              None => {
                sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs.insert((i, k), forward_tmp_part_func_set_mat);
              },
            }
            match sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs.get_mut(&(j, l)) {
              Some(part_func_set_mat) => {
                *part_func_set_mat = backward_tmp_part_func_set_mat;
              },
              None => {
                sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs.insert((j, l), backward_tmp_part_func_set_mat);
              },
            }
          }, None => {},
        }
      }
    }
  }
  let leftmost_pos_pair = (T::zero(), T::zero());
  let rightmost_pos_pair = (seq_len_pair.0 - T::one(), seq_len_pair.1 - T::one());
  let mut part_funcs = TmpPartFuncs::new();
  part_funcs.part_func_4_align = 0.;
  sta_part_func_mats.forward_part_func_set_mat_4_external_loop.insert(leftmost_pos_pair, part_funcs.clone());
  sta_part_func_mats.backward_part_func_set_mat_4_external_loop.insert(rightmost_pos_pair, part_funcs);
  for i in range(T::zero(), seq_len_pair.0 - T::one()) {
    let long_i = i.to_usize().unwrap();
    let insert_score = sta_fe_params.insert_scores[long_i];
    for j in range(T::zero(), seq_len_pair.1 - T::one()) {
      let pos_pair = (i, j);
      if pos_pair == (T::zero(), T::zero()) {continue;}
      if !is_min_gap_ok(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
      let mut part_funcs = TmpPartFuncs::new();
      let mut sum = NEG_INFINITY;
      match backward_pos_pair_mat_set.get(&pos_pair) {
        Some(backward_pos_pair_mat) => {
          for &(k, l) in backward_pos_pair_mat {
            let pos_pair_2 = (k - T::one(), l - T::one());
            let pos_quadruple = (k, i, l, j);
            let is_begin = pos_pair_2 == leftmost_pos_pair;
            match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
                  Some(part_funcs) => {
                    let score = part_funcs.part_func_4_align + part_func + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
                    logsumexp(&mut sum, score);
                    let score = part_funcs.part_func_4_insert + part_func + INSERT_2_MATCH_SCORE;
                    logsumexp(&mut sum, score);
                    let score = part_funcs.part_func_4_insert_2 + part_func + INSERT_2_MATCH_SCORE;
                    logsumexp(&mut sum, score);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      if i > T::zero() && j > T::zero() {
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        let pos_pair_2 = (i - T::one(), j - T::one());
        let is_begin = pos_pair_2 == leftmost_pos_pair;
        match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + ba_score + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert + ba_score + INSERT_2_MATCH_SCORE;
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert_2 + ba_score + INSERT_2_MATCH_SCORE;
            logsumexp(&mut sum, score);
          }, None => {},
        }
        part_funcs.part_func_4_align = sum;
      }
      if i > T::zero() {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i - T::one(), j);
        let is_begin = pos_pair_2 == leftmost_pos_pair;
        match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
            logsumexp(&mut sum, score);
          }, None => {},
        }
        let sum = sum + insert_score;
        part_funcs.part_func_4_insert = sum;
      }
      if j > T::zero() {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i, j - T::one());
        let is_begin = pos_pair_2 == leftmost_pos_pair;
        match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
            logsumexp(&mut sum, score);
          }, None => {},
        }
        let long_j = j.to_usize().unwrap();
        let insert_score_2 = sta_fe_params.insert_scores_2[long_j];
        let sum = sum + insert_score_2;
        part_funcs.part_func_4_insert_2 = sum;
      }
      sta_part_func_mats.forward_part_func_set_mat_4_external_loop.insert(pos_pair, part_funcs);
    }
  }
  let mut final_sum = NEG_INFINITY;
  let ref part_funcs = sta_part_func_mats.forward_part_func_set_mat_4_external_loop[&(seq_len_pair.0 - T::from_usize(2).unwrap(), seq_len_pair.1 - T::from_usize(2).unwrap())];
  logsumexp(&mut final_sum, part_funcs.part_func_4_align);
  logsumexp(&mut final_sum, part_funcs.part_func_4_insert);
  logsumexp(&mut final_sum, part_funcs.part_func_4_insert_2);
  for i in range(T::one(), seq_len_pair.0).rev() {
    let long_i = i.to_usize().unwrap();
    let insert_score = sta_fe_params.insert_scores[long_i];
    for j in range(T::one(), seq_len_pair.1).rev() {
      let pos_pair = (i, j);
      if pos_pair == (seq_len_pair.0 - T::one(), seq_len_pair.1 - T::one()) {continue;}
      if !is_min_gap_ok(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
      let mut part_funcs = TmpPartFuncs::new();
      let mut sum = NEG_INFINITY;
      match forward_pos_pair_mat_set.get(&pos_pair) {
        Some(forward_pos_pair_mat) => {
          for &(k, l) in forward_pos_pair_mat {
            let pos_pair_2 = (k + T::one(), l + T::one());
            let pos_quadruple = (i, k, j, l);
            let is_end = pos_pair_2 == rightmost_pos_pair;
            match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
                  Some(part_funcs) => {
                    let score = part_funcs.part_func_4_align + part_func + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                    logsumexp(&mut sum, score);
                    let score = part_funcs.part_func_4_insert + part_func + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                    logsumexp(&mut sum, score);
                    let score = part_funcs.part_func_4_insert_2 + part_func + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                    logsumexp(&mut sum, score);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      if i < seq_len_pair.0 - T::one() && j < seq_len_pair.1 - T::one() {
        let pos_pair_2 = (i + T::one(), j + T::one());
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        let is_end = pos_pair_2 == rightmost_pos_pair;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + ba_score + if is_end {0.} else {MATCH_2_MATCH_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert + ba_score + if is_end {0.} else {MATCH_2_INSERT_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert_2 + ba_score + if is_end {0.} else {MATCH_2_INSERT_SCORE};
            logsumexp(&mut sum, score);
          }, None => {},
        }
        part_funcs.part_func_4_align = sum;
      }
      if i < seq_len_pair.0 - T::one() {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i + T::one(), j);
        let is_end = pos_pair_2 == rightmost_pos_pair;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert + if is_end {0.} else {INSERT_EXTEND_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert_2 + if is_end {0.} else {INSERT_SWITCH_SCORE};
            logsumexp(&mut sum, score);
          }, None => {},
        }
        part_funcs.part_func_4_insert = sum + insert_score;
      }
      if j < seq_len_pair.1 - T::one() {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i, j + T::one());
        let is_end = pos_pair_2 == rightmost_pos_pair;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert + if is_end {0.} else {INSERT_SWITCH_SCORE};
            logsumexp(&mut sum, score);
            let score = part_funcs.part_func_4_insert_2 + if is_end {0.} else {INSERT_EXTEND_SCORE};
            logsumexp(&mut sum, score);
          }, None => {},
        }
        let long_j = j.to_usize().unwrap();
        let insert_score_2 = sta_fe_params.insert_scores_2[long_j];
        part_funcs.part_func_4_insert_2 = sum + insert_score_2;
      }
      sta_part_func_mats.backward_part_func_set_mat_4_external_loop.insert(pos_pair, part_funcs);
    }
  }
  (sta_part_func_mats, final_sum)
}

pub fn get_tmp_part_func_set_mat<T>(seq_len_pair: &PosPair<T>, sta_fe_params: &StaFeParams<T>, max_gap_num_4_il: T, pos_quadruple: &PosQuadruple<T>, sta_part_func_mats: &StaPartFuncMats<T>, is_forward: bool, forward_pos_pair_mat_set: &PosPairMatSet<T>, backward_pos_pair_mat_set: &PosPairMatSet<T>,) -> (TmpPartFuncSetMat<T>, PartFunc, PartFunc)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let pseudo_pos_quadruple = (T::zero(), seq_len_pair.0 - T::one(), T::zero(), seq_len_pair.1 - T::one());
  let mut cache_is_used = false;
  let &(i, j, k, l) = pos_quadruple;
  let leftmost_pos_pair = if is_forward {(i, k)} else {(i + T::one(), k + T::one())};
  let rightmost_pos_pair = if is_forward {(j - T::one(), l - T::one())} else {(j, l)};
  let tmp_part_func_set_mats_with_pos_pairs = if is_forward{&sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs} else {&sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs};
  let mut tmp_part_func_set_mat = match tmp_part_func_set_mats_with_pos_pairs.get(&if is_forward {leftmost_pos_pair} else {rightmost_pos_pair}) {
    Some(cache) => {
      cache_is_used = true;
      cache.clone()
    }, None => {
      TmpPartFuncSetMat::<T>::new()
    },
  };
  let iter: Poss<T> = if is_forward {range(i, j).collect()} else {range_inclusive(i + T::one(), j).rev().collect()};
  let iter_2: Poss<T> = if is_forward {range(k, l).collect()} else {range_inclusive(k + T::one(), l).rev().collect()};
  let pos_pair_mat_set = if is_forward {
    backward_pos_pair_mat_set
  } else {
    forward_pos_pair_mat_set
  };
  for &u in iter.iter() {
    let long_u = u.to_usize().unwrap();
    let insert_score = sta_fe_params.insert_scores[long_u];
    for &v in iter_2.iter() {
      let pos_pair = (u, v);
      if !is_min_gap_ok(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
      if cache_is_used && tmp_part_func_set_mat.contains_key(&pos_pair) {
        continue;
      }
      let mut tmp_part_func_sets = TmpPartFuncSets::new();
      if (is_forward && u == i && v == k) || (!is_forward && u == j && v == l) {
        tmp_part_func_sets.part_funcs_on_sa.part_func_4_align = 0.;
        tmp_part_func_set_mat.insert(pos_pair, tmp_part_func_sets);
        continue;
      }
      let long_v = v.to_usize().unwrap();
      let mut sum_on_sa = NEG_INFINITY;
      let mut sum_4_ml = sum_on_sa;
      let mut sum_4_first_bpas_on_mls = sum_on_sa;
      let mut tmp_sum = sum_on_sa;
      // For alignments.
      match pos_pair_mat_set.get(&pos_pair) {
        Some(pos_pair_mat) => {
          for &(m, n) in pos_pair_mat {
            if is_forward {
              if !(i < m && k < n) {continue;}
            } else {
              if !(m < j && n < l) {continue;}
            }
            let pos_pair_2 = if is_forward {(m - T::one(), n - T::one())} else {(m + T::one(), n + T::one())};
            let pos_quadruple_2 = if is_forward {(m, u, n, v)} else {(u, m, v, n)};
            match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
              Some(&part_func) => {
                let is_begin = pos_pair_2 == if is_forward {leftmost_pos_pair} else {rightmost_pos_pair};
                match tmp_part_func_set_mat.get(&pos_pair_2) {
                  Some(part_func_sets) => {
                    let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                    let score = part_funcs.part_func_4_align + part_func + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
                    logsumexp(&mut sum_4_ml, score);
                    let score = part_funcs.part_func_4_insert + part_func + if is_forward {
                      INSERT_2_MATCH_SCORE
                    } else {
                      MATCH_2_INSERT_SCORE
                    };
                    logsumexp(&mut sum_4_ml, score);
                    let score = part_funcs.part_func_4_insert_2 + part_func + if is_forward {
                      INSERT_2_MATCH_SCORE
                    } else {
                      MATCH_2_INSERT_SCORE
                    };
                    logsumexp(&mut sum_4_ml, score);
                    let ref part_funcs = part_func_sets.part_funcs_on_sa;
                    let score = part_funcs.part_func_4_align + part_func + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
                    logsumexp(&mut sum_4_first_bpas_on_mls, score);
                    let score = part_funcs.part_func_4_insert + part_func + if is_forward {
                      INSERT_2_MATCH_SCORE
                    } else {
                      MATCH_2_INSERT_SCORE
                    };
                    logsumexp(&mut sum_4_first_bpas_on_mls, score);
                    let score = part_funcs.part_func_4_insert_2 + part_func + if is_forward {
                      INSERT_2_MATCH_SCORE
                    } else {
                      MATCH_2_INSERT_SCORE
                    };
                    logsumexp(&mut sum_4_first_bpas_on_mls, score);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      let pos_pair_2 = if is_forward {(u - T::one(), v - T::one())} else {(u + T::one(), v + T::one())};
      let is_begin = if is_forward {pos_pair_2 == leftmost_pos_pair} else {pos_pair_2 == rightmost_pos_pair};
      let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
      match tmp_part_func_set_mat.get(&pos_pair_2) {
        Some(part_func_sets) => {
          let ref part_funcs = part_func_sets.part_funcs_4_ml;
          let score = part_funcs.part_func_4_align + ba_score + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
          logsumexp(&mut sum_4_ml, score);
          let score = part_funcs.part_func_4_insert + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          logsumexp(&mut sum_4_ml, score);
          let score = part_funcs.part_func_4_insert_2 + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          logsumexp(&mut sum_4_ml, score);
          let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
          let score = part_funcs.part_func_4_align + ba_score + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let score = part_funcs.part_func_4_insert + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let score = part_funcs.part_func_4_insert_2 + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let ref part_funcs = part_func_sets.part_funcs_on_sa;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
          logsumexp(&mut sum_on_sa, score);
          let score = part_funcs.part_func_4_insert + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          logsumexp(&mut sum_on_sa, score);
          let score = part_funcs.part_func_4_insert_2 + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          logsumexp(&mut sum_on_sa, score);
        }, None => {},
      }
      tmp_part_func_sets.part_funcs_4_ml.part_func_4_align = sum_4_ml;
      logsumexp(&mut tmp_sum, sum_4_ml);
      tmp_part_func_sets.part_funcs_4_first_bpas_on_mls.part_func_4_align = sum_4_first_bpas_on_mls;
      logsumexp(&mut tmp_sum, sum_4_first_bpas_on_mls);
      tmp_part_func_sets.part_funcs_4_bpas_on_mls.part_func_4_align = tmp_sum;
      let sum_on_sa = sum_on_sa + ba_score;
      tmp_part_func_sets.part_funcs_on_sa.part_func_4_align = sum_on_sa;
      logsumexp(&mut tmp_sum, sum_on_sa);
      tmp_part_func_sets.part_funcs_on_mls.part_func_4_align = tmp_sum;
      // For inserts.
      let mut sum_on_sa = NEG_INFINITY;
      let mut sum_4_ml = sum_on_sa;
      let mut sum_4_first_bpas_on_mls = sum_on_sa;
      let mut tmp_sum = sum_on_sa;
      let pos_pair_2 = if is_forward {(u - T::one(), v)} else {(u + T::one(), v)};
      let is_begin = pos_pair_2 == if is_forward {leftmost_pos_pair} else {rightmost_pos_pair};
      match tmp_part_func_set_mat.get(&pos_pair_2) {
        Some(part_func_sets) => {
          let ref part_funcs = part_func_sets.part_funcs_4_ml;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          logsumexp(&mut sum_4_ml, score);
          let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
          logsumexp(&mut sum_4_ml, score);
          let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
          logsumexp(&mut sum_4_ml, score);
          let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let ref part_funcs = part_func_sets.part_funcs_on_sa;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          logsumexp(&mut sum_on_sa, score);
          let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
          logsumexp(&mut sum_on_sa, score);
          let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
          logsumexp(&mut sum_on_sa, score);
        }, None => {},
      }
      let sum_4_ml = sum_4_ml + insert_score;
      tmp_part_func_sets.part_funcs_4_ml.part_func_4_insert = sum_4_ml;
      logsumexp(&mut tmp_sum, sum_4_ml);
      let sum_4_first_bpas_on_mls = sum_4_first_bpas_on_mls + insert_score;
      tmp_part_func_sets.part_funcs_4_first_bpas_on_mls.part_func_4_insert = sum_4_first_bpas_on_mls;
      logsumexp(&mut tmp_sum, sum_4_first_bpas_on_mls);
      tmp_part_func_sets.part_funcs_4_bpas_on_mls.part_func_4_insert = tmp_sum;
      let sum_on_sa = sum_on_sa + insert_score;
      tmp_part_func_sets.part_funcs_on_sa.part_func_4_insert = sum_on_sa;
      logsumexp(&mut tmp_sum, sum_on_sa);
      tmp_part_func_sets.part_funcs_on_mls.part_func_4_insert = tmp_sum;
      // For inserts on the other side.
      let mut sum_on_sa = NEG_INFINITY;
      let mut sum_4_ml = sum_on_sa;
      let mut sum_4_first_bpas_on_mls = sum_on_sa;
      let mut tmp_sum = sum_on_sa;
      let pos_pair_2 = if is_forward {(u, v - T::one())} else {(u, v + T::one())};
      let is_begin = pos_pair_2 == if is_forward {leftmost_pos_pair} else {rightmost_pos_pair};
      let insert_score = sta_fe_params.insert_scores_2[long_v];
      match tmp_part_func_set_mat.get(&pos_pair_2) {
        Some(part_func_sets) => {
          let ref part_funcs = part_func_sets.part_funcs_4_ml;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          logsumexp(&mut sum_4_ml, score);
          let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
          logsumexp(&mut sum_4_ml, score);
          let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
          logsumexp(&mut sum_4_ml, score);
          let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
          logsumexp(&mut sum_4_first_bpas_on_mls, score);
          let ref part_funcs = part_func_sets.part_funcs_on_sa;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          logsumexp(&mut sum_on_sa, score);
          let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
          logsumexp(&mut sum_on_sa, score);
          let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
          logsumexp(&mut sum_on_sa, score);
        }, None => {},
      }
      let sum_4_ml = sum_4_ml + insert_score;
      tmp_part_func_sets.part_funcs_4_ml.part_func_4_insert_2 = sum_4_ml;
      logsumexp(&mut tmp_sum, sum_4_ml);
      let sum_4_first_bpas_on_ml = sum_4_first_bpas_on_mls + insert_score;
      tmp_part_func_sets.part_funcs_4_first_bpas_on_mls.part_func_4_insert_2 = sum_4_first_bpas_on_ml;
      logsumexp(&mut tmp_sum, sum_4_first_bpas_on_mls);
      tmp_part_func_sets.part_funcs_4_bpas_on_mls.part_func_4_insert_2 = tmp_sum;
      let sum_on_sa = sum_on_sa + insert_score;
      tmp_part_func_sets.part_funcs_on_sa.part_func_4_insert_2 = sum_on_sa;
      logsumexp(&mut tmp_sum, sum_on_sa);
      tmp_part_func_sets.part_funcs_on_mls.part_func_4_insert_2 = tmp_sum;
      tmp_part_func_set_mat.insert(pos_pair, tmp_part_func_sets);
    }
  }
  let mut final_sum_on_sa = NEG_INFINITY;
  let mut final_sum_4_ml = final_sum_on_sa;
  if is_forward {
    let ref part_func_sets = tmp_part_func_set_mat[&rightmost_pos_pair];
    let ref part_funcs = part_func_sets.part_funcs_on_sa;
    logsumexp(&mut final_sum_on_sa, part_funcs.part_func_4_align);
    logsumexp(&mut final_sum_on_sa, part_funcs.part_func_4_insert);
    logsumexp(&mut final_sum_on_sa, part_funcs.part_func_4_insert_2);
    let ref part_funcs = part_func_sets.part_funcs_4_ml;
    logsumexp(&mut final_sum_4_ml, part_funcs.part_func_4_align);
    logsumexp(&mut final_sum_4_ml, part_funcs.part_func_4_insert);
    logsumexp(&mut final_sum_4_ml, part_funcs.part_func_4_insert_2);
  }
  (tmp_part_func_set_mat, final_sum_on_sa, final_sum_4_ml)
}

pub fn get_sta_prob_mats<T>(seq_pair: &SeqPair, seq_len_pair: &PosPair<T>, sta_fe_params: &StaFeParams<T>, max_bp_span_pair: &PosPair<T>, max_gap_num: T, max_gap_num_4_il: T, sta_part_func_mats: &StaPartFuncMats<T>, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair<T>, produces_access_probs: bool, global_part_func: PartFunc, pos_quadruple_mat: &PosQuadrupleMat<T>,) -> StaProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let pseudo_pos_quadruple = (T::zero(), seq_len_pair.0 - T::one(), T::zero(), seq_len_pair.1 - T::one());
  let mut sta_outside_part_func_4d_mat_4_bpas = PartFunc4dMat::<T>::default();
  let mut sta_prob_mats = StaProbMats::<T>::new(&seq_len_pair);
  for substr_len_1 in range_inclusive(T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(), max_bp_span_pair.0).rev() {
    for substr_len_2 in range_inclusive(T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(), max_bp_span_pair.1).rev() {
      for &(i, j, k, l) in pos_quadruple_mat {
        if j - i + T::one() != substr_len_1 || l - k + T::one() != substr_len_2 {continue;}
        let pos_quadruple = (i, j, k, l);
        match sta_part_func_mats.part_func_4d_mat_4_bpas.get(&pos_quadruple) {
          Some(&part_func_4_bpa) => {
            let (long_i, long_j, long_k, long_l) = (i.to_usize().unwrap(), j.to_usize().unwrap(), k.to_usize().unwrap(), l.to_usize().unwrap());
            let pos_pair = (i, j);
            let pos_pair_2 = (k, l);
            let prob_coeff = part_func_4_bpa - global_part_func;
            let mut sum = NEG_INFINITY;
            let mut forward_term = sum;
            let mut backward_term = sum;
            let is_begin = i == T::one() && k == T::one();
            let is_end = j == seq_len_pair.0 - T::from_usize(2).unwrap() && l == seq_len_pair.1 - T::from_usize(2).unwrap();
            match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&(i - T::one(), k - T::one())) {
              Some(part_funcs) => {
                let term = part_funcs.part_func_4_align + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
                logsumexp(&mut forward_term, term);
                let term = part_funcs.part_func_4_insert + if is_begin {INIT_MATCH_SCORE} else {INSERT_2_MATCH_SCORE};
                logsumexp(&mut forward_term, term);
                let term = part_funcs.part_func_4_insert_2 + if is_begin {INIT_MATCH_SCORE} else {INSERT_2_MATCH_SCORE};
                logsumexp(&mut forward_term, term);
              }, None => {},
            }
            match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&(j + T::one(), l + T::one())) {
              Some(part_funcs) => {
                let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                logsumexp(&mut backward_term, term);
                let term = part_funcs.part_func_4_insert + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                logsumexp(&mut backward_term, term);
                let term = part_funcs.part_func_4_insert_2 + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                logsumexp(&mut backward_term, term);
              }, None => {},
            }
            let part_func_4_el = forward_term + backward_term;
            if part_func_4_el > NEG_INFINITY {
              let coefficient = sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els[&pos_quadruple] - part_func_4_bpa;
              sum = coefficient + part_func_4_el;
              if produces_access_probs {
                let bpap_4_el = prob_coeff + sum;
                match sta_prob_mats.bpp_mat_pair_4_el.0.get_mut(&pos_pair) {
                  Some(bpp_4_el) => {
                    logsumexp(bpp_4_el, bpap_4_el);
                  }, None => {
                    sta_prob_mats.bpp_mat_pair_4_el.0.insert(pos_pair, bpap_4_el);
                  },
                }
                match sta_prob_mats.bpp_mat_pair_4_el.1.get_mut(&pos_pair_2) {
                  Some(bpp_4_el) => {
                    logsumexp(bpp_4_el, bpap_4_el);
                  }, None => {
                    sta_prob_mats.bpp_mat_pair_4_el.1.insert(pos_pair_2, bpap_4_el);
                  },
                }
              }
            }
            for &(m, n, o, p) in pos_quadruple_mat {
              if !(m < i && j < n) || !(o < k && l < p) {continue;}
              let (long_m, long_n, long_o, long_p ) = (m.to_usize().unwrap(), n.to_usize().unwrap(), o.to_usize().unwrap(), p.to_usize().unwrap());
              if long_n - long_j - 1 + long_i - long_m - 1 > MAX_2_LOOP_LEN {continue;}
              if long_p - long_l - 1 + long_k - long_o - 1 > MAX_2_LOOP_LEN {continue;}
              let pos_quadruple_2 = (m, n, o, p);
              match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                Some(&part_func) => {
                  let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs[&(m, o)];
                  let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs[&(n, p)];
                  let mut forward_term = NEG_INFINITY;
                  let mut backward_term = forward_term;
                  match forward_tmp_part_func_set_mat.get(&(i - T::one(), k - T::one())) {
                    Some(part_func_sets) => {
                      let ref part_funcs = part_func_sets.part_funcs_on_sa;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                      let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                      let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                    }, None => {},
                  }
                  match backward_tmp_part_func_set_mat.get(&(j + T::one(), l + T::one())) {
                    Some(part_func_sets) => {
                      let ref part_funcs = part_func_sets.part_funcs_on_sa;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut backward_term, term);
                      let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term, term);
                      let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term, term);
                    }, None => {},
                  }
                  let part_func_4_2l = forward_term + backward_term;
                  if part_func_4_2l > NEG_INFINITY {
                    let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                    let twoloop_fe = ss_free_energy_mat_set_pair.0.twoloop_fe_4d_mat[&(m, n, i, j)];
                    let twoloop_fe_2 = ss_free_energy_mat_set_pair.1.twoloop_fe_4d_mat[&(o, p, k, l)];
                    let coefficient = bpa_score + twoloop_fe + twoloop_fe_2 + part_func;
                    let part_func_4_2l = coefficient + part_func_4_2l;
                    logsumexp(&mut sum, part_func_4_2l);
                    if produces_access_probs {
                      let bpap_4_2l = prob_coeff + part_func_4_2l;
                      match sta_prob_mats.access_bpp_mat_pair_4_2l.0.get_mut(&pos_pair) {
                        Some(bpp_4_2l) => {
                          logsumexp(bpp_4_2l, bpap_4_2l);
                        }, None => {
                          sta_prob_mats.access_bpp_mat_pair_4_2l.0.insert(pos_pair, bpap_4_2l);
                        },
                      }
                      match sta_prob_mats.access_bpp_mat_pair_4_2l.1.get_mut(&pos_pair_2) {
                        Some(bpp_4_2l) => {
                          logsumexp(bpp_4_2l, bpap_4_2l);
                        }, None => {
                          sta_prob_mats.access_bpp_mat_pair_4_2l.1.insert(pos_pair_2, bpap_4_2l);
                        },
                      }
                      for q in long_m + 1 .. long_i {
                        logsumexp(&mut sta_prob_mats.upp_mat_pair_4_2l.0[q], bpap_4_2l);
                      }
                      for q in long_j + 1 .. long_n {
                        logsumexp(&mut sta_prob_mats.upp_mat_pair_4_2l.0[q], bpap_4_2l);
                      }
                      for q in long_o + 1 .. long_k {
                        logsumexp(&mut sta_prob_mats.upp_mat_pair_4_2l.1[q], bpap_4_2l);
                      }
                      for q in long_l + 1 .. long_p {
                        logsumexp(&mut sta_prob_mats.upp_mat_pair_4_2l.1[q], bpap_4_2l);
                      }
                    }
                  }
                }, None => {},
              }
            }
            let part_func_ratio = sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple] - part_func_4_bpa;
            for &(m, n, o, p) in pos_quadruple_mat {
              if !(m < i && j < n) || !(o < k && l < p) {continue;}
              let pos_quadruple_2 = (m, n, o, p);
              match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                Some(&part_func_4_bpa_2) => {
                  let (long_m, long_n, long_o, long_p) = (m.to_usize().unwrap(), n.to_usize().unwrap(), o.to_usize().unwrap(), p.to_usize().unwrap());
                  let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs[&(m, o)];
                  let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs[&(n, p)];
                  let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                  let mut forward_term = NEG_INFINITY;
                  let mut forward_term_2 = forward_term;
                  let mut backward_term = forward_term;
                  let mut backward_term_2 = forward_term;
                  match forward_tmp_part_func_set_mat.get(&(i - T::one(), k - T::one())) {
                    Some(part_func_sets) => {
                      let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                      let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                      let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term, term);
                      let ref part_funcs = part_func_sets.part_funcs_on_sa;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut forward_term_2, term);
                      let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term_2, term);
                      let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                      logsumexp(&mut forward_term_2, term);
                    }, None => {},
                  }
                  match backward_tmp_part_func_set_mat.get(&(j + T::one(), l + T::one())) {
                    Some(part_func_sets) => {
                      let ref part_funcs = part_func_sets.part_funcs_on_mls;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut backward_term, term);
                      let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term, term);
                      let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term, term);
                      let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                      let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                      logsumexp(&mut backward_term_2, term);
                      let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term_2, term);
                      let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                      logsumexp(&mut backward_term_2, term);
                    }, None => {},
                  }
                  let mut part_func_4_ml = forward_term + backward_term;
                  logsumexp(&mut part_func_4_ml, forward_term_2 + backward_term_2);
                  if part_func_4_ml > NEG_INFINITY {
                    let multi_loop_closing_basepairing_fe = get_ml_closing_basepairing_fe(
                      seq_pair.0,
                      &(long_m, long_n),
                    );
                    let multi_loop_closing_basepairing_fe_2 =
                      get_ml_closing_basepairing_fe(
                        seq_pair.1,
                        &(long_o, long_p),
                      );
                    let coefficient = part_func_ratio + bpa_score + multi_loop_closing_basepairing_fe + multi_loop_closing_basepairing_fe_2 + part_func_4_bpa_2;
                    let part_func_4_ml = coefficient + part_func_4_ml;
                    logsumexp(&mut sum, part_func_4_ml);
                    if produces_access_probs {
                      let bpap_4_ml = prob_coeff + part_func_4_ml;
                      match sta_prob_mats.access_bpp_mat_pair_4_ml.0.get_mut(&pos_pair) {
                        Some(bpp_4_ml) => {
                          logsumexp(bpp_4_ml, bpap_4_ml);
                        }, None => {
                          sta_prob_mats.access_bpp_mat_pair_4_ml.0.insert(pos_pair, bpap_4_ml);
                        },
                      }
                      match sta_prob_mats.access_bpp_mat_pair_4_ml.1.get_mut(&pos_pair_2) {
                        Some(bpp_4_ml) => {
                          logsumexp(bpp_4_ml, bpap_4_ml);
                        }, None => {
                          sta_prob_mats.access_bpp_mat_pair_4_ml.1.insert(pos_pair_2, bpap_4_ml);
                        },
                      }
                    }
                  }
                }, None => {},
              }
            }
            if sum > NEG_INFINITY {
              sta_outside_part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
              let bpap = expf(prob_coeff + sum);
              debug_assert!(0. <= bpap && bpap <= 1.);
              match sta_prob_mats.bpp_mat_pair.0.get_mut(&(i, j)) {
                Some(bpp) => {
                  *bpp += bpap;
                }, None => {
                  sta_prob_mats.bpp_mat_pair.0.insert((i, j), bpap);
                },
              }
              match sta_prob_mats.bpp_mat_pair.1.get_mut(&(k, l)) {
                Some(bpp) => {
                  *bpp += bpap;
                }, None => {
                  sta_prob_mats.bpp_mat_pair.1.insert((k, l), bpap);
                },
              }
              sta_prob_mats.upp_mat_pair.0[long_i] -= bpap;
              sta_prob_mats.upp_mat_pair.0[long_j] -= bpap;
              sta_prob_mats.upp_mat_pair.1[long_k] -= bpap;
              sta_prob_mats.upp_mat_pair.1[long_l] -= bpap;
            }
          }, None => {},
        }
      }
    }
  }
  if produces_access_probs {
    let leftmost_pos_pair = (T::zero(), T::zero());
    let rightmost_pos_pair = (seq_len_pair.0 - T::one(), seq_len_pair.1 - T::one());
    for u in range(T::zero(), seq_len_pair.0 - T::one()) {
      let long_u = u.to_usize().unwrap();
      let insert_score = sta_fe_params.insert_scores[long_u];
      for v in range(T::zero(), seq_len_pair.1 - T::one()) {
        if u == T::zero() && v == T::zero() {continue;}
        let pos_pair = (u, v);
        if !is_min_gap_ok(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
        let long_v = v.to_usize().unwrap();
        let insert_score_2 = sta_fe_params.insert_scores_2[long_v];
        let pos_pair_4_ba = (u - T::one(), v - T::one());
        let pos_pair_4_gap_1 = (u - T::one(), v);
        let pos_pair_4_gap_2 = (u, v - T::one());
        let pos_pair_2 = (u + T::one(), v + T::one());
        let is_end = pos_pair_2 == rightmost_pos_pair;
        let mut backward_term_4_align = NEG_INFINITY;
        let mut backward_term_4_insert = backward_term_4_align;
        let mut backward_term_4_insert_2 = backward_term_4_align;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            if u > T::zero() && v > T::zero() {
              let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
              logsumexp(&mut backward_term_4_align, term);
              let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
              logsumexp(&mut backward_term_4_align, term);
              let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
              logsumexp(&mut backward_term_4_align, term);
            }
            if u > T::zero() {
              let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
              logsumexp(&mut backward_term_4_insert, term);
              let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
              logsumexp(&mut backward_term_4_insert, term);
              let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
              logsumexp(&mut backward_term_4_insert, term);
            }
            if v > T::zero() {
              let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
              logsumexp(&mut backward_term_4_insert_2, term);
              let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
              logsumexp(&mut backward_term_4_insert_2, term);
              let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
              logsumexp(&mut backward_term_4_insert_2, term);
            }
          }, None => {},
        }
        if u > T::zero() && v > T::zero() {
          match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_4_ba) {
            Some(part_funcs) => {
              let mut forward_term = NEG_INFINITY;
              let is_begin = pos_pair_4_ba == leftmost_pos_pair;
              let term = part_funcs.part_func_4_align + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
              logsumexp(&mut forward_term, term);
              let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
              logsumexp(&mut forward_term, term);
              let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
              logsumexp(&mut forward_term, term);
              let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
              let bap_4_el = ba_score + forward_term + backward_term_4_align - global_part_func;
              logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.0[long_u], bap_4_el);
              logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.1[long_v], bap_4_el);
            }, None => {},
          }
        }
        if u > T::zero() {
          match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_4_gap_1) {
            Some(part_funcs) => {
              let mut forward_term = NEG_INFINITY;
              let is_begin = pos_pair_4_gap_1 == leftmost_pos_pair;
              let term = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
              logsumexp(&mut forward_term, term);
              let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
              logsumexp(&mut forward_term, term);
              let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
              logsumexp(&mut forward_term, term);
              let upp_4_el = insert_score + forward_term + backward_term_4_insert - global_part_func;
              logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.0[long_u], upp_4_el);
            }, None => {},
          }
        }
        if v > T::zero() {
          match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_4_gap_2) {
            Some(part_funcs) => {
              let mut forward_term = NEG_INFINITY;
              let is_begin = pos_pair_4_gap_2 == leftmost_pos_pair;
              let term = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
              logsumexp(&mut forward_term, term);
              let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
              logsumexp(&mut forward_term, term);
              let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
              logsumexp(&mut forward_term, term);
              let upp_4_el = insert_score_2 + forward_term + backward_term_4_insert_2 - global_part_func;
              logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.1[long_v], upp_4_el);
            }, None => {},
          }
        }
        if !is_min_gap_ok(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
        if !(u > T::zero() && v > T::zero()) {continue;}
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        for &(i, j, k, l) in pos_quadruple_mat {
          if !(i < u && u < j) || !(k < v && v < l) {
            continue;
          }
          let pos_quadruple = (i, j, k, l);
          match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple) {
            Some(&part_func_4_bpa) => {
              let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
              let prob_coeff = part_func_4_bpa - global_part_func + bpa_score;
              let hl_fe = ss_free_energy_mat_set_pair.0.hl_fe_mat[&(i, j)];
              let hl_fe_2 = ss_free_energy_mat_set_pair.1.hl_fe_mat[&(k, l)];
              let (long_i, long_j, long_k, long_l) = (i.to_usize().unwrap(), j.to_usize().unwrap(), k.to_usize().unwrap(), l.to_usize().unwrap());
              let multi_loop_closing_basepairing_fe = get_ml_closing_basepairing_fe(
                seq_pair.0,
                &(long_i, long_j),
              );
              let multi_loop_closing_basepairing_fe_2 =
                get_ml_closing_basepairing_fe(
                  seq_pair.1,
                  &(long_k, long_l),
                );
              let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs[&(i, k)];
              let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs[&(j, l)];
              let rightmost_pos_pair = (j, l);
              let is_end = pos_pair_2 == rightmost_pos_pair;
              let mut backward_term_4_align_on_sa = NEG_INFINITY;
              let mut backward_term_4_insert_on_sa = backward_term_4_align_on_sa;
              let mut backward_term_4_insert_2_on_sa = backward_term_4_align_on_sa;
              let mut backward_term_4_align_4_ml = backward_term_4_align_on_sa;
              let mut backward_term_4_insert_4_ml = backward_term_4_align_on_sa;
              let mut backward_term_4_insert_2_4_ml = backward_term_4_align_on_sa;
              let mut backward_term_4_align_4_bpas_on_mls = backward_term_4_align_on_sa;
              let mut backward_term_4_insert_4_bpas_on_mls = backward_term_4_align_on_sa;
              let mut backward_term_4_insert_2_4_bpas_on_mls = backward_term_4_align_on_sa;
              let mut backward_term_4_align_on_mls = backward_term_4_align_on_sa;
              let mut backward_term_4_insert_on_mls = backward_term_4_align_on_sa;
              let mut backward_term_4_insert_2_on_mls = backward_term_4_align_on_sa;
              match backward_tmp_part_func_set_mat.get(&pos_pair_2) {
                Some(part_func_sets) => {
                  let ref part_funcs = part_func_sets.part_funcs_on_sa;
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_align_on_sa, term);
                  let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_on_sa, term);
                  let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_on_sa, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_on_sa, term);
                  let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_on_sa, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_on_sa, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_2_on_sa, term);
                  let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_on_sa, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_on_sa, term);
                  let ref part_funcs = part_func_sets.part_funcs_4_ml;
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_align_4_ml, term);
                  let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_4_ml, term);
                  let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_4_ml, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_4_ml, term);
                  let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_4_ml, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_4_ml, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_2_4_ml, term);
                  let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_4_ml, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_4_ml, term);
                  let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_align_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_2_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_4_bpas_on_mls, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_4_bpas_on_mls, term);
                  let ref part_funcs = part_func_sets.part_funcs_on_mls;
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_align_on_mls, term);
                  let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_on_mls, term);
                  let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut backward_term_4_align_on_mls, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_on_mls, term);
                  let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_on_mls, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_on_mls, term);
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                  logsumexp(&mut backward_term_4_insert_2_on_mls, term);
                  let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_on_mls, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                  logsumexp(&mut backward_term_4_insert_2_on_mls, term);
                }, None => {},
              }
              let prob_coeff_4_hl = prob_coeff + hl_fe + hl_fe_2;
              let prob_coeff_4_ml = prob_coeff + multi_loop_closing_basepairing_fe + multi_loop_closing_basepairing_fe_2;
              match forward_tmp_part_func_set_mat.get(&pos_pair_4_ba) {
                Some(part_func_sets) => {
                  let ref part_funcs = part_func_sets.part_funcs_on_sa;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let bap_4_hl = prob_coeff_4_hl + ba_score + forward_term + backward_term_4_align_on_sa;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.0[long_u], bap_4_hl);
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.1[long_v], bap_4_hl);
                  let bap_4_ml = prob_coeff_4_ml + ba_score + forward_term + backward_term_4_align_4_ml;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], bap_4_ml);
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], bap_4_ml);
                  let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let bap_4_ml = prob_coeff_4_ml + ba_score + forward_term + backward_term_4_align_4_bpas_on_mls;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], bap_4_ml);
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], bap_4_ml);
                  let ref part_funcs = part_func_sets.part_funcs_4_ml;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let bap_4_ml = prob_coeff_4_ml + ba_score + forward_term + backward_term_4_align_on_mls;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], bap_4_ml);
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], bap_4_ml);
                }, None => {},
              }
              match forward_tmp_part_func_set_mat.get(&pos_pair_4_gap_1) {
                Some(part_func_sets) => {
                  let ref part_funcs = part_func_sets.part_funcs_on_sa;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let upp_4_hl = prob_coeff_4_hl + insert_score + forward_term + backward_term_4_insert_on_sa;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.0[long_u], upp_4_hl);
                  let upp_4_ml = prob_coeff_4_ml + insert_score + forward_term + backward_term_4_insert_4_ml;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], upp_4_ml);
                  let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let upp_4_ml = prob_coeff_4_ml + insert_score + forward_term + backward_term_4_insert_4_bpas_on_mls;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], upp_4_ml);
                  let ref part_funcs = part_func_sets.part_funcs_4_ml;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let upp_4_ml = prob_coeff_4_ml + insert_score + forward_term + backward_term_4_insert_on_mls;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], upp_4_ml);
                }, None => {},
              }
              match forward_tmp_part_func_set_mat.get(&pos_pair_4_gap_2) {
                Some(part_func_sets) => {
                  let ref part_funcs = part_func_sets.part_funcs_on_sa;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                  logsumexp(&mut forward_term, term);
                  let upp_4_hl = prob_coeff_4_hl + insert_score_2 + forward_term + backward_term_4_insert_2_on_sa;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.1[long_v], upp_4_hl);
                  let upp_4_ml = prob_coeff_4_ml + insert_score_2 + forward_term + backward_term_4_insert_2_4_ml;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], upp_4_ml);
                  let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                  logsumexp(&mut forward_term, term);
                  let upp_4_ml = prob_coeff_4_ml + insert_score_2 + forward_term + backward_term_4_insert_2_4_bpas_on_mls;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], upp_4_ml);
                  let ref part_funcs = part_func_sets.part_funcs_4_ml;
                  let mut forward_term = NEG_INFINITY;
                  let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                  logsumexp(&mut forward_term, term);
                  let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                  logsumexp(&mut forward_term, term);
                  let upp_4_ml = prob_coeff_4_ml + insert_score_2 + forward_term + backward_term_4_insert_2_on_mls;
                  logsumexp(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], upp_4_ml);
                }, None => {},
              }
            }, None => {},
          }
        }
      }
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_2l.0.values_mut() {
      *bpp = expf(*bpp)
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_2l.1.values_mut() {
      *bpp = expf(*bpp)
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_ml.0.values_mut() {
      *bpp = expf(*bpp)
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_ml.1.values_mut() {
      *bpp = expf(*bpp)
    }
    for bpp in sta_prob_mats.bpp_mat_pair_4_el.0.values_mut() {
      *bpp = expf(*bpp)
    }
    for bpp in sta_prob_mats.bpp_mat_pair_4_el.1.values_mut() {
      *bpp = expf(*bpp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_hl.0.iter_mut() {
      *upp = expf(*upp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_hl.1.iter_mut() {
      *upp = expf(*upp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_ml.0.iter_mut() {
      *upp = expf(*upp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_ml.1.iter_mut() {
      *upp = expf(*upp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_2l.0.iter_mut() {
      *upp = expf(*upp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_2l.1.iter_mut() {
      *upp = expf(*upp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_el.0.iter_mut() {
      *upp = expf(*upp)
    }
    for upp in sta_prob_mats.upp_mat_pair_4_el.1.iter_mut() {
      *upp = expf(*upp)
    }
  }
  sta_prob_mats
}

pub fn pct_of_prob_mats<T>(prob_mats_with_rna_id_pairs: &StaProbMatsWithRnaIdPairs<T>, rna_id: RnaId, num_of_rnas: usize, bpp_mat: &SparseProbMat<T>, upp_mat_len: usize, produces_access_probs: bool) -> PctStaProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let weight = 1. / (num_of_rnas - 1) as Prob;
  let mut pct_prob_mats = PctStaProbMats::new(upp_mat_len);
  pct_prob_mats.bpp_mat_on_ss = bpp_mat.clone();
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id == rna_id_2 {continue;}
    let rna_id_pair = if rna_id < rna_id_2 {(rna_id, rna_id_2)} else {(rna_id_2, rna_id)};
    let ref ref_2_prob_mats = prob_mats_with_rna_id_pairs[&rna_id_pair];
    let ref_2_bpp_mat = if rna_id < rna_id_2 {
      &ref_2_prob_mats.bpp_mat_pair.0
    } else {
      &ref_2_prob_mats.bpp_mat_pair.1
    };
    for (pos_pair, &bpp) in ref_2_bpp_mat.iter() {
      let weighted_bpp = weight * bpp;
      match pct_prob_mats.bpp_mat.get_mut(pos_pair) {
        Some(bpp) => {
          *bpp += weighted_bpp;
        },
        None => {
          pct_prob_mats.bpp_mat.insert(*pos_pair, weighted_bpp);
        },
      }
    }
    let ref_2_upp_mat = if rna_id < rna_id_2 {
      &ref_2_prob_mats.upp_mat_pair.0
    } else {
      &ref_2_prob_mats.upp_mat_pair.1
    };
    for (i, &upp) in ref_2_upp_mat.iter().enumerate() {
      let weighted_upp = weight * upp;
      pct_prob_mats.upp_mat[i] += weighted_upp;
    }
    if produces_access_probs {
      let ref_2_access_bpp_mat_4_2l = if rna_id < rna_id_2 {
        &ref_2_prob_mats.access_bpp_mat_pair_4_2l.0
      } else {
        &ref_2_prob_mats.access_bpp_mat_pair_4_2l.1
      };
      for (pos_pair, &bpp) in ref_2_access_bpp_mat_4_2l.iter() {
        let weighted_bpp = weight * bpp;
        match pct_prob_mats.access_bpp_mat_4_2l.get_mut(pos_pair) {
          Some(bpp_4_2l) => {
            *bpp_4_2l += weighted_bpp;
          },
          None => {
            pct_prob_mats.access_bpp_mat_4_2l.insert(*pos_pair, weighted_bpp);
          },
        }
      }
      let ref_2_access_bpp_mat_4_ml = if rna_id < rna_id_2 {
        &ref_2_prob_mats.access_bpp_mat_pair_4_ml.0
      } else {
        &ref_2_prob_mats.access_bpp_mat_pair_4_ml.1
      };
      for (pos_pair, &bpp) in ref_2_access_bpp_mat_4_ml.iter() {
        let weighted_bpp = weight * bpp;
        match pct_prob_mats.access_bpp_mat_4_ml.get_mut(pos_pair) {
          Some(bpp_4_ml) => {
            *bpp_4_ml += weighted_bpp;
          },
          None => {
            pct_prob_mats.access_bpp_mat_4_ml.insert(*pos_pair, weighted_bpp);
          },
        }
      }
      let ref_2_bpp_mat_4_el = if rna_id < rna_id_2 {
        &ref_2_prob_mats.bpp_mat_pair_4_el.0
      } else {
        &ref_2_prob_mats.bpp_mat_pair_4_el.1
      };
      for (pos_pair, &bpp) in ref_2_bpp_mat_4_el.iter() {
        let weighted_bpp = weight * bpp;
        match pct_prob_mats.bpp_mat_4_el.get_mut(pos_pair) {
          Some(bpp_4_el) => {
            *bpp_4_el += weighted_bpp;
          },
          None => {
            pct_prob_mats.bpp_mat_4_el.insert(*pos_pair, weighted_bpp);
          },
        }
      }
      let ref_2_upp_mat = if rna_id < rna_id_2 {
        &ref_2_prob_mats.upp_mat_pair_4_hl.0
      } else {
        &ref_2_prob_mats.upp_mat_pair_4_hl.1
      };
      for (i, &upp) in ref_2_upp_mat.iter().enumerate() {
        let weighted_upp = weight * upp;
        pct_prob_mats.upp_mat_4_hl[i] += weighted_upp;
      }
      let ref_2_upp_mat = if rna_id < rna_id_2 {
        &ref_2_prob_mats.upp_mat_pair_4_2l.0
      } else {
        &ref_2_prob_mats.upp_mat_pair_4_2l.1
      };
      for (i, &upp) in ref_2_upp_mat.iter().enumerate() {
        let weighted_upp = weight * upp;
        pct_prob_mats.upp_mat_4_2l[i] += weighted_upp;
      }
      let ref_2_upp_mat = if rna_id < rna_id_2 {
        &ref_2_prob_mats.upp_mat_pair_4_ml.0
      } else {
        &ref_2_prob_mats.upp_mat_pair_4_ml.1
      };
      for (i, &upp) in ref_2_upp_mat.iter().enumerate() {
        let weighted_upp = weight * upp;
        pct_prob_mats.upp_mat_4_ml[i] += weighted_upp;
      }
      let ref_2_upp_mat = if rna_id < rna_id_2 {
        &ref_2_prob_mats.upp_mat_pair_4_el.0
      } else {
        &ref_2_prob_mats.upp_mat_pair_4_el.1
      };
      for (i, &upp) in ref_2_upp_mat.iter().enumerate() {
        let weighted_upp = weight * upp;
        pct_prob_mats.upp_mat_4_el[i] += weighted_upp;
      }
    }
  }
  for (&(i, j), &bpp) in bpp_mat {
    let pos_pair = (i + T::one(), j + T::one());
    if !pct_prob_mats.bpp_mat.contains_key(&pos_pair) {
      pct_prob_mats.bpp_mat.insert(pos_pair, bpp);
    }
  }
  if produces_access_probs {
    for (i, &upp) in pct_prob_mats.upp_mat_4_hl.iter().enumerate() {
      pct_prob_mats.max_upp_mat[i] = upp;
    }
    for (i, &upp) in pct_prob_mats.upp_mat_4_2l.iter().enumerate() {
      let old_upp = pct_prob_mats.max_upp_mat[i];
      if upp > old_upp {
        pct_prob_mats.max_upp_mat[i] = upp;
      }
    }
    for (i, &upp) in pct_prob_mats.upp_mat_4_ml.iter().enumerate() {
      let old_upp = pct_prob_mats.max_upp_mat[i];
      if upp > old_upp {
        pct_prob_mats.max_upp_mat[i] = upp;
      }
    }
    for (i, &upp) in pct_prob_mats.upp_mat_4_el.iter().enumerate() {
      let old_upp = pct_prob_mats.max_upp_mat[i];
      if upp > old_upp {
        pct_prob_mats.max_upp_mat[i] = upp;
      }
    }
    for (pos_pair, &bpp) in pct_prob_mats.access_bpp_mat_4_2l.iter() {
      pct_prob_mats.max_bpp_mat.insert(*pos_pair, bpp);
    }
    for (pos_pair, &bpp) in pct_prob_mats.access_bpp_mat_4_ml.iter() {
      match pct_prob_mats.max_bpp_mat.get_mut(pos_pair) {
        Some(old_bpp) => {
          if bpp > *old_bpp {
            *old_bpp = bpp;
          }
        }, None => {
          pct_prob_mats.max_bpp_mat.insert(*pos_pair, bpp);
        }
      }
    }
    for (pos_pair, &bpp) in pct_prob_mats.bpp_mat_4_el.iter() {
      match pct_prob_mats.max_bpp_mat.get_mut(pos_pair) {
        Some(old_bpp) => {
          if bpp > *old_bpp {
            *old_bpp = bpp;
          }
        }, None => {
          pct_prob_mats.max_bpp_mat.insert(*pos_pair, bpp);
        }
      }
    }
  }
  pct_prob_mats
}

pub fn get_max_bp_span<T>(sparse_bpp_mat: &SparseProbMat<T>) -> T
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  sparse_bpp_mat.iter().map(|(pos_pair, _)| {pos_pair.1 - pos_pair.0 + T::one()}).max().unwrap()
}

pub fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
}

pub fn get_seq_len_diff<T>(pos_quadruple: &PosQuadruple<T>) -> T
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let seq_len_pair = (pos_quadruple.1 + T::one() - pos_quadruple.0, pos_quadruple.3 + T::one() - pos_quadruple.2);
  max(seq_len_pair.0, seq_len_pair.1) - min(seq_len_pair.0, seq_len_pair.1)
}

pub fn remove_small_bpps_from_bpp_mat<T>(sparse_bpp_mat: &SparseProbMat<T>, min_bpp: Prob) -> SparseProbMat<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  sparse_bpp_mat.iter().filter(|(_, &bpp)| {bpp >= min_bpp}).map(|(&(i, j), &bpp)| {((i + T::one(), j + T::one()), bpp)}).collect()
}

fn is_min_gap_ok<T>(pos_pair: &PosPair<T>, pos_quadruple: &PosQuadruple<T>, max_gap_num: T) -> bool
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let min_gap_num_1 = get_seq_len_diff::<T>(&(pos_quadruple.0, pos_pair.0, pos_quadruple.2, pos_pair.1));
  let min_gap_num_2 = get_seq_len_diff::<T>(&(pos_pair.0, pos_quadruple.1, pos_pair.1, pos_quadruple.3));
  if min_gap_num_1 <= max_gap_num && min_gap_num_2 <= max_gap_num {
    true
  } else {
    false
  }
}

pub fn consprob<T>(thread_pool: &mut Pool, fasta_records: &FastaRecords, min_bpp: Prob, offset_4_max_gap_num: T, produces_access_probs: bool) -> ProbMatSets<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let num_of_fasta_records = fasta_records.len();
  let mut bpp_mats = vec![SparseProbMat::<T>::new(); num_of_fasta_records];
  let mut sparse_bpp_mats = vec![SparseProbMat::<T>::new(); num_of_fasta_records];
  let mut max_bp_spans = vec![T::zero(); num_of_fasta_records];
  let mut ss_free_energy_mat_sets = vec![SsFreeEnergyMats::<T>::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (bpp_mat, sparse_bpp_mat, max_bp_span, fasta_record, ss_free_energy_mats) in multizip((bpp_mats.iter_mut(), sparse_bpp_mats.iter_mut(), max_bp_spans.iter_mut(), fasta_records.iter(), ss_free_energy_mat_sets.iter_mut())) {
      let seq_len = fasta_record.seq.len();
      scope.execute(move || {
        let (obtained_bpp_mat, obtained_ss_free_energy_mats) = mccaskill_algo(&fasta_record.seq[1 .. seq_len - 1], false);
        *bpp_mat = obtained_bpp_mat;
        *ss_free_energy_mats = sparsify(&obtained_ss_free_energy_mats, bpp_mat, min_bpp);
        *sparse_bpp_mat = remove_small_bpps_from_bpp_mat::<T>(bpp_mat, min_bpp);
        *max_bp_span = get_max_bp_span::<T>(sparse_bpp_mat);
      });
    }
  });
  let mut prob_mats_with_rna_id_pairs = StaProbMatsWithRnaIdPairs::<T>::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      prob_mats_with_rna_id_pairs.insert(rna_id_pair, StaProbMats::<T>::origin());
    }
  }
  thread_pool.scoped(|scope| {
    for (rna_id_pair, prob_mats) in prob_mats_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
      let max_gap_num = offset_4_max_gap_num + T::from_usize(max(seq_len_pair.0, seq_len_pair.1) - min(seq_len_pair.0, seq_len_pair.1)).unwrap();
      let max_gap_num_4_il = max(min(max_gap_num, T::from_usize(MAX_GAP_NUM_4_IL).unwrap()), T::from_usize(MIN_GAP_NUM_4_IL).unwrap());
      let max_bp_span_pair = (max_bp_spans[rna_id_pair.0], max_bp_spans[rna_id_pair.1]);
      let bpp_mat_pair = (&sparse_bpp_mats[rna_id_pair.0], &sparse_bpp_mats[rna_id_pair.1]);
      let ss_free_energy_mat_set_pair = (&ss_free_energy_mat_sets[rna_id_pair.0], &ss_free_energy_mat_sets[rna_id_pair.1]);
      let seq_len_pair = (
        T::from_usize(seq_pair.0.len()).unwrap(),
        T::from_usize(seq_pair.1.len()).unwrap(),
      );
      let pseudo_pos_quadruple = (
        T::zero(),
        seq_len_pair.0 - T::one(),
        T::zero(),
        seq_len_pair.1 - T::one(),
      );
      let (forward_pos_pair_mat_set, backward_pos_pair_mat_set, pos_quadruple_mat) = get_sparse_pos_sets(&bpp_mat_pair, max_gap_num_4_il, &pseudo_pos_quadruple);
      let max = max(max_gap_num, max_gap_num_4_il);
      scope.execute(move || {
        let sta_fe_params = StaFeParams::<T>::new(&seq_pair, &seq_len_pair, max, &pos_quadruple_mat);
        *prob_mats = io_algo_4_prob_mats::<T>(&seq_pair, &seq_len_pair, &sta_fe_params, &max_bp_span_pair, max_gap_num, max_gap_num_4_il, &ss_free_energy_mat_set_pair, produces_access_probs, &forward_pos_pair_mat_set, &backward_pos_pair_mat_set, &pos_quadruple_mat,);
      });
    }
  });
  let mut prob_mat_sets = vec![PctStaProbMats::<T>::origin(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (rna_id, prob_mats, bpp_mat) in multizip((0 .. num_of_fasta_records, prob_mat_sets.iter_mut(), bpp_mats.iter_mut())) {
      let ref ref_2_prob_mats_with_rna_id_pairs = prob_mats_with_rna_id_pairs;
      let seq_len = fasta_records[rna_id].seq.len();
      scope.execute(move || {
        *prob_mats = pct_of_prob_mats::<T>(ref_2_prob_mats_with_rna_id_pairs, rna_id, num_of_fasta_records, bpp_mat, seq_len, produces_access_probs);
      });
    }
  });
  prob_mat_sets
}

pub fn sparsify<T>(ss_free_energy_mats: &SsFreeEnergyMats<T>, bpp_mat: &SparseProbMat<T>, min_bpp: Prob) -> SsFreeEnergyMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let mut new_ss_free_energy_mats = SsFreeEnergyMats::new();
  new_ss_free_energy_mats.hl_fe_mat = ss_free_energy_mats.hl_fe_mat.iter().filter(|(pos_pair, _)| {bpp_mat[pos_pair] >= min_bpp}).map(|(&(i, j), &free_energy)| {((i + T::one(), j + T::one()), free_energy)}).collect();
  new_ss_free_energy_mats.twoloop_fe_4d_mat = ss_free_energy_mats.twoloop_fe_4d_mat.iter().filter(|(&(i, j, k, l), _)| {bpp_mat[&(i, j)] >= min_bpp && bpp_mat[&(k, l)] >= min_bpp}).map(|(&(i, j, k, l), &free_energy)| {((i + T::one(), j + T::one(), k + T::one(), l + T::one()), free_energy)}).collect();
  new_ss_free_energy_mats
}

pub fn write_prob_mat_sets<T>(output_dir_path: &Path, prob_mat_sets: &ProbMatSets<T>, produces_access_probs: bool)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Display + Ord,
{
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let bpp_mat_file_path = output_dir_path.join(BPP_MAT_FILE_NAME);
  let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).unwrap());
  let mut buf_4_writer_2_bpp_mat_file = String::new();
  for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in prob_mats.bpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - T::one(), j - T::one(), bpp));
    }
    buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
  let upp_mat_file_path = output_dir_path.join(UPP_MAT_FILE_NAME);
  let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).unwrap());
  let mut buf_4_writer_2_upp_mat_file = String::new();
  for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
    let seq_len = prob_mats.upp_mat.len();
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (i, &upp) in prob_mats.upp_mat.iter().enumerate() {
      if i == 0 || i == seq_len - 1 {continue;}
      buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
    }
    buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
  let bpp_mat_file_path = output_dir_path.join(BPP_MAT_ON_SS_FILE_NAME);
  let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).unwrap());
  let mut buf_4_writer_2_bpp_mat_file = String::new();
  for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in prob_mats.bpp_mat_on_ss.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i, j, bpp));
    }
    buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
  if produces_access_probs {
    let bpp_mat_file_path = output_dir_path.join(MAX_BPP_MAT_FILE_NAME);
    let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).unwrap());
    let mut buf_4_writer_2_bpp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (&(i, j), &bpp) in prob_mats.max_bpp_mat.iter() {
        buf_4_rna_id.push_str(&format!("{},{},{} ", i - T::one(), j - T::one(), bpp));
      }
      buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
    let bpp_mat_file_path = output_dir_path.join(ACCESS_BPP_MAT_ON_2L_FILE_NAME);
    let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).unwrap());
    let mut buf_4_writer_2_bpp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (&(i, j), &bpp) in prob_mats.access_bpp_mat_4_2l.iter() {
        buf_4_rna_id.push_str(&format!("{},{},{} ", i - T::one(), j - T::one(), bpp));
      }
      buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
    let bpp_mat_file_path = output_dir_path.join(ACCESS_BPP_MAT_ON_ML_FILE_NAME);
    let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).unwrap());
    let mut buf_4_writer_2_bpp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (&(i, j), &bpp) in prob_mats.access_bpp_mat_4_ml.iter() {
        buf_4_rna_id.push_str(&format!("{},{},{} ", i - T::one(), j - T::one(), bpp));
      }
      buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
    let bpp_mat_file_path = output_dir_path.join(BPP_MAT_ON_EL_FILE_NAME);
    let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).unwrap());
    let mut buf_4_writer_2_bpp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (&(i, j), &bpp) in prob_mats.bpp_mat_4_el.iter() {
        buf_4_rna_id.push_str(&format!("{},{},{} ", i - T::one(), j - T::one(), bpp));
      }
      buf_4_writer_2_bpp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_bpp_mat_file.write_all(buf_4_writer_2_bpp_mat_file.as_bytes());
    let upp_mat_file_path = output_dir_path.join(MAX_UPP_MAT_FILE_NAME);
    let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).unwrap());
    let mut buf_4_writer_2_upp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let seq_len = prob_mats.max_upp_mat.len();
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (i, &upp) in prob_mats.max_upp_mat.iter().enumerate() {
        if i == 0 || i == seq_len - 1 {continue;}
        buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
      }
      buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
    let upp_mat_file_path = output_dir_path.join(UPP_MAT_ON_HL_FILE_NAME);
    let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).unwrap());
    let mut buf_4_writer_2_upp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let seq_len = prob_mats.upp_mat_4_hl.len();
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (i, &upp) in prob_mats.upp_mat_4_hl.iter().enumerate() {
        if i == 0 || i == seq_len - 1 {continue;}
        buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
      }
      buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
    let upp_mat_file_path = output_dir_path.join(UPP_MAT_ON_2L_FILE_NAME);
    let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).unwrap());
    let mut buf_4_writer_2_upp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let seq_len = prob_mats.upp_mat_4_2l.len();
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (i, &upp) in prob_mats.upp_mat_4_2l.iter().enumerate() {
        if i == 0 || i == seq_len - 1 {continue;}
        buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
      }
      buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
    let upp_mat_file_path = output_dir_path.join(UPP_MAT_ON_ML_FILE_NAME);
    let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).unwrap());
    let mut buf_4_writer_2_upp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let seq_len = prob_mats.upp_mat_4_ml.len();
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (i, &upp) in prob_mats.upp_mat_4_ml.iter().enumerate() {
        if i == 0 || i == seq_len - 1 {continue;}
        buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
      }
      buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
    let upp_mat_file_path = output_dir_path.join(UPP_MAT_ON_EL_FILE_NAME);
    let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).unwrap());
    let mut buf_4_writer_2_upp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let seq_len = prob_mats.upp_mat_4_el.len();
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (i, &upp) in prob_mats.upp_mat_4_el.iter().enumerate() {
        if i == 0 || i == seq_len - 1 {continue;}
        buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
      }
      buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
  }
}

pub fn get_sparse_pos_sets<T>(bpp_mat_pair: &ProbMatPair<T>, max_gap_num_4_il: T, pseudo_pos_quadruple: &PosQuadruple<T>) -> (PosPairMatSet<T>, PosPairMatSet<T>, PosQuadrupleMat<T>)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let mut forward_pos_pair_mat_set = PosPairMatSet::<T>::default();
  let mut backward_pos_pair_mat_set = PosPairMatSet::<T>::default();
  let mut pos_quadruple_mat = PosQuadrupleMat::<T>::default();
  for pos_pair in bpp_mat_pair.0.keys() {
    for pos_pair_2 in bpp_mat_pair.1.keys() {
      let forward_pos_pair = (pos_pair.0, pos_pair_2.0);
      let backward_pos_pair = (pos_pair.1, pos_pair_2.1);
      if !is_min_gap_ok(&forward_pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {
        continue;
      }
      if !is_min_gap_ok(&backward_pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {
        continue;
      }
      if !forward_pos_pair_mat_set.contains_key(&forward_pos_pair) {
        forward_pos_pair_mat_set.insert(forward_pos_pair, PosPairMat::<T>::default());
      }
      forward_pos_pair_mat_set.get_mut(&forward_pos_pair).unwrap().insert(backward_pos_pair);
      if !backward_pos_pair_mat_set.contains_key(&backward_pos_pair) {
        backward_pos_pair_mat_set.insert(backward_pos_pair, PosPairMat::<T>::default());
      }
      backward_pos_pair_mat_set.get_mut(&backward_pos_pair).unwrap().insert(forward_pos_pair);
      pos_quadruple_mat.insert((pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1));
    }
  }
  (forward_pos_pair_mat_set, backward_pos_pair_mat_set, pos_quadruple_mat)
}
