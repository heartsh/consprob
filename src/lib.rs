extern crate rna_algos;
extern crate getopts;
extern crate scoped_threadpool;
extern crate itertools;
extern crate bio;
extern crate num_cpus;
extern crate hashbrown;

pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use self::scoped_threadpool::Pool;
pub use itertools::multizip;
pub use std::path::Path;
pub use bio::io::fasta::Reader;
pub use std::io::prelude::*;
pub use std::io::BufWriter;
pub use std::fs::File;
pub use std::fs::create_dir;
pub use std::cmp::Ord;
pub use std::marker::{Sync, Send};
pub use hashbrown::HashSet;
pub use rna_algos::utils::*;
pub use rna_algos::mccaskill_algo::*;
pub use rna_algos::durbin_algo::*;
pub use rna_algos::compiled_seq_align_params::*;

pub type PosQuadrupleMat<T> = HashSet<PosQuadruple<T>>;
pub type PosPairMatSet<T> = HashMap<PosPair<T>, PosPairMat<T>>;
pub type PosPairMat<T> = HashSet<PosPair<T>>;
pub type Prob4dMat<T> = HashMap<PosQuadruple<T>, Prob>;
pub type PartFunc4dMat<T> = HashMap<PosQuadruple<T>, PartFunc>;
pub type TmpPartFuncSetMat<T> = HashMap<PosPair<T>, TmpPartFuncs>;

#[derive(Clone)]
pub struct TmpPartFuncs {
  pub part_func_on_sa: PartFunc,
  pub part_func_on_sa_4_ml: PartFunc,
  pub part_func_4_ml: PartFunc,
  pub part_func_4_first_bpas_on_mls: PartFunc,
  pub part_func_4_bpas_on_mls: PartFunc,
  pub part_func_on_mls: PartFunc,
}

#[derive(Clone)]
pub struct StaPartFuncMats<T> {
  pub part_func_4d_mat_4_bpas: PartFunc4dMat<T>,
  pub part_func_4d_mat_4_bpas_accessible_on_els: PartFunc4dMat<T>,
  pub part_func_4d_mat_4_bpas_accessible_on_mls: PartFunc4dMat<T>,
  pub forward_part_func_mat_4_external_loop: SparsePartFuncMat<T>,
  pub backward_part_func_mat_4_external_loop: SparsePartFuncMat<T>,
  pub forward_part_func_mat_4_external_loop_decode: SparsePartFuncMat<T>,
  pub backward_part_func_mat_4_external_loop_decode: SparsePartFuncMat<T>,
  pub forward_tmp_part_func_set_mats_with_pos_pairs: TmpPartFuncSetMatsWithPosPairs<T>,
  pub backward_tmp_part_func_set_mats_with_pos_pairs: TmpPartFuncSetMatsWithPosPairs<T>,
  pub forward_tmp_part_func_set_mats_with_pos_pairs_decode: TmpPartFuncSetMatsWithPosPairs<T>,
  pub backward_tmp_part_func_set_mats_with_pos_pairs_decode: TmpPartFuncSetMatsWithPosPairs<T>,
}
pub type FreeEnergyMat = Vec<FreeEnergies>;

pub struct StaFeParams<T> {
  pub ba_score_mat: SparseFreeEnergyMat<T>,
  pub bpa_score_mat: FreeEnergy4dMat<T>,
  pub insert_scores_range: FreeEnergyMat,
  pub insert_scores_range_4_el: FreeEnergyMat,
  pub insert_scores_range_4_ml: FreeEnergyMat,
  pub insert_scores_range_2: FreeEnergyMat,
  pub insert_scores_range_4_el_2: FreeEnergyMat,
  pub insert_scores_range_4_ml_2: FreeEnergyMat,
}

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
  pub upp_mat_pair_4_hl: ProbSetPair,
  pub upp_mat_pair_4_2l: ProbSetPair,
  pub upp_mat_pair_4_ml: ProbSetPair,
  pub upp_mat_pair_4_el: ProbSetPair,
  pub bpp_mat_pair_2: ProbSetPair,
  pub basepair_align_prob_mat: Prob4dMat<T>,
  pub loop_align_prob_mat: SparseProbMat<T>,
}

#[derive(Clone)]
pub struct PctStaProbMats<T> {
  pub bpp_mat: SparseProbMat<T>,
  pub upp_mat_4_hl: Probs,
  pub upp_mat_4_2l: Probs,
  pub upp_mat_4_ml: Probs,
  pub upp_mat_4_el: Probs,
  pub bpp_mat_2: Probs,
}

pub type SparseProbMatPair<T> = (SparseProbMat<T>, SparseProbMat<T>);
pub type ProbSetPair = (Probs, Probs);
pub type ProbMatSets<T> = Vec<PctStaProbMats<T>>;
pub type StaProbMatsWithRnaIdPairs<T> = HashMap<RnaIdPair, StaProbMats<T>>;
pub type ProbSeqPair<'a> = (&'a Probs, &'a Probs);
pub type Poss<T> = Vec<T>;
pub type TmpPartFuncSetMatsWithPosPairs<T> = HashMap<PosPair<T>, TmpPartFuncSetMat<T>>;

#[derive(Clone)]
pub struct AlignProbMats<T> {
  pub loop_align_prob_mat: SparseProbMat<T>,
  pub basepair_align_prob_mat: Prob4dMat<T>,
  pub align_prob_mat: SparseProbMat<T>,
}

pub type AlignProbMatSetsWithRnaIdPairs<T> = HashMap<RnaIdPair, AlignProbMats<T>>;
pub type SparseProbMatsWithRnaIdPairs<T> = HashMap<RnaIdPair, SparseProbMat<T>>;
pub type PosQuadrupleMatWithLenPairs<T> = HashMap<PosPair<T>, PosPairMat<T>>;
pub type SparsePoss<T> = HashSet<T>;
pub type SparsePosSets<T> = HashMap<T, SparsePoss<T>>;

impl<T: Hash + ToPrimitive + Clone> StaProbMats<T> {
  pub fn origin() -> StaProbMats<T> {
    let prob_mat_pair = (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    let prob_set_pair = (Vec::new(), Vec::new());
    StaProbMats {
      bpp_mat_pair: prob_mat_pair.clone(),
      upp_mat_pair_4_hl: prob_set_pair.clone(),
      upp_mat_pair_4_2l: prob_set_pair.clone(),
      upp_mat_pair_4_ml: prob_set_pair.clone(),
      upp_mat_pair_4_el: prob_set_pair.clone(),
      bpp_mat_pair_2: prob_set_pair,
      basepair_align_prob_mat: Prob4dMat::<T>::default(),
      loop_align_prob_mat: SparseProbMat::<T>::default(),
    }
  }

  pub fn new(seq_len_pair: &PosPair<T>) -> StaProbMats<T> {
    let prob_mat_pair = (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    let prob_set_pair = (vec![NEG_INFINITY; seq_len_pair.0.to_usize().unwrap()], vec![NEG_INFINITY; seq_len_pair.1.to_usize().unwrap()]);
    StaProbMats {
      bpp_mat_pair: prob_mat_pair.clone(),
      upp_mat_pair_4_hl: prob_set_pair.clone(),
      upp_mat_pair_4_2l: prob_set_pair.clone(),
      upp_mat_pair_4_ml: prob_set_pair.clone(),
      upp_mat_pair_4_el: prob_set_pair.clone(),
      bpp_mat_pair_2: prob_set_pair,
      basepair_align_prob_mat: Prob4dMat::<T>::default(),
      loop_align_prob_mat: SparseProbMat::<T>::default(),
    }
  }
}

impl<T: Hash + Clone> PctStaProbMats<T> {
  pub fn origin() -> PctStaProbMats<T> {
    let prob_mat = SparseProbMat::<T>::default();
    let probs = Vec::new();
    PctStaProbMats {
      bpp_mat: prob_mat,
      upp_mat_4_hl: probs.clone(),
      upp_mat_4_2l: probs.clone(),
      upp_mat_4_ml: probs.clone(),
      upp_mat_4_el: probs.clone(),
      bpp_mat_2: probs,
    }
  }

  pub fn new(seq_len: usize) -> PctStaProbMats<T> {
    let prob_mat = SparseProbMat::<T>::default();
    let probs = vec![0.; seq_len as usize];
    PctStaProbMats {
      bpp_mat: prob_mat,
      upp_mat_4_hl: probs.clone(),
      upp_mat_4_2l: probs.clone(),
      upp_mat_4_ml: probs.clone(),
      upp_mat_4_el: probs.clone(),
      bpp_mat_2: probs,
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
      insert_scores_range: insert_scores.clone(),
      insert_scores_range_4_el: insert_scores.clone(),
      insert_scores_range_4_ml: insert_scores.clone(),
      insert_scores_range_2: insert_scores.clone(),
      insert_scores_range_4_el_2: insert_scores.clone(),
      insert_scores_range_4_ml_2: insert_scores,
    }
  }

  pub fn new(seq_pair: &SeqPair, seq_len_pair: &PosPair<T>, pos_quadruple_mat: &PosQuadrupleMat<T>, align_prob_mat: &SparseProbMat<T>, uses_contra_model: bool, align_feature_score_sets: &AlignFeatureCountSets,) -> StaFeParams<T> {
    let mut sta_fe_params = StaFeParams::<T>::origin();
    let mat = vec![vec![NEG_INFINITY; seq_len_pair.0.to_usize().unwrap()]; seq_len_pair.0.to_usize().unwrap()];
    sta_fe_params.insert_scores_range = mat.clone();
    sta_fe_params.insert_scores_range_4_el = mat.clone();
    sta_fe_params.insert_scores_range_4_ml = mat;
    let mat = vec![vec![NEG_INFINITY; seq_len_pair.1.to_usize().unwrap()]; seq_len_pair.1.to_usize().unwrap()];
    sta_fe_params.insert_scores_range_2 = mat.clone();
    sta_fe_params.insert_scores_range_4_el_2 = mat.clone();
    sta_fe_params.insert_scores_range_4_ml_2 = mat.clone();
    for i in range(T::one(), seq_len_pair.1 - T::one()) {
      let long_i = i.to_usize().unwrap();
      let base = seq_pair.1[long_i];
      let mut sum = align_feature_score_sets.insert_counts[base];
      let mut sum_4_el = sum + if !uses_contra_model {0.} else {CONTRA_EL_UNPAIRED_FE};
      let mut sum_4_ml = sum + if !uses_contra_model {0.} else {CONTRA_ML_UNPAIRED_FE};
      sta_fe_params.insert_scores_range_2[long_i][long_i] = sum;
      sta_fe_params.insert_scores_range_4_el_2[long_i][long_i] = sum_4_el;
      sta_fe_params.insert_scores_range_4_ml_2[long_i][long_i] = sum_4_ml;
      for j in range(i + T::one(), seq_len_pair.1 - T::one()) {
        let long_j = j.to_usize().unwrap();
        let base = seq_pair.1[long_j];
        let term = align_feature_score_sets.insert_counts[base] + align_feature_score_sets.insert_extend_count;
        sum += term;
        sum_4_el += term + if !uses_contra_model {0.} else {CONTRA_EL_UNPAIRED_FE};
        sum_4_ml += term + if !uses_contra_model {0.} else {CONTRA_ML_UNPAIRED_FE};
        sta_fe_params.insert_scores_range_2[long_i][long_j] = sum;
        sta_fe_params.insert_scores_range_4_el_2[long_i][long_j] = sum_4_el;
        sta_fe_params.insert_scores_range_4_ml_2[long_i][long_j] = sum_4_ml;
      }
    }
    for i in range(T::one(), seq_len_pair.0 - T::one()) {
      let long_i = i.to_usize().unwrap();
      let base = seq_pair.0[long_i];
      let term = align_feature_score_sets.insert_counts[base];
      let mut sum = term;
      let mut sum_4_el = sum + if !uses_contra_model {0.} else {CONTRA_EL_UNPAIRED_FE};
      let mut sum_4_ml = sum + if !uses_contra_model {0.} else {CONTRA_ML_UNPAIRED_FE};
      sta_fe_params.insert_scores_range[long_i][long_i] = sum;
      sta_fe_params.insert_scores_range_4_el[long_i][long_i] = sum_4_el;
      sta_fe_params.insert_scores_range_4_ml[long_i][long_i] = sum_4_ml;
      for j in range(i + T::one(), seq_len_pair.0 - T::one()) {
        let long_j = j.to_usize().unwrap();
        let base = seq_pair.0[long_j];
        let term = align_feature_score_sets.insert_counts[base] + align_feature_score_sets.insert_extend_count;
        sum += term;
        sum_4_el += term + if !uses_contra_model {0.} else {CONTRA_EL_UNPAIRED_FE};
        sum_4_ml += term + if !uses_contra_model {0.} else {CONTRA_ML_UNPAIRED_FE};
        sta_fe_params.insert_scores_range[long_i][long_j] = sum;
        sta_fe_params.insert_scores_range_4_el[long_i][long_j] = sum_4_el;
        sta_fe_params.insert_scores_range_4_ml[long_i][long_j] = sum_4_ml;
      }
    }
    for pos_pair in align_prob_mat.keys() {
      let &(i, j) = pos_pair;
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let base_pair = (seq_pair.0[long_i], seq_pair.1[long_j]);
      sta_fe_params.ba_score_mat.insert(*pos_pair, align_feature_score_sets.align_count_mat[base_pair.0][base_pair.1]);
    }
    for &(i, j, k, l) in pos_quadruple_mat {
      let (long_i, long_j, long_k, long_l) = (i.to_usize().unwrap(), j.to_usize().unwrap(), k.to_usize().unwrap(), l.to_usize().unwrap());
      let pos_quadruple = (i, j, k, l);
      let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
      let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
      sta_fe_params.bpa_score_mat.insert(pos_quadruple, align_feature_score_sets.align_count_mat[base_pair.0][base_pair_2.0] + align_feature_score_sets.align_count_mat[base_pair.1][base_pair_2.1]);
    }
    sta_fe_params
  }
}

impl TmpPartFuncs {
  pub fn new() -> TmpPartFuncs {
    TmpPartFuncs {
      part_func_on_sa: NEG_INFINITY,
      part_func_on_sa_4_ml: NEG_INFINITY,
      part_func_4_ml: NEG_INFINITY,
      part_func_4_first_bpas_on_mls: NEG_INFINITY,
      part_func_4_bpas_on_mls: NEG_INFINITY,
      part_func_on_mls: NEG_INFINITY,
    }
  }
}

impl<T: Hash + Clone> StaPartFuncMats<T> {
  pub fn new() -> StaPartFuncMats<T> {
    let part_func_4d_mat = PartFunc4dMat::<T>::default();
    let part_func_mat = SparsePartFuncMat::<T>::default();
    let tmp_part_func_set_mats_with_pos_pairs = TmpPartFuncSetMatsWithPosPairs::<T>::default();
    StaPartFuncMats {
      part_func_4d_mat_4_bpas: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_els: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_mls: part_func_4d_mat,
      forward_part_func_mat_4_external_loop: part_func_mat.clone(),
      backward_part_func_mat_4_external_loop: part_func_mat.clone(),
      forward_part_func_mat_4_external_loop_decode: part_func_mat.clone(),
      backward_part_func_mat_4_external_loop_decode: part_func_mat,
      forward_tmp_part_func_set_mats_with_pos_pairs: tmp_part_func_set_mats_with_pos_pairs.clone(),
      backward_tmp_part_func_set_mats_with_pos_pairs: tmp_part_func_set_mats_with_pos_pairs.clone(),
      forward_tmp_part_func_set_mats_with_pos_pairs_decode: tmp_part_func_set_mats_with_pos_pairs.clone(),
      backward_tmp_part_func_set_mats_with_pos_pairs_decode: tmp_part_func_set_mats_with_pos_pairs,
    }
  }
}

impl<T: Hash + Clone> AlignProbMats<T> {
  pub fn new() -> AlignProbMats<T> {
    AlignProbMats {
      loop_align_prob_mat: SparseProbMat::<T>::default(),
      basepair_align_prob_mat: Prob4dMat::<T>::default(),
      align_prob_mat: SparseProbMat::<T>::default(),
    }
  }
}

pub const DEFAULT_MIN_BPP: Prob = 0.01;
pub const DEFAULT_MIN_ALIGN_PROB: Prob = 0.01;
pub const BPP_MAT_FILE_NAME: &'static str = "bpp_mats.dat";
pub const UPP_MAT_ON_HL_FILE_NAME: &'static str = "upp_mats_on_hl.dat";
pub const UPP_MAT_ON_2L_FILE_NAME: &'static str = "upp_mats_on_2l.dat";
pub const UPP_MAT_ON_ML_FILE_NAME: &'static str = "upp_mats_on_ml.dat";
pub const UPP_MAT_ON_EL_FILE_NAME: &'static str = "upp_mats_on_el.dat";
pub const BPP_MAT_FILE_NAME_2: &'static str = "bpp_mats_2.dat";
pub const BASEPAIR_ALIGN_PROB_MAT_FILE_NAME: &'static str = "basepair_align_prob_mat.dat";
pub const LOOP_ALIGN_PROB_MAT_FILE_NAME: &'static str = "loop_align_prob_mat.dat";
pub const ALIGN_PROB_MAT_FILE_NAME: &'static str = "align_prob_mat.dat";
pub const README_FILE_NAME: &str = "README.md";
pub const README_CONTENTS: &str = "# bpp_mats.dat\n
This file contains average probabilistic consistency based on posterior nucleotide pair-matching probabilities. You can treat this average probabilistic consistency like conventional nucleotide base-pairing probabilities. Nucleotide positions are indexed starting from zero.\n\n
# bpp_mats_2.dat\n
This file contains average probabilistic consistency per nucleotide. This average probabilistic consistency is obtained by marginalizing each nucleotide for average probabilistic consistency in \"bpp_mats.dat.\"\n\n
# upp_mats_on_x.dat\n
This file type contains average probabilistic consistency per nucleotide. This average probabilistic consistency is for nucleotide unpairing and under the structural context \"x.\" \"hl,\" \"2l,\" \"ml,\" \"el\" stand for hairpin loops, 2-loops, multi-loops, external loops, respectively.\n\n
# basepair_align_prob_mat.dat\n
This file contains posterior nucleotide pair-matching probabilities.\n\n
# loop_align_prob_mat.dat\n
This file contains posterior nucleotide loop-matching probabilities.\n\n
align_prob_mat.dat\n
This file contains posterior nucleotide matching probabilities.";
pub const INSERT_2_MATCH_SCORE: FreeEnergy = MATCH_2_INSERT_SCORE;

pub fn io_algo_4_prob_mats<T>(seq_len_pair: &PosPair<T>, sta_fe_params: &StaFeParams<T>, max_bp_span_pair: &PosPair<T>, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair<T>, produces_struct_profs: bool, forward_pos_pair_mat_set: &PosPairMatSet<T>, backward_pos_pair_mat_set: &PosPairMatSet<T>, uses_contra_model: bool, pos_quadruple_mat_with_len_pairs: &PosQuadrupleMatWithLenPairs<T>, produces_align_probs: bool, matchable_pos_sets_1: &SparsePosSets<T>, matchable_pos_sets_2: &SparsePosSets<T>, align_feature_score_sets: &AlignFeatureCountSets,) -> StaProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let (sta_part_func_mats, global_part_func) = 
    get_sta_inside_part_func_mats::<T>(seq_len_pair, sta_fe_params, max_bp_span_pair, ss_free_energy_mat_set_pair, forward_pos_pair_mat_set, backward_pos_pair_mat_set, uses_contra_model, pos_quadruple_mat_with_len_pairs, matchable_pos_sets_1, matchable_pos_sets_2, align_feature_score_sets,);
  get_sta_prob_mats::<T>(seq_len_pair, sta_fe_params, max_bp_span_pair, &sta_part_func_mats, ss_free_energy_mat_set_pair, produces_struct_profs, global_part_func, uses_contra_model, pos_quadruple_mat_with_len_pairs, produces_align_probs, forward_pos_pair_mat_set, backward_pos_pair_mat_set, matchable_pos_sets_1, matchable_pos_sets_2, align_feature_score_sets,)
}

pub fn get_sta_inside_part_func_mats<T>(seq_len_pair: &PosPair<T>, sta_fe_params: &StaFeParams<T>, max_bp_span_pair: &PosPair<T>, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair<T>, forward_pos_pair_mat_set: &PosPairMatSet<T>, backward_pos_pair_mat_set: &PosPairMatSet<T>, uses_contra_model: bool, pos_quadruple_mat_with_len_pairs: &PosQuadrupleMatWithLenPairs<T>, matchable_pos_sets_1: &SparsePosSets<T>, matchable_pos_sets_2: &SparsePosSets<T>, align_feature_score_sets: &AlignFeatureCountSets,) -> (StaPartFuncMats<T>, PartFunc)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let mut sta_part_func_mats = StaPartFuncMats::<T>::new();
  for substr_len_1 in range_inclusive(
    T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(),
    max_bp_span_pair.0,
  ) {
    for substr_len_2 in range_inclusive(
      T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(),
      max_bp_span_pair.1,
    ) {
      match pos_quadruple_mat_with_len_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let (j, l) = (i + substr_len_1 - T::one(), k + substr_len_2 - T::one());
            let (long_i, long_j, long_k, long_l) = (i.to_usize().unwrap(), j.to_usize().unwrap(), k.to_usize().unwrap(), l.to_usize().unwrap());
            let pos_quadruple = (i, j, k, l);
            let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
            let (part_func_on_sa, part_func_4_ml) = get_tmp_part_func_set_mats::<T>(sta_fe_params, &pos_quadruple, &mut sta_part_func_mats, true, forward_pos_pair_mat_set, uses_contra_model, matchable_pos_sets_1, matchable_pos_sets_2, align_feature_score_sets);
            let _ = get_tmp_part_func_set_mats::<T>(sta_fe_params, &pos_quadruple, &mut sta_part_func_mats, false, backward_pos_pair_mat_set, uses_contra_model, matchable_pos_sets_1, matchable_pos_sets_2, align_feature_score_sets);
            let mut sum = NEG_INFINITY;
            if !uses_contra_model || (uses_contra_model && long_j - long_i - 1 <= CONTRA_MAX_LOOP_LEN && long_l - long_k - 1 <= CONTRA_MAX_LOOP_LEN) {
              let score = bpa_score + ss_free_energy_mat_set_pair.0.hl_fe_mat[&(i, j)] + ss_free_energy_mat_set_pair.1.hl_fe_mat[&(k, l)] + part_func_on_sa;
              logsumexp(&mut sum, score);
            }
            let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs_decode[&(i, k)];
            let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs_decode[&(j, l)];
            let min = T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap();
            let min_len_pair = (
              if substr_len_1 <= min + T::from_usize(MAX_2_LOOP_LEN + 2).unwrap() {
                min
              } else {
                substr_len_1 - T::from_usize(MAX_2_LOOP_LEN + 2).unwrap()
              },
              if substr_len_2 <= min + T::from_usize(MAX_2_LOOP_LEN + 2).unwrap() {
                min
              } else {
                substr_len_2 - T::from_usize(MAX_2_LOOP_LEN + 2).unwrap()
              },
              );
            for substr_len_3 in range(
              min_len_pair.0,
              substr_len_1 - T::one(),
            ) {
              for substr_len_4 in range(
                min_len_pair.1,
                substr_len_2 - T::one(),
              ) {
                match pos_quadruple_mat_with_len_pairs.get(&(substr_len_3, substr_len_4)) {
                  Some(pos_pairs_2) => {
                    for &(m, o) in pos_pairs_2 {
                      let (n, p) = (m + substr_len_3 - T::one(), o + substr_len_4 - T::one());
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
                          let pos_pair_2 = (m - T::one(), o - T::one());
                          match forward_tmp_part_func_set_mat.get(&pos_pair_2) {
                            Some(part_funcs) => {
                              let part_func = part_funcs.part_func_on_sa;
                              let term = part_func;
                              logsumexp(&mut forward_term, term);
                            }, None => {},
                          }
                          let pos_pair_2 = (n + T::one(), p + T::one());
                          match backward_tmp_part_func_set_mat.get(&pos_pair_2) {
                            Some(part_funcs) => {
                              let part_func = part_funcs.part_func_on_sa;
                              let term = part_func;
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
                  }, None => {},
                }
              }
            }
            let multi_loop_closing_basepairing_fe = ss_free_energy_mat_set_pair.0.ml_closing_bp_fe_mat[&(i, j)];
            let multi_loop_closing_basepairing_fe_2 = ss_free_energy_mat_set_pair.1.ml_closing_bp_fe_mat[&(k, l)];
            let score = bpa_score + multi_loop_closing_basepairing_fe + multi_loop_closing_basepairing_fe_2 + part_func_4_ml;
            logsumexp(&mut sum, score);
            if sum > NEG_INFINITY {
              sta_part_func_mats.part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
              let accessible_basepairing_fe = ss_free_energy_mat_set_pair.0.accessible_bp_fe_mat[&(i, j)];
              let accessible_basepairing_fe_2 = ss_free_energy_mat_set_pair.1.accessible_bp_fe_mat[&(k, l)];
              sum += accessible_basepairing_fe + accessible_basepairing_fe_2;
              sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.insert(pos_quadruple, sum + if !uses_contra_model {0.} else {2. * CONTRA_EL_PAIRED_FE});
              sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, sum + 2. * if !uses_contra_model {COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE} else {CONTRA_ML_PAIRED_FE});
            }
          }
        }, None => {},
      }
    }
  }
  let leftmost_pos_pair = (T::zero(), T::zero());
  let rightmost_pos_pair = (seq_len_pair.0 - T::one(), seq_len_pair.1 - T::one());
  sta_part_func_mats.forward_part_func_mat_4_external_loop.insert(leftmost_pos_pair, 0.);
  sta_part_func_mats.backward_part_func_mat_4_external_loop.insert(rightmost_pos_pair, 0.);
  for i in range(T::zero(), seq_len_pair.0 - T::one()) {
    for j in range(T::zero(), seq_len_pair.1 - T::one()) {
      let pos_pair = (i, j);
      if pos_pair == leftmost_pos_pair {continue;}
      let mut sum = NEG_INFINITY;
      match forward_pos_pair_mat_set.get(&pos_pair) {
        Some(forward_pos_pair_mat) => {
          for &(k, l) in forward_pos_pair_mat {
            let pos_pair_2 = (k - T::one(), l - T::one());
            let pos_quadruple = (k, i, l, j);
            match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_part_func_mats.forward_part_func_mat_4_external_loop_decode.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    let score = part_func_2 + part_func;
                    logsumexp(&mut sum, score);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      if i > T::zero() && j > T::zero() {
        match sta_fe_params.ba_score_mat.get(&pos_pair) {
          Some(&ba_score) => {
            let mut sum_2 = NEG_INFINITY;
            let ba_score = ba_score + if !uses_contra_model {0.} else {2. * CONTRA_EL_UNPAIRED_FE};
            let pos_pair_2 = (i - T::one(), j - T::one());
            let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
            let is_begin = pos_pair_2 == leftmost_pos_pair;
            match sta_part_func_mats.forward_part_func_mat_4_external_loop.get(&pos_pair_2) {
              Some(&part_func) => {
                let score = part_func + if is_begin {align_feature_score_sets.init_match_count} else {align_feature_score_sets.match_2_match_count};
                logsumexp(&mut sum_2, score);
              }, None => {},
            }
            match matchable_pos_sets_1.get(&pos_pair_2.0) {
              Some(matchable_poss) => {
                for &x in matchable_poss {
                  if x >= pos_pair_2.1 {continue;}
                  let pos_pair_3 = (pos_pair_2.0, x);
                  match sta_part_func_mats.forward_part_func_mat_4_external_loop.get(&pos_pair_3) {
                    Some(&part_func) => {
                      let long_x = x.to_usize().unwrap();
                      let is_begin = pos_pair_3 == leftmost_pos_pair;
                      let insert_score_range = sta_fe_params.insert_scores_range_4_el_2[long_x + 1][long_pos_pair_2.1] + if is_begin {align_feature_score_sets.init_insert_count} else {align_feature_score_sets.match_2_insert_count};
                      let score = part_func + insert_score_range + align_feature_score_sets.match_2_insert_count;
                      logsumexp(&mut sum_2, score);
                    }, None => {},
                  }
                }
              }, None => {},
            }
            match matchable_pos_sets_2.get(&pos_pair_2.1) {
              Some(matchable_poss) => {
                for &x in matchable_poss {
                  if x >= pos_pair_2.0 {continue;}
                  let pos_pair_3 = (x, pos_pair_2.1);
                  match sta_part_func_mats.forward_part_func_mat_4_external_loop.get(&pos_pair_3) {
                    Some(&part_func) => {
                      let long_x = x.to_usize().unwrap();
                      let is_begin = pos_pair_3 == leftmost_pos_pair;
                      let insert_score_range = sta_fe_params.insert_scores_range_4_el[long_x + 1][long_pos_pair_2.0] + if is_begin {align_feature_score_sets.init_insert_count} else {align_feature_score_sets.match_2_insert_count};
                      let score = part_func + insert_score_range + align_feature_score_sets.match_2_insert_count;
                      logsumexp(&mut sum_2, score);
                    }, None => {},
                  }
                }
              }, None => {},
            }
            if sum_2 > NEG_INFINITY {
              sta_part_func_mats.forward_part_func_mat_4_external_loop_decode.insert(pos_pair_2, sum_2);
            }
            let term = sum_2 + ba_score;
            logsumexp(&mut sum, term);
            if sum > NEG_INFINITY {
              sta_part_func_mats.forward_part_func_mat_4_external_loop.insert(pos_pair, sum);
            }
          }, None => {},
        }
      }
    }
  }
  let mut final_sum = NEG_INFINITY;
  let pos_pair_2 = (rightmost_pos_pair.0 - T::one(), rightmost_pos_pair.1 - T::one());
  let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
  match sta_part_func_mats.forward_part_func_mat_4_external_loop.get(&pos_pair_2) {
    Some(&part_func) => {
      logsumexp(&mut final_sum, part_func);
    }, None => {},
  }
  match matchable_pos_sets_1.get(&pos_pair_2.0) {
    Some(matchable_poss) => {
      for &x in matchable_poss {
        if x >= pos_pair_2.1 {continue;}
        match sta_part_func_mats.forward_part_func_mat_4_external_loop.get(&(pos_pair_2.0, x)) {
          Some(&part_func) => {
            let long_x = x.to_usize().unwrap();
            let insert_score_range = sta_fe_params.insert_scores_range_4_el_2[long_x + 1][long_pos_pair_2.1] + align_feature_score_sets.match_2_insert_count;
            let score = part_func + insert_score_range;
            logsumexp(&mut final_sum, score);
          }, None => {},
        }
      }
    }, None => {},
  }
  match matchable_pos_sets_2.get(&pos_pair_2.1) {
    Some(matchable_poss) => {
      for &x in matchable_poss {
        if x >= pos_pair_2.0 {continue;}
        match sta_part_func_mats.forward_part_func_mat_4_external_loop.get(&(x, pos_pair_2.1)) {
          Some(&part_func) => {
            let long_x = x.to_usize().unwrap();
            let insert_score_range = sta_fe_params.insert_scores_range_4_el[long_x + 1][long_pos_pair_2.0] + align_feature_score_sets.match_2_insert_count;
            let score = part_func + insert_score_range;
            logsumexp(&mut final_sum, score);
          }, None => {},
        }
      }
    }, None => {},
  }
  for i in range(T::one(), seq_len_pair.0).rev() {
    for j in range(T::one(), seq_len_pair.1).rev() {
      let pos_pair = (i, j);
      if pos_pair == rightmost_pos_pair {continue;}
      let mut sum = NEG_INFINITY;
      match backward_pos_pair_mat_set.get(&pos_pair) {
        Some(backward_pos_pair_mat) => {
          for &(k, l) in backward_pos_pair_mat {
            let pos_pair_2 = (k + T::one(), l + T::one());
            let pos_quadruple = (i, k, j, l);
            match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
              Some(&part_func) => {
                match sta_part_func_mats.backward_part_func_mat_4_external_loop_decode.get(&pos_pair_2) {
                  Some(&part_func_2) => {
                    let score = part_func_2 + part_func;
                    logsumexp(&mut sum, score);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      if i < seq_len_pair.0 - T::one() && j < seq_len_pair.1 - T::one() {
        match sta_fe_params.ba_score_mat.get(&pos_pair) {
          Some(&ba_score) => {
            let mut sum_2 = NEG_INFINITY;
            let ba_score = ba_score + if !uses_contra_model {0.} else {2. * CONTRA_EL_UNPAIRED_FE};
            let pos_pair_2 = (i + T::one(), j + T::one());
            let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
            let is_end = pos_pair_2 == rightmost_pos_pair;
            match sta_part_func_mats.backward_part_func_mat_4_external_loop.get(&pos_pair_2) {
              Some(&part_func) => {
                let score = part_func + if is_end {0.} else {align_feature_score_sets.match_2_match_count};
                logsumexp(&mut sum_2, score);
              }, None => {},
            }
            match matchable_pos_sets_1.get(&pos_pair_2.0) {
              Some(matchable_poss) => {
                for &x in matchable_poss {
                  if x <= pos_pair_2.1 {continue;}
                  let pos_pair_3 = (pos_pair_2.0, x);
                  match sta_part_func_mats.backward_part_func_mat_4_external_loop.get(&pos_pair_3) {
                    Some(&part_func) => {
                      let long_x = x.to_usize().unwrap();
                      let is_end = pos_pair_3 == rightmost_pos_pair;
                      let insert_score_range = sta_fe_params.insert_scores_range_4_el_2[long_pos_pair_2.1][long_x - 1] + if is_end {0.} else {align_feature_score_sets.match_2_insert_count};
                      let score = part_func + insert_score_range + align_feature_score_sets.match_2_insert_count;
                      logsumexp(&mut sum_2, score);
                    }, None => {},
                  }
                }
              }, None => {},
            }
            match matchable_pos_sets_2.get(&pos_pair_2.1) {
              Some(matchable_poss) => {
                for &x in matchable_poss {
                  if x <= pos_pair_2.0 {continue;}
                  let pos_pair_3 = (x, pos_pair_2.1);
                  match sta_part_func_mats.backward_part_func_mat_4_external_loop.get(&pos_pair_3) {
                    Some(&part_func) => {
                      let long_x = x.to_usize().unwrap();
                      let is_end = pos_pair_3 == rightmost_pos_pair;
                      let insert_score_range = sta_fe_params.insert_scores_range_4_el[long_pos_pair_2.0][long_x - 1] + if is_end {0.} else {align_feature_score_sets.match_2_insert_count};
                      let score = part_func + insert_score_range + align_feature_score_sets.match_2_insert_count;
                      logsumexp(&mut sum_2, score);
                    }, None => {},
                  }
                }
              }, None => {},
            }
            if sum_2 > NEG_INFINITY {
              sta_part_func_mats.backward_part_func_mat_4_external_loop_decode.insert(pos_pair_2, sum_2);
            }
            let term = sum_2 + ba_score;
            logsumexp(&mut sum, term);
            if sum > NEG_INFINITY {
              sta_part_func_mats.backward_part_func_mat_4_external_loop.insert(pos_pair, sum);
            }
          }, None => {},
        }
      }
    }
  }
  (sta_part_func_mats, final_sum)
}

pub fn get_tmp_part_func_set_mats<T>(sta_fe_params: &StaFeParams<T>, pos_quadruple: &PosQuadruple<T>, sta_part_func_mats: &mut StaPartFuncMats<T>, is_forward: bool, pos_pair_mat_set: &PosPairMatSet<T>, uses_contra_model: bool, matchable_pos_sets_1: &SparsePosSets<T>, matchable_pos_sets_2: &SparsePosSets<T>, align_feature_score_sets: &AlignFeatureCountSets,) -> (PartFunc, PartFunc)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let &(i, j, k, l) = pos_quadruple;
  let leftmost_pos_pair = if is_forward {(i, k)} else {(i + T::one(), k + T::one())};
  let rightmost_pos_pair = if is_forward {(j - T::one(), l - T::one())} else {(j, l)};
  let tmp_part_func_set_mats_with_pos_pairs = if is_forward{&mut sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs} else {&mut sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs};
  let tmp_part_func_set_mats_with_pos_pairs_decode = if is_forward{&mut sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs_decode} else {&mut sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs_decode};
  if !tmp_part_func_set_mats_with_pos_pairs.contains_key(&if is_forward {
    leftmost_pos_pair
  } else {
    rightmost_pos_pair
  }) {
    tmp_part_func_set_mats_with_pos_pairs.insert(if is_forward {
      leftmost_pos_pair
      } else {
        rightmost_pos_pair
      }, TmpPartFuncSetMat::<T>::new()
    );
  }
  if !tmp_part_func_set_mats_with_pos_pairs_decode.contains_key(&if is_forward {
    leftmost_pos_pair
  } else {
    rightmost_pos_pair
  }) {
    tmp_part_func_set_mats_with_pos_pairs_decode.insert(if is_forward {
      leftmost_pos_pair
      } else {
        rightmost_pos_pair
      }, TmpPartFuncSetMat::<T>::new()
    );
  }
  let ref mut tmp_part_func_set_mat = tmp_part_func_set_mats_with_pos_pairs.get_mut(&if is_forward {
    leftmost_pos_pair
  } else {
    rightmost_pos_pair
  }).unwrap();
  let ref mut tmp_part_func_set_mat_decode = tmp_part_func_set_mats_with_pos_pairs_decode.get_mut(&if is_forward {
    leftmost_pos_pair
  } else {
    rightmost_pos_pair
  }).unwrap();
  let iter: Poss<T> = if is_forward {range(i, j).collect()} else {range_inclusive(i + T::one(), j).rev().collect()};
  let iter_2: Poss<T> = if is_forward {range(k, l).collect()} else {range_inclusive(k + T::one(), l).rev().collect()};
  for &u in iter.iter() {
    for &v in iter_2.iter() {
      let pos_pair = (u, v);
      if tmp_part_func_set_mat.contains_key(&pos_pair) {
        continue;
      }
      let mut tmp_part_funcs = TmpPartFuncs::new();
      if (is_forward && u == i && v == k) || (!is_forward && u == j && v == l) {
        tmp_part_funcs.part_func_on_sa = 0.;
        tmp_part_funcs.part_func_on_sa_4_ml = 0.;
        tmp_part_funcs.part_func_on_mls = 0.;
        tmp_part_func_set_mat.insert(pos_pair, tmp_part_funcs);
        continue;
      }
      let mut sum_4_ml = NEG_INFINITY;
      let mut sum_4_first_bpas_on_mls = sum_4_ml;
      let mut tmp_sum = sum_4_ml;
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
                match tmp_part_func_set_mat_decode.get(&pos_pair_2) {
                  Some(part_funcs) => {
                    let part_func_2 = part_funcs.part_func_4_bpas_on_mls;
                    let score = part_func_2 + part_func;
                    logsumexp(&mut sum_4_ml, score);
                    let part_func_2 = if !uses_contra_model {part_funcs.part_func_on_sa} else {part_funcs.part_func_on_sa_4_ml};
                    let score = part_func_2 + part_func;
                    logsumexp(&mut sum_4_first_bpas_on_mls, score);
                  }, None => {},
                }
              }, None => {},
            }
          }
        }, None => {},
      }
      let pos_pair_2 = if is_forward {(u - T::one(), v - T::one())} else {(u + T::one(), v + T::one())};
      let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
      match sta_fe_params.ba_score_mat.get(&pos_pair) {
        Some(&ba_score) => {
          let mut tmp_part_funcs_decode = TmpPartFuncs::new();
          let mut sum_on_sa_decode = NEG_INFINITY;
          let mut sum_on_sa_4_ml_decode = sum_on_sa_decode;
          let mut sum_4_ml_decode = sum_on_sa_decode;
          let mut sum_4_first_bpas_on_mls_decode = sum_on_sa_decode;
          let mut tmp_sum_decode = sum_on_sa_decode;
          let ba_score_4_ml = ba_score + if !uses_contra_model {0.} else {2. * CONTRA_ML_UNPAIRED_FE};
          match tmp_part_func_set_mat.get(&pos_pair_2) {
            Some(part_funcs) => {
              let part_func = part_funcs.part_func_4_ml;
              let score = part_func + align_feature_score_sets.match_2_match_count;
              logsumexp(&mut sum_4_ml_decode, score);
              let part_func = part_funcs.part_func_4_first_bpas_on_mls;
              let score = part_func + align_feature_score_sets.match_2_match_count;
              logsumexp(&mut sum_4_first_bpas_on_mls_decode, score);
              let part_func = part_funcs.part_func_on_sa;
              let score = part_func + align_feature_score_sets.match_2_match_count;
              logsumexp(&mut sum_on_sa_decode, score);
              if uses_contra_model {
                let part_func = part_funcs.part_func_on_sa_4_ml;
                let score = part_func + align_feature_score_sets.match_2_match_count;
                logsumexp(&mut sum_on_sa_4_ml_decode, score);
              }
            }, None => {},
          }
          match matchable_pos_sets_1.get(&pos_pair_2.0) {
            Some(matchable_poss) => {
              for &x in matchable_poss {
                if is_forward && x >= pos_pair_2.1 || (!is_forward && x <= pos_pair_2.1) {continue;}
                match tmp_part_func_set_mat.get(&(pos_pair_2.0, x)) {
                  Some(part_funcs) => {
                    let long_x = x.to_usize().unwrap();
                    let insert_score_range = if is_forward {
                      sta_fe_params.insert_scores_range_2[long_x + 1][long_pos_pair_2.1]
                    } else {
                      sta_fe_params.insert_scores_range_2[long_pos_pair_2.1][long_x - 1]
                    } + align_feature_score_sets.match_2_insert_count;
                    let insert_score_range_4_ml = if is_forward {
                      sta_fe_params.insert_scores_range_4_ml_2[long_x + 1][long_pos_pair_2.1]
                    } else {
                      sta_fe_params.insert_scores_range_4_ml_2[long_pos_pair_2.1][long_x - 1]
                    } + align_feature_score_sets.match_2_insert_count;
                    let part_func = part_funcs.part_func_4_ml;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
                    logsumexp(&mut sum_4_ml_decode, score);
                    let part_func = part_funcs.part_func_4_first_bpas_on_mls;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
                    logsumexp(&mut sum_4_first_bpas_on_mls_decode, score);
                    let part_func = part_funcs.part_func_on_sa;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range;
                    logsumexp(&mut sum_on_sa_decode, score);
                    if uses_contra_model {
                      let part_func = part_funcs.part_func_on_sa_4_ml;
                      let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
                      logsumexp(&mut sum_on_sa_4_ml_decode, score);
                    }
                  }, None => {},
                }
              }
            }, None => {},
          }
          match matchable_pos_sets_2.get(&pos_pair_2.1) {
            Some(matchable_poss) => {
              for &x in matchable_poss {
                if is_forward && x >= pos_pair_2.0 || (!is_forward && x <= pos_pair_2.0) {continue;}
                match tmp_part_func_set_mat.get(&(x, pos_pair_2.1)) {
                  Some(part_funcs) => {
                    let long_x = x.to_usize().unwrap();
                    let insert_score_range = if is_forward {
                      sta_fe_params.insert_scores_range[long_x + 1][long_pos_pair_2.0]
                    } else {
                      sta_fe_params.insert_scores_range[long_pos_pair_2.0][long_x - 1]
                    } + align_feature_score_sets.match_2_insert_count;
                    let insert_score_range_4_ml = if is_forward {
                      sta_fe_params.insert_scores_range_4_ml[long_x + 1][long_pos_pair_2.0]
                    } else {
                      sta_fe_params.insert_scores_range_4_ml[long_pos_pair_2.0][long_x - 1]
                    } + align_feature_score_sets.match_2_insert_count;
                    let part_func = part_funcs.part_func_4_ml;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
                    logsumexp(&mut sum_4_ml_decode, score);
                    let part_func = part_funcs.part_func_4_first_bpas_on_mls;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
                    logsumexp(&mut sum_4_first_bpas_on_mls_decode, score);
                    let part_func = part_funcs.part_func_on_sa;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range;
                    logsumexp(&mut sum_on_sa_decode, score);
                    if uses_contra_model {
                      let part_func = part_funcs.part_func_on_sa_4_ml;
                      let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
                      logsumexp(&mut sum_on_sa_4_ml_decode, score);
                    }
                  }, None => {},
                }
              }
            }, None => {},
          }
          tmp_part_funcs_decode.part_func_4_ml = sum_4_ml_decode;
          logsumexp(&mut tmp_sum_decode, sum_4_ml_decode);
          tmp_part_funcs_decode.part_func_4_first_bpas_on_mls = sum_4_first_bpas_on_mls_decode;
          logsumexp(&mut tmp_sum_decode, sum_4_first_bpas_on_mls_decode);
          tmp_part_funcs_decode.part_func_4_bpas_on_mls = tmp_sum_decode;
          tmp_part_funcs_decode.part_func_on_sa = sum_on_sa_decode;
          if uses_contra_model {
            tmp_part_funcs_decode.part_func_on_sa_4_ml = sum_on_sa_4_ml_decode;
            logsumexp(&mut tmp_sum_decode, sum_on_sa_4_ml_decode);
          } else {
            logsumexp(&mut tmp_sum_decode, sum_on_sa_decode);
          }
          tmp_part_funcs_decode.part_func_on_mls = tmp_sum_decode;
          if !is_empty(&tmp_part_funcs_decode) {
            tmp_part_func_set_mat_decode.insert(pos_pair_2, tmp_part_funcs_decode);
          }
          let term = sum_4_ml_decode + ba_score_4_ml;
          logsumexp(&mut sum_4_ml, term);
          tmp_part_funcs.part_func_4_ml = sum_4_ml;
          logsumexp(&mut tmp_sum, sum_4_ml);
          let term = sum_4_first_bpas_on_mls_decode + ba_score_4_ml;
          logsumexp(&mut sum_4_first_bpas_on_mls, term);
          tmp_part_funcs.part_func_4_first_bpas_on_mls = sum_4_first_bpas_on_mls;
          logsumexp(&mut tmp_sum, sum_4_first_bpas_on_mls);
          tmp_part_funcs.part_func_4_bpas_on_mls = tmp_sum;
          let sum_on_sa = sum_on_sa_decode + ba_score;
          tmp_part_funcs.part_func_on_sa = sum_on_sa;
          if uses_contra_model {
            let sum_on_sa_4_ml = sum_on_sa_4_ml_decode + ba_score_4_ml;
            tmp_part_funcs.part_func_on_sa_4_ml = sum_on_sa_4_ml;
            logsumexp(&mut tmp_sum, sum_on_sa_4_ml);
          } else {
            logsumexp(&mut tmp_sum, sum_on_sa);
          }
          tmp_part_funcs.part_func_on_mls = tmp_sum;
          if !is_empty(&tmp_part_funcs) {
            tmp_part_func_set_mat.insert(pos_pair, tmp_part_funcs);
          }
        }, None => {},
      }
    }
  }
  let mut final_sum_on_sa = NEG_INFINITY;
  let mut final_sum_4_ml = final_sum_on_sa;
  if is_forward {
    let pos_pair_2 = rightmost_pos_pair;
    let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
    match tmp_part_func_set_mat.get(&pos_pair_2) {
      Some(part_funcs) => {
        let part_func = part_funcs.part_func_4_ml;
        let score = part_func + align_feature_score_sets.match_2_match_count;
        logsumexp(&mut final_sum_4_ml, score);
        let part_func = part_funcs.part_func_on_sa;
        let score = part_func + align_feature_score_sets.match_2_match_count;
        logsumexp(&mut final_sum_on_sa, score);
      }, None => {},
    }
    match matchable_pos_sets_1.get(&pos_pair_2.0) {
      Some(matchable_poss) => {
        for &x in matchable_poss {
          if x >= pos_pair_2.1 {continue;}
          match tmp_part_func_set_mat.get(&(pos_pair_2.0, x)) {
            Some(part_funcs) => {
              let part_func = part_funcs.part_func_4_ml;
              let long_x = x.to_usize().unwrap();
              let insert_score_range = sta_fe_params.insert_scores_range_2[long_x + 1][long_pos_pair_2.1] + align_feature_score_sets.match_2_insert_count;
              let insert_score_range_4_ml = sta_fe_params.insert_scores_range_4_ml_2[long_x + 1][long_pos_pair_2.1] + align_feature_score_sets.match_2_insert_count;
              let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
              logsumexp(&mut final_sum_4_ml, score);
              let part_func = part_funcs.part_func_on_sa;
              let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range;
              logsumexp(&mut final_sum_on_sa, score);
            }, None => {},
          }
        }
      }, None => {},
    }
    match matchable_pos_sets_2.get(&pos_pair_2.1) {
      Some(matchable_poss) => {
        for &x in matchable_poss {
          if x >= pos_pair_2.0 {continue;}
          match tmp_part_func_set_mat.get(&(x, pos_pair_2.1)) {
            Some(part_funcs) => {
              let part_func = part_funcs.part_func_4_ml;
              let long_x = x.to_usize().unwrap();
              let insert_score_range = sta_fe_params.insert_scores_range[long_x + 1][long_pos_pair_2.0] + align_feature_score_sets.match_2_insert_count;
              let insert_score_range_4_ml = sta_fe_params.insert_scores_range_4_ml[long_x + 1][long_pos_pair_2.0] + align_feature_score_sets.match_2_insert_count;
              let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range_4_ml;
              logsumexp(&mut final_sum_4_ml, score);
              let part_func = part_funcs.part_func_on_sa;
              let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range;
              logsumexp(&mut final_sum_on_sa, score);
            }, None => {},
          }
        }
      }, None => {},
    }
  }
  (final_sum_on_sa, final_sum_4_ml)
}

pub fn get_part_func_mats_2loop<T>(
  pos_quadruple: &PosQuadruple<T>,
  sta_part_func_mats: &StaPartFuncMats<T>,
  is_forward: bool,
  pos_pair_mat_set: &PosPairMatSet<T>,
  ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair<T>,
  sta_fe_params: &StaFeParams<T>,
  matchable_pos_sets_1: &SparsePosSets<T>,
  matchable_pos_sets_2: &SparsePosSets<T>,
  align_feature_score_sets: &AlignFeatureCountSets,
) -> (SparsePartFuncMat<T>, SparsePartFuncMat<T>)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let &(i, j, k, l) = pos_quadruple;
  let leftmost_pos_pair = if is_forward {
    (i, k)
  } else {
    (i + T::one(), k + T::one())
  };
  let rightmost_pos_pair = if is_forward {
    (j - T::one(), l - T::one())
  } else {
    (j, l)
  };
  let tmp_part_func_set_mats_with_pos_pairs_decode = if is_forward {
    &sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs_decode
  } else {
    &sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs_decode
  };
  let ref tmp_part_func_set_mat_decode = tmp_part_func_set_mats_with_pos_pairs_decode[&if is_forward {leftmost_pos_pair} else {rightmost_pos_pair}];
  let iter: Poss<T> = if is_forward {
    range(i, j).collect()
  } else {
    range_inclusive(i + T::one(), j).rev().collect()
  };
  let iter_2: Poss<T> = if is_forward {
    range(k, l).collect()
  } else {
    range_inclusive(k + T::one(), l).rev().collect()
  };
  let mut part_func_mat_4_2loop = SparsePartFuncMat::<T>::default();
  let mut part_func_mat_4_2loop_decode = part_func_mat_4_2loop.clone();
  for &u in iter.iter() {
    for &v in iter_2.iter() {
      let pos_pair = (u, v);
      if (is_forward && u == i && v == k) || (!is_forward && u == j && v == l) {
        continue;
      }
      let mut sum_4_2loop = NEG_INFINITY;
      // For alignments.
      match pos_pair_mat_set.get(&pos_pair) {
        Some(pos_pair_mat) => {
          for &(m, n) in pos_pair_mat {
            if is_forward {
              if !(i < m && k < n) {continue;}
            } else {
              if !(m < j && n < l) {continue;}
            }
            let pos_pair_2 = if is_forward {
              (m - T::one(), n - T::one())
            } else {
              (m + T::one(), n + T::one())
            };
            let pos_quadruple_2 = if is_forward {
              (m, u, n, v)
            } else {
              (u, m, v, n)
            };
            if pos_quadruple_2.0 - i - T::one() + j - pos_quadruple_2.1 - T::one() > T::from_usize(MAX_2_LOOP_LEN).unwrap() {
              continue;
            }
            if pos_quadruple_2.2 - k - T::one() + l - pos_quadruple_2.3 - T::one() > T::from_usize(MAX_2_LOOP_LEN).unwrap() {
              continue;
            }
            match sta_part_func_mats
              .part_func_4d_mat_4_bpas
              .get(&pos_quadruple_2)
            {
              Some(&part_func) => {
                match tmp_part_func_set_mat_decode.get(&pos_pair_2) {
                  Some(part_funcs) => {
                    let part_func_2 = part_funcs.part_func_on_sa;
                    let twoloop_score = ss_free_energy_mat_set_pair.0.twoloop_fe_4d_mat[&(i, j, pos_quadruple_2.0, pos_quadruple_2.1)];
                    let twoloop_score_2 = ss_free_energy_mat_set_pair.1.twoloop_fe_4d_mat[&(k, l, pos_quadruple_2.2, pos_quadruple_2.3)];
                    let score = part_func_2 + part_func + twoloop_score + twoloop_score_2;
                    logsumexp(&mut sum_4_2loop, score);
                  },
                  None => {},
                }
              },
              None => {},
            }
          }
        }, None => {},
      }
      let pos_pair_2 = if is_forward {
        (u - T::one(), v - T::one())
      } else {
        (u + T::one(), v + T::one())
      };
      let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
      match sta_fe_params.ba_score_mat.get(&pos_pair) {
        Some(&ba_score) => {
          let mut sum_4_2loop_decode = NEG_INFINITY;
          match part_func_mat_4_2loop.get(&pos_pair_2) {
            Some(&part_func) => {
              let score = part_func + align_feature_score_sets.match_2_match_count;
              logsumexp(&mut sum_4_2loop_decode, score);
            }, None => {},
          }
          match matchable_pos_sets_1.get(&pos_pair_2.0) {
            Some(matchable_poss) => {
              for &x in matchable_poss {
                if is_forward && x >= pos_pair_2.1 || (!is_forward && x <= pos_pair_2.1) {continue;}
                match part_func_mat_4_2loop.get(&(pos_pair_2.0, x)) {
                  Some(&part_func) => {
                    let long_x = x.to_usize().unwrap();
                    let insert_score_range = if is_forward {
                      sta_fe_params.insert_scores_range_2[long_x + 1][long_pos_pair_2.1]
                    } else {
                      sta_fe_params.insert_scores_range_2[long_pos_pair_2.1][long_x - 1]
                    } + align_feature_score_sets.match_2_insert_count;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range;
                    logsumexp(&mut sum_4_2loop_decode, score);
                  }, None => {},
                }
              }
            }, None => {},
          }
          match matchable_pos_sets_2.get(&pos_pair_2.1) {
            Some(matchable_poss) => {
              for &x in matchable_poss {
                if is_forward && x >= pos_pair_2.0 || (!is_forward && x <= pos_pair_2.0) {continue;}
                match part_func_mat_4_2loop.get(&(x, pos_pair_2.1)) {
                  Some(&part_func) => {
                    let long_x = x.to_usize().unwrap();
                    let insert_score_range = if is_forward {
                      sta_fe_params.insert_scores_range[long_x + 1][long_pos_pair_2.0]
                    } else {
                      sta_fe_params.insert_scores_range[long_pos_pair_2.0][long_x - 1]
                    } + align_feature_score_sets.match_2_insert_count;
                    let score = part_func + align_feature_score_sets.match_2_insert_count + insert_score_range;
                    logsumexp(&mut sum_4_2loop_decode, score);
                  }, None => {},
                }
              }
            }, None => {},
          }
          if sum_4_2loop_decode > NEG_INFINITY {
            part_func_mat_4_2loop_decode.insert(pos_pair_2, sum_4_2loop_decode);
          }
          let term = sum_4_2loop_decode + ba_score;
          logsumexp(&mut sum_4_2loop, term);
          if sum_4_2loop > NEG_INFINITY {
            part_func_mat_4_2loop.insert(pos_pair, sum_4_2loop);
          }
        }, None => {},
      }
    }
  }
  (part_func_mat_4_2loop, part_func_mat_4_2loop_decode)
}

pub fn is_empty(tmp_part_funcs: &TmpPartFuncs) -> bool {
  tmp_part_funcs.part_func_on_sa == NEG_INFINITY &&
  tmp_part_funcs.part_func_4_ml == NEG_INFINITY &&
  tmp_part_funcs.part_func_4_first_bpas_on_mls == NEG_INFINITY
}

pub fn get_sta_prob_mats<T>(seq_len_pair: &PosPair<T>, sta_fe_params: &StaFeParams<T>, max_bp_span_pair: &PosPair<T>, sta_part_func_mats: &StaPartFuncMats<T>, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair<T>, produces_struct_profs: bool, global_part_func: PartFunc, uses_contra_model: bool, pos_quadruple_mat_with_len_pairs: &PosQuadrupleMatWithLenPairs<T>, produces_align_probs: bool, forward_pos_pair_mat_set: &PosPairMatSet<T>, backward_pos_pair_mat_set: &PosPairMatSet<T>, matchable_pos_sets_1: &SparsePosSets<T>, matchable_pos_sets_2: &SparsePosSets<T>, align_feature_score_sets: &AlignFeatureCountSets,) -> StaProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let mut sta_outside_part_func_4d_mat_4_bpas = PartFunc4dMat::<T>::default();
  let mut sta_prob_mats = StaProbMats::<T>::new(&seq_len_pair);
  let mut prob_coeff_mat_4_ml = PartFunc4dMat::<T>::default();
  let mut prob_coeff_mat_4_ml_2 = prob_coeff_mat_4_ml.clone();
  let leftmost_pos_pair = (T::zero(), T::zero());
  let rightmost_pos_pair = (seq_len_pair.0 - T::one(), seq_len_pair.1 - T::one());
  for substr_len_1 in range_inclusive(
    T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(),
    max_bp_span_pair.0,
  )
  .rev()
  {
    for substr_len_2 in range_inclusive(
      T::from_usize(MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL).unwrap(),
      max_bp_span_pair.1,
    )
    .rev()
    {
      match pos_quadruple_mat_with_len_pairs.get(&(substr_len_1, substr_len_2)) {
        Some(pos_pairs) => {
          for &(i, k) in pos_pairs {
            let (j, l) = (i + substr_len_1 - T::one(), k + substr_len_2 - T::one());
            let pos_quadruple = (i, j, k, l);
            match sta_part_func_mats.part_func_4d_mat_4_bpas.get(&pos_quadruple) {
              Some(&part_func_4_bpa) => {
                let (long_i, long_j, long_k, long_l) = (i.to_usize().unwrap(), j.to_usize().unwrap(), k.to_usize().unwrap(), l.to_usize().unwrap());
                let prob_coeff = part_func_4_bpa - global_part_func;
                let mut sum = NEG_INFINITY;
                let mut forward_term = sum;
                let mut backward_term = sum;
                let pos_pair_2 = (i - T::one(), k - T::one());
                match sta_part_func_mats.forward_part_func_mat_4_external_loop_decode.get(&pos_pair_2) {
                  Some(&part_func) => {
                    logsumexp(&mut forward_term, part_func);
                  }, None => {},
                }
                let pos_pair_2 = (j + T::one(), l + T::one());
                match sta_part_func_mats.backward_part_func_mat_4_external_loop_decode.get(&pos_pair_2) {
                  Some(&part_func) => {
                    logsumexp(&mut backward_term, part_func);
                  }, None => {},
                }
                let part_func_4_el = forward_term + backward_term;
                if part_func_4_el > NEG_INFINITY {
                  let coefficient = sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els[&pos_quadruple] - part_func_4_bpa;
                  sum = coefficient + part_func_4_el;
                }
                for substr_len_3 in range_inclusive(
                  substr_len_1 + T::from_usize(2).unwrap(),
                  (substr_len_1 + T::from_usize(MAX_2_LOOP_LEN + 2).unwrap()).min(max_bp_span_pair.0),
                ) {
                  for substr_len_4 in range_inclusive(
                    substr_len_2 + T::from_usize(2).unwrap(),
                    (substr_len_2 + T::from_usize(MAX_2_LOOP_LEN + 2).unwrap()).min(max_bp_span_pair.1),
                  ) {
                    match pos_quadruple_mat_with_len_pairs.get(&(substr_len_3, substr_len_4)) {
                      Some(pos_pairs_2) => {
                        for &(m, o) in pos_pairs_2 {
                          let (n, p) = (m + substr_len_3 - T::one(), o + substr_len_4 - T::one());
                          if !(m < i && j < n) || !(o < k && l < p) {continue;}
                          let (long_m, long_n, long_o, long_p ) = (m.to_usize().unwrap(), n.to_usize().unwrap(), o.to_usize().unwrap(), p.to_usize().unwrap());
                          if long_n - long_j - 1 + long_i - long_m - 1 > MAX_2_LOOP_LEN {continue;}
                          if long_p - long_l - 1 + long_k - long_o - 1 > MAX_2_LOOP_LEN {continue;}
                          let pos_quadruple_2 = (m, n, o, p);
                          match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                            Some(&part_func) => {
                              let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs_decode[&(m, o)];
                              let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs_decode[&(n, p)];
                              let mut forward_term = NEG_INFINITY;
                              let mut backward_term = forward_term;
                              let pos_pair_2 = (i - T::one(), k - T::one());
                              match forward_tmp_part_func_set_mat.get(&pos_pair_2) {
                                Some(part_funcs) => {
                                  let part_func = part_funcs.part_func_on_sa;
                                  logsumexp(&mut forward_term, part_func);
                                }, None => {},
                              }
                              let pos_pair_2 = (j + T::one(), l + T::one());
                              match backward_tmp_part_func_set_mat.get(&pos_pair_2) {
                                Some(part_funcs) => {
                                  let part_func = part_funcs.part_func_on_sa;
                                  logsumexp(&mut backward_term, part_func);
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
                                if produces_struct_profs {
                                  let bpap_4_2l = prob_coeff + part_func_4_2l;
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
                      }, None => {},
                    }
                  }
                }
                let part_func_ratio = sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple] - part_func_4_bpa;
                for (pos_pair, forward_tmp_part_func_set_mat) in &sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs_decode {
                  let &(u, v) = pos_pair;
                  if !(u < i && v < k) {continue;}
                  let pos_quadruple_2 = (u, j, v, l);
                  let mut forward_term = NEG_INFINITY;
                  let mut forward_term_2 = forward_term;
                  let pos_pair_2 = (i - T::one(), k - T::one());
                  match forward_tmp_part_func_set_mat.get(&pos_pair_2) {
                    Some(part_funcs) => {
                      let part_func = part_funcs.part_func_4_bpas_on_mls;
                      logsumexp(&mut forward_term, part_func);
                      let part_func = if !uses_contra_model {part_funcs.part_func_on_sa} else {part_funcs.part_func_on_sa_4_ml};
                      logsumexp(&mut forward_term_2, part_func);
                    }, None => {},
                  }
                  let mut part_func_4_ml = NEG_INFINITY;
                  match prob_coeff_mat_4_ml.get(&pos_quadruple_2) {
                    Some(prob_coeff_ml) => {
                      let prob_coeff_ml = prob_coeff_ml + part_func_ratio;
                      let term = prob_coeff_ml + forward_term;
                      logsumexp(&mut part_func_4_ml, term);
                    }, None => {},
                  }
                  match prob_coeff_mat_4_ml_2.get(&pos_quadruple_2) {
                    Some(prob_coeff_ml_2) => {
                      let prob_coeff_ml_2 = prob_coeff_ml_2 + part_func_ratio;
                      let term = prob_coeff_ml_2 + forward_term_2;
                      logsumexp(&mut part_func_4_ml, term);
                    }, None => {},
                  }
                  if part_func_4_ml > NEG_INFINITY {
                    logsumexp(&mut sum, part_func_4_ml);
                  }
                }
                if sum > NEG_INFINITY {
                  sta_outside_part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
                  let bpap = prob_coeff + sum;
                  if produces_align_probs {
                    sta_prob_mats.basepair_align_prob_mat.insert(pos_quadruple, bpap);
                  }
                  debug_assert!(NEG_INFINITY <= bpap && bpap <= 0.);
                  match sta_prob_mats.bpp_mat_pair.0.get_mut(&(i, j)) {
                    Some(bpp) => {
                      logsumexp(bpp, bpap);
                    }, None => {
                      sta_prob_mats.bpp_mat_pair.0.insert((i, j), bpap);
                    },
                  }
                  match sta_prob_mats.bpp_mat_pair.1.get_mut(&(k, l)) {
                    Some(bpp) => {
                      logsumexp(bpp, bpap);
                    }, None => {
                      sta_prob_mats.bpp_mat_pair.1.insert((k, l), bpap);
                    },
                  }
                  if produces_struct_profs {
                    logsumexp(&mut sta_prob_mats.bpp_mat_pair_2.0[long_i], bpap);
                    logsumexp(&mut sta_prob_mats.bpp_mat_pair_2.0[long_j], bpap);
                    logsumexp(&mut sta_prob_mats.bpp_mat_pair_2.1[long_k], bpap);
                    logsumexp(&mut sta_prob_mats.bpp_mat_pair_2.1[long_l], bpap);
                  }
                  let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
                  let multi_loop_closing_basepairing_fe = ss_free_energy_mat_set_pair.0.ml_closing_bp_fe_mat[&(i, j)];
                  let multi_loop_closing_basepairing_fe_2 = ss_free_energy_mat_set_pair.1.ml_closing_bp_fe_mat[&(k, l)];
                  let coefficient = sum + bpa_score
                    + multi_loop_closing_basepairing_fe
                    + multi_loop_closing_basepairing_fe_2;
                  let ref backward_tmp_part_func_set_mat = sta_part_func_mats
                    .backward_tmp_part_func_set_mats_with_pos_pairs_decode[&(j, l)];
                  for pos_pair in sta_fe_params.ba_score_mat.keys() {
                    let &(u, v) = pos_pair;
                    if !(i < u && u < j && k < v && v < l) {continue;}
                    let mut backward_term = NEG_INFINITY;
                    let mut backward_term_2 = backward_term;
                    let pos_pair_2 = (u + T::one(), v + T::one());
                    match backward_tmp_part_func_set_mat.get(&pos_pair_2) {
                      Some(part_funcs) => {
                        let part_func = part_funcs.part_func_on_mls;
                        logsumexp(&mut backward_term, part_func);
                        let part_func = part_funcs.part_func_4_bpas_on_mls;
                        logsumexp(&mut backward_term_2, part_func);
                      }, None => {},
                    }
                    let pos_quadruple_2 = (i, u, k, v);
                    let prob_coeff_4_ml = coefficient + backward_term;
                    match prob_coeff_mat_4_ml.get_mut(&pos_quadruple_2) {
                      Some(x) => {
                        logsumexp(x, prob_coeff_4_ml);
                      }, None => {
                        prob_coeff_mat_4_ml.insert(pos_quadruple_2, prob_coeff_4_ml);
                      },
                    }
                    let prob_coeff_4_ml_2 = coefficient + backward_term_2;
                    match prob_coeff_mat_4_ml_2.get_mut(&pos_quadruple_2) {
                      Some(x) => {
                        logsumexp(x, prob_coeff_4_ml_2);
                      }, None => {
                        prob_coeff_mat_4_ml_2.insert(pos_quadruple_2, prob_coeff_4_ml_2);
                      },
                    }
                  }
                }
              }, None => {},
            }
          }
        }, None => {},
      }
    }
  }
  for bpp in sta_prob_mats.bpp_mat_pair.0.values_mut() {
    *bpp = expf(*bpp);
  }
  for bpp in sta_prob_mats.bpp_mat_pair.1.values_mut() {
    *bpp = expf(*bpp);
  }
  if produces_struct_profs || produces_align_probs {
    let mut upp_mat_pair_4_el_range = (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    let mut upp_mat_pair_4_hl_range = (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    for u in range(T::zero(), seq_len_pair.0 - T::one()) {
      let long_u = u.to_usize().unwrap();
      for v in range(T::zero(), seq_len_pair.1 - T::one()) {
        let pos_pair = (u, v);
        let long_v = v.to_usize().unwrap();
        let pos_pair_2 = (u + T::one(), v + T::one());
        let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
        if u > T::zero() && v > T::zero() {
          let pos_pair_4_ba = (u - T::one(), v - T::one());
          match sta_fe_params.ba_score_mat.get(&pos_pair) {
            Some(&ba_score) => {
              let mut forward_term = NEG_INFINITY;
              match sta_part_func_mats.forward_part_func_mat_4_external_loop_decode.get(&pos_pair_4_ba) {
                Some(&part_func) => {
                  logsumexp(&mut forward_term, part_func);
                }, None => {},
              }
              let mut backward_term = NEG_INFINITY;
              match sta_part_func_mats.backward_part_func_mat_4_external_loop_decode.get(&pos_pair_2) {
                Some(&part_func) => {
                  logsumexp(&mut backward_term, part_func);
                }, None => {},
              }
              let ba_score = ba_score + if !uses_contra_model {0.} else {2. * CONTRA_EL_UNPAIRED_FE};
              let bap_4_el = ba_score + forward_term + backward_term - global_part_func;
              if produces_struct_profs {
                logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.0[long_u], bap_4_el);
                logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.1[long_v], bap_4_el);
              }
              if produces_align_probs {
                match sta_prob_mats.loop_align_prob_mat.get_mut(&pos_pair) {
                  Some(loop_align_prob) => {
                    logsumexp(loop_align_prob, bap_4_el);
                  }
                  None => {
                    sta_prob_mats.loop_align_prob_mat.insert(pos_pair, bap_4_el);
                  }
                }
              }
            }, None => {},
          }
        }
        if produces_struct_profs {
          match sta_part_func_mats.forward_part_func_mat_4_external_loop.get(&pos_pair) {
            Some(&part_func) => {
              let is_begin = pos_pair == leftmost_pos_pair;
              let forward_term = part_func + if is_begin {align_feature_score_sets.init_insert_count} else {align_feature_score_sets.match_2_insert_count};
              match matchable_pos_sets_1.get(&pos_pair_2.0) {
                Some(matchable_poss) => {
                  for &x in matchable_poss {
                    if x <= pos_pair_2.1 {continue;}
                    let pos_pair_3 = (pos_pair_2.0, x);
                    match sta_part_func_mats.backward_part_func_mat_4_external_loop.get(&pos_pair_3) {
                      Some(&backward_term) => {
                        let long_x = x.to_usize().unwrap();
                        let is_end = pos_pair_3 == rightmost_pos_pair;
                        let insert_score_range = sta_fe_params.insert_scores_range_4_el_2[long_pos_pair_2.1][long_x - 1] + if is_end {0.} else {align_feature_score_sets.match_2_insert_count};
                        let upp_4_el = forward_term + insert_score_range + backward_term - global_part_func;
                        let pos_pair_4 = (pos_pair_2.1, x - T::one());
                        match upp_mat_pair_4_el_range.1.get_mut(&pos_pair_4) {
                          Some(upp) => {
                            logsumexp(upp, upp_4_el);
                          }, None => {
                            upp_mat_pair_4_el_range.1.insert(pos_pair_4, upp_4_el);
                          },
                        }
                      }, None => {},
                    }
                  }
                }, None => {},
              }
              match matchable_pos_sets_2.get(&pos_pair_2.1) {
                Some(matchable_poss) => {
                  for &x in matchable_poss {
                    if x <= pos_pair_2.0 {continue;}
                    let pos_pair_3 = (x, pos_pair_2.1);
                    match sta_part_func_mats.backward_part_func_mat_4_external_loop.get(&pos_pair_3) {
                      Some(&backward_term) => {
                        let long_x = x.to_usize().unwrap();
                        let is_end = pos_pair_3 == rightmost_pos_pair;
                        let insert_score_range = sta_fe_params.insert_scores_range_4_el[long_pos_pair_2.0][long_x - 1] + if is_end {0.} else {align_feature_score_sets.match_2_insert_count};
                        let upp_4_el = forward_term + insert_score_range + backward_term - global_part_func;
                        let pos_pair_4 = (pos_pair_2.0, x - T::one());
                        match upp_mat_pair_4_el_range.0.get_mut(&pos_pair_4) {
                          Some(upp) => {
                            logsumexp(upp, upp_4_el);
                          }, None => {
                            upp_mat_pair_4_el_range.0.insert(pos_pair_4, upp_4_el);
                          },
                        }
                      }, None => {},
                    }
                  }
                }, None => {},
              }
            }, None => {},
          }
        }
      }
    }
    for (pos_quadruple, &part_func_4_bpa) in &sta_outside_part_func_4d_mat_4_bpas {
      let (i, j, k, l) = *pos_quadruple;
      let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
      let prob_coeff = part_func_4_bpa - global_part_func + bpa_score;
      let hl_fe = match ss_free_energy_mat_set_pair.0.hl_fe_mat.get(&(i, j)) {
        Some(&hl_fe) => {
          hl_fe
        }, None => {
          NEG_INFINITY
        },
      };
      let hl_fe_2 = match ss_free_energy_mat_set_pair.1.hl_fe_mat.get(&(k, l)) {
        Some(&hl_fe_2) => {
          hl_fe_2
        }, None => {
          NEG_INFINITY
        },
      };
      let prob_coeff_4_hl = prob_coeff + hl_fe + hl_fe_2;
      let multi_loop_closing_basepairing_score = ss_free_energy_mat_set_pair.0.ml_closing_bp_fe_mat[&(i, j)];
      let multi_loop_closing_basepairing_score_2 = ss_free_energy_mat_set_pair.1.ml_closing_bp_fe_mat[&(k, l)];
      let prob_coeff_4_ml = prob_coeff
        + multi_loop_closing_basepairing_score
        + multi_loop_closing_basepairing_score_2;
      let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs[&(i, k)];
      let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs[&(j, l)];
      let ref forward_tmp_part_func_set_mat_decode = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_pairs_decode[&(i, k)];
      let ref backward_tmp_part_func_set_mat_decode = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_pairs_decode[&(j, l)];
      let (_, forward_part_func_mat_4_2loop_decode) =
        if produces_align_probs {
          get_part_func_mats_2loop(&pos_quadruple, &sta_part_func_mats, true, forward_pos_pair_mat_set, ss_free_energy_mat_set_pair, sta_fe_params, matchable_pos_sets_1, matchable_pos_sets_2, align_feature_score_sets)
        } else {
          (SparsePartFuncMat::<T>::default(), SparsePartFuncMat::<T>::default())
        };
      let (_, backward_part_func_mat_4_2loop_decode) =
        if produces_align_probs {
          get_part_func_mats_2loop(&pos_quadruple, &sta_part_func_mats, false, backward_pos_pair_mat_set, ss_free_energy_mat_set_pair, sta_fe_params, matchable_pos_sets_1, matchable_pos_sets_2, align_feature_score_sets)
        } else {
          (SparsePartFuncMat::<T>::default(), SparsePartFuncMat::<T>::default())
        };
      for u in range(i, j) {
        let long_u = u.to_usize().unwrap();
        for v in range(k, l) {
          let pos_pair = (u, v);
          let long_v = v.to_usize().unwrap();
          let pos_pair_2 = (u + T::one(), v + T::one());
          let long_pos_pair_2 = (pos_pair_2.0.to_usize().unwrap(), pos_pair_2.1.to_usize().unwrap());
          match sta_fe_params.ba_score_mat.get(&pos_pair) {
            Some(&ba_score) => {
              let mut backward_term_4_align_on_sa = NEG_INFINITY;
              let mut backward_term_4_align_4_ml = backward_term_4_align_on_sa;
              let mut backward_term_4_align_4_bpas_on_mls = backward_term_4_align_on_sa;
              let mut backward_term_4_align_on_mls = backward_term_4_align_on_sa;
              let mut backward_term_4_align_4_2loop = backward_term_4_align_on_sa;
              match backward_tmp_part_func_set_mat_decode.get(&pos_pair_2) {
                Some(part_funcs) => {
                  let part_func = part_funcs.part_func_on_sa;
                  logsumexp(&mut backward_term_4_align_on_sa, part_func);
                  if produces_align_probs {
                    let part_func = part_funcs.part_func_4_ml;
                    logsumexp(&mut backward_term_4_align_4_ml, part_func);
                    let part_func = part_funcs.part_func_4_bpas_on_mls;
                    logsumexp(&mut backward_term_4_align_4_bpas_on_mls, part_func);
                    let part_func = part_funcs.part_func_on_mls;
                    logsumexp(&mut backward_term_4_align_on_mls, part_func);
                  }
                }, None => {},
              }
              if produces_align_probs {
                match backward_part_func_mat_4_2loop_decode.get(&pos_pair_2) {
                  Some(&part_func_4_2loop) => {
                    logsumexp(&mut backward_term_4_align_4_2loop, part_func_4_2loop);
                  }, None => {},
                }
              }
              let pos_pair_4_ba = (u - T::one(), v - T::one());
              let mut bap_4_hl = NEG_INFINITY;
              let mut bap_4_ml = bap_4_hl;
              let mut bap_4_2l = bap_4_hl;
              let loop_align_score_ml = ba_score + if !uses_contra_model {0.} else {2. * CONTRA_ML_UNPAIRED_FE};
              match forward_tmp_part_func_set_mat_decode.get(&pos_pair_4_ba) {
                Some(part_funcs) => {
                  let part_func = part_funcs.part_func_on_sa;
                  let term = prob_coeff_4_hl + ba_score + part_func + backward_term_4_align_on_sa;
                  logsumexp(&mut bap_4_hl, term);
                  if produces_align_probs {
                    let part_func = if !uses_contra_model {part_funcs.part_func_on_sa} else {part_funcs.part_func_on_sa_4_ml};
                    let term = prob_coeff_4_ml + loop_align_score_ml + part_func + backward_term_4_align_4_ml;
                    logsumexp(&mut bap_4_ml, term);
                    let part_func = part_funcs.part_func_4_first_bpas_on_mls;
                    let term = prob_coeff_4_ml
                      + loop_align_score_ml
                      + part_func
                      + backward_term_4_align_4_bpas_on_mls;
                    logsumexp(&mut bap_4_ml, term);
                    let part_func = part_funcs.part_func_4_ml;
                    let term = prob_coeff_4_ml + loop_align_score_ml + part_func + backward_term_4_align_on_mls;
                    logsumexp(&mut bap_4_ml, term);
                    let part_func = part_funcs.part_func_on_sa;
                    let term = prob_coeff + ba_score + part_func + backward_term_4_align_4_2loop;
                    logsumexp(&mut bap_4_2l, term);
                  }
                }, None => {},
              }
              if produces_struct_profs {
                logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.0[long_u], bap_4_hl);
                logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.1[long_v], bap_4_hl);
              }
              if produces_align_probs {
                match forward_part_func_mat_4_2loop_decode.get(&pos_pair_4_ba) {
                  Some(&part_func_4_2loop) => {
                    let term = prob_coeff + ba_score + part_func_4_2loop + backward_term_4_align_on_sa;
                    logsumexp(&mut bap_4_2l, term);
                  }, None => {},
                }
                let mut prob = NEG_INFINITY;
                logsumexp(&mut prob, bap_4_hl);
                logsumexp(&mut prob, bap_4_ml);
                logsumexp(&mut prob, bap_4_2l);
                match sta_prob_mats.loop_align_prob_mat.get_mut(&pos_pair) {
                  Some(loop_align_prob) => {
                    logsumexp(loop_align_prob, prob);
                  }
                  None => {
                    sta_prob_mats.loop_align_prob_mat.insert(pos_pair, prob);
                  }
                }
              }
            }, None => {},
          }
          if produces_struct_profs {
            match forward_tmp_part_func_set_mat.get(&pos_pair) {
              Some(part_funcs) => {
                let part_func = part_funcs.part_func_on_sa;
                let forward_term = part_func + align_feature_score_sets.match_2_insert_count;
                match matchable_pos_sets_1.get(&pos_pair_2.0) {
                  Some(matchable_poss) => {
                    for &x in matchable_poss {
                      if x <= pos_pair_2.1 {continue;}
                      let pos_pair_3 = (pos_pair_2.0, x);
                      match backward_tmp_part_func_set_mat.get(&pos_pair_3) {
                        Some(part_funcs) => {
                          let backward_term = part_funcs.part_func_on_sa;
                          let long_x = x.to_usize().unwrap();
                          let insert_score_range = sta_fe_params.insert_scores_range_4_el_2[long_pos_pair_2.1][long_x - 1] + align_feature_score_sets.match_2_insert_count;
                          let upp_4_hl = prob_coeff_4_hl + forward_term + insert_score_range + backward_term;
                          let pos_pair_4 = (pos_pair_2.1, x - T::one());
                          match upp_mat_pair_4_hl_range.1.get_mut(&pos_pair_4) {
                            Some(upp) => {
                              logsumexp(upp, upp_4_hl);
                            }, None => {
                              upp_mat_pair_4_hl_range.1.insert(pos_pair_4, upp_4_hl);
                            },
                          }
                        }, None => {},
                      }
                    }
                  }, None => {},
                }
                match matchable_pos_sets_2.get(&pos_pair_2.1) {
                  Some(matchable_poss) => {
                    for &x in matchable_poss {
                      if x <= pos_pair_2.0 {continue;}
                      let pos_pair_3 = (x, pos_pair_2.1);
                      match backward_tmp_part_func_set_mat.get(&pos_pair_3) {
                        Some(part_funcs) => {
                          let backward_term = part_funcs.part_func_on_sa;
                          let long_x = x.to_usize().unwrap();
                          let insert_score_range = sta_fe_params.insert_scores_range_4_el[long_pos_pair_2.0][long_x - 1] + align_feature_score_sets.match_2_insert_count;
                          let upp_4_hl = prob_coeff_4_hl + forward_term + insert_score_range + backward_term;
                          let pos_pair_4 = (pos_pair_2.0, x - T::one());
                          match upp_mat_pair_4_hl_range.0.get_mut(&pos_pair_4) {
                            Some(upp) => {
                              logsumexp(upp, upp_4_hl);
                            }, None => {
                              upp_mat_pair_4_hl_range.0.insert(pos_pair_4, upp_4_hl);
                            },
                          }
                        }, None => {},
                      }
                    }
                  }, None => {},
                }
              }, None => {},
            }
          }
        }
      }
    }
    if produces_struct_profs {
      for (pos_pair, &upp_4_el) in &upp_mat_pair_4_el_range.0 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.0[long_i], upp_4_el);
        }
      }
      for (pos_pair, &upp_4_el) in &upp_mat_pair_4_el_range.1 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(&mut sta_prob_mats.upp_mat_pair_4_el.1[long_i], upp_4_el);
        }
      }
      for (pos_pair, &upp_4_hl) in &upp_mat_pair_4_hl_range.0 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.0[long_i], upp_4_hl);
        }
      }
      for (pos_pair, &upp_4_hl) in &upp_mat_pair_4_hl_range.1 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(&mut sta_prob_mats.upp_mat_pair_4_hl.1[long_i], upp_4_hl);
        }
      }
      for (i, upp) in sta_prob_mats.upp_mat_pair_4_ml.0.iter_mut().enumerate() {
        let mut sum = sta_prob_mats.upp_mat_pair_4_hl.0[i];
        logsumexp(&mut sum, sta_prob_mats.upp_mat_pair_4_2l.0[i]);
        logsumexp(&mut sum, sta_prob_mats.upp_mat_pair_4_el.0[i]);
        logsumexp(&mut sum, sta_prob_mats.bpp_mat_pair_2.0[i]);
        *upp = 1. - expf(sum);
      }
      for (i, upp) in sta_prob_mats.upp_mat_pair_4_ml.1.iter_mut().enumerate() {
        let mut sum = sta_prob_mats.upp_mat_pair_4_hl.1[i];
        logsumexp(&mut sum, sta_prob_mats.upp_mat_pair_4_2l.1[i]);
        logsumexp(&mut sum, sta_prob_mats.upp_mat_pair_4_el.1[i]);
        logsumexp(&mut sum, sta_prob_mats.bpp_mat_pair_2.1[i]);
        *upp = 1. - expf(sum);
      }
      for upp in sta_prob_mats.upp_mat_pair_4_hl.0.iter_mut() {
        *upp = expf(*upp);
      }
      for upp in sta_prob_mats.upp_mat_pair_4_hl.1.iter_mut() {
        *upp = expf(*upp);
      }
      for upp in sta_prob_mats.upp_mat_pair_4_2l.0.iter_mut() {
        *upp = expf(*upp);
      }
      for upp in sta_prob_mats.upp_mat_pair_4_2l.1.iter_mut() {
        *upp = expf(*upp);
      }
      for upp in sta_prob_mats.upp_mat_pair_4_el.0.iter_mut() {
        *upp = expf(*upp);
      }
      for upp in sta_prob_mats.upp_mat_pair_4_el.1.iter_mut() {
        *upp = expf(*upp);
      }
      for bpp in sta_prob_mats.bpp_mat_pair_2.0.iter_mut() {
        *bpp = expf(*bpp);
      }
      for bpp in sta_prob_mats.bpp_mat_pair_2.1.iter_mut() {
        *bpp = expf(*bpp);
      }
    }
    if produces_align_probs {
      for loop_align_prob in sta_prob_mats.loop_align_prob_mat.values_mut() {
        *loop_align_prob = expf(*loop_align_prob);
      }
      for basepair_align_prob in sta_prob_mats.basepair_align_prob_mat.values_mut() {
        *basepair_align_prob = expf(*basepair_align_prob);
      }
    }
  }
  sta_prob_mats
}

pub fn pct_of_prob_mats<T>(prob_mats_with_rna_id_pairs: &StaProbMatsWithRnaIdPairs<T>, rna_id: RnaId, num_of_rnas: usize, upp_mat_len: usize, produces_struct_profs: bool) -> PctStaProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let weight = 1. / (num_of_rnas - 1) as Prob;
  let mut pct_prob_mats = PctStaProbMats::new(upp_mat_len);
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
    let ref_2_bpp_mat_2 = if rna_id < rna_id_2 {
      &ref_2_prob_mats.bpp_mat_pair_2.0
    } else {
      &ref_2_prob_mats.bpp_mat_pair_2.1
    };
    for (i, &bpp) in ref_2_bpp_mat_2.iter().enumerate() {
      let weighted_bpp = weight * bpp;
      pct_prob_mats.bpp_mat_2[i] += weighted_bpp;
    }
    if produces_struct_profs {
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
  pct_prob_mats
}

pub fn get_max_bp_span<T>(sparse_bpp_mat: &SparseProbMat<T>) -> T
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let max_bp_span = sparse_bpp_mat.iter().map(|(pos_pair, _)| {pos_pair.1 - pos_pair.0 + T::one()}).max();
  match max_bp_span {
    Some(max_bp_span) => {max_bp_span},
    None => {T::zero()},
  }
}

pub fn sparsify_bpp_mat<T>(sparse_bpp_mat: &SparseProbMat<T>, min_bpp: Prob) -> SparseProbMat<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  sparse_bpp_mat.iter().filter(|(_, &bpp)| {bpp >= min_bpp}).map(|(&(i, j), &bpp)| {((i + T::one(), j + T::one()), bpp)}).collect()
}

pub fn sparsify_align_prob_mat<T>(align_prob_mat: &ProbMat, min_align_prob: Prob) -> SparseProbMat<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let mut sparse_align_prob_mat = SparseProbMat::<T>::default();
  for (i, align_probs) in align_prob_mat.iter().enumerate() {
    let i = T::from_usize(i).unwrap();
    for (j, &align_prob) in align_probs.iter().enumerate() {
      if align_prob >= min_align_prob {
        let j = T::from_usize(j).unwrap();
        sparse_align_prob_mat.insert((i, j), align_prob);
      }
    }
  }
  sparse_align_prob_mat
}

pub fn consprob<T>(thread_pool: &mut Pool, fasta_records: &FastaRecords, min_bpp: Prob, min_align_prob: Prob, produces_struct_profs: bool, uses_contra_model: bool, produces_align_probs: bool, align_feature_score_sets: &AlignFeatureCountSets,) -> (ProbMatSets<T>, AlignProbMatSetsWithRnaIdPairs<T>)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let num_of_fasta_records = fasta_records.len();
  let mut sparse_bpp_mats = vec![SparseProbMat::<T>::new(); num_of_fasta_records];
  let mut max_bp_spans = vec![T::zero(); num_of_fasta_records];
  let mut ss_free_energy_mat_sets = vec![SsFreeEnergyMats::<T>::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (sparse_bpp_mat, max_bp_span, fasta_record, ss_free_energy_mats) in multizip((sparse_bpp_mats.iter_mut(), max_bp_spans.iter_mut(), fasta_records.iter(), ss_free_energy_mat_sets.iter_mut())) {
      let seq_len = fasta_record.seq.len();
      scope.execute(move || {
        let (bpp_mat, obtained_ss_free_energy_mats) = mccaskill_algo(&fasta_record.seq[1 .. seq_len - 1], false, false, &StructFeatureCountSets::new(0.));
        *ss_free_energy_mats = sparsify(&obtained_ss_free_energy_mats, &bpp_mat, min_bpp);
        *sparse_bpp_mat = sparsify_bpp_mat::<T>(&bpp_mat, min_bpp);
        *max_bp_span = get_max_bp_span::<T>(sparse_bpp_mat);
      });
    }
  });
  let mut prob_mats_with_rna_id_pairs = StaProbMatsWithRnaIdPairs::<T>::default();
  let mut align_prob_mats_with_rna_id_pairs = SparseProbMatsWithRnaIdPairs::<T>::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      prob_mats_with_rna_id_pairs.insert(rna_id_pair, StaProbMats::<T>::origin());
      align_prob_mats_with_rna_id_pairs.insert(rna_id_pair, SparseProbMat::<T>::default());
    }
  }
  thread_pool.scoped(|scope| {
    for (rna_id_pair, align_prob_mat) in align_prob_mats_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      scope.execute(move || {
        *align_prob_mat = sparsify_align_prob_mat(&durbin_algo(&seq_pair, align_feature_score_sets,), min_align_prob);
      });
    }
  });
  thread_pool.scoped(|scope| {
    for (rna_id_pair, prob_mats) in prob_mats_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      let max_bp_span_pair = (max_bp_spans[rna_id_pair.0], max_bp_spans[rna_id_pair.1]);
      let bpp_mat_pair = (&sparse_bpp_mats[rna_id_pair.0], &sparse_bpp_mats[rna_id_pair.1]);
      let ss_free_energy_mat_set_pair = (&ss_free_energy_mat_sets[rna_id_pair.0], &ss_free_energy_mat_sets[rna_id_pair.1]);
      let seq_len_pair = (
        T::from_usize(seq_pair.0.len()).unwrap(),
        T::from_usize(seq_pair.1.len()).unwrap(),
      );
      let ref align_prob_mat = align_prob_mats_with_rna_id_pairs[rna_id_pair];
      let (forward_pos_pair_mat_set, backward_pos_pair_mat_set, pos_quadruple_mat, pos_quadruple_mat_with_len_pairs, matchable_pos_sets_1, matchable_pos_sets_2) = get_sparse_pos_sets(&bpp_mat_pair, align_prob_mat, &seq_len_pair);
      scope.execute(move || {
        let sta_fe_params = StaFeParams::<T>::new(&seq_pair, &seq_len_pair, &pos_quadruple_mat, align_prob_mat, uses_contra_model, align_feature_score_sets);
        *prob_mats = io_algo_4_prob_mats::<T>(&seq_len_pair, &sta_fe_params, &max_bp_span_pair, &ss_free_energy_mat_set_pair, produces_struct_profs, &forward_pos_pair_mat_set, &backward_pos_pair_mat_set, uses_contra_model, &pos_quadruple_mat_with_len_pairs, produces_align_probs, &matchable_pos_sets_1, &matchable_pos_sets_2, align_feature_score_sets,);
      });
    }
  });
  let mut prob_mat_sets = vec![PctStaProbMats::<T>::origin(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (rna_id, prob_mats) in prob_mat_sets.iter_mut().enumerate() {
      let ref ref_2_prob_mats_with_rna_id_pairs = prob_mats_with_rna_id_pairs;
      let seq_len = fasta_records[rna_id].seq.len();
      scope.execute(move || {
        *prob_mats = pct_of_prob_mats::<T>(ref_2_prob_mats_with_rna_id_pairs, rna_id, num_of_fasta_records, seq_len, produces_struct_profs);
      });
    }
  });
  let mut align_prob_mat_sets_with_rna_id_pairs = AlignProbMatSetsWithRnaIdPairs::<T>::default();
  if produces_align_probs {
    for rna_id_1 in 0 .. num_of_fasta_records {
      for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
        let rna_id_pair = (rna_id_1, rna_id_2);
        let ref prob_mats = prob_mats_with_rna_id_pairs[&rna_id_pair];
        let mut align_prob_mats = AlignProbMats::<T>::new();
        align_prob_mats.loop_align_prob_mat = prob_mats.loop_align_prob_mat.clone();
        align_prob_mats.basepair_align_prob_mat = prob_mats.basepair_align_prob_mat.clone();
        align_prob_mats.align_prob_mat = align_prob_mats.loop_align_prob_mat.clone();
        for (pos_quadruple, &bpap) in &align_prob_mats.basepair_align_prob_mat {
          let pos_pair = (pos_quadruple.0, pos_quadruple.2);
          match align_prob_mats.align_prob_mat.get_mut(&pos_pair) {
            Some(bap) => {
              *bap += bpap;
            }, None => {
              align_prob_mats.align_prob_mat.insert(pos_pair, bpap);
            }
          }
          let pos_pair = (pos_quadruple.1, pos_quadruple.3);
          match align_prob_mats.align_prob_mat.get_mut(&pos_pair) {
            Some(bap) => {
              *bap += bpap;
            }, None => {
              align_prob_mats.align_prob_mat.insert(pos_pair, bpap);
            }
          }
        }
        align_prob_mat_sets_with_rna_id_pairs.insert(rna_id_pair, align_prob_mats);
      }
    }
  }
  (prob_mat_sets, align_prob_mat_sets_with_rna_id_pairs)
}

pub fn sparsify<T>(ss_free_energy_mats: &SsFreeEnergyMats<T>, bpp_mat: &SparseProbMat<T>, min_bpp: Prob) -> SsFreeEnergyMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let mut new_ss_free_energy_mats = SsFreeEnergyMats::new();
  new_ss_free_energy_mats.hl_fe_mat = ss_free_energy_mats.hl_fe_mat.iter().filter(|(pos_pair, _)| {bpp_mat[pos_pair] >= min_bpp}).map(|(&(i, j), &free_energy)| {((i + T::one(), j + T::one()), free_energy)}).collect();
  new_ss_free_energy_mats.twoloop_fe_4d_mat = ss_free_energy_mats.twoloop_fe_4d_mat.iter().filter(|(&(i, j, k, l), _)| {bpp_mat[&(i, j)] >= min_bpp && bpp_mat[&(k, l)] >= min_bpp}).map(|(&(i, j, k, l), &free_energy)| {((i + T::one(), j + T::one(), k + T::one(), l + T::one()), free_energy)}).collect();
  new_ss_free_energy_mats.ml_closing_bp_fe_mat = ss_free_energy_mats.ml_closing_bp_fe_mat.iter().filter(|(pos_pair, _)| {bpp_mat[pos_pair] >= min_bpp}).map(|(&(i, j), &free_energy)| {((i + T::one(), j + T::one()), free_energy)}).collect();
  new_ss_free_energy_mats.accessible_bp_fe_mat = ss_free_energy_mats.accessible_bp_fe_mat.iter().filter(|(pos_pair, _)| {bpp_mat[pos_pair] >= min_bpp}).map(|(&(i, j), &free_energy)| {((i + T::one(), j + T::one()), free_energy)}).collect();
  new_ss_free_energy_mats
}

pub fn write_prob_mat_sets<T>(
  output_dir_path: &Path,
  prob_mat_sets: &ProbMatSets<T>,
  produces_struct_profs: bool,
  align_prob_mat_sets_with_rna_id_pairs: &AlignProbMatSetsWithRnaIdPairs<T>,
  produces_align_probs: bool,
) where
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
  if produces_struct_profs {
    let bpp_mat_file_path_2 = output_dir_path.join(BPP_MAT_FILE_NAME_2);
    let mut writer_2_bpp_mat_file_2 = BufWriter::new(File::create(bpp_mat_file_path_2).unwrap());
    let mut buf_4_writer_2_bpp_mat_file_2 = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let seq_len = prob_mats.bpp_mat_2.len();
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (i, &bpp) in prob_mats.bpp_mat_2.iter().enumerate() {
        if i == 0 || i == seq_len - 1 {continue;}
        buf_4_rna_id.push_str(&format!("{},{} ", i - 1, bpp));
      }
      buf_4_writer_2_bpp_mat_file_2.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_bpp_mat_file_2.write_all(buf_4_writer_2_bpp_mat_file_2.as_bytes());
    let upp_mat_file_path = output_dir_path.join(UPP_MAT_ON_HL_FILE_NAME);
    let mut writer_2_upp_mat_file = BufWriter::new(File::create(upp_mat_file_path).unwrap());
    let mut buf_4_writer_2_upp_mat_file = String::new();
    for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
      let seq_len = prob_mats.upp_mat_4_hl.len();
      let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
      for (i, &upp) in prob_mats.upp_mat_4_hl.iter().enumerate() {
        if i == 0 || i == seq_len - 1 {
          continue;
        }
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
        if i == 0 || i == seq_len - 1 {
          continue;
        }
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
        if i == 0 || i == seq_len - 1 {
          continue;
        }
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
        if i == 0 || i == seq_len - 1 {
          continue;
        }
        buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
      }
      buf_4_writer_2_upp_mat_file.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_upp_mat_file.write_all(buf_4_writer_2_upp_mat_file.as_bytes());
  }
  if produces_align_probs {
    let loop_align_prob_mat_file_path = output_dir_path.join(LOOP_ALIGN_PROB_MAT_FILE_NAME);
    let basepair_align_prob_mat_file_path = output_dir_path.join(BASEPAIR_ALIGN_PROB_MAT_FILE_NAME);
    let align_prob_mat_file_path = output_dir_path.join(ALIGN_PROB_MAT_FILE_NAME);
    let mut writer_2_loop_align_prob_mat_file = BufWriter::new(File::create(loop_align_prob_mat_file_path).unwrap());
    let mut writer_2_basepair_align_prob_mat_file = BufWriter::new(File::create(basepair_align_prob_mat_file_path).unwrap());
    let mut writer_2_align_prob_mat_file = BufWriter::new(File::create(align_prob_mat_file_path).unwrap());
    let mut buf_4_writer_2_loop_align_prob_mat_file = String::new();
    let mut buf_4_writer_2_basepair_align_prob_mat_file = String::new();
    let mut buf_4_writer_2_align_prob_mat_file = String::new();
    for (rna_id_pair, align_prob_mats) in align_prob_mat_sets_with_rna_id_pairs.iter() {
      let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
      for (&(i, j), &loop_align_prob) in align_prob_mats.loop_align_prob_mat.iter() {
        buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - T::one(), j - T::one(), loop_align_prob));
      }
      buf_4_writer_2_loop_align_prob_mat_file.push_str(&buf_4_rna_id_pair);
      let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
      for (&(i, j, k, l), &basepair_align_prob) in align_prob_mats.basepair_align_prob_mat.iter() {
        buf_4_rna_id_pair.push_str(&format!("{},{},{},{},{} ", i - T::one(), j - T::one(), k - T::one(), l - T::one(), basepair_align_prob));
      }
      buf_4_writer_2_basepair_align_prob_mat_file.push_str(&buf_4_rna_id_pair);
      let mut buf_4_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
      for (&(i, j), &align_prob) in align_prob_mats.align_prob_mat.iter() {
        buf_4_rna_id_pair.push_str(&format!("{},{},{} ", i - T::one(), j - T::one(), align_prob));
      }
      buf_4_writer_2_align_prob_mat_file.push_str(&buf_4_rna_id_pair);
    }
    let _ = writer_2_loop_align_prob_mat_file.write_all(buf_4_writer_2_loop_align_prob_mat_file.as_bytes());
    let _ = writer_2_basepair_align_prob_mat_file.write_all(buf_4_writer_2_basepair_align_prob_mat_file.as_bytes());
    let _ = writer_2_align_prob_mat_file.write_all(buf_4_writer_2_align_prob_mat_file.as_bytes());
  }
}

pub fn get_sparse_pos_sets<T>(bpp_mat_pair: &ProbMatPair<T>, align_prob_mat: &SparseProbMat<T>, seq_len_pair: &(T, T)) -> (PosPairMatSet<T>, PosPairMatSet<T>, PosQuadrupleMat<T>, PosQuadrupleMatWithLenPairs<T>, SparsePosSets<T>, SparsePosSets<T>)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let mut forward_pos_pair_mat_set = PosPairMatSet::<T>::default();
  let mut backward_pos_pair_mat_set = PosPairMatSet::<T>::default();
  let mut pos_quadruple_mat = PosQuadrupleMat::<T>::default();
  let mut pos_quadruple_mat_with_len_pairs = PosQuadrupleMatWithLenPairs::<T>::default();
  let mut matchable_pos_sets_1 = SparsePosSets::<T>::default();
  let mut matchable_pos_sets_2 = SparsePosSets::<T>::default();
  for pos_pair in bpp_mat_pair.0.keys() {
    for pos_pair_2 in bpp_mat_pair.1.keys() {
      let forward_pos_pair = (pos_pair.0, pos_pair_2.0);
      let backward_pos_pair = (pos_pair.1, pos_pair_2.1);
      if !align_prob_mat.contains_key(&forward_pos_pair) {continue;}
      if !align_prob_mat.contains_key(&backward_pos_pair) {continue;}
      match forward_pos_pair_mat_set.get_mut(&backward_pos_pair) {
        Some(forward_pos_pair_mat) => {
          forward_pos_pair_mat.insert(forward_pos_pair);
        }, None => {
          let mut empty_mat = PosPairMat::<T>::default();
          empty_mat.insert(forward_pos_pair);
          forward_pos_pair_mat_set.insert(backward_pos_pair, empty_mat);
        },
      }
      match backward_pos_pair_mat_set.get_mut(&forward_pos_pair) {
        Some(backward_pos_pair_mat) => {
          backward_pos_pair_mat.insert(backward_pos_pair);
        }, None => {
          let mut empty_mat = PosPairMat::<T>::default();
          empty_mat.insert(backward_pos_pair);
          backward_pos_pair_mat_set.insert(forward_pos_pair, empty_mat);
        },
      }
      pos_quadruple_mat.insert((pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1));
      let substr_len_pair = (pos_pair.1 - pos_pair.0 + T::one(), pos_pair_2.1 - pos_pair_2.0 + T::one());
      let left_pos_pair = (pos_pair.0, pos_pair_2.0);
      match pos_quadruple_mat_with_len_pairs.get_mut(&substr_len_pair) {
        Some(pos_pairs) => {
          pos_pairs.insert(left_pos_pair);
        }, None => {
          let mut empty_mat = PosPairMat::<T>::default();
          empty_mat.insert(left_pos_pair);
          pos_quadruple_mat_with_len_pairs.insert(substr_len_pair, empty_mat);
        },
      }
    }
  }
  for pos_pair in align_prob_mat.keys() {
    match matchable_pos_sets_1.get_mut(&pos_pair.0) {
      Some(matchable_poss) => {
        matchable_poss.insert(pos_pair.1);
      }, None => {
        let mut matchable_poss = SparsePoss::<T>::default();
        matchable_poss.insert(pos_pair.1);
        matchable_pos_sets_1.insert(pos_pair.0, matchable_poss);
      },
    }
  }
  for pos_pair in align_prob_mat.keys() {
    match matchable_pos_sets_2.get_mut(&pos_pair.1) {
      Some(matchable_poss) => {
        matchable_poss.insert(pos_pair.0);
      }, None => {
        let mut matchable_poss = SparsePoss::<T>::default();
        matchable_poss.insert(pos_pair.0);
        matchable_pos_sets_2.insert(pos_pair.1, matchable_poss);
      },
    }
  }
  let mut matchable_poss = SparsePoss::<T>::default();
  matchable_poss.insert(T::zero());
  matchable_pos_sets_1.insert(T::zero(), matchable_poss.clone());
  matchable_pos_sets_2.insert(T::zero(), matchable_poss);
  let mut matchable_poss = SparsePoss::<T>::default();
  matchable_poss.insert(seq_len_pair.1 - T::one());
  matchable_pos_sets_1.insert(seq_len_pair.0 - T::one(), matchable_poss);
  let mut matchable_poss = SparsePoss::<T>::default();
  matchable_poss.insert(seq_len_pair.0 - T::one());
  matchable_pos_sets_2.insert(seq_len_pair.1 - T::one(), matchable_poss);
  (forward_pos_pair_mat_set, backward_pos_pair_mat_set, pos_quadruple_mat, pos_quadruple_mat_with_len_pairs, matchable_pos_sets_1, matchable_pos_sets_2)
}

pub fn write_readme(output_dir_path: &Path, readme_contents: &String) {
  let readme_file_path = output_dir_path.join(README_FILE_NAME);
  let mut writer_2_readme_file = BufWriter::new(File::create(readme_file_path).unwrap());
  let buf_4_writer_2_readme_file = String::from(readme_contents);
  let _ = writer_2_readme_file.write_all(buf_4_writer_2_readme_file.as_bytes());
}
