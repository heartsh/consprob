extern crate rna_algos;
#[macro_use]
extern crate lazy_static;
extern crate getopts;
extern crate scoped_threadpool;
extern crate itertools;
extern crate bio;
extern crate num_cpus;

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

pub type Prob4dMat = HashMap<PosQuadruple, Prob>;
type PartFunc4dMat = HashMap<PosQuadruple, PartFunc>;
pub type TmpPartFuncSetMat = HashMap<PosPair, TmpPartFuncSets>;
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
pub type TmpPartFuncSetMat4El = HashMap<PosPair, TmpPartFuncs>;
#[derive(Clone)]
pub struct StaPartFuncMats {
  pub part_func_4d_mat_4_bpas: PartFunc4dMat,
  pub part_func_4d_mat_4_bpas_accessible_on_els: PartFunc4dMat,
  pub part_func_4d_mat_4_bpas_accessible_on_mls: PartFunc4dMat,
  pub forward_part_func_set_mat_4_external_loop: TmpPartFuncSetMat4El,
  pub backward_part_func_set_mat_4_external_loop: TmpPartFuncSetMat4El,
  pub forward_tmp_part_func_set_mats_with_pos_quadruples: TmpPartFuncSetMatsWithPosQuadruples,
  pub backward_tmp_part_func_set_mats_with_pos_quadruples: TmpPartFuncSetMatsWithPosQuadruples,
}
pub type TmpPartFuncSetMatsWithPosQuadruples = HashMap<PosQuadruple, TmpPartFuncSetMat>;
pub struct StaFeParams {
  pub ba_score_mat: SparseFreeEnergyMat,
  pub bpa_score_mat: FreeEnergy4dMat,
  pub insert_scores: FreeEnergies,
  pub insert_scores_2: FreeEnergies,
}
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type StaFeParamSetsWithRnaIdPairs = HashMap<RnaIdPair, StaFeParams>;
pub type Prob4dMatsWithRnaIdPairs = HashMap<RnaIdPair, Prob4dMat>;
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
pub type SparseFreeEnergyMat = HashMap<PosPair, FreeEnergy>;
pub type PosPairsWithPosPairs = HashMap<PosPair, PosPair>;
pub type BoolsWithPosPairs = HashMap<PosPair, bool>;
pub type ProbMatPair<'a> =  (&'a SparseProbMat, &'a SparseProbMat);
pub type SsFreeEnergyMatSetPair<'a> =  (&'a SsFreeEnergyMats, &'a SsFreeEnergyMats);
pub type NumOfThreads = u32;
pub type FreeEnergySetPair<'a> = (&'a FreeEnergies, &'a FreeEnergies);
pub struct StaProbMats {
  pub bpp_mat_pair: SparseProbMatPair,
  pub access_bpp_mat_pair_4_2l: SparseProbMatPair,
  pub access_bpp_mat_pair_4_ml: SparseProbMatPair,
  pub bpp_mat_pair_4_el: SparseProbMatPair,
  pub upp_mat_pair: ProbSetPair,
  pub upp_mat_pair_4_hl: ProbSetPair,
  pub upp_mat_pair_4_2l: ProbSetPair,
  pub upp_mat_pair_4_ml: ProbSetPair,
  pub upp_mat_pair_4_el: ProbSetPair,
}
#[derive(Clone)]
pub struct PctStaProbMats {
  pub max_bpp_mat: SparseProbMat,
  pub bpp_mat: SparseProbMat,
  pub access_bpp_mat_4_2l: SparseProbMat,
  pub access_bpp_mat_4_ml: SparseProbMat,
  pub bpp_mat_4_el: SparseProbMat,
  pub bpp_mat_on_ss: SparseProbMat,
  pub max_upp_mat: Probs,
  pub upp_mat: Probs,
  pub upp_mat_4_hl: Probs,
  pub upp_mat_4_2l: Probs,
  pub upp_mat_4_ml: Probs,
  pub upp_mat_4_el: Probs,
}
pub type SparseProbMatPair = (SparseProbMat, SparseProbMat);
pub type ProbSetPair = (Probs, Probs);
pub type ProbMatSets = Vec<PctStaProbMats>;
pub type StaProbMatsWithRnaIdPairs = HashMap<RnaIdPair, StaProbMats>;
pub type ProbSeqPair<'a> = (&'a Probs, &'a Probs);
pub type Poss = Vec<Pos>;
pub type TmpPartFuncSetMatsWithPosPairs = HashMap<PosPair, TmpPartFuncSetMat>;

impl StaProbMats {
  pub fn origin() -> StaProbMats {
    let prob_mat_pair = (SparseProbMat::default(), SparseProbMat::default());
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
  pub fn new(seq_len_pair: &(u16, u16)) -> StaProbMats {
    let prob_mat_pair = (SparseProbMat::default(), SparseProbMat::default());
    let prob_set_pair = (vec![NEG_INFINITY; seq_len_pair.0 as usize], vec![NEG_INFINITY; seq_len_pair.1 as usize]);
    StaProbMats {
      bpp_mat_pair: prob_mat_pair.clone(),
      access_bpp_mat_pair_4_2l: prob_mat_pair.clone(),
      access_bpp_mat_pair_4_ml: prob_mat_pair.clone(),
      bpp_mat_pair_4_el: prob_mat_pair,
      upp_mat_pair: (vec![1.; seq_len_pair.0 as usize], vec![1.; seq_len_pair.1 as usize]),
      upp_mat_pair_4_hl: prob_set_pair.clone(),
      upp_mat_pair_4_2l: prob_set_pair.clone(),
      upp_mat_pair_4_ml: prob_set_pair.clone(),
      upp_mat_pair_4_el: prob_set_pair,
    }
  }
}

impl PctStaProbMats {
  pub fn origin() -> PctStaProbMats {
    let prob_mat = SparseProbMat::default();
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
  pub fn new(seq_len: usize) -> PctStaProbMats {
    let prob_mat = SparseProbMat::default();
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

impl StaFeParams {
  pub fn origin() -> StaFeParams {
    let sparse_free_energy_mat = SparseFreeEnergyMat::default();
    let free_energy_4d_mat = FreeEnergy4dMat::default();
    let insert_scores = Vec::new();
    StaFeParams {
      ba_score_mat: sparse_free_energy_mat.clone(),
      bpa_score_mat: free_energy_4d_mat.clone(),
      insert_scores: insert_scores.clone(),
      insert_scores_2: insert_scores,
    }
  }

  pub fn new(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), max_bp_span_pair: &PosPair, max_gap_num: Pos, max_gap_num_4_il: Pos, bpp_mat_pair: &ProbMatPair, uses_bpps: bool) -> StaFeParams {
    let max = max(max_gap_num, max_gap_num_4_il);
    let mut sta_fe_params = StaFeParams::origin();
    sta_fe_params.insert_scores = vec![0.; seq_len_pair.0];
    sta_fe_params.insert_scores_2 = vec![0.; seq_len_pair.1];
    for j in 1 .. seq_len_pair.1 - 1 {
      let base = seq_pair.1[j];
      sta_fe_params.insert_scores_2[j] = INSERT_SCORES[base];
    }
    let seq_len_pair = (seq_len_pair.0 as Pos, seq_len_pair.1 as Pos);
    let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
    for i in 1 .. seq_len_pair.0 - 1 {
      let long_i = i as usize;
      let base = seq_pair.0[long_i];
      sta_fe_params.insert_scores[long_i] = INSERT_SCORES[base];
      for j in 1 .. seq_len_pair.1 - 1 {
        let long_j = j as usize;
        let pos_pair = (i, j);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max) {continue;}
        let base_pair = (base, seq_pair.1[long_j]);
        sta_fe_params.ba_score_mat.insert(pos_pair, BA_SCORE_MAT[&base_pair] + RIBOSUM_BA_SCORE_MAT[&base_pair]);
      }
      let upper_j = if i + max_bp_span_pair.0 >= seq_len_pair.0 - 1 {seq_len_pair.0 - 1} else {i + max_bp_span_pair.0};
      for j in i + 1 .. upper_j {
        let long_j = j as usize;
        let pos_pair = (i, j);
        let base_pair = (base, seq_pair.0[long_j]);
        if !is_canonical(&base_pair) {continue;}
        if !bpp_mat_pair.0.contains_key(&pos_pair) {continue;}
        let bpp = if uses_bpps {bpp_mat_pair.0[&pos_pair].ln()} else {0.};
        for k in 1 .. seq_len_pair.1 - 1 {
          if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
          let long_k = k as usize;
          let base_2 = seq_pair.1[long_k];
          let align_score = BA_SCORE_MAT[&(base, base_2)];
          let upper_l = if k + max_bp_span_pair.1 >= seq_len_pair.1 - 1 {seq_len_pair.1 - 1} else {k + max_bp_span_pair.1};
          for l in k + 1 .. upper_l {
            if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
            let pos_pair_2 = (k, l);
            let pos_quadruple = (pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1);
            let long_l = l as usize;
            let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
            if !is_canonical(&base_pair_2) {continue;}
            if !(is_min_gap_ok_2(&pos_quadruple, max_gap_num_4_il) && bpp_mat_pair.1.contains_key(&pos_pair_2)) {continue;}
            let bpp_2 = if uses_bpps {bpp_mat_pair.1[&pos_pair_2].ln()} else {0.};
            sta_fe_params.bpa_score_mat.insert(pos_quadruple, align_score + BA_SCORE_MAT[&(base_pair.1, base_pair_2.1)] + bpp + bpp_2 + RIBOSUM_BPA_SCORE_MAT[&(base_pair, base_pair_2)]);
          }
        }
      }
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

impl StaPartFuncMats {
  pub fn new() -> StaPartFuncMats {
    let part_func_4d_mat = PartFunc4dMat::default();
    let part_func_set_mat = TmpPartFuncSetMat4El::new();
    let tmp_part_func_set_mats_with_pos_quadruples = TmpPartFuncSetMatsWithPosQuadruples::default();
    StaPartFuncMats {
      part_func_4d_mat_4_bpas: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_els: part_func_4d_mat.clone(),
      part_func_4d_mat_4_bpas_accessible_on_mls: part_func_4d_mat,
      forward_part_func_set_mat_4_external_loop: part_func_set_mat.clone(),
      backward_part_func_set_mat_4_external_loop: part_func_set_mat,
      forward_tmp_part_func_set_mats_with_pos_quadruples: tmp_part_func_set_mats_with_pos_quadruples.clone(),
      backward_tmp_part_func_set_mats_with_pos_quadruples: tmp_part_func_set_mats_with_pos_quadruples,
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

pub const MAX_GAP_NUM_4_IL: Pos = 100;
pub const MIN_GAP_NUM_4_IL: Pos = 5;
pub const DEFAULT_MIN_BPP: Prob = 0.005;
pub const DEFAULT_OFFSET_4_MAX_GAP_NUM: Pos = 5;
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

pub fn io_algo_4_prob_mats(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &PosPair, max_gap_num: Pos, max_gap_num_4_il: Pos, bpp_mat_pair: &ProbMatPair, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair, uses_bpps: bool, produces_access_probs: bool) -> StaProbMats {
  let (sta_part_func_mats, global_part_func) = get_sta_inside_part_func_mats(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, max_gap_num_4_il, bpp_mat_pair, ss_free_energy_mat_set_pair, false, uses_bpps);
  get_sta_prob_mats(seq_pair, seq_len_pair, sta_fe_params, max_bp_span_pair, max_gap_num, max_gap_num_4_il, &sta_part_func_mats, bpp_mat_pair, ss_free_energy_mat_set_pair, uses_bpps, produces_access_probs, global_part_func)
}

pub fn get_sta_inside_part_func_mats(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &PosPair, max_gap_num: Pos, max_gap_num_4_il: Pos, bpp_mat_pair: &ProbMatPair, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair, is_viterbi: bool, uses_bpps: bool) -> (StaPartFuncMats, PartFunc) {
  let seq_len_pair = (seq_len_pair.0 as Pos, seq_len_pair.1 as Pos);
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut sta_part_func_mats = StaPartFuncMats::new();
  let mut forward_tmp_part_func_set_mats_with_pos_pairs = TmpPartFuncSetMatsWithPosPairs::default();
  let mut backward_tmp_part_func_set_mats_with_pos_pairs = forward_tmp_part_func_set_mats_with_pos_pairs.clone();
  for substr_len_1 in MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.0 + 1 {
    for i in 1 .. seq_len_pair.0 - substr_len_1 {
      let j = i + substr_len_1 - 1;
      let (long_i, long_j) = (i as usize, j as usize);
      let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
      if !is_canonical(&base_pair) {continue;}
      if !bpp_mat_pair.0.contains_key(&(i, j)) {continue;}
      let invert_base_pair = invert_bp(&base_pair);
      let invert_stacking_base_pair = invert_bp(&(seq_pair.0[long_i + 1], seq_pair.0[long_j - 1]));
      let ml_tm_delta_fe = if uses_bpps {0.} else {ML_TM_DELTA_FES[invert_base_pair.0][invert_base_pair.1][invert_stacking_base_pair.0][invert_stacking_base_pair.1]};
      let stacking_bp = (seq_pair.0[long_i - 1], seq_pair.0[long_j + 1]);
      let ml_tm_or_de_delta_fe = if uses_bpps {0.} else {
          if i > 1 && j < seq_len_pair.0 - 2 {
          ML_TM_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.0][stacking_bp.1]
        } else if i > 1 {
          FIVE_PRIME_DE_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.0]
        } else if j < seq_len_pair.0 - 2 {
          THREE_PRIME_DE_DELTA_FES[base_pair.0][base_pair.1][stacking_bp.1]
        } else {
          0.
        }
      };
      let au_or_gu_end_penalty_delta_fe = if uses_bpps {0.} else {
        if is_au_or_gu(&base_pair) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
      };
      let hl_fe = if uses_bpps {0.} else {ss_free_energy_mat_set_pair.0.hl_fe_mat[&(i, j)]};
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
              let (forward_tmp_part_func_set_mat, part_func_on_sa, part_func_4_ml) = get_tmp_part_func_set_mat(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_part_func_mats, bpp_mat_pair, is_viterbi, true, &forward_tmp_part_func_set_mats_with_pos_pairs);
              let (backward_tmp_part_func_set_mat, _, _) = get_tmp_part_func_set_mat(&seq_len_pair, sta_fe_params, max_gap_num_4_il, &pos_quadruple, &sta_part_func_mats, bpp_mat_pair, is_viterbi, false, &backward_tmp_part_func_set_mats_with_pos_pairs);
              let mut sum = NEG_INFINITY;
              let score = bpa_score + hl_fe + if uses_bpps {0.} else {ss_free_energy_mat_set_pair.1.hl_fe_mat[&(k, l)]} + part_func_on_sa;
              sumormax(&mut sum, score, is_viterbi);
              for m in i + 1 .. j {
                let long_m = m as usize;
                for n in m + 1 .. j {
                  let long_n = n as usize;
                  let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                  if !is_canonical(&base_pair_3) {continue;}
                  if !bpp_mat_pair.0.contains_key(&(m, n)) {continue;}
                  if long_m - long_i - 1 + long_j - long_n - 1 > MAX_2_LOOP_LEN {continue;}
                  let twoloop_fe = if uses_bpps {0.} else {ss_free_energy_mat_set_pair.0.twoloop_fe_4d_mat[&(i, j, m, n)]};
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
                      match sta_part_func_mats.part_func_4d_mat_4_bpas.get(&pos_quadruple_2) {
                        Some(&part_func) => {
                          let mut forward_term = NEG_INFINITY;
                          let mut backward_term = forward_term;
                          match forward_tmp_part_func_set_mat.get(&(m - 1, o - 1)) {
                            Some(part_func_sets) => {
                              let ref part_funcs = part_func_sets.part_funcs_on_sa;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                            }, None => {},
                          }
                          match backward_tmp_part_func_set_mat.get(&(n + 1, p + 1)) {
                            Some(part_func_sets) => {
                              let ref part_funcs = part_func_sets.part_funcs_on_sa;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                            }, None => {},
                          }
                          let part_func_4_2l = forward_term + backward_term;
                          let twoloop_fe_2 = if uses_bpps {0.} else {ss_free_energy_mat_set_pair.1.twoloop_fe_4d_mat[&(k, l, o, p)]};
                          let coefficient = bpa_score + twoloop_fe + twoloop_fe_2 + part_func;
                          sumormax(&mut sum, coefficient + part_func_4_2l, is_viterbi);
                        }, None => {},
                      }
                    }
                  }
                }
              }
              let au_or_gu_end_penalty_delta_fe_2 = if uses_bpps {0.} else {
                if is_au_or_gu(&base_pair_2) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
              };
              let invert_base_pair_2 = invert_bp(&base_pair_2);
              let invert_stacking_base_pair_2 = invert_bp(&(seq_pair.1[long_k + 1], seq_pair.1[long_l - 1]));
              let ml_tm_delta_fe_2 = if uses_bpps {0.} else {ML_TM_DELTA_FES[invert_base_pair_2.0][invert_base_pair_2.1][invert_stacking_base_pair_2.0][invert_stacking_base_pair_2.1]};
              let score = bpa_score + if uses_bpps {0.} else {
                2. * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe + ml_tm_delta_fe_2 + au_or_gu_end_penalty_delta_fe + au_or_gu_end_penalty_delta_fe_2
              } + part_func_4_ml;
              sumormax(&mut sum, score, is_viterbi);
              if sum > NEG_INFINITY {
                sta_part_func_mats.part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
                let stacking_bp_2 = (seq_pair.1[long_k - 1], seq_pair.1[long_l + 1]);
                let ml_tm_or_de_delta_fe_2 = if uses_bpps {0.} else {
                  if k > 1 && l < seq_len_pair.1 - 2 {
                    ML_TM_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.0][stacking_bp_2.1]
                  } else if k > 1 {
                    FIVE_PRIME_DE_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.0]
                  } else if l < seq_len_pair.1 - 2 {
                    THREE_PRIME_DE_DELTA_FES[base_pair_2.0][base_pair_2.1][stacking_bp_2.1]
                  } else {
                    0.
                  }
                };
                sum += ml_tm_or_de_delta_fe + ml_tm_or_de_delta_fe_2 + au_or_gu_end_penalty_delta_fe + au_or_gu_end_penalty_delta_fe_2;
                sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.insert(pos_quadruple, sum);
                sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls.insert(pos_quadruple, sum + 2. * COEFFICIENT_4_TERM_OF_NUM_OF_BRANCHING_HELICES_ON_INIT_ML_DELTA_FE);
              }
              sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_quadruples.insert(pos_quadruple, forward_tmp_part_func_set_mat.clone());
              sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_quadruples.insert(pos_quadruple, backward_tmp_part_func_set_mat.clone());
              match forward_tmp_part_func_set_mats_with_pos_pairs.get_mut(&(i, k)) {
                Some(part_func_set_mat) => {
                  *part_func_set_mat = forward_tmp_part_func_set_mat;
                },
                None => {
                  forward_tmp_part_func_set_mats_with_pos_pairs.insert((i, k), forward_tmp_part_func_set_mat);
                },
              }
              match backward_tmp_part_func_set_mats_with_pos_pairs.get_mut(&(j, l)) {
                Some(part_func_set_mat) => {
                  *part_func_set_mat = backward_tmp_part_func_set_mat;
                },
                None => {
                  backward_tmp_part_func_set_mats_with_pos_pairs.insert((j, l), backward_tmp_part_func_set_mat);
                },
              }
            }, None => {},
          }
        }
      }
    }
  }
  let leftmost_pos_pair = (0, 0);
  let rightmost_pos_pair = (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
  let mut part_funcs = TmpPartFuncs::new();
  part_funcs.part_func_4_align = 0.;
  sta_part_func_mats.forward_part_func_set_mat_4_external_loop.insert(leftmost_pos_pair, part_funcs.clone());
  sta_part_func_mats.backward_part_func_set_mat_4_external_loop.insert(rightmost_pos_pair, part_funcs);
  for i in 0 .. seq_len_pair.0 - 1 {
    let long_i = i as usize;
    let insert_score = sta_fe_params.insert_scores[long_i];
    for j in 0 .. seq_len_pair.1 - 1 {
      let pos_pair = (i, j);
      if pos_pair == (0, 0) {continue;}
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
      let mut part_funcs = TmpPartFuncs::new();
      let mut sum = NEG_INFINITY;
      for k in 1 .. i {
        if !bpp_mat_pair.0.contains_key(&(k, i)) {continue;}
        for l in 1 .. j {
          let pos_pair_2 = (k - 1, l - 1);
          if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num) {continue;}
          let pos_quadruple = (k, i, l, j);
          let is_begin = pos_pair_2 == leftmost_pos_pair;
          match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
            Some(&part_func) => {
              match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
                Some(part_funcs) => {
                  let score = part_funcs.part_func_4_align + part_func + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
                  sumormax(&mut sum, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert + part_func + INSERT_2_MATCH_SCORE;
                  sumormax(&mut sum, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert_2 + part_func + INSERT_2_MATCH_SCORE;
                  sumormax(&mut sum, score, is_viterbi);
                }, None => {},
              }
            }, None => {},
          }
        }
      }
      if i > 0 && j > 0 {
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        let pos_pair_2 = (i - 1, j - 1);
        let is_begin = pos_pair_2 == leftmost_pos_pair;
        match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + ba_score + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert + ba_score + INSERT_2_MATCH_SCORE;
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert_2 + ba_score + INSERT_2_MATCH_SCORE;
            sumormax(&mut sum, score, is_viterbi);
          }, None => {},
        }
        part_funcs.part_func_4_align = sum;
      }
      if i > 0 {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i - 1, j);
        let is_begin = pos_pair_2 == leftmost_pos_pair;
        match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
            sumormax(&mut sum, score, is_viterbi);
          }, None => {},
        }
        let sum = sum + insert_score;
        part_funcs.part_func_4_insert = sum;
      }
      if j > 0 {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i, j - 1);
        let is_begin = pos_pair_2 == leftmost_pos_pair;
        match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
            sumormax(&mut sum, score, is_viterbi);
          }, None => {},
        }
        let long_j = j as usize;
        let insert_score_2 = sta_fe_params.insert_scores_2[long_j];
        let sum = sum + insert_score_2;
        part_funcs.part_func_4_insert_2 = sum;
      }
      sta_part_func_mats.forward_part_func_set_mat_4_external_loop.insert(pos_pair, part_funcs);
    }
  }
  let mut final_sum = NEG_INFINITY;
  let ref part_funcs = sta_part_func_mats.forward_part_func_set_mat_4_external_loop[&(seq_len_pair.0 - 2, seq_len_pair.1 - 2)];
  sumormax(&mut final_sum, part_funcs.part_func_4_align, is_viterbi);
  sumormax(&mut final_sum, part_funcs.part_func_4_insert, is_viterbi);
  sumormax(&mut final_sum, part_funcs.part_func_4_insert_2, is_viterbi);
  for i in (1 .. seq_len_pair.0).rev() {
    let long_i = i as usize;
    let insert_score = sta_fe_params.insert_scores[long_i];
    for j in (1 .. seq_len_pair.1).rev() {
      let pos_pair = (i, j);
      if pos_pair == (seq_len_pair.0 - 1, seq_len_pair.1 - 1) {continue;}
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
      let mut part_funcs = TmpPartFuncs::new();
      let mut sum = NEG_INFINITY;
      for k in i + 1 .. seq_len_pair.0 - 1 {
        if !bpp_mat_pair.0.contains_key(&(i, k)) {continue;}
        for l in j + 1 .. seq_len_pair.1 - 1 {
          let pos_pair_2 = (k + 1, l + 1);
          if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num) {continue;}
          let pos_quadruple = (i, k, j, l);
          let is_end = pos_pair_2 == rightmost_pos_pair;
          match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_els.get(&pos_quadruple) {
            Some(&part_func) => {
              match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
                Some(part_funcs) => {
                  let score = part_funcs.part_func_4_align + part_func + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                  sumormax(&mut sum, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert + part_func + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                  sumormax(&mut sum, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert_2 + part_func + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                  sumormax(&mut sum, score, is_viterbi);
                }, None => {},
              }
            }, None => {},
          }
        }
      }
      if i < seq_len_pair.0 - 1 && j < seq_len_pair.1 - 1 {
        let pos_pair_2 = (i + 1, j + 1);
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        let is_end = pos_pair_2 == rightmost_pos_pair;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + ba_score + if is_end {0.} else {MATCH_2_MATCH_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert + ba_score + if is_end {0.} else {MATCH_2_INSERT_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert_2 + ba_score + if is_end {0.} else {MATCH_2_INSERT_SCORE};
            sumormax(&mut sum, score, is_viterbi);
          }, None => {},
        }
        part_funcs.part_func_4_align = sum;
      }
      if i < seq_len_pair.0 - 1 {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i + 1, j);
        let is_end = pos_pair_2 == rightmost_pos_pair;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert + if is_end {0.} else {INSERT_EXTEND_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert_2 + if is_end {0.} else {INSERT_SWITCH_SCORE};
            sumormax(&mut sum, score, is_viterbi);
          }, None => {},
        }
        part_funcs.part_func_4_insert = sum + insert_score;
      }
      if j < seq_len_pair.1 - 1 {
        sum = NEG_INFINITY;
        let pos_pair_2 = (i, j + 1);
        let is_end = pos_pair_2 == rightmost_pos_pair;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            let score = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert + if is_end {0.} else {INSERT_SWITCH_SCORE};
            sumormax(&mut sum, score, is_viterbi);
            let score = part_funcs.part_func_4_insert_2 + if is_end {0.} else {INSERT_EXTEND_SCORE};
            sumormax(&mut sum, score, is_viterbi);
          }, None => {},
        }
        let long_j = j as usize;
        let insert_score_2 = sta_fe_params.insert_scores_2[long_j];
        part_funcs.part_func_4_insert_2 = sum + insert_score_2;
      }
      sta_part_func_mats.backward_part_func_set_mat_4_external_loop.insert(pos_pair, part_funcs);
    }
  }
  (sta_part_func_mats, final_sum)
}

pub fn get_tmp_part_func_set_mat(seq_len_pair: &PosPair, sta_fe_params: &StaFeParams, max_gap_num_4_il: Pos, pos_quadruple: &PosQuadruple, sta_part_func_mats: &StaPartFuncMats, bpp_mat_pair: &ProbMatPair, is_viterbi: bool, is_forward: bool, tmp_part_func_set_mats_with_pos_pairs: &TmpPartFuncSetMatsWithPosPairs) -> (TmpPartFuncSetMat, PartFunc, PartFunc) {
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut cache_is_used = false;
  let &(i, j, k, l) = pos_quadruple;
  let leftmost_pos_pair = if is_forward {(i, k)} else {(i + 1, k + 1)};
  let rightmost_pos_pair = if is_forward {(j - 1, l - 1)} else {(j, l)};
  let mut tmp_part_func_set_mat = match tmp_part_func_set_mats_with_pos_pairs.get(&if is_forward {leftmost_pos_pair} else {rightmost_pos_pair}) {
    Some(cache) => {
      cache_is_used = true;
      cache.clone()
    }, None => {
      TmpPartFuncSetMat::new()
    },
  };
  let iter: Poss = if is_forward {(i .. j).collect()} else {(i + 1 .. j + 1).rev().collect()};
  let iter_2: Poss = if is_forward {(k .. l).collect()} else {(k + 1 .. l + 1).rev().collect()};
  for &u in iter.iter() {
    let long_u = u as usize;
    let insert_score = sta_fe_params.insert_scores[long_u];
    for &v in iter_2.iter() {
      let pos_pair = (u, v);
      if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
      if cache_is_used && tmp_part_func_set_mat.contains_key(&pos_pair) {
        continue;
      }
      let mut tmp_part_func_sets = TmpPartFuncSets::new();
      if (is_forward && u == i && v == k) || (!is_forward && u == j && v == l) {
        tmp_part_func_sets.part_funcs_on_sa.part_func_4_align = 0.;
        tmp_part_func_set_mat.insert(pos_pair, tmp_part_func_sets);
        continue;
      }
      let long_v = v as usize;
      let mut sum_on_sa = NEG_INFINITY;
      let mut sum_4_ml = sum_on_sa;
      let mut sum_4_first_bpas_on_mls = sum_on_sa;
      let mut tmp_sum = sum_on_sa;
      // For alignments.
      for m in if is_forward {i + 1 .. u} else {u + 1 .. j} {
        if !bpp_mat_pair.0.contains_key(&if is_forward{(m, u)} else {(u, m)}) {continue;}
        for n in if is_forward {k + 1 .. v} else {v + 1 .. l} {
          let pos_pair_2 = if is_forward {(m - 1, n - 1)} else {(m + 1, n + 1)};
          if !is_min_gap_ok_1(&pos_pair_2, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
          let pos_quadruple_2 = if is_forward {(m, u, n, v)} else {(u, m, v, n)};
          match sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls.get(&pos_quadruple_2) {
            Some(&part_func) => {
              let is_begin = pos_pair_2 == if is_forward {leftmost_pos_pair} else {rightmost_pos_pair};
              match tmp_part_func_set_mat.get(&pos_pair_2) {
                Some(part_func_sets) => {
                  let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                  let score = part_funcs.part_func_4_align + part_func + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
                  sumormax(&mut sum_4_ml, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert + part_func + if is_forward {
                    INSERT_2_MATCH_SCORE
                  } else {
                    MATCH_2_INSERT_SCORE
                  };
                  sumormax(&mut sum_4_ml, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert_2 + part_func + if is_forward {
                    INSERT_2_MATCH_SCORE
                  } else {
                    MATCH_2_INSERT_SCORE
                  };
                  sumormax(&mut sum_4_ml, score, is_viterbi);
                  let ref part_funcs = part_func_sets.part_funcs_on_sa;
                  let score = part_funcs.part_func_4_align + part_func + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
                  sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert + part_func + if is_forward {
                    INSERT_2_MATCH_SCORE
                  } else {
                    MATCH_2_INSERT_SCORE
                  };
                  sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
                  let score = part_funcs.part_func_4_insert_2 + part_func + if is_forward {
                    INSERT_2_MATCH_SCORE
                  } else {
                    MATCH_2_INSERT_SCORE
                  };
                  sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
                }, None => {},
              }
            }, None => {},
          }
        }
      }
      let pos_pair_2 = if is_forward {(u - 1, v - 1)} else {(u + 1, v + 1)};
      let is_begin = if is_forward {pos_pair_2 == leftmost_pos_pair} else {pos_pair_2 == rightmost_pos_pair};
      let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
      match tmp_part_func_set_mat.get(&pos_pair_2) {
        Some(part_func_sets) => {
          let ref part_funcs = part_func_sets.part_funcs_4_ml;
          let score = part_funcs.part_func_4_align + ba_score + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
          let score = part_funcs.part_func_4_align + ba_score + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + ba_score + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let ref part_funcs = part_func_sets.part_funcs_on_sa;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {MATCH_2_MATCH_SCORE};
          sumormax(&mut sum_on_sa, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          sumormax(&mut sum_on_sa, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + if is_forward {
            INSERT_2_MATCH_SCORE
          } else {
            MATCH_2_INSERT_SCORE
          };
          sumormax(&mut sum_on_sa, score, is_viterbi);
        }, None => {},
      }
      tmp_part_func_sets.part_funcs_4_ml.part_func_4_align = sum_4_ml;
      sumormax(&mut tmp_sum, sum_4_ml, is_viterbi);
      tmp_part_func_sets.part_funcs_4_first_bpas_on_mls.part_func_4_align = sum_4_first_bpas_on_mls;
      sumormax(&mut tmp_sum, sum_4_first_bpas_on_mls, is_viterbi);
      tmp_part_func_sets.part_funcs_4_bpas_on_mls.part_func_4_align = tmp_sum;
      let sum_on_sa = sum_on_sa + ba_score;
      tmp_part_func_sets.part_funcs_on_sa.part_func_4_align = sum_on_sa;
      sumormax(&mut tmp_sum, sum_on_sa, is_viterbi);
      tmp_part_func_sets.part_funcs_on_mls.part_func_4_align = tmp_sum;
      // For inserts.
      let mut sum_on_sa = NEG_INFINITY;
      let mut sum_4_ml = sum_on_sa;
      let mut sum_4_first_bpas_on_mls = sum_on_sa;
      let mut tmp_sum = sum_on_sa;
      let pos_pair_2 = if is_forward {(u - 1, v)} else {(u + 1, v)};
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
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let ref part_funcs = part_func_sets.part_funcs_on_sa;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          sumormax(&mut sum_on_sa, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
          sumormax(&mut sum_on_sa, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
          sumormax(&mut sum_on_sa, score, is_viterbi);
        }, None => {},
      }
      let sum_4_ml = sum_4_ml + insert_score;
      tmp_part_func_sets.part_funcs_4_ml.part_func_4_insert = sum_4_ml;
      sumormax(&mut tmp_sum, sum_4_ml, is_viterbi);
      let sum_4_first_bpas_on_mls = sum_4_first_bpas_on_mls + insert_score;
      tmp_part_func_sets.part_funcs_4_first_bpas_on_mls.part_func_4_insert = sum_4_first_bpas_on_mls;
      sumormax(&mut tmp_sum, sum_4_first_bpas_on_mls, is_viterbi);
      tmp_part_func_sets.part_funcs_4_bpas_on_mls.part_func_4_insert = tmp_sum;
      let sum_on_sa = sum_on_sa + insert_score;
      tmp_part_func_sets.part_funcs_on_sa.part_func_4_insert = sum_on_sa;
      sumormax(&mut tmp_sum, sum_on_sa, is_viterbi);
      tmp_part_func_sets.part_funcs_on_mls.part_func_4_insert = tmp_sum;
      // For inserts on the other side.
      let mut sum_on_sa = NEG_INFINITY;
      let mut sum_4_ml = sum_on_sa;
      let mut sum_4_first_bpas_on_mls = sum_on_sa;
      let mut tmp_sum = sum_on_sa;
      let pos_pair_2 = if is_forward {(u, v - 1)} else {(u, v + 1)};
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
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
          sumormax(&mut sum_4_ml, score, is_viterbi);
          let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
          sumormax(&mut sum_4_first_bpas_on_mls, score, is_viterbi);
          let ref part_funcs = part_func_sets.part_funcs_on_sa;
          let score = part_funcs.part_func_4_align + if !is_forward && is_begin {0.} else {
            if is_forward {
              MATCH_2_INSERT_SCORE
            } else {
              INSERT_2_MATCH_SCORE
            }
          };
          sumormax(&mut sum_on_sa, score, is_viterbi);
          let score = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
          sumormax(&mut sum_on_sa, score, is_viterbi);
          let score = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
          sumormax(&mut sum_on_sa, score, is_viterbi);
        }, None => {},
      }
      let sum_4_ml = sum_4_ml + insert_score;
      tmp_part_func_sets.part_funcs_4_ml.part_func_4_insert_2 = sum_4_ml;
      sumormax(&mut tmp_sum, sum_4_ml, is_viterbi);
      let sum_4_first_bpas_on_ml = sum_4_first_bpas_on_mls + insert_score;
      tmp_part_func_sets.part_funcs_4_first_bpas_on_mls.part_func_4_insert_2 = sum_4_first_bpas_on_ml;
      sumormax(&mut tmp_sum, sum_4_first_bpas_on_mls, is_viterbi);
      tmp_part_func_sets.part_funcs_4_bpas_on_mls.part_func_4_insert_2 = tmp_sum;
      let sum_on_sa = sum_on_sa + insert_score;
      tmp_part_func_sets.part_funcs_on_sa.part_func_4_insert_2 = sum_on_sa;
      sumormax(&mut tmp_sum, sum_on_sa, is_viterbi);
      tmp_part_func_sets.part_funcs_on_mls.part_func_4_insert_2 = tmp_sum;
      tmp_part_func_set_mat.insert(pos_pair, tmp_part_func_sets);
    }
  }
  let mut final_sum_on_sa = NEG_INFINITY;
  let mut final_sum_4_ml = final_sum_on_sa;
  if is_forward {
    let ref part_func_sets = tmp_part_func_set_mat[&rightmost_pos_pair];
    let ref part_funcs = part_func_sets.part_funcs_on_sa;
    sumormax(&mut final_sum_on_sa, part_funcs.part_func_4_align, is_viterbi);
    sumormax(&mut final_sum_on_sa, part_funcs.part_func_4_insert, is_viterbi);
    sumormax(&mut final_sum_on_sa, part_funcs.part_func_4_insert_2, is_viterbi);
    let ref part_funcs = part_func_sets.part_funcs_4_ml;
    sumormax(&mut final_sum_4_ml, part_funcs.part_func_4_align, is_viterbi);
    sumormax(&mut final_sum_4_ml, part_funcs.part_func_4_insert, is_viterbi);
    sumormax(&mut final_sum_4_ml, part_funcs.part_func_4_insert_2, is_viterbi);
  }
  (tmp_part_func_set_mat, final_sum_on_sa, final_sum_4_ml)
}

pub fn get_sta_prob_mats(seq_pair: &SeqPair, seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_bp_span_pair: &PosPair, max_gap_num: Pos, max_gap_num_4_il: Pos, sta_part_func_mats: &StaPartFuncMats, bpp_mat_pair: &ProbMatPair, ss_free_energy_mat_set_pair: &SsFreeEnergyMatSetPair, uses_bpps: bool, produces_access_probs: bool, global_part_func: PartFunc) -> StaProbMats {
  let is_viterbi = false;
  let seq_len_pair = (seq_len_pair.0 as Pos, seq_len_pair.1 as Pos);
  let pseudo_pos_quadruple = (0, seq_len_pair.0 - 1, 0, seq_len_pair.1 - 1);
  let mut sta_outside_part_func_4d_mat_4_bpas = PartFunc4dMat::default();
  let mut sta_prob_mats = StaProbMats::new(&seq_len_pair);
  for substr_len_1 in (MIN_SPAN_OF_INDEX_PAIR_CLOSING_HL as Pos .. max_bp_span_pair.0 + 1).rev() {
    for i in 1 .. seq_len_pair.0 - substr_len_1 {
      let j = i + substr_len_1 - 1;
      let (long_i, long_j) = (i as usize, j as usize);
      let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
      if !is_canonical(&base_pair) {continue;}
      let pos_pair = (i, j);
      if !bpp_mat_pair.0.contains_key(&pos_pair) {continue;}
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
          let pos_pair_2 = (k, l);
          let pos_quadruple = (i, j, k, l);
          match sta_part_func_mats.part_func_4d_mat_4_bpas.get(&pos_quadruple) {
            Some(&part_func_4_bpa) => {
              let prob_coeff = part_func_4_bpa - global_part_func;
              let mut sum = NEG_INFINITY;
              let mut forward_term = sum;
              let mut backward_term = sum;
              let is_begin = i == 1 && k == 1;
              let is_end = j == seq_len_pair.0 - 2 && l == seq_len_pair.1 - 2;
              match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&(i - 1, k - 1)) {
                Some(part_funcs) => {
                  let term = part_funcs.part_func_4_align + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
                  sumormax(&mut forward_term, term, is_viterbi);
                  let term = part_funcs.part_func_4_insert + if is_begin {INIT_MATCH_SCORE} else {INSERT_2_MATCH_SCORE};
                  sumormax(&mut forward_term, term, is_viterbi);
                  let term = part_funcs.part_func_4_insert_2 + if is_begin {INIT_MATCH_SCORE} else {INSERT_2_MATCH_SCORE};
                  sumormax(&mut forward_term, term, is_viterbi);
                }, None => {},
              }
              match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&(j + 1, l + 1)) {
                Some(part_funcs) => {
                  let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                  sumormax(&mut backward_term, term, is_viterbi);
                  let term = part_funcs.part_func_4_insert + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                  sumormax(&mut backward_term, term, is_viterbi);
                  let term = part_funcs.part_func_4_insert_2 + if is_end {0.} else {MATCH_2_INSERT_SCORE};
                  sumormax(&mut backward_term, term, is_viterbi);
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
                      sumormax(bpp_4_el, bpap_4_el, is_viterbi);
                    }, None => {
                      sta_prob_mats.bpp_mat_pair_4_el.0.insert(pos_pair, bpap_4_el);
                    },
                  }
                  match sta_prob_mats.bpp_mat_pair_4_el.1.get_mut(&pos_pair_2) {
                    Some(bpp_4_el) => {
                      sumormax(bpp_4_el, bpap_4_el, is_viterbi);
                    }, None => {
                      sta_prob_mats.bpp_mat_pair_4_el.1.insert(pos_pair_2, bpap_4_el);
                    },
                  }
                }
              }
              for m in 1 .. i {
                let long_m = m as usize;
                for n in j + 1 .. seq_len_pair.0 - 1 {
                  let long_n = n as usize;
                  let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                  if !is_canonical(&base_pair_3) {continue;}
                  let pos_pair_3 = (m, n);
                  if long_n - long_j - 1 + long_i - long_m - 1 > MAX_2_LOOP_LEN {continue;}
                  if !bpp_mat_pair.0.contains_key(&pos_pair_3) {continue;}
                  let twoloop_fe = if uses_bpps {0.} else {ss_free_energy_mat_set_pair.0.twoloop_fe_4d_mat[&(m, n, i, j)]};
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
                          let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_quadruples[&pos_quadruple_2];
                          let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_quadruples[&pos_quadruple_2];
                          let mut forward_term = NEG_INFINITY;
                          let mut backward_term = forward_term;
                          match forward_tmp_part_func_set_mat.get(&(i - 1, k - 1)) {
                            Some(part_func_sets) => {
                              let ref part_funcs = part_func_sets.part_funcs_on_sa;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                            }, None => {},
                          }
                          match backward_tmp_part_func_set_mat.get(&(j + 1, l + 1)) {
                            Some(part_func_sets) => {
                              let ref part_funcs = part_func_sets.part_funcs_on_sa;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                            }, None => {},
                          }
                          let part_func_4_2l = forward_term + backward_term;
                          if part_func_4_2l > NEG_INFINITY {
                            let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                            let twoloop_fe_2 = if uses_bpps {0.} else {ss_free_energy_mat_set_pair.1.twoloop_fe_4d_mat[&(o, p, k, l)]};
                            let coefficient = bpa_score + twoloop_fe + twoloop_fe_2 + part_func;
                            let part_func_4_2l = coefficient + part_func_4_2l;
                            sumormax(&mut sum, part_func_4_2l, is_viterbi);
                            if produces_access_probs {
                              let bpap_4_2l = prob_coeff + part_func_4_2l;
                              match sta_prob_mats.access_bpp_mat_pair_4_2l.0.get_mut(&pos_pair) {
                                Some(bpp_4_2l) => {
                                  sumormax(bpp_4_2l, bpap_4_2l, is_viterbi);
                                }, None => {
                                  sta_prob_mats.access_bpp_mat_pair_4_2l.0.insert(pos_pair, bpap_4_2l);
                                },
                              }
                              match sta_prob_mats.access_bpp_mat_pair_4_2l.1.get_mut(&pos_pair_2) {
                                Some(bpp_4_2l) => {
                                  sumormax(bpp_4_2l, bpap_4_2l, is_viterbi);
                                }, None => {
                                  sta_prob_mats.access_bpp_mat_pair_4_2l.1.insert(pos_pair_2, bpap_4_2l);
                                },
                              }
                              for q in long_m + 1 .. long_i {
                                sumormax(&mut sta_prob_mats.upp_mat_pair_4_2l.0[q], bpap_4_2l, is_viterbi);
                              }
                              for q in long_j + 1 .. long_n {
                                sumormax(&mut sta_prob_mats.upp_mat_pair_4_2l.0[q], bpap_4_2l, is_viterbi);
                              }
                              for q in long_o + 1 .. long_k {
                                sumormax(&mut sta_prob_mats.upp_mat_pair_4_2l.1[q], bpap_4_2l, is_viterbi);
                              }
                              for q in long_l + 1 .. long_p {
                                sumormax(&mut sta_prob_mats.upp_mat_pair_4_2l.1[q], bpap_4_2l, is_viterbi);
                              }
                            }
                          }
                        }, None => {},
                      }
                    }
                  }
                }
              }
              let part_func_ratio = sta_part_func_mats.part_func_4d_mat_4_bpas_accessible_on_mls[&pos_quadruple] - part_func_4_bpa;
              for m in 1 .. i {
                let long_m = m as usize;
                for n in j + 1 .. seq_len_pair.0 - 1 {
                  let long_n = n as usize;
                  let base_pair_3 = (seq_pair.0[long_m], seq_pair.0[long_n]);
                  if !is_canonical(&base_pair_3) {continue;}
                  let pos_pair_3 = (m, n);
                  if !bpp_mat_pair.0.contains_key(&pos_pair_3) {continue;}
                  let base_pair = (seq_pair.0[long_m], seq_pair.0[long_n]);
                  let invert_base_pair = invert_bp(&base_pair);
                  let invert_stacking_bp = invert_bp(&(seq_pair.0[long_m + 1], seq_pair.0[long_n - 1]));
                  let ml_tm_delta_fe = if uses_bpps {0.} else {
                    ML_TM_DELTA_FES[invert_base_pair.0][invert_base_pair.1][invert_stacking_bp.0][invert_stacking_bp.1]
                  };
                  let au_or_gu_end_penalty_delta_fe = if uses_bpps {0.} else {
                    if is_au_or_gu(&base_pair) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
                  };
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
                          let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_quadruples[&pos_quadruple_2];
                          let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_quadruples[&pos_quadruple_2];
                          let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple_2];
                          let ml_tm_delta_fe_2 = if uses_bpps {0.} else {
                            ML_TM_DELTA_FES[invert_base_pair_2.0][invert_base_pair_2.1][invert_stacking_bp_2.0][invert_stacking_bp_2.1]
                          };
                          let au_or_gu_end_penalty_delta_fe_2 = if uses_bpps {0.} else {
                            if is_au_or_gu(&base_pair_2) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
                          };
                          let mut forward_term = NEG_INFINITY;
                          let mut forward_term_2 = forward_term;
                          let mut backward_term = forward_term;
                          let mut backward_term_2 = forward_term;
                          match forward_tmp_part_func_set_mat.get(&(i - 1, k - 1)) {
                            Some(part_func_sets) => {
                              let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term, term, is_viterbi);
                              let ref part_funcs = part_func_sets.part_funcs_on_sa;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut forward_term_2, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term_2, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                              sumormax(&mut forward_term_2, term, is_viterbi);
                            }, None => {},
                          }
                          match backward_tmp_part_func_set_mat.get(&(j + 1, l + 1)) {
                            Some(part_func_sets) => {
                              let ref part_funcs = part_func_sets.part_funcs_on_mls;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term, term, is_viterbi);
                              let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                              let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                              sumormax(&mut backward_term_2, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term_2, term, is_viterbi);
                              let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                              sumormax(&mut backward_term_2, term, is_viterbi);
                            }, None => {},
                          }
                          let mut part_func_4_ml = forward_term + backward_term;
                          sumormax(&mut part_func_4_ml, forward_term_2 + backward_term_2, is_viterbi);
                          if part_func_4_ml > NEG_INFINITY {
                            let coefficient = part_func_ratio + bpa_score + if uses_bpps {0.} else {2. * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe + ml_tm_delta_fe_2 + au_or_gu_end_penalty_delta_fe + au_or_gu_end_penalty_delta_fe_2} + part_func_4_bpa_2;
                            let part_func_4_ml = coefficient + part_func_4_ml;
                            sumormax(&mut sum, part_func_4_ml, is_viterbi);
                            if produces_access_probs {
                              let bpap_4_ml = prob_coeff + part_func_4_ml;
                              match sta_prob_mats.access_bpp_mat_pair_4_ml.0.get_mut(&pos_pair) {
                                Some(bpp_4_ml) => {
                                  sumormax(bpp_4_ml, bpap_4_ml, is_viterbi);
                                }, None => {
                                  sta_prob_mats.access_bpp_mat_pair_4_ml.0.insert(pos_pair, bpap_4_ml);
                                },
                              }
                              match sta_prob_mats.access_bpp_mat_pair_4_ml.1.get_mut(&pos_pair_2) {
                                Some(bpp_4_ml) => {
                                  sumormax(bpp_4_ml, bpap_4_ml, is_viterbi);
                                }, None => {
                                  sta_prob_mats.access_bpp_mat_pair_4_ml.1.insert(pos_pair_2, bpap_4_ml);
                                },
                              }
                            }
                          }
                        }, None => {},
                      }
                    }
                  }
                }
              }
              if sum > NEG_INFINITY {
                sta_outside_part_func_4d_mat_4_bpas.insert(pos_quadruple, sum);
                let bpap = (prob_coeff + sum).exp();
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
  }
  if produces_access_probs {
    let leftmost_pos_pair = (0, 0);
    let rightmost_pos_pair = (seq_len_pair.0 - 1, seq_len_pair.1 - 1);
    for u in 0 .. seq_len_pair.0 - 1 {
      let long_u = u as usize;
      let insert_score = sta_fe_params.insert_scores[long_u];
      for v in 0 .. seq_len_pair.1 - 1 {
        if u == 0 && v == 0 {continue;}
        let pos_pair = (u, v);
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num) {continue;}
        let long_v = v as usize;
        let insert_score_2 = sta_fe_params.insert_scores_2[long_v];
        let pos_pair_4_ba = (u - 1, v - 1);
        let pos_pair_4_gap_1 = (u - 1, v);
        let pos_pair_4_gap_2 = (u, v - 1);
        let pos_pair_2 = (u + 1, v + 1);
        let is_end = pos_pair_2 == rightmost_pos_pair;
        let mut backward_term_4_align = NEG_INFINITY;
        let mut backward_term_4_insert = backward_term_4_align;
        let mut backward_term_4_insert_2 = backward_term_4_align;
        match sta_part_func_mats.backward_part_func_set_mat_4_external_loop.get(&pos_pair_2) {
          Some(part_funcs) => {
            if u > 0 && v > 0 {
              let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
              sumormax(&mut backward_term_4_align, term, is_viterbi);
              let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
              sumormax(&mut backward_term_4_align, term, is_viterbi);
              let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
              sumormax(&mut backward_term_4_align, term, is_viterbi);
            }
            if u > 0 {
              let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
              sumormax(&mut backward_term_4_insert, term, is_viterbi);
              let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
              sumormax(&mut backward_term_4_insert, term, is_viterbi);
              let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
              sumormax(&mut backward_term_4_insert, term, is_viterbi);
            }
            if v > 0 {
              let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
              sumormax(&mut backward_term_4_insert_2, term, is_viterbi);
              let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
              sumormax(&mut backward_term_4_insert_2, term, is_viterbi);
              let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
              sumormax(&mut backward_term_4_insert_2, term, is_viterbi);
            }
          }, None => {},
        }
        if u > 0 && v > 0 {
          match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_4_ba) {
            Some(part_funcs) => {
              let mut forward_term = NEG_INFINITY;
              let is_begin = pos_pair_4_ba == leftmost_pos_pair;
              let term = part_funcs.part_func_4_align + if is_begin {INIT_MATCH_SCORE} else {MATCH_2_MATCH_SCORE};
              sumormax(&mut forward_term, term, is_viterbi);
              let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
              sumormax(&mut forward_term, term, is_viterbi);
              let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
              sumormax(&mut forward_term, term, is_viterbi);
              let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
              let bap_4_el = ba_score + forward_term + backward_term_4_align - global_part_func;
              sumormax(&mut sta_prob_mats.upp_mat_pair_4_el.0[long_u], bap_4_el, is_viterbi);
              sumormax(&mut sta_prob_mats.upp_mat_pair_4_el.1[long_v], bap_4_el, is_viterbi);
            }, None => {},
          }
        }
        if u > 0 {
          match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_4_gap_1) {
            Some(part_funcs) => {
              let mut forward_term = NEG_INFINITY;
              let is_begin = pos_pair_4_gap_1 == leftmost_pos_pair;
              let term = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
              sumormax(&mut forward_term, term, is_viterbi);
              let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
              sumormax(&mut forward_term, term, is_viterbi);
              let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
              sumormax(&mut forward_term, term, is_viterbi);
              let upp_4_el = insert_score + forward_term + backward_term_4_insert - global_part_func;
              sumormax(&mut sta_prob_mats.upp_mat_pair_4_el.0[long_u], upp_4_el, is_viterbi);
            }, None => {},
          }
        }
        if v > 0 {
          match sta_part_func_mats.forward_part_func_set_mat_4_external_loop.get(&pos_pair_4_gap_2) {
            Some(part_funcs) => {
              let mut forward_term = NEG_INFINITY;
              let is_begin = pos_pair_4_gap_2 == leftmost_pos_pair;
              let term = part_funcs.part_func_4_align + if is_begin {INIT_INSERT_SCORE} else {MATCH_2_INSERT_SCORE};
              sumormax(&mut forward_term, term, is_viterbi);
              let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
              sumormax(&mut forward_term, term, is_viterbi);
              let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
              sumormax(&mut forward_term, term, is_viterbi);
              let upp_4_el = insert_score_2 + forward_term + backward_term_4_insert_2 - global_part_func;
              sumormax(&mut sta_prob_mats.upp_mat_pair_4_el.1[long_v], upp_4_el, is_viterbi);
            }, None => {},
          }
        }
        if !is_min_gap_ok_1(&pos_pair, &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
        if !(u > 0 && v > 0) {continue;}
        let ba_score = sta_fe_params.ba_score_mat[&pos_pair];
        for i in 1 .. u {
          for j in u + 1 .. seq_len_pair.0 - 1 {
            let (long_i, long_j) = (i as usize, j as usize);
            let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
            if !is_canonical(&base_pair) {continue;}
            if !bpp_mat_pair.0.contains_key(&(i, j)) {continue;}
            let hl_fe = if uses_bpps {0.} else {ss_free_energy_mat_set_pair.0.hl_fe_mat[&(i, j)]};
            let invert_base_pair = invert_bp(&base_pair);
            let invert_stacking_bp = invert_bp(&(seq_pair.0[long_i + 1], seq_pair.0[long_j - 1]));
            let ml_tm_delta_fe = if uses_bpps {0.} else {
              ML_TM_DELTA_FES[invert_base_pair.0][invert_base_pair.1][invert_stacking_bp.0][invert_stacking_bp.1]
            };
            let au_or_gu_end_penalty_delta_fe = if uses_bpps {0.} else{
              if is_au_or_gu(&base_pair) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
            };
            for k in 1 .. v {
              if !is_min_gap_ok_1(&(i, k), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
              for l in v + 1 .. seq_len_pair.1 - 1 {
                let (long_k, long_l) = (k as usize, l as usize);
                if !is_min_gap_ok_1(&(j, l), &pseudo_pos_quadruple, max_gap_num_4_il) {continue;}
                let base_pair_2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
                if !is_canonical(&base_pair_2) {continue;}
                if !bpp_mat_pair.1.contains_key(&(k, l)) {continue;}
                let hl_fe_2 = if uses_bpps {0.} else {ss_free_energy_mat_set_pair.1.hl_fe_mat[&(k, l)]};
                let invert_base_pair_2 = invert_bp(&base_pair_2);
                let invert_stacking_bp_2 = invert_bp(&(seq_pair.1[long_k + 1], seq_pair.1[long_l - 1]));
                let ml_tm_delta_fe_2 = if uses_bpps {0.} else {
                  ML_TM_DELTA_FES[invert_base_pair_2.0][invert_base_pair_2.1][invert_stacking_bp_2.0][invert_stacking_bp_2.1]
                };
                let au_or_gu_end_penalty_delta_fe_2 = if uses_bpps {0.} else {
                  if is_au_or_gu(&base_pair_2) {HELIX_AU_OR_GU_END_PENALTY_DELTA_FE} else {0.}
                };
                let pos_quadruple = (i, j, k, l);
                match sta_outside_part_func_4d_mat_4_bpas.get(&pos_quadruple) {
                  Some(&part_func_4_bpa) => {
                    let bpa_score = sta_fe_params.bpa_score_mat[&pos_quadruple];
                    let prob_coeff = part_func_4_bpa - global_part_func + bpa_score;
                    let ref forward_tmp_part_func_set_mat = sta_part_func_mats.forward_tmp_part_func_set_mats_with_pos_quadruples[&pos_quadruple];
                    let ref backward_tmp_part_func_set_mat = sta_part_func_mats.backward_tmp_part_func_set_mats_with_pos_quadruples[&pos_quadruple];
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
                        sumormax(&mut backward_term_4_align_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_2_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_2_on_sa, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_2_on_sa, term, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_ml;
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_align_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_2_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_2_4_ml, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_2_4_ml, term, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_bpas_on_mls;
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_align_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_2_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_2_4_bpas_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_2_4_bpas_on_mls, term, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_on_mls;
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {MATCH_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_align_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + MATCH_2_INSERT_SCORE;
                        sumormax(&mut backward_term_4_align_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_align + if is_end {0.} else {INSERT_2_MATCH_SCORE};
                        sumormax(&mut backward_term_4_insert_2_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                        sumormax(&mut backward_term_4_insert_2_on_mls, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                        sumormax(&mut backward_term_4_insert_2_on_mls, term, is_viterbi);
                      }, None => {},
                    }
                    let prob_coeff_4_hl = prob_coeff + hl_fe + hl_fe_2;
                    let prob_coeff_4_ml = prob_coeff + if uses_bpps {0.} else {
                      2. * CONST_4_INIT_ML_DELTA_FE + ml_tm_delta_fe + ml_tm_delta_fe_2 + au_or_gu_end_penalty_delta_fe + au_or_gu_end_penalty_delta_fe_2
                    };
                    match forward_tmp_part_func_set_mat.get(&pos_pair_4_ba) {
                      Some(part_func_sets) => {
                        let ref part_funcs = part_func_sets.part_funcs_on_sa;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let bap_4_hl = prob_coeff_4_hl + ba_score + forward_term + backward_term_4_align_on_sa;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_hl.0[long_u], bap_4_hl, is_viterbi);
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_hl.1[long_v], bap_4_hl, is_viterbi);
                        let bap_4_ml = prob_coeff_4_ml + ba_score + forward_term + backward_term_4_align_4_ml;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], bap_4_ml, is_viterbi);
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], bap_4_ml, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let bap_4_ml = prob_coeff_4_ml + ba_score + forward_term + backward_term_4_align_4_bpas_on_mls;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], bap_4_ml, is_viterbi);
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], bap_4_ml, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_ml;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_2_MATCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let bap_4_ml = prob_coeff_4_ml + ba_score + forward_term + backward_term_4_align_on_mls;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], bap_4_ml, is_viterbi);
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], bap_4_ml, is_viterbi);
                      }, None => {},
                    }
                    match forward_tmp_part_func_set_mat.get(&pos_pair_4_gap_1) {
                      Some(part_func_sets) => {
                        let ref part_funcs = part_func_sets.part_funcs_on_sa;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let upp_4_hl = prob_coeff_4_hl + insert_score + forward_term + backward_term_4_insert_on_sa;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_hl.0[long_u], upp_4_hl, is_viterbi);
                        let upp_4_ml = prob_coeff_4_ml + insert_score + forward_term + backward_term_4_insert_4_ml;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], upp_4_ml, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let upp_4_ml = prob_coeff_4_ml + insert_score + forward_term + backward_term_4_insert_4_bpas_on_mls;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], upp_4_ml, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_ml;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_EXTEND_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_SWITCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let upp_4_ml = prob_coeff_4_ml + insert_score + forward_term + backward_term_4_insert_on_mls;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.0[long_u], upp_4_ml, is_viterbi);
                      }, None => {},
                    }
                    match forward_tmp_part_func_set_mat.get(&pos_pair_4_gap_2) {
                      Some(part_func_sets) => {
                        let ref part_funcs = part_func_sets.part_funcs_on_sa;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let upp_4_hl = prob_coeff_4_hl + insert_score_2 + forward_term + backward_term_4_insert_2_on_sa;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_hl.1[long_v], upp_4_hl, is_viterbi);
                        let upp_4_ml = prob_coeff_4_ml + insert_score_2 + forward_term + backward_term_4_insert_2_4_ml;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], upp_4_ml, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_first_bpas_on_mls;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let upp_4_ml = prob_coeff_4_ml + insert_score_2 + forward_term + backward_term_4_insert_2_4_bpas_on_mls;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], upp_4_ml, is_viterbi);
                        let ref part_funcs = part_func_sets.part_funcs_4_ml;
                        let mut forward_term = NEG_INFINITY;
                        let term = part_funcs.part_func_4_align + MATCH_2_INSERT_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert + INSERT_SWITCH_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let term = part_funcs.part_func_4_insert_2 + INSERT_EXTEND_SCORE;
                        sumormax(&mut forward_term, term, is_viterbi);
                        let upp_4_ml = prob_coeff_4_ml + insert_score_2 + forward_term + backward_term_4_insert_2_on_mls;
                        sumormax(&mut sta_prob_mats.upp_mat_pair_4_ml.1[long_v], upp_4_ml, is_viterbi);
                      }, None => {},
                    }
                  }, None => {},
                }
              }
            }
          }
        }
      }
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_2l.0.values_mut() {
      *bpp = bpp.exp()
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_2l.1.values_mut() {
      *bpp = bpp.exp()
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_ml.0.values_mut() {
      *bpp = bpp.exp()
    }
    for bpp in sta_prob_mats.access_bpp_mat_pair_4_ml.1.values_mut() {
      *bpp = bpp.exp()
    }
    for bpp in sta_prob_mats.bpp_mat_pair_4_el.0.values_mut() {
      *bpp = bpp.exp()
    }
    for bpp in sta_prob_mats.bpp_mat_pair_4_el.1.values_mut() {
      *bpp = bpp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_hl.0.iter_mut() {
      *upp = upp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_hl.1.iter_mut() {
      *upp = upp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_ml.0.iter_mut() {
      *upp = upp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_ml.1.iter_mut() {
      *upp = upp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_2l.0.iter_mut() {
      *upp = upp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_2l.1.iter_mut() {
      *upp = upp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_el.0.iter_mut() {
      *upp = upp.exp()
    }
    for upp in sta_prob_mats.upp_mat_pair_4_el.1.iter_mut() {
      *upp = upp.exp()
    }
  }
  sta_prob_mats
}

pub fn pct_of_prob_mats(prob_mats_with_rna_id_pairs: &StaProbMatsWithRnaIdPairs, rna_id: RnaId, num_of_rnas: usize, bpp_mats: &SsProbMats, upp_mat_len: usize, produces_access_probs: bool) -> PctStaProbMats {
  let weight = 1. / (num_of_rnas - 1) as Prob;
  let mut pct_prob_mats = PctStaProbMats::new(upp_mat_len);
  pct_prob_mats.bpp_mat_on_ss = bpp_mats.bpp_mat.clone();
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
  for (&(i, j), &bpp) in &bpp_mats.bpp_mat {
    let pos_pair = (i + 1, j + 1);
    if !pct_prob_mats.bpp_mat.contains_key(&pos_pair) {
      pct_prob_mats.bpp_mat.insert(pos_pair, bpp);
    }
  }
  if produces_access_probs {
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
    for (&(i, j), &bpp) in &bpp_mats.access_bpp_mat_4_2l {
      let pos_pair = (i + 1, j + 1);
      if !pct_prob_mats.access_bpp_mat_4_2l.contains_key(&pos_pair) {
        pct_prob_mats.access_bpp_mat_4_2l.insert(pos_pair, bpp);
      }
    }
    for (&(i, j), &bpp) in &bpp_mats.access_bpp_mat_4_ml {
      let pos_pair = (i + 1, j + 1);
      if !pct_prob_mats.access_bpp_mat_4_ml.contains_key(&pos_pair) {
        pct_prob_mats.access_bpp_mat_4_ml.insert(pos_pair, bpp);
      }
    }
    for (&(i, j), &bpp) in &bpp_mats.bpp_mat_4_el {
      let pos_pair = (i + 1, j + 1);
      if !pct_prob_mats.bpp_mat_4_el.contains_key(&pos_pair) {
        pct_prob_mats.bpp_mat_4_el.insert(pos_pair, bpp);
      }
    }
  }
  pct_prob_mats
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

pub fn remove_small_bpps_from_bpp_mat(sparse_bpp_mat: &SparseProbMat, min_bpp: Prob) -> SparseProbMat {
  sparse_bpp_mat.iter().filter(|(_, &bpp)| {bpp >= min_bpp}).map(|(&(i, j), &bpp)| {((i + 1, j + 1), bpp)}).collect()
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

pub fn sumormax(sum: &mut FreeEnergy, new_term: FreeEnergy, is_viterbi: bool) {
  if new_term > NEG_INFINITY {
    *sum = if !sum.is_finite() {
     new_term
    } else {
      let max = sum.max(new_term);
      if is_viterbi {
        max
      } else {
        max + ((if *sum == max {new_term - max} else {*sum - max}).exp() + 1.).ln()
      }
    };
  }
}

pub fn consprob(thread_pool: &mut Pool, fasta_records: &FastaRecords, min_bpp: Prob, offset_4_max_gap_num: Pos, uses_bpps: bool, produces_access_probs: bool) -> ProbMatSets {
  let num_of_fasta_records = fasta_records.len();
  let mut bpp_mat_sets = vec![SsProbMats::new(); num_of_fasta_records];
  let mut sparse_bpp_mats = vec![SparseProbMat::new(); num_of_fasta_records];
  let mut max_bp_spans = vec![0; num_of_fasta_records];
  let mut ss_free_energy_mat_sets = vec![SsFreeEnergyMats::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (bpp_mats, sparse_bpp_mat, max_bp_span, fasta_record, ss_free_energy_mats) in multizip((bpp_mat_sets.iter_mut(), sparse_bpp_mats.iter_mut(), max_bp_spans.iter_mut(), fasta_records.iter(), ss_free_energy_mat_sets.iter_mut())) {
      let seq_len = fasta_record.seq.len();
      scope.execute(move || {
        let (obtained_bpp_mats, obtained_ss_free_energy_mats) = mccaskill_algo(&fasta_record.seq[1 .. seq_len - 1]);
        *bpp_mats = obtained_bpp_mats;
        *ss_free_energy_mats = obtained_ss_free_energy_mats;
        ss_free_energy_mats.sparsify(&bpp_mats.bpp_mat, min_bpp);
        *sparse_bpp_mat = remove_small_bpps_from_bpp_mat(&bpp_mats.bpp_mat, min_bpp);
        *max_bp_span = get_max_bp_span(sparse_bpp_mat);
      });
    }
  });
  let mut prob_mats_with_rna_id_pairs = StaProbMatsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      prob_mats_with_rna_id_pairs.insert(rna_id_pair, StaProbMats::origin());
    }
  }
  thread_pool.scoped(|scope| {
    for (rna_id_pair, prob_mats) in prob_mats_with_rna_id_pairs.iter_mut() {
      let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
      let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
      let max_gap_num = offset_4_max_gap_num + (max(seq_len_pair.0, seq_len_pair.1) - min(seq_len_pair.0, seq_len_pair.1)) as Pos;
      let max_gap_num_4_il = max(min(max_gap_num, MAX_GAP_NUM_4_IL), MIN_GAP_NUM_4_IL);
      let max_bp_span_pair = (max_bp_spans[rna_id_pair.0], max_bp_spans[rna_id_pair.1]);
      let bpp_mat_pair = (&sparse_bpp_mats[rna_id_pair.0], &sparse_bpp_mats[rna_id_pair.1]);
      let ss_free_energy_mat_set_pair = (&ss_free_energy_mat_sets[rna_id_pair.0], &ss_free_energy_mat_sets[rna_id_pair.1]);
      scope.execute(move || {
        let sta_fe_params = StaFeParams::new(&seq_pair, &seq_len_pair, &max_bp_span_pair, max_gap_num, max_gap_num_4_il, &bpp_mat_pair, uses_bpps);
        *prob_mats = io_algo_4_prob_mats(&seq_pair, &seq_len_pair, &sta_fe_params, &max_bp_span_pair, max_gap_num, max_gap_num_4_il, &bpp_mat_pair, &ss_free_energy_mat_set_pair, uses_bpps, produces_access_probs);
      });
    }
  });
  let mut prob_mat_sets = vec![PctStaProbMats::origin(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (rna_id, prob_mats, bpp_mats) in multizip((0 .. num_of_fasta_records, prob_mat_sets.iter_mut(), bpp_mat_sets.iter_mut())) {
      let ref ref_2_prob_mats_with_rna_id_pairs = prob_mats_with_rna_id_pairs;
      let seq_len = fasta_records[rna_id].seq.len();
      scope.execute(move || {
        *prob_mats = pct_of_prob_mats(ref_2_prob_mats_with_rna_id_pairs, rna_id, num_of_fasta_records, bpp_mats, seq_len, produces_access_probs);
      });
    }
  });
  prob_mat_sets
}

pub fn write_prob_mat_sets(output_dir_path: &Path, prob_mat_sets: &ProbMatSets, produces_access_probs: bool) {
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let bpp_mat_file_path = output_dir_path.join(BPP_MAT_FILE_NAME);
  let mut writer_2_bpp_mat_file = BufWriter::new(File::create(bpp_mat_file_path).unwrap());
  let mut buf_4_writer_2_bpp_mat_file = String::new();
  for (rna_id, prob_mats) in prob_mat_sets.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in prob_mats.bpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
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
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
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
        buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
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
        buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
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
        buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
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
