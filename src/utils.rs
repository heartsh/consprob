pub use rna_algos::utils::*;
pub use rna_algos::mccaskill_algo::*;
use std::f64::consts::LOG2_E;

pub type BaScoreMat = FxHashMap<BasePair, FreeEnergy>;
pub type BpaScoreMat = FxHashMap<(BasePair, BasePair), FreeEnergy>;

// pub const SEQ_ALPHABET: [Base; NUM_OF_BASES] = [A, C, G, U];
pub const PSEUDO_BASE: Base = U + 1 as Base;
lazy_static! {
  pub static ref RIBOSUM_85_60_BA_SCORE_MAT: BaScoreMat = {
    [
      (AA, 2.22), (AC, -1.86), (AG, -1.46), (AU, -1.39),
      (CA, -1.86), (CC, 1.16), (CG, -2.48), (CU, -1.05),
      (GA, -1.46), (GC, -2.48), (GG, 1.03), (GU, -1.74),
      (UA, -1.39), (UC, -1.05), (UG, -1.74), (UU, 1.65),
    ].iter().map(|(base_pair, ba_score)| {(*base_pair, ba_score / LOG2_E)}).collect()
  };
  pub static ref EXP_RIBOSUM_85_60_BA_SCORE_MAT: BaScoreMat = {RIBOSUM_85_60_BA_SCORE_MAT.iter().map(|(base_pair, &ba_score)| {(*base_pair, ba_score.exp())}).collect()};
  pub static ref RIBOSUM_85_60_BPA_SCORE_MAT: BpaScoreMat = {
    [
      ((AU, AU), 4.49), ((AU, CG), 1.67), ((AU, GC), 2.70), ((AU, GU), 0.59), ((AU, UA), 1.61), ((AU, UG), -0.51),
      ((CG, AU), 1.67), ((CG, CG), 5.36), ((CG, GC), 2.11), ((CG, GU), -0.27), ((CG, UA), 2.75), ((CG, UG), 1.32),
      ((GC, AU), 2.70), ((GC, CG), 2.11), ((GC, GC), 5.62), ((GC, GU), 1.21), ((GC, UA), 1.6), ((GC, UG), -0.08),
      ((GU, AU), 0.59), ((GU, CG), -0.27), ((GU, GC), 1.21), ((GU, GU), 3.47), ((GU, UA), -0.57), ((GU, UG), -2.09),
      ((UA, AU), 1.61), ((UA, CG), 2.75), ((UA, GC), 1.6), ((UA, GU), -0.57), ((UA, UA), 4.97), ((UA, UG), 1.14),
      ((UG, AU), -0.51), ((UG, CG), 1.32), ((UG, GC), -0.08), ((UG, GU), -2.09), ((UG, UA), 1.14), ((UG, UG), 3.36),
    ].iter().map(|(base_quadruple, bpa_score)| {(*base_quadruple, bpa_score / LOG2_E)}).collect()
  };
  pub static ref EXP_RIBOSUM_85_60_BPA_SCORE_MAT: BpaScoreMat = {RIBOSUM_85_60_BPA_SCORE_MAT.iter().map(|(base_quadruple, &bpa_score)| {(*base_quadruple, bpa_score.exp())}).collect()};
}
