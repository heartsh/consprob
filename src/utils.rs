pub use rna_algos::utils::*;
pub use rna_algos::mccaskill_algo::*;
use std::f32::consts::LOG2_E;

pub type BaScoreMat = HashMap<BasePair, FreeEnergy>;
pub type BpaScoreMat = HashMap<(BasePair, BasePair), FreeEnergy>;
pub type InsertScores = [FreeEnergy; NUM_OF_BASES];

pub const PSEUDO_BASE: Base = U + 1 as Base;
pub const INSERT_SCORES: InsertScores = [-0.002521927159, -0.08313891561, -0.07443970653, -0.01290054598];
pub const INIT_MATCH_SCORE: FreeEnergy = 0.3959924457;
pub const INIT_INSERT_SCORE: FreeEnergy = -0.3488104904;
pub const MATCH_2_MATCH_SCORE: FreeEnergy = 2.50575671;
pub const MATCH_2_INSERT_SCORE: FreeEnergy = 0.1970448791;
pub const INSERT_2_MATCH_SCORE: FreeEnergy = MATCH_2_INSERT_SCORE;
pub const INSERT_EXTEND_SCORE: FreeEnergy = 1.014026583;
pub const INSERT_SWITCH_SCORE: FreeEnergy = -7.346968782;
lazy_static! {
  pub static ref BA_SCORE_MAT: BaScoreMat = {
    [
      (AA, 0.5256508867), (AC, -0.40906402), (AG, -0.2502759109), (AU, -0.3252306723),
      (CA, -0.40906402), (CC, 0.6665219366), (CG, -0.3289391181), (CU, -0.1326088918),
      (GA, -0.2502759109), (GC, -0.3289391181), (GG, 0.6684676551), (GU, -0.3565888168),
      (UA, -0.3252306723), (UC, -0.1326088918), (UG, -0.3565888168), (UU, 0.459052045),
    ].iter().map(|(base_pair, ba_score)| {(*base_pair, *ba_score)}).collect()
  };
  pub static ref RIBOSUM_BA_SCORE_MAT: BaScoreMat = {
    [
      (AA, 2.22), (AC, -1.86), (AG, -1.46), (AU, -1.39),
      (CA, -1.86), (CC, 1.16), (CG, -2.48), (CU, -1.05),
      (GA, -1.46), (GC, -2.48), (GG, 1.03), (GU, -1.74),
      (UA, -1.39), (UC, -1.05), (UG, -1.74), (UU, 1.65),
    ].iter().map(|(base_pair, ba_score)| {(*base_pair, ba_score / LOG2_E)}).collect()
  };
  pub static ref RIBOSUM_BPA_SCORE_MAT: BpaScoreMat = {
    [
      ((AU, AU), 4.49), ((AU, CG), 1.67), ((AU, GC), 2.70), ((AU, GU), 0.59), ((AU, UA), 1.61), ((AU, UG), -0.51),
      ((CG, AU), 1.67), ((CG, CG), 5.36), ((CG, GC), 2.11), ((CG, GU), -0.27), ((CG, UA), 2.75), ((CG, UG), 1.32),
      ((GC, AU), 2.70), ((GC, CG), 2.11), ((GC, GC), 5.62), ((GC, GU), 1.21), ((GC, UA), 1.6), ((GC, UG), -0.08),
      ((GU, AU), 0.59), ((GU, CG), -0.27), ((GU, GC), 1.21), ((GU, GU), 3.47), ((GU, UA), -0.57), ((GU, UG), -2.09),
      ((UA, AU), 1.61), ((UA, CG), 2.75), ((UA, GC), 1.6), ((UA, GU), -0.57), ((UA, UA), 4.97), ((UA, UG), 1.14),
      ((UG, AU), -0.51), ((UG, CG), 1.32), ((UG, GC), -0.08), ((UG, GU), -2.09), ((UG, UA), 1.14), ((UG, UG), 3.36),
    ].iter().map(|(base_quadruple, bpa_score)| {(*base_quadruple, bpa_score / LOG2_E)}).collect()
  };
}
