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
}
