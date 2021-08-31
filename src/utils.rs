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
      ((AA, AA), -2.49), ((AA, AC), -7.04), ((AA, AG), -8.24), ((AA, AU), -4.32), ((AA, CA), -8.84), ((AA, CC), -14.37), ((AA, CG), -4.68), ((AA, CU), -12.64), ((AA, GA), -6.86), ((AA, GC), -5.03), ((AA, GG), -8.39), ((AA, GU), -5.84), ((AA, UA), -4.01), ((AA, UC), -11.32), ((AA, UG), -6.16), ((AA, UU), -9.05),
      ((AC, AA), -7.04), ((AC, AC), -2.11), ((AC, AG), -8.89), ((AC, AU), -2.04), ((AC, CA), -9.37), ((AC, CC), -9.08), ((AC, CG), -5.86), ((AC, CU), -10.45), ((AC, GA), -9.73), ((AC, GC), -3.81), ((AC, GG), -11.05), ((AC, GU), -4.72), ((AC, UA), -5.33), ((AC, UC), -8.67), ((AC, UG), -6.93), ((AC, UU), -7.83),
      ((AG, AA), -8.24), ((AG, AC), -8.89), ((AG, AG), -0.80), ((AG, AU), -5.13), ((AG, CA), -10.41), ((AG, CC), -14.53), ((AG, CG), -4.57), ((AG, CU), -10.14), ((AG, GA), -8.61), ((AG, GC), -5.77), ((AG, GG), -5.38), ((AG, GU), -6.60), ((AG, UA), -5.43), ((AG, UC), -8.87), ((AG, UG), -5.94), ((AG, UU), -11.07),
      ((AU, AA), -4.32), ((AU, AC), -2.04), ((AU, AG), -5.13), ((AU, AU), 4.49), ((AU, CA), -5.56), ((AU, CC), -6.71), ((AU, CG), 1.67), ((AU, CU), -5.17), ((AU, GA), -5.33), ((AU, GC), 2.70), ((AU, GG), -5.61), ((AU, GU), 0.59), ((AU, UA), 1.61), ((AU, UC), -4.81), ((AU, UG), -0.51), ((AU, UU), -2.98),
      ((CA, AA), -8.84), ((CA, AC), -9.37), ((CA, AG), -10.41), ((CA, AU), -5.56), ((CA, CA), -5.13), ((CA, CC), -10.45), ((CA, CG), -3.57), ((CA, CU), -8.49), ((CA, GA), -7.98), ((CA, GC), -5.95), ((CA, GG), -11.36), ((CA, GU), -7.93), ((CA, UA), -2.42), ((CA, UC), -7.08), ((CA, UG), -5.63), ((CA, UU), -8.39),
      ((CC, AA), -14.37), ((CC, AC), -9.08), ((CC, AG), -14.53), ((CC, AU), -6.71), ((CC, CA), -10.45), ((CC, CC), -3.59), ((CC, CG), -5.71), ((CC, CU), -5.77), ((CC, GA), -12.43), ((CC, GC), -3.70), ((CC, GG), -12.58), ((CC, GU), -7.88), ((CC, UA), -6.88), ((CC, UC), -7.40), ((CC, UG), -8.41), ((CC, UU), -5.41),
      ((CG, AA), -4.68), ((CG, AC), -5.86), ((CG, AG), -4.57), ((CG, AU), 1.67), ((CG, CA), -3.57), ((CG, CC), -5.71), ((CG, CG), 5.36), ((CG, CU), -4.96), ((CG, GA), -6.00), ((CG, GC), 2.11), ((CG, GG), -4.66), ((CG, GU), -0.27), ((CG, UA), 2.75), ((CG, UC), -4.91), ((CG, UG), 1.32), ((CG, UU), -3.67),
      ((CU, AA), -12.64), ((CU, AC), -10.45), ((CU, AG), -10.14), ((CU, AU), -5.17), ((CU, CA), -8.49), ((CU, CC), -5.77), ((CU, CG), -4.96), ((CU, CU), -2.28), ((CU, GA), -7.71), ((CU, GC), -5.84), ((CU, GG), -13.69), ((CU, GU), -5.61), ((CU, UA), -4.72), ((CU, UC), -3.83), ((CU, UG), -7.36), ((CU, UU), -5.21),
      ((GA, AA), -6.86), ((GA, AC), -9.73), ((GA, AG), -8.61), ((GA, AU), -5.33), ((GA, CA), -7.98), ((GA, CC), -12.43), ((GA, CG), -6.00), ((GA, CU), -7.71), ((GA, GA), -1.05), ((GA, GC), -4.88), ((GA, GG), -8.67), ((GA, GU), -6.10), ((GA, UA), -5.85), ((GA, UC), -6.63), ((GA, UG), -7.55), ((GA, UU), -11.54),
      ((GC, AA), -5.03), ((GC, AC), -3.81), ((GC, AG), -5.77), ((GC, AU), 2.70), ((GC, CA), -5.95), ((GC, CC), -3.70), ((GC, CG), 2.11), ((GC, CU), -5.84), ((GC, GA), -4.88), ((GC, GC), 5.62), ((GC, GG), -4.13), ((GC, GU), 1.21), ((GC, UA), 1.60), ((GC, UC), -4.49), ((GC, UG), -0.08), ((GC, UU), -3.90),
      ((GG, AA), -8.39), ((GG, AC), -11.05), ((GG, AG), -5.38), ((GG, AU), -5.61), ((GG, CA), -11.36), ((GG, CC), -12.58), ((GG, CG), -4.66), ((GG, CU), -13.69), ((GG, GA), -8.67), ((GG, GC), -4.13), ((GG, GG), -1.98), ((GG, GU), -5.77), ((GG, UA), -5.75), ((GG, UC), -12.01), ((GG, UG), -4.27), ((GG, UU), -10.79),
      ((GU, AA), -5.84), ((GU, AC), -4.72), ((GU, AG), -6.60), ((GU, AU), 0.59), ((GU, CA), -7.93), ((GU, CC), -7.88), ((GU, CG), -0.27), ((GU, CU), -5.61), ((GU, GA), -6.10), ((GU, GC), 1.21), ((GU, GG), -5.77), ((GU, GU), 3.47), ((GU, UA), -0.57), ((GU, UC), -5.30), ((GU, UG), -2.09), ((GU, UU), -4.45),
      ((UA, AA), -4.01), ((UA, AC), -5.33), ((UA, AG), -5.43), ((UA, AU), 1.61), ((UA, CA), -2.42), ((UA, CC), -6.88), ((UA, CG), 2.75), ((UA, CU), -4.72), ((UA, GA), -5.85), ((UA, GC), 1.60), ((UA, GG), -5.75), ((UA, GU), -0.57), ((UA, UA), 4.97), ((UA, UC), -2.98), ((UA, UG), 1.14), ((UA, UU), -3.39),
      ((UC, AA), -11.32), ((UC, AC), -8.67), ((UC, AG), -8.87), ((UC, AU), -4.81), ((UC, CA), -7.08), ((UC, CC), -7.40), ((UC, CG), -4.91), ((UC, CU), -3.83), ((UC, GA), -6.63), ((UC, GC), -4.49), ((UC, GG), -12.01), ((UC, GU), -5.30), ((UC, UA), -2.98), ((UC, UC), -3.21), ((UC, UG), -4.76), ((UC, UU), -5.97),
      ((UG, AA), -6.16), ((UG, AC), -6.93), ((UG, AG), -5.94), ((UG, AU), -0.51), ((UG, CA), -5.63), ((UG, CC), -8.41), ((UG, CG), 1.32), ((UG, CU), -7.36), ((UG, GA), -7.55), ((UG, GC), -0.08), ((UG, GG), -4.27), ((UG, GU), -2.09), ((UG, UA), 1.14), ((UG, UC), -4.76), ((UG, UG), 3.36), ((UG, UU), -4.28),
      ((UU, AA), -9.05), ((UU, AC), -7.83), ((UU, AG), -11.07), ((UU, AU), -2.98), ((UU, CA), -8.39), ((UU, CC), -5.41), ((UU, CG), -3.67), ((UU, CU), -5.21), ((UU, GA), -11.54), ((UU, GC), -3.90), ((UU, GG), -10.79), ((UU, GU), -4.45), ((UU, UA), -3.39), ((UU, UC), -5.97), ((UU, UG), -4.28), ((UU, UU), -0.02),
    ].iter().map(|(base_quadruple, bpa_score)| {(*base_quadruple, bpa_score / LOG2_E)}).collect()
  };
}
