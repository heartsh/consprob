pub use rna_algos::utils::*;
pub use rna_algos::mccaskill_algo::*;

type LbapsWithBasePairs = HashMap<BasePair, LogProb, Hasher>;
type LogGapProbsWithBases = HashMap<Base, LogProb, Hasher>;
pub type BaseQuadruple = (Base, Base, Base, Base);
type LbpapsWithBaseQuadruples = HashMap<BaseQuadruple, LogProb, Hasher>;
type LbppsWithBasePairs = HashMap<BasePair, LogProb, Hasher>;
type LnbppsWithBases = HashMap<Base, LogProb, Hasher>;
#[derive(Debug)]
pub struct StemParams {
  pub lbaps_with_base_pairs: LbapsWithBasePairs,
  pub logps_with_bases: LogGapProbsWithBases,
  pub legps_with_bases: LogGapProbsWithBases,
  pub lbpaps_with_base_quadruples: LbpapsWithBaseQuadruples,
  pub lbpps_with_base_pairs: LbppsWithBasePairs,
  pub lnbpps_with_bases: LnbppsWithBases,
}

pub const SEQ_ALPHABET: [Base; 4] = [A, C, G, U];
pub const PSEUDO_BASE: Base = '$' as Base;
lazy_static! {
  pub static ref BA_ALPHABET: Vec<BasePair> = {
    let mut ba_alphabet = Vec::new();
    for (i, &base_1) in SEQ_ALPHABET.iter().enumerate() {
      for &base_2 in &SEQ_ALPHABET[i ..] {
        ba_alphabet.push((base_1, base_2));
      }
    }
    ba_alphabet
  };
  pub static ref BP_ALPHABET: HashMap<BasePair, bool, Hasher> = {
    [AU, CG, GC, GU, UA, UG].iter().map(|base_pair| {(*base_pair, true)}).collect()
  };
  pub static ref BPA_ALPHABET: Vec<BaseQuadruple> = {
    let mut bpa_alphabet = Vec::new();
    for base_pair_1 in BP_ALPHABET.keys() {
      for base_pair_2 in BP_ALPHABET.keys() {
        if base_pair_1 <= base_pair_2 {
          bpa_alphabet.push((base_pair_1.0, base_pair_1.1, base_pair_2.0, base_pair_2.1));
        }
      }
    }
    bpa_alphabet
  };
}

impl StemParams {
  pub fn new() -> StemParams {
    let lbaps_with_base_pairs = BA_ALPHABET.iter().map(|base_pair| {(*base_pair, NEG_INFINITY)}).collect();
    let lgps_with_bases: LogGapProbsWithBases = SEQ_ALPHABET.iter().map(|&base| {(base, NEG_INFINITY)}).collect();
    let lbpaps_with_base_quadruples = BPA_ALPHABET.iter().map(|base_quadruple| {(*base_quadruple, NEG_INFINITY)}).collect();
    StemParams {
      lbaps_with_base_pairs: lbaps_with_base_pairs,
      logps_with_bases: lgps_with_bases.clone(),
      legps_with_bases: lgps_with_bases,
      lbpaps_with_base_quadruples: lbpaps_with_base_quadruples,
      lbpps_with_base_pairs: BP_ALPHABET.keys().map(|base_pair| {(*base_pair, NEG_INFINITY)}).collect(),
      lnbpps_with_bases: SEQ_ALPHABET.iter().map(|&base| {(base, NEG_INFINITY)}).collect(),
    }
  }
}
