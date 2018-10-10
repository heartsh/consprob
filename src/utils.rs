pub use rna_algos::utils::*;

type LbapsWithBasePairs = HashMap<BasePair, LogProb, Hasher>;
type LogGapProbsWithBases = HashMap<Base, LogProb, Hasher>;
pub type BaseQuadruple = (Base, Base, Base, Base);
type LbpapsWithBaseQuadruples = HashMap<BaseQuadruple, LogProb, Hasher>;
type LogGpProbsWithBasePairs = LbapsWithBasePairs;
pub type BaseTriple = (Base, Base, Base);
type LogSgProbsWithBaseTriples = HashMap<BaseTriple, LogProb, Hasher>;
type LogBaseProbsWithBases = HashMap<Base, LogProb, Hasher>;
#[derive(Debug)]
pub struct StemParams {
  pub lbaps_with_base_pairs: LbapsWithBasePairs,
  pub logps_with_bases: LogGapProbsWithBases,
  pub legps_with_bases: LogGapProbsWithBases,
  pub lbpaps_with_base_quadruples_1: LbpapsWithBaseQuadruples,
  pub lbpaps_with_base_quadruples_2: LbpapsWithBaseQuadruples,
  pub logpps_with_base_pairs: LogGpProbsWithBasePairs,
  pub legpps_with_base_pairs: LogGpProbsWithBasePairs,
  pub llgps_with_base_triples: LogSgProbsWithBaseTriples,
  pub lrgps_with_base_triples: LogSgProbsWithBaseTriples,
  pub lbps_with_bases: LogBaseProbsWithBases,
}

pub const SEQ_ALPHABET: [Base; 4] = [A, C, G, U];
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
  pub static ref BPA_ALPHABET_1: Vec<BaseQuadruple> = {
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
  pub static ref BPA_ALPHABET_2: Vec<BaseQuadruple> = {
    let mut bpa_alphabet = Vec::new();
    for &(base_1, base_2) in BP_ALPHABET.keys() {
      for &base_3 in &SEQ_ALPHABET {
        for &base_4 in &SEQ_ALPHABET {
          bpa_alphabet.push((base_1, base_2, base_3, base_4));
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
    let lbpaps_with_base_quadruples_1 = BPA_ALPHABET_1.iter().map(|base_quadruple| {(*base_quadruple, NEG_INFINITY)}).collect();
    let lbpaps_with_base_quadruples_2 = BPA_ALPHABET_2.iter().map(|base_quadruple| {(*base_quadruple, NEG_INFINITY)}).collect();
    let lgpps_with_base_pairs: LogGpProbsWithBasePairs = BP_ALPHABET.keys().map(|base_pair| {(*base_pair, NEG_INFINITY)}).collect();
    let mut lsgps_with_base_triples = LogSgProbsWithBaseTriples::default();
    for &(base_1, base_2) in BP_ALPHABET.keys() {
      for &base_3 in &SEQ_ALPHABET {
        lsgps_with_base_triples.insert((base_1, base_2, base_3), NEG_INFINITY);
      }
    }
    StemParams {
      lbaps_with_base_pairs: lbaps_with_base_pairs,
      logps_with_bases: lgps_with_bases.clone(),
      legps_with_bases: lgps_with_bases,
      lbpaps_with_base_quadruples_1: lbpaps_with_base_quadruples_1,
      lbpaps_with_base_quadruples_2: lbpaps_with_base_quadruples_2,
      logpps_with_base_pairs: lgpps_with_base_pairs.clone(),
      legpps_with_base_pairs: lgpps_with_base_pairs,
      llgps_with_base_triples: lsgps_with_base_triples.clone(),
      lrgps_with_base_triples: lsgps_with_base_triples,
      lbps_with_bases: SEQ_ALPHABET.iter().map(|&base| {(base, NEG_INFINITY)}).collect(),
    }
  }
}
