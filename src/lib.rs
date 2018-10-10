extern crate rna_algos;
extern crate itertools;
#[macro_use]
extern crate lazy_static;
extern crate getopts;

pub mod utils;
pub mod stem_params;

use std::process::{Command, Output};
pub use std::str::from_utf8_unchecked;
pub use getopts::Options;
pub use utils::*;
pub use stem_params::*;

pub type PosTriple = (Pos, Pos, Pos);
pub type PosQuadruple = (Pos, Pos, Pos, Pos);
pub type Prob1dMat = HashMap<Pos, Prob, Hasher>;
pub type SparseProbMat = HashMap<PosPair, Prob, Hasher>;
pub type Prob3dMat = HashMap<PosTriple, Prob, Hasher>;
pub type Prob4dMat = HashMap<PosQuadruple, Prob, Hasher>;
pub struct Stapmt {
  pub base_align_prob_mat: ProbMat,
  pub opening_gap_prob_mat_1: Probs,
  pub ogp_mat_2: Probs,
  pub extending_gap_prob_mat_1: Probs,
  pub egp_mat_2: Probs,
  pub base_pair_align_prob_mat_1: Prob4dMat,
  pub bpap_mat_2: Prob4dMat,
  pub bpap_mat_3: Prob4dMat,
  pub opening_gap_pair_prob_mat_1: SparseProbMat,
  pub ogpp_mat_2: SparseProbMat,
  pub extending_gap_pair_prob_mat_1: SparseProbMat,
  pub egpp_mat_2: SparseProbMat,
  pub left_gap_prob_mat_1: Prob3dMat,
  pub lgp_mat_2: Prob3dMat,
  pub right_gap_prob_mat_1: Prob3dMat,
  pub rgp_mat_2: Prob3dMat,
  pub base_pairing_prob_mat_1: SparseProbMat,
  pub bpp_mat_2: SparseProbMat,
  pub not_bpp_mat_1: Probs,
  pub nbpp_mat_2: Probs,
}
pub type SparseLogProbMat = HashMap<PosPair, LogProb, Hasher>;
pub type LogProb3dMat = HashMap<PosTriple, LogProb, Hasher>;
pub type LogProb4dMat = HashMap<PosQuadruple, LogProb, Hasher>;
#[derive(Clone)]
pub struct LogStapmt {
  pub lbap_mat: LogProbMat,
  pub logp_mat_1: LogProbs,
  pub logp_mat_2: LogProbs,
  pub legp_mat_1: LogProbs,
  pub legp_mat_2: LogProbs,
  pub lbpap_mat_1: LogProb4dMat,
  pub lbpap_mat_2: LogProb4dMat,
  pub lbpap_mat_3: LogProb4dMat,
  pub logpp_mat_1: SparseLogProbMat,
  pub logpp_mat_2: SparseLogProbMat,
  pub legpp_mat_1: SparseLogProbMat,
  pub legpp_mat_2: SparseLogProbMat,
  pub llgp_mat_1: LogProb3dMat,
  pub llgp_mat_2: LogProb3dMat,
  pub lrgp_mat_1: LogProb3dMat,
  pub lrgp_mat_2: LogProb3dMat,
  pub lbpp_mat_1: SparseLogProbMat,
  pub lbpp_mat_2: SparseLogProbMat,
  pub lnbpp_mat_1: LogProbs,
  pub lnbpp_mat_2: LogProbs,
}
type SparseLogPpfMat = HashMap<PosPair, LogPf, Hasher>;
type LogPpf4dMat = HashMap<PosQuadruple, LogPf, Hasher>;
pub struct LogStaPpf4dMats {
  pub log_ppf_mat_4_bpas_1: LogPpf4dMat,
  pub log_ppf_mat_4_bpas_2: LogPpf4dMat,
  pub log_ppf_mat_4_bpas_3: LogPpf4dMat,
  pub log_ppf_mat_4_ogps_1: LogPpf4dMat,
  pub log_ppf_mat_4_ogps_2: LogPpf4dMat,
  pub log_ppf_mat_4_egps_1: LogPpf4dMat,
  pub log_ppf_mat_4_egps_2: LogPpf4dMat,
  pub log_ppf_mat_4_lgs_1: LogPpf4dMat,
  pub log_ppf_mat_4_lgs_2: LogPpf4dMat,
  pub log_ppf_mat_4_rgs_1: LogPpf4dMat,
  pub log_ppf_mat_4_rgs_2: LogPpf4dMat,
}
struct LogStaPpfMats {
  pub log_ppf_mat_1: SparseLogPpfMat,
  pub log_ppf_mat_2: SparseLogPpfMat,
  pub log_ppf_mat_3: SparseLogPpfMat,
  pub log_ppf_mat_4_bas: SparseLogPpfMat,
  pub log_ppf_mat_4_ogs_1: SparseLogPpfMat,
  pub log_ppf_mat_4_ogs_2: SparseLogPpfMat,
  pub log_ppf_mat_4_egs_1: SparseLogPpfMat,
  pub log_ppf_mat_4_egs_2: SparseLogPpfMat,
  pub log_ppf_mat_4_gaps_1: SparseLogPpfMat,
  pub log_ppf_mat_4_gaps_2: SparseLogPpfMat,
}
pub type StaFreeEnergy = LogProb;
pub struct StaFeParams {
  pub lstapmt: LogStapmt,
  pub lstapmt_on_random_assump: LogStapmt,
}
pub type RnaId = usize;
pub type RnaIdPair = (RnaId, RnaId);
pub type StaFeParamSetsWithRnaIdPairs = HashMap<RnaIdPair, StaFeParams, Hasher>;
type StapmtsWithRnaIdPairs = HashMap<RnaIdPair, Stapmt, Hasher>;
pub type LstapmtsWithRnaIdPairs = HashMap<RnaIdPair, LogStapmt, Hasher>;
pub type LogProbMats = Vec<SparseLogProbMat>;
pub type LogProbSeqs = Vec<LogProbs>;
type SeqsOfEpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosPairs = HashMap<PosPair, ExpPartOfTerm4LogProb, Hasher>;
type SeqsOfEpsOfTerms4LogProbsWithPosTriples = HashMap<PosTriple, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosTriples = HashMap<PosTriple, ExpPartOfTerm4LogProb, Hasher>;
type SeqsOfEpsOfTerms4LogProbsWithPosQuadruples = HashMap<PosQuadruple, EpsOfTerms4LogProb, Hasher>;
type EpsOfTerms4LogProbsWithPosQuadruples = HashMap<PosQuadruple, ExpPartOfTerm4LogProb, Hasher>;
type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
pub struct FastaRecord {
  pub fasta_id: FastaId,
  pub seq: Seq,
}
pub type FastaRecords = Vec<FastaRecord>;

const PARASOR_COMMAND: &'static str = "ParasoR";

impl LogStapmt {
  pub fn new(seq_len_pair: &(usize, usize)) -> LogStapmt {
    let log_prob_seq_pair = (vec![NEG_INFINITY; seq_len_pair.0 + 2], vec![NEG_INFINITY; seq_len_pair.1 + 2]);
    let log_prob_mat = SparseLogProbMat::default();
    let log_prob_3d_mat = LogProb3dMat::default();
    let log_prob_4d_mat = LogProb4dMat::default();
    LogStapmt {
      lbap_mat: vec![vec![NEG_INFINITY; seq_len_pair.1 + 2]; seq_len_pair.0 + 2],
      logp_mat_1: log_prob_seq_pair.0.clone(),
      logp_mat_2: log_prob_seq_pair.1.clone(),
      legp_mat_1: log_prob_seq_pair.0.clone(),
      legp_mat_2: log_prob_seq_pair.1.clone(),
      lbpap_mat_1: log_prob_4d_mat.clone(),
      lbpap_mat_2: log_prob_4d_mat.clone(),
      lbpap_mat_3: log_prob_4d_mat,
      logpp_mat_1: log_prob_mat.clone(),
      logpp_mat_2: log_prob_mat.clone(),
      legpp_mat_1: log_prob_mat.clone(),
      legpp_mat_2: log_prob_mat.clone(),
      llgp_mat_1: log_prob_3d_mat.clone(),
      llgp_mat_2: log_prob_3d_mat.clone(),
      lrgp_mat_1: log_prob_3d_mat.clone(),
      lrgp_mat_2: log_prob_3d_mat,
      lbpp_mat_1: log_prob_mat.clone(),
      lbpp_mat_2: log_prob_mat,
      lnbpp_mat_1: log_prob_seq_pair.0,
      lnbpp_mat_2: log_prob_seq_pair.1,
    }
  }
  pub fn origin() -> LogStapmt {
    let log_probs = LogProbs::new();
    let log_prob_mat = SparseLogProbMat::default();
    let log_prob_3d_mat = LogProb3dMat::default();
    let log_prob_4d_mat = LogProb4dMat::default();
    LogStapmt {
      lbap_mat: LogProbMat::new(),
      logp_mat_1: log_probs.clone(),
      logp_mat_2: log_probs.clone(),
      legp_mat_1: log_probs.clone(),
      legp_mat_2: log_probs.clone(),
      lbpap_mat_1: log_prob_4d_mat.clone(),
      lbpap_mat_2: log_prob_4d_mat.clone(),
      lbpap_mat_3: log_prob_4d_mat,
      logpp_mat_1: log_prob_mat.clone(),
      logpp_mat_2: log_prob_mat.clone(),
      legpp_mat_1: log_prob_mat.clone(),
      legpp_mat_2: log_prob_mat.clone(),
      llgp_mat_1: log_prob_3d_mat.clone(),
      llgp_mat_2: log_prob_3d_mat.clone(),
      lrgp_mat_1: log_prob_3d_mat.clone(),
      lrgp_mat_2: log_prob_3d_mat,
      lbpp_mat_1: log_prob_mat.clone(),
      lbpp_mat_2: log_prob_mat,
      lnbpp_mat_1: log_probs.clone(),
      lnbpp_mat_2: log_probs,
    }
  }
}

impl StaFeParams {
  pub fn origin() -> StaFeParams {
    let lstapmt = LogStapmt::origin();
    StaFeParams {
      lstapmt: lstapmt.clone(),
      lstapmt_on_random_assump: lstapmt,
    }
  }
  pub fn new(rna_id_pair: &RnaIdPair, fasta_records: &FastaRecords, gap_num: usize, lbpp_mats: &LogProbMats, lnbpp_mats: &LogProbSeqs) -> StaFeParams {
    let seq_pair = (&fasta_records[rna_id_pair.0].seq[..], &fasta_records[rna_id_pair.1].seq[..]);
    let seq_len_pair = (seq_pair.0.len(), seq_pair.1.len());
    let max_gap_num = get_seq_len_diff(&seq_len_pair) + gap_num;
    let lstapmt = LogStapmt::new(&seq_len_pair);
    let mut sta_fe_params = StaFeParams {
      lstapmt: lstapmt.clone(),
      lstapmt_on_random_assump: lstapmt,
    };
    let lbpp_mat_pair = (&lbpp_mats[rna_id_pair.0], &lbpp_mats[rna_id_pair.1]);
    let pseudo_pos_quadruple = (0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1);
    sta_fe_params.lstapmt.lbpap_mat_1.insert(pseudo_pos_quadruple, 0.);
    sta_fe_params.lstapmt_on_random_assump.lbpap_mat_1.insert(pseudo_pos_quadruple, 0.);
    let pseudo_pos_pair = (0, seq_len_pair.0 + 1); 
    let lbpp = lbpp_mat_pair.0[&pseudo_pos_pair];
    sta_fe_params.lstapmt.lbpp_mat_1.insert(pseudo_pos_pair, lbpp);
    sta_fe_params.lstapmt_on_random_assump.lbpp_mat_1.insert(pseudo_pos_pair, lbpp);
    let pseudo_pos_pair = (0, seq_len_pair.1 + 1); 
    let lbpp = lbpp_mat_pair.1[&pseudo_pos_pair];
    sta_fe_params.lstapmt.lbpp_mat_2.insert(pseudo_pos_pair, lbpp);
    sta_fe_params.lstapmt_on_random_assump.lbpp_mat_2.insert(pseudo_pos_pair, lbpp);
    let lnbpp_mat_pair = (&lnbpp_mats[rna_id_pair.0], &lnbpp_mats[rna_id_pair.1]);
    for i in 0 .. seq_len_pair.0 {
      let base = seq_pair.0[i];
      let lbp = STEM_PARAMS.lbps_with_bases[&base];
      sta_fe_params.lstapmt.logp_mat_1[i + 1] = STEM_PARAMS.logps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.logp_mat_1[i + 1] = lbp;
      sta_fe_params.lstapmt.legp_mat_1[i + 1] = STEM_PARAMS.legps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.legp_mat_1[i + 1] = lbp;
      sta_fe_params.lstapmt.lnbpp_mat_1[i + 1] = lnbpp_mat_pair.0[i + 1];
      sta_fe_params.lstapmt_on_random_assump.lnbpp_mat_1[i + 1] = lbp;
      for j in 0 .. seq_len_pair.1 {
        let pos_pair = (i + 1, j + 1);
        let base_pair = (base, seq_pair.1[j]);
        sta_fe_params.lstapmt.lbap_mat[pos_pair.0][pos_pair.1] = STEM_PARAMS.lbaps_with_base_pairs[&base_pair];
        sta_fe_params.lstapmt_on_random_assump.lbap_mat[pos_pair.0][pos_pair.1] = STEM_PARAMS.lbps_with_bases[&base_pair.0] + STEM_PARAMS.lbps_with_bases[&base_pair.1];
      }
      for j in i + 1 .. seq_len_pair.0 {
        let pos_pair = (i + 1, j + 1);
        let base_pair = (base, seq_pair.0[j]);
        if STEM_PARAMS.logpps_with_base_pairs.contains_key(&base_pair) && lbpp_mat_pair.0.contains_key(&pos_pair) {
          let sum_of_lbp_pair = STEM_PARAMS.lbps_with_bases[&base_pair.0] + STEM_PARAMS.lbps_with_bases[&base_pair.1];
          sta_fe_params.lstapmt.logpp_mat_1.insert(pos_pair, STEM_PARAMS.logpps_with_base_pairs[&base_pair]);
          sta_fe_params.lstapmt_on_random_assump.logpp_mat_1.insert(pos_pair, sum_of_lbp_pair);
          sta_fe_params.lstapmt.legpp_mat_1.insert(pos_pair, STEM_PARAMS.legpps_with_base_pairs[&base_pair]);
          sta_fe_params.lstapmt_on_random_assump.legpp_mat_1.insert(pos_pair, sum_of_lbp_pair);
          sta_fe_params.lstapmt.lbpp_mat_1.insert(pos_pair, lbpp_mat_pair.0[&pos_pair]);
          sta_fe_params.lstapmt_on_random_assump.lbpp_mat_1.insert(pos_pair, sum_of_lbp_pair);
        }
        for k in 0 .. seq_len_pair.1 {
          let pos_triple = (pos_pair.0, pos_pair.1, k + 1);
          let base_triple = (base_pair.0, base_pair.1, seq_pair.1[k]);
          if STEM_PARAMS.llgps_with_base_triples.contains_key(&base_triple) && lbpp_mat_pair.0.contains_key(&pos_pair) {
            let sum_of_lbp_triple = STEM_PARAMS.lbps_with_bases[&base_triple.0] + STEM_PARAMS.lbps_with_bases[&base_triple.1] + STEM_PARAMS.lbps_with_bases[&base_triple.2];
            sta_fe_params.lstapmt.llgp_mat_1.insert(pos_triple, STEM_PARAMS.llgps_with_base_triples[&base_triple]);
            sta_fe_params.lstapmt_on_random_assump.llgp_mat_1.insert(pos_triple, sum_of_lbp_triple);
            sta_fe_params.lstapmt.lrgp_mat_1.insert(pos_triple, STEM_PARAMS.lrgps_with_base_triples[&base_triple]);
            sta_fe_params.lstapmt_on_random_assump.lrgp_mat_1.insert(pos_triple, sum_of_lbp_triple);
          }
          for l in k + 1 .. seq_len_pair.1 {
            let pos_pair_2 = (pos_triple.2, l + 1);
            let pos_quadruple = (pos_pair.0, pos_pair.1, pos_pair_2.0, pos_pair_2.1);
            if get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
            let base_quadruple_1 = (base_triple.0, base_triple.1, base_triple.2, seq_pair.1[l]);
            let base_quadruple_2 = (base_quadruple_1.2, base_quadruple_1.3, base_quadruple_1.0, base_quadruple_1.1);
            let sum_of_lbp_quadruple = STEM_PARAMS.lbps_with_bases[&base_quadruple_1.0] + STEM_PARAMS.lbps_with_bases[&base_quadruple_1.1] + STEM_PARAMS.lbps_with_bases[&base_quadruple_1.2] + STEM_PARAMS.lbps_with_bases[&base_quadruple_1.3];
            if STEM_PARAMS.lbpaps_with_base_quadruples_1.contains_key(&base_quadruple_1) && lbpp_mat_pair.0.contains_key(&pos_pair) && lbpp_mat_pair.1.contains_key(&pos_pair_2) {
              sta_fe_params.lstapmt.lbpap_mat_1.insert(pos_quadruple, STEM_PARAMS.lbpaps_with_base_quadruples_1[&base_quadruple_1]);
              sta_fe_params.lstapmt_on_random_assump.lbpap_mat_1.insert(pos_quadruple, sum_of_lbp_quadruple);
            }
            if STEM_PARAMS.lbpaps_with_base_quadruples_2.contains_key(&base_quadruple_2) && lbpp_mat_pair.1.contains_key(&pos_pair_2) {
              sta_fe_params.lstapmt.lbpap_mat_2.insert(pos_quadruple, STEM_PARAMS.lbpaps_with_base_quadruples_2[&base_quadruple_2]);
              sta_fe_params.lstapmt_on_random_assump.lbpap_mat_2.insert(pos_quadruple, sum_of_lbp_quadruple);
            }
            if STEM_PARAMS.lbpaps_with_base_quadruples_2.contains_key(&base_quadruple_1) && lbpp_mat_pair.0.contains_key(&pos_pair) {
              sta_fe_params.lstapmt.lbpap_mat_3.insert(pos_quadruple, STEM_PARAMS.lbpaps_with_base_quadruples_2[&base_quadruple_1]);
              sta_fe_params.lstapmt_on_random_assump.lbpap_mat_3.insert(pos_quadruple, sum_of_lbp_quadruple);
            }
          }
        }
      }
    }
    for i in 0 .. seq_len_pair.1 {
      let base = seq_pair.1[i];
      sta_fe_params.lstapmt.logp_mat_2[i + 1] = STEM_PARAMS.logps_with_bases[&base];
      let lbp = STEM_PARAMS.lbps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.logp_mat_2[i + 1] = lbp;
      sta_fe_params.lstapmt.legp_mat_2[i + 1] = STEM_PARAMS.legps_with_bases[&base];
      sta_fe_params.lstapmt_on_random_assump.legp_mat_2[i + 1] = lbp;
      sta_fe_params.lstapmt.lnbpp_mat_2[i + 1] = lnbpp_mat_pair.1[i + 1];
      sta_fe_params.lstapmt_on_random_assump.lnbpp_mat_2[i + 1] = lbp;
      for j in i + 1 .. seq_len_pair.1 {
        let pos_pair = (i + 1, j + 1);
        let base_pair = (seq_pair.1[i], seq_pair.1[j]);
        if STEM_PARAMS.logpps_with_base_pairs.contains_key(&base_pair) && lbpp_mat_pair.1.contains_key(&pos_pair) {
          let sum_of_lbp_pair = STEM_PARAMS.lbps_with_bases[&base_pair.0] + STEM_PARAMS.lbps_with_bases[&base_pair.1];
          sta_fe_params.lstapmt.logpp_mat_2.insert(pos_pair, STEM_PARAMS.logpps_with_base_pairs[&base_pair]);
          sta_fe_params.lstapmt_on_random_assump.logpp_mat_2.insert(pos_pair, sum_of_lbp_pair);
          sta_fe_params.lstapmt.legpp_mat_2.insert(pos_pair, STEM_PARAMS.legpps_with_base_pairs[&base_pair]);
          sta_fe_params.lstapmt_on_random_assump.legpp_mat_2.insert(pos_pair, sum_of_lbp_pair);
          sta_fe_params.lstapmt.lbpp_mat_2.insert(pos_pair, lbpp_mat_pair.1[&pos_pair]);
          sta_fe_params.lstapmt_on_random_assump.lbpp_mat_2.insert(pos_pair, sum_of_lbp_pair);
        }
        for k in 0 .. seq_len_pair.0 {
          let pos_triple = (k + 1, pos_pair.0, pos_pair.1);
          let base_triple = (seq_pair.0[k], base_pair.0, base_pair.1);
          if STEM_PARAMS.llgps_with_base_triples.contains_key(&base_triple) && lbpp_mat_pair.1.contains_key(&pos_pair) {
            let sum_of_lbp_triple = STEM_PARAMS.lbps_with_bases[&base_triple.0] + STEM_PARAMS.lbps_with_bases[&base_triple.1] + STEM_PARAMS.lbps_with_bases[&base_triple.2];
            sta_fe_params.lstapmt.llgp_mat_2.insert(pos_triple, STEM_PARAMS.llgps_with_base_triples[&base_triple]);
            sta_fe_params.lstapmt_on_random_assump.llgp_mat_2.insert(pos_triple, sum_of_lbp_triple);
            sta_fe_params.lstapmt.lrgp_mat_2.insert(pos_triple, STEM_PARAMS.lrgps_with_base_triples[&base_triple]);
            sta_fe_params.lstapmt_on_random_assump.lrgp_mat_2.insert(pos_triple, sum_of_lbp_triple);
          }
        }
      }
    }
    sta_fe_params
  }
}

impl LogStaPpf4dMats {
  pub fn new() -> LogStaPpf4dMats {
    let log_ppf_mat = LogPpf4dMat::default();
    LogStaPpf4dMats {
      log_ppf_mat_4_bpas_1: log_ppf_mat.clone(),
      log_ppf_mat_4_bpas_2: log_ppf_mat.clone(),
      log_ppf_mat_4_bpas_3: log_ppf_mat.clone(),
      log_ppf_mat_4_ogps_1: log_ppf_mat.clone(),
      log_ppf_mat_4_ogps_2: log_ppf_mat.clone(),
      log_ppf_mat_4_egps_1: log_ppf_mat.clone(),
      log_ppf_mat_4_egps_2: log_ppf_mat.clone(),
      log_ppf_mat_4_lgs_1: log_ppf_mat.clone(),
      log_ppf_mat_4_lgs_2: log_ppf_mat.clone(),
      log_ppf_mat_4_rgs_1: log_ppf_mat.clone(),
      log_ppf_mat_4_rgs_2: log_ppf_mat,
    }
  }
}

impl LogStaPpfMats {
  pub fn new() -> LogStaPpfMats {
    let log_ppf_mat = SparseLogPpfMat::default();
    LogStaPpfMats {
      log_ppf_mat_1: log_ppf_mat.clone(),
      log_ppf_mat_2: log_ppf_mat.clone(),
      log_ppf_mat_3: log_ppf_mat.clone(),
      log_ppf_mat_4_bas: log_ppf_mat.clone(),
      log_ppf_mat_4_ogs_1: log_ppf_mat.clone(),
      log_ppf_mat_4_ogs_2: log_ppf_mat.clone(),
      log_ppf_mat_4_egs_1: log_ppf_mat.clone(),
      log_ppf_mat_4_egs_2: log_ppf_mat.clone(),
      log_ppf_mat_4_gaps_1: log_ppf_mat.clone(),
      log_ppf_mat_4_gaps_2: log_ppf_mat,
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
pub fn io_algo_4_rna_stapmt(seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_gap_num: usize) -> Stapmt {
  let lstapmt = io_algo_4_rna_lstapmt(seq_len_pair, sta_fe_params, max_gap_num);
  get_stapmt(&lstapmt)
}

#[inline]
pub fn io_algo_4_rna_lstapmt(seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_gap_num: usize) -> LogStapmt {
  let log_sta_ppf_mats = get_log_sta_ppf_4d_mats(seq_len_pair, sta_fe_params, max_gap_num);
  get_lstapmt(seq_len_pair, sta_fe_params, &log_sta_ppf_mats, max_gap_num)
}

#[inline]
fn get_stapmt(lstapmt: &LogStapmt) -> Stapmt {
  Stapmt {
    base_align_prob_mat: lstapmt.lbap_mat.iter().map(|lbaps| {lbaps.iter().map(|&lbap| lbap.exp()).collect()}).collect(),
    opening_gap_prob_mat_1: lstapmt.logp_mat_1.iter().map(|&logp| logp.exp()).collect(),
    ogp_mat_2: lstapmt.logp_mat_2.iter().map(|&logp| logp.exp()).collect(),
    extending_gap_prob_mat_1: lstapmt.legp_mat_1.iter().map(|&legp| legp.exp()).collect(),
    egp_mat_2: lstapmt.legp_mat_2.iter().map(|&legp| legp.exp()).collect(),
    base_pair_align_prob_mat_1: lstapmt.lbpap_mat_1.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect(),
    bpap_mat_2: lstapmt.lbpap_mat_2.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect(),
    bpap_mat_3: lstapmt.lbpap_mat_3.iter().map(|(pos_quadruple, &lbpap)| (*pos_quadruple, lbpap.exp())).collect(),
    opening_gap_pair_prob_mat_1: lstapmt.logpp_mat_1.iter().map(|(pos_pair, &logpp)| (*pos_pair, logpp.exp())).collect(),
    ogpp_mat_2: lstapmt.logpp_mat_2.iter().map(|(pos_pair, &logpp)| (*pos_pair, logpp.exp())).collect(),
    extending_gap_pair_prob_mat_1: lstapmt.legpp_mat_1.iter().map(|(pos_pair, &legpp)| (*pos_pair, legpp.exp())).collect(),
    egpp_mat_2: lstapmt.legpp_mat_2.iter().map(|(pos_pair, &legpp)| (*pos_pair, legpp.exp())).collect(),
    left_gap_prob_mat_1: lstapmt.llgp_mat_1.iter().map(|(pos_triple, &llgp)| (*pos_triple, llgp.exp())).collect(),
    lgp_mat_2: lstapmt.llgp_mat_2.iter().map(|(pos_triple, &llgp)| (*pos_triple, llgp.exp())).collect(),
    right_gap_prob_mat_1: lstapmt.lrgp_mat_1.iter().map(|(pos_triple, &lrgp)| (*pos_triple, lrgp.exp())).collect(),
    rgp_mat_2: lstapmt.lrgp_mat_2.iter().map(|(pos_triple, &lrgp)| (*pos_triple, lrgp.exp())).collect(),
    base_pairing_prob_mat_1: lstapmt.lbpp_mat_1.iter().map(|(pos_pair, &lbpp)| (*pos_pair, lbpp.exp())).collect(),
    bpp_mat_2: lstapmt.lbpp_mat_2.iter().map(|(pos_pair, &lbpp)| (*pos_pair, lbpp.exp())).collect(),
    not_bpp_mat_1: lstapmt.lnbpp_mat_1.iter().map(|&lnbpp| lnbpp.exp()).collect(),
    nbpp_mat_2: lstapmt.lnbpp_mat_2.iter().map(|&lnbpp| lnbpp.exp()).collect(),
  }
}

#[inline]
pub fn get_log_sta_ppf_4d_mats(seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, max_gap_num: usize) -> LogStaPpf4dMats {
  let mut log_sta_ppf_mats = LogStaPpf4dMats::new();
  for substr_len_1 in 2 .. seq_len_pair.0 + 3 {
    for substr_len_2 in 2 .. seq_len_pair.1 + 3 {
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_gap_num {continue;}
      for i in 0 .. seq_len_pair.0 + 3 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 0 .. seq_len_pair.1 + 3 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !(sta_fe_params.lstapmt.lbpp_mat_1.contains_key(&(i, j)) || sta_fe_params.lstapmt.lbpp_mat_2.contains_key(&(k, l))) {continue;}
          let log_sta_pf = get_log_sta_forward_ppf_mats(&pos_quadruple, sta_fe_params, &log_sta_ppf_mats, max_gap_num).log_ppf_mat_1[&(j - 1, l - 1)];
          if !log_sta_pf.is_finite() {continue;}
          if sta_fe_params.lstapmt.lbpap_mat_1.contains_key(&pos_quadruple) {
            log_sta_ppf_mats.log_ppf_mat_4_bpas_1.insert(pos_quadruple, get_bpa_lor_1(&pos_quadruple, sta_fe_params) + log_sta_pf);
          }
          if sta_fe_params.lstapmt.lbpap_mat_2.contains_key(&pos_quadruple) {
            log_sta_ppf_mats.log_ppf_mat_4_bpas_2.insert(pos_quadruple, get_bpa_lor_2(&pos_quadruple, sta_fe_params) + log_sta_pf);
          }
          if sta_fe_params.lstapmt.lbpap_mat_3.contains_key(&pos_quadruple) {
            log_sta_ppf_mats.log_ppf_mat_4_bpas_3.insert(pos_quadruple, get_bpa_lor_3(&pos_quadruple, sta_fe_params) + log_sta_pf);
          }
          let pos_pair = (i, j);
          let pos_quadruple = (i, j, k + 1, l - 1);
          if sta_fe_params.lstapmt.logpp_mat_1.contains_key(&pos_pair) {
            log_sta_ppf_mats.log_ppf_mat_4_ogps_1.insert(pos_quadruple, get_ogp_lor_1(&pos_pair, sta_fe_params) + log_sta_pf);
          }
          let pos_pair = (k, l);
          let pos_quadruple = (i + 1, j - 1, k, l);
          if sta_fe_params.lstapmt.logpp_mat_2.contains_key(&pos_pair) {
            log_sta_ppf_mats.log_ppf_mat_4_ogps_2.insert(pos_quadruple, get_ogp_lor_2(&pos_pair, sta_fe_params) + log_sta_pf);
          }
          let pos_pair = (i, j);
          let pos_quadruple = (i, j, k + 1, l - 1);
          if sta_fe_params.lstapmt.legpp_mat_1.contains_key(&pos_pair) {
            log_sta_ppf_mats.log_ppf_mat_4_egps_1.insert(pos_quadruple, get_egp_lor_1(&pos_pair, sta_fe_params) + log_sta_pf);
          }
          let pos_pair = (k, l);
          let pos_quadruple = (i + 1, j - 1, k, l);
          if sta_fe_params.lstapmt.legpp_mat_2.contains_key(&pos_pair) {
            log_sta_ppf_mats.log_ppf_mat_4_egps_2.insert(pos_quadruple, get_egp_lor_2(&pos_pair, sta_fe_params) + log_sta_pf);
          }
          let pos_triple = (i, j, l);
          let pos_quadruple = (i, j, k + 1, l);
          if sta_fe_params.lstapmt.llgp_mat_1.contains_key(&pos_triple) {
            log_sta_ppf_mats.log_ppf_mat_4_lgs_1.insert(pos_quadruple, get_lg_lor_1(&pos_triple, sta_fe_params) + log_sta_pf);
          }
          let pos_triple = (j, k, l);
          let pos_quadruple = (i + 1, j, k, l);
          if sta_fe_params.lstapmt.llgp_mat_2.contains_key(&pos_triple) {
            log_sta_ppf_mats.log_ppf_mat_4_lgs_2.insert(pos_quadruple, get_lg_lor_2(&pos_triple, sta_fe_params) + log_sta_pf);
          }
          let pos_triple = (i, j, k);
          let pos_quadruple = (i, j, k, l - 1);
          if sta_fe_params.lstapmt.lrgp_mat_1.contains_key(&pos_triple) {
            log_sta_ppf_mats.log_ppf_mat_4_rgs_1.insert(pos_quadruple, get_rg_lor_1(&pos_triple, sta_fe_params) + log_sta_pf);
          }
          let pos_triple = (i, k, l);
          let pos_quadruple = (i, j - 1, k, l);
          if sta_fe_params.lstapmt.lrgp_mat_2.contains_key(&pos_triple) {
            log_sta_ppf_mats.log_ppf_mat_4_rgs_2.insert(pos_quadruple, get_rg_lor_2(&pos_triple, sta_fe_params) + log_sta_pf);
          }
        }
      }
    }
  }
  log_sta_ppf_mats
}

#[inline]
fn get_log_sta_forward_ppf_mats(pos_quadruple: &PosQuadruple, sta_fe_params: &StaFeParams, log_sta_ppf_4d_mats: &LogStaPpf4dMats, max_gap_num: usize) -> LogStaPpfMats {
  let &(i, j, k, l) = pos_quadruple;
  let mut log_sta_ppf_mats = LogStaPpfMats::new();
  for n in i .. j {
    for p in k .. l {
      let pos_pair_1 = (n, p);
      if n == i && p == k {
        log_sta_ppf_mats.log_ppf_mat_4_bas.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_1.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_2.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_3.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
      } else {
        if get_min_gap_num(&(i, n, k, p)) > max_gap_num {continue;}
        if n == i || p == k {
          log_sta_ppf_mats.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let pos_pair_2 = (n - 1, p - 1);
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair_2] + get_ba_log_odds_ratio(&pos_pair_1, sta_fe_params);
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num || !(sta_fe_params.lstapmt.lbpp_mat_1.contains_key(&(m, n)) || sta_fe_params.lstapmt.lbpp_mat_2.contains_key(&(o, p))) {continue;}
              let pos_pair = (m - 1, o - 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_2[&pos_pair];
              let log_sta_pf_3 = log_sta_ppf_mats.log_ppf_mat_3[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_1 + log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_1 + log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_2[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_3.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_1 + log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_3[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_2 + log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_3 + log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if n == i {
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(i, n - 1, k, p)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (n - 1, p);
            max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_2[&pos_pair_2] + get_og_lor_1(n, sta_fe_params);
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (m - 1, o - 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_2[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_2 + log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_1 + log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if p == k {
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(i, n, k, p - 1)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (n, p - 1);
            max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_3[&pos_pair_2] + get_og_lor_2(p, sta_fe_params);
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (m - 1, o - 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_3[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_2 + log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_1 + log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if n == i {
          log_sta_ppf_mats.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(i, n - 1, k, p)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (n - 1, p);
            max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_gaps_1[&pos_pair_2] + get_eg_lor_1(n, sta_fe_params);
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (m - 1, o - 1);
              if log_sta_ppf_4d_mats.log_ppf_mat_4_egps_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_gaps_1[&pos_pair] + log_sta_ppf_4d_mats.log_ppf_mat_4_egps_1[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_egs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if p == k {
          log_sta_ppf_mats.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(i, n, k, p - 1)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (n, p - 1);
            max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_gaps_2[&pos_pair_2] + get_eg_lor_2(p, sta_fe_params);
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (m - 1, o - 1);
              if log_sta_ppf_4d_mats.log_ppf_mat_4_egps_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_gaps_2[&pos_pair] + log_sta_ppf_4d_mats.log_ppf_mat_4_egps_2[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_egs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if n == i {
          log_sta_ppf_mats.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let pos_pair_2 = (n - 1, p);
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_egs_1[&pos_pair_1];
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          if get_min_gap_num(&(i, n - 1, k, p)) <= max_gap_num && max_gap_num >= 1 {
            let og_lor = get_og_lor_1(n, sta_fe_params);
            let ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (m - 1, o - 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair];
              let log_sta_pf_3 = log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_1 + log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1.contains_key(&pos_quadruple) {
                let log_sta_pf_4_ogp = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1[&pos_quadruple];
                let ep_of_term_4_log_pf = log_sta_pf_2 + log_sta_pf_4_ogp;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
                let ep_of_term_4_log_pf = log_sta_pf_3 + log_sta_pf_4_ogp;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if p == k {
          log_sta_ppf_mats.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let pos_pair_2 = (n, p - 1);
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_egs_2[&pos_pair_1];
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          if get_min_gap_num(&(i, n, k, p - 1)) <= max_gap_num && max_gap_num >= 1 {
            let og_lor = get_og_lor_2(p, sta_fe_params);
            let ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair_2] + og_lor;
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          for m in i + 1 .. n {
            for o in k + 1 .. p {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (m - 1, o - 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair];
              let log_sta_pf_3 = log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_pf_1 + log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2[&pos_quadruple];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2.contains_key(&pos_quadruple) {
                let log_sta_pf_4_ogp = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2[&pos_quadruple];
                let ep_of_term_4_log_pf = log_sta_pf_2 + log_sta_pf_4_ogp;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
                let ep_of_term_4_log_pf = log_sta_pf_3 + log_sta_pf_4_ogp;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mats.log_ppf_mat_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_gaps_1[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mats.log_ppf_mat_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_gaps_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mats.log_ppf_mat_3.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
      }
    }
  }
  log_sta_ppf_mats
}

#[inline]
fn get_ba_log_odds_ratio(pos_pair: &PosPair, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.lbap_mat[pos_pair.0][pos_pair.1] - sta_fe_params.lstapmt_on_random_assump.lbap_mat[pos_pair.0][pos_pair.1]
}

#[inline]
fn get_og_lor_1(pos: Pos, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.logp_mat_1[pos] - sta_fe_params.lstapmt_on_random_assump.logp_mat_1[pos]
}

#[inline]
fn get_og_lor_2(pos: Pos, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.logp_mat_2[pos] - sta_fe_params.lstapmt_on_random_assump.logp_mat_2[pos]
}

#[inline]
fn get_eg_lor_1(pos: Pos, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.legp_mat_1[pos] - sta_fe_params.lstapmt_on_random_assump.legp_mat_1[pos]
}

#[inline]
fn get_eg_lor_2(pos: Pos, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.legp_mat_2[pos] - sta_fe_params.lstapmt_on_random_assump.legp_mat_2[pos]
}

#[inline]
fn get_bp_lor_1(pos_pair: &PosPair, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.lbpp_mat_1[&pos_pair] - sta_fe_params.lstapmt_on_random_assump.lbpp_mat_1[&pos_pair]
}

#[inline]
fn get_bp_lor_2(pos_pair: &PosPair, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.lbpp_mat_2[&pos_pair] - sta_fe_params.lstapmt_on_random_assump.lbpp_mat_2[&pos_pair]
}

#[inline]
fn get_nbp_lor_1(pos: Pos, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.lnbpp_mat_1[pos] - sta_fe_params.lstapmt_on_random_assump.lnbpp_mat_1[pos]
}

#[inline]
fn get_nbp_lor_2(pos: Pos, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.lnbpp_mat_2[pos] - sta_fe_params.lstapmt_on_random_assump.lnbpp_mat_2[pos]
}

#[inline]
fn get_bpa_lor_1(pos_quadruple: &PosQuadruple, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  let bp_pos_pair_1 = (pos_quadruple.0, pos_quadruple.1);
  let bp_pos_pair_2 = (pos_quadruple.2, pos_quadruple.3);
  sta_fe_params.lstapmt.lbpap_mat_1[pos_quadruple]
  - sta_fe_params.lstapmt_on_random_assump.lbpap_mat_1[pos_quadruple]
  + get_bp_lor_1(&bp_pos_pair_1, sta_fe_params)
  + get_bp_lor_2(&bp_pos_pair_2, sta_fe_params)
}

#[inline]
fn get_bpa_lor_2(pos_quadruple: &PosQuadruple, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  let nbp_pos_pair = (pos_quadruple.0, pos_quadruple.1);
  let bp_pos_pair = (pos_quadruple.2, pos_quadruple.3);
  sta_fe_params.lstapmt.lbpap_mat_2[pos_quadruple]
  - sta_fe_params.lstapmt_on_random_assump.lbpap_mat_2[pos_quadruple]
  + get_nbp_lor_1(nbp_pos_pair.0, sta_fe_params)
  + get_nbp_lor_1(nbp_pos_pair.1, sta_fe_params)
  + get_bp_lor_2(&bp_pos_pair, sta_fe_params)
}

#[inline]
fn get_bpa_lor_3(pos_quadruple: &PosQuadruple, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  let bp_pos_pair = (pos_quadruple.0, pos_quadruple.1);
  let nbp_pos_pair = (pos_quadruple.2, pos_quadruple.3);
  sta_fe_params.lstapmt.lbpap_mat_3[pos_quadruple]
  - sta_fe_params.lstapmt_on_random_assump.lbpap_mat_3[pos_quadruple]
  + get_bp_lor_1(&bp_pos_pair, sta_fe_params)
  + get_nbp_lor_2(nbp_pos_pair.0, sta_fe_params)
  + get_nbp_lor_2(nbp_pos_pair.1, sta_fe_params)
}

#[inline]
fn get_ogp_lor_1(pos_pair: &PosPair, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.logpp_mat_1[pos_pair]
  - sta_fe_params.lstapmt_on_random_assump.logpp_mat_1[pos_pair]
  + get_bp_lor_1(&pos_pair, sta_fe_params)
}

#[inline]
fn get_ogp_lor_2(pos_pair: &PosPair, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.logpp_mat_2[pos_pair]
  - sta_fe_params.lstapmt_on_random_assump.logpp_mat_2[pos_pair]
  + get_bp_lor_2(&pos_pair, sta_fe_params)
}

#[inline]
fn get_egp_lor_1(pos_pair: &PosPair, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.legpp_mat_1[pos_pair]
  - sta_fe_params.lstapmt_on_random_assump.legpp_mat_1[pos_pair]
  + get_bp_lor_1(&pos_pair, sta_fe_params)
}

#[inline]
fn get_egp_lor_2(pos_pair: &PosPair, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  sta_fe_params.lstapmt.legpp_mat_2[pos_pair]
  - sta_fe_params.lstapmt_on_random_assump.legpp_mat_2[pos_pair]
  + get_bp_lor_2(&pos_pair, sta_fe_params)
}

#[inline]
fn get_lg_lor_1(pos_triple: &PosTriple, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  let bp_pos_pair = (pos_triple.0, pos_triple.1);
  let nbp_pos = pos_triple.2;
  sta_fe_params.lstapmt.llgp_mat_1[pos_triple]
  - sta_fe_params.lstapmt_on_random_assump.llgp_mat_1[pos_triple]
  + get_bp_lor_1(&bp_pos_pair, sta_fe_params)
  + get_nbp_lor_2(nbp_pos, sta_fe_params)
}

#[inline]
fn get_lg_lor_2(pos_triple: &PosTriple, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  let bp_pos_pair = (pos_triple.1, pos_triple.2);
  let nbp_pos = pos_triple.0;
  sta_fe_params.lstapmt.llgp_mat_2[pos_triple]
  - sta_fe_params.lstapmt_on_random_assump.llgp_mat_2[pos_triple]
  + get_bp_lor_2(&bp_pos_pair, sta_fe_params)
  + get_nbp_lor_1(nbp_pos, sta_fe_params)
}

#[inline]
fn get_rg_lor_1(pos_triple: &PosTriple, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  let bp_pos_pair = (pos_triple.0, pos_triple.1);
  let nbp_pos = pos_triple.2;
  sta_fe_params.lstapmt.lrgp_mat_1[pos_triple]
  - sta_fe_params.lstapmt_on_random_assump.lrgp_mat_1[pos_triple]
  + get_bp_lor_1(&bp_pos_pair, sta_fe_params)
  + get_nbp_lor_2(nbp_pos, sta_fe_params)
}

#[inline]
fn get_rg_lor_2(pos_triple: &PosTriple, sta_fe_params: &StaFeParams) -> StaFreeEnergy {
  let bp_pos_pair = (pos_triple.1, pos_triple.2);
  let nbp_pos = pos_triple.0;
  sta_fe_params.lstapmt.lrgp_mat_2[pos_triple]
  - sta_fe_params.lstapmt_on_random_assump.lrgp_mat_2[pos_triple]
  + get_bp_lor_2(&bp_pos_pair, sta_fe_params)
  + get_nbp_lor_1(nbp_pos, sta_fe_params)
}

#[inline]
fn get_lstapmt(seq_len_pair: &(usize, usize), sta_fe_params: &StaFeParams, log_sta_ppf_4d_mats: &LogStaPpf4dMats, max_gap_num: usize) -> LogStapmt {
  let mut lstapmt = LogStapmt::new(seq_len_pair);
  let mut mat_of_seqs_of_eps_of_terms_4_log_probs = vec![vec![EpsOfTerms4LogProb::new(); seq_len_pair.1 + 2]; seq_len_pair.0 + 2];
  let mut mat_of_min_eps_of_terms_4_log_probs = vec![vec![INFINITY; seq_len_pair.1 + 2]; seq_len_pair.0 + 2];
  let mut seqs_of_eps_of_terms_4_logps_1 = vec![EpsOfTerms4LogProb::new(); seq_len_pair.0 + 2];
  let mut min_eps_of_terms_4_logps_1 = vec![INFINITY; seq_len_pair.0 + 2];
  let mut seqs_of_eps_of_terms_4_logps_2 = vec![EpsOfTerms4LogProb::new(); seq_len_pair.1 + 2];
  let mut min_eps_of_terms_4_logps_2 = vec![INFINITY; seq_len_pair.1 + 2];
  let mut seqs_of_eps_of_terms_4_legps_1 = seqs_of_eps_of_terms_4_logps_1.clone();
  let mut min_eps_of_terms_4_legps_1 = min_eps_of_terms_4_logps_1.clone();
  let mut seqs_of_eps_of_terms_4_legps_2 = seqs_of_eps_of_terms_4_logps_2.clone();
  let mut min_eps_of_terms_4_legps_2 = min_eps_of_terms_4_logps_2.clone();
  let mut seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_1 = SeqsOfEpsOfTerms4LogProbsWithPosQuadruples::default();
  let mut min_eps_of_terms_4_lbpaps_with_pos_quadruples_1 = EpsOfTerms4LogProbsWithPosQuadruples::default();
  let mut seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_2 = seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_1.clone();
  let mut min_eps_of_terms_4_lbpaps_with_pos_quadruples_2 = min_eps_of_terms_4_lbpaps_with_pos_quadruples_1.clone();
  let mut seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_3 = seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_1.clone();
  let mut min_eps_of_terms_4_lbpaps_with_pos_quadruples_3 = min_eps_of_terms_4_lbpaps_with_pos_quadruples_1.clone();
  let mut seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1 = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut min_eps_of_terms_4_logpps_with_pos_pairs_1 = EpsOfTerms4LogProbsWithPosPairs::default();
  let mut seqs_of_eps_of_terms_4_logpps_with_pos_pairs_2 = seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1.clone();
  let mut min_eps_of_terms_4_logpps_with_pos_pairs_2 = min_eps_of_terms_4_logpps_with_pos_pairs_1.clone();
  let mut seqs_of_eps_of_terms_4_legpps_with_pos_pairs_1 = seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1.clone();
  let mut min_eps_of_terms_4_legpps_with_pos_pairs_1 = min_eps_of_terms_4_logpps_with_pos_pairs_1.clone();
  let mut seqs_of_eps_of_terms_4_legpps_with_pos_pairs_2 = seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1.clone();
  let mut min_eps_of_terms_4_legpps_with_pos_pairs_2 = min_eps_of_terms_4_logpps_with_pos_pairs_1.clone();
  let mut seqs_of_eps_of_terms_4_llgps_with_pos_triples_1 = SeqsOfEpsOfTerms4LogProbsWithPosTriples::default();
  let mut min_eps_of_terms_4_llgps_with_pos_triples_1 = EpsOfTerms4LogProbsWithPosTriples::default();
  let mut seqs_of_eps_of_terms_4_llgps_with_pos_triples_2 = seqs_of_eps_of_terms_4_llgps_with_pos_triples_1.clone();
  let mut min_eps_of_terms_4_llgps_with_pos_triples_2 = min_eps_of_terms_4_llgps_with_pos_triples_1.clone();
  let mut seqs_of_eps_of_terms_4_lrgps_with_pos_triples_1 = seqs_of_eps_of_terms_4_llgps_with_pos_triples_1.clone();
  let mut min_eps_of_terms_4_lrgps_with_pos_triples_1 = min_eps_of_terms_4_llgps_with_pos_triples_1.clone();
  let mut seqs_of_eps_of_terms_4_lrgps_with_pos_triples_2 = seqs_of_eps_of_terms_4_llgps_with_pos_triples_1.clone();
  let mut min_eps_of_terms_4_lrgps_with_pos_triples_2 = min_eps_of_terms_4_llgps_with_pos_triples_1.clone();
  let mut seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1 = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut min_eps_of_terms_4_lbpps_with_pos_pairs_1 = EpsOfTerms4LogProbsWithPosPairs::default();
  let mut seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2 = seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.clone();
  let mut min_eps_of_terms_4_lbpps_with_pos_pairs_2 = min_eps_of_terms_4_lbpps_with_pos_pairs_1.clone();
  let mut seqs_of_eps_of_terms_4_lnbpps_1 = vec![EpsOfTerms4LogProb::new(); seq_len_pair.0 + 2];
  let mut min_eps_of_terms_4_lnbpps_1 = vec![INFINITY; seq_len_pair.0 + 2];
  let mut seqs_of_eps_of_terms_4_lnbpps_2 = vec![EpsOfTerms4LogProb::new(); seq_len_pair.1 + 2];
  let mut min_eps_of_terms_4_lnbpps_2 = vec![INFINITY; seq_len_pair.1 + 2];
  let pseudo_pos_quadruple = (0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1);
  for pos_pair in sta_fe_params.lstapmt.lbpp_mat_1.keys() {
    seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.insert(*pos_pair, EpsOfTerms4LogProb::new());
    min_eps_of_terms_4_lbpps_with_pos_pairs_1.insert(*pos_pair, INFINITY);
  }
  for pos_pair in sta_fe_params.lstapmt.lbpp_mat_2.keys() {
    seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2.insert(*pos_pair, EpsOfTerms4LogProb::new());
    min_eps_of_terms_4_lbpps_with_pos_pairs_2.insert(*pos_pair, INFINITY);
  }
  for substr_len_1 in (3 .. seq_len_pair.0 + 3).rev() {
    for substr_len_2 in (3 .. seq_len_pair.1 + 3).rev() {
      if max(substr_len_1, substr_len_2) - min(substr_len_1, substr_len_2) > max_gap_num {continue;}
      for i in 0 .. seq_len_pair.0 + 3 - substr_len_1 {
        let j = i + substr_len_1 - 1;
        for k in 0 .. seq_len_pair.1 + 3 - substr_len_2 {
          let l = k + substr_len_2 - 1;
          let pos_quadruple = (i, j, k, l);
          if !(sta_fe_params.lstapmt.lbpp_mat_1.contains_key(&(i, j)) || sta_fe_params.lstapmt.lbpp_mat_2.contains_key(&(k, l))) {continue;}
          if pos_quadruple == pseudo_pos_quadruple {
            lstapmt.lbpap_mat_1.insert(pos_quadruple, 0.);
          } else {
            if min_eps_of_terms_4_lbpaps_with_pos_quadruples_1.contains_key(&pos_quadruple) {
              let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpaps_with_pos_quadruples_1[&pos_quadruple];
              if min_ep_of_term_4_log_prob.is_finite() {
                let lbpap = logsumexp(&seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_1[&pos_quadruple][..], min_ep_of_term_4_log_prob);
                if lbpap.is_finite() {
                  lstapmt.lbpap_mat_1.insert(pos_quadruple, lbpap);
                  let pos_pair_1 = (i, j);
                  let pos_pair_2 = (k, l);
                  if min_eps_of_terms_4_lbpps_with_pos_pairs_1.contains_key(&pos_pair_1) {
                    seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair_1).expect("Failed to get an element of a hash map.").push(lbpap);
                    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair_1).expect("Failed to get an element of a hash map.");
                    if lbpap < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = lbpap;
                    }
                  }
                  if min_eps_of_terms_4_lbpps_with_pos_pairs_2.contains_key(&pos_pair_2) {
                    seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair_2).expect("Failed to get an element of a hash map.").push(lbpap);
                    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair_2).expect("Failed to get an element of a hash map.");
                    if lbpap < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = lbpap;
                    }
                  }
                }
              }
            }
            if min_eps_of_terms_4_lbpaps_with_pos_quadruples_2.contains_key(&pos_quadruple) {
              let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpaps_with_pos_quadruples_2[&pos_quadruple];
              if min_ep_of_term_4_log_prob.is_finite() {
                let lbpap = logsumexp(&seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_2[&pos_quadruple][..], min_ep_of_term_4_log_prob);
                if lbpap.is_finite() {
                  lstapmt.lbpap_mat_2.insert(pos_quadruple, lbpap);
                  let pos_pair = (k, l);
                  seqs_of_eps_of_terms_4_lnbpps_1[i].push(lbpap);
                  let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[i];
                  if lbpap < min_ep_of_term_4_log_prob {
                    min_eps_of_terms_4_lnbpps_1[i] = lbpap;
                  }
                  seqs_of_eps_of_terms_4_lnbpps_1[j].push(lbpap);
                  let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[j];
                  if lbpap < min_ep_of_term_4_log_prob {
                    min_eps_of_terms_4_lnbpps_1[j] = lbpap;
                  }
                  if min_eps_of_terms_4_lbpps_with_pos_pairs_2.contains_key(&pos_pair) {
                    seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(lbpap);
                    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                    if lbpap < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = lbpap;
                    }
                  }
                }
              }
            }
            if min_eps_of_terms_4_lbpaps_with_pos_quadruples_3.contains_key(&pos_quadruple) {
              let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpaps_with_pos_quadruples_3[&pos_quadruple];
              if min_ep_of_term_4_log_prob.is_finite() {
                let lbpap = logsumexp(&seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_3[&pos_quadruple][..], min_ep_of_term_4_log_prob);
                if lbpap.is_finite() {
                  lstapmt.lbpap_mat_3.insert(pos_quadruple, lbpap);
                  let pos_pair = (i, j);
                  if min_eps_of_terms_4_lbpps_with_pos_pairs_1.contains_key(&pos_pair) {
                    seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(lbpap);
                    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                    if lbpap < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = lbpap;
                    }
                  }
                  seqs_of_eps_of_terms_4_lnbpps_2[k].push(lbpap);
                  let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[k];
                  if lbpap < min_ep_of_term_4_log_prob {
                    min_eps_of_terms_4_lnbpps_2[k] = lbpap;
                  }
                  seqs_of_eps_of_terms_4_lnbpps_2[l].push(lbpap);
                  let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[l];
                  if lbpap < min_ep_of_term_4_log_prob {
                    min_eps_of_terms_4_lnbpps_2[l] = lbpap;
                  }
                }
              }
            }
          }
          if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1.contains_key(&pos_quadruple) && lstapmt.lbpap_mat_1.contains_key(&pos_quadruple) {
            let log_sta_pf_4_bpa = log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1[&pos_quadruple];
            let lbpap = lstapmt.lbpap_mat_1[&pos_quadruple];
            let log_sta_forward_ppf_mats = get_log_sta_forward_ppf_mats(&pos_quadruple, sta_fe_params, &log_sta_ppf_4d_mats, max_gap_num);
            let log_sta_backward_ppf_mats = get_log_sta_backward_ppf_mats(&pos_quadruple, sta_fe_params, &log_sta_ppf_4d_mats, max_gap_num);
            for n in i + 1 .. j {
              for p in k + 1 .. l {
                if get_min_gap_num(&(i, n, k, p)) > max_gap_num || get_min_gap_num(&(n, j, p, l)) > max_gap_num {continue;}
                let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(n - 1, p - 1)] + get_ba_log_odds_ratio(&(n, p), sta_fe_params) + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                if ep_of_term_4_log_prob.is_finite() {
                  mat_of_seqs_of_eps_of_terms_4_log_probs[n][p].push(ep_of_term_4_log_prob);
                  let ref mut min_ep_of_term_4_log_prob = mat_of_min_eps_of_terms_4_log_probs[n][p];
                  if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                    *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
                }
                if get_min_gap_num(&(i, n - 1, k, p)) <= max_gap_num && max_gap_num >= 1 {
                  let og_lor = get_og_lor_1(n, sta_fe_params);
                  let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_2[&(n - 1, p)] + og_lor + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                  if ep_of_term_4_log_prob.is_finite() {
                    seqs_of_eps_of_terms_4_logps_1[n].push(ep_of_term_4_log_prob);
                    let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_logps_1[n];
                    if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                    }
                  }
                  let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(n - 1, p)] + og_lor + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                  if ep_of_term_4_log_prob.is_finite() {
                    seqs_of_eps_of_terms_4_logps_1[n].push(ep_of_term_4_log_prob);
                    let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_logps_1[n];
                    if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                    }
                  }
                }
                if get_min_gap_num(&(i, n, k, p - 1)) <= max_gap_num && max_gap_num >= 1 {
                  let og_lor = get_og_lor_2(p, sta_fe_params);
                  let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_3[&(n, p - 1)] + og_lor + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                  if ep_of_term_4_log_prob.is_finite() {
                    seqs_of_eps_of_terms_4_logps_2[p].push(ep_of_term_4_log_prob);
                    let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_logps_2[p];
                    if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                    }
                  }
                  let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(n, p - 1)] + og_lor + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                  if ep_of_term_4_log_prob.is_finite() {
                    seqs_of_eps_of_terms_4_logps_2[p].push(ep_of_term_4_log_prob);
                    let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_logps_2[p];
                    if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                    }
                  }
                }
                if get_min_gap_num(&(i, n - 1, k, p)) <= max_gap_num && max_gap_num >= 1 {
                  let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_4_gaps_1[&(n - 1, p)] + get_eg_lor_1(n, sta_fe_params) + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                  if ep_of_term_4_log_prob.is_finite() {
                    seqs_of_eps_of_terms_4_legps_1[n].push(ep_of_term_4_log_prob);
                    let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_legps_1[n];
                    if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                    }
                  }
                }
                if get_min_gap_num(&(i, n, k, p - 1)) <= max_gap_num && max_gap_num >= 1 {
                  let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_4_gaps_2[&(n, p - 1)] + get_eg_lor_2(p, sta_fe_params) + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                  if ep_of_term_4_log_prob.is_finite() {
                    seqs_of_eps_of_terms_4_legps_2[p].push(ep_of_term_4_log_prob);
                    let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_legps_2[p];
                    if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                      *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                    }
                  }
                }
              }
            }
            for m in i + 1 .. j {
              for n in m + 1 .. j {
                for o in k + 1 .. l {
                  for p in o + 1 .. l {
                    let pos_quadruple = (m, n, o, p);
                    if get_min_gap_num(&(i, m, k, o)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num || get_min_gap_num(&(n, j, p, l)) > max_gap_num || !(sta_fe_params.lstapmt.lbpp_mat_1.contains_key(&(m, n)) || sta_fe_params.lstapmt.lbpp_mat_2.contains_key(&(o, p))) {continue;}
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1.contains_key(&pos_quadruple) {
                      if !min_eps_of_terms_4_lbpaps_with_pos_quadruples_1.contains_key(&pos_quadruple) {
                        seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_1.insert(pos_quadruple, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_lbpaps_with_pos_quadruples_1.insert(pos_quadruple, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_1.get_mut(&pos_quadruple).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpaps_with_pos_quadruples_1.get_mut(&pos_quadruple).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_2.contains_key(&pos_quadruple) {
                      if !min_eps_of_terms_4_lbpaps_with_pos_quadruples_2.contains_key(&pos_quadruple) {
                        seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_2.insert(pos_quadruple, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_lbpaps_with_pos_quadruples_2.insert(pos_quadruple, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_2[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_2.get_mut(&pos_quadruple).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpaps_with_pos_quadruples_2.get_mut(&pos_quadruple).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_3.contains_key(&pos_quadruple) {
                      if !min_eps_of_terms_4_lbpaps_with_pos_quadruples_3.contains_key(&pos_quadruple) {
                        seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_3.insert(pos_quadruple, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_lbpaps_with_pos_quadruples_3.insert(pos_quadruple, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_3[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_lbpaps_with_pos_quadruples_3.get_mut(&pos_quadruple).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpaps_with_pos_quadruples_3.get_mut(&pos_quadruple).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1.contains_key(&pos_quadruple) {
                      let pos_pair = (m, n);
                      if !min_eps_of_terms_4_logpps_with_pos_pairs_1.contains_key(&pos_pair) {
                        seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1.insert(pos_pair, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_logpps_with_pos_pairs_1.insert(pos_pair, INFINITY);
                      }
                      let log_sta_pf_4_ogp = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1[&pos_quadruple];
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_2[&(m - 1, o - 1)] + log_sta_pf_4_ogp + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(m - 1, o - 1)] + log_sta_pf_4_ogp + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2.contains_key(&pos_quadruple) {
                      let pos_pair = (o, p);
                      if !min_eps_of_terms_4_logpps_with_pos_pairs_2.contains_key(&pos_pair) {
                        seqs_of_eps_of_terms_4_logpps_with_pos_pairs_2.insert(pos_pair, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_logpps_with_pos_pairs_2.insert(pos_pair, INFINITY);
                      }
                      let log_sta_pf_4_ogp = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2[&pos_quadruple];
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_3[&(m - 1, o - 1)] + log_sta_pf_4_ogp + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_logpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(m - 1, o - 1)] + log_sta_pf_4_ogp + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_logpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_egps_1.contains_key(&pos_quadruple) {
                      let pos_pair = (m, n);
                      if !min_eps_of_terms_4_legpps_with_pos_pairs_1.contains_key(&pos_pair) {
                        seqs_of_eps_of_terms_4_legpps_with_pos_pairs_1.insert(pos_pair, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_legpps_with_pos_pairs_1.insert(pos_pair, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_4_gaps_1[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_egps_1[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_legpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_legpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_egps_2.contains_key(&pos_quadruple) {
                      let pos_pair = (o, p);
                      if !min_eps_of_terms_4_legpps_with_pos_pairs_2.contains_key(&pos_pair) {
                        seqs_of_eps_of_terms_4_legpps_with_pos_pairs_2.insert(pos_pair, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_legpps_with_pos_pairs_2.insert(pos_pair, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_4_gaps_2[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_egps_2[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_4_gaps_2[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_legpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_legpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1.contains_key(&pos_quadruple) {
                      let pos_triple = (m, n, p);
                      if !min_eps_of_terms_4_llgps_with_pos_triples_1.contains_key(&pos_triple) {
                        seqs_of_eps_of_terms_4_llgps_with_pos_triples_1.insert(pos_triple, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_llgps_with_pos_triples_1.insert(pos_triple, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_2[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_llgps_with_pos_triples_1.get_mut(&pos_triple).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_llgps_with_pos_triples_1.get_mut(&pos_triple).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2.contains_key(&pos_quadruple) {
                      let pos_triple = (n, o, p);
                      if !min_eps_of_terms_4_llgps_with_pos_triples_2.contains_key(&pos_triple) {
                        seqs_of_eps_of_terms_4_llgps_with_pos_triples_2.insert(pos_triple, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_llgps_with_pos_triples_2.insert(pos_triple, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_3[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_1[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_llgps_with_pos_triples_2.get_mut(&pos_triple).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_llgps_with_pos_triples_2.get_mut(&pos_triple).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1.contains_key(&pos_quadruple) {
                      let pos_triple = (m, n, o);
                      if !min_eps_of_terms_4_lrgps_with_pos_triples_1.contains_key(&pos_triple) {
                        seqs_of_eps_of_terms_4_lrgps_with_pos_triples_1.insert(pos_triple, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_lrgps_with_pos_triples_1.insert(pos_triple, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_2[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_lrgps_with_pos_triples_1.get_mut(&pos_triple).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lrgps_with_pos_triples_1.get_mut(&pos_triple).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                    if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2.contains_key(&pos_quadruple) {
                      let pos_triple = (m, o, p);
                      if !min_eps_of_terms_4_lrgps_with_pos_triples_2.contains_key(&pos_triple) {
                        seqs_of_eps_of_terms_4_lrgps_with_pos_triples_2.insert(pos_triple, EpsOfTerms4LogProb::new());
                        min_eps_of_terms_4_lrgps_with_pos_triples_2.insert(pos_triple, INFINITY);
                      }
                      let ep_of_term_4_log_prob = lbpap + log_sta_forward_ppf_mats.log_ppf_mat_1[&(m - 1, o - 1)] + log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2[&pos_quadruple] + log_sta_backward_ppf_mats.log_ppf_mat_3[&(n + 1, p + 1)] - log_sta_pf_4_bpa;
                      if ep_of_term_4_log_prob.is_finite() {
                        seqs_of_eps_of_terms_4_lrgps_with_pos_triples_2.get_mut(&pos_triple).expect("Failed to get an element of a hash map.").push(ep_of_term_4_log_prob);
                        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lrgps_with_pos_triples_2.get_mut(&pos_triple).expect("Failed to get an element of a hash map.");
                        if ep_of_term_4_log_prob < *min_ep_of_term_4_log_prob {
                          *min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  for i in 0 .. seq_len_pair.0 + 2 {
    for j in 0 .. seq_len_pair.1 + 2 {
      let min_ep_of_term_4_log_prob = mat_of_min_eps_of_terms_4_log_probs[i][j];
      let lbap = if min_ep_of_term_4_log_prob.is_finite() {
        logsumexp(&mat_of_seqs_of_eps_of_terms_4_log_probs[i][j][..], min_ep_of_term_4_log_prob)
      } else {
        NEG_INFINITY
      };
      if lbap.is_finite() {
        lstapmt.lbap_mat[i][j] = lbap;
        seqs_of_eps_of_terms_4_lnbpps_1[i].push(lbap);
        let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[i];
        if lbap < *min_ep_of_term_4_log_prob {
          *min_ep_of_term_4_log_prob = lbap;
        }
        seqs_of_eps_of_terms_4_lnbpps_2[j].push(lbap);
        let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[j];
        if lbap < *min_ep_of_term_4_log_prob {
          *min_ep_of_term_4_log_prob = lbap;
        }
      }
    }
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logps_1[i];
    let logp = if min_ep_of_term_4_log_prob.is_finite() {
      logsumexp(&seqs_of_eps_of_terms_4_logps_1[i][..], min_ep_of_term_4_log_prob)
    } else {
      NEG_INFINITY
    };
    if logp.is_finite() {
      lstapmt.logp_mat_1[i] = logp;
      seqs_of_eps_of_terms_4_lnbpps_1[i].push(logp);
      let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[i];
      if logp < min_ep_of_term_4_log_prob {
        min_eps_of_terms_4_lnbpps_1[i]= logp;
      }
    }
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_legps_1[i];
    let legp = if min_ep_of_term_4_log_prob.is_finite() {
      logsumexp(&seqs_of_eps_of_terms_4_legps_1[i][..], min_ep_of_term_4_log_prob)
    } else {
      NEG_INFINITY
    };
    if legp.is_finite() {
      lstapmt.legp_mat_1[i] = legp;
      seqs_of_eps_of_terms_4_lnbpps_1[i].push(legp);
      let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[i];
      if legp < min_ep_of_term_4_log_prob {
        min_eps_of_terms_4_lnbpps_1[i] = legp;
      }
    }
  }
  for i in 0 .. seq_len_pair.1 + 2 {
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logps_2[i];
    let logp = if min_ep_of_term_4_log_prob.is_finite() {
      logsumexp(&seqs_of_eps_of_terms_4_logps_2[i][..], min_ep_of_term_4_log_prob)
    } else {
      NEG_INFINITY
    };
    if logp.is_finite() {
      lstapmt.logp_mat_2[i] = logp;
      seqs_of_eps_of_terms_4_lnbpps_2[i].push(logp);
      let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[i];
      if logp < min_ep_of_term_4_log_prob {
        min_eps_of_terms_4_lnbpps_2[i] = logp;
      }
    }
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_legps_2[i];
    let legp = if min_ep_of_term_4_log_prob.is_finite() {
      logsumexp(&seqs_of_eps_of_terms_4_legps_2[i][..], min_ep_of_term_4_log_prob)
    } else {
      NEG_INFINITY
    };
    if legp.is_finite() {
      lstapmt.legp_mat_2[i] = legp;
      seqs_of_eps_of_terms_4_lnbpps_2[i].push(legp);
      let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[i];
      if legp < *min_ep_of_term_4_log_prob {
        *min_ep_of_term_4_log_prob = legp;
      }
    }
  }
  for i in 0 .. seq_len_pair.0 + 1 {
    for j in i + 1 .. seq_len_pair.0 + 2 {
      let pos_pair = (i, j);
      for k in 0 .. seq_len_pair.1 + 2 {
        let pos_triple = (i, j, k);
        if min_eps_of_terms_4_llgps_with_pos_triples_1.contains_key(&pos_triple) {
          let min_ep_of_term_4_log_prob = min_eps_of_terms_4_llgps_with_pos_triples_1[&pos_triple];
          if min_ep_of_term_4_log_prob.is_finite() {
            let llgp = logsumexp(&seqs_of_eps_of_terms_4_llgps_with_pos_triples_1[&pos_triple][..], min_ep_of_term_4_log_prob);
            if llgp.is_finite() {
              lstapmt.llgp_mat_1.insert(pos_triple, llgp);
              seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(llgp);
              let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
              if llgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = llgp;
              }
              seqs_of_eps_of_terms_4_lnbpps_2[k].push(llgp);
              let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[k];
              if llgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = llgp;
              }
            }
          }
        }
        if min_eps_of_terms_4_lrgps_with_pos_triples_1.contains_key(&pos_triple) {
          let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lrgps_with_pos_triples_1[&pos_triple];
          if min_ep_of_term_4_log_prob.is_finite() {
            let lrgp = logsumexp(&seqs_of_eps_of_terms_4_lrgps_with_pos_triples_1[&pos_triple][..], min_ep_of_term_4_log_prob);
            if lrgp.is_finite() {
              lstapmt.lrgp_mat_1.insert(pos_triple, lrgp);
              seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(lrgp);
              let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
              if lrgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = lrgp;
              }
              seqs_of_eps_of_terms_4_lnbpps_2[k].push(lrgp);
              let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[k];
              if lrgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = lrgp;
              }
            }
          }
        }
      }
      if min_eps_of_terms_4_logpps_with_pos_pairs_1.contains_key(&pos_pair) {
        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logpps_with_pos_pairs_1[&pos_pair];
        if min_ep_of_term_4_log_prob.is_finite() {
          let logpp = logsumexp(&seqs_of_eps_of_terms_4_logpps_with_pos_pairs_1[&pos_pair][..], min_ep_of_term_4_log_prob);
          if logpp.is_finite() {
            lstapmt.logpp_mat_1.insert(pos_pair, logpp);
            seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(logpp);
            let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
            if logpp < *min_ep_of_term_4_log_prob {
              *min_ep_of_term_4_log_prob = logpp;
            }
          }
        }
      }
      if min_eps_of_terms_4_legpps_with_pos_pairs_1.contains_key(&pos_pair) {
        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_legpps_with_pos_pairs_1[&pos_pair];
        if min_ep_of_term_4_log_prob.is_finite() {
          let legpp = logsumexp(&seqs_of_eps_of_terms_4_legpps_with_pos_pairs_1[&pos_pair][..], min_ep_of_term_4_log_prob);
          if legpp.is_finite() {
            lstapmt.legpp_mat_1.insert(pos_pair, legpp);
            seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(legpp);
            let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_1.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
            if legpp < *min_ep_of_term_4_log_prob {
              *min_ep_of_term_4_log_prob = legpp;
            }
          }
        }
      }
    }
  }
  for i in 0 .. seq_len_pair.1 + 1 {
    for j in i + 1 .. seq_len_pair.1 + 2 {
      let pos_pair = (i, j);
      for k in 0 .. seq_len_pair.0 + 2 {
        let pos_triple = (k, i, j);
        if min_eps_of_terms_4_llgps_with_pos_triples_2.contains_key(&pos_triple) {
          let min_ep_of_term_4_log_prob = min_eps_of_terms_4_llgps_with_pos_triples_2[&pos_triple];
          if min_ep_of_term_4_log_prob.is_finite() {
            let llgp = logsumexp(&seqs_of_eps_of_terms_4_llgps_with_pos_triples_2[&pos_triple][..], min_ep_of_term_4_log_prob);
            if llgp.is_finite() {
              lstapmt.llgp_mat_2.insert(pos_triple, llgp);
              seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(llgp);
              let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
              if llgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = llgp;
              }
              seqs_of_eps_of_terms_4_lnbpps_1[k].push(llgp);
              let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[k];
              if llgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = llgp;
              }
            }
          }
        }
        if min_eps_of_terms_4_lrgps_with_pos_triples_2.contains_key(&pos_triple) {
          let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lrgps_with_pos_triples_2[&pos_triple];
          if min_ep_of_term_4_log_prob.is_finite() {
            let lrgp = logsumexp(&seqs_of_eps_of_terms_4_lrgps_with_pos_triples_2[&pos_triple][..], min_ep_of_term_4_log_prob);
            if lrgp.is_finite() {
              lstapmt.lrgp_mat_2.insert(pos_triple, lrgp);
              seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(lrgp);
              let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
              if lrgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = lrgp;
              }
              seqs_of_eps_of_terms_4_lnbpps_1[k].push(lrgp);
              let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[k];
              if lrgp < *min_ep_of_term_4_log_prob {
                *min_ep_of_term_4_log_prob = lrgp;
              }
            }
          }
        }
      }
      let pos_pair = (i, j);
      if min_eps_of_terms_4_logpps_with_pos_pairs_2.contains_key(&pos_pair) {
        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_logpps_with_pos_pairs_2[&pos_pair];
        if min_ep_of_term_4_log_prob.is_finite() {
          let logpp = logsumexp(&seqs_of_eps_of_terms_4_logpps_with_pos_pairs_2[&pos_pair][..], min_ep_of_term_4_log_prob);
          if logpp.is_finite() {
            lstapmt.logpp_mat_2.insert(pos_pair, logpp);
            seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(logpp);
            let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
            if logpp < *min_ep_of_term_4_log_prob {
              *min_ep_of_term_4_log_prob = logpp;
            }
          }
        }
      }
      if min_eps_of_terms_4_legpps_with_pos_pairs_2.contains_key(&pos_pair) {
        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_legpps_with_pos_pairs_2[&pos_pair];
        if min_ep_of_term_4_log_prob.is_finite() {
          let legpp = logsumexp(&seqs_of_eps_of_terms_4_legpps_with_pos_pairs_2[&pos_pair][..], min_ep_of_term_4_log_prob);
          if legpp.is_finite() {
            lstapmt.legpp_mat_2.insert(pos_pair, legpp);
            seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.").push(legpp);
            let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_2.get_mut(&pos_pair).expect("Failed to get an element of a hash map.");
            if legpp < *min_ep_of_term_4_log_prob {
              *min_ep_of_term_4_log_prob = legpp;
            }
          }
        }
      }
    }
  }
  for i in 0 .. seq_len_pair.0 + 2 {
    for j in i + 1 .. seq_len_pair.0 + 2 {
      let pos_pair = (i, j);
      if min_eps_of_terms_4_lbpps_with_pos_pairs_1.contains_key(&pos_pair) {
        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_1[&pos_pair];
        if min_ep_of_term_4_log_prob.is_finite() {
          lstapmt.lbpp_mat_1.insert(pos_pair, logsumexp(&seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_1[&pos_pair][..], min_ep_of_term_4_log_prob));
        }
      }
    }
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_1[i];
    if min_ep_of_term_4_log_prob.is_finite() {
      lstapmt.lnbpp_mat_1[i] = logsumexp(&seqs_of_eps_of_terms_4_lnbpps_1[i][..], min_ep_of_term_4_log_prob);
    }
  }
  for i in 0 .. seq_len_pair.1 + 2 {
    for j in i + 1 .. seq_len_pair.1 + 2 {
      let pos_pair = (i, j);
      if min_eps_of_terms_4_lbpps_with_pos_pairs_2.contains_key(&pos_pair) {
        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lbpps_with_pos_pairs_2[&pos_pair];
        if min_ep_of_term_4_log_prob.is_finite() {
          lstapmt.lbpp_mat_2.insert(pos_pair, logsumexp(&seqs_of_eps_of_terms_4_lbpps_with_pos_pairs_2[&pos_pair][..], min_ep_of_term_4_log_prob));
        }
      }
    }
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_lnbpps_2[i];
    if min_ep_of_term_4_log_prob.is_finite() {
      lstapmt.lnbpp_mat_2[i] = logsumexp(&seqs_of_eps_of_terms_4_lnbpps_2[i][..], min_ep_of_term_4_log_prob);
    }
  }
  lstapmt
}

#[inline]
fn get_log_sta_backward_ppf_mats(pos_quadruple: &PosQuadruple, sta_fe_params: &StaFeParams, log_sta_ppf_4d_mats: &LogStaPpf4dMats, max_gap_num: usize) -> LogStaPpfMats {
  let &(i, j, k, l) = pos_quadruple;
  let mut log_sta_ppf_mats = LogStaPpfMats::new();
  for m in (i + 1 .. j + 1).rev() {
    for o in (k + 1 .. l + 1).rev() {
      let pos_pair_1 = (m, o);
      if m == j && o == l {
        log_sta_ppf_mats.log_ppf_mat_4_bas.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_1.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_2.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_3.insert(pos_pair_1, 0.);
        log_sta_ppf_mats.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        log_sta_ppf_mats.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
      } else {
        if get_min_gap_num(&(m, j, o, l)) > max_gap_num {continue;}
        if m == j || o == l {
          log_sta_ppf_mats.log_ppf_mat_4_bas.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let pos_pair_2 = (m + 1, o + 1);
          let mut max_ep_of_term_4_log_pf = get_ba_log_odds_ratio(&pos_pair_1, sta_fe_params) + log_sta_ppf_mats.log_ppf_mat_1[&pos_pair_2];
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(n, j, p, l)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num || !(sta_fe_params.lstapmt.lbpp_mat_1.contains_key(&(m, n)) || sta_fe_params.lstapmt.lbpp_mat_2.contains_key(&(o, p))) {continue;}
              let pos_pair = (n + 1, p + 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_2[&pos_pair];
              let log_sta_pf_3 = log_sta_ppf_mats.log_ppf_mat_3[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_1[&pos_quadruple] + log_sta_pf_1;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_2[&pos_quadruple] + log_sta_pf_1;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_3.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_bpas_3[&pos_quadruple] + log_sta_pf_1;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_1[&pos_quadruple] + log_sta_pf_2;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_rgs_2[&pos_quadruple] + log_sta_pf_3;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_bas.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if m == j {
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(m + 1, j, o, l)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (m + 1, o);
            max_ep_of_term_4_log_pf = get_og_lor_1(m, sta_fe_params) + log_sta_ppf_mats.log_ppf_mat_2[&pos_pair_2];
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(n, j, p, l)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (n + 1, p + 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_2[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1[&pos_quadruple] + log_sta_pf_2;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1[&pos_quadruple] + log_sta_pf_1;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if o == l {
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(m, j, o + 1, l)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (m, o + 1);
            max_ep_of_term_4_log_pf = get_og_lor_2(o, sta_fe_params) + log_sta_ppf_mats.log_ppf_mat_3[&pos_pair_2];
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(n, j, p, l)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (n + 1, p + 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_3[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2[&pos_quadruple] + log_sta_pf_2;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2[&pos_quadruple] + log_sta_pf_1;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if m == j {
          log_sta_ppf_mats.log_ppf_mat_4_egs_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(m + 1, j, o, l)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (m + 1, o);
            max_ep_of_term_4_log_pf = get_eg_lor_1(m, sta_fe_params) + log_sta_ppf_mats.log_ppf_mat_4_gaps_1[&pos_pair_2];
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(n, j, p, l)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (n + 1, p + 1);
              if log_sta_ppf_4d_mats.log_ppf_mat_4_egps_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_egps_1[&pos_quadruple] + log_sta_ppf_mats.log_ppf_mat_4_gaps_1[&pos_pair];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_egs_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if o == l {
          log_sta_ppf_mats.log_ppf_mat_4_egs_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let mut max_ep_of_term_4_log_pf;
          if get_min_gap_num(&(m, j, o + 1, l)) <= max_gap_num && max_gap_num >= 1 {
            let pos_pair_2 = (m, o + 1);
            max_ep_of_term_4_log_pf = get_eg_lor_2(o, sta_fe_params) + log_sta_ppf_mats.log_ppf_mat_4_gaps_2[&pos_pair_2];
            if max_ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
            }
          } else {
            max_ep_of_term_4_log_pf = NEG_INFINITY;
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(n, j, p, l)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (n + 1, p + 1);
              if log_sta_ppf_4d_mats.log_ppf_mat_4_egps_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_egps_2[&pos_quadruple] + log_sta_ppf_mats.log_ppf_mat_4_gaps_2[&pos_pair];
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_egs_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if m == j {
          log_sta_ppf_mats.log_ppf_mat_4_gaps_1.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let pos_pair_2 = (m + 1, o);
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_egs_1[&pos_pair_1];
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          if get_min_gap_num(&(m + 1, j, o, l)) <= max_gap_num && max_gap_num >= 1 {
            let og_lor = get_og_lor_1(m, sta_fe_params);
            let ep_of_term_4_log_pf = og_lor + log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_2];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = og_lor + log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair_2];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(n, j, p, l)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (n + 1, p + 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair];
              let log_sta_pf_3 = log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_1[&pos_quadruple] + log_sta_pf_1;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1.contains_key(&pos_quadruple) {
                let log_sta_pf_4_ogp = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_1[&pos_quadruple];
                let ep_of_term_4_log_pf = log_sta_pf_4_ogp + log_sta_pf_2;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
                let ep_of_term_4_log_pf = log_sta_pf_4_ogp + log_sta_pf_3;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_gaps_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        if o == l {
          log_sta_ppf_mats.log_ppf_mat_4_gaps_2.insert(pos_pair_1, NEG_INFINITY);
        } else {
          let mut eps_of_terms_4_log_pf = EpsOfTerms4LogPf::new();
          let pos_pair_2 = (m, o + 1);
          let mut max_ep_of_term_4_log_pf = log_sta_ppf_mats.log_ppf_mat_4_egs_2[&pos_pair_1];
          if max_ep_of_term_4_log_pf.is_finite() {
            eps_of_terms_4_log_pf.push(max_ep_of_term_4_log_pf);
          }
          if get_min_gap_num(&(m, j, o + 1, l)) <= max_gap_num && max_gap_num >= 1 {
            let og_lor = get_og_lor_2(o, sta_fe_params);
            let ep_of_term_4_log_pf = og_lor + log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_2];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
            let ep_of_term_4_log_pf = og_lor + log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair_2];
            if ep_of_term_4_log_pf.is_finite() {
              eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
              if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
              }
            }
          }
          for n in m + 1 .. j {
            for p in o + 1 .. l {
              let pos_quadruple = (m, n, o, p);
              if get_min_gap_num(&(n, j, p, l)) > max_gap_num || get_min_gap_num(&pos_quadruple) > max_gap_num {continue;}
              let pos_pair = (n + 1, p + 1);
              let log_sta_pf_1 = log_sta_ppf_mats.log_ppf_mat_1[&pos_pair];
              let log_sta_pf_2 = log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair];
              let log_sta_pf_3 = log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair];
              if log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2.contains_key(&pos_quadruple) {
                let ep_of_term_4_log_pf = log_sta_ppf_4d_mats.log_ppf_mat_4_lgs_2[&pos_quadruple] + log_sta_pf_1;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
              if log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2.contains_key(&pos_quadruple) {
                let log_sta_pf_4_ogp = log_sta_ppf_4d_mats.log_ppf_mat_4_ogps_2[&pos_quadruple];
                let ep_of_term_4_log_pf = log_sta_pf_4_ogp + log_sta_pf_2;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
                let ep_of_term_4_log_pf = log_sta_pf_4_ogp + log_sta_pf_3;
                if ep_of_term_4_log_pf.is_finite() {
                  eps_of_terms_4_log_pf.push(ep_of_term_4_log_pf);
                  if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
                    max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
                  }
                }
              }
            }
          }
          log_sta_ppf_mats.log_ppf_mat_4_gaps_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
            logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
          } else {
            NEG_INFINITY
          });
        }
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mats.log_ppf_mat_1.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_2[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_gaps_1[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mats.log_ppf_mat_2.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
        let eps_of_terms_4_log_pf = [
          log_sta_ppf_mats.log_ppf_mat_4_bas[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_ogs_1[&pos_pair_1],
          log_sta_ppf_mats.log_ppf_mat_4_gaps_2[&pos_pair_1],
        ];
        let mut max_ep_of_term_4_log_pf = NEG_INFINITY;
        for &ep_of_term_4_log_pf in &eps_of_terms_4_log_pf {
          if ep_of_term_4_log_pf > max_ep_of_term_4_log_pf {
            max_ep_of_term_4_log_pf = ep_of_term_4_log_pf;
          }
        }
        log_sta_ppf_mats.log_ppf_mat_3.insert(pos_pair_1, if max_ep_of_term_4_log_pf.is_finite() {
          logsumexp(&eps_of_terms_4_log_pf[..], max_ep_of_term_4_log_pf)
        } else {
          NEG_INFINITY
        });
      }
    }
  }
  log_sta_ppf_mats
}

#[inline]
pub fn prob_cons_transformation_of_lstapmt(lstapmts_with_rna_id_pairs: &LstapmtsWithRnaIdPairs, rna_id_pair: &RnaIdPair, num_of_rnas: usize) -> LogStapmt {
  let mut lstapmt = lstapmts_with_rna_id_pairs[rna_id_pair].clone();
  let log_coefficient = -((num_of_rnas - 1) as Prob).ln();
  for (pos_quadruple, lbpap) in lstapmt.lbpap_mat_1.iter_mut() {
    let (i, j, k, l) = *pos_quadruple;
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut min_ep_of_term_4_log_prob = *lbpap;
    if min_ep_of_term_4_log_prob.is_finite() {
      eps_of_terms_4_log_prob.push(min_ep_of_term_4_log_prob);
    }
    for rna_id in 0 .. num_of_rnas {
      if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
      for (&(m, n, o, p), &lbpap_1) in lstapmts_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].lbpap_mat_1.iter() {
        if (rna_id_pair.0 < rna_id && m != i && n != j) || (rna_id_pair.0 > rna_id && o != i && p != j) {continue;}
        for (&(m, n, o, p), &lbpap_2) in lstapmts_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].lbpap_mat_1.iter() {
          if (rna_id_pair.1 < rna_id && m != k && n != l) || (rna_id_pair.1 > rna_id && o != k && p != l) {continue;}
          let ep_of_term_4_log_prob = lbpap_1 + lbpap_2;
          if ep_of_term_4_log_prob.is_finite() {
            eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
            if ep_of_term_4_log_prob < min_ep_of_term_4_log_prob {
              min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
            }
          }
        }
      }
    }
    if eps_of_terms_4_log_prob.len() > 0 {
      *lbpap = log_coefficient + logsumexp(&eps_of_terms_4_log_prob[..], min_ep_of_term_4_log_prob);
    }
  }
  for (i, lbaps) in lstapmt.lbap_mat.iter_mut().enumerate() {
    for (j, lbap) in lbaps.iter_mut().enumerate() {
      let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
      let mut min_ep_of_term_4_log_prob = *lbap;
      if min_ep_of_term_4_log_prob.is_finite() {
        eps_of_terms_4_log_prob.push(min_ep_of_term_4_log_prob);
      }
      for rna_id in 0 .. num_of_rnas {
        if rna_id == rna_id_pair.0 || rna_id == rna_id_pair.1 {continue;}
        for (k, lbaps_1) in lstapmts_with_rna_id_pairs[& if rna_id_pair.0 < rna_id {(rna_id_pair.0, rna_id)} else {(rna_id, rna_id_pair.0)}].lbap_mat.iter().enumerate() {
          for (l, lbap_1) in lbaps_1.iter().enumerate() {
            if (rna_id_pair.0 < rna_id && k != i) || (rna_id_pair.0 > rna_id && l != i) {continue;}
            for (k, lbaps_2) in lstapmts_with_rna_id_pairs[& if rna_id_pair.1 < rna_id {(rna_id_pair.1, rna_id)} else {(rna_id, rna_id_pair.1)}].lbap_mat.iter().enumerate() {
              for (l, lbap_2) in lbaps_2.iter().enumerate() {
                if (rna_id_pair.1 < rna_id && k != j) || (rna_id_pair.1 > rna_id && l != j) {continue;}
                let ep_of_term_4_log_prob = lbap_1 + lbap_2;
                if ep_of_term_4_log_prob.is_finite() {
                  eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
                  if ep_of_term_4_log_prob < min_ep_of_term_4_log_prob {
                    min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
                  }
                }
              }
            }
          }
        }
      }
      if eps_of_terms_4_log_prob.len() > 0 {
        *lbap = log_coefficient + logsumexp(&eps_of_terms_4_log_prob[..], min_ep_of_term_4_log_prob);
      }
    }
  }
  lstapmt
}

#[inline]
pub fn pct_of_lbpp_mat(lstapmts_with_rna_id_pairs: &LstapmtsWithRnaIdPairs, rna_id: RnaId, num_of_rnas: usize, lbpp_mat: &SparseLogProbMat) -> SparseLogProbMat {
  let log_coefficient = -(num_of_rnas as Prob).ln();
  let mut seqs_of_eps_of_terms_4_log_probs_with_pos_pairs = SeqsOfEpsOfTerms4LogProbsWithPosPairs::default();
  let mut min_eps_of_terms_4_log_probs_with_pos_pairs = EpsOfTerms4LogProbsWithPosPairs::default();
  for (pos_pair, &lbpp) in lbpp_mat.iter() {
    seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, vec![lbpp]);
    min_eps_of_terms_4_log_probs_with_pos_pairs.insert(*pos_pair, lbpp);
  }
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id == rna_id_2 {continue;}
    let rna_id_pair = if rna_id < rna_id_2 {(rna_id, rna_id_2)} else {(rna_id_2, rna_id)};
    let ref ref_2_lstapmt = lstapmts_with_rna_id_pairs[&rna_id_pair];
    for (pos_pair, &lbpp) in if rna_id < rna_id_2 {ref_2_lstapmt.lbpp_mat_1.iter()} else {ref_2_lstapmt.lbpp_mat_2.iter()} {
      if lbpp.is_finite() {
        let eps_of_terms_4_log_prob = seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
        eps_of_terms_4_log_prob.push(lbpp);
        let min_ep_of_term_4_log_prob = min_eps_of_terms_4_log_probs_with_pos_pairs.get_mut(pos_pair).expect("Failed to get an element of a hash map.");
        if lbpp < *min_ep_of_term_4_log_prob {
          *min_ep_of_term_4_log_prob = lbpp;
        }
      }
    }
  }
  let mut lbpp_mat = lbpp_mat.clone();
  for (pos_pair, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs_with_pos_pairs.iter() {
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_log_probs_with_pos_pairs[pos_pair];
    let lbpp = log_coefficient + logsumexp(eps_of_terms_4_log_prob, min_ep_of_term_4_log_prob);
    *lbpp_mat.get_mut(pos_pair).expect("Failed to get an element of a hash map.") = lbpp;
  }
  lbpp_mat
}

#[inline]
pub fn pct_of_lnbpp_mat(lstapmts_with_rna_id_pairs: &LstapmtsWithRnaIdPairs, rna_id: RnaId, num_of_rnas: usize, lnbpp_mat: &LogProbs, seq_len: usize) -> LogProbs {
  let log_coefficient = -(num_of_rnas as Prob).ln();
  let mut seqs_of_eps_of_terms_4_log_probs = vec![EpsOfTerms4LogProb::new(); seq_len + 2];
  let mut min_eps_of_terms_4_log_probs = vec![NEG_INFINITY; seq_len + 2];
  for i in 0 .. seq_len + 2 {
    let lnbpp = lnbpp_mat[i];
    seqs_of_eps_of_terms_4_log_probs[i].push(lnbpp);
    min_eps_of_terms_4_log_probs[i] = lnbpp;
  }
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id == rna_id_2 {continue;}
    let rna_id_pair = if rna_id < rna_id_2 {(rna_id, rna_id_2)} else {(rna_id_2, rna_id)};
    let ref ref_2_lstapmt = lstapmts_with_rna_id_pairs[&rna_id_pair];
    for (i, &lnbpp) in if rna_id < rna_id_2 {ref_2_lstapmt.lnbpp_mat_1.iter().enumerate()} else {ref_2_lstapmt.lnbpp_mat_2.iter().enumerate()} {
      if lnbpp.is_finite() {
        let ref mut eps_of_terms_4_log_prob = seqs_of_eps_of_terms_4_log_probs[i];
        eps_of_terms_4_log_prob.push(lnbpp);
        let ref mut min_ep_of_term_4_log_prob = min_eps_of_terms_4_log_probs[i];
        if lnbpp < *min_ep_of_term_4_log_prob {
          *min_ep_of_term_4_log_prob = lnbpp;
        }
      }
    }
  }
  let mut lnbpp_mat = lnbpp_mat.clone();
  for (i, eps_of_terms_4_log_prob) in seqs_of_eps_of_terms_4_log_probs.iter().enumerate() {
    let min_ep_of_term_4_log_prob = min_eps_of_terms_4_log_probs[i];
    lnbpp_mat[i] = log_coefficient + logsumexp(eps_of_terms_4_log_prob, min_ep_of_term_4_log_prob);
  }
  lnbpp_mat
}

#[inline]
pub fn update_lstapmts_with_rna_id_pairs(lstapmts_with_rna_id_pairs: &mut LstapmtsWithRnaIdPairs, lbpp_mats: &LogProbMats, lnbpp_mats: &LogProbSeqs) {
  for (&(rna_id_1, rna_id_2), lstapmt) in lstapmts_with_rna_id_pairs.iter_mut() {
    lstapmt.lbpp_mat_1 = lbpp_mats[rna_id_1].clone();
    lstapmt.lbpp_mat_2 = lbpp_mats[rna_id_2].clone();
    lstapmt.lnbpp_mat_1 = lnbpp_mats[rna_id_1].clone();
    lstapmt.lnbpp_mat_2 = lnbpp_mats[rna_id_2].clone();
  }
}

#[inline]
pub fn get_stapmts_with_rna_id_pairs(lstapmts_with_rna_id_pairs: &LstapmtsWithRnaIdPairs) -> StapmtsWithRnaIdPairs {
  let mut stapmts_with_rna_id_pairs = StapmtsWithRnaIdPairs::default();
  for (rna_id_pair, lstapmt) in lstapmts_with_rna_id_pairs.iter() {
    stapmts_with_rna_id_pairs.insert(*rna_id_pair, get_stapmt(lstapmt));
  }
  stapmts_with_rna_id_pairs
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
pub fn get_lbpp_mat(seq: SeqSlice, max_bp_span: usize) -> SparseLogProbMat {
  let arg = max_bp_span.to_string();
  let args = unsafe {["-f", from_utf8_unchecked(seq), "--pre", "--constraint", &arg]};
  let parasor_output = unsafe {String::from_utf8_unchecked(run_command(PARASOR_COMMAND, &args, "Failed to run the ParasoR program.").stdout)};
  let seq_len = seq.len();
  let mut lbpp_mat = SparseLogProbMat::default();
  lbpp_mat.insert((0, seq_len + 1), 0.);
  for line in parasor_output.lines().filter(|line| {!line.starts_with("#")}) {
    let strings = line.trim().splitn(4, '\t').collect::<Vec<&str>>();
    let bpp = strings[3].parse::<Prob>().expect("Failed to parse a string.");
    if bpp > 0. {
      lbpp_mat.insert((
        strings[1].parse().expect("Failed to parse a string."),
        strings[2].parse().expect("Failed to parse a string."),
      ), bpp.ln());
    }
  }
  lbpp_mat
}

#[inline]
pub fn get_lnbpp_mat(lbpp_mat: &SparseLogProbMat, seq_len: usize) -> LogProbs {
  let mut lnbpp_mat = vec![NEG_INFINITY; seq_len + 2];
  /* lnbpp_mat[0] = NEG_INFINITY;
  lnbpp_mat[seq_len + 1] = NEG_INFINITY; */
  for i in 0 .. seq_len + 2 {
    let mut eps_of_terms_4_log_prob = EpsOfTerms4LogProb::new();
    let mut min_ep_of_term_4_log_prob = INFINITY;
    for j in 0 .. seq_len + 2 {
      if j == i {continue;}
      let pos_pair = if j < i {(j, i)} else {(i, j)};
      if lbpp_mat.contains_key(&pos_pair) {
        let ep_of_term_4_log_prob = lbpp_mat[&pos_pair];
        if ep_of_term_4_log_prob.is_finite() {
          eps_of_terms_4_log_prob.push(ep_of_term_4_log_prob);
          if ep_of_term_4_log_prob < min_ep_of_term_4_log_prob {
            min_ep_of_term_4_log_prob = ep_of_term_4_log_prob;
          }
        }
      }
    }
    lnbpp_mat[i] = (1. - logsumexp(&eps_of_terms_4_log_prob[..], min_ep_of_term_4_log_prob).exp()).ln();
  }
  lnbpp_mat
}

#[inline]
pub fn remove_little_lbpps_from_lbpp_mat(lbpp_mat: &SparseLogProbMat, min_lbpp: LogProb) -> SparseLogProbMat {
  lbpp_mat.iter().filter(|&(_, &lbpp)| {lbpp >= min_lbpp}).map(|(pos_pair, &lbpp)| {(*pos_pair, lbpp)}).collect::<SparseLogProbMat>()
}

#[inline]
pub fn run_command(command: &str, args: &[&str], expect: &str) -> Output {
  Command::new(command).args(args).output().expect(expect)
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
fn get_min_gap_num(pos_quadruple: &PosQuadruple) -> usize {
  let substr_len_pair = ((pos_quadruple.1 as isize - pos_quadruple.0 as isize).abs(), (pos_quadruple.3 as isize - pos_quadruple.2 as isize).abs());
  (substr_len_pair.0 - substr_len_pair.1).abs() as usize
}
