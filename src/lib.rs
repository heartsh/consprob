extern crate bio;
extern crate getopts;
extern crate hashbrown;
extern crate itertools;
extern crate ndarray;
extern crate num_cpus;
extern crate rna_algos;
extern crate scoped_threadpool;

pub use bio::io::fasta::Reader;
pub use getopts::Options;
pub use hashbrown::HashSet;
pub use itertools::multizip;
pub use ndarray::prelude::*;
pub use rna_algos::compiled_align_scores::*;
pub use rna_algos::durbin_algo::*;
pub use rna_algos::mccaskill_algo::*;
pub use rna_algos::utils::*;
pub use scoped_threadpool::Pool;
pub use std::cmp::Ord;
pub use std::fs::create_dir;
pub use std::fs::File;
pub use std::io::prelude::*;
pub use std::io::BufWriter;
pub use std::marker::{Send, Sync};
pub use std::path::Path;
pub use std::str::from_utf8_unchecked;

pub type ContextProfs = Array2<Prob>;
pub type ContextProf = Array1<Prob>;
pub type ContextProfSetPair = (ContextProfs, ContextProfs);
pub type PosQuadMat<T> = HashSet<PosQuad<T>>;
pub type PosPairMatSet<T> = HashMap<PosPair<T>, PosPairMat<T>>;
pub type PosPairMat<T> = HashSet<PosPair<T>>;
pub type ProbMat4d<T> = HashMap<PosQuad<T>, Prob>;
pub type SumMat4d<T> = HashMap<PosQuad<T>, Sum>;
pub type LoopSumsMat<T> = HashMap<PosPair<T>, LoopSums>;

#[derive(Clone)]
pub struct LoopSums {
  pub sum_seqalign: Sum,
  pub sum_seqalign_multibranch: Sum,
  pub sum_multibranch: Sum,
  pub sum_1st_pairmatches: Sum,
  pub sum_1ormore_pairmatches: Sum,
  pub sum_0ormore_pairmatches: Sum,
}

#[derive(Clone)]
pub struct AlignfoldSums<T> {
  pub sums_close: SumMat4d<T>,
  pub sums_accessible_external: SumMat4d<T>,
  pub sums_accessible_multibranch: SumMat4d<T>,
  pub forward_sums_external: SparseSumMat<T>,
  pub backward_sums_external: SparseSumMat<T>,
  pub forward_sums_external2: SparseSumMat<T>,
  pub backward_sums_external2: SparseSumMat<T>,
  pub forward_sums_hashed_poss: LoopSumsHashedPoss<T>,
  pub backward_sums_hashed_poss: LoopSumsHashedPoss<T>,
  pub forward_sums_hashed_poss2: LoopSumsHashedPoss<T>,
  pub backward_sums_hashed_poss2: LoopSumsHashedPoss<T>,
}

pub type ScoreMat = Vec<Vec<Score>>;

pub struct AlignfoldScores<T> {
  pub loopmatch_scores: SparseScoreMat<T>,
  pub pairmatch_scores: ScoreMat4d<T>,
  pub range_insert_scores: ScoreMat,
  pub range_insert_scores2: ScoreMat,
}

pub type ParamSetsHashedIds<T> = HashMap<RnaIdPair, AlignfoldScores<T>>;
pub type Probs4dHashedIds<T> = HashMap<RnaIdPair, ProbMat4d<T>>;
pub type ProbMats<T> = Vec<SparseProbMat<T>>;
pub type Prob1dMats = Vec<Probs>;
pub type Arg = String;
pub type Args = Vec<Arg>;
pub type FastaId = String;
pub type SeqPair<'a> = (SeqSlice<'a>, SeqSlice<'a>);
pub type SeqSlices<'a> = Vec<SeqSlice<'a>>;
pub type ScorePair = (Score, Score);
pub type SparseScoreMat<T> = HashMap<PosPair<T>, Score>;
pub type PosPairsHashedPoss<T> = HashMap<PosPair<T>, PosPair<T>>;
pub type BoolsHashedPoss<T> = HashMap<PosPair<T>, bool>;
pub type ProbMatPair<'a, T> = (&'a SparseProbMat<T>, &'a SparseProbMat<T>);
pub type FoldScoreSetPair<'a, T> = (&'a FoldScores<T>, &'a FoldScores<T>);
pub type NumThreads = u32;

pub struct AlignfoldProbMats<T> {
  pub basepair_probs_pair: SparseProbMatPair<T>,
  pub context_profs_pair: ContextProfSetPair,
  pub pairmatch_probs: ProbMat4d<T>,
  pub loopmatch_probs: SparseProbMat<T>,
  pub match_probs: SparseProbMat<T>,
}

#[derive(Clone)]
pub struct AlignfoldProbMatsAvg<T> {
  pub basepair_probs: SparseProbMat<T>,
  pub context_profs: ContextProfs,
}

pub type SparseProbMatPair<T> = (SparseProbMat<T>, SparseProbMat<T>);
pub type ProbMatSetsAvg<T> = Vec<AlignfoldProbMatsAvg<T>>;
pub type AlignfoldProbsHashedIds<T> = HashMap<RnaIdPair, AlignfoldProbMats<T>>;
pub type Poss<T> = Vec<T>;
pub type LoopSumsHashedPoss<T> = HashMap<PosPair<T>, LoopSumsMat<T>>;

#[derive(Clone)]
pub struct MatchProbMats<T> {
  pub loopmatch_probs: SparseProbMat<T>,
  pub pairmatch_probs: ProbMat4d<T>,
  pub match_probs: SparseProbMat<T>,
}

pub type MatchProbsHashedIds<T> = HashMap<RnaIdPair, MatchProbMats<T>>;
pub type SparseProbsHashedIds<T> = HashMap<RnaIdPair, SparseProbMat<T>>;
pub type PosQuadsHashedLens<T> = HashMap<PosPair<T>, PosPairMat<T>>;
pub type SparsePoss<T> = HashSet<T>;
pub type SparsePosSets<T> = HashMap<T, SparsePoss<T>>;

pub type OutputsSparsePossGetter<T> = (
  PosPairMatSet<T>,
  PosPairMatSet<T>,
  PosQuadMat<T>,
  PosQuadsHashedLens<T>,
  SparsePosSets<T>,
  SparsePosSets<T>,
);

pub type InputsAlignfoldProbsGetter<'a, T> = (
  &'a PosPair<T>,
  &'a AlignfoldScores<T>,
  &'a PosPair<T>,
  &'a AlignfoldSums<T>,
  &'a FoldScoreSetPair<'a, T>,
  bool,
  Sum,
  &'a PosQuadsHashedLens<T>,
  bool,
  &'a PosPairMatSet<T>,
  &'a PosPairMatSet<T>,
  &'a SparsePosSets<T>,
  &'a SparsePosSets<T>,
  &'a AlignScores,
);

pub type InputsConsprobCore<'a, T> = (
  &'a PosPair<T>,
  &'a AlignfoldScores<T>,
  &'a PosPair<T>,
  &'a FoldScoreSetPair<'a, T>,
  bool,
  &'a PosPairMatSet<T>,
  &'a PosPairMatSet<T>,
  &'a PosQuadsHashedLens<T>,
  bool,
  &'a SparsePosSets<T>,
  &'a SparsePosSets<T>,
  &'a AlignScores,
);

pub type Inputs2loopSumsGetter<'a, T> = (
  &'a PosQuad<T>,
  &'a AlignfoldSums<T>,
  bool,
  &'a PosPairMatSet<T>,
  &'a FoldScoreSetPair<'a, T>,
  &'a AlignfoldScores<T>,
  &'a SparsePosSets<T>,
  &'a SparsePosSets<T>,
  &'a AlignScores,
);

pub type InputsLoopSumsGetter<'a, T> = (
  &'a AlignfoldScores<T>,
  &'a PosQuad<T>,
  &'a mut AlignfoldSums<T>,
  bool,
  &'a PosPairMatSet<T>,
  &'a SparsePosSets<T>,
  &'a SparsePosSets<T>,
  &'a AlignScores,
);

pub type InputsInsideSumsGetter<'a, T> = (
  &'a PosPair<T>,
  &'a AlignfoldScores<T>,
  &'a PosPair<T>,
  &'a FoldScoreSetPair<'a, T>,
  &'a PosPairMatSet<T>,
  &'a PosPairMatSet<T>,
  &'a PosQuadsHashedLens<T>,
  &'a SparsePosSets<T>,
  &'a SparsePosSets<T>,
  &'a AlignScores,
);

impl<T: HashIndex> AlignfoldProbMats<T> {
  pub fn origin() -> AlignfoldProbMats<T> {
    AlignfoldProbMats {
      basepair_probs_pair: (SparseProbMat::<T>::default(), SparseProbMat::<T>::default()),
      context_profs_pair: (ContextProfs::default((0, 0)), ContextProfs::default((0, 0))),
      pairmatch_probs: ProbMat4d::<T>::default(),
      loopmatch_probs: SparseProbMat::<T>::default(),
      match_probs: SparseProbMat::<T>::default(),
    }
  }

  pub fn new(seq_len_pair: &(usize, usize)) -> AlignfoldProbMats<T> {
    let neg_infs_pair = (
      NEG_INFINITY * ContextProfs::ones((seq_len_pair.0, NUM_CONTEXTS)),
      NEG_INFINITY * ContextProfs::ones((seq_len_pair.1, NUM_CONTEXTS)),
    );
    AlignfoldProbMats {
      basepair_probs_pair: (SparseProbMat::<T>::default(), SparseProbMat::<T>::default()),
      context_profs_pair: neg_infs_pair,
      pairmatch_probs: ProbMat4d::<T>::default(),
      loopmatch_probs: SparseProbMat::<T>::default(),
      match_probs: SparseProbMat::<T>::default(),
    }
  }
}

impl<T: HashIndex> AlignfoldProbMatsAvg<T> {
  pub fn origin() -> AlignfoldProbMatsAvg<T> {
    AlignfoldProbMatsAvg {
      basepair_probs: SparseProbMat::<T>::default(),
      context_profs: ContextProfs::default((0, 0)),
    }
  }

  pub fn new(seq_len: usize) -> AlignfoldProbMatsAvg<T> {
    AlignfoldProbMatsAvg {
      basepair_probs: SparseProbMat::<T>::default(),
      context_profs: ContextProfs::zeros((seq_len, NUM_CONTEXTS)),
    }
  }
}

impl<T: HashIndex> AlignfoldScores<T> {
  pub fn origin() -> AlignfoldScores<T> {
    let scores = SparseScoreMat::<T>::default();
    let score_mat = ScoreMat4d::<T>::default();
    let insert_scores = Vec::new();
    AlignfoldScores {
      loopmatch_scores: scores,
      pairmatch_scores: score_mat,
      range_insert_scores: insert_scores.clone(),
      range_insert_scores2: insert_scores,
    }
  }

  pub fn new(
    seq_pair: &SeqPair,
    pos_quads: &PosQuadMat<T>,
    match_probs: &SparseProbMat<T>,
    align_scores: &AlignScores,
  ) -> AlignfoldScores<T> {
    let seq_len_pair = (
      T::from_usize(seq_pair.0.len()).unwrap(),
      T::from_usize(seq_pair.1.len()).unwrap(),
    );
    let mut alignfold_scores = AlignfoldScores::<T>::origin();
    let mat = vec![
      vec![NEG_INFINITY; seq_len_pair.0.to_usize().unwrap()];
      seq_len_pair.0.to_usize().unwrap()
    ];
    alignfold_scores.range_insert_scores = mat;
    let mat = vec![
      vec![NEG_INFINITY; seq_len_pair.1.to_usize().unwrap()];
      seq_len_pair.1.to_usize().unwrap()
    ];
    alignfold_scores.range_insert_scores2 = mat;
    for i in range(T::one(), seq_len_pair.1 - T::one()) {
      let long_i = i.to_usize().unwrap();
      let base = seq_pair.1[long_i];
      let mut sum = align_scores.insert_scores[base];
      alignfold_scores.range_insert_scores2[long_i][long_i] = sum;
      for j in range(i + T::one(), seq_len_pair.1 - T::one()) {
        let long_j = j.to_usize().unwrap();
        let base = seq_pair.1[long_j];
        let term = align_scores.insert_scores[base] + align_scores.insert_extend_score;
        sum += term;
        alignfold_scores.range_insert_scores2[long_i][long_j] = sum;
      }
    }
    for i in range(T::one(), seq_len_pair.0 - T::one()) {
      let long_i = i.to_usize().unwrap();
      let base = seq_pair.0[long_i];
      let term = align_scores.insert_scores[base];
      let mut sum = term;
      alignfold_scores.range_insert_scores[long_i][long_i] = sum;
      for j in range(i + T::one(), seq_len_pair.0 - T::one()) {
        let long_j = j.to_usize().unwrap();
        let base = seq_pair.0[long_j];
        let term = align_scores.insert_scores[base] + align_scores.insert_extend_score;
        sum += term;
        alignfold_scores.range_insert_scores[long_i][long_j] = sum;
      }
    }
    for pos_pair in match_probs.keys() {
      let &(i, j) = pos_pair;
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let base_pair = (seq_pair.0[long_i], seq_pair.1[long_j]);
      alignfold_scores.loopmatch_scores.insert(
        *pos_pair,
        align_scores.match_scores[base_pair.0][base_pair.1],
      );
    }
    for &(i, j, k, l) in pos_quads {
      let (long_i, long_j, long_k, long_l) = (
        i.to_usize().unwrap(),
        j.to_usize().unwrap(),
        k.to_usize().unwrap(),
        l.to_usize().unwrap(),
      );
      let pos_quad = (i, j, k, l);
      let base_pair = (seq_pair.0[long_i], seq_pair.0[long_j]);
      let base_pair2 = (seq_pair.1[long_k], seq_pair.1[long_l]);
      alignfold_scores.pairmatch_scores.insert(
        pos_quad,
        align_scores.match_scores[base_pair.0][base_pair2.0]
          + align_scores.match_scores[base_pair.1][base_pair2.1],
      );
    }
    alignfold_scores
  }
}

impl Default for LoopSums {
  fn default() -> Self {
    Self::new()
  }
}

impl LoopSums {
  pub fn new() -> LoopSums {
    LoopSums {
      sum_seqalign: NEG_INFINITY,
      sum_seqalign_multibranch: NEG_INFINITY,
      sum_multibranch: NEG_INFINITY,
      sum_1st_pairmatches: NEG_INFINITY,
      sum_1ormore_pairmatches: NEG_INFINITY,
      sum_0ormore_pairmatches: NEG_INFINITY,
    }
  }
}

impl<T: HashIndex> Default for AlignfoldSums<T> {
  fn default() -> Self {
    Self::new()
  }
}

impl<T: HashIndex> AlignfoldSums<T> {
  pub fn new() -> AlignfoldSums<T> {
    let sum_mat = SumMat4d::<T>::default();
    let sums = SparseSumMat::<T>::default();
    let loop_sums_hashed_poss = LoopSumsHashedPoss::<T>::default();
    AlignfoldSums {
      sums_close: sum_mat.clone(),
      sums_accessible_external: sum_mat.clone(),
      sums_accessible_multibranch: sum_mat,
      forward_sums_external: sums.clone(),
      backward_sums_external: sums.clone(),
      forward_sums_external2: sums.clone(),
      backward_sums_external2: sums,
      forward_sums_hashed_poss: loop_sums_hashed_poss.clone(),
      backward_sums_hashed_poss: loop_sums_hashed_poss.clone(),
      forward_sums_hashed_poss2: loop_sums_hashed_poss.clone(),
      backward_sums_hashed_poss2: loop_sums_hashed_poss,
    }
  }
}

impl<T> Default for MatchProbMats<T> {
  fn default() -> Self {
    Self::new()
  }
}

impl<T> MatchProbMats<T> {
  pub fn new() -> MatchProbMats<T> {
    MatchProbMats {
      loopmatch_probs: SparseProbMat::<T>::default(),
      pairmatch_probs: ProbMat4d::<T>::default(),
      match_probs: SparseProbMat::<T>::default(),
    }
  }
}

pub const DEFAULT_MIN_BASEPAIR_PROB: Prob = 0.01;
pub const DEFAULT_MIN_MATCH_PROB: Prob = 0.01;
pub const BASEPAIR_PROBS_FILE: &str = "basepair_probs.dat";
pub const UNPAIR_PROBS_FILE_HAIRPIN: &str = "unpair_probs_hairpin.dat";
pub const UNPAIR_PROBS_FILE_BULGE: &str = "unpair_probs_bulge.dat";
pub const UNPAIR_PROBS_FILE_INTERIOR: &str = "unpair_probs_interior.dat";
pub const UNPAIR_PROBS_FILE_MULTIBRANCH: &str = "unpair_probs_multibranch.dat";
pub const UNPAIR_PROBS_FILE_EXTERNAL: &str = "unpair_probs_external.dat";
pub const BASEPAIR_PROBS_FILE2: &str = "basepair_probs2.dat";
pub const PAIRMATCH_PROBS_FILE: &str = "pairmatch_probs.dat";
pub const LOOPMATCH_PROBS_FILE: &str = "loopmatch_probs.dat";
pub const MATCH_PROBS_FILE: &str = "match_probs.dat";
pub const README_FILE: &str = "README.md";
pub const README_CONTENTS: &str = "# basepair_probs.dat\n
This file contains average probabilistic consistency based on posterior nucleotide pair-matching probabilities. You can treat this average probabilistic consistency like conventional nucleotide base-pairing probabilities. Nucleotide positions are indexed starting from zero.\n\n
# basepair_probs2.dat\n
This file contains average probabilistic consistency per nucleotide. This average probabilistic consistency is obtained by marginalizing each nucleotide for average probabilistic consistency in \"basepair_probs.dat.\"\n\n
# unpair_probs_x.dat\n
This file type contains average probabilistic consistency per nucleotide. This average probabilistic consistency is for nucleotide unpairing and under the structural context \"x.\" \"hairpin,\" \"bulge,\" \"interior,\" \"multibranch,\" \"external\" stand for hairpin loops, bulge loops, interior loops, multi-loops, external loops, respectively.\n\n
# pairmatch_probs.dat\n
This file contains posterior nucleotide pair-matching probabilities.\n\n
# loopmatch_probs.dat\n
This file contains posterior nucleotide loop-matching probabilities.\n\n
match_probs.dat\n
This file contains posterior nucleotide matching probabilities.";
pub const INSERT2MATCH_SCORE: Score = MATCH2INSERT_SCORE;
pub const NUM_CONTEXTS: usize = 6;
pub const CONTEXT_INDEX_HAIRPIN: usize = 0;
pub const CONTEXT_INDEX_BULGE: usize = 1;
pub const CONTEXT_INDEX_INTERIOR: usize = 2;
pub const CONTEXT_INDEX_EXTERNAL: usize = 3;
pub const CONTEXT_INDEX_BASEPAIR: usize = 4;
pub const CONTEXT_INDEX_MULTIBRANCH: usize = 5;
pub const EXAMPLE_FASTA_FILE_PATH: &str = "assets/sampled_trnas.fa";
pub const EPSILON: Prob = 0.00_1;
pub const PROB_BOUND_LOWER: Prob = -EPSILON;
pub const PROB_BOUND_UPPER: Prob = 1. + EPSILON;

pub fn consprob_core<T>(inputs: InputsConsprobCore<T>) -> AlignfoldProbMats<T>
where
  T: HashIndex,
{
  let (
    seq_len_pair,
    alignfold_scores,
    max_basepair_span_pair,
    fold_scores_pair,
    produces_context_profs,
    forward_pos_pairs,
    backward_pos_pairs,
    pos_quads_hashed_lens,
    produces_match_probs,
    matchable_poss,
    matchable_poss2,
    align_scores,
  ) = inputs;
  let (alignfold_sums, global_sum) = get_alignfold_sums::<T>((
    seq_len_pair,
    alignfold_scores,
    max_basepair_span_pair,
    fold_scores_pair,
    forward_pos_pairs,
    backward_pos_pairs,
    pos_quads_hashed_lens,
    matchable_poss,
    matchable_poss2,
    align_scores,
  ));
  get_alignfold_probs::<T>((
    seq_len_pair,
    alignfold_scores,
    max_basepair_span_pair,
    &alignfold_sums,
    fold_scores_pair,
    produces_context_profs,
    global_sum,
    pos_quads_hashed_lens,
    produces_match_probs,
    forward_pos_pairs,
    backward_pos_pairs,
    matchable_poss,
    matchable_poss2,
    align_scores,
  ))
}

pub fn get_alignfold_sums<T>(inputs: InputsInsideSumsGetter<T>) -> (AlignfoldSums<T>, Sum)
where
  T: HashIndex,
{
  let (
    seq_len_pair,
    alignfold_scores,
    max_basepair_span_pair,
    fold_scores_pair,
    forward_pos_pairs,
    backward_pos_pairs,
    pos_quads_hashed_lens,
    matchable_poss,
    matchable_poss2,
    align_scores,
  ) = inputs;
  let mut alignfold_sums = AlignfoldSums::<T>::new();
  for substr_len in range_inclusive(
    T::from_usize(MIN_SPAN_HAIRPIN_CLOSE).unwrap(),
    max_basepair_span_pair.0,
  ) {
    for substr_len2 in range_inclusive(
      T::from_usize(MIN_SPAN_HAIRPIN_CLOSE).unwrap(),
      max_basepair_span_pair.1,
    ) {
      if let Some(pos_pairs) = pos_quads_hashed_lens.get(&(substr_len, substr_len2)) {
        for &(i, k) in pos_pairs {
          let (j, l) = (i + substr_len - T::one(), k + substr_len2 - T::one());
          let (long_i, long_j, long_k, long_l) = (
            i.to_usize().unwrap(),
            j.to_usize().unwrap(),
            k.to_usize().unwrap(),
            l.to_usize().unwrap(),
          );
          let pos_quad = (i, j, k, l);
          let pairmatch_score = alignfold_scores.pairmatch_scores[&pos_quad];
          let computes_forward_sums = true;
          let (sum_seqalign, sum_multibranch) = get_loop_sums::<T>((
            alignfold_scores,
            &pos_quad,
            &mut alignfold_sums,
            computes_forward_sums,
            forward_pos_pairs,
            matchable_poss,
            matchable_poss2,
            align_scores,
          ));
          let computes_forward_sums = false;
          let _ = get_loop_sums::<T>((
            alignfold_scores,
            &pos_quad,
            &mut alignfold_sums,
            computes_forward_sums,
            backward_pos_pairs,
            matchable_poss,
            matchable_poss2,
            align_scores,
          ));
          let mut sum = NEG_INFINITY;
          let score = pairmatch_score
            + fold_scores_pair.0.hairpin_scores[&(i, j)]
            + fold_scores_pair.1.hairpin_scores[&(k, l)]
            + sum_seqalign;
          logsumexp(&mut sum, score);
          let forward_sums = &alignfold_sums.forward_sums_hashed_poss2[&(i, k)];
          let backward_sums = &alignfold_sums.backward_sums_hashed_poss2[&(j, l)];
          let min = T::from_usize(MIN_SPAN_HAIRPIN_CLOSE).unwrap();
          let min_len_pair = (
            if substr_len <= min + T::from_usize(MAX_2LOOP_LEN + 2).unwrap() {
              min
            } else {
              substr_len - T::from_usize(MAX_2LOOP_LEN + 2).unwrap()
            },
            if substr_len2 <= min + T::from_usize(MAX_2LOOP_LEN + 2).unwrap() {
              min
            } else {
              substr_len2 - T::from_usize(MAX_2LOOP_LEN + 2).unwrap()
            },
          );
          for substr_len3 in range(min_len_pair.0, substr_len - T::one()) {
            for substr_len4 in range(min_len_pair.1, substr_len2 - T::one()) {
              if let Some(pos_pairs2) = pos_quads_hashed_lens.get(&(substr_len3, substr_len4)) {
                for &(m, o) in pos_pairs2 {
                  let (n, p) = (m + substr_len3 - T::one(), o + substr_len4 - T::one());
                  if !(i < m && n < j && k < o && p < l) {
                    continue;
                  }
                  let (long_m, long_n, long_o, long_p) = (
                    m.to_usize().unwrap(),
                    n.to_usize().unwrap(),
                    o.to_usize().unwrap(),
                    p.to_usize().unwrap(),
                  );
                  if long_m - long_i - 1 + long_j - long_n - 1 > MAX_2LOOP_LEN {
                    continue;
                  }
                  if long_o - long_k - 1 + long_l - long_p - 1 > MAX_2LOOP_LEN {
                    continue;
                  }
                  let pos_quad2 = (m, n, o, p);
                  if let Some(&x) = alignfold_sums.sums_close.get(&pos_quad2) {
                    let mut forward_term = NEG_INFINITY;
                    let mut backward_term = forward_term;
                    let pos_pair2 = (m - T::one(), o - T::one());
                    if let Some(x) = forward_sums.get(&pos_pair2) {
                      logsumexp(&mut forward_term, x.sum_seqalign);
                    }
                    let pos_pair2 = (n + T::one(), p + T::one());
                    if let Some(x) = backward_sums.get(&pos_pair2) {
                      logsumexp(&mut backward_term, x.sum_seqalign);
                    }
                    let sum_2loop = forward_term + backward_term;
                    let twoloop_score = fold_scores_pair.0.twoloop_scores[&(i, j, m, n)];
                    let twoloop_score2 = fold_scores_pair.1.twoloop_scores[&(k, l, o, p)];
                    let x = pairmatch_score + twoloop_score + twoloop_score2 + x;
                    logsumexp(&mut sum, x + sum_2loop);
                  }
                }
              }
            }
          }
          let multibranch_close_score = fold_scores_pair.0.multibranch_close_scores[&(i, j)];
          let multibranch_close_score2 = fold_scores_pair.1.multibranch_close_scores[&(k, l)];
          let score =
            pairmatch_score + multibranch_close_score + multibranch_close_score2 + sum_multibranch;
          logsumexp(&mut sum, score);
          if sum > NEG_INFINITY {
            alignfold_sums.sums_close.insert(pos_quad, sum);
            let accessible_score = fold_scores_pair.0.accessible_scores[&(i, j)];
            let accessible_score2 = fold_scores_pair.1.accessible_scores[&(k, l)];
            sum += accessible_score + accessible_score2;
            alignfold_sums
              .sums_accessible_external
              .insert(pos_quad, sum);
            alignfold_sums
              .sums_accessible_multibranch
              .insert(pos_quad, sum + 2. * COEFF_NUM_BRANCHES);
          }
        }
      }
    }
  }
  let leftmost_pos_pair = (T::zero(), T::zero());
  let rightmost_pos_pair = (seq_len_pair.0 - T::one(), seq_len_pair.1 - T::one());
  alignfold_sums
    .forward_sums_external
    .insert(leftmost_pos_pair, 0.);
  alignfold_sums
    .backward_sums_external
    .insert(rightmost_pos_pair, 0.);
  for i in range(T::zero(), seq_len_pair.0 - T::one()) {
    for j in range(T::zero(), seq_len_pair.1 - T::one()) {
      let pos_pair = (i, j);
      if pos_pair == leftmost_pos_pair {
        continue;
      }
      let mut sum = NEG_INFINITY;
      if let Some(x) = forward_pos_pairs.get(&pos_pair) {
        for &(k, l) in x {
          let pos_pair2 = (k - T::one(), l - T::one());
          let pos_quad = (k, i, l, j);
          if let Some(&x) = alignfold_sums.sums_accessible_external.get(&pos_quad) {
            if let Some(&y) = alignfold_sums.forward_sums_external2.get(&pos_pair2) {
              let y = x + y;
              logsumexp(&mut sum, y);
            }
          }
        }
      }
      if i > T::zero() && j > T::zero() {
        if let Some(&loopmatch_score) = alignfold_scores.loopmatch_scores.get(&pos_pair) {
          let mut sum2 = NEG_INFINITY;
          let pos_pair2 = (i - T::one(), j - T::one());
          let long_pos_pair2 = (
            pos_pair2.0.to_usize().unwrap(),
            pos_pair2.1.to_usize().unwrap(),
          );
          let begins_sum = pos_pair2 == leftmost_pos_pair;
          if let Some(&x) = alignfold_sums.forward_sums_external.get(&pos_pair2) {
            let x = x
              + if begins_sum {
                align_scores.init_match_score
              } else {
                align_scores.match2match_score
              };
            logsumexp(&mut sum2, x);
          }
          if let Some(x) = matchable_poss.get(&pos_pair2.0) {
            for &x in x {
              if x >= pos_pair2.1 {
                continue;
              }
              let pos_pair3 = (pos_pair2.0, x);
              if let Some(&y) = alignfold_sums.forward_sums_external.get(&pos_pair3) {
                let long_x = x.to_usize().unwrap();
                let begins_sum = pos_pair3 == leftmost_pos_pair;
                let z = alignfold_scores.range_insert_scores2[long_x + 1][long_pos_pair2.1]
                  + if begins_sum {
                    align_scores.init_insert_score
                  } else {
                    align_scores.match2insert_score
                  };
                let z = y + z + align_scores.match2insert_score;
                logsumexp(&mut sum2, z);
              }
            }
          }
          if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
            for &x in x {
              if x >= pos_pair2.0 {
                continue;
              }
              let pos_pair3 = (x, pos_pair2.1);
              if let Some(&y) = alignfold_sums.forward_sums_external.get(&pos_pair3) {
                let long_x = x.to_usize().unwrap();
                let begins_sum = pos_pair3 == leftmost_pos_pair;
                let z = alignfold_scores.range_insert_scores[long_x + 1][long_pos_pair2.0]
                  + if begins_sum {
                    align_scores.init_insert_score
                  } else {
                    align_scores.match2insert_score
                  };
                let z = y + z + align_scores.match2insert_score;
                logsumexp(&mut sum2, z);
              }
            }
          }
          if sum2 > NEG_INFINITY {
            alignfold_sums
              .forward_sums_external2
              .insert(pos_pair2, sum2);
          }
          let term = sum2 + loopmatch_score;
          logsumexp(&mut sum, term);
          if sum > NEG_INFINITY {
            alignfold_sums.forward_sums_external.insert(pos_pair, sum);
          }
        }
      }
    }
  }
  let mut final_sum = NEG_INFINITY;
  let pos_pair2 = (
    rightmost_pos_pair.0 - T::one(),
    rightmost_pos_pair.1 - T::one(),
  );
  let long_pos_pair2 = (
    pos_pair2.0.to_usize().unwrap(),
    pos_pair2.1.to_usize().unwrap(),
  );
  if let Some(&x) = alignfold_sums.forward_sums_external.get(&pos_pair2) {
    logsumexp(&mut final_sum, x);
  }
  if let Some(x) = matchable_poss.get(&pos_pair2.0) {
    for &x in x {
      if x >= pos_pair2.1 {
        continue;
      }
      if let Some(&y) = alignfold_sums.forward_sums_external.get(&(pos_pair2.0, x)) {
        let long_x = x.to_usize().unwrap();
        let z = alignfold_scores.range_insert_scores2[long_x + 1][long_pos_pair2.1]
          + align_scores.match2insert_score;
        let z = y + z;
        logsumexp(&mut final_sum, z);
      }
    }
  }
  if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
    for &x in x {
      if x >= pos_pair2.0 {
        continue;
      }
      if let Some(&y) = alignfold_sums.forward_sums_external.get(&(x, pos_pair2.1)) {
        let long_x = x.to_usize().unwrap();
        let z = alignfold_scores.range_insert_scores[long_x + 1][long_pos_pair2.0]
          + align_scores.match2insert_score;
        let z = y + z;
        logsumexp(&mut final_sum, z);
      }
    }
  }
  for i in range(T::one(), seq_len_pair.0).rev() {
    for j in range(T::one(), seq_len_pair.1).rev() {
      let pos_pair = (i, j);
      if pos_pair == rightmost_pos_pair {
        continue;
      }
      let mut sum = NEG_INFINITY;
      if let Some(x) = backward_pos_pairs.get(&pos_pair) {
        for &(k, l) in x {
          let pos_pair2 = (k + T::one(), l + T::one());
          let pos_quad = (i, k, j, l);
          if let Some(&x) = alignfold_sums.sums_accessible_external.get(&pos_quad) {
            if let Some(&y) = alignfold_sums.backward_sums_external2.get(&pos_pair2) {
              let y = x + y;
              logsumexp(&mut sum, y);
            }
          }
        }
      }
      if i < seq_len_pair.0 - T::one() && j < seq_len_pair.1 - T::one() {
        if let Some(&loopmatch_score) = alignfold_scores.loopmatch_scores.get(&pos_pair) {
          let mut sum2 = NEG_INFINITY;
          let pos_pair2 = (i + T::one(), j + T::one());
          let long_pos_pair2 = (
            pos_pair2.0.to_usize().unwrap(),
            pos_pair2.1.to_usize().unwrap(),
          );
          let ends_sum = pos_pair2 == rightmost_pos_pair;
          if let Some(&x) = alignfold_sums.backward_sums_external.get(&pos_pair2) {
            let x = x
              + if ends_sum {
                0.
              } else {
                align_scores.match2match_score
              };
            logsumexp(&mut sum2, x);
          }
          if let Some(x) = matchable_poss.get(&pos_pair2.0) {
            for &x in x {
              if x <= pos_pair2.1 {
                continue;
              }
              let pos_pair3 = (pos_pair2.0, x);
              if let Some(&y) = alignfold_sums.backward_sums_external.get(&pos_pair3) {
                let long_x = x.to_usize().unwrap();
                let ends_sum = pos_pair3 == rightmost_pos_pair;
                let z = alignfold_scores.range_insert_scores2[long_pos_pair2.1][long_x - 1]
                  + if ends_sum {
                    0.
                  } else {
                    align_scores.match2insert_score
                  };
                let z = y + z + align_scores.match2insert_score;
                logsumexp(&mut sum2, z);
              }
            }
          }
          if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
            for &x in x {
              if x <= pos_pair2.0 {
                continue;
              }
              let pos_pair3 = (x, pos_pair2.1);
              if let Some(&y) = alignfold_sums.backward_sums_external.get(&pos_pair3) {
                let long_x = x.to_usize().unwrap();
                let ends_sum = pos_pair3 == rightmost_pos_pair;
                let z = alignfold_scores.range_insert_scores[long_pos_pair2.0][long_x - 1]
                  + if ends_sum {
                    0.
                  } else {
                    align_scores.match2insert_score
                  };
                let z = y + z + align_scores.match2insert_score;
                logsumexp(&mut sum2, z);
              }
            }
          }
          if sum2 > NEG_INFINITY {
            alignfold_sums
              .backward_sums_external2
              .insert(pos_pair2, sum2);
          }
          let term = sum2 + loopmatch_score;
          logsumexp(&mut sum, term);
          if sum > NEG_INFINITY {
            alignfold_sums.backward_sums_external.insert(pos_pair, sum);
          }
        }
      }
    }
  }
  (alignfold_sums, final_sum)
}

pub fn get_loop_sums<T>(inputs: InputsLoopSumsGetter<T>) -> (Sum, Sum)
where
  T: HashIndex,
{
  let (
    alignfold_scores,
    pos_quad,
    alignfold_sums,
    computes_forward_sums,
    pos_pairs,
    matchable_poss,
    matchable_poss2,
    align_scores,
  ) = inputs;
  let &(i, j, k, l) = pos_quad;
  let leftmost_pos_pair = if computes_forward_sums {
    (i, k)
  } else {
    (i + T::one(), k + T::one())
  };
  let rightmost_pos_pair = if computes_forward_sums {
    (j - T::one(), l - T::one())
  } else {
    (j, l)
  };
  let sums_hashed_poss = if computes_forward_sums {
    &mut alignfold_sums.forward_sums_hashed_poss
  } else {
    &mut alignfold_sums.backward_sums_hashed_poss
  };
  let sums_hashed_poss2 = if computes_forward_sums {
    &mut alignfold_sums.forward_sums_hashed_poss2
  } else {
    &mut alignfold_sums.backward_sums_hashed_poss2
  };
  if !sums_hashed_poss.contains_key(&if computes_forward_sums {
    leftmost_pos_pair
  } else {
    rightmost_pos_pair
  }) {
    sums_hashed_poss.insert(
      if computes_forward_sums {
        leftmost_pos_pair
      } else {
        rightmost_pos_pair
      },
      LoopSumsMat::<T>::new(),
    );
  }
  if !sums_hashed_poss2.contains_key(&if computes_forward_sums {
    leftmost_pos_pair
  } else {
    rightmost_pos_pair
  }) {
    sums_hashed_poss2.insert(
      if computes_forward_sums {
        leftmost_pos_pair
      } else {
        rightmost_pos_pair
      },
      LoopSumsMat::<T>::new(),
    );
  }
  let sums_mat = &mut sums_hashed_poss
    .get_mut(&if computes_forward_sums {
      leftmost_pos_pair
    } else {
      rightmost_pos_pair
    })
    .unwrap();
  let sums_mat2 = &mut sums_hashed_poss2
    .get_mut(&if computes_forward_sums {
      leftmost_pos_pair
    } else {
      rightmost_pos_pair
    })
    .unwrap();
  let iter: Poss<T> = if computes_forward_sums {
    range(i, j).collect()
  } else {
    range_inclusive(i + T::one(), j).rev().collect()
  };
  let iter2: Poss<T> = if computes_forward_sums {
    range(k, l).collect()
  } else {
    range_inclusive(k + T::one(), l).rev().collect()
  };
  for &u in iter.iter() {
    for &v in iter2.iter() {
      let pos_pair = (u, v);
      if sums_mat.contains_key(&pos_pair) {
        continue;
      }
      let mut sums = LoopSums::new();
      if (computes_forward_sums && u == i && v == k) || (!computes_forward_sums && u == j && v == l)
      {
        sums.sum_seqalign = 0.;
        sums.sum_0ormore_pairmatches = 0.;
        sums_mat.insert(pos_pair, sums);
        continue;
      }
      let mut sum_multibranch = NEG_INFINITY;
      let mut sum_1st_pairmatches = sum_multibranch;
      let mut sum = sum_multibranch;
      if let Some(pos_pairs) = pos_pairs.get(&pos_pair) {
        for &(m, n) in pos_pairs {
          if computes_forward_sums {
            if !(i < m && k < n) {
              continue;
            }
          } else if !(m < j && n < l) {
            continue;
          }
          let pos_pair2 = if computes_forward_sums {
            (m - T::one(), n - T::one())
          } else {
            (m + T::one(), n + T::one())
          };
          let pos_quad2 = if computes_forward_sums {
            (m, u, n, v)
          } else {
            (u, m, v, n)
          };
          if let Some(&x) = alignfold_sums.sums_accessible_multibranch.get(&pos_quad2) {
            if let Some(y) = sums_mat2.get(&pos_pair2) {
              let z = x + y.sum_1ormore_pairmatches;
              logsumexp(&mut sum_multibranch, z);
              let z = x + y.sum_seqalign;
              logsumexp(&mut sum_1st_pairmatches, z);
            }
          }
        }
      }
      let pos_pair2 = if computes_forward_sums {
        (u - T::one(), v - T::one())
      } else {
        (u + T::one(), v + T::one())
      };
      let long_pos_pair2 = (
        pos_pair2.0.to_usize().unwrap(),
        pos_pair2.1.to_usize().unwrap(),
      );
      if let Some(&loopmatch_score) = alignfold_scores.loopmatch_scores.get(&pos_pair) {
        let mut sums2 = LoopSums::new();
        let mut sum_seqalign2 = NEG_INFINITY;
        let mut sum_multibranch2 = sum_seqalign2;
        let mut sum_1st_pairmatches2 = sum_seqalign2;
        let mut sum2 = sum_seqalign2;
        if let Some(x) = sums_mat.get(&pos_pair2) {
          let y = x.sum_multibranch + align_scores.match2match_score;
          logsumexp(&mut sum_multibranch2, y);
          let y = x.sum_1st_pairmatches + align_scores.match2match_score;
          logsumexp(&mut sum_1st_pairmatches2, y);
          let y = x.sum_seqalign + align_scores.match2match_score;
          logsumexp(&mut sum_seqalign2, y);
        }
        if let Some(x) = matchable_poss.get(&pos_pair2.0) {
          for &x in x {
            if computes_forward_sums && x >= pos_pair2.1
              || (!computes_forward_sums && x <= pos_pair2.1)
            {
              continue;
            }
            if let Some(y) = sums_mat.get(&(pos_pair2.0, x)) {
              let long_x = x.to_usize().unwrap();
              let z = if computes_forward_sums {
                alignfold_scores.range_insert_scores2[long_x + 1][long_pos_pair2.1]
              } else {
                alignfold_scores.range_insert_scores2[long_pos_pair2.1][long_x - 1]
              } + align_scores.match2insert_score;
              let x = y.sum_multibranch + align_scores.match2insert_score + z;
              logsumexp(&mut sum_multibranch2, x);
              let x = y.sum_1st_pairmatches + align_scores.match2insert_score + z;
              logsumexp(&mut sum_1st_pairmatches2, x);
              let x = y.sum_seqalign + align_scores.match2insert_score + z;
              logsumexp(&mut sum_seqalign2, x);
            }
          }
        }
        if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
          for &x in x {
            if computes_forward_sums && x >= pos_pair2.0
              || (!computes_forward_sums && x <= pos_pair2.0)
            {
              continue;
            }
            if let Some(y) = sums_mat.get(&(x, pos_pair2.1)) {
              let long_x = x.to_usize().unwrap();
              let z = if computes_forward_sums {
                alignfold_scores.range_insert_scores[long_x + 1][long_pos_pair2.0]
              } else {
                alignfold_scores.range_insert_scores[long_pos_pair2.0][long_x - 1]
              } + align_scores.match2insert_score;
              let x = y.sum_multibranch + align_scores.match2insert_score + z;
              logsumexp(&mut sum_multibranch2, x);
              let x = y.sum_1st_pairmatches + align_scores.match2insert_score + z;
              logsumexp(&mut sum_1st_pairmatches2, x);
              let x = y.sum_seqalign + align_scores.match2insert_score + z;
              logsumexp(&mut sum_seqalign2, x);
            }
          }
        }
        sums2.sum_multibranch = sum_multibranch2;
        logsumexp(&mut sum2, sum_multibranch2);
        sums2.sum_1st_pairmatches = sum_1st_pairmatches2;
        logsumexp(&mut sum2, sum_1st_pairmatches2);
        sums2.sum_1ormore_pairmatches = sum2;
        sums2.sum_seqalign = sum_seqalign2;
        logsumexp(&mut sum2, sum_seqalign2);
        sums2.sum_0ormore_pairmatches = sum2;
        if has_valid_sums(&sums2) {
          sums_mat2.insert(pos_pair2, sums2);
        }
        let term = sum_multibranch2 + loopmatch_score;
        logsumexp(&mut sum_multibranch, term);
        sums.sum_multibranch = sum_multibranch;
        logsumexp(&mut sum, sum_multibranch);
        let term = sum_1st_pairmatches2 + loopmatch_score;
        logsumexp(&mut sum_1st_pairmatches, term);
        sums.sum_1st_pairmatches = sum_1st_pairmatches;
        logsumexp(&mut sum, sum_1st_pairmatches);
        sums.sum_1ormore_pairmatches = sum;
        let sum_seqalign = sum_seqalign2 + loopmatch_score;
        sums.sum_seqalign = sum_seqalign;
        logsumexp(&mut sum, sum_seqalign);
        sums.sum_0ormore_pairmatches = sum;
        if has_valid_sums(&sums) {
          sums_mat.insert(pos_pair, sums);
        }
      }
    }
  }
  let mut final_sum_seqalign = NEG_INFINITY;
  let mut final_sum_multibranch = final_sum_seqalign;
  if computes_forward_sums {
    let pos_pair2 = rightmost_pos_pair;
    let long_pos_pair2 = (
      pos_pair2.0.to_usize().unwrap(),
      pos_pair2.1.to_usize().unwrap(),
    );
    if let Some(x) = sums_mat.get(&pos_pair2) {
      let y = x.sum_multibranch + align_scores.match2match_score;
      logsumexp(&mut final_sum_multibranch, y);
      let y = x.sum_seqalign + align_scores.match2match_score;
      logsumexp(&mut final_sum_seqalign, y);
    }
    if let Some(x) = matchable_poss.get(&pos_pair2.0) {
      for &x in x {
        if x >= pos_pair2.1 {
          continue;
        }
        if let Some(y) = sums_mat.get(&(pos_pair2.0, x)) {
          let long_x = x.to_usize().unwrap();
          let z = alignfold_scores.range_insert_scores2[long_x + 1][long_pos_pair2.1]
            + align_scores.match2insert_score;
          let x = y.sum_multibranch + align_scores.match2insert_score + z;
          logsumexp(&mut final_sum_multibranch, x);
          let x = y.sum_seqalign + align_scores.match2insert_score + z;
          logsumexp(&mut final_sum_seqalign, x);
        }
      }
    }
    if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
      for &x in x {
        if x >= pos_pair2.0 {
          continue;
        }
        if let Some(y) = sums_mat.get(&(x, pos_pair2.1)) {
          let long_x = x.to_usize().unwrap();
          let z = alignfold_scores.range_insert_scores[long_x + 1][long_pos_pair2.0]
            + align_scores.match2insert_score;
          let x = y.sum_multibranch + align_scores.match2insert_score + z;
          logsumexp(&mut final_sum_multibranch, x);
          let x = y.sum_seqalign + align_scores.match2insert_score + z;
          logsumexp(&mut final_sum_seqalign, x);
        }
      }
    }
  }
  (final_sum_seqalign, final_sum_multibranch)
}

pub fn get_2loop_sums<T>(inputs: Inputs2loopSumsGetter<T>) -> (SparseSumMat<T>, SparseSumMat<T>)
where
  T: HashIndex,
{
  let (
    pos_quad,
    alignfold_sums,
    computes_forward_sums,
    pos_pairs,
    fold_scores_pair,
    alignfold_scores,
    matchable_poss,
    matchable_poss2,
    align_scores,
  ) = inputs;
  let &(i, j, k, l) = pos_quad;
  let leftmost_pos_pair = if computes_forward_sums {
    (i, k)
  } else {
    (i + T::one(), k + T::one())
  };
  let rightmost_pos_pair = if computes_forward_sums {
    (j - T::one(), l - T::one())
  } else {
    (j, l)
  };
  let sums_hashed_poss = if computes_forward_sums {
    &alignfold_sums.forward_sums_hashed_poss2
  } else {
    &alignfold_sums.backward_sums_hashed_poss2
  };
  let sums = &sums_hashed_poss[&if computes_forward_sums {
    leftmost_pos_pair
  } else {
    rightmost_pos_pair
  }];
  let iter: Poss<T> = if computes_forward_sums {
    range(i, j).collect()
  } else {
    range_inclusive(i + T::one(), j).rev().collect()
  };
  let iter2: Poss<T> = if computes_forward_sums {
    range(k, l).collect()
  } else {
    range_inclusive(k + T::one(), l).rev().collect()
  };
  let mut sum_mat = SparseSumMat::<T>::default();
  let mut sum_mat2 = sum_mat.clone();
  for &u in iter.iter() {
    for &v in iter2.iter() {
      let pos_pair = (u, v);
      if (computes_forward_sums && u == i && v == k) || (!computes_forward_sums && u == j && v == l)
      {
        continue;
      }
      let mut sum = NEG_INFINITY;
      if let Some(x) = pos_pairs.get(&pos_pair) {
        for &(m, n) in x {
          if computes_forward_sums {
            if !(i < m && k < n) {
              continue;
            }
          } else if !(m < j && n < l) {
            continue;
          }
          let pos_pair2 = if computes_forward_sums {
            (m - T::one(), n - T::one())
          } else {
            (m + T::one(), n + T::one())
          };
          let pos_quad2 = if computes_forward_sums {
            (m, u, n, v)
          } else {
            (u, m, v, n)
          };
          if pos_quad2.0 - i - T::one() + j - pos_quad2.1 - T::one()
            > T::from_usize(MAX_2LOOP_LEN).unwrap()
          {
            continue;
          }
          if pos_quad2.2 - k - T::one() + l - pos_quad2.3 - T::one()
            > T::from_usize(MAX_2LOOP_LEN).unwrap()
          {
            continue;
          }
          if let Some(&x) = alignfold_sums.sums_close.get(&pos_quad2) {
            if let Some(y) = sums.get(&pos_pair2) {
              let twoloop_score =
                fold_scores_pair.0.twoloop_scores[&(i, j, pos_quad2.0, pos_quad2.1)];
              let twoloop_score2 =
                fold_scores_pair.1.twoloop_scores[&(k, l, pos_quad2.2, pos_quad2.3)];
              let y = x + y.sum_seqalign + twoloop_score + twoloop_score2;
              logsumexp(&mut sum, y);
            }
          }
        }
      }
      let pos_pair2 = if computes_forward_sums {
        (u - T::one(), v - T::one())
      } else {
        (u + T::one(), v + T::one())
      };
      let long_pos_pair2 = (
        pos_pair2.0.to_usize().unwrap(),
        pos_pair2.1.to_usize().unwrap(),
      );
      if let Some(&loopmatch_score) = alignfold_scores.loopmatch_scores.get(&pos_pair) {
        let mut sum2 = NEG_INFINITY;
        if let Some(&x) = sum_mat.get(&pos_pair2) {
          let x = x + align_scores.match2match_score;
          logsumexp(&mut sum2, x);
        }
        if let Some(x) = matchable_poss.get(&pos_pair2.0) {
          for &x in x {
            if computes_forward_sums && x >= pos_pair2.1
              || (!computes_forward_sums && x <= pos_pair2.1)
            {
              continue;
            }
            if let Some(&y) = sum_mat.get(&(pos_pair2.0, x)) {
              let long_x = x.to_usize().unwrap();
              let z = if computes_forward_sums {
                alignfold_scores.range_insert_scores2[long_x + 1][long_pos_pair2.1]
              } else {
                alignfold_scores.range_insert_scores2[long_pos_pair2.1][long_x - 1]
              } + align_scores.match2insert_score;
              let z = y + align_scores.match2insert_score + z;
              logsumexp(&mut sum2, z);
            }
          }
        }
        if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
          for &x in x {
            if computes_forward_sums && x >= pos_pair2.0
              || (!computes_forward_sums && x <= pos_pair2.0)
            {
              continue;
            }
            if let Some(&y) = sum_mat.get(&(x, pos_pair2.1)) {
              let long_x = x.to_usize().unwrap();
              let z = if computes_forward_sums {
                alignfold_scores.range_insert_scores[long_x + 1][long_pos_pair2.0]
              } else {
                alignfold_scores.range_insert_scores[long_pos_pair2.0][long_x - 1]
              } + align_scores.match2insert_score;
              let z = y + align_scores.match2insert_score + z;
              logsumexp(&mut sum2, z);
            }
          }
        }
        if sum2 > NEG_INFINITY {
          sum_mat2.insert(pos_pair2, sum2);
        }
        let term = sum2 + loopmatch_score;
        logsumexp(&mut sum, term);
        if sum > NEG_INFINITY {
          sum_mat.insert(pos_pair, sum);
        }
      }
    }
  }
  (sum_mat, sum_mat2)
}

pub fn has_valid_sums(x: &LoopSums) -> bool {
  x.sum_seqalign > NEG_INFINITY
    || x.sum_multibranch > NEG_INFINITY
    || x.sum_1st_pairmatches > NEG_INFINITY
}

pub fn get_alignfold_probs<T>(inputs: InputsAlignfoldProbsGetter<T>) -> AlignfoldProbMats<T>
where
  T: HashIndex,
{
  let (
    seq_len_pair,
    alignfold_scores,
    max_basepair_span_pair,
    alignfold_sums,
    fold_scores_pair,
    produces_context_profs,
    global_sum,
    pos_quads_hashed_lens,
    produces_match_probs,
    forward_pos_pairs,
    backward_pos_pairs,
    matchable_poss,
    matchable_poss2,
    align_scores,
  ) = inputs;
  let mut alignfold_outside_sums = SumMat4d::<T>::default();
  let mut alignfold_probs = AlignfoldProbMats::<T>::new(&(
    seq_len_pair.0.to_usize().unwrap(),
    seq_len_pair.1.to_usize().unwrap(),
  ));
  let mut prob_coeffs_multibranch = SumMat4d::<T>::default();
  let mut prob_coeffs_multibranch2 = prob_coeffs_multibranch.clone();
  let leftmost_pos_pair = (T::zero(), T::zero());
  let rightmost_pos_pair = (seq_len_pair.0 - T::one(), seq_len_pair.1 - T::one());
  for substr_len in range_inclusive(
    T::from_usize(MIN_SPAN_HAIRPIN_CLOSE).unwrap(),
    max_basepair_span_pair.0,
  )
  .rev()
  {
    for substr_len2 in range_inclusive(
      T::from_usize(MIN_SPAN_HAIRPIN_CLOSE).unwrap(),
      max_basepair_span_pair.1,
    )
    .rev()
    {
      if let Some(pos_pairs) = pos_quads_hashed_lens.get(&(substr_len, substr_len2)) {
        for &(i, k) in pos_pairs {
          let (j, l) = (i + substr_len - T::one(), k + substr_len2 - T::one());
          let pos_quad = (i, j, k, l);
          if let Some(&sum_close) = alignfold_sums.sums_close.get(&pos_quad) {
            let (long_i, long_j, long_k, long_l) = (
              i.to_usize().unwrap(),
              j.to_usize().unwrap(),
              k.to_usize().unwrap(),
              l.to_usize().unwrap(),
            );
            let prob_coeff = sum_close - global_sum;
            let mut sum = NEG_INFINITY;
            let mut forward_term = sum;
            let mut backward_term = sum;
            let pos_pair2 = (i - T::one(), k - T::one());
            if let Some(&x) = alignfold_sums.forward_sums_external2.get(&pos_pair2) {
              logsumexp(&mut forward_term, x);
            }
            let pos_pair2 = (j + T::one(), l + T::one());
            if let Some(&x) = alignfold_sums.backward_sums_external2.get(&pos_pair2) {
              logsumexp(&mut backward_term, x);
            }
            let sum_external = forward_term + backward_term;
            if sum_external > NEG_INFINITY {
              let coeff = alignfold_sums.sums_accessible_external[&pos_quad] - sum_close;
              sum = coeff + sum_external;
            }
            for substr_len3 in range_inclusive(
              substr_len + T::from_usize(2).unwrap(),
              (substr_len + T::from_usize(MAX_2LOOP_LEN + 2).unwrap())
                .min(max_basepair_span_pair.0),
            ) {
              for substr_len4 in range_inclusive(
                substr_len2 + T::from_usize(2).unwrap(),
                (substr_len2 + T::from_usize(MAX_2LOOP_LEN + 2).unwrap())
                  .min(max_basepair_span_pair.1),
              ) {
                if let Some(pos_pairs2) = pos_quads_hashed_lens.get(&(substr_len3, substr_len4)) {
                  for &(m, o) in pos_pairs2 {
                    let (n, p) = (m + substr_len3 - T::one(), o + substr_len4 - T::one());
                    if !(m < i && j < n && o < k && l < p) {
                      continue;
                    }
                    let (long_m, long_n, long_o, long_p) = (
                      m.to_usize().unwrap(),
                      n.to_usize().unwrap(),
                      o.to_usize().unwrap(),
                      p.to_usize().unwrap(),
                    );
                    if long_n - long_j - 1 + long_i - long_m - 1 > MAX_2LOOP_LEN {
                      continue;
                    }
                    if long_p - long_l - 1 + long_k - long_o - 1 > MAX_2LOOP_LEN {
                      continue;
                    }
                    let pos_quad2 = (m, n, o, p);
                    if let Some(&outside_sum) = alignfold_outside_sums.get(&pos_quad2) {
                      let forward_sums = &alignfold_sums.forward_sums_hashed_poss2[&(m, o)];
                      let backward_sums = &alignfold_sums.backward_sums_hashed_poss2[&(n, p)];
                      let mut forward_term = NEG_INFINITY;
                      let mut backward_term = forward_term;
                      let pos_pair2 = (i - T::one(), k - T::one());
                      if let Some(x) = forward_sums.get(&pos_pair2) {
                        logsumexp(&mut forward_term, x.sum_seqalign);
                      }
                      let pos_pair2 = (j + T::one(), l + T::one());
                      if let Some(x) = backward_sums.get(&pos_pair2) {
                        logsumexp(&mut backward_term, x.sum_seqalign);
                      }
                      let sum_2loop = forward_term + backward_term;
                      if sum_2loop > NEG_INFINITY {
                        let pairmatch_score = alignfold_scores.pairmatch_scores[&pos_quad2];
                        let twoloop_score = fold_scores_pair.0.twoloop_scores[&(m, n, i, j)];
                        let twoloop_score2 = fold_scores_pair.1.twoloop_scores[&(o, p, k, l)];
                        let coeff = pairmatch_score + twoloop_score + twoloop_score2 + outside_sum;
                        let sum_2loop = coeff + sum_2loop;
                        logsumexp(&mut sum, sum_2loop);
                        if produces_context_profs {
                          let loop_len_pair = (long_i - long_m - 1, long_n - long_j - 1);
                          let found_bulge_loop = (loop_len_pair.0 == 0) ^ (loop_len_pair.1 == 0);
                          let found_interior_loop = loop_len_pair.0 > 0 && loop_len_pair.1 > 0;
                          let pairmatch_prob_2loop = prob_coeff + sum_2loop;
                          for q in long_m + 1..long_i {
                            if found_bulge_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.0[(q, CONTEXT_INDEX_BULGE)],
                                pairmatch_prob_2loop,
                              );
                            } else if found_interior_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.0
                                  [(q, CONTEXT_INDEX_INTERIOR)],
                                pairmatch_prob_2loop,
                              );
                            }
                          }
                          for q in long_j + 1..long_n {
                            if found_bulge_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.0[(q, CONTEXT_INDEX_BULGE)],
                                pairmatch_prob_2loop,
                              );
                            } else if found_interior_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.0
                                  [(q, CONTEXT_INDEX_INTERIOR)],
                                pairmatch_prob_2loop,
                              );
                            }
                          }
                          let loop_len_pair = (long_k - long_o - 1, long_p - long_l - 1);
                          let found_bulge_loop = (loop_len_pair.0 == 0) ^ (loop_len_pair.1 == 0);
                          let found_interior_loop = loop_len_pair.0 > 0 && loop_len_pair.1 > 0;
                          for q in long_o + 1..long_k {
                            if found_bulge_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.1[(q, CONTEXT_INDEX_BULGE)],
                                pairmatch_prob_2loop,
                              );
                            } else if found_interior_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.1
                                  [(q, CONTEXT_INDEX_INTERIOR)],
                                pairmatch_prob_2loop,
                              );
                            }
                          }
                          for q in long_l + 1..long_p {
                            if found_bulge_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.1[(q, CONTEXT_INDEX_BULGE)],
                                pairmatch_prob_2loop,
                              );
                            } else if found_interior_loop {
                              logsumexp(
                                &mut alignfold_probs.context_profs_pair.1
                                  [(q, CONTEXT_INDEX_INTERIOR)],
                                pairmatch_prob_2loop,
                              );
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            let sum_ratio = alignfold_sums.sums_accessible_multibranch[&pos_quad] - sum_close;
            for (&(u, v), x) in &alignfold_sums.forward_sums_hashed_poss2 {
              if !(u < i && v < k) {
                continue;
              }
              let pos_quad2 = (u, j, v, l);
              let mut forward_term = NEG_INFINITY;
              let mut forward_term2 = forward_term;
              let pos_pair2 = (i - T::one(), k - T::one());
              if let Some(x) = x.get(&pos_pair2) {
                logsumexp(&mut forward_term, x.sum_1ormore_pairmatches);
                logsumexp(&mut forward_term2, x.sum_seqalign);
              }
              let mut sum_multibranch = NEG_INFINITY;
              if let Some(x) = prob_coeffs_multibranch.get(&pos_quad2) {
                let x = x + sum_ratio + forward_term;
                logsumexp(&mut sum_multibranch, x);
              }
              if let Some(x) = prob_coeffs_multibranch2.get(&pos_quad2) {
                let x = x + sum_ratio + forward_term2;
                logsumexp(&mut sum_multibranch, x);
              }
              if sum_multibranch > NEG_INFINITY {
                logsumexp(&mut sum, sum_multibranch);
              }
            }
            if sum > NEG_INFINITY {
              alignfold_outside_sums.insert(pos_quad, sum);
              let pairmatch_prob = prob_coeff + sum;
              if produces_match_probs {
                alignfold_probs
                  .pairmatch_probs
                  .insert(pos_quad, pairmatch_prob);
                match alignfold_probs.match_probs.get_mut(&(i, k)) {
                  Some(x) => {
                    logsumexp(x, pairmatch_prob);
                  }
                  None => {
                    alignfold_probs.match_probs.insert((i, k), pairmatch_prob);
                  }
                }
                match alignfold_probs.match_probs.get_mut(&(j, l)) {
                  Some(x) => {
                    logsumexp(x, pairmatch_prob);
                  }
                  None => {
                    alignfold_probs.match_probs.insert((j, l), pairmatch_prob);
                  }
                }
              }
              match alignfold_probs.basepair_probs_pair.0.get_mut(&(i, j)) {
                Some(x) => {
                  logsumexp(x, pairmatch_prob);
                }
                None => {
                  alignfold_probs
                    .basepair_probs_pair
                    .0
                    .insert((i, j), pairmatch_prob);
                }
              }
              match alignfold_probs.basepair_probs_pair.1.get_mut(&(k, l)) {
                Some(x) => {
                  logsumexp(x, pairmatch_prob);
                }
                None => {
                  alignfold_probs
                    .basepair_probs_pair
                    .1
                    .insert((k, l), pairmatch_prob);
                }
              }
              if produces_context_profs {
                logsumexp(
                  &mut alignfold_probs.context_profs_pair.0[(long_i, CONTEXT_INDEX_BASEPAIR)],
                  pairmatch_prob,
                );
                logsumexp(
                  &mut alignfold_probs.context_profs_pair.0[(long_j, CONTEXT_INDEX_BASEPAIR)],
                  pairmatch_prob,
                );
                logsumexp(
                  &mut alignfold_probs.context_profs_pair.1[(long_k, CONTEXT_INDEX_BASEPAIR)],
                  pairmatch_prob,
                );
                logsumexp(
                  &mut alignfold_probs.context_profs_pair.1[(long_l, CONTEXT_INDEX_BASEPAIR)],
                  pairmatch_prob,
                );
              }
              let pairmatch_score = alignfold_scores.pairmatch_scores[&pos_quad];
              let multibranch_close_score = fold_scores_pair.0.multibranch_close_scores[&(i, j)];
              let multibranch_close_score2 = fold_scores_pair.1.multibranch_close_scores[&(k, l)];
              let coeff =
                sum + pairmatch_score + multibranch_close_score + multibranch_close_score2;
              let backward_sums = &alignfold_sums.backward_sums_hashed_poss2[&(j, l)];
              for pos_pair in alignfold_scores.loopmatch_scores.keys() {
                let &(u, v) = pos_pair;
                if !(i < u && u < j && k < v && v < l) {
                  continue;
                }
                let mut backward_term = NEG_INFINITY;
                let mut backward_term2 = backward_term;
                let pos_pair2 = (u + T::one(), v + T::one());
                if let Some(x) = backward_sums.get(&pos_pair2) {
                  logsumexp(&mut backward_term, x.sum_0ormore_pairmatches);
                  logsumexp(&mut backward_term2, x.sum_1ormore_pairmatches);
                }
                let pos_quad2 = (i, u, k, v);
                let x = coeff + backward_term;
                match prob_coeffs_multibranch.get_mut(&pos_quad2) {
                  Some(y) => {
                    logsumexp(y, x);
                  }
                  None => {
                    prob_coeffs_multibranch.insert(pos_quad2, x);
                  }
                }
                let x = coeff + backward_term2;
                match prob_coeffs_multibranch2.get_mut(&pos_quad2) {
                  Some(y) => {
                    logsumexp(y, x);
                  }
                  None => {
                    prob_coeffs_multibranch2.insert(pos_quad2, x);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  for x in alignfold_probs.basepair_probs_pair.0.values_mut() {
    *x = expf(*x);
  }
  for x in alignfold_probs.basepair_probs_pair.1.values_mut() {
    *x = expf(*x);
  }
  if produces_context_profs || produces_match_probs {
    let mut unpair_probs_range_external =
      (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    let mut unpair_probs_range_hairpin =
      (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
    for u in range(T::zero(), seq_len_pair.0 - T::one()) {
      let long_u = u.to_usize().unwrap();
      for v in range(T::zero(), seq_len_pair.1 - T::one()) {
        let pos_pair = (u, v);
        let long_v = v.to_usize().unwrap();
        let pos_pair2 = (u + T::one(), v + T::one());
        let long_pos_pair2 = (
          pos_pair2.0.to_usize().unwrap(),
          pos_pair2.1.to_usize().unwrap(),
        );
        if u > T::zero() && v > T::zero() {
          let pos_pair_loopmatch = (u - T::one(), v - T::one());
          if let Some(&loopmatch_score) = alignfold_scores.loopmatch_scores.get(&pos_pair) {
            let mut forward_term = NEG_INFINITY;
            if let Some(&x) = alignfold_sums
              .forward_sums_external2
              .get(&pos_pair_loopmatch)
            {
              logsumexp(&mut forward_term, x);
            }
            let mut backward_term = NEG_INFINITY;
            if let Some(&x) = alignfold_sums.backward_sums_external2.get(&pos_pair2) {
              logsumexp(&mut backward_term, x);
            }
            let loopmatch_prob_external =
              loopmatch_score + forward_term + backward_term - global_sum;
            if produces_context_profs {
              logsumexp(
                &mut alignfold_probs.context_profs_pair.0[(long_u, CONTEXT_INDEX_EXTERNAL)],
                loopmatch_prob_external,
              );
              logsumexp(
                &mut alignfold_probs.context_profs_pair.1[(long_v, CONTEXT_INDEX_EXTERNAL)],
                loopmatch_prob_external,
              );
            }
            if produces_match_probs {
              match alignfold_probs.loopmatch_probs.get_mut(&pos_pair) {
                Some(x) => {
                  logsumexp(x, loopmatch_prob_external);
                }
                None => {
                  alignfold_probs
                    .loopmatch_probs
                    .insert(pos_pair, loopmatch_prob_external);
                }
              }
            }
          }
        }
        if produces_context_profs {
          if let Some(&sum) = alignfold_sums.forward_sums_external.get(&pos_pair) {
            let begins_sum = pos_pair == leftmost_pos_pair;
            let forward_term = sum
              + if begins_sum {
                align_scores.init_insert_score
              } else {
                align_scores.match2insert_score
              };
            if let Some(x) = matchable_poss.get(&pos_pair2.0) {
              for &x in x {
                if x <= pos_pair2.1 {
                  continue;
                }
                let pos_pair3 = (pos_pair2.0, x);
                if let Some(&y) = alignfold_sums.backward_sums_external.get(&pos_pair3) {
                  let long_x = x.to_usize().unwrap();
                  let ends_sum = pos_pair3 == rightmost_pos_pair;
                  let z = alignfold_scores.range_insert_scores2[long_pos_pair2.1][long_x - 1]
                    + if ends_sum {
                      0.
                    } else {
                      align_scores.match2insert_score
                    };
                  let z = forward_term + y + z - global_sum;
                  let pos_pair4 = (pos_pair2.1, x - T::one());
                  match unpair_probs_range_external.1.get_mut(&pos_pair4) {
                    Some(x) => {
                      logsumexp(x, z);
                    }
                    None => {
                      unpair_probs_range_external.1.insert(pos_pair4, z);
                    }
                  }
                }
              }
            }
            if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
              for &x in x {
                if x <= pos_pair2.0 {
                  continue;
                }
                let pos_pair3 = (x, pos_pair2.1);
                if let Some(&y) = alignfold_sums.backward_sums_external.get(&pos_pair3) {
                  let long_x = x.to_usize().unwrap();
                  let ends_sum = pos_pair3 == rightmost_pos_pair;
                  let z = alignfold_scores.range_insert_scores[long_pos_pair2.0][long_x - 1]
                    + if ends_sum {
                      0.
                    } else {
                      align_scores.match2insert_score
                    };
                  let z = forward_term + y + z - global_sum;
                  let pos_pair4 = (pos_pair2.0, x - T::one());
                  match unpair_probs_range_external.0.get_mut(&pos_pair4) {
                    Some(x) => {
                      logsumexp(x, z);
                    }
                    None => {
                      unpair_probs_range_external.0.insert(pos_pair4, z);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    for (pos_quad, &outside_sum) in &alignfold_outside_sums {
      let (i, j, k, l) = *pos_quad;
      let pairmatch_score = alignfold_scores.pairmatch_scores[pos_quad];
      let prob_coeff = outside_sum - global_sum + pairmatch_score;
      let hairpin_score = match fold_scores_pair.0.hairpin_scores.get(&(i, j)) {
        Some(&x) => x,
        None => NEG_INFINITY,
      };
      let hairpin_score2 = match fold_scores_pair.1.hairpin_scores.get(&(k, l)) {
        Some(&x) => x,
        None => NEG_INFINITY,
      };
      let prob_coeff_hairpin = prob_coeff + hairpin_score + hairpin_score2;
      let multibranch_close_score = fold_scores_pair.0.multibranch_close_scores[&(i, j)];
      let multibranch_close_score2 = fold_scores_pair.1.multibranch_close_scores[&(k, l)];
      let prob_coeff_multibranch = prob_coeff + multibranch_close_score + multibranch_close_score2;
      let forward_sums = &alignfold_sums.forward_sums_hashed_poss[&(i, k)];
      let backward_sums = &alignfold_sums.backward_sums_hashed_poss[&(j, l)];
      let forward_sums2 = &alignfold_sums.forward_sums_hashed_poss2[&(i, k)];
      let backward_sums2 = &alignfold_sums.backward_sums_hashed_poss2[&(j, l)];
      let (_, forward_sums_2loop2) = if produces_match_probs {
        let computes_forward_sums = true;
        get_2loop_sums((
          pos_quad,
          alignfold_sums,
          computes_forward_sums,
          forward_pos_pairs,
          fold_scores_pair,
          alignfold_scores,
          matchable_poss,
          matchable_poss2,
          align_scores,
        ))
      } else {
        (SparseSumMat::<T>::default(), SparseSumMat::<T>::default())
      };
      let (_, backward_sums_2loop2) = if produces_match_probs {
        let computes_forward_sums = false;
        get_2loop_sums((
          pos_quad,
          alignfold_sums,
          computes_forward_sums,
          backward_pos_pairs,
          fold_scores_pair,
          alignfold_scores,
          matchable_poss,
          matchable_poss2,
          align_scores,
        ))
      } else {
        (SparseSumMat::<T>::default(), SparseSumMat::<T>::default())
      };
      for u in range(i, j) {
        let long_u = u.to_usize().unwrap();
        for v in range(k, l) {
          let pos_pair = (u, v);
          let long_v = v.to_usize().unwrap();
          let pos_pair2 = (u + T::one(), v + T::one());
          let long_pos_pair2 = (
            pos_pair2.0.to_usize().unwrap(),
            pos_pair2.1.to_usize().unwrap(),
          );
          if let Some(&loopmatch_score) = alignfold_scores.loopmatch_scores.get(&pos_pair) {
            let mut backward_term_match_seqalign = NEG_INFINITY;
            let mut backward_term_match_multibranch = backward_term_match_seqalign;
            let mut backward_term_match_1ormore = backward_term_match_seqalign;
            let mut backward_term_match_0ormore = backward_term_match_seqalign;
            let mut backward_term_match_2loop = backward_term_match_seqalign;
            if let Some(x) = backward_sums2.get(&pos_pair2) {
              logsumexp(&mut backward_term_match_seqalign, x.sum_seqalign);
              if produces_match_probs {
                logsumexp(&mut backward_term_match_multibranch, x.sum_multibranch);
                logsumexp(&mut backward_term_match_1ormore, x.sum_1ormore_pairmatches);
                logsumexp(&mut backward_term_match_0ormore, x.sum_0ormore_pairmatches);
              }
            }
            if produces_match_probs {
              if let Some(&x) = backward_sums_2loop2.get(&pos_pair2) {
                logsumexp(&mut backward_term_match_2loop, x);
              }
            }
            let pos_pair_loopmatch = (u - T::one(), v - T::one());
            let mut loopmatch_prob_hairpin = NEG_INFINITY;
            let mut loopmatch_prob_multibranch = loopmatch_prob_hairpin;
            let mut loopmatch_prob_2loop = loopmatch_prob_hairpin;
            if let Some(x) = forward_sums2.get(&pos_pair_loopmatch) {
              let y = prob_coeff_hairpin
                + loopmatch_score
                + x.sum_seqalign
                + backward_term_match_seqalign;
              logsumexp(&mut loopmatch_prob_hairpin, y);
              if produces_match_probs {
                let y = prob_coeff_multibranch
                  + loopmatch_score
                  + x.sum_seqalign
                  + backward_term_match_multibranch;
                logsumexp(&mut loopmatch_prob_multibranch, y);
                let y = prob_coeff_multibranch
                  + loopmatch_score
                  + x.sum_1st_pairmatches
                  + backward_term_match_1ormore;
                logsumexp(&mut loopmatch_prob_multibranch, y);
                let y = prob_coeff_multibranch
                  + loopmatch_score
                  + x.sum_multibranch
                  + backward_term_match_0ormore;
                logsumexp(&mut loopmatch_prob_multibranch, y);
                let y = prob_coeff + loopmatch_score + x.sum_seqalign + backward_term_match_2loop;
                logsumexp(&mut loopmatch_prob_2loop, y);
              }
            }
            if produces_context_profs {
              logsumexp(
                &mut alignfold_probs.context_profs_pair.0[(long_u, CONTEXT_INDEX_HAIRPIN)],
                loopmatch_prob_hairpin,
              );
              logsumexp(
                &mut alignfold_probs.context_profs_pair.1[(long_v, CONTEXT_INDEX_HAIRPIN)],
                loopmatch_prob_hairpin,
              );
            }
            if produces_match_probs {
              if let Some(&x) = forward_sums_2loop2.get(&pos_pair_loopmatch) {
                let x = prob_coeff + loopmatch_score + x + backward_term_match_seqalign;
                logsumexp(&mut loopmatch_prob_2loop, x);
              }
              let mut prob = NEG_INFINITY;
              logsumexp(&mut prob, loopmatch_prob_hairpin);
              logsumexp(&mut prob, loopmatch_prob_multibranch);
              logsumexp(&mut prob, loopmatch_prob_2loop);
              match alignfold_probs.loopmatch_probs.get_mut(&pos_pair) {
                Some(x) => {
                  logsumexp(x, prob);
                }
                None => {
                  alignfold_probs.loopmatch_probs.insert(pos_pair, prob);
                }
              }
            }
          }
          if produces_context_profs {
            if let Some(sums) = forward_sums.get(&pos_pair) {
              let forward_term = sums.sum_seqalign + align_scores.match2insert_score;
              if let Some(x) = matchable_poss.get(&pos_pair2.0) {
                for &x in x {
                  if x <= pos_pair2.1 {
                    continue;
                  }
                  let pos_pair3 = (pos_pair2.0, x);
                  if let Some(y) = backward_sums.get(&pos_pair3) {
                    let long_x = x.to_usize().unwrap();
                    let z = alignfold_scores.range_insert_scores2[long_pos_pair2.1][long_x - 1]
                      + align_scores.match2insert_score;
                    let z = prob_coeff_hairpin + forward_term + y.sum_seqalign + z;
                    let pos_pair4 = (pos_pair2.1, x - T::one());
                    match unpair_probs_range_hairpin.1.get_mut(&pos_pair4) {
                      Some(x) => {
                        logsumexp(x, z);
                      }
                      None => {
                        unpair_probs_range_hairpin.1.insert(pos_pair4, z);
                      }
                    }
                  }
                }
              }
              if let Some(x) = matchable_poss2.get(&pos_pair2.1) {
                for &x in x {
                  if x <= pos_pair2.0 {
                    continue;
                  }
                  let pos_pair3 = (x, pos_pair2.1);
                  if let Some(y) = backward_sums.get(&pos_pair3) {
                    let long_x = x.to_usize().unwrap();
                    let z = alignfold_scores.range_insert_scores[long_pos_pair2.0][long_x - 1]
                      + align_scores.match2insert_score;
                    let z = prob_coeff_hairpin + forward_term + z + y.sum_seqalign;
                    let pos_pair4 = (pos_pair2.0, x - T::one());
                    match unpair_probs_range_hairpin.0.get_mut(&pos_pair4) {
                      Some(x) => {
                        logsumexp(x, z);
                      }
                      None => {
                        unpair_probs_range_hairpin.0.insert(pos_pair4, z);
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
    if produces_context_profs {
      for (pos_pair, &x) in &unpair_probs_range_external.0 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(
            &mut alignfold_probs.context_profs_pair.0[(long_i, CONTEXT_INDEX_EXTERNAL)],
            x,
          );
        }
      }
      for (pos_pair, &x) in &unpair_probs_range_external.1 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(
            &mut alignfold_probs.context_profs_pair.1[(long_i, CONTEXT_INDEX_EXTERNAL)],
            x,
          );
        }
      }
      for (pos_pair, &x) in &unpair_probs_range_hairpin.0 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(
            &mut alignfold_probs.context_profs_pair.0[(long_i, CONTEXT_INDEX_HAIRPIN)],
            x,
          );
        }
      }
      for (pos_pair, &x) in &unpair_probs_range_hairpin.1 {
        for i in range_inclusive(pos_pair.0, pos_pair.1) {
          let long_i = i.to_usize().unwrap();
          logsumexp(
            &mut alignfold_probs.context_profs_pair.1[(long_i, CONTEXT_INDEX_HAIRPIN)],
            x,
          );
        }
      }
      alignfold_probs
        .context_profs_pair
        .0
        .slice_mut(s![.., ..CONTEXT_INDEX_MULTIBRANCH])
        .mapv_inplace(expf);
      let fold = 1.
        - alignfold_probs
          .context_profs_pair
          .0
          .slice_mut(s![.., ..CONTEXT_INDEX_MULTIBRANCH])
          .sum_axis(Axis(1));
      alignfold_probs
        .context_profs_pair
        .0
        .slice_mut(s![.., CONTEXT_INDEX_MULTIBRANCH])
        .assign(&fold);
      alignfold_probs
        .context_profs_pair
        .1
        .slice_mut(s![.., ..CONTEXT_INDEX_MULTIBRANCH])
        .mapv_inplace(expf);
      let fold = 1.
        - alignfold_probs
          .context_profs_pair
          .1
          .slice_mut(s![.., ..CONTEXT_INDEX_MULTIBRANCH])
          .sum_axis(Axis(1));
      alignfold_probs
        .context_profs_pair
        .1
        .slice_mut(s![.., CONTEXT_INDEX_MULTIBRANCH])
        .assign(&fold);
    }
    if produces_match_probs {
      for (pos_pair, x) in alignfold_probs.loopmatch_probs.iter_mut() {
        match alignfold_probs.match_probs.get_mut(pos_pair) {
          Some(y) => {
            logsumexp(y, *x);
            *y = expf(*y);
          }
          None => {
            alignfold_probs.match_probs.insert(*pos_pair, expf(*x));
          }
        }
        *x = expf(*x);
      }
      for x in alignfold_probs.pairmatch_probs.values_mut() {
        *x = expf(*x);
      }
    }
  }
  alignfold_probs
}

pub fn pair_probs2avg_probs<T>(
  alignfold_probs_hashed_ids: &AlignfoldProbsHashedIds<T>,
  rna_id: RnaId,
  num_rnas: usize,
  unpair_probs_len: usize,
  produces_context_profs: bool,
) -> AlignfoldProbMatsAvg<T>
where
  T: HashIndex,
{
  let weight = 1. / (num_rnas - 1) as Prob;
  let mut alignfold_probs_avg = AlignfoldProbMatsAvg::new(unpair_probs_len);
  for rna_id2 in 0..num_rnas {
    if rna_id == rna_id2 {
      continue;
    }
    let rna_id_pair = if rna_id < rna_id2 {
      (rna_id, rna_id2)
    } else {
      (rna_id2, rna_id)
    };
    let alignfold_probs = &alignfold_probs_hashed_ids[&rna_id_pair];
    let basepair_probs = if rna_id < rna_id2 {
      &alignfold_probs.basepair_probs_pair.0
    } else {
      &alignfold_probs.basepair_probs_pair.1
    };
    for (x, &y) in basepair_probs.iter() {
      let y = weight * y;
      match alignfold_probs_avg.basepair_probs.get_mut(x) {
        Some(x) => {
          *x += y;
        }
        None => {
          alignfold_probs_avg.basepair_probs.insert(*x, y);
        }
      }
    }
    if produces_context_profs {
      let context_profs = if rna_id < rna_id2 {
        &alignfold_probs.context_profs_pair.0
      } else {
        &alignfold_probs.context_profs_pair.1
      };
      alignfold_probs_avg.context_profs =
        alignfold_probs_avg.context_profs + weight * context_profs;
    }
  }
  alignfold_probs_avg
}

pub fn get_max_basepair_span<T>(basepair_probs: &SparseProbMat<T>) -> T
where
  T: HashIndex,
{
  let max_basepair_span = basepair_probs.keys().map(|x| x.1 - x.0 + T::one()).max();
  match max_basepair_span {
    Some(x) => x,
    None => T::zero(),
  }
}

pub fn filter_basepair_probs<T>(
  basepair_probs: &SparseProbMat<T>,
  min_basepair_prob: Prob,
) -> SparseProbMat<T>
where
  T: HashIndex,
{
  basepair_probs
    .iter()
    .filter(|(_, &x)| x >= min_basepair_prob)
    .map(|(x, &y)| ((x.0 + T::one(), x.1 + T::one()), y))
    .collect()
}

pub fn filter_match_probs<T>(match_probs: &ProbMat, min_match_prob: Prob) -> SparseProbMat<T>
where
  T: HashIndex,
{
  let mut match_probs_sparse = SparseProbMat::<T>::default();
  for (i, x) in match_probs.iter().enumerate() {
    let i = T::from_usize(i).unwrap();
    for (j, &x) in x.iter().enumerate() {
      if x >= min_match_prob {
        let j = T::from_usize(j).unwrap();
        match_probs_sparse.insert((i, j), x);
      }
    }
  }
  match_probs_sparse
}

pub fn consprob<'a, T>(
  thread_pool: &mut Pool,
  seqs: &SeqSlices<'a>,
  min_basepair_prob: Prob,
  min_match_prob: Prob,
  produces_context_profs: bool,
  produces_match_probs: bool,
  align_scores: &AlignScores,
) -> (ProbMatSetsAvg<T>, MatchProbsHashedIds<T>)
where
  T: HashIndex,
{
  let num_seqs = seqs.len();
  let mut basepair_prob_mats = vec![SparseProbMat::<T>::new(); num_seqs];
  let mut max_basepair_spans = vec![T::zero(); num_seqs];
  let mut fold_score_sets = vec![FoldScores::<T>::new(); num_seqs];
  let uses_contra_model = false;
  let allows_short_hairpins = false;
  thread_pool.scoped(|scope| {
    for (x, y, z, a) in multizip((
      basepair_prob_mats.iter_mut(),
      max_basepair_spans.iter_mut(),
      seqs.iter(),
      fold_score_sets.iter_mut(),
    )) {
      let b = z.len();
      scope.execute(move || {
        let (c, d) = mccaskill_algo(
          &z[1..b - 1],
          uses_contra_model,
          allows_short_hairpins,
          &FoldScoreSets::new(0.),
        );
        *x = filter_basepair_probs::<T>(&c, min_basepair_prob);
        *a = filter_fold_scores(&d, x);
        *y = get_max_basepair_span::<T>(x);
      });
    }
  });
  let mut alignfold_probs_hashed_ids = AlignfoldProbsHashedIds::<T>::default();
  let mut match_probs_hashed_ids = SparseProbsHashedIds::<T>::default();
  for x in 0..num_seqs {
    for y in x + 1..num_seqs {
      let y = (x, y);
      alignfold_probs_hashed_ids.insert(y, AlignfoldProbMats::<T>::origin());
      match_probs_hashed_ids.insert(y, SparseProbMat::<T>::default());
    }
  }
  thread_pool.scoped(|x| {
    for (y, z) in match_probs_hashed_ids.iter_mut() {
      let y = (seqs[y.0], seqs[y.1]);
      x.execute(move || {
        *z = filter_match_probs(&durbin_algo(&y, align_scores), min_match_prob);
      });
    }
  });
  thread_pool.scoped(|x| {
    for (y, z) in alignfold_probs_hashed_ids.iter_mut() {
      let seq_pair = (seqs[y.0], seqs[y.1]);
      let max_basepair_span_pair = (max_basepair_spans[y.0], max_basepair_spans[y.1]);
      let basepair_probs_pair = (&basepair_prob_mats[y.0], &basepair_prob_mats[y.1]);
      let fold_scores_pair = (&fold_score_sets[y.0], &fold_score_sets[y.1]);
      let seq_len_pair = (
        T::from_usize(seq_pair.0.len()).unwrap(),
        T::from_usize(seq_pair.1.len()).unwrap(),
      );
      let match_probs = &match_probs_hashed_ids[y];
      let (
        forward_pos_pairs,
        backward_pos_pairs,
        pos_quad_mat,
        pos_quads_hashed_lens,
        matchable_poss,
        matchable_poss2,
      ) = get_sparse_poss(&basepair_probs_pair, match_probs, &seq_len_pair);
      x.execute(move || {
        let alignfold_scores =
          AlignfoldScores::<T>::new(&seq_pair, &pos_quad_mat, match_probs, align_scores);
        *z = consprob_core::<T>((
          &seq_len_pair,
          &alignfold_scores,
          &max_basepair_span_pair,
          &fold_scores_pair,
          produces_context_profs,
          &forward_pos_pairs,
          &backward_pos_pairs,
          &pos_quads_hashed_lens,
          produces_match_probs,
          &matchable_poss,
          &matchable_poss2,
          align_scores,
        ));
      });
    }
  });
  let mut alignfold_prob_mats_avg = vec![AlignfoldProbMatsAvg::<T>::origin(); num_seqs];
  thread_pool.scoped(|x| {
    let y = &alignfold_probs_hashed_ids;
    for (z, a) in alignfold_prob_mats_avg.iter_mut().enumerate() {
      let b = seqs[z].len();
      x.execute(move || {
        *a = pair_probs2avg_probs::<T>(y, z, num_seqs, b, produces_context_profs);
      });
    }
  });
  let mut match_probs_hashed_ids = MatchProbsHashedIds::<T>::default();
  if produces_match_probs {
    for x in 0..num_seqs {
      for y in x + 1..num_seqs {
        let y = (x, y);
        let z = &alignfold_probs_hashed_ids[&y];
        let mut a = MatchProbMats::<T>::new();
        a.loopmatch_probs = z.loopmatch_probs.clone();
        a.pairmatch_probs = z.pairmatch_probs.clone();
        a.match_probs = z.match_probs.clone();
        match_probs_hashed_ids.insert(y, a);
      }
    }
  }
  (alignfold_prob_mats_avg, match_probs_hashed_ids)
}

pub fn filter_fold_scores<T>(
  fold_scores_dense: &FoldScores<T>,
  basepair_prob_mat: &SparseProbMat<T>,
) -> FoldScores<T>
where
  T: HashIndex,
{
  let mut fold_scores = FoldScores::new();
  fold_scores.hairpin_scores = fold_scores_dense
    .hairpin_scores
    .iter()
    .map(|(x, &y)| ((x.0 + T::one(), x.1 + T::one()), y))
    .filter(|(x, _)| basepair_prob_mat.contains_key(x))
    .collect();
  fold_scores.twoloop_scores = fold_scores_dense
    .twoloop_scores
    .iter()
    .map(|(x, &y)| {
      (
        (
          x.0 + T::one(),
          x.1 + T::one(),
          x.2 + T::one(),
          x.3 + T::one(),
        ),
        y,
      )
    })
    .filter(|(x, _)| {
      basepair_prob_mat.contains_key(&(x.0, x.1)) && basepair_prob_mat.contains_key(&(x.2, x.3))
    })
    .collect();
  fold_scores.multibranch_close_scores = fold_scores_dense
    .multibranch_close_scores
    .iter()
    .map(|(x, &y)| ((x.0 + T::one(), x.1 + T::one()), y))
    .filter(|(x, _)| basepair_prob_mat.contains_key(x))
    .collect();
  fold_scores.accessible_scores = fold_scores_dense
    .accessible_scores
    .iter()
    .map(|(x, &y)| ((x.0 + T::one(), x.1 + T::one()), y))
    .filter(|(x, _)| basepair_prob_mat.contains_key(x))
    .collect();
  fold_scores
}

pub fn write_alignfold_prob_mats<T>(
  output_dir_path: &Path,
  alignfold_prob_mats_avg: &ProbMatSetsAvg<T>,
  match_probs_hashed_ids: &MatchProbsHashedIds<T>,
  produces_context_profs: bool,
  produces_match_probs: bool,
) where
  T: HashIndex,
{
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let basepair_probs_file_path = output_dir_path.join(BASEPAIR_PROBS_FILE);
  let mut writer_basepair_prob = BufWriter::new(File::create(basepair_probs_file_path).unwrap());
  let mut buf_basepair_prob = String::new();
  for (rna_id, alignfold_probs_avg) in alignfold_prob_mats_avg.iter().enumerate() {
    let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
    for (x, &y) in alignfold_probs_avg.basepair_probs.iter() {
      buf_rna_id.push_str(&format!("{},{},{} ", x.0 - T::one(), x.1 - T::one(), y));
    }
    buf_basepair_prob.push_str(&buf_rna_id);
  }
  let _ = writer_basepair_prob.write_all(buf_basepair_prob.as_bytes());
  if produces_context_profs {
    let basepair_prob_file_path2 = output_dir_path.join(BASEPAIR_PROBS_FILE2);
    let mut writer_basepair_prob2 = BufWriter::new(File::create(basepair_prob_file_path2).unwrap());
    let unpair_prob_file_path = output_dir_path.join(UNPAIR_PROBS_FILE_HAIRPIN);
    let mut writer_unpair_prob_hairpin =
      BufWriter::new(File::create(unpair_prob_file_path).unwrap());
    let unpair_prob_file_path = output_dir_path.join(UNPAIR_PROBS_FILE_BULGE);
    let mut writer_unpair_prob_bulge = BufWriter::new(File::create(unpair_prob_file_path).unwrap());
    let unpair_prob_file_path = output_dir_path.join(UNPAIR_PROBS_FILE_INTERIOR);
    let mut writer_unpair_prob_interior =
      BufWriter::new(File::create(unpair_prob_file_path).unwrap());
    let unpair_prob_file_path = output_dir_path.join(UNPAIR_PROBS_FILE_MULTIBRANCH);
    let mut writer_unpair_prob_multibranch =
      BufWriter::new(File::create(unpair_prob_file_path).unwrap());
    let unpair_prob_file_path = output_dir_path.join(UNPAIR_PROBS_FILE_EXTERNAL);
    let mut writer_unpair_prob_external =
      BufWriter::new(File::create(unpair_prob_file_path).unwrap());
    let mut buf_basepair_prob2 = String::new();
    let mut buf_unpair_prob_hairpin = String::new();
    let mut buf_unpair_prob_bulge = buf_unpair_prob_hairpin.clone();
    let mut buf_unpair_prob_interior = buf_unpair_prob_hairpin.clone();
    let mut buf_unpair_prob_external = buf_unpair_prob_hairpin.clone();
    let mut buf_unpair_prob_multibranch = buf_unpair_prob_hairpin.clone();
    for (rna_id, alignfold_probs_avg) in alignfold_prob_mats_avg.iter().enumerate() {
      let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
      let slice = alignfold_probs_avg
        .context_profs
        .slice(s![.., CONTEXT_INDEX_BASEPAIR]);
      let slice_len = slice.len();
      for (i, &x) in slice.iter().enumerate() {
        if i == 0 || i == slice_len - 1 {
          continue;
        }
        buf_rna_id.push_str(&format!("{},{} ", i - 1, x));
      }
      buf_basepair_prob2.push_str(&buf_rna_id);
      let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
      let slice = alignfold_probs_avg
        .context_profs
        .slice(s![.., CONTEXT_INDEX_HAIRPIN]);
      let slice_len = slice.len();
      for (i, &x) in slice.iter().enumerate() {
        if i == 0 || i == slice_len - 1 {
          continue;
        }
        buf_rna_id.push_str(&format!("{},{} ", i - 1, x));
      }
      buf_unpair_prob_hairpin.push_str(&buf_rna_id);
      let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
      let slice = alignfold_probs_avg
        .context_profs
        .slice(s![.., CONTEXT_INDEX_BULGE]);
      let slice_len = slice.len();
      for (i, &x) in slice.iter().enumerate() {
        if i == 0 || i == slice_len - 1 {
          continue;
        }
        buf_rna_id.push_str(&format!("{},{} ", i - 1, x));
      }
      buf_unpair_prob_bulge.push_str(&buf_rna_id);
      let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
      let slice = alignfold_probs_avg
        .context_profs
        .slice(s![.., CONTEXT_INDEX_INTERIOR]);
      let slice_len = slice.len();
      for (i, &x) in slice.iter().enumerate() {
        if i == 0 || i == slice_len - 1 {
          continue;
        }
        buf_rna_id.push_str(&format!("{},{} ", i - 1, x));
      }
      buf_unpair_prob_interior.push_str(&buf_rna_id);
      let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
      let slice = alignfold_probs_avg
        .context_profs
        .slice(s![.., CONTEXT_INDEX_MULTIBRANCH]);
      let slice_len = slice.len();
      for (i, &x) in slice.iter().enumerate() {
        if i == 0 || i == slice_len - 1 {
          continue;
        }
        buf_rna_id.push_str(&format!("{},{} ", i - 1, x));
      }
      buf_unpair_prob_multibranch.push_str(&buf_rna_id);
      let mut buf_rna_id = format!("\n\n>{}\n", rna_id);
      let slice = alignfold_probs_avg
        .context_profs
        .slice(s![.., CONTEXT_INDEX_EXTERNAL]);
      let slice_len = slice.len();
      for (i, &x) in slice.iter().enumerate() {
        if i == 0 || i == slice_len - 1 {
          continue;
        }
        buf_rna_id.push_str(&format!("{},{} ", i - 1, x));
      }
      buf_unpair_prob_external.push_str(&buf_rna_id);
    }
    let _ = writer_basepair_prob2.write_all(buf_basepair_prob2.as_bytes());
    let _ = writer_unpair_prob_hairpin.write_all(buf_unpair_prob_hairpin.as_bytes());
    let _ = writer_unpair_prob_bulge.write_all(buf_unpair_prob_bulge.as_bytes());
    let _ = writer_unpair_prob_interior.write_all(buf_unpair_prob_interior.as_bytes());
    let _ = writer_unpair_prob_multibranch.write_all(buf_unpair_prob_multibranch.as_bytes());
    let _ = writer_unpair_prob_external.write_all(buf_unpair_prob_external.as_bytes());
  }
  if produces_match_probs {
    let loopmatch_probs_file = output_dir_path.join(LOOPMATCH_PROBS_FILE);
    let pairmatch_probs_file = output_dir_path.join(PAIRMATCH_PROBS_FILE);
    let match_probs_file = output_dir_path.join(MATCH_PROBS_FILE);
    let mut writer_loopmatch_prob = BufWriter::new(File::create(loopmatch_probs_file).unwrap());
    let mut writer_pairmatch_prob = BufWriter::new(File::create(pairmatch_probs_file).unwrap());
    let mut writer_match_prob = BufWriter::new(File::create(match_probs_file).unwrap());
    let mut buf_loopmatch_prob = String::new();
    let mut buf_pairmatch_prob = String::new();
    let mut buf_match_prob = String::new();
    for (rna_id_pair, match_probs) in match_probs_hashed_ids.iter() {
      let mut buf_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
      for (x, &y) in match_probs.loopmatch_probs.iter() {
        buf_rna_id_pair.push_str(&format!("{},{},{} ", x.0 - T::one(), x.1 - T::one(), y));
      }
      buf_loopmatch_prob.push_str(&buf_rna_id_pair);
      let mut buf_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
      for (x, &y) in match_probs.pairmatch_probs.iter() {
        buf_rna_id_pair.push_str(&format!(
          "{},{},{},{},{} ",
          x.0 - T::one(),
          x.1 - T::one(),
          x.2 - T::one(),
          x.3 - T::one(),
          y
        ));
      }
      buf_pairmatch_prob.push_str(&buf_rna_id_pair);
      let mut buf_rna_id_pair = format!("\n\n>{},{}\n", rna_id_pair.0, rna_id_pair.1);
      for (x, &y) in match_probs.match_probs.iter() {
        buf_rna_id_pair.push_str(&format!("{},{},{} ", x.0 - T::one(), x.1 - T::one(), y));
      }
      buf_match_prob.push_str(&buf_rna_id_pair);
    }
    let _ = writer_loopmatch_prob.write_all(buf_loopmatch_prob.as_bytes());
    let _ = writer_pairmatch_prob.write_all(buf_pairmatch_prob.as_bytes());
    let _ = writer_match_prob.write_all(buf_match_prob.as_bytes());
  }
}

pub fn get_sparse_poss<T>(
  basepair_probs_pair: &ProbMatPair<T>,
  match_probs: &SparseProbMat<T>,
  seq_len_pair: &(T, T),
) -> OutputsSparsePossGetter<T>
where
  T: HashIndex,
{
  let mut forward_pos_pairs = PosPairMatSet::<T>::default();
  let mut backward_pos_pairs = PosPairMatSet::<T>::default();
  let mut pos_quads = PosQuadMat::<T>::default();
  let mut pos_quads_hashed_lens = PosQuadsHashedLens::<T>::default();
  let mut matchable_poss = SparsePosSets::<T>::default();
  let mut matchable_poss2 = SparsePosSets::<T>::default();
  for pos_pair in basepair_probs_pair.0.keys() {
    for pos_pair2 in basepair_probs_pair.1.keys() {
      let forward_pos_pair = (pos_pair.0, pos_pair2.0);
      let backward_pos_pair = (pos_pair.1, pos_pair2.1);
      if !match_probs.contains_key(&forward_pos_pair) {
        continue;
      }
      if !match_probs.contains_key(&backward_pos_pair) {
        continue;
      }
      match forward_pos_pairs.get_mut(&backward_pos_pair) {
        Some(x) => {
          x.insert(forward_pos_pair);
        }
        None => {
          let mut x = PosPairMat::<T>::default();
          x.insert(forward_pos_pair);
          forward_pos_pairs.insert(backward_pos_pair, x);
        }
      }
      match backward_pos_pairs.get_mut(&forward_pos_pair) {
        Some(x) => {
          x.insert(backward_pos_pair);
        }
        None => {
          let mut x = PosPairMat::<T>::default();
          x.insert(backward_pos_pair);
          backward_pos_pairs.insert(forward_pos_pair, x);
        }
      }
      pos_quads.insert((pos_pair.0, pos_pair.1, pos_pair2.0, pos_pair2.1));
      let substr_len_pair = (
        pos_pair.1 - pos_pair.0 + T::one(),
        pos_pair2.1 - pos_pair2.0 + T::one(),
      );
      let left_pos_pair = (pos_pair.0, pos_pair2.0);
      match pos_quads_hashed_lens.get_mut(&substr_len_pair) {
        Some(x) => {
          x.insert(left_pos_pair);
        }
        None => {
          let mut x = PosPairMat::<T>::default();
          x.insert(left_pos_pair);
          pos_quads_hashed_lens.insert(substr_len_pair, x);
        }
      }
    }
  }
  for x in match_probs.keys() {
    match matchable_poss.get_mut(&x.0) {
      Some(y) => {
        y.insert(x.1);
      }
      None => {
        let mut y = SparsePoss::<T>::default();
        y.insert(x.1);
        matchable_poss.insert(x.0, y);
      }
    }
    match matchable_poss2.get_mut(&x.1) {
      Some(y) => {
        y.insert(x.0);
      }
      None => {
        let mut y = SparsePoss::<T>::default();
        y.insert(x.0);
        matchable_poss2.insert(x.1, y);
      }
    }
  }
  let mut matchable_poss_leftmost = SparsePoss::<T>::default();
  matchable_poss_leftmost.insert(T::zero());
  matchable_poss.insert(T::zero(), matchable_poss_leftmost.clone());
  matchable_poss2.insert(T::zero(), matchable_poss_leftmost);
  let mut matchable_poss_rightmost = SparsePoss::<T>::default();
  matchable_poss_rightmost.insert(seq_len_pair.1 - T::one());
  matchable_poss.insert(seq_len_pair.0 - T::one(), matchable_poss_rightmost);
  let mut matchable_poss_rightmost = SparsePoss::<T>::default();
  matchable_poss_rightmost.insert(seq_len_pair.0 - T::one());
  matchable_poss2.insert(seq_len_pair.1 - T::one(), matchable_poss_rightmost);
  (
    forward_pos_pairs,
    backward_pos_pairs,
    pos_quads,
    pos_quads_hashed_lens,
    matchable_poss,
    matchable_poss2,
  )
}

pub fn write_readme(output_dir_path: &Path, readme_contents: &String) {
  let readme_file_path = output_dir_path.join(README_FILE);
  let mut writer = BufWriter::new(File::create(readme_file_path).unwrap());
  let buf = String::from(readme_contents);
  let _ = writer.write_all(buf.as_bytes());
}
