extern crate consprob;

use consprob::*;

#[test]
fn test_consprob() {
  let fasta_file_reader = Reader::from_file(Path::new(&EXAMPLE_FASTA_FILE_PATH)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for x in fasta_file_reader.records() {
    let x = x.unwrap();
    let mut y = bytes2seq(x.seq());
    y.insert(0, PSEUDO_BASE);
    y.push(PSEUDO_BASE);
    let z = y.len();
    if z > max_seq_len {
      max_seq_len = z;
    }
    fasta_records.push(FastaRecord::new(String::from(x.id()), y));
  }
  let mut align_scores = AlignScores::new(0.);
  align_scores.transfer();
  let seqs = fasta_records.iter().map(|x| &x.seq[..]).collect();
  let num_threads = num_cpus::get() as NumThreads;
  let mut thread_pool = Pool::new(num_threads);
  let produces_struct_profs = true;
  let produces_match_probs = true;
  let (alignfold_prob_mats_avg, match_probs_hashed_ids) = consprob::<u8>(
    &mut thread_pool,
    &seqs,
    DEFAULT_MIN_BASEPAIR_PROB,
    DEFAULT_MIN_MATCH_PROB,
    produces_struct_profs,
    produces_match_probs,
    &align_scores,
  );
  for alignfold_probs_avg in &alignfold_prob_mats_avg {
    for x in alignfold_probs_avg.basepair_probs.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(x));
    }
    for x in alignfold_probs_avg.context_profs.iter() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(x));
    }
  }
  for match_probs in match_probs_hashed_ids.values() {
    for x in match_probs.loopmatch_probs.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(x));
    }
    for x in match_probs.pairmatch_probs.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(x));
    }
    for x in match_probs.match_probs.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(x));
    }
  }
}
