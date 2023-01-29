extern crate consprob;

use consprob::*;

#[test]
fn test_consprob() {
  let fasta_file_reader = Reader::from_file(Path::new(&EXAMPLE_FASTA_FILE_PATH)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = convert(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    let seq_len = seq.len();
    if seq_len > max_seq_len {
      max_seq_len = seq_len;
    }
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut align_feature_score_sets = AlignFeatureCountSets::new(0.);
  align_feature_score_sets.transfer();
  let seqs = fasta_records.iter().map(|x| &x.seq[..]).collect();
  let num_of_threads = num_cpus::get() as NumOfThreads;
  let mut thread_pool = Pool::new(num_of_threads);
  let (prob_mat_sets, align_prob_mat_sets_with_rna_id_pairs) = consprob::<u8>(
    &mut thread_pool,
    &seqs,
    DEFAULT_MIN_BPP,
    DEFAULT_MIN_ALIGN_PROB,
    true,
    true,
    &align_feature_score_sets,
  );
  for prob_mats in &prob_mat_sets {
    for &bpp in prob_mats.bpp_mat.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&bpp));
    }
    for &prob in prob_mats.contexts.iter() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&prob));
    }
  }
  for align_prob_mats in align_prob_mat_sets_with_rna_id_pairs.values() {
    for &align_prob in align_prob_mats.loop_align_prob_mat.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&align_prob));
    }
    for &align_prob in align_prob_mats.basepair_align_prob_mat.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&align_prob));
    }
    for &align_prob in align_prob_mats.align_prob_mat.values() {
      assert!((PROB_BOUND_LOWER..PROB_BOUND_UPPER).contains(&align_prob));
    }
  }
}
