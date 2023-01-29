extern crate consprob;
extern crate criterion;

use consprob::*;
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_consprob(criterion: &mut Criterion) {
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
  criterion.bench_function("consprob::<u8>", |b| {
    b.iter(|| {
      let _ = consprob::<u8>(
        &mut thread_pool,
        &seqs,
        DEFAULT_MIN_BPP,
        DEFAULT_MIN_ALIGN_PROB,
        true,
        true,
        &align_feature_score_sets,
      );
    });
  });
}

criterion_group!(benches, bench_consprob);
criterion_main!(benches);
