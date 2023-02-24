extern crate consprob;
extern crate criterion;

use consprob::*;
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_consprob(criterion: &mut Criterion) {
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
  criterion.bench_function("consprob::<u8>", |x| {
    x.iter(|| {
      let _ = consprob::<u8>(
        &mut thread_pool,
        &seqs,
        DEFAULT_MIN_BASEPAIR_PROB,
        DEFAULT_MIN_MATCH_PROB,
        produces_struct_profs,
        produces_match_probs,
        &align_scores,
      );
    });
  });
}

criterion_group!(benches, bench_consprob);
criterion_main!(benches);
