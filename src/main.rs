extern crate consprob;

use consprob::*;
use std::env;

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt(
    "i",
    "input_file_path",
    "An input FASTA file path containing RNA sequences to predict probabilities",
    "STR",
  );
  opts.reqopt("o", "output_dir_path", "An output directory path", "STR");
  opts.optopt(
    "",
    "min_base_pair_prob",
    &format!(
      "A minimum base-pairing probability (Use {} by default)",
      DEFAULT_MIN_BPP
    ),
    "FLOAT",
  );
  opts.optopt(
    "",
    "min_align_prob",
    &format!(
      "A minimum aligning probability (Use {} by default)",
      DEFAULT_MIN_ALIGN_PROB
    ),
    "FLOAT",
  );
  opts.optopt(
    "t",
    "num_of_threads",
    "The number of threads in multithreading (Use all the threads of this computer by default)",
    "UINT",
  );
  opts.optflag(
    "s",
    "produce_struct_profs",
    &format!("Also compute RNA structural context profiles"),
  );
  opts.optflag(
    "a",
    "produce_align_probs",
    &format!("Also compute nucleotide alignment probabilities"),
  );
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1..]) {
    Ok(opt) => opt,
    Err(failure) => {
      print_program_usage(&program_name, &opts);
      panic!("{}", failure.to_string())
    }
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").unwrap();
  let input_file_path = Path::new(&input_file_path);
  let min_bpp = if matches.opt_present("min_base_pair_prob") {
    matches
      .opt_str("min_base_pair_prob")
      .unwrap()
      .parse()
      .unwrap()
  } else {
    DEFAULT_MIN_BPP
  };
  let min_align_prob = if matches.opt_present("min_align_prob") {
    matches.opt_str("min_align_prob").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIN_ALIGN_PROB
  };
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let produce_struct_profs = matches.opt_present("s");
  let produce_align_probs = matches.opt_present("a");
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
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
  let mut thread_pool = Pool::new(num_of_threads);
  if max_seq_len <= u8::MAX as usize {
    let (prob_mat_sets, align_prob_mat_sets_with_rna_id_pairs) = consprob::<u8>(
      &mut thread_pool,
      &fasta_records,
      min_bpp,
      min_align_prob,
      produce_struct_profs,
      produce_align_probs,
      &align_feature_score_sets,
    );
    write_prob_mat_sets(
      &output_dir_path,
      &prob_mat_sets,
      produce_struct_profs,
      &align_prob_mat_sets_with_rna_id_pairs,
      produce_align_probs,
    );
  } else {
    let (prob_mat_sets, align_prob_mat_sets_with_rna_id_pairs) = consprob::<u16>(
      &mut thread_pool,
      &fasta_records,
      min_bpp,
      min_align_prob,
      produce_struct_profs,
      produce_align_probs,
      &align_feature_score_sets,
    );
    write_prob_mat_sets(
      &output_dir_path,
      &prob_mat_sets,
      produce_struct_profs,
      &align_prob_mat_sets_with_rna_id_pairs,
      produce_align_probs,
    );
  }
  write_readme(output_dir_path, &String::from(README_CONTENTS));
}
