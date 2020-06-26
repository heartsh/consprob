extern crate consprob;

use consprob::*;
use std::env;

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences to predict probabilities", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "offset_4_max_gap_num", &format!("An offset for maximum numbers of gaps (Uses {} by default)", DEFAULT_OFFSET_4_MAX_GAP_NUM), "UINT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
  opts.optflag("u", "is_posterior_model", "Uses posterior model to score secondary structures (Not recommended due to poor accuracy)");
  opts.optflag("a", "produces_access_probs", &format!("Also compute accessible probabilities (only for Turner model)"));
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").unwrap();
  let input_file_path = Path::new(&input_file_path);
  let min_bpp = if matches.opt_present("min_base_pair_prob") {
    matches.opt_str("min_base_pair_prob").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIN_BPP
  };
  let offset_4_max_gap_num = if matches.opt_present("offset_4_max_gap_num") {
    matches.opt_str("offset_4_max_gap_num").unwrap().parse().unwrap()
  } else {
    DEFAULT_OFFSET_4_MAX_GAP_NUM
  };
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let is_posterior_model = matches.opt_present("u");
  let produces_access_probs = matches.opt_present("a") & !is_posterior_model;
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = convert(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut thread_pool = Pool::new(num_of_threads);
  let prob_mat_sets = consprob(&mut thread_pool, &fasta_records, min_bpp, offset_4_max_gap_num, is_posterior_model, produces_access_probs);
  write_prob_mat_sets(&output_dir_path, &prob_mat_sets, produces_access_probs);
}
