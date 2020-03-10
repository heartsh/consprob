extern crate phyloprob;
extern crate bio;
extern crate num_cpus;

use phyloprob::*;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::fs::create_dir;

const BPP_MAT_ON_STA_FILE_NAME: &'static str = "bpp_mats_on_sta.dat";
const UPP_MAT_ON_STA_FILE_NAME: &'static str = "upp_mats_on_sta.dat";
const VERSION: &'static str = "0.1.1";

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "opening_gap_penalty", &format!("An opening-gap penalty (Uses {} by default)", DEFAULT_OPENING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "extending_gap_penalty", &format!("An extending-gap penalty (Uses {} by default)", DEFAULT_EXTENDING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "offset_4_max_gap_num", &format!("An offset for maximum numbers of gaps (Uses {} by default)", DEFAULT_OFFSET_4_MAX_GAP_NUM), "UINT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of the threads of this computer by default)", "UINT");
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
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
  let opening_gap_penalty = if matches.opt_present("opening_gap_penalty") {
    matches.opt_str("opening_gap_penalty").unwrap().parse().unwrap()
  } else {
    DEFAULT_OPENING_GAP_PENALTY
  };
  let extending_gap_penalty = if matches.opt_present("extending_gap_penalty") {
    matches.opt_str("extending_gap_penalty").unwrap().parse().unwrap()
  } else {
    DEFAULT_EXTENDING_GAP_PENALTY
  };
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
  let (bpp_mats, upp_mats) = phyloprob(&mut thread_pool, &fasta_records, opening_gap_penalty, extending_gap_penalty, min_bpp, offset_4_max_gap_num);
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let output_file_header = format!(" in this file = \"{}\".\n; The values of the parameters used to the matrices are as follows.\n; \"opening_gap_penalty\" = {}, \"extending_gap_penalty\" = {}, \"min_bpp\" = {}, \"offset_4_max_gap_num\" = {}, \"num_of_threads\" = {}.", input_file_path.display(), opening_gap_penalty, extending_gap_penalty, min_bpp, offset_4_max_gap_num, num_of_threads);
  let bpp_mat_on_sta_file_path = output_dir_path.join(BPP_MAT_ON_STA_FILE_NAME);
  let upp_mat_on_sta_file_path = output_dir_path.join(UPP_MAT_ON_STA_FILE_NAME);
  let mut writer_2_bpp_mat_on_sta_file = BufWriter::new(File::create(bpp_mat_on_sta_file_path).unwrap());
  let mut buf_4_writer_2_bpp_mat_on_sta_file = format!("; The version {} of the PhyloProb program.\n; The path to the input file in order to compute the average base-pairing probability matrices on structural alignment in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of the input RNA sequence correpsonding to the row. The next row to the row is with the average base-pairing probability matrix on structural alignment of the sequence.";
  for (rna_id, bpp_mat) in bpp_mats.iter().enumerate() {
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (&(i, j), &bpp) in bpp_mat.iter() {
      buf_4_rna_id.push_str(&format!("{},{},{} ", i - 1, j - 1, bpp));
    }
    buf_4_writer_2_bpp_mat_on_sta_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_bpp_mat_on_sta_file.write_all(buf_4_writer_2_bpp_mat_on_sta_file.as_bytes());
  let mut writer_2_upp_mat_on_sta_file = BufWriter::new(File::create(upp_mat_on_sta_file_path).unwrap());
  let mut buf_4_writer_2_upp_mat_on_sta_file = format!("; The version {} of the PhyloProb program.\n; The path to the input file in order to compute the average unpairing probability matrices on structural alignment in this file", VERSION) + &output_file_header + "\n; Each row beginning with \">\" is with the ID of an RNA sequence. The next row to the row is with the average unpairing probability matrix on structural alignment of the sequence.";
  for (rna_id, upp_mat) in upp_mats.iter().enumerate() {
    let seq_len = fasta_records[rna_id].seq.len();
    let mut buf_4_rna_id = format!("\n\n>{}\n", rna_id);
    for (i, &upp) in upp_mat.iter().enumerate() {
      if i == 0 || i == seq_len - 1 {continue;}
      buf_4_rna_id.push_str(&format!("{},{} ", i - 1, upp));
    }
    buf_4_writer_2_upp_mat_on_sta_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_upp_mat_on_sta_file.write_all(buf_4_writer_2_upp_mat_on_sta_file.as_bytes());
}
