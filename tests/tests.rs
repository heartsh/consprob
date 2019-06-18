extern crate rnafamprob;
extern crate time;

use time::precise_time_s;
use std::process::{Command, Output};

#[test]
fn test_rnafam_program() {
  let input_fasta_file_path = "assets/sampled_trnas.fa";
  let args = ["-i", input_fasta_file_path, "-o", "assets/sampled_trnas"];
  let begin = precise_time_s();
  run_command("target/release/rnafamprob", &args, "Failed to run the RNAfamProb program.");
  let elapsed_time = precise_time_s() - begin;
  println!("The elapsed time for computing the posterior probabilities on structural alignment given each pair of RNA sequences in the FASTA file \"{}\" = {} [s].", input_fasta_file_path, elapsed_time);
}

#[inline]
pub fn run_command(command: &str, args: &[&str], expect: &str) -> Output {
  Command::new(command).args(args).output().expect(expect)
}
