#! /usr/bin/env python

import utils
from Bio import SeqIO
from Bio import AlignIO
import os

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  rna_file_path = asset_dir_path + "/sampled_trnas.fa"
  records = [record for record in SeqIO.parse(rna_file_path, "fasta")]
  record_ids = [record.id for record in records]
  input_sa_file_path = asset_dir_path + "/trna.sth"
  sa = AlignIO.read(input_sa_file_path, "stockholm")
  records_other_than_sampled_rnas_are_removed = False
  while not records_other_than_sampled_rnas_are_removed:
    num_of_records = len(sa)
    for (i, record) in sa.enumerate():
      if record.id in record_ids:
        sa
  print(sa)
  output_sa_file_path = asset_dir_path + "/sampled_trnas/ref_sa.aln"
  AlignIO.write(sa, output_sa_file_path, "clustal")
  # sa_len = sa.get_alignment_length()
  # for i in 0 .. sa_len:
    # for j in 0 .. num_of_records:

if __name__ == "__main__":
  main()
