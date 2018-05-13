#! /usr/bin/env python

import utils
from Bio import SeqIO
from Bio import AlignIO
import numpy
import seaborn
from matplotlib import pyplot
import os
from clint.textui import colored

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  utils.init_matplotlib()
  rna_file_path = asset_dir_path + "/sampled_trnas.fa"
  records = [record for record in SeqIO.parse(rna_file_path, "fasta")]
  seq_lens = [len(record.seq) for record in records]
  num_of_records = len(records)
  sa_file_path = asset_dir_path + "/sampled_trnas/sa_with_mafft_ginsi.aln"
  sa = AlignIO.read(sa_file_path, "clustal")
  sa_len = sa.get_alignment_length()
  nums_of_gaps_in_front_of_chars = utils.get_nums_of_gaps_in_front_of_chars(sa, num_of_records, sa_len)
  bap_mat_on_sa_file_path = asset_dir_path + "/sampled_trnas/bap_mats_on_sa.dat"
  bap_mats_on_sa = utils.get_bap_mats(bap_mat_on_sa_file_path, seq_lens)
  uabp_seq_pairs_on_sa = utils.get_unaligned_base_prob_seq_pairs(bap_mats_on_sa, num_of_records, seq_lens)
  print(colored.black("(A)"))
  utils.print_color_coded_sa(sa, bap_mats_on_sa, uabp_seq_pairs_on_sa, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  bap_mat_on_sta_file_path = asset_dir_path + "/sampled_trnas/bap_mats_on_sta.dat"
  bap_mats_on_sta = utils.get_bap_mats(bap_mat_on_sta_file_path, seq_lens)
  uabp_seq_pairs_on_sta = utils.get_unaligned_base_prob_seq_pairs(bap_mats_on_sta, num_of_records, seq_lens)
  print(colored.black("(B)"))
  utils.print_color_coded_sa(sa, bap_mats_on_sta, uabp_seq_pairs_on_sta, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  sa_file_path = asset_dir_path + "/sampled_trnas/sa_with_mafft_xinsi.aln"
  sa = AlignIO.read(sa_file_path, "clustal")
  sa_len = sa.get_alignment_length()
  nums_of_gaps_in_front_of_chars = utils.get_nums_of_gaps_in_front_of_chars(sa, num_of_records, sa_len)
  print(colored.black("(C)"))
  utils.print_color_coded_sa(sa, bap_mats_on_sa, uabp_seq_pairs_on_sa, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  print(colored.black("(D)"))
  utils.print_color_coded_sa(sa, bap_mats_on_sta, uabp_seq_pairs_on_sta, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  sa_file_path = asset_dir_path + "/sampled_trnas/ref_sa.aln"
  sa = AlignIO.read(sa_file_path, "clustal")
  sa_len = sa.get_alignment_length()
  nums_of_gaps_in_front_of_chars = utils.get_nums_of_gaps_in_front_of_chars(sa, num_of_records, sa_len)
  print(colored.black("(E)"))
  utils.print_color_coded_sa(sa, bap_mats_on_sa, uabp_seq_pairs_on_sa, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  print(colored.black("(F)"))
  utils.print_color_coded_sa(sa, bap_mats_on_sta, uabp_seq_pairs_on_sta, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  print(
    colored.red("p >= %.0e, " % (0.1 ** 6))
    + colored.yellow("p >= %.0e, " % (0.1 ** 9))
    + colored.green("p >= %.0e, " % (0.1 ** 12))
    + colored.cyan("p >= %.0e, " % (0.1 ** 15))
    + colored.blue("p < %.0e" % (0.1 ** 15))
    + colored.black(" where p is the product of the BAPs\nfor the base pairs and unaligned base probabilities for the bases aligned with base gaps in a\ncolor-coded column")
  )

if __name__ == "__main__":
  main()
