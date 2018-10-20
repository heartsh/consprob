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
  css_file_path = asset_dir_path + "/sampled_trnas/css_from_sa_with_mafft_ginsi.sth"
  css_string = utils.get_css_string(css_file_path)
  css = utils.get_ss(css_string)
  nums_of_gaps_in_front_of_chars = utils.get_nums_of_gaps_in_front_of_chars(sa, num_of_records, sa_len)
  bpap_mat_file_path = asset_dir_path + "/sampled_trnas/bpap_mats_1.dat"
  bpap_mats = utils.get_bpap_mats(bpap_mat_file_path, seq_lens)
  print(colored.black("(A)"))
  utils.print_color_coded_css_with_sa(css, css_string, sa, bpap_mats, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  sa_file_path = asset_dir_path + "/sampled_trnas/sa_with_mafft_xinsi.aln"
  sa = AlignIO.read(sa_file_path, "clustal")
  sa_len = sa.get_alignment_length()
  css_file_path = asset_dir_path + "/sampled_trnas/css_from_sa_with_mafft_xinsi.sth"
  css_string = utils.get_css_string(css_file_path)
  css = utils.get_ss(css_string)
  nums_of_gaps_in_front_of_chars = utils.get_nums_of_gaps_in_front_of_chars(sa, num_of_records, sa_len)
  print(colored.black("(B)"))
  utils.print_color_coded_css_with_sa(css, css_string, sa, bpap_mats, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  sa_file_path = asset_dir_path + "/sampled_trnas/ref_sa.aln"
  sa = AlignIO.read(sa_file_path, "clustal")
  sa_len = sa.get_alignment_length()
  css_file_path = asset_dir_path + "/sampled_trnas/ref_css.sth"
  css_string = utils.get_css_string(css_file_path)
  css = utils.get_ss(css_string)
  nums_of_gaps_in_front_of_chars = utils.get_nums_of_gaps_in_front_of_chars(sa, num_of_records, sa_len)
  print(colored.black("(C)"))
  utils.print_color_coded_css_with_sa(css, css_string, sa, bpap_mats, nums_of_gaps_in_front_of_chars, num_of_records, sa_len)
  print(
    colored.red("p >= %.0e, " % (0.1 ** 66))
    + colored.yellow("p >= %.0e, " % (0.1 ** 74))
    + colored.green("p >= %.0e, " % (0.1 ** 82))
    + colored.cyan("p >= %.0e, " % (0.1 ** 90))
    + colored.blue("p < %.0e" % (0.1 ** 90))
    + colored.black(" where p is the product of the STAPs\nfor the base quadruples in a color-coded column-pair")
  )

if __name__ == "__main__":
  main()
