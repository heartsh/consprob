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
  ss_file_path = asset_dir_path + "/sampled_trnas/sss.dat"
  ss_strings = utils.get_ss_strings(ss_file_path)
  sss = utils.get_sss(ss_strings)
  bpp_mat_on_ss_file_path = asset_dir_path + "/sampled_trnas/bpp_mats_on_ss.dat"
  bpp_mats_on_ss = utils.get_bpp_mats(bpp_mat_on_ss_file_path, seq_lens)
  print(colored.black("(A)"))
  utils.print_color_coded_sss(sss, ss_strings, bpp_mats_on_ss, records, num_of_records)
  bpp_mat_on_sta_file_path = asset_dir_path + "/sampled_trnas/bpp_mats_on_sta.dat"
  bpp_mats_on_sta = utils.get_bpp_mats(bpp_mat_on_sta_file_path, seq_lens)
  print(colored.black("(B)"))
  utils.print_color_coded_sss(sss, ss_strings, bpp_mats_on_sta, records, num_of_records)
  print(
    colored.red("p >= %.2e, " % (0.5))
    + colored.yellow("p >= %.2e, " % (0.5 ** 2))
    + colored.green("p >= %.2e, " % (0.5 ** 3))
    + colored.cyan("p >= %.2e, " % (0.5 ** 4))
    + colored.blue("p < %.2e" % (0.5 ** 4))
    + colored.black(" where p is a\n posterior BPP for a base pair")
  )

if __name__ == "__main__":
  main()
