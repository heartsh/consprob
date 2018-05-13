#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  utils.init_matplotlib()
  rna_file_path = asset_dir_path + "/sampled_trnas.fa"
  records = [record for record in SeqIO.parse(rna_file_path, "fasta")]
  seq_lens = [len(record.seq) for record in records]
  num_of_records = len(records)
  fig_size_in_inch = 8.3
  (_, axes) = pyplot.subplots(nrows = 2, ncols = num_of_records, figsize = (fig_size_in_inch, fig_size_in_inch / 3))
  bpp_mat_on_ss_file_path = asset_dir_path + "/sampled_trnas/bpp_mats_on_ss.dat"
  bpp_mats_on_ss = utils.get_bpp_mats(bpp_mat_on_ss_file_path, seq_lens)
  bpp_mat_on_sta_file_path = asset_dir_path + "/sampled_trnas/bpp_mats_on_sta.dat"
  bpp_mats_on_sta = utils.get_bpp_mats(bpp_mat_on_sta_file_path, seq_lens)
  for rna_id in range(0, num_of_records):
    seaborn.heatmap(bpp_mats_on_ss[rna_id], ax = axes[0][rna_id], xticklabels = False, yticklabels = False, cbar = False)
    seaborn.heatmap(bpp_mats_on_sta[rna_id], ax = axes[1][rna_id], xticklabels = False, yticklabels = False, cbar = False)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/bpp_mats_on_ss_and_bpp_mats_on_sta.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
