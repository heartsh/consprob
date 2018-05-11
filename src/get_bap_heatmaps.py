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
  (_, axes) = pyplot.subplots(nrows = num_of_records, ncols = num_of_records, figsize = (fig_size_in_inch, fig_size_in_inch))
  bap_mat_on_sa_file_path = asset_dir_path + "/sampled_trnas/bap_mats_on_sa.dat"
  bap_mats_on_sa = utils.get_bap_mats(bap_mat_on_sa_file_path, seq_lens)
  bap_mat_on_sta_file_path = asset_dir_path + "/sampled_trnas/bap_mats_on_sta.dat"
  bap_mats_on_sta = utils.get_bap_mats(bap_mat_on_sta_file_path, seq_lens)
  for rna_id_1 in range(0, num_of_records):
    axes[num_of_records - 1 - rna_id_1][rna_id_1].set_visible(False)
    for rna_id_2 in range(rna_id_1 + 1, num_of_records):
      seaborn.heatmap(bap_mats_on_sa[(rna_id_1, rna_id_2)], ax = axes[num_of_records - 1 - rna_id_2][rna_id_1], xticklabels = False, yticklabels = False, cbar = False)
      seaborn.heatmap(bap_mats_on_sta[(rna_id_1, rna_id_2)].transpose(), ax = axes[num_of_records - 1 - rna_id_1][rna_id_2], xticklabels = False, yticklabels = False, cbar = False)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/bap_mats_on_sa_and_bap_mats_on_sta.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
