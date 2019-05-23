import os
import matplotlib
from matplotlib import pylab
import numpy
from math import log
from clint.textui import colored
from itertools import combinations

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)

def init_matplotlib(): 
  params = {
    "legend.fontsize": "x-large",
    "axes.labelsize": "x-large",
    "axes.titlesize":"x-large",
    "xtick.labelsize":"x-large",
    "ytick.labelsize":"x-large"
  }
  pylab.rcParams.update(params)
  matplotlib.rcParams['ps.fonttype'] = 42

def get_nums_of_gaps_in_front_of_chars(sa, num_of_records, sa_len):
  nums_of_gaps_in_front_of_chars = numpy.zeros((num_of_records, sa_len), dtype = int)
  for i in range(0, num_of_records):
    for j in range(0, sa_len):
      if j == 0:
        if sa[i].seq[j] == "-":
          nums_of_gaps_in_front_of_chars[i][j] += 1
      else:
        if sa[i].seq[j] == "-":
          nums_of_gaps_in_front_of_chars[i][j] = nums_of_gaps_in_front_of_chars[i][j - 1] + 1
        else:
          nums_of_gaps_in_front_of_chars[i][j] = nums_of_gaps_in_front_of_chars[i][j - 1]
  return nums_of_gaps_in_front_of_chars

def get_css_string(css_file_path):
  css_file = open(css_file_path)
  line = css_file.readlines()[-2]
  css_string = line.split()[2]
  return css_string

def get_ss(ss_string):
  ss = []
  stack = []
  for (i, char) in enumerate(ss_string):
    if char == "(":
      stack.append(i)
    elif char == ")":
      ss.insert(0, (stack.pop(), i))
  return ss

def print_color_coded_css_with_sa(css, css_string, sa, bpap_mats, nums_of_gaps_in_front_of_chars, num_of_records, sa_len):
  color_coded_css_with_sa = [list(map(colored.black, sa[i].seq.upper())) for i in range(0, num_of_records)]
  color_coded_css_with_sa.append(list(map(colored.black, css_string)))
  combination_num = len(list(combinations(range(0, num_of_records), 2)))
  for (i, j) in css:
    mean_bpap = 0
    for k in range(0, num_of_records):
      pos_without_gaps_1 = i - nums_of_gaps_in_front_of_chars[k][i]
      pos_without_gaps_2 = j - nums_of_gaps_in_front_of_chars[k][j]
      for l in range(k + 1, num_of_records):
        pos_without_gaps_3 = i - nums_of_gaps_in_front_of_chars[l][i]
        pos_without_gaps_4 = j - nums_of_gaps_in_front_of_chars[l][j]
        bpap = bpap_mats[(k, l)][pos_without_gaps_1][pos_without_gaps_2][pos_without_gaps_3][pos_without_gaps_4]
        if bpap > 0:
          mean_bpap += bpap
    mean_bpap /= combination_num
    for k in range(0, num_of_records + 1):
      for l in (i, j):
        char = colored.clean(color_coded_css_with_sa[k][l])
        color_coded_char = colored.blue(char)
        if mean_bpap >= 0.25:
          color_coded_char = colored.red(char)
        elif mean_bpap >= 0.25 ** 2:
          color_coded_char = colored.yellow(char)
        elif mean_bpap >= 0.25 ** 3:
          color_coded_char = colored.green(char)
        elif mean_bpap >= 0.25 ** 4:
          color_coded_char = colored.cyan(char)
        color_coded_css_with_sa[k][l] = color_coded_char
  for string in color_coded_css_with_sa:
    color_coded_string = ""
    for char in string:
      color_coded_string += char
    print(color_coded_string)

def get_bpap_mats(bpap_mat_file_path, seq_lens):
  bpap_mats = {}
  bpap_mat_file = open(bpap_mat_file_path)
  lines = bpap_mat_file.readlines()
  lines = [line for line in lines if line[0].isdigit() or line[0].startswith(">")]
  num_of_lines = len(lines)
  for i in range(0, num_of_lines - 1, 2):
    tail = lines[i][1 :].split(",")
    rna_id_pair = (int(tail[0]), int(tail[1]))
    seq_len_pair = (seq_lens[rna_id_pair[0]], seq_lens[rna_id_pair[1]])
    bpap_mat = numpy.zeros((seq_len_pair[0],seq_len_pair[0], seq_len_pair[1],  seq_len_pair[1]))
    for string in lines[i + 1].strip().split(" "):
      substrings = string.split(",")
      (j, k, l, m, bpap) = (int(substrings[0]), int(substrings[1]), int(substrings[2]), int(substrings[3]), float(substrings[4]))
      bpap_mat[j, k, l, m] = bpap
    bpap_mats[rna_id_pair] = bpap_mat
  return bpap_mats

def get_bpp_mats(bpp_mat_file_path, seq_lens):
  bpp_mats = {}
  bpp_mat_file = open(bpp_mat_file_path)
  lines = bpp_mat_file.readlines()
  lines = [line for line in lines if line[0].isdigit() or line[0].startswith(">")]
  num_of_lines = len(lines)
  for i in range(0, num_of_lines - 1, 2):
    rna_id = int(lines[i][1 :])
    seq_len = seq_lens[rna_id]
    bpp_mat = numpy.zeros((seq_len, seq_len))
    for string in lines[i + 1].strip().split(" "):
      substrings = string.split(",")
      (j, k, bpp) = (int(substrings[0]), int(substrings[1]), float(substrings[2]))
      bpp_mat[j, k] = bpp
    bpp_mats[rna_id] = bpp_mat
  return bpp_mats

def print_color_coded_sss(sss, ss_strings, bpp_mats, records, num_of_records):
  color_coded_seqs = [list(map(colored.black, record.seq)) for record in records]
  color_coded_sss = [list(map(colored.black, ss_string)) for ss_string in ss_strings]
  for i in range(0, num_of_records):
    for (j, k) in sss[i]:
      bpp = bpp_mats[i][j, k]
      for l in (j, k):
        char_pair = (colored.clean(color_coded_seqs[i][l]), colored.clean(color_coded_sss[i][l]))
        color_coded_char_pair = (colored.blue(char_pair[0]), colored.blue(char_pair[1]))
        if bpp >= 0.5:
          color_coded_char_pair = (colored.red(char_pair[0]), colored.red(char_pair[1]))
        elif bpp >= 0.5 ** 2:
          color_coded_char_pair = (colored.yellow(char_pair[0]), colored.yellow(char_pair[1]))
        elif bpp >= 0.5 ** 3:
          color_coded_char_pair = (colored.green(char_pair[0]), colored.green(char_pair[1]))
        elif bpp >= 0.5 ** 4:
          color_coded_char_pair = (colored.cyan(char_pair[0]), colored.cyan(char_pair[1]))
        color_coded_seqs[i][l] = color_coded_char_pair[0]
        color_coded_sss[i][l] = color_coded_char_pair[1]
  for seq, ss in zip(color_coded_seqs, color_coded_sss):
    color_coded_seq = ""
    for char in seq:
      color_coded_seq += char
    color_coded_ss = ""
    for char in ss:
      color_coded_ss += char
    print(color_coded_seq)
    print(color_coded_ss)

def get_ss_strings(ss_file_path):
  ss_strings = []
  ss_file = open(ss_file_path)
  lines = ss_file.readlines()
  num_of_lines = len(lines)
  for i in range(0, num_of_lines - 1, 7):
    ss_string = lines[i + 5].split()[0]
    ss_strings.append(ss_string)
  return ss_strings

def get_sss(ss_strings):
  return list(map(get_ss, ss_strings))
