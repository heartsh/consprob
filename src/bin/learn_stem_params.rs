extern crate stem;

use stem::*;
use std::env;
use std::path::Path;
use std::io::prelude::*;
use std::io::{BufWriter, BufReader};
use std::fs::File;

type Col = Pos;
type ColPair = (Col, Col);
type ColPairs = Vec<ColPair>;
type ColPairSeqsWithColPairs = HashMap<ColPair, ColPairs, Hasher>;
type Css = ColPairSeqsWithColPairs;
type Sa = Vec<Seq>;
struct Sta {
  pub cons_second_struct: Css,
  pub seq_align: Sa,
}
type Stas = Vec<Sta>;
type StaEventCount = Prob;
type BaCounts = HashMap<BasePair, StaEventCount, Hasher>;
type GapCounts = HashMap<Base, StaEventCount, Hasher>;
type BpaCounts = HashMap<BaseQuadruple, StaEventCount, Hasher>;
type BpCounts = HashMap<BasePair, StaEventCount, Hasher>;
type NbpCounts = HashMap<Base, StaEventCount, Hasher>;
struct StaEventCountSets {
  pub ba_counts: BaCounts,
  pub og_counts: GapCounts,
  pub eg_counts: GapCounts,
  pub bpa_counts: BpaCounts,
  pub bp_counts: BpCounts,
  pub nbp_counts: NbpCounts,
}

const CSS_REPLACEMENT_PATTERN_CHARS: &'static str = "[{<>}],_-:";
const CSS_REPLACEMENT_DEST_CHARS: &'static str = "((()))....";
const BP_LEFT_BRACKET: char = '(';
const BP_RIGHT_BRACKET: char = ')';
const GAP: Base = '-' as Base;
const DEFAULT_STA_EVENT_PSEUDO_COUNT: StaEventCount = 0.;
const DEFAULT_MAX_INDEL_FREQ_RELATIVE_2_ALL_PROF_FREQS: Prob = 0.2;
const DEFAULT_MIN_BP_FREQ_RELATIVE_2_ALL_PROF_FREQS: Prob = 0.4;

impl Sta {
  pub fn new() -> Sta {
    Sta {
      cons_second_struct: Css::default(),
      seq_align: Sa::new(),
    }
  }
}

impl StaEventCountSets {
  pub fn new(sta_event_pseudo_count: StaEventCount) -> StaEventCountSets {
    let ba_counts = BA_ALPHABET.iter().map(|base_pair| {(*base_pair, sta_event_pseudo_count)}).collect();
    let gap_counts: GapCounts = SEQ_ALPHABET.iter().map(|&base| {(base, sta_event_pseudo_count)}).collect();
    let bpa_counts = BPA_ALPHABET.iter().map(|base_quadruple| {(*base_quadruple, sta_event_pseudo_count)}).collect();
    StaEventCountSets {
      ba_counts: ba_counts,
      og_counts: gap_counts.clone(),
      eg_counts: gap_counts,
      bpa_counts: bpa_counts,
      bp_counts: BP_ALPHABET.keys().map(|base_pair| {(*base_pair, sta_event_pseudo_count)}).collect(),
      nbp_counts: SEQ_ALPHABET.iter().map(|&base| {(base, sta_event_pseudo_count)}).collect(),
    }
  }
}

fn main() {
  let args = env::args().collect::<Args>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input Stockholm file containing RNA structural alignments", "STR");
  opts.reqopt("o", "output_file_path", "The path to an output file", "STR");
  opts.optopt("p", "struct_align_event_pseudo_count", "An structural-alignment event pseudo-count", "UINT");
  opts.optopt("", "max_indel_freq_relative_2_all_prof_freqs", "A maximum indel-frequency relative to all profile frequencies", "FLOAT");
  opts.optopt("", "min_base_pair_freq_relative_2_all_prof_freqs", "A minimum base-pairing-frequency relative to all profile frequencies", "FLOAT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_file_path = opts.opt_str("i").expect("Failed to get the path to an input Stockholm file containing RNA structural alignments from command arguments.");
  let input_file_path = Path::new(&input_file_path);
  let output_file_path = opts.opt_str("o").expect("Failed to get the path to an output file from command arguments.");
  let output_file_path = Path::new(&output_file_path);
  let sta_event_pseudo_count = if opts.opt_present("struct_align_event_pseudo_count") {
    opts.opt_str("struct_align_event_pseudo_count").expect("Failed to get a structural-alignment event pseudo-count from command arguments.").parse().expect("Failed to parse a structural-alignment event pseudo-count.")
  } else {
    DEFAULT_STA_EVENT_PSEUDO_COUNT
  };
  let min_bp_freq_relative_2_all_prof_freqs = if opts.opt_present("min_base_pair_freq_relative_2_all_prof_freqs") {
    opts.opt_str("min_base_pair_freq_relative_2_all_prof_freqs").expect("Failed to get a minimum base-pairing-frequency relative to all profile frequencies from command arguments.").parse().expect("Failed to parse a minimum base-pairing-frequency relative to all profile frequencies.")
  } else {
    DEFAULT_MIN_BP_FREQ_RELATIVE_2_ALL_PROF_FREQS
  };
  let max_indel_freq_relative_2_all_prof_freqs = if opts.opt_present("max_indel_freq_relative_2_all_prof_freqs") {
    opts.opt_str("max_indel_freq_relative_2_all_prof_freqs").expect("Failed to get a maximum indel-frequency relative to all profile frequencies from command arguments.").parse().expect("Failed to parse a maximum indel-frequency relative to all profile frequencies.")
  } else {
    DEFAULT_MAX_INDEL_FREQ_RELATIVE_2_ALL_PROF_FREQS
  };
  let reader_2_sta_file = BufReader::new(File::open(input_file_path).expect("Failed to read an input file."));
  let mut stas = Stas::new();
  let mut sta = Sta::new();
  let mut is_sta_valid = true;
  for line in reader_2_sta_file.lines() {
    let line = line.expect("Failed to get a line.");
    if line.is_empty() || line == "# STOCKHOLM 1.0" {
      continue;
    } else if line == "//" {
      if is_sta_valid {
        stas.push(sta);
      }
      sta = Sta::new();
      is_sta_valid = true;
    } else if line.starts_with("#=GC SS_cons") {
      let sta_prof = line.split_whitespace().skip(2).next().expect("Failed to split a string.");
      let bp_freq_relative_2_all_prof_freqs = sta_prof.chars().filter(|&prof_char| prof_char == '(' || prof_char == ')').count() as Prob / sta_prof.len() as Prob;
      let indel_freq_relative_2_all_prof_freqs = sta_prof.chars().filter(|&prof_char| prof_char == '.').count() as Prob / sta_prof.len() as Prob;
      is_sta_valid = (bp_freq_relative_2_all_prof_freqs >= min_bp_freq_relative_2_all_prof_freqs) & (indel_freq_relative_2_all_prof_freqs <= max_indel_freq_relative_2_all_prof_freqs) & is_sta_valid;
      sta.cons_second_struct = get_css(&convert_css_str(&(String::from("(") + sta_prof + ")")));
    } else {
      if is_sta_valid {
        let mut sa_mem = vec![PSEUDO_BASE];
        let temp_sa_mem = line.split_whitespace().skip(1).next().expect("Failed to split a string.").as_bytes();
        is_sta_valid = temp_sa_mem.iter().fold(true, |acc, &base| {acc & (is_rna_base(base) || base == GAP)}) & is_sta_valid;
        if is_sta_valid {
          sa_mem.extend(temp_sa_mem);
          sa_mem.push(PSEUDO_BASE);
          sta.seq_align.push(sa_mem);
        }
      }
    }
  }
  let mut sta_event_count_sets = StaEventCountSets::new(sta_event_pseudo_count);
  for sta in &stas {
    let num_of_sa_rows = sta.seq_align.len();
    let num_of_sa_cols = sta.seq_align[0].len();
    let mut col_pair_stack = vec![(0, num_of_sa_cols - 1)];
    while col_pair_stack.len() > 0 {
      let col_pair_1 = col_pair_stack.pop().expect("Failed to pop a vector.");
      let (i, j) = col_pair_1;
      match sta.cons_second_struct.get(&col_pair_1) {
        Some(col_pairs) => {
          let mut start = i + 1;
          for &(k, l) in col_pairs {
            col_pair_stack.push((k, l));
            for m in start .. k {
              for n in 0 .. num_of_sa_rows {
                let base_1 = sta.seq_align[n][m];
                let front_base_1 = sta.seq_align[n][m - 1];
                let behind_base_1 = sta.seq_align[n][m + 1];
                if base_1 != GAP {
                  *sta_event_count_sets.nbp_counts.get_mut(&base_1).expect("Failed to get an element from a hash map.") += 1.;
                }
                for o in n + 1 .. num_of_sa_rows {
                  let base_2 = sta.seq_align[o][m];
                  let front_base_2 = sta.seq_align[o][m - 1];
                  let behind_base_2 = sta.seq_align[o][m + 1];
                  let base_pair = get_ordered_base_pair(&(base_1, base_2));
                  let base_not_being_gap = base_pair.1;
                  if base_1 != GAP && base_2 != GAP {
                    *sta_event_count_sets.ba_counts.get_mut(&base_pair).expect("Failed to get an element from a hash map.") += 1.;
                  } else if base_1 != GAP {
                    if front_base_2 == GAP && behind_base_2 == GAP {
                      *sta_event_count_sets.eg_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                    } else {
                      *sta_event_count_sets.og_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                    }
                  } else if base_2 != GAP {
                    if front_base_1 == GAP && behind_base_1 == GAP {
                      *sta_event_count_sets.eg_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                    } else {
                      *sta_event_count_sets.og_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                    }
                  }
                }
              }
            }
            for m in 0 .. num_of_sa_rows {
              let base_1 = sta.seq_align[m][k];
              let base_2 = sta.seq_align[m][l];
              if base_1 != GAP && base_2 != GAP && sta_event_count_sets.bp_counts.contains_key(&(base_1, base_2)) {
                *sta_event_count_sets.bp_counts.get_mut(&(base_1, base_2)).expect("Failed to get an element from a hash map.") += 1.;
              } else if base_1 != GAP {
                *sta_event_count_sets.nbp_counts.get_mut(&base_1).expect("Failed to get an element from a hash map.") += 1.;
              } else if base_2 != GAP {
                *sta_event_count_sets.nbp_counts.get_mut(&base_2).expect("Failed to get an element from a hash map.") += 1.;
              }
              for n in m + 1 .. num_of_sa_rows {
                let base_3 = sta.seq_align[n][k];
                let base_4 = sta.seq_align[n][l];
                let base_quadruple = get_ordered_base_quadruple(&(base_1, base_2, base_3, base_4));
                if sta_event_count_sets.bpa_counts.contains_key(&base_quadruple) {
                  *sta_event_count_sets.bpa_counts.get_mut(&base_quadruple).expect("Failed to get an element from a hash map.") += 1.;
                }
              }
            }
            start = l + 1;
          }
          for m in start .. j {
            for n in 0 .. num_of_sa_rows {
              let base_1 = sta.seq_align[n][m];
              let front_base_1 = sta.seq_align[n][m - 1];
              let behind_base_1 = sta.seq_align[n][m + 1];
              if base_1 != GAP {
                *sta_event_count_sets.nbp_counts.get_mut(&base_1).expect("Failed to get an element from a hash map.") += 1.;
              }
              for o in n + 1 .. num_of_sa_rows {
                let base_2 = sta.seq_align[o][m];
                let front_base_2 = sta.seq_align[o][m - 1];
                let behind_base_2 = sta.seq_align[o][m + 1];
                let base_pair = get_ordered_base_pair(&(base_1, base_2));
                let base_not_being_gap = base_pair.1;
                if base_1 != GAP && base_2 != GAP {
                  *sta_event_count_sets.ba_counts.get_mut(&base_pair).expect("Failed to get an element from a hash map.") += 1.;
                } else if base_1 != GAP {
                  if front_base_2 == GAP && behind_base_2 == GAP {
                    *sta_event_count_sets.eg_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  } else {
                    *sta_event_count_sets.og_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  }
                } else if base_2 != GAP {
                  if front_base_1 == GAP && behind_base_1 == GAP {
                    *sta_event_count_sets.eg_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  } else {
                    *sta_event_count_sets.og_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  }
                }
              }
            }
          }
        }, None => {
          for m in i + 1 .. j {
            for n in 0 .. num_of_sa_rows {
              let base_1 = sta.seq_align[n][m];
              let front_base_1 = sta.seq_align[n][m - 1];
              let behind_base_1 = sta.seq_align[n][m + 1];
              if base_1 != GAP {
                *sta_event_count_sets.nbp_counts.get_mut(&base_1).expect("Failed to get an element from a hash map.") += 1.;
              }
              for o in n + 1 .. num_of_sa_rows {
                let base_2 = sta.seq_align[o][m];
                let front_base_2 = sta.seq_align[o][m - 1];
                let behind_base_2 = sta.seq_align[o][m + 1];
                let base_pair = get_ordered_base_pair(&(base_1, base_2));
                let base_not_being_gap = base_pair.1;
                if base_1 != GAP && base_2 != GAP {
                  *sta_event_count_sets.ba_counts.get_mut(&base_pair).expect("Failed to get an element from a hash map.") += 1.;
                } else if base_1 != GAP {
                  if front_base_2 == GAP && behind_base_2 == GAP {
                    *sta_event_count_sets.eg_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  } else {
                    *sta_event_count_sets.og_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  }
                } else if base_2 != GAP {
                  if front_base_1 == GAP && behind_base_1 == GAP {
                    *sta_event_count_sets.eg_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  } else {
                    *sta_event_count_sets.og_counts.get_mut(&base_not_being_gap).expect("Failed to get an element from a hash map.") += 1.;
                  }
                }
              }
            }
          }
        },
      }
    }
  }
  let mut stem_params = StemParams::new();
  let mut count_of_all_sta_events = 0.;
  count_of_all_sta_events += sta_event_count_sets.ba_counts.values().fold(0., |acc, &count| {acc + count});
  count_of_all_sta_events += sta_event_count_sets.og_counts.values().fold(0., |acc, &count| {acc + count});
  count_of_all_sta_events += sta_event_count_sets.eg_counts.values().fold(0., |acc, &count| {acc + count});
  count_of_all_sta_events += sta_event_count_sets.bpa_counts.values().fold(0., |acc, &count| {acc + count});
  stem_params.lbaps_with_base_pairs = sta_event_count_sets.ba_counts.iter().map(|(&key, &val)| {(key, (val / count_of_all_sta_events).ln())}).collect();
  stem_params.logps_with_bases = sta_event_count_sets.og_counts.iter().map(|(&key, &val)| {(key, (val / count_of_all_sta_events).ln())}).collect();
  stem_params.legps_with_bases = sta_event_count_sets.eg_counts.iter().map(|(&key, &val)| {(key, (val / count_of_all_sta_events).ln())}).collect();
  stem_params.lbpaps_with_base_quadruples = sta_event_count_sets.bpa_counts.iter().map(|(&key, &val)| {(key, (val / count_of_all_sta_events).ln())}).collect();
  let count_of_all_bps_and_nbps = sta_event_count_sets.bp_counts.values().fold(0., |acc, &count| {acc + count}) + sta_event_count_sets.nbp_counts.values().fold(0., |acc, &count| {acc + count});
  stem_params.lbpps_with_base_pairs = sta_event_count_sets.bp_counts.iter().map(|(&key, &val)| {(key, (val / count_of_all_bps_and_nbps).ln())}).collect();
  stem_params.lnbpps_with_bases = sta_event_count_sets.nbp_counts.iter().map(|(&key, &val)| {(key, (val / count_of_all_bps_and_nbps).ln())}).collect();
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_output_file = String::from("use utils::*;\nlazy_static! {\n  pub static ref STEM_PARAMS: StemParams = {\n    StemParams {\n      lbaps_with_base_pairs: [");
  for (base_pair, &lbap) in &stem_params.lbaps_with_base_pairs {
    buf_4_writer_2_output_file += &format!("(({}, {}), {:e}), ", base_pair.0 as char, base_pair.1 as char, lbap);
    if base_pair.0 != base_pair.1 {
      buf_4_writer_2_output_file += &format!("(({}, {}), {:e}), ", base_pair.1 as char, base_pair.0 as char, lbap);
    }
  }
  buf_4_writer_2_output_file += "].iter().cloned().collect(),\n      logps_with_bases: [";
  for (&base, &logp) in &stem_params.logps_with_bases {
    buf_4_writer_2_output_file += &format!("({}, {:e}), ", base as char, logp);
  }
  buf_4_writer_2_output_file += "].iter().cloned().collect(),\n      legps_with_bases: [";
  for (&base, &legp) in &stem_params.legps_with_bases {
    buf_4_writer_2_output_file += &format!("({}, {:e}), ", base as char, legp);
  }
  buf_4_writer_2_output_file += "].iter().cloned().collect(),\n      lbpaps_with_base_quadruples: [";
  for (base_quadruple, &lbpap) in &stem_params.lbpaps_with_base_quadruples {
    buf_4_writer_2_output_file += &format!("(({}, {}, {}, {}), {:e}), ", base_quadruple.0 as char, base_quadruple.1 as char, base_quadruple.2 as char, base_quadruple.3 as char, lbpap);
    if (base_quadruple.0, base_quadruple.1) != (base_quadruple.2, base_quadruple.3) {
      buf_4_writer_2_output_file += &format!("(({}, {}, {}, {}), {:e}), ", base_quadruple.2 as char, base_quadruple.3 as char, base_quadruple.0 as char, base_quadruple.1 as char, lbpap);
    }
  }
  buf_4_writer_2_output_file += "].iter().cloned().collect(),\n      lbpps_with_base_pairs: [";
  for (base_pair, &lbpp) in &stem_params.lbpps_with_base_pairs {
    buf_4_writer_2_output_file += &format!("(({}, {}), {:e}), ", base_pair.0 as char, base_pair.1 as char, lbpp);
  }
  buf_4_writer_2_output_file += "].iter().cloned().collect(),\n      lnbpps_with_bases: [";
  for (&base, &lnbpp) in &stem_params.lnbpps_with_bases {
    buf_4_writer_2_output_file += &format!("({}, {:e}), ", base as char, lnbpp);
  }
  buf_4_writer_2_output_file += "].iter().cloned().collect(),\n    }\n  };\n}";
  let _ = writer_2_output_file.write_all(buf_4_writer_2_output_file.as_bytes());
}

#[inline]
fn convert_css_str(old_css_str: &str) -> String {
  let mut new_css_str = String::from(old_css_str);
  for (pattern_char, dest_char) in CSS_REPLACEMENT_PATTERN_CHARS.chars().zip(CSS_REPLACEMENT_DEST_CHARS.chars()) {
    new_css_str = new_css_str.replace(pattern_char, &dest_char.to_string());
  }
  new_css_str
}

#[inline]
fn get_css(css_str: &str) -> Css {
  let mut flat_css: HashMap<ColPair, bool, Hasher> = HashMap::default();
  let mut col_stack = Vec::new();
  for (i, css_char) in css_str.chars().enumerate() {
    if css_char == BP_LEFT_BRACKET {
      col_stack.push(i);
    } else if css_char == BP_RIGHT_BRACKET {
      flat_css.insert((col_stack.pop().expect("Failed to pop a vector."), i), true);
    }
  }
  let mut css = Css::default();
  let css_str_len = css_str.len();
  let mut col_pair_stack = vec![(0, css_str_len - 1)];
  while col_pair_stack.len() > 0 {
    let col_pair_1 = col_pair_stack.pop().expect("Failed to pop a vector.");
    let (i, j) = col_pair_1;
    let substr_len_1 = j - i + 1;
    let mut temp_col_pair_stack = vec![];
    for substr_len_2 in (2 .. substr_len_1 - 1).rev() {
      for k in i + 1 .. j + 1 - substr_len_2 {
        let l = k + substr_len_2 - 1;
        let col_pair_2 = (k, l);
        match flat_css.get(&col_pair_2) {
          Some(_) => {
            let mut is_base_pair_accessible = true;
            for &(m, n) in &temp_col_pair_stack {
              if m < k && l < n {
                is_base_pair_accessible = false;
                break;
              }
            }
            if is_base_pair_accessible {
              temp_col_pair_stack.push(col_pair_2);
            }
          }, None => {},
        }
      }
    }
    if temp_col_pair_stack.len() > 0 {
      temp_col_pair_stack.sort_unstable();
      css.insert(col_pair_1, temp_col_pair_stack.clone());
      col_pair_stack.extend(&temp_col_pair_stack);
      temp_col_pair_stack.clear();
    }
  }
  css
}

#[inline]
fn get_ordered_base_pair(base_pair: &BasePair) -> BasePair {
  if base_pair.0 < base_pair.1 {*base_pair} else {(base_pair.1, base_pair.0)}
}

#[inline]
fn get_ordered_base_quadruple(base_quadruple: &BaseQuadruple) -> BaseQuadruple {
  let base_pair_1 = (base_quadruple.0, base_quadruple.1);
  let base_pair_2 = (base_quadruple.2, base_quadruple.3);
  let can_base_pair_bp_1 = BP_ALPHABET.contains_key(&base_pair_1);
  let can_base_pair_bp_2 = BP_ALPHABET.contains_key(&base_pair_2);
  let ordered_base_quadruple = (base_quadruple.2, base_quadruple.3, base_quadruple.0, base_quadruple.1);
  if can_base_pair_bp_1 && can_base_pair_bp_2 {
    if base_pair_1 < base_pair_2 {
      *base_quadruple
    } else {
      ordered_base_quadruple
    }
  } else if can_base_pair_bp_1 {
    *base_quadruple
  } else if can_base_pair_bp_2 {
    ordered_base_quadruple
  } else {
    if base_pair_1 < base_pair_2 {
      *base_quadruple
    } else {
      ordered_base_quadruple
    }
  }
}
