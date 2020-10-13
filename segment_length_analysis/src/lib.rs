use needletail::fastx;
use std::collections::VecDeque;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

pub mod hash_functions;
pub mod qgram_iterator;
use crate::qgram_iterator::{Canonical, HashFunction, QGrams};

#[derive(PartialEq, Debug)]
pub struct Segment {
    pub length: usize,
    pub start: usize,
    pub minimizer: u64,
}

/// Find the rightmost minimum in the active window.
pub fn rightmost_min(window: &VecDeque<u64>) -> (u64, usize) {
    let current_min = *window.iter().min().unwrap();
    for (i, val) in window.iter().enumerate().rev() {
        if *val == current_min {
            return (current_min, i);
        }
    }
    panic!("Could not find a minimum. This is impossible for non-empty windows.");
}

/// Compute segment length from two start positions.
pub fn segment_length(current_start: Option<usize>, last_start: Option<usize>) -> usize {
    match (current_start, last_start) {
        (Some(current_start), Some(last_start)) => current_start - last_start,
        _ => panic!(format!(
            "invalid length computation between {:?} and {:?}",
            current_start, last_start
        )),
    }
}

/// Iterate over the q-gram iterator and compute compressed winnowing segments.
/// Return a vector of segments, each defined by its start position, length (nr of
/// window positions with the same minimizer), and the minimizer value itself.
pub fn segment_length_analysis_by_pos<I>(q_grams: I, w: usize) -> Vec<Segment>
where
    I: Iterator<Item = Option<u64>>,
{
    let mut segments = Vec::with_capacity(100);
    let mut window = VecDeque::with_capacity(w);
    let mut total_q_grams = 0;

    let mut current_min = <u64>::max_value();
    let mut current_min_pos_qgrams = None;
    let mut current_segment_start = None;

    for (q_gram_pos, q_gram) in q_grams.enumerate() {
        total_q_grams += 1; // keep track of this to compute the last window's length

        if let Some(q_gram) = q_gram {
            if window.len() < w {
                // Window is not yet full.
                window.push_back(q_gram);
                if q_gram <= current_min {
                    // prefer rightmost minimum, i.e. use <= instead of < to allow extension of window with the same value
                    current_min = q_gram;
                    current_min_pos_qgrams = Some(q_gram_pos);
                }
                if let None = current_segment_start {
                    // new segment is started, its first entry is the first added q_gram
                    current_segment_start = Some(q_gram_pos);
                }
            } else {
                // Window is full. Update window content.
                let _pushed_out_q_gram = window.pop_front();
                window.push_back(q_gram);

                if current_min_pos_qgrams == Some(q_gram_pos - w) {
                    // The old minimizer has been pushed out. Find a new minimum, create a new segment.
                    // Note: the definitely ends a segment, since same minhash values extend the segment (see below)

                    // Compute the new segment start position i.e.
                    // The old min [#] was pushed out and [*] becomes the new min.
                    // The position of the first q-gram in the window that hits the new min.
                    // Example:
                    // pos: 0    56
                    //      #=====*======
                    //   w: \____/>
                    // The first segments starts at 0. After one shift,
                    // the new minimizer at pos 6 is pushed in.
                    // The new segments starts a position q_gram_pos (6) - w (6) + 1 = 1
                    // Its length is the difference between the new start_positions (1) and the new one (0): 1 - 0 = 1
                    let last_segment_start = current_segment_start;
                    current_segment_start = Some(q_gram_pos - w + 1);
                    segments.push(Segment {
                        length: segment_length(current_segment_start, last_segment_start),
                        start: last_segment_start.unwrap(),
                        minimizer: current_min,
                    });

                    // Search for the minimum of the current window.
                    // If the minimizer occurs several times, use the rightmost.
                    let (minimizer, pos_in_w) = rightmost_min(&window);

                    current_min = minimizer;
                    // Compute q_gram_position of the minimizer from the position in the window
                    //
                    // pos: 0    56
                    //      =====*=======
                    // w:    \____/
                    // q_gram_pos: 6, startpoint of window: q_gram_pos - w + 1
                    // new min pos in q-grams = window_start (6 - 6 + 1 = 1) + position in window (4)
                    // => 1 + 4 = 5
                    // Note that pos_in_w is an index and hence can be at most w-1.
                    current_min_pos_qgrams = Some(q_gram_pos - w + 1 + pos_in_w);
                } else if q_gram == current_min {
                    // Enlongate segment.
                    // The pushed in minimizer is equal to the present minimizer.
                    // No new segment will be opened, only the minimizer position will be updated.
                    current_min_pos_qgrams = Some(q_gram_pos);
                } else if q_gram < current_min {
                    // Smaller minimizer pushed in.
                    // Finish segment, open new one. (See above for example.)
                    let last_segment_start = current_segment_start;
                    current_segment_start = Some(q_gram_pos - w + 1);
                    segments.push(Segment {
                        length: segment_length(current_segment_start, last_segment_start),
                        start: last_segment_start.unwrap(),
                        minimizer: current_min,
                    });

                    current_min = q_gram;
                    current_min_pos_qgrams = Some(q_gram_pos);
                }
            }
        } else {
            // q_gram is None

            if window.len() == w {
                // Window was full. Finalize.
                let last_segment_start = current_segment_start;
                current_segment_start = Some(q_gram_pos - w + 1);
                segments.push(Segment {
                    length: segment_length(current_segment_start, last_segment_start),
                    start: last_segment_start.unwrap(),
                    minimizer: current_min,
                });
            }

            // Either way, drop the remainder of the window.
            window.clear();
            current_min = <u64>::max_value();
            current_segment_start = None;
            current_min_pos_qgrams = None;
        }
    }

    // Handle the last window, if a window is still open but the q-grams ran out.
    if window.len() == w {
        let last_segment_start = current_segment_start;
        current_segment_start = Some(total_q_grams - w + 1);
        segments.push(Segment {
            length: segment_length(current_segment_start, last_segment_start),
            start: last_segment_start.unwrap(),
            minimizer: current_min,
        });
    }

    segments
}


/// Compute empirical distribution of segment lengths.
/// Open the given FASTA file, read in its content, winnow it
/// and write the counts (nr of times a certain segment length
/// wasd observed) to the output file.
/// Additionally a file containing the minimizers themselves and
/// how often they were observed is written.
pub fn count_segment_lengths(
    fasta_path: &str,
    seg_len_path: &str,
    minimizer_path: &str,
    q: usize,
    w: usize,
    canonical: Canonical,
    hash_function: HashFunction,
    canonisize_qgrams: bool,
    write_minimizers: bool,
) {
    let mut counts = counter::Counter::new();
    let mut minimizers = counter::Counter::new();
    fastx::fastx_file(&fasta_path[..], |seq| {

        if seq.seq.len() > q {
            let q_grams = QGrams::new(&seq.seq, q, canonical.clone(), hash_function.clone(), canonisize_qgrams);

            let segments = segment_length_analysis_by_pos(q_grams, w);
            counts += segments.iter().map(|x| x.length);
            if write_minimizers {
                minimizers += segments.iter().map(|x| x.minimizer);
            };
        }
    })
    .expect("failed to iterate through FASTA file.");

    // Handle output for segment lengths
    let mut file = match File::create(seg_len_path) {
        Ok(file) => {
            eprintln!("Writing output to {}", seg_len_path);
            BufWriter::new(file)
        }
        Err(e) => panic!("Could not open file: {}", e),
    };

    let mut counts: Vec<(&usize, &usize)> = counts.iter().collect();
    counts.sort();

    for (key, count) in counts {
        match file.write(format!("{} {}\n", key, count).as_bytes()) {
            Ok(_) => (),
            Err(e) => panic!("Writing line to file failed: {}", e),
        }
    }
    
    // Handle output for minhash tally file
    let mut distr_file = match File::create(minimizer_path) {
        Ok(distr_file) => {
            eprintln!("Writing minimizer counts to {}", minimizer_path);
            BufWriter::new(distr_file)
        }
        Err(e) => panic!("Could not open output file for distribution: {}", e),
    };

    if write_minimizers {
        let mut minimizers: Vec<(&u64, &usize)> = minimizers.iter().collect();
        minimizers.sort();
    
        for (key, count) in minimizers {
            match distr_file.write(format!("{} {}\n", key, count).as_bytes()) {
                Ok(_) => (),
                Err(e) => panic!("Writing line to file failed: {}", e),
            }
        }
    }
}

#[cfg(test)]
mod tests;
mod unwieldy_tests;
