use super::hash_functions::HlinParams;
use super::qgram_iterator::{Canonical, HashFunction, QGrams};

use super::*;


#[test]
fn test_qgrams_iupac() {
    let seq = b"TGACTNTGACT";
    /*
    TGACTNTGACT
    forward:    1110000111XX1110000111 -> 111000, 100001, 000111, None, 111000, 100001, 000111
    rev comp:                          -> 110100, 101101, 001011, None, 110100, 101101, 001011
    min canonical:                        110100, 100001, 000111, None, 110100, 100001, 000111
    max canonical:                        111000, 101101, 001011, None, 111000, 101101, 001011
    */
    let qgrams = QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(
        qgrams,
        vec![
            Some(0b_111000),
            Some(0b_100001),
            Some(0b_000111),
            None,
            Some(0b_111000),
            Some(0b_100001),
            Some(0b_000111)
        ]
    );

    let qgrams =
        QGrams::new(seq, 3, Canonical::Min, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(
        qgrams,
        vec![
            Some(0b_110100),
            Some(0b_100001),
            Some(0b_000111),
            None,
            Some(0b_110100),
            Some(0b_100001),
            Some(0b_000111)
        ]
    );

    let qgrams =
        QGrams::new(seq, 3, Canonical::Max, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(
        qgrams,
        vec![
            Some(0b_111000),
            Some(0b_101101),
            Some(0b_001011),
            None,
            Some(0b_111000),
            Some(0b_101101),
            Some(0b_001011)
        ]
    );
}

#[test]
fn test_edgecases() {
    // these should panic due to q < seq.len()
    let seq = b"";
    let qgrams = std::panic::catch_unwind(|| {
        QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>()
    });
    assert!(qgrams.is_err());

    let seq = b"ACGT";
    let qgrams = std::panic::catch_unwind(|| {
        QGrams::new(seq, 5, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>()
    });
    assert!(qgrams.is_err());

    // This should onyl return None
    let seq = b"NNNNN";
    let qgrams = QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![None]);
}

#[test]
fn test_leading_iupac() {
    let seq = b"NTGA";
    let qgrams = QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![None, Some(0b_111000)]);
}

#[test]
fn test_trailing_iupac() {
    let seq = b"TGAN";
    let qgrams = QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![Some(0b_111000), None]);
}

#[test]
fn test_infix_iupac() {
    let seq = b"NAANAAN";
    let qgrams = QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![None]);

    let seq = b"NAANAANTTTNAAN";
    let qgrams = QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![None, Some(0b_111111), None]);

    let seq = b"TTTT!!!TTTT";
    let qgrams = QGrams::new(seq, 4, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![Some(0b_11111111), None, Some(0b_11111111)]);
}

#[test]
fn test_capitalization() {
    let seq_lower_case = b"acgttgca";
    let qgrams_lower_case = QGrams::new(seq_lower_case, 3, Canonical::No, HashFunction::No, false)
        .collect::<Vec<Option<u64>>>();
    let seq_upper_case = b"ACGTTGCA";
    let qgrams_upper_case = QGrams::new(seq_upper_case, 3, Canonical::No, HashFunction::No, false)
        .collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams_lower_case, qgrams_upper_case);
}

#[test]
fn test_qgram_generation_small() {
    let seq = b"ACGT";
    /*
    ACGT
    forward:    00011011
    rev comp:   00011011
    canonical:  00011011
     */
    let qgrams = QGrams::new(seq, 4, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![Some(0b_00_01_10_11)],);

    let qgrams =
        QGrams::new(seq, 4, Canonical::Min, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![Some(0b_00_01_10_11)],);

    let qgrams =
        QGrams::new(seq, 4, Canonical::Max, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(qgrams, vec![Some(0b_00_01_10_11)],);

    let seq = b"TGACT";
    /*
    TGACT  (rc: AGTCA)
    forward:    1110000111 -> 111000, 100001, 000111
    rev comp:              -> 110100, 101101, 001011
    canonical:                110100, 100001, 000111
     */
    let qgrams = QGrams::new(seq, 3, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(
        qgrams,
        vec![Some(0b_111000), Some(0b_100001), Some(0b_000111)],
    );

    let qgrams =
        QGrams::new(seq, 3, Canonical::Min, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(
        qgrams,
        vec![Some(0b_110100), Some(0b_100001), Some(0b_000111)],
    );

    let qgrams =
        QGrams::new(seq, 3, Canonical::Max, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(
        qgrams,
        vec![Some(0b_111000), Some(0b_101101), Some(0b_001011)],
    );
}

#[test]
fn test_qgrams_small_all_variants() {
    let seq = b"GACGTTTGGA";
    // 10_00_01_10_11_11_11_10_10_00
    // length 9

    let expected_no_no = vec![
        Some(0b_10_00_01_10_11_11_11),
        Some(0b_00_01_10_11_11_11_10),
        Some(0b_01_10_11_11_11_10_10),
        Some(0b_10_11_11_11_10_10_00),
    ];
    let qgrams = QGrams::new(seq, 7, Canonical::No, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_no_no, qgrams);

    let expected_min_no = vec![
        Some(0b_00_00_00_01_10_11_01),
        Some(0b_00_01_10_11_11_11_10),
        Some(0b_01_01_00_00_00_01_10),
        Some(0b_10_11_11_11_10_10_00),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::Min, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_min_no, qgrams);

    let expected_max_no = vec![
        Some(0b_10_00_01_10_11_11_11),
        Some(0b_01_00_00_00_01_10_11),
        Some(0b_01_10_11_11_11_10_10),
        Some(0b_11_01_01_00_00_00_01),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::Max, HashFunction::No, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_max_no, qgrams);

    let expected_no_swap = vec![
        Some(0b_0111111_1000011),
        Some(0b_1111110_0001101),
        Some(0b_1111010_0110111),
        Some(0b_1101000_1011111),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::No, HashFunction::WordSwap, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_no_swap, qgrams);

    let expected_min_swap = vec![
        Some(0b_0111111_1000011),
        Some(0b_0011011_0100000),
        Some(0b_0000110_0101000),
        Some(0b_0000001_1101010),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::Min, HashFunction::WordSwap, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_min_swap, qgrams);

    let expected_max_swap = vec![
        Some(0b_1101101_0000000),
        Some(0b_1111110_0001101),
        Some(0b_1111010_0110111),
        Some(0b_1101000_1011111),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::Max, HashFunction::WordSwap, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_max_swap, qgrams);

    // params randomly chosen
    let params = HlinParams::with_params(
        42,
        42, // unused here
        129809215886441603082368137654056488438_u128,
        249051908688722022154549480275671269000_u128,
    );
    let expected_no_df = vec![
        Some(5431125397752018329),
        Some(10554260287900737422),
        Some(15787591469578849677),
        Some(2606460887913369201),
    ];
    let qgrams = QGrams::new(seq, 7, Canonical::No, HashFunction::Hlin(params.clone()), false)
        .collect::<Vec<Option<u64>>>();
    assert_eq!(expected_no_df, qgrams);

    let expected_min_df = vec![
        Some(5431125397752018329),
        Some(10206995939169286940),
        Some(3186624425311059255),
        Some(2606460887913369201),
    ];
    let qgrams = QGrams::new(seq, 7, Canonical::Min, HashFunction::Hlin(params.clone()), false)
        .collect::<Vec<Option<u64>>>();
    assert_eq!(expected_min_df, qgrams);

    let expected_max_df = vec![
        Some(5767794533620661932),
        Some(10554260287900737422),
        Some(15787591469578849677),
        Some(17822716449235860937),
    ];
    let qgrams = QGrams::new(seq, 7, Canonical::Max, HashFunction::Hlin(params.clone()), false)
        .collect::<Vec<Option<u64>>>();
    assert_eq!(expected_max_df, qgrams);

    let expected_no_fmix = vec![
        Some(7526011328732148264),
        Some(3162698631890770771),
        Some(4195643915253645620),
        Some(453378680581374969),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::No, HashFunction::Fmix, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_no_fmix, qgrams);

    let expected_min_fmix = vec![
        Some(7526011328732148264),
        Some(3162698631890770771),
        Some(4195643915253645620),
        Some(453378680581374969),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::Min, HashFunction::Fmix, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_min_fmix, qgrams);

    let expected_max_fmix = vec![
        Some(16996655880070524895),
        Some(13352430448427503177),
        Some(11095935145862035941),
        Some(13051068276403513641),
    ];
    let qgrams =
        QGrams::new(seq, 7, Canonical::Max, HashFunction::Fmix, false).collect::<Vec<Option<u64>>>();
    assert_eq!(expected_max_fmix, qgrams);
}

#[test]
fn test_rightmost_min() {
    let window = std::collections::VecDeque::from(vec![0, 1, 2, 3, 4, 5]);
    assert_eq!(rightmost_min(&window), (0, 0));

    let window = std::collections::VecDeque::from(vec![5, 4, 3, 2, 1, 0]);
    assert_eq!(rightmost_min(&window), (0, 5));

    let window = std::collections::VecDeque::from(vec![2, 1, 3, 0, 4, 5]);
    assert_eq!(rightmost_min(&window), (0, 3));

    let window = std::collections::VecDeque::from(vec![0, 0, 0, 0, 0, 0]);
    assert_eq!(rightmost_min(&window), (0, 5));
}

#[test]
fn test_windows_small() {
    let q_grams = vec![
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_01),
        Some(0b_00_00_00_00_00_01_00),
        Some(0b_00_00_00_00_01_00_00),
        Some(0b_00_00_00_01_00_00_00),
        Some(0b_00_00_01_00_00_00_00),
        Some(0b_00_01_00_00_00_00_00),
    ];
    // 0, 0, 1, 4, 16, 64, 256, 1024
    // \____________/
    //       \_______________/
    //              \_______________/
    let lengths = segment_length_analysis_by_pos(q_grams.iter().cloned(), 5);
    let expected_lengths = vec![
        Segment {
            length: 2,
            start: 0,
            minimizer: 0,
        }, // first window (pos 0-5)
        Segment {
            length: 1,
            start: 2,
            minimizer: 0b_01,
        }, // second window (pos 2-6)
        Segment {
            length: 1,
            start: 3,
            minimizer: 0b_01_00,
        }, // third window (pos 3-7)
    ];
    assert_eq!(lengths, expected_lengths);

    let q_grams = vec![
        Some(0b_11_11_11_11_11_11_11),
        Some(0b_11_11_11_11_11_11_01),
        Some(0b_11_11_11_11_11_01_11),
        Some(0b_11_11_11_11_01_11_11),
        Some(0b_11_11_11_01_11_11_11),
        Some(0b_11_11_01_11_11_11_11),
        Some(0b_11_01_11_11_11_11_11),
        Some(0b_01_11_11_11_11_11_11),
    ];
    // 16383, 16381, 16375, 16351, 16255, 15871, 14335, 8191
    // \________________________________/
    //        \________________________________/
    //               \________________________________/
    //                     \________________________________/
    let lengths = segment_length_analysis_by_pos(q_grams.iter().cloned(), 5);
    let expected_lengths = vec![
        Segment {
            length: 1,
            start: 0,
            minimizer: 16255,
        },
        Segment {
            length: 1,
            start: 1,
            minimizer: 15871,
        },
        Segment {
            length: 1,
            start: 2,
            minimizer: 14335,
        },
        Segment {
            length: 1,
            start: 3,
            minimizer: 8191,
        },
    ];
    assert_eq!(lengths, expected_lengths);
}

#[test]
fn test_windows_with_subsequences() {
    let q_grams = vec![
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_01),
        Some(0b_00_00_00_00_00_01_00),
        Some(0b_00_00_00_00_01_00_00),
        Some(0b_00_00_00_01_00_00_00),
        Some(0b_00_00_01_00_00_00_00),
        Some(0b_00_01_00_00_00_00_00),
        None,
        Some(0b_11_11_11_11_11_11_11),
        Some(0b_11_11_11_11_11_11_01),
        Some(0b_11_11_11_11_11_01_11),
        Some(0b_11_11_11_11_01_11_11),
        Some(0b_11_11_11_01_11_11_11),
        Some(0b_11_11_01_11_11_11_11),
        Some(0b_11_01_11_11_11_11_11),
        Some(0b_01_11_11_11_11_11_11),
    ];
    let expected_lengths = vec![
        Segment {
            length: 2,
            start: 0,
            minimizer: 0b_0,
        },
        Segment {
            length: 1,
            start: 2,
            minimizer: 0b_1,
        },
        Segment {
            length: 1,
            start: 3,
            minimizer: 0b_01_00,
        },
        // this is where the None hits
        Segment {
            length: 1,
            start: 9,
            minimizer: 0b_11_11_11_01_11_11_11,
        },
        Segment {
            length: 1,
            start: 10,
            minimizer: 0b_11_11_01_11_11_11_11,
        },
        Segment {
            length: 1,
            start: 11,
            minimizer: 0b_11_01_11_11_11_11_11,
        },
        Segment {
            length: 1,
            start: 12,
            minimizer: 0b_01_11_11_11_11_11_11,
        },
    ];
    let lengths = segment_length_analysis_by_pos(q_grams.iter().cloned(), 5);
    assert_eq!(lengths, expected_lengths);
}

#[test]
fn test_subsequence_edge_cases() {
    // None at the front. Should be skipped.
    let q_grams = vec![
        None,
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_01),
        Some(0b_00_00_00_00_00_01_00),
        Some(0b_00_00_00_00_01_00_00),
        Some(0b_00_00_00_01_00_00_00),
        Some(0b_00_00_01_00_00_00_00),
        Some(0b_00_01_00_00_00_00_00),
    ];
    let expected_lengths = vec![
        Segment {
            length: 2,
            start: 1,
            minimizer: 0,
        }, // first window (pos 0-5)
        Segment {
            length: 1,
            start: 3,
            minimizer: 0b_01,
        }, // second window (pos 2-6)
        Segment {
            length: 1,
            start: 4,
            minimizer: 0b_01_00,
        }, // third window (pos 3-7)
    ];
    let lengths = segment_length_analysis_by_pos(q_grams.iter().cloned(), 5);
    assert_eq!(lengths, expected_lengths);

    // None in the back. Should terminate early.
    let q_grams = vec![
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_01),
        Some(0b_00_00_00_00_00_01_00),
        Some(0b_00_00_00_00_01_00_00),
        Some(0b_00_00_00_01_00_00_00),
        Some(0b_00_00_01_00_00_00_00),
        Some(0b_00_01_00_00_00_00_00),
        None,
    ];
    let expected_lengths = vec![
        Segment {
            length: 2,
            start: 0,
            minimizer: 0,
        }, // first window (pos 0-5)
        Segment {
            length: 1,
            start: 2,
            minimizer: 0b_01,
        }, // second window (pos 2-6)
        Segment {
            length: 1,
            start: 3,
            minimizer: 0b_01_00,
        }, // third window (pos 3-7)
    ];
    let lengths = segment_length_analysis_by_pos(q_grams.iter().cloned(), 5);
    assert_eq!(lengths, expected_lengths);

    // One None to rule them all (out)
    let q_grams = vec![
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_00),
        Some(0b_00_00_00_00_00_00_01),
        Some(0b_00_00_00_00_00_01_00),
        None,
        Some(0b_00_00_00_00_01_00_00),
        Some(0b_00_00_00_01_00_00_00),
        Some(0b_00_00_01_00_00_00_00),
        Some(0b_00_01_00_00_00_00_00),
    ];
    let expected_lengths: Vec<Segment> = vec![];
    let lengths = segment_length_analysis_by_pos(q_grams.iter().cloned(), 5);
    assert_eq!(lengths, expected_lengths);
}

#[test]
fn test_windows_big() {
    /*
    import random
    x = [random.randint(0,1000) for _ in range(100)]

    lengths = []
    window_length = 0
    window_start = 0
    last_min = None

    for i in range(len(x) - 5 + 1):
        m = min(x[i:i+5])
        if last_min is None:
            last_min = m
            window_length = 1
        elif m == last_min:
            window_length += 1
        else:
            lengths.append((window_length, window_start, last_min))
            last_min = m
            window_length = 1
            window_start = i

    This misses the last window!!
    print(lengths)


    for (l, s, m) in lengths:
        print(f"Segment{{\n   length: {l},\n    start: {s},\n}},")
    */

    // w = 5
    let q_grams = vec![
        273, 508, 360, 985, 320, 733, 470, 852, 158, 63, 23, 381, 692, 498, 143, 739, 811, 840,
        286, 209, 917, 328, 748, 663, 570, 410, 689, 148, 848, 170, 523, 323, 589, 62, 305, 302,
        937, 153, 332, 114, 250, 214, 923, 508, 701, 290, 315, 294, 491, 364, 231, 929, 917, 621,
        861, 282, 278, 708, 660, 19, 36, 127, 213, 403, 202, 583, 808, 522, 903, 18, 537, 890, 937,
        263, 886, 585, 564, 898, 748, 391, 179, 872, 330, 351, 127, 999, 988, 567, 936, 513, 253,
        756, 625, 20, 366, 558, 876, 313, 293, 773,
    ];
    let expected_lengths = vec![
        Segment {
            length: 1,
            start: 0,
            minimizer: 273,
        },
        Segment {
            length: 3,
            start: 1,
            minimizer: 320,
        },
        Segment {
            length: 1,
            start: 4,
            minimizer: 158,
        },
        Segment {
            length: 1,
            start: 5,
            minimizer: 63,
        },
        Segment {
            length: 5,
            start: 6,
            minimizer: 23,
        },
        Segment {
            length: 4,
            start: 11,
            minimizer: 143,
        },
        Segment {
            length: 5,
            start: 15,
            minimizer: 209,
        },
        Segment {
            length: 2,
            start: 20,
            minimizer: 328,
        },
        Segment {
            length: 1,
            start: 22,
            minimizer: 410,
        },
        Segment {
            length: 5,
            start: 23,
            minimizer: 148,
        },
        Segment {
            length: 1,
            start: 28,
            minimizer: 170,
        },
        Segment {
            length: 5,
            start: 29,
            minimizer: 62,
        },
        Segment {
            length: 1,
            start: 34,
            minimizer: 153,
        },
        Segment {
            length: 5,
            start: 35,
            minimizer: 114,
        },
        Segment {
            length: 2,
            start: 40,
            minimizer: 214,
        },
        Segment {
            length: 4,
            start: 42,
            minimizer: 290,
        },
        Segment {
            length: 5,
            start: 46,
            minimizer: 231,
        },
        Segment {
            length: 1,
            start: 51,
            minimizer: 282,
        },
        Segment {
            length: 3,
            start: 52,
            minimizer: 278,
        },
        Segment {
            length: 5,
            start: 55,
            minimizer: 19,
        },
        Segment {
            length: 1,
            start: 60,
            minimizer: 36,
        },
        Segment {
            length: 1,
            start: 61,
            minimizer: 127,
        },
        Segment {
            length: 3,
            start: 62,
            minimizer: 202,
        },
        Segment {
            length: 5,
            start: 65,
            minimizer: 18,
        },
        Segment {
            length: 4,
            start: 70,
            minimizer: 263,
        },
        Segment {
            length: 1,
            start: 74,
            minimizer: 564,
        },
        Segment {
            length: 1,
            start: 75,
            minimizer: 391,
        },
        Segment {
            length: 4,
            start: 76,
            minimizer: 179,
        },
        Segment {
            length: 5,
            start: 80,
            minimizer: 127,
        },
        Segment {
            length: 1,
            start: 85,
            minimizer: 513,
        },
        Segment {
            length: 3,
            start: 86,
            minimizer: 253,
        },
        Segment {
            length: 5,
            start: 89,
            minimizer: 20,
        },
        Segment {
            length: 2,
            start: 94,
            minimizer: 293,
        },
    ];

    let lengths = segment_length_analysis_by_pos(q_grams.iter().cloned().map(|x| Some(x)), 5);
    assert_eq!(lengths, expected_lengths);
}

#[test]
fn test_swap_even_q() {
    let qgrams = vec![
        0b_00000000_00000000_00000000_00000000_11111111_11111111_11111111_11111111,
        0b_11111111_11111111_11111111_11111111_00000000_00000000_00000000_00000000,
        0b_11111111_00000000_11111111_00000000_11111111_00000000_11111111_00000000,
        0b_10101010_10101010_10101010_10101010_01010101_01010101_01010101_01010101,
        0b_10000110_00001000_10000011_10110010_00000011_01110000_01101110_01000111,
    ];

    let expected_swapped = vec![
        0b_11111111_11111111_11111111_11111111_00000000_00000000_00000000_00000000,
        0b_00000000_00000000_00000000_00000000_11111111_11111111_11111111_11111111,
        0b_11111111_00000000_11111111_00000000_11111111_00000000_11111111_00000000,
        0b_01010101_01010101_01010101_01010101_10101010_10101010_10101010_10101010,
        0b_00000011_01110000_01101110_01000111_10000110_00001000_10000011_10110010,
    ];

    for (qgram, expected) in qgrams.iter().zip(expected_swapped.iter()) {
        assert_eq!(hash_functions::swap_words_q(*qgram, 32), *expected);
    }
}

#[test]
fn test_swap_uneven_q() {
    let qgrams = vec![
        0b_0000000000000000000000000000000_0000000000000000000000000000000,
        0b_1111111111111111111111111111111_1111111111111111111111111111111,
        0b_1000000000000000000000000000001_1000000000000000000000000000001,
        0b_0000000000000000000000000000001_1111111111111111111111111111111,
        0b_1111111111111111111111111111110_0000000000000000000000000000000,
        0b_1111110000000011111111000000001_1111111000000001111111100000000,
        0b_1111111100000000111111110000000_0111111110000000011111111000000,
        0b_1010101010101010101010101010100_1010101010101010101010101010101,
        0b_1010101010101010101010101010101_0010101010101010101010101010101,
        0b_0001100000100010000011101100100_0000011011100000110111001000111,
        0b_1000011000001000100000111011001_0000000110111000001101110010001,
    ];

    let expected_swapped = vec![
        0b_0000000000000000000000000000000_0000000000000000000000000000000,
        0b_1111111111111111111111111111111_1111111111111111111111111111111,
        0b_1000000000000000000000000000001_1000000000000000000000000000001,
        0b_1111111111111111111111111111111_0000000000000000000000000000001,
        0b_0000000000000000000000000000000_1111111111111111111111111111110,
        0b_1111111000000001111111100000000_1111110000000011111111000000001,
        0b_0111111110000000011111111000000_1111111100000000111111110000000,
        0b_1010101010101010101010101010101_1010101010101010101010101010100,
        0b_0010101010101010101010101010101_1010101010101010101010101010101,
        0b_0000011011100000110111001000111_0001100000100010000011101100100,
        0b_0000000110111000001101110010001_1000011000001000100000111011001,
    ];

    for (qgram, expected) in qgrams.iter().zip(expected_swapped.iter()) {
        assert_eq!(hash_functions::swap_words_q(*qgram, 31), *expected);
    }

    let qgrams = vec![
        0b_00011100111101000_00110000000111011,
        0b_01000000110101100_10001100101110110,
        0b_00101000010110110_00111011000101001,
        0b_01011111110111111_00000000010101111,
        0b_01010101100011110_10000001010101010,
    ];

    let expected_swapped = vec![
        0b_00110000000111011_00011100111101000,
        0b_10001100101110110_01000000110101100,
        0b_00111011000101001_00101000010110110,
        0b_00000000010101111_01011111110111111,
        0b_10000001010101010_01010101100011110,
    ];

    for (qgram, expected) in qgrams.iter().zip(expected_swapped.iter()) {
        assert_eq!(hash_functions::swap_words_q(*qgram, 17), *expected);
    }
}
