use mod_exp;
use rand;
use rand::Rng;
use serde::{Deserialize, Serialize};
use murmur3;
use std::mem;

#[derive(Clone, Serialize, Deserialize)]
pub struct InvMultParams {
    pub a: u64,
}


impl InvMultParams {
    pub fn new() -> Self {
        // restrict the random value to be uneven and "large" so that the hash function scatters well
        let mut rng = rand::thread_rng();
        let mut mult = rng.gen_range(1_000_000, 18446744073709551615); // 2^{64}-1
        if mult % 2 == 0 {
            mult += 1;
        }
        InvMultParams { a: mult }
    }

    pub fn with_params(a: u64) -> Self {
        InvMultParams { a }
    }
}


#[derive(Clone, Serialize, Deserialize)]
pub struct HlinParams {
    a64: u64,
    b64: u64,
    a128: u128,
    b128: u128,
}

impl HlinParams {
    pub fn new() -> Self {
        HlinParams {
            a64: rand::random::<u64>(),
            b64: rand::random::<u64>(),
            a128: rand::random::<u128>(),
            b128: rand::random::<u128>(),
        }
    }

    pub fn with_params(a64: u64, b64: u64, a128: u128, b128: u128) -> Self {
        HlinParams {
            a64,
            b64,
            a128,
            b128,
        }
    }
}


/// Dietzfelbinger's 2-independant H^{lin} hash function for 32 bits.
pub fn df_32(x: u32, params: &HlinParams) -> u32 {
    ((params.a64.wrapping_mul(x as u64).wrapping_add(params.b64)) >> 32) as u32
}

/// Same as above, but adapted for 64 bit values
pub fn df_64(x: u64, params: &HlinParams) -> u64 {
    ((params
        .a128
        .wrapping_mul(x as u128)
        .wrapping_add(params.b128))
        >> 64) as u64
}

/// Hash mixing from A. Appleby's MurmurHash3
/// https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
pub fn fmix_64(y: u64) -> u64 {
    let mut x = y;
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^= x >> 33;

    x
}

const LOW: u64 = 0b_00000000_00000000_00000000_00000000_11111111_11111111_11111111_11111111;
const HIGH: u64 = 0b_11111111_11111111_11111111_11111111_00000000_00000000_00000000_00000000;

/// Swap low and high words of a u64 integer.
pub fn swap_words(hash: u64) -> u64 {
    let low = hash & LOW;
    let high = hash & HIGH;
    (low << 32) | (high >> 32)
}

/// Swap low and high bits of a u64 integer.
pub fn swap_words_q(hash: u64, q: u64) -> u64 {
    let (low_bits, high_bits) = if q == 32 {
        (LOW, HIGH)
    } else {
        (
            2_u64.pow(q as u32) - 1,
            (2_u64.pow(2 * q as u32) - 1) ^ (2_u64.pow(q as u32) - 1)
        )
    };
    let low = hash & low_bits;
    let high = hash & high_bits;
    (low << q) | (high >> q)
}


/// Invertible multiplicative hash function H^{\text{im}} (cf. Section 4.4)
/// Compute an (address, fingerprint) tuple from a u64 key
///
/// m: (uneven) multiplier (for universal hashing)
/// n: hash table size (# slots)
/// 
/// NOTE: size of key universe u needs to be a power of 2
/// m needs to be uneven
pub fn inv_mult_hash(key: u64, m: u64, n: u64) -> (u64, u64) {
    let hv = key.wrapping_mul(m);
    let p = hv % n;
    let f = hv / n;
    (p, f)
}


pub fn simplified_inv_mult_hash_32(key: u32, a: u32) -> u32 {
    key.wrapping_mul(a)
}

pub fn simplified_inv_mult_hash_64(key: u64, a: u64) -> u64 {
    key.wrapping_mul(a)
}

const U64_BITMASK: u128 = 0b_11111111_11111111_11111111_11111111_11111111_11111111_11111111_11111111;

/// 64-bit Murmurhash3 using its Rust implementation of Stuart Small
/// https://crates.io/crates/murmur3
pub fn mmh3_64(key: u64, seed: u32) -> u64 {
    let key = unsafe {
        mem::transmute::<u64, [u8; 8]>(key)
    };
    (murmur3::murmur3_x64_128(&mut key.as_ref(), seed).unwrap() & U64_BITMASK) as u64
}



#[derive(Clone, Serialize, Deserialize)]
pub struct Mmh3Params {
    pub seed: u32,
}

impl Mmh3Params {
    pub fn new() -> Self {
        let mut rng = rand::thread_rng();
        let seed = rng.gen_range(0, 4294967295); // 2^{32}-1
        Mmh3Params { seed }
    }

    pub fn with_params(seed: u32) -> Self {
        Mmh3Params { seed }
    }
}
