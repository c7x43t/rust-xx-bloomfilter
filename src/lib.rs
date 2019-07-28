//! Bloom filter for Rust - forked from https://github.com/jedisct1/rust-bloom-filter
//!
//! This is a simple but fast Bloom filter implementation, that requires only
//! 2 hash functions, generated with XXHash64 using randomized keys.

extern crate bit_vec;
extern crate rand;
extern crate twox_hash;

use bit_vec::BitVec;
use std::cmp;
use std::f64;
use std::hash::{Hash, Hasher};
use twox_hash::XxHash64;

/// Bloom filter structure
pub struct Bloom {
    bitmap: BitVec,
    bitmap_size: u64,
    k: u64,
    xx: (XxHash64, XxHash64),
}

impl Bloom {
    /// Create a new bloom filter structure.
    /// bitmap_size is the size in bytes (not bits) that will be allocated in memory
    /// items_count is an estimation of the maximum number of items to store.
    pub fn new(bitmap_size: usize, items_count: usize) -> Self {
        assert!(bitmap_size > 0 && items_count > 0);
        let bitmap_size = (bitmap_size as u64) * 8u64;
        let k = Self::optimal_k_num(bitmap_size, items_count);
        let bitmap = BitVec::from_elem(bitmap_size as usize, false);
        let xx = (Self::xx_new(), Self::xx_new());
        Self {
            bitmap,
            bitmap_size,
            k,
            xx,
        }
    }

    /// Create a new bloom filter structure.
    /// items_count is an estimation of the maximum number of items to store.
    /// fp_p is the wanted rate of false positives, in ]0.0, 1.0[
    /// ```
    /// extern crate xx_bloomfilter;
    /// extern crate rand;
    ///
    /// use xx_bloomfilter::Bloom;
    ///
    /// let mut bloom = Bloom::new_with_rate(1_000_000, 1e-6);
    /// let item: u64 = rand::random();
    /// assert_eq!(false, bloom.check_and_add(&item));
    /// assert_eq!(true, bloom.check(&item));
    /// // Clear all values
    /// bloom.clear();
    /// assert_eq!(false, bloom.check_and_add(&item));
    /// ```
    pub fn new_with_rate(items_count: usize, fp_p: f64) -> Self {
        let bitmap_size = Self::compute_bitmap_size(items_count, fp_p);
        Bloom::new(bitmap_size, items_count)
    }

    pub fn from_existing_struct(other: &Bloom) -> Self {
//        let xx = [
//            XxHash64::with_seed(other.xx()[0].seed),
//            XxHash64::with_seed(other.xx()[1].seed)
//        ];
        Self {
            bitmap: BitVec::from_bytes(other.bitmap().to_bytes().as_slice()),
            bitmap_size: other.bitmap_size,
            k: other.k,
            xx: other.xx(),
        }
    }


    /// Create a bloom filter structure with an existing state.
    /// The state is assumed to be retrieved from an existing bloom filter.
    pub fn from_existing(
        bitmap: &[u8],
        bitmap_size: u64,
        k: u64,
        xx: (XxHash64, XxHash64),
    ) -> Self {
        Self {
            bitmap: BitVec::from_bytes(bitmap),
            bitmap_size,
            k,
            xx,
        }
    }

    /// Compute a recommended bitmap size for items_count items
    /// and a fp_p rate of false positives.
    /// fp_p has to be within the ]0.0, 1.0[ range.
    pub fn compute_bitmap_size(items_count: usize, fp_p: f64) -> usize {
        assert!(items_count > 0);
        assert!(fp_p > 0.0 && fp_p < 1.0);
        let log2 = f64::consts::LN_2;
        let log2_2 = log2 * log2;
        ((items_count as f64) * f64::ln(fp_p) / (-8.0 * log2_2)).ceil() as usize
    }

    fn bit_offset(&self, hashes: (u64, u64), i_k: u64) -> usize {
        (self.double_hash(hashes, i_k) % self.bitmap_size) as usize
    }

    /// Record the presence of an item.
    pub fn add<T: Hash>(&mut self, item: &T) {
        let hashes = self.hashes(item);
//        let offsets = (0..self.k).map(|k| );
        for i_k in 0..self.k {
            let bit_offset = self.bit_offset(hashes, i_k);
            self.bitmap.set(bit_offset, true);
        }
    }

    /// Check if an item is present in the set.
    /// There can be false positives, but no false negatives.
    pub fn check<T: Hash>(&self, item: &T) -> bool {
        let hashes = self.hashes(item);
        for i_k in 0..self.k {
            let bit_offset = self.bit_offset(hashes, i_k);
            if !self.bitmap.get(bit_offset).unwrap_or_else(|| panic!("bit_offset {} not in bitmap!", bit_offset)) {
                return false;
            }
        }
        true
    }

    /// Record the presence of an item in the set,
    /// and return the previous state of this item.
    pub fn check_and_add<T: Hash>(&mut self, item: &T) -> bool {
        let hashes = self.hashes(item);
        let mut found = true;
        for i_k in 0..self.k {
            let bit_offset = self.bit_offset(hashes, i_k);
            if !self.bitmap.get(bit_offset).unwrap_or_else(|| panic!("bit_offset {} not in bitmap!", bit_offset)) {
                found = false;
                self.bitmap.set(bit_offset, true);
            }
        }
        found
    }

    /// Return the bitmap
    pub fn bitmap(&self) -> BitVec {
        self.bitmap.clone()
    }

    /// Return the number of bits in the filter
    pub fn number_of_bits(&self) -> u64 {
        self.bitmap_size
    }

    /// Return the number of hash functions used for `check` and `set`
    pub fn number_of_hash_functions(&self) -> u64 {
        self.k
    }

    pub fn xx(&self) -> (XxHash64, XxHash64) {
        self.xx
    }

    fn optimal_k_num(bitmap_size: u64, items_count: usize) -> u64 {
        let m = bitmap_size as f64;
        let n = items_count as f64;
        let k = (m / n * f64::ln(2.0f64)).ceil() as u64;
        cmp::max(k, 1)
    }

    fn hash1<T: Hash>(&self, t: &T) -> u64 {
        let mut hasher = self.xx().0;
        t.hash(&mut hasher);
        hasher.finish()
    }

    fn hash2<T: Hash>(&self, t: &T) -> u64 {
        let mut hasher = self.xx().1;
        t.hash(&mut hasher);
        hasher.finish()
    }

    fn hashes<T: Hash>(&self, t: &T) -> (u64, u64) {
        (self.hash1(t), self.hash2(t))
    }

    fn double_hash(&self, hashes: (u64, u64), i_k: u64) -> u64 {
        hashes.0.wrapping_add(i_k.wrapping_mul(hashes.1)) % self.bitmap_size
    }

    /// Clear all of the bits in the filter, removing all keys from the set
    pub fn clear(&mut self) {
        self.bitmap.clear()
    }

    fn xx_new() -> XxHash64 {
        XxHash64::with_seed(rand::random())
    }
}

#[test]
fn bloom_test_add() {
    let mut bloom = Bloom::new(100, 10);
    let key: u64 = rand::random();
    assert_eq!(bloom.check(&key), false);
    bloom.add(&key);
    assert_eq!(bloom.check(&key), true);
}

#[test]
fn bloom_test_check_and_add() {
    let mut bloom = Bloom::new(100, 10);
    let key: u64 = rand::random();
    assert_eq!(bloom.check_and_add(&key), false);
    assert_eq!(bloom.check_and_add(&key), true);
}

#[test]
fn bloom_test_clear() {
    let mut bloom = Bloom::new(100, 10);
    let key: u64 = rand::random();
    bloom.add(&key);
    assert_eq!(bloom.check(&key), true);
    bloom.clear();
    assert_eq!(bloom.check(&key), false);
}

#[test]
fn bloom_test_load() {
    let mut original = Bloom::new(100, 10);
    let key: u64 = rand::random();
    original.add(&key);
    assert_eq!(original.check(&key), true);

    let cloned = Bloom::from_existing(
        &original.bitmap().to_bytes(),
        original.number_of_bits(),
        original.number_of_hash_functions(),
        original.xx(),
    );
    assert_eq!(cloned.check(&key), true);
}

#[test]
fn bloom_test_load_struct() {
    let mut original = Bloom::new(100, 10);
    let key: u64 = rand::random();
    original.add(&key);
    assert_eq!(original.check(&key), true);

    let cloned = Bloom::from_existing_struct(&original);
    assert_eq!(cloned.check(&key), true);
}
