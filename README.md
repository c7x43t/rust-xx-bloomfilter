## xx-bloomfilter

Hard fork of https://github.com/jedisct1/rust-bloom-filter. Reworked the most of the internals to make tbe algorithm
cleaner and more efficient. Uses extremly fast [XxHash64](https://docs.rs/twox-hash/1.4.2/twox_hash/struct.XxHash64.html) for hashing.

### Usage

*In your Cargo.toml*

```toml
[dependencies]
xx-bloomfilter = "0.10.0"
```

Initialize with expected number of items and a wanted false positive rates

```rust
extern crate xx_bloomfilter;
extern crate rand;

use xx_bloomfilter::Bloom;

fn main () {

    let mut bloom = Bloom::new_with_rate(1_000_000, 1e-6);
    let item: u64 = rand::random();

    assert_eq!(false, bloom.check_and_add(&item));
    assert_eq!(true, bloom.check(&item));

    // Clear all values
    bloom.clear();

    assert_eq!(false, bloom.check_and_add(&item));
}
```


