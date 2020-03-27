//! Coalescent process as described in [Coalescent Theory](https://en.wikipedia.org/wiki/Coalescent_theory)

pub use coalescent::*;

pub mod coalescent;

pub mod traits;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
