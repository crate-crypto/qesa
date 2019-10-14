#![feature(test)]

mod copy;
mod ipa;
mod lmpa;
mod math_utils;
mod matrix;
pub mod qesa_zk;
mod transcript;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
