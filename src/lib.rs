#![allow(clippy::type_complexity)]
#![allow(unstable_name_collisions)]
#![allow(dead_code)] // TODO: remove it later
#![allow(non_snake_case)] // UPPER_CASE is used for ease of compatibility with Nova documentation

pub mod commitment;
pub mod constants;
pub mod digest;
pub mod fft;
pub mod gadgets;
pub mod ivc;
pub mod main_gate;
pub mod nifs;
pub mod plonk;
pub mod polynomial;
pub mod poseidon;
pub mod sps;
pub mod table;
pub mod traits;
pub mod util;

pub mod error;

pub use halo2_proofs::{
    self, halo2curves,
    halo2curves::{ff, group},
};
