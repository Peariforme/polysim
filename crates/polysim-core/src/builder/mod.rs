//! Polymer chain builders.
//!
//! Each builder takes a parsed [`BigSmiles`](bigsmiles::BigSmiles) and a
//! [`BuildStrategy`] that controls chain length, then produces one or more
//! [`PolymerChain`](crate::PolymerChain) instances.

pub mod branched;
pub mod linear;
pub mod strategy;

pub use strategy::BuildStrategy;
