//! # polysim-core
//!
//! Polymer structure generator and physical property simulator built on top of
//! [BigSMILES](https://olsenlabmit.github.io/BigSMILES/docs/line_notation.html).
//!
//! ## Overview
//!
//! polysim-core turns a **BigSMILES** string into one or more concrete polymer
//! chains (represented as SMILES) and computes physical/chemical properties on
//! them.
//!
//! The typical workflow is:
//!
//! 1. **Parse** a BigSMILES string with [`parse`].
//! 2. **Build** a chain with one of the builders in [`builder`].
//! 3. **Compute** properties via [`properties`].
//!
//! ## Quick start
//!
//! ```rust
//! use polysim_core::{parse, builder::{linear::LinearBuilder, BuildStrategy}};
//!
//! // Polyethylene — 10 repeat units
//! let bs = parse("{[]CC[]}").unwrap();
//! let chain = LinearBuilder::new(bs, BuildStrategy::ByRepeatCount(10))
//!     .homopolymer()
//!     .unwrap();
//!
//! assert_eq!(chain.repeat_count, 10);
//! assert_eq!(chain.smiles, "CCCCCCCCCCCCCCCCCCCC");
//! ```

pub mod builder;
pub mod error;
pub mod polymer;
pub mod properties;

pub use bigsmiles::{parse, BigSmiles};
pub use builder::BuildStrategy;
pub use error::PolySimError;
pub use polymer::PolymerChain;
