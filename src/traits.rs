use crate::ff::PrimeFieldBits;
use halo2_proofs::arithmetic::CurveAffine;

/// A trait alias for CurveAffine types with a Base field that supports bit decomposition
pub trait BDCurveAffine: CurveAffine<Base: PrimeFieldBits, ScalarExt: PrimeFieldBits> {}

// Implement it for all types that satisfy the constraint
impl<C: CurveAffine<Base: PrimeFieldBits, ScalarExt: PrimeFieldBits>> BDCurveAffine for C {}
