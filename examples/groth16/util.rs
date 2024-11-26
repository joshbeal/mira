use std::ops::MulAssign;

use halo2_proofs::halo2curves::bn256;
use mira::ff::Field;
use rand_core::OsRng;

/// A verification key in the Groth16 SNARK.
#[derive(Clone, Debug, PartialEq)]
pub struct Bn256VerifyingKey {
    /// The `alpha * G`, where `G` is the generator of `G1`.
    pub alpha_g1: bn256::G1Affine,
    /// The `alpha * H`, where `H` is the generator of `G2`.
    pub beta_g2: bn256::G2Affine,
    /// The `gamma * H`, where `H` is the generator of `G2`.
    pub gamma_g2: bn256::G2Affine,
    /// The `delta * H`, where `H` is the generator of `G2`.
    pub delta_g2: bn256::G2Affine,
    /// The `gamma^{-1} * (beta * a_i + alpha * b_i + c_i) * H`, where `H` is the generator of `G1`.
    pub gamma_abc_g1: Vec<bn256::G1Affine>,
}

/// A proof in the Groth16 SNARK.
#[derive(Clone, Debug, PartialEq)]
pub struct Bn256Proof {
    /// The `A` element in `G1`.
    pub a: bn256::G1Affine,
    /// The `B` element in `G2`.
    pub b: bn256::G2Affine,
    /// The `C` element in `G1`.
    pub c: bn256::G1Affine,
}

pub fn generate_random_cross_terms(nproofs: usize) -> Vec<Vec<bn256::Gt>> {
    (0..nproofs)
        .map(|_| {
            let mut a = bn256::G1::generator();
            let ka = bn256::Fr::random(&mut OsRng);
            a.mul_assign(ka);

            let mut b = bn256::G2::generator();
            let kb = bn256::Fr::random(&mut OsRng);
            b.mul_assign(kb);

            vec![bn256::pairing(
                &bn256::G1Affine::from(a),
                &bn256::G2Affine::from(b),
            )]
        })
        .collect::<Vec<_>>()
}
