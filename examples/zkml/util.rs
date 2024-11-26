use std::ops::MulAssign;

use halo2_proofs::halo2curves::bn256;
use mira::ff::Field;
use rand_core::OsRng;

pub fn generate_random_g1_elems(nproofs: usize, k: usize) -> Vec<Vec<bn256::G1Affine>> {
    (0..nproofs)
        .map(|_| {
            (0..k)
                .map(|_| {
                    let mut a = bn256::G1::generator();
                    let ka = bn256::Fr::random(&mut OsRng);
                    a.mul_assign(ka);
                    bn256::G1Affine::from(a)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

pub fn generate_random_g2_elems(nproofs: usize, k: usize) -> Vec<Vec<bn256::G2Affine>> {
    (0..nproofs)
        .map(|_| {
            (0..k)
                .map(|_| {
                    let mut b = bn256::G2::generator();
                    let kb = bn256::Fr::random(&mut OsRng);
                    b.mul_assign(kb);
                    bn256::G2Affine::from(b)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

pub fn generate_random_cross_terms(nproofs: usize, k: usize) -> Vec<Vec<bn256::Gt>> {
    (0..nproofs)
        .map(|_| {
            (0..k)
                .map(|_| {
                    let mut a = bn256::G1::generator();
                    let ka = bn256::Fr::random(&mut OsRng);
                    a.mul_assign(ka);

                    let mut b = bn256::G2::generator();
                    let kb = bn256::Fr::random(&mut OsRng);
                    b.mul_assign(kb);

                    bn256::pairing(&bn256::G1Affine::from(a), &bn256::G2Affine::from(b))
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}
