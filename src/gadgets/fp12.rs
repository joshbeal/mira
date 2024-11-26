use rand_core::{OsRng, RngCore};
use std::cmp;

use halo2_proofs::{
    arithmetic::CurveAffine,
    circuit::{Chip, Value},
    plonk::Error,
};

use crate::{
    ff::{Field, PrimeField, PrimeFieldBits},
    main_gate::{AssignedValue, MainGate, MainGateConfig, RegionCtx},
};

/// Represent Fp12 point as FqPoint with degree = 12
/// `Fp12 = Fp2[w] / (w^6 - u - xi)`
/// This implementation assumes p = 3 (mod 4) in order for the polynomial u^2 + 1 to
/// be irreducible over Fp; i.e., in order for -1 to not be a square (quadratic residue) in Fp
/// This means we store an Fp12 point as `\sum_{i = 0}^6 (a_{i0} + a_{i1} * u) * w^i`
/// This is encoded in an FqPoint of degree 12 as `(a_{00}, ..., a_{50}, a_{01}, ..., a_{51})`

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Tuple12<C: CurveAffine> {
    pub elements: [C::Base; 12],
}

impl<C: CurveAffine> Default for Tuple12<C> {
    fn default() -> Self {
        Self {
            elements: [C::Base::ZERO; 12],
        }
    }
}

impl<C: CurveAffine> Tuple12<C> {
    pub fn one() -> Self {
        let mut result = [C::Base::ZERO; 12];
        result[0] = C::Base::ONE;
        Self { elements: result }
    }

    pub fn add(&self, other: &Tuple12<C>) -> Self {
        let elements = self
            .elements
            .iter()
            .zip(other.elements.iter())
            .map(|(&a, &b)| a + b)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        Self { elements }
    }

    pub fn negate(&self) -> Self {
        let elements = self
            .elements
            .iter()
            .map(|&a| -a)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        Self { elements }
    }

    pub fn mul<const XI_0: i64>(&self, other: &Tuple12<C>) -> Self {
        let mut a0b0_coeffs = [C::Base::ZERO; 11];
        let mut a0b1_coeffs = [C::Base::ZERO; 11];
        let mut a1b0_coeffs = [C::Base::ZERO; 11];
        let mut a1b1_coeffs = [C::Base::ZERO; 11];

        for i in 0..6 {
            for j in 0..6 {
                let coeff00 = self.elements[i] * other.elements[j];
                let coeff01 = self.elements[i] * other.elements[j + 6];
                let coeff10 = self.elements[i + 6] * other.elements[j];
                let coeff11 = self.elements[i + 6] * other.elements[j + 6];

                if i + j < 11 {
                    a0b0_coeffs[i + j] += coeff00;
                    a0b1_coeffs[i + j] += coeff01;
                    a1b0_coeffs[i + j] += coeff10;
                    a1b1_coeffs[i + j] += coeff11;
                }
            }
        }

        let mut a0b0_minus_a1b1 = [C::Base::ZERO; 11];
        let mut a0b1_plus_a1b0 = [C::Base::ZERO; 11];

        for i in 0..11 {
            a0b0_minus_a1b1[i] = a0b0_coeffs[i] - a1b1_coeffs[i];
            a0b1_plus_a1b0[i] = a0b1_coeffs[i] + a1b0_coeffs[i];
        }

        let mut result = [C::Base::ZERO; 12];
        let xi_0_f = C::Base::from(XI_0 as u64);

        for i in 0..6 {
            if i < 5 {
                result[i] =
                    xi_0_f * a0b0_minus_a1b1[i + 6] + a0b0_minus_a1b1[i] - a0b1_plus_a1b0[i + 6];
            } else {
                result[i] = a0b0_minus_a1b1[i];
            }
        }

        for i in 0..6 {
            if i < 5 {
                result[i + 6] =
                    a0b1_plus_a1b0[i] + a0b0_minus_a1b1[i + 6] + xi_0_f * a0b1_plus_a1b0[i + 6];
            } else {
                result[i + 6] = a0b1_plus_a1b0[i];
            }
        }

        Self { elements: result }
    }

    pub fn scalar_mul<F: PrimeFieldBits, const XI_0: i64>(&self, scalar: &F) -> Self {
        let bits = scalar.to_le_bits();

        let split_len = std::cmp::min(bits.len(), (F::NUM_BITS - 2) as usize);
        let (incomplete_bits, complete_bits) = bits.split_at(split_len);

        let mut acc = if incomplete_bits[0] {
            self.clone()
        } else {
            Self::one()
        };

        let mut p = self.mul::<XI_0>(self);

        for bit in incomplete_bits.iter().skip(1) {
            if *bit {
                acc = acc.mul::<XI_0>(&p);
            }
            p = p.mul::<XI_0>(&p);
        }

        for bit in complete_bits.iter() {
            if *bit {
                acc = acc.mul::<XI_0>(&p);
            }
            p = p.mul::<XI_0>(&p);
        }

        acc
    }

    pub fn generator() -> Self {
        gt_generator()
    }

    // pub fn random(mut rng: impl RngCore) -> Self {
    //     let elements = std::array::from_fn(|_| C::Base::random(&mut rng));
    //     Self { elements }
    // }

    // pub fn random_vartime() -> Self {
    //     let elements = std::array::from_fn(|_| C::Base::random(&mut OsRng));
    //     Self { elements }
    // }

    pub fn random<F: PrimeFieldBits, const XI_0: i64>(mut rng: impl RngCore) -> Self {
        let generator = Self::generator();
        let random_scalar = F::random(&mut rng);
        generator.scalar_mul::<F, XI_0>(&random_scalar)
    }

    pub fn random_vartime<F: PrimeFieldBits, const XI_0: i64>() -> Self {
        let generator = Self::generator();
        let random_scalar = F::random(&mut OsRng);
        generator.scalar_mul::<F, XI_0>(&random_scalar)
    }
}

// Known generator point for BN254 GT
fn gt_generator<C: CurveAffine>() -> Tuple12<C> {
    Tuple12 {
        elements: [
            C::Base::from_str_vartime(
                "8493334370784016972005089913588211327688223499729897951716206968320726508021",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "20049218015652006197026173611347504489508678646783216776320737476707192559881",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "6565798094314091391201231504228224566495939541538094766881371862976727043038",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "12145052038566888241256672223106590273978429515702193755778990643425246950730",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "634997487638609332803583491743335852620873788902390365055086820718589720118",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "6223602427219597392892794664899549544171383137467762280768257680446283161705",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "3758435817766288188804561253838670030762970764366672594784247447067868088068",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "18059168546148152671857026372711724379319778306792011146784665080987064164612",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "14656606573936501743457633041048024656612227301473084805627390748872617280984",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "17918828665069491344039743589118342552553375221610735811112289083834142789347",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "19455424343576886430889849773367397946457449073528455097210946839000147698372",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "7484542354754424633621663080190936924481536615300815203692506276894207018007",
            )
            .unwrap(),
        ],
    }
}

#[derive(Clone, Debug)]
pub struct AssignedTuple12<C: CurveAffine> {
    pub(crate) elements: Vec<AssignedValue<C::Base>>,
}

impl<C: CurveAffine> AssignedTuple12<C> {
    pub fn new(elements: Vec<AssignedValue<C::Base>>) -> Self {
        assert_eq!(
            elements.len(),
            12,
            "AssignedTuple12 must contain exactly 12 elements."
        );
        Self { elements }
    }

    pub fn get(&self, index: usize) -> &AssignedValue<C::Base> {
        &self.elements[index]
    }

    pub fn set(&mut self, index: usize, value: AssignedValue<C::Base>) {
        self.elements[index] = value;
    }
}

impl<C: CurveAffine> AssignedTuple12<C> {
    pub fn elements(&self) -> [&AssignedValue<C::Base>; 12] {
        self.elements.iter().collect::<Vec<_>>().try_into().unwrap()
    }

    pub fn to_tuple12(&self) -> Option<Tuple12<C>> {
        let mut elements = [C::Base::ZERO; 12];
        for (i, e) in self.elements.iter().enumerate().take(12) {
            let elem = e.value().copied().unwrap();
            elements[i] = elem?;
        }
        Some(Tuple12 { elements })
    }
}

pub struct Fp12Chip<C: CurveAffine<Base = F>, F: PrimeFieldBits, const T: usize, const XI_0: i64> {
    main_gate: MainGate<C::Base, T>,
}

impl<C: CurveAffine<Base = F>, F: PrimeFieldBits, const T: usize, const XI_0: i64>
    Fp12Chip<C, F, T, XI_0>
{
    pub fn new(config: MainGateConfig<T>) -> Self {
        let main_gate = MainGate::new(config);
        Self { main_gate }
    }

    pub fn assign_tuple<AN: Into<String>>(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        annotation: impl Fn() -> AN,
        elems: Option<[C::Base; 12]>,
    ) -> Result<AssignedTuple12<C>, Error> {
        let mut assigned_values = Vec::with_capacity(12);

        for (i, &elem) in elems.unwrap_or([C::Base::ZERO; 12]).iter().enumerate() {
            let assigned_value = ctx.assign_advice(
                || format!("{}.c{}", annotation().into(), i),
                self.main_gate.config().state[i % 4],
                Value::known(elem),
            )?;
            assigned_values.push(assigned_value);
            if i % 4 == 3 {
                ctx.next();
            }
        }

        Ok(AssignedTuple12::new(assigned_values))
    }

    pub fn one(&self, ctx: &mut RegionCtx<'_, C::Base>) -> Result<AssignedTuple12<C>, Error> {
        let elems = {
            let mut id = [C::Base::ZERO; 12];
            id[0] = C::Base::ONE;
            id
        };

        let mut assigned_values = Vec::with_capacity(12);

        for (i, &elem) in elems.iter().enumerate() {
            let assigned_value = ctx.assign_advice(
                || format!("id.c{}", i),
                self.main_gate.config().state[i % 4],
                Value::known(elem),
            )?;
            assigned_values.push(assigned_value);
            if i % 4 == 3 {
                ctx.next();
            }
        }

        for (i, val) in assigned_values.iter().enumerate() {
            if i == 0 {
                self.main_gate
                    .assert_equal_const(ctx, val.clone(), C::Base::ONE)?;
            } else {
                self.main_gate
                    .assert_equal_const(ctx, val.clone(), C::Base::ZERO)?;
            }
        }

        ctx.next();

        Ok(AssignedTuple12::new(assigned_values))
    }

    pub fn negate(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        p: &AssignedTuple12<C>,
    ) -> Result<AssignedTuple12<C>, Error> {
        let mut negated_elements = Vec::with_capacity(12);

        for i in 0..12 {
            let element = &p.elements[i];
            let negated_value: Value<C::Base> = -element.value().copied();
            let negated_element = self.main_gate.apply(
                ctx,
                (
                    Some(vec![C::Base::ONE]),
                    None,
                    Some(vec![element.clone().into()]),
                ),
                None,
                (C::Base::ONE, negated_value.into()),
            )?;
            negated_elements.push(negated_element);
        }

        Ok(AssignedTuple12::new(negated_elements))
    }

    pub fn add(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        p: &AssignedTuple12<C>,
        q: &AssignedTuple12<C>,
    ) -> Result<AssignedTuple12<C>, Error> {
        let mut sum_elements = Vec::with_capacity(12);

        for i in 0..12 {
            let sum = self.main_gate.add(ctx, &p.elements[i], &q.elements[i])?;
            sum_elements.push(sum);
        }

        Ok(AssignedTuple12::new(sum_elements))
    }

    pub fn mul(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        p: &AssignedTuple12<C>,
        q: &AssignedTuple12<C>,
    ) -> Result<AssignedTuple12<C>, Error> {
        let mut a0b0_coeffs: Vec<AssignedValue<C::Base>> = Vec::with_capacity(11);
        let mut a0b1_coeffs: Vec<AssignedValue<C::Base>> = Vec::with_capacity(11);
        let mut a1b0_coeffs: Vec<AssignedValue<C::Base>> = Vec::with_capacity(11);
        let mut a1b1_coeffs: Vec<AssignedValue<C::Base>> = Vec::with_capacity(11);

        for i in 0..6 {
            for j in 0..6 {
                let coeff00 = self.main_gate.mul(ctx, &p.elements[i], &q.elements[j])?;
                let coeff01 = self
                    .main_gate
                    .mul(ctx, &p.elements[i], &q.elements[j + 6])?;
                let coeff10 = self
                    .main_gate
                    .mul(ctx, &p.elements[i + 6], &q.elements[j])?;
                let coeff11 = self
                    .main_gate
                    .mul(ctx, &p.elements[i + 6], &q.elements[j + 6])?;

                if i + j < a0b0_coeffs.len() {
                    a0b0_coeffs[i + j] = self.main_gate.add(ctx, &a0b0_coeffs[i + j], &coeff00)?;
                    a0b1_coeffs[i + j] = self.main_gate.add(ctx, &a0b1_coeffs[i + j], &coeff01)?;
                    a1b0_coeffs[i + j] = self.main_gate.add(ctx, &a1b0_coeffs[i + j], &coeff10)?;
                    a1b1_coeffs[i + j] = self.main_gate.add(ctx, &a1b1_coeffs[i + j], &coeff11)?;
                } else {
                    a0b0_coeffs.push(coeff00);
                    a0b1_coeffs.push(coeff01);
                    a1b0_coeffs.push(coeff10);
                    a1b1_coeffs.push(coeff11);
                }
            }
        }

        let mut a0b0_minus_a1b1 = Vec::with_capacity(11);
        let mut a0b1_plus_a1b0 = Vec::with_capacity(11);

        for i in 0..11 {
            let a0b0_minus_a1b1_entry =
                self.main_gate.sub(ctx, &a0b0_coeffs[i], &a1b1_coeffs[i])?;
            let a0b1_plus_a1b0_entry = self.main_gate.add(ctx, &a0b1_coeffs[i], &a1b0_coeffs[i])?;

            a0b0_minus_a1b1.push(a0b0_minus_a1b1_entry);
            a0b1_plus_a1b0.push(a0b1_plus_a1b0_entry);
        }

        let mut out_coeffs = Vec::with_capacity(12);

        for i in 0..6 {
            if i < 5 {
                // Multiply by XI_0 and then add
                let xi_0_f = C::Base::from(XI_0 as u64);
                let mul_result =
                    self.main_gate
                        .mul_by_const(ctx, &a0b0_minus_a1b1[i + 6], xi_0_f)?;
                let mut coeff = self.main_gate.add(ctx, &mul_result, &a0b0_minus_a1b1[i])?;
                coeff = self.main_gate.sub(ctx, &coeff, &a0b1_plus_a1b0[i + 6])?;
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b0_minus_a1b1[i].clone());
            }
        }

        for i in 0..6 {
            if i < 5 {
                let mut coeff =
                    self.main_gate
                        .add(ctx, &a0b1_plus_a1b0[i], &a0b0_minus_a1b1[i + 6])?;
                let xi_0_f = C::Base::from(XI_0 as u64);
                let mul_result =
                    self.main_gate
                        .mul_by_const(ctx, &a0b1_plus_a1b0[i + 6], xi_0_f)?;
                coeff = self.main_gate.add(ctx, &coeff, &mul_result)?;
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b1_plus_a1b0[i].clone());
            }
        }

        Ok(AssignedTuple12::new(out_coeffs))
    }

    pub fn scalar_mul(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        p0: &AssignedTuple12<C>,
        scalar_bits: &[AssignedValue<C::Base>],
    ) -> Result<AssignedTuple12<C>, Error> {
        // Ensure that scalar_bits has a meaningful length
        let split_len = cmp::min(scalar_bits.len(), (C::Base::NUM_BITS - 2) as usize);
        let (incomplete_bits, complete_bits) = scalar_bits.split_at(split_len);

        // Initialize acc to p0 and perform the first doubling
        let id = self.one(ctx)?;
        let mut acc = self.conditional_select(ctx, &p0.clone(), &id, &scalar_bits[0])?;
        let mut p = self.mul(ctx, p0, p0)?;

        // Process the incomplete bits
        for bit in incomplete_bits.iter().skip(1) {
            let sum = self.mul(ctx, &acc, &p)?;
            acc = self.conditional_select(ctx, &sum, &acc, bit)?;
            p = self.mul(ctx, &p, &p)?;
        }

        // Process the remaining bits
        for bit in complete_bits {
            let sum = self.mul(ctx, &acc, &p)?;
            acc = self.conditional_select(ctx, &sum, &acc, bit)?;
            p = self.mul(ctx, &p, &p)?;
        }

        Ok(acc)
    }

    pub fn conditional_select(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        lhs: &AssignedTuple12<C>,
        rhs: &AssignedTuple12<C>,
        condition: &AssignedValue<C::Base>,
    ) -> Result<AssignedTuple12<C>, Error> {
        let mut result_elements = Vec::with_capacity(12);

        for i in 0..12 {
            let selected = self.main_gate.conditional_select(
                ctx,
                &lhs.elements[i],
                &rhs.elements[i],
                condition,
            )?;
            result_elements.push(selected);
        }

        Ok(AssignedTuple12::new(result_elements))
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use halo2_proofs::{
        circuit::{Layouter, SimpleFloorPlanner},
        plonk::{Circuit, Column, ConstraintSystem, Instance},
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use tracing_test::traced_test;

    use super::*;
    use crate::{
        create_and_verify_proof,
        ff::Field,
        halo2curves::{
            bn256,
            pasta::{pallas, EqAffine, Fp, Fq},
        },
        run_mock_prover_test,
        util::fe_to_fe_safe,
    };

    const T: usize = 5;
    const XI_0: i64 = 9;

    #[derive(Clone, Debug)]
    struct TestCircuitConfig {
        config: MainGateConfig<T>,
        instance: Column<Instance>,
    }

    struct TestCircuit<C: CurveAffine<Base = F>, F: PrimeFieldBits> {
        a: Tuple12<C>,
        b: Tuple12<C>,
        lambda: C::Scalar,
        test_case: usize, // 0: add, 1: mul, 2: scalar mul
    }

    impl<C: CurveAffine<Base = F>, F: PrimeFieldBits> TestCircuit<C, F> {
        fn new(a: Tuple12<C>, b: Tuple12<C>, lambda: C::Scalar, test_case: usize) -> Self {
            Self {
                a,
                b,
                lambda,
                test_case,
            }
        }
    }

    impl<C: CurveAffine<Base = F>, F: PrimeFieldBits> Circuit<C::Base> for TestCircuit<C, F> {
        type Config = TestCircuitConfig;
        type FloorPlanner = SimpleFloorPlanner;

        fn without_witnesses(&self) -> Self {
            TestCircuit::new(Tuple12::default(), Tuple12::default(), C::Scalar::ZERO, 0)
        }

        fn configure(meta: &mut ConstraintSystem<C::Base>) -> Self::Config {
            let instance = meta.instance_column();
            meta.enable_equality(instance);
            let config = MainGate::configure(meta);
            Self::Config { config, instance }
        }

        fn synthesize(
            &self,
            config: Self::Config,
            mut layouter: impl Layouter<C::Base>,
        ) -> Result<(), Error> {
            let fp12_chip = Fp12Chip::<C, F, T, XI_0>::new(config.config);
            let output = layouter.assign_region(
                || "fp12 test circuit",
                |region| {
                    let ctx = &mut RegionCtx::new(region, 0);
                    let mut a_elements = Vec::with_capacity(12);
                    for (i, &elem) in self.a.elements.iter().enumerate() {
                        let assigned_value = ctx.assign_advice(
                            || format!("a.c{}", i),
                            fp12_chip.main_gate.config().state[i % 4],
                            Value::known(elem),
                        )?;
                        a_elements.push(assigned_value);
                        if i % 4 == 3 && i < 11 {
                            ctx.next();
                        }
                    }
                    let a = AssignedTuple12 {
                        elements: a_elements,
                    };

                    if self.test_case == 2 {
                        let lambda: C::Base = fe_to_fe_safe(&self.lambda).unwrap();
                        let bit_len =
                            NonZeroUsize::new(lambda.to_le_bits().len()).expect("Non Zero");
                        let lambda = ctx.assign_advice(
                            || "lambda",
                            fp12_chip.main_gate.config().state[4],
                            Value::known(lambda),
                        )?;
                        ctx.next();
                        let bits = fp12_chip.main_gate.le_num_to_bits(ctx, lambda, bit_len)?;
                        fp12_chip.scalar_mul(ctx, &a, &bits)
                    } else {
                        ctx.next();

                        let mut b_elements = Vec::with_capacity(12);
                        for (i, &elem) in self.b.elements.iter().enumerate() {
                            let assigned_value = ctx.assign_advice(
                                || format!("b.c{}", i),
                                fp12_chip.main_gate.config().state[i % 4],
                                Value::known(elem),
                            )?;
                            b_elements.push(assigned_value);
                            if i % 4 == 3 && i < 11 {
                                ctx.next();
                            }
                        }
                        let b = AssignedTuple12 {
                            elements: b_elements,
                        };
                        ctx.next();
                        if self.test_case == 0 {
                            fp12_chip.add(ctx, &a, &b)
                        } else {
                            fp12_chip.mul(ctx, &a, &b)
                        }
                    }
                },
            )?;
            for i in 0..12 {
                layouter.constrain_instance(output.elements[i].cell(), config.instance, i)?;
            }
            Ok(())
        }
    }

    #[traced_test]
    #[test]
    fn test_fp12_scalar_mul() {
        let p: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let lambda = Fq::from(41);
        let r1 = p.scalar_mul::<Fq, XI_0>(&lambda);
        let mut r2 = p.clone();
        for _ in 0..40 {
            r2 = r2.mul::<XI_0>(&p);
        }
        for i in 0..12 {
            assert_eq!(
                r1.elements[i], r2.elements[i],
                "Element {} does not match",
                i
            );
        }
    }

    #[traced_test]
    #[test]
    #[ignore = "cause it takes over a minute to run"]
    fn test_fp12_add_prove() {
        println!("-----running Fp12 add Circuit-----");
        let p: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let q: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let lambda = Fq::from(1);
        let r = p.add(&q);
        let circuit = TestCircuit::new(p, q, lambda, 0);
        let public_inputs: &[&[Fp]] = &[&r.elements];

        let K: u32 = 14;
        create_and_verify_proof!(IPA, K, circuit, public_inputs, EqAffine);
        println!("-----Fp12 add circuit works fine-----");
    }

    #[traced_test]
    #[test]
    #[ignore = "cause it takes over a minute to run"]
    fn test_fp12_mul_prove() {
        println!("-----running Fp12 mul Circuit-----");
        let p: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let q: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let lambda = Fq::from(1);
        let r = p.mul::<XI_0>(&q);
        let circuit = TestCircuit::new(p, q, lambda, 1);
        let public_inputs: &[&[Fp]] = &[&r.elements];

        let K: u32 = 14;
        create_and_verify_proof!(IPA, K, circuit, public_inputs, EqAffine);
        println!("-----Fp12 mul circuit works fine-----");
    }

    #[traced_test]
    #[test]
    #[ignore = "cause it takes over a minute to run"]
    fn test_fp12_scalar_mul_prove() {
        println!("-----running Fp12 scalar mul Circuit-----");
        let p: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let lambda = Fq::random(&mut OsRng);
        let r = p.scalar_mul::<Fq, XI_0>(&lambda);
        let circuit = TestCircuit::new(p, Tuple12::default(), lambda, 2);
        let public_inputs: &[&[Fp]] = &[&r.elements];

        let K: u32 = 18;
        create_and_verify_proof!(IPA, K, circuit, public_inputs, EqAffine);

        println!("-----Fp12 scalar mul circuit works fine-----");
    }

    #[test]
    fn test_fp12_scalar_mul_mock() {
        let K: u32 = 18;
        let p: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let q: Tuple12<pallas::Affine> = Tuple12::random_vartime::<Fq, XI_0>();
        let lambda = Fq::random(&mut OsRng);
        let r = p.scalar_mul::<Fq, XI_0>(&lambda);
        let circuit = TestCircuit::new(p, q, lambda, 2);
        let public_inputs: Vec<Vec<Fp>> = vec![r.elements.to_vec()];

        run_mock_prover_test!(K, circuit, public_inputs);
    }

    fn coeffs(x: bn256::Gt) -> Vec<bn256::Fq> {
        vec![
            x.0.c0.c0.c0,
            x.0.c1.c0.c0,
            x.0.c0.c1.c0,
            x.0.c1.c1.c0,
            x.0.c0.c2.c0,
            x.0.c1.c2.c0,
            x.0.c0.c0.c1,
            x.0.c1.c0.c1,
            x.0.c0.c1.c1,
            x.0.c1.c1.c1,
            x.0.c0.c2.c1,
            x.0.c1.c2.c1,
        ]
    }

    #[test]
    fn test_impl_equivalence() {
        use std::ops::MulAssign;

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            let mut a = bn256::G1::generator();
            let ka = bn256::Fr::random(&mut rng);
            a.mul_assign(ka);

            let mut b = bn256::G2::generator();
            let kb = bn256::Fr::random(&mut rng);
            b.mul_assign(kb);

            let e = bn256::pairing(&bn256::G1Affine::from(a), &bn256::G2Affine::from(b));

            let p: Tuple12<bn256::G1Affine> = Tuple12::<bn256::G1Affine> {
                elements: coeffs(e).try_into().unwrap(),
            };
            let r1 = p.mul::<XI_0>(&p);

            let x = e + e;
            let r2: Tuple12<bn256::G1Affine> = Tuple12::<bn256::G1Affine> {
                elements: coeffs(x).try_into().unwrap(),
            };

            assert_eq!(r1.elements, r2.elements);
        }
    }

    #[test]
    fn test_generator() {
        // use num_bigint::BigUint;

        let a = bn256::G1::generator();
        let b = bn256::G2::generator();
        let e = bn256::pairing(&bn256::G1Affine::from(a), &bn256::G2Affine::from(b));
        let p = Tuple12::<bn256::G1Affine> {
            elements: coeffs(e).try_into().unwrap(),
        };

        // p.elements.iter().for_each(|elem| {
        //     let repr = elem.to_repr();
        //     let decimal_str = BigUint::from_bytes_le(&repr).to_string();
        //     println!("{}", decimal_str);
        //     let r = <bn256::G1Affine as CurveAffine>::Base::from_str_vartime(&decimal_str).unwrap();
        //     assert_eq!(*elem, r);
        // });

        let q = gt_generator::<bn256::G1Affine>();
        assert_eq!(p.elements, q.elements);
    }
}
