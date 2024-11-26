use rand_core::OsRng;
use subtle::CtOption;

use halo2_proofs::{
    arithmetic::CurveAffine,
    circuit::{Chip, Value},
    plonk::Error,
};

use crate::{
    ff::{Field, PrimeFieldBits},
    main_gate::{AssignedValue, MainGate, MainGateConfig, RegionCtx},
};

/// Represent Fp2 point as `FieldVector` with degree = 2
/// `Fp2 = Fp[u] / (u^2 + 1)`
/// This implementation assumes p = 3 (mod 4) in order for the polynomial u^2 + 1 to be irreducible over Fp; i.e., in order for -1 to not be a square (quadratic residue) in Fp
/// This means we store an Fp2 point as `a_0 + a_1 * u` where `a_0, a_1 in Fp`

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Tuple2<C: CurveAffine> {
    pub c0: C::Base,
    pub c1: C::Base,
}

impl<C: CurveAffine> Default for Tuple2<C> {
    fn default() -> Self {
        Self {
            c0: C::Base::ZERO,
            c1: C::Base::ZERO,
        }
    }
}

impl<C: CurveAffine> Tuple2<C> {
    pub fn new(c0: C::Base, c1: C::Base) -> Self {
        Self { c0, c1 }
    }

    pub fn from(value: C::Base) -> Self {
        Self {
            c0: value,
            c1: C::Base::ZERO,
        }
    }

    /// Add two elements
    pub fn add(&self, other: &Tuple2<C>) -> Self {
        let c0 = self.c0 + other.c0;
        let c1 = self.c0 + other.c1;
        Self { c0, c1 }
    }

    /// Subtract two elements
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }

    /// Multiply two elements
    pub fn mul(&self, other: &Tuple2<C>) -> Self {
        // (a_0 + a_1 * u) * (b_0 + b_1 * u) = (a_0 b_0 - a_1 b_1) + (a_0 b_1 + a_1 b_0) * u
        let a0b0_minus_a1b1 = self.c0 * other.c0 - self.c1 * other.c1;
        let a0b1_plus_a1b0 = self.c0 * other.c1 + self.c1 * other.c0;
        Self {
            c0: a0b0_minus_a1b1,
            c1: a0b1_plus_a1b0,
        }
    }

    /// Square this element
    pub fn double(&self) -> Self {
        Self::mul(self, self)
    }

    /// Negate this element
    pub fn negate(&self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
        }
    }

    /// Invert this element
    // pub fn invert(&self) -> Option<Self> {
    //     let norm = self.c0.square() + self.c1.square();
    //     norm.invert().map(|norm_inv| Self {
    //         c0: self.c0 * norm_inv,
    //         c1: -self.c1 * norm_inv,
    //     })
    // }
    pub fn invert(&self) -> Option<Self> {
        let norm = self.c0.square() + self.c1.square();
        CtOption::new(norm, !norm.is_zero())
            .and_then(|norm| {
                norm.invert().map(|norm_inv| Self {
                    c0: self.c0 * norm_inv,
                    c1: -self.c1 * norm_inv,
                })
            })
            .into()
    }

    /// Generate random element
    pub fn random_vartime() -> Self {
        let c0 = C::Base::random(&mut OsRng);
        let c1 = C::Base::random(&mut OsRng);
        Self { c0, c1 }
    }
}

#[derive(Clone, Debug)]
pub struct AssignedTuple2<C: CurveAffine> {
    pub(crate) c0: AssignedValue<C::Base>,
    pub(crate) c1: AssignedValue<C::Base>,
}

impl<C: CurveAffine> AssignedTuple2<C> {
    pub fn elements(&self) -> (&AssignedValue<C::Base>, &AssignedValue<C::Base>) {
        (&self.c0, &self.c1)
    }

    pub fn elements_values(&self) -> Option<(C::Base, C::Base)> {
        let c0 = self.c0.value().copied().unwrap();
        let c1 = self.c1.value().copied().unwrap();

        Some((c0?, c1?))
    }
}

pub struct Fp2Chip<C: CurveAffine<Base = F>, F: PrimeFieldBits, const T: usize> {
    main_gate: MainGate<C::Base, T>,
}

impl<C: CurveAffine<Base = F>, F: PrimeFieldBits, const T: usize> Fp2Chip<C, F, T> {
    pub fn new(config: MainGateConfig<T>) -> Self {
        let main_gate = MainGate::new(config);
        Self { main_gate }
    }

    pub fn assign_tuple<AN: Into<String>>(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        annotation: impl Fn() -> AN,
        elems: Option<(C::Base, C::Base)>,
    ) -> Result<AssignedTuple2<C>, Error> {
        let c0 = ctx.assign_advice(
            || format!("{}.c0", annotation().into()),
            self.main_gate.config().state[0],
            Value::known(elems.map_or(C::Base::ZERO, |c| c.0)),
        )?;
        let c1 = ctx.assign_advice(
            || format!("{}.c1", annotation().into()),
            self.main_gate.config().state[1],
            Value::known(elems.map_or(C::Base::ZERO, |c| c.1)),
        )?;
        ctx.next();

        Ok(AssignedTuple2 { c0, c1 })
    }

    pub fn negate(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        p: &AssignedTuple2<C>,
    ) -> Result<AssignedTuple2<C>, Error> {
        let c0 = &p.c0;
        let c1 = &p.c1;
        let c0_minus_val: Value<C::Base> = -c0.value().copied();
        let c1_minus_val: Value<C::Base> = -c1.value().copied();
        let c0 = self.main_gate.apply(
            ctx,
            (Some(vec![C::Base::ONE]), None, Some(vec![c0.into()])),
            None,
            (C::Base::ONE, c0_minus_val.into()),
        )?;
        let c1 = self.main_gate.apply(
            ctx,
            (Some(vec![C::Base::ONE]), None, Some(vec![c1.into()])),
            None,
            (C::Base::ONE, c1_minus_val.into()),
        )?;
        Ok(AssignedTuple2 { c0, c1 })
    }

    pub fn add(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        p: &AssignedTuple2<C>,
        q: &AssignedTuple2<C>,
    ) -> Result<AssignedTuple2<C>, Error> {
        let c0 = self.main_gate.add(ctx, &p.c0, &q.c0)?;
        let c1 = self.main_gate.add(ctx, &p.c1, &q.c1)?;
        Ok(AssignedTuple2 { c0, c1 })
    }

    fn mul(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        p: &AssignedTuple2<C>,
        q: &AssignedTuple2<C>,
    ) -> Result<AssignedTuple2<C>, Error> {
        let x00 = self.main_gate.mul(ctx, &p.c0, &q.c0)?;
        let x11 = self.main_gate.mul(ctx, &p.c1, &q.c1)?;
        let c0 = self.main_gate.sub(ctx, &x00, &x11)?;

        let x01 = self.main_gate.mul(ctx, &p.c0, &q.c1)?;
        let x10 = self.main_gate.mul(ctx, &p.c1, &q.c0)?;
        let c1 = self.main_gate.add(ctx, &x01, &x10)?;

        Ok(AssignedTuple2 { c0, c1 })
    }

    pub fn conditional_select(
        &self,
        ctx: &mut RegionCtx<'_, C::Base>,
        lhs: &AssignedTuple2<C>,
        rhs: &AssignedTuple2<C>,
        condition: &AssignedValue<C::Base>,
    ) -> Result<AssignedTuple2<C>, Error> {
        Ok(AssignedTuple2 {
            c0: self
                .main_gate
                .conditional_select(ctx, &lhs.c0, &rhs.c0, condition)?,
            c1: self
                .main_gate
                .conditional_select(ctx, &lhs.c1, &rhs.c1, condition)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use halo2_proofs::{
        circuit::{Layouter, SimpleFloorPlanner},
        plonk::{Circuit, Column, ConstraintSystem, Instance},
    };
    use tracing_test::traced_test;

    use super::*;
    use crate::{
        create_and_verify_proof,
        halo2curves::pasta::{pallas, EqAffine, Fp},
        run_mock_prover_test,
    };

    const T: usize = 4;
    #[derive(Clone, Debug)]
    struct TestCircuitConfig {
        config: MainGateConfig<T>,
        instance: Column<Instance>,
    }

    struct TestCircuit<C: CurveAffine<Base = F>, F: PrimeFieldBits> {
        a: Tuple2<C>,
        b: Tuple2<C>,
        test_case: usize, // 0: add, 1: mul
    }
    impl<C: CurveAffine<Base = F>, F: PrimeFieldBits> TestCircuit<C, F> {
        fn new(a: Tuple2<C>, b: Tuple2<C>, test_case: usize) -> Self {
            Self { a, b, test_case }
        }
    }

    impl<C: CurveAffine<Base = F>, F: PrimeFieldBits> Circuit<C::Base> for TestCircuit<C, F> {
        type Config = TestCircuitConfig;
        type FloorPlanner = SimpleFloorPlanner;

        fn without_witnesses(&self) -> Self {
            TestCircuit::new(Tuple2::default(), Tuple2::default(), 0)
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
            let fp2_chip = Fp2Chip::<C, F, T>::new(config.config);
            let output = layouter.assign_region(
                || "ecc test circuit",
                |region| {
                    let ctx = &mut RegionCtx::new(region, 0);
                    let a0 = ctx.assign_advice(
                        || "a.c0",
                        fp2_chip.main_gate.config().state[0],
                        Value::known(self.a.c0),
                    )?;
                    let a1 = ctx.assign_advice(
                        || "a.c1",
                        fp2_chip.main_gate.config().state[1],
                        Value::known(self.a.c1),
                    )?;
                    let a = AssignedTuple2 { c0: a0, c1: a1 };
                    if self.test_case == 0 {
                        let b0 = ctx.assign_advice(
                            || "b.c0",
                            fp2_chip.main_gate.config().state[2],
                            Value::known(self.b.c0),
                        )?;
                        let b1 = ctx.assign_advice(
                            || "b.c1",
                            fp2_chip.main_gate.config().state[3],
                            Value::known(self.b.c1),
                        )?;
                        let b = AssignedTuple2 { c0: b0, c1: b1 };
                        ctx.next();
                        fp2_chip.add(ctx, &a, &b)
                    } else {
                        let b0 = ctx.assign_advice(
                            || "b.c0",
                            fp2_chip.main_gate.config().state[2],
                            Value::known(self.b.c0),
                        )?;
                        let b1 = ctx.assign_advice(
                            || "b.c1",
                            fp2_chip.main_gate.config().state[3],
                            Value::known(self.b.c1),
                        )?;
                        let b = AssignedTuple2 { c0: b0, c1: b1 };
                        ctx.next();
                        fp2_chip.mul(ctx, &a, &b)
                    }
                },
            )?;
            layouter.constrain_instance(output.c0.cell(), config.instance, 0)?;
            layouter.constrain_instance(output.c1.cell(), config.instance, 1)?;
            Ok(())
        }
    }

    #[traced_test]
    #[test]
    #[ignore = "cause it takes over a minute to run"]
    fn test_fp2_op() {
        println!("-----running Fp2 Circuit-----");
        let p: Tuple2<pallas::Affine> = Tuple2::random_vartime();
        let q: Tuple2<pallas::Affine> = Tuple2::random_vartime();
        let r = p.mul(&q);
        let circuit = TestCircuit::new(p, q, 1);
        let public_inputs: &[&[Fp]] = &[&[r.c0, r.c1]];

        let K: u32 = 14;
        create_and_verify_proof!(IPA, K, circuit, public_inputs, EqAffine);
        println!("-----Fp2 circuit works fine-----");
    }

    #[test]
    fn test_fp2_mock() {
        let K: u32 = 14;
        let p: Tuple2<pallas::Affine> = Tuple2::random_vartime();
        let q: Tuple2<pallas::Affine> = Tuple2::random_vartime();
        let r = p.mul(&q);
        let circuit = TestCircuit::new(p, q, 1);
        let public_inputs = vec![vec![r.c0, r.c1]];
        run_mock_prover_test!(K, circuit, public_inputs);
    }
}
