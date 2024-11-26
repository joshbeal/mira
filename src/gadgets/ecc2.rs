use rand_core::{OsRng, RngCore};

use crate::{
    ff::{PrimeField, PrimeFieldBits},
    gadgets::fp2::Tuple2,
    main_gate::{AssignedValue, MainGate, MainGateConfig, RegionCtx},
};
use halo2_proofs::{
    arithmetic::{CurveAffine, Field},
    circuit::{Chip, Value},
    plonk::Error,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct G2Point<C: CurveAffine<ScalarExt = F>, F: PrimeFieldBits> {
    pub x: Tuple2<C>,
    pub y: Tuple2<C>,
    pub is_inf: bool,
}

impl<C: CurveAffine<ScalarExt = F>, F: PrimeFieldBits> Default for G2Point<C, F>
where
    C::ScalarExt: PrimeFieldBits,
{
    fn default() -> Self {
        Self {
            x: Tuple2::default(),
            y: Tuple2::default(),
            is_inf: true,
        }
    }
}

impl<C: CurveAffine<ScalarExt = F>, F: PrimeFieldBits> G2Point<C, F>
where
    C::ScalarExt: PrimeFieldBits,
{
    pub fn negate(&self) -> Self {
        if self.is_inf {
            return self.clone();
        }
        Self {
            x: self.x.clone(),
            y: self.y.negate(),
            is_inf: self.is_inf,
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        if self.is_inf {
            return other.clone();
        }
        if other.is_inf {
            return self.clone();
        }
        if self.x == other.x {
            if self.y == other.y {
                return self.double();
            } else {
                return Self::default();
            }
        }
        let lambda = other
            .y
            .sub(&self.y)
            .mul(&other.x.sub(&self.x).invert().unwrap());
        let x3 = lambda.double().sub(&self.x).sub(&other.x);
        let y3 = lambda.mul(&self.x.sub(&x3)).sub(&self.y);
        Self {
            x: x3,
            y: y3,
            is_inf: false,
        }
    }

    pub fn double(&self) -> Self {
        if self.is_inf {
            return Self::default();
        }
        if self.y == Tuple2::default() {
            return Self::default();
        }

        // Compute lambda = 3 * x^2 / (2 * y)
        let x_squared = self.x.double();
        let three_x_squared = x_squared.mul(&Tuple2::from(C::Base::from(3)));
        let two_y = self.y.mul(&Tuple2::from(C::Base::from(2)));
        let lambda = three_x_squared.mul(&two_y.invert().unwrap());

        // Compute x3 = lambda^2 - 2x
        let x3 = lambda
            .double()
            .sub(&self.x.mul(&Tuple2::from(C::Base::from(2))));

        // Compute y3 = lambda * (x - x3) - y
        let y3 = lambda.mul(&self.x.sub(&x3)).sub(&self.y);

        Self {
            x: x3,
            y: y3,
            is_inf: false,
        }
    }

    pub fn scalar_mul(&self, scalar: &C::ScalarExt) -> Self {
        let bits = scalar.to_le_bits();
        let split_len = std::cmp::min(bits.len(), (C::ScalarExt::NUM_BITS - 2) as usize);
        let (incomplete_bits, complete_bits) = bits.split_at(split_len);

        // (1) Assume p is not infinity
        let mut acc = self.clone(); // Initial accumulator set to the point itself
        let mut double_p = self.double(); // Double the point to get 2P

        // Handle incomplete bits
        for bit in incomplete_bits.iter().skip(1) {
            let sum = acc.add(&double_p); // Sum current acc and 2P
            acc = if *bit { sum } else { acc }; // Conditional selection based on bit
            double_p = double_p.double(); // Double the point again for the next iteration
        }

        // Correct if the first bit is 0
        if !bits[0] {
            let neg_p = self.negate(); // Negate the original point
            acc = acc.add(&neg_p); // Subtract the original point from the accumulator
        }

        // (2) Handle the case where `self` is infinity
        if self.is_inf {
            return Self::default(); // Return point at infinity
        }

        // (3) Handle complete bits
        for bit in complete_bits.iter() {
            let sum = acc.add(&double_p); // Add current acc and double_p
            acc = if *bit { sum } else { acc }; // Conditional selection based on bit
            double_p = double_p.double(); // Double the point again for the next iteration
        }

        acc // Return the final accumulator
    }

    pub fn random(rng: impl RngCore) -> Self {
        let generator = g2_generator();
        let random_scalar = C::ScalarExt::random(rng);
        generator.scalar_mul(&random_scalar)
    }

    pub fn random_vartime() -> Self {
        let generator = g2_generator();
        let random_scalar = C::ScalarExt::random(&mut OsRng);
        generator.scalar_mul(&random_scalar)
    }
}

// Known generator point for BN254 G2
fn g2_generator<C: CurveAffine<ScalarExt = F>, F: PrimeFieldBits>() -> G2Point<C, F> {
    G2Point {
        x: Tuple2::new(
            C::Base::from_str_vartime(
                "10857046999023057135944570762232829481370756359578518086990519993285655852781",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "11559732032986387107991004021392285783925812861821192530917403151452391805634",
            )
            .unwrap(),
        ),
        y: Tuple2::new(
            C::Base::from_str_vartime(
                "8495653923123431417604973247489272438418190587263600148770280649306958101930",
            )
            .unwrap(),
            C::Base::from_str_vartime(
                "4082367875863433681332203403145435568316851327593401208105741076214120093531",
            )
            .unwrap(),
        ),
        is_inf: false,
    }
}

// G2 point representation using G1's base field
#[derive(Clone, Debug)]
pub struct AssignedG2Point<C: CurveAffine> {
    pub(crate) x: (AssignedValue<C::Base>, AssignedValue<C::Base>),
    pub(crate) y: (AssignedValue<C::Base>, AssignedValue<C::Base>),
}

impl<C: CurveAffine<Base = F>, F: Field> AssignedG2Point<C>
where
    C::ScalarExt: PrimeFieldBits,
{
    pub fn coordinates(
        &self,
    ) -> (
        &(AssignedValue<F>, AssignedValue<F>),
        &(AssignedValue<F>, AssignedValue<F>),
    ) {
        (&self.x, &self.y)
    }

    pub fn coordinates_values(&self) -> Option<(Tuple2<C>, Tuple2<C>)> {
        let x0 = self.x.0.value().copied().unwrap();
        let x1 = self.x.1.value().copied().unwrap();
        let y0 = self.y.0.value().copied().unwrap();
        let y1 = self.y.1.value().copied().unwrap();

        Some((Tuple2::new(x0?, x1?), Tuple2::new(y0?, y1?)))
    }

    pub fn to_curve(&self) -> Option<G2Point<C, C::ScalarExt>> {
        let x0 = self.x.0.value().copied().unwrap();
        let x1 = self.x.1.value().copied().unwrap();
        let y0 = self.y.0.value().copied().unwrap();
        let y1 = self.y.1.value().copied().unwrap();

        Some(G2Point {
            x: Tuple2::new(x0?, x1?),
            y: Tuple2::new(y0?, y1?),
            is_inf: false, // TODO(jbeal): Add check for infinity point
        })
    }
}

pub struct G2EccChip<C: CurveAffine<Base = F>, F: PrimeField, const T: usize> {
    main_gate: MainGate<C::Base, T>,
}

impl<C: CurveAffine<Base = F>, F: PrimeField, const T: usize> G2EccChip<C, F, T>
where
    C::ScalarExt: PrimeFieldBits,
{
    pub fn new(config: MainGateConfig<T>) -> Self {
        let main_gate = MainGate::new(config);
        Self { main_gate }
    }

    pub fn assign_g2_point<AN: Into<String>>(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        annotation: impl Fn() -> AN,
        coords: Option<((F, F), (F, F))>,
    ) -> Result<AssignedG2Point<C>, Error> {
        let x0 = ctx.assign_advice(
            || format!("{}.x0", annotation().into()),
            self.main_gate.config().state[0],
            Value::known(coords.map_or(F::ZERO, |c| c.0 .0)),
        )?;
        let x1 = ctx.assign_advice(
            || format!("{}.x1", annotation().into()),
            self.main_gate.config().state[1],
            Value::known(coords.map_or(F::ZERO, |c| c.0 .1)),
        )?;
        let y0 = ctx.assign_advice(
            || format!("{}.y0", annotation().into()),
            self.main_gate.config().state[2],
            Value::known(coords.map_or(F::ZERO, |c| c.1 .0)),
        )?;
        let y1 = ctx.assign_advice(
            || format!("{}.y1", annotation().into()),
            self.main_gate.config().state[3],
            Value::known(coords.map_or(F::ZERO, |c| c.1 .1)),
        )?;
        ctx.next();

        Ok(AssignedG2Point {
            x: (x0, x1),
            y: (y0, y1),
        })
    }

    // G2 point addition
    pub fn add_g2(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        p: &AssignedG2Point<C>,
        q: &AssignedG2Point<C>,
    ) -> Result<AssignedG2Point<C>, Error> {
        // Check if either point is the point at infinity
        let p_is_infinity = self.fq2_is_zero(ctx, (&p.x.0, &p.x.1))?;
        let q_is_infinity = self.fq2_is_zero(ctx, (&q.x.0, &q.x.1))?;

        // Compute x_equal and y_equal to check if points are the same (P == Q)
        let x_equal = self.fq2_is_equal(ctx, (&p.x.0, &p.x.1), (&q.x.0, &q.x.1))?;
        let y_equal = self.fq2_is_equal(ctx, (&p.y.0, &p.y.1), (&q.y.0, &q.y.1))?;
        let points_equal = self.main_gate.mul(ctx, &x_equal, &y_equal)?;

        // Compute lambda = (y2 - y1) / (x2 - x1)
        let y2_minus_y1 = self.fq2_sub(ctx, (&q.y.0, &q.y.1), (&p.y.0, &p.y.1))?;
        let x2_minus_x1 = self.fq2_sub(ctx, (&q.x.0, &q.x.1), (&p.x.0, &p.x.1))?;
        let x2_minus_x1_inv = self.fq2_inv(ctx, (&x2_minus_x1.0, &x2_minus_x1.1))?;
        let lambda = self.fq2_mul(
            ctx,
            (&y2_minus_y1.0, &y2_minus_y1.1),
            (&x2_minus_x1_inv.0, &x2_minus_x1_inv.1),
        )?;

        // Compute x3 = lambda^2 - x1 - x2
        let lambda_squared = self.fq2_mul(ctx, (&lambda.0, &lambda.1), (&lambda.0, &lambda.1))?;
        let lambda_squared_minus_x1 = self.fq2_sub(
            ctx,
            (&lambda_squared.0, &lambda_squared.1),
            (&p.x.0, &p.x.1),
        )?;
        let x3 = self.fq2_sub(
            ctx,
            (&lambda_squared_minus_x1.0, &lambda_squared_minus_x1.1),
            (&q.x.0, &q.x.1),
        )?;

        // Compute y3 = lambda * (x1 - x3) - y1
        let x1_minus_x3 = self.fq2_sub(ctx, (&p.x.0, &p.x.1), (&x3.0, &x3.1))?;
        let lambda_times_x1_minus_x3 = self.fq2_mul(
            ctx,
            (&lambda.0, &lambda.1),
            (&x1_minus_x3.0, &x1_minus_x3.1),
        )?;
        let y3 = self.fq2_sub(
            ctx,
            (&lambda_times_x1_minus_x3.0, &lambda_times_x1_minus_x3.1),
            (&p.y.0, &p.y.1),
        )?;

        let normal_result = AssignedG2Point { x: x3, y: y3 };
        let zero_point = self.zero_g2(ctx)?;

        // Special case logic: Handle P == -Q (opposite points) -> return the point at infinity
        let result_if_not_opposite =
            self.conditional_select_g2(ctx, &zero_point, &normal_result, &x_equal)?;

        // Special case logic: Handle P == Q (point doubling)
        let doubled_p = self.double_g2(ctx, p)?;
        let result_if_not_equal =
            self.conditional_select_g2(ctx, &doubled_p, &result_if_not_opposite, &points_equal)?;

        // Special case logic: Handle P being the point at infinity
        let result_if_not_p_infinity =
            self.conditional_select_g2(ctx, q, &result_if_not_equal, &p_is_infinity)?;

        // Special case logic: Handle Q being the point at infinity
        let result_if_not_q_infinity =
            self.conditional_select_g2(ctx, p, &result_if_not_p_infinity, &q_is_infinity)?;

        // Final result
        Ok(result_if_not_q_infinity)
    }

    // G2 point doubling
    pub fn double_g2(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        p: &AssignedG2Point<C>,
    ) -> Result<AssignedG2Point<C>, Error> {
        // Check if the point is the point at infinity
        let is_infinity = self.fq2_is_zero(ctx, (&p.x.0, &p.x.1))?;

        // Check if y == 0
        let y_is_zero = self.fq2_is_zero(ctx, (&p.y.0, &p.y.1))?;

        // Compute lambda = 3 * x^2 / (2 * y)
        let x_squared = self.fq2_mul(ctx, (&p.x.0, &p.x.1), (&p.x.0, &p.x.1))?;

        // Assign the scalar value 3 for the first element, and 0 for the second element
        let three = self.main_gate.assign_value(ctx, Value::known(F::from(3)))?;
        let zero = self.main_gate.assign_value(ctx, Value::known(F::ZERO))?;

        // Multiply (3, 0) with the Fq2 value of x^2
        let three_x_squared = self.fq2_mul(ctx, (&three, &zero), (&x_squared.0, &x_squared.1))?;

        // Assign the scalar value 2 for the first element, and 0 for the second element
        let two = self.main_gate.assign_value(ctx, Value::known(F::from(2)))?;
        let zero = self.main_gate.assign_value(ctx, Value::known(F::ZERO))?;

        // Multiply (2, 0) with the Fq2 value of y
        let two_y = self.fq2_mul(ctx, (&two, &zero), (&p.y.0, &p.y.1))?;
        let two_y_inv = self.fq2_inv_or_zero(ctx, (&two_y.0, &two_y.1))?;

        let lambda = self.fq2_mul(
            ctx,
            (&three_x_squared.0, &three_x_squared.1),
            (&two_y_inv.0, &two_y_inv.1),
        )?;

        // Compute x3 = lambda^2 - 2x
        let lambda_squared = self.fq2_mul(ctx, (&lambda.0, &lambda.1), (&lambda.0, &lambda.1))?;
        let two_x = self.fq2_mul(ctx, (&two, &zero), (&p.x.0, &p.x.1))?;
        let x3 = self.fq2_sub(
            ctx,
            (&lambda_squared.0, &lambda_squared.1),
            (&two_x.0, &two_x.1),
        )?;

        // Compute y3 = lambda * (x - x3) - y
        let x_minus_x3 = self.fq2_sub(ctx, (&p.x.0, &p.x.1), (&x3.0, &x3.1))?;
        let lambda_times_diff =
            self.fq2_mul(ctx, (&lambda.0, &lambda.1), (&x_minus_x3.0, &x_minus_x3.1))?;
        let y3 = self.fq2_sub(
            ctx,
            (&lambda_times_diff.0, &lambda_times_diff.1),
            (&p.y.0, &p.y.1),
        )?;

        let result = AssignedG2Point { x: x3, y: y3 };
        let infinity_point = self.zero_g2(ctx)?;

        // Handle special cases
        let result_if_not_zero_y =
            self.conditional_select_g2(ctx, &infinity_point, &result, &y_is_zero)?;
        let final_result =
            self.conditional_select_g2(ctx, &infinity_point, &result_if_not_zero_y, &is_infinity)?;

        Ok(final_result)
    }

    pub fn scalar_mul(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        p: &AssignedG2Point<C>,
        scalar_bits: &[AssignedValue<F>],
    ) -> Result<AssignedG2Point<C>, Error> {
        let split_len = std::cmp::min(scalar_bits.len(), (F::NUM_BITS - 2) as usize);
        let (incomplete_bits, complete_bits) = scalar_bits.split_at(split_len);

        // (1) Assume p is not infinity

        // Assume first bit of scalar_bits is 1 for now
        let mut acc = p.clone();
        let mut double_p = self.double_g2(ctx, p)?;

        // Handle incomplete bits
        for bit in incomplete_bits.iter().skip(1) {
            let sum = self.add_g2(ctx, &acc, &double_p)?;

            acc = AssignedG2Point {
                x: (
                    self.main_gate
                        .conditional_select(ctx, &sum.x.0, &acc.x.0, bit)?,
                    self.main_gate
                        .conditional_select(ctx, &sum.x.1, &acc.x.1, bit)?,
                ),
                y: (
                    self.main_gate
                        .conditional_select(ctx, &sum.y.0, &acc.y.0, bit)?,
                    self.main_gate
                        .conditional_select(ctx, &sum.y.1, &acc.y.1, bit)?,
                ),
            };
            double_p = self.double_g2(ctx, &double_p)?;
        }

        // Correct if first bit is 0
        acc = {
            let acc_minus_initial = {
                let neg_p = self.negate_g2(ctx, p)?;
                self.add_g2(ctx, &acc, &neg_p)?
            };
            AssignedG2Point {
                x: (
                    self.main_gate.conditional_select(
                        ctx,
                        &acc.x.0,
                        &acc_minus_initial.x.0,
                        &scalar_bits[0],
                    )?,
                    self.main_gate.conditional_select(
                        ctx,
                        &acc.x.1,
                        &acc_minus_initial.x.1,
                        &scalar_bits[0],
                    )?,
                ),
                y: (
                    self.main_gate.conditional_select(
                        ctx,
                        &acc.y.0,
                        &acc_minus_initial.y.0,
                        &scalar_bits[0],
                    )?,
                    self.main_gate.conditional_select(
                        ctx,
                        &acc.y.1,
                        &acc_minus_initial.y.1,
                        &scalar_bits[0],
                    )?,
                ),
            }
        };

        // (2) Handle the case where p is infinity
        let infinity_point = self.zero_g2(ctx)?;
        let is_p_infinity = self.is_infinity_g2(ctx, p)?;
        acc = AssignedG2Point {
            x: (
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.x.0,
                    &acc.x.0,
                    &is_p_infinity,
                )?,
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.x.1,
                    &acc.x.1,
                    &is_p_infinity,
                )?,
            ),
            y: (
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.y.0,
                    &acc.y.0,
                    &is_p_infinity,
                )?,
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.y.1,
                    &acc.y.1,
                    &is_p_infinity,
                )?,
            ),
        };

        double_p = AssignedG2Point {
            x: (
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.x.0,
                    &double_p.x.0,
                    &is_p_infinity,
                )?,
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.x.1,
                    &double_p.x.1,
                    &is_p_infinity,
                )?,
            ),
            y: (
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.y.0,
                    &double_p.y.0,
                    &is_p_infinity,
                )?,
                self.main_gate.conditional_select(
                    ctx,
                    &infinity_point.y.1,
                    &double_p.y.1,
                    &is_p_infinity,
                )?,
            ),
        };

        // (3) Handle complete bits
        for bit in complete_bits.iter() {
            let sum = self.add_g2(ctx, &acc, &double_p)?;
            acc = AssignedG2Point {
                x: (
                    self.main_gate
                        .conditional_select(ctx, &sum.x.0, &acc.x.0, bit)?,
                    self.main_gate
                        .conditional_select(ctx, &sum.x.1, &acc.x.1, bit)?,
                ),
                y: (
                    self.main_gate
                        .conditional_select(ctx, &sum.y.0, &acc.y.0, bit)?,
                    self.main_gate
                        .conditional_select(ctx, &sum.y.1, &acc.y.1, bit)?,
                ),
            };
            double_p = self.double_g2(ctx, &double_p)?;
        }

        Ok(acc)
    }

    // Helper method to negate a G2 point
    fn negate_g2(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        p: &AssignedG2Point<C>,
    ) -> Result<AssignedG2Point<C>, Error> {
        let zero = self.main_gate.assign_value(ctx, Value::known(F::ZERO))?;
        let negated_y0 = self.main_gate.sub(ctx, &zero, &p.y.0)?;
        let negated_y1 = self.main_gate.sub(ctx, &zero, &p.y.1)?;
        Ok(AssignedG2Point {
            x: p.x.clone(),
            y: (negated_y0, negated_y1),
        })
    }

    // Helper method to check if a G2 point is infinity
    fn is_infinity_g2(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        p: &AssignedG2Point<C>,
    ) -> Result<AssignedValue<F>, Error> {
        let x_is_zero = self.fq2_is_zero(ctx, (&p.x.0, &p.x.1))?;
        let y_is_zero = self.fq2_is_zero(ctx, (&p.y.0, &p.y.1))?;
        self.main_gate.mul(ctx, &x_is_zero, &y_is_zero)
    }

    // Helper method to create a zero G2 point
    fn zero_g2(&self, ctx: &mut RegionCtx<'_, F>) -> Result<AssignedG2Point<C>, Error> {
        let zero = self.main_gate.assign_value(ctx, Value::known(F::ZERO))?;
        Ok(AssignedG2Point {
            x: (zero.clone(), zero.clone()),
            y: (zero.clone(), zero),
        })
    }

    // Helper method for conditional selection of G2 points
    pub fn conditional_select_g2(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: &AssignedG2Point<C>,
        b: &AssignedG2Point<C>,
        condition: &AssignedValue<F>,
    ) -> Result<AssignedG2Point<C>, Error> {
        let x0 = self
            .main_gate
            .conditional_select(ctx, &a.x.0, &b.x.0, condition)?;
        let x1 = self
            .main_gate
            .conditional_select(ctx, &a.x.1, &b.x.1, condition)?;
        let y0 = self
            .main_gate
            .conditional_select(ctx, &a.y.0, &b.y.0, condition)?;
        let y1 = self
            .main_gate
            .conditional_select(ctx, &a.y.1, &b.y.1, condition)?;
        Ok(AssignedG2Point {
            x: (x0, x1),
            y: (y0, y1),
        })
    }

    // Helper method to check if two Fq2 elements are equal
    fn fq2_is_equal(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: (&AssignedValue<F>, &AssignedValue<F>),
        b: (&AssignedValue<F>, &AssignedValue<F>),
    ) -> Result<AssignedValue<F>, Error> {
        // Check if first elements are equal
        let first_equal = self.main_gate.is_equal_term(ctx, a.0, b.0)?;

        // Check if second elements are equal
        let second_equal = self.main_gate.is_equal_term(ctx, a.1, b.1)?;

        // Both elements must be equal for Fq2 elements to be equal
        self.main_gate.mul(ctx, &first_equal, &second_equal)
    }

    // Helper method for Fq2 addition
    fn fq2_add(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: (&AssignedValue<F>, &AssignedValue<F>),
        b: (&AssignedValue<F>, &AssignedValue<F>),
    ) -> Result<(AssignedValue<F>, AssignedValue<F>), Error> {
        let c0 = self.main_gate.add(ctx, a.0, b.0)?;
        let c1 = self.main_gate.add(ctx, a.1, b.1)?;
        Ok((c0, c1))
    }

    // Helper method for Fq2 subtraction
    fn fq2_sub(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: (&AssignedValue<F>, &AssignedValue<F>),
        b: (&AssignedValue<F>, &AssignedValue<F>),
    ) -> Result<(AssignedValue<F>, AssignedValue<F>), Error> {
        let c0 = self.main_gate.sub(ctx, a.0, b.0)?;
        let c1 = self.main_gate.sub(ctx, a.1, b.1)?;
        Ok((c0, c1))
    }

    // Helper method for Fq2 multiplication
    fn fq2_mul(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: (&AssignedValue<F>, &AssignedValue<F>),
        b: (&AssignedValue<F>, &AssignedValue<F>),
    ) -> Result<(AssignedValue<F>, AssignedValue<F>), Error> {
        let ac = self.main_gate.mul(ctx, a.0, b.0)?;
        let bd = self.main_gate.mul(ctx, a.1, b.1)?;
        let a_plus_b = self.main_gate.add(ctx, a.0, a.1)?;
        let c_plus_d = self.main_gate.add(ctx, b.0, b.1)?;
        let ad_plus_bc = self.main_gate.mul(ctx, &a_plus_b, &c_plus_d)?;
        let ac_plus_bd = self.main_gate.add(ctx, &ac, &bd)?;
        let c0 = self.main_gate.sub(ctx, &ac, &bd)?;
        let c1 = self.main_gate.sub(ctx, &ad_plus_bc, &ac_plus_bd)?;
        Ok((c0, c1))
    }

    // Helper method for Fq2 inversion
    fn fq2_inv(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: (&AssignedValue<F>, &AssignedValue<F>),
    ) -> Result<(AssignedValue<F>, AssignedValue<F>), Error> {
        // 1. Compute a^2 + b^2
        let a_squared = self.main_gate.square(ctx, a.0)?;
        let b_squared = self.main_gate.square(ctx, a.1)?;
        let norm = self.main_gate.add(ctx, &a_squared, &b_squared)?;

        // 2. Compute the inverse of the norm
        let (_is_zero, norm_inv) = self.main_gate.invert_with_flag(ctx, norm)?;

        // 3. Compute (a - bi) * norm_inv
        let a_norm = self.main_gate.mul(ctx, a.0, &norm_inv)?;
        let b_norm = self.main_gate.mul(ctx, a.1, &norm_inv)?;

        // Negate b_norm
        let zero = self.main_gate.assign_value(ctx, Value::known(F::ZERO))?;
        let b_norm_neg = self.main_gate.sub(ctx, &zero, &b_norm)?;

        Ok((a_norm, b_norm_neg))
    }

    // Helper method to check if an Fq2 element is zero
    fn fq2_is_zero(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: (&AssignedValue<F>, &AssignedValue<F>),
    ) -> Result<AssignedValue<F>, Error> {
        let a_is_zero = self.main_gate.is_zero_term(ctx, a.0.clone())?;
        let b_is_zero = self.main_gate.is_zero_term(ctx, a.1.clone())?;
        self.main_gate.mul(ctx, &a_is_zero, &b_is_zero)
    }

    // Fq2 inversion with zero check
    pub fn fq2_inv_or_zero(
        &self,
        ctx: &mut RegionCtx<'_, F>,
        a: (&AssignedValue<F>, &AssignedValue<F>),
    ) -> Result<(AssignedValue<F>, AssignedValue<F>), Error> {
        let is_zero = self.fq2_is_zero(ctx, a)?;
        let inv = self.fq2_inv(ctx, a)?;

        // Conditionally select between (0, 0) and the computed inverse
        let zero = self.main_gate.assign_value(ctx, Value::known(F::ZERO))?;
        let x = self
            .main_gate
            .conditional_select(ctx, &zero, &inv.0, &is_zero)?;
        let y = self
            .main_gate
            .conditional_select(ctx, &zero, &inv.1, &is_zero)?;

        Ok((x, y))
    }
}

// Test module
#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::{
        ff::Field,
        halo2curves::bn256::{Fr, G1Affine},
        run_mock_prover_test,
        util::fe_to_fe_safe,
    };
    use halo2_proofs::{
        circuit::{Layouter, SimpleFloorPlanner},
        plonk::{Circuit, Column, ConstraintSystem, Error, Instance},
    };
    use rand_core::OsRng;

    #[derive(Clone, Debug)]
    struct TestG2CircuitConfig {
        config: MainGateConfig<5>,
        instance: Column<Instance>,
    }

    struct TestG2Circuit<
        C: CurveAffine<Base = F, ScalarExt = E>,
        F: PrimeFieldBits,
        E: PrimeFieldBits,
    > {
        a: G2Point<C, E>,
        b: G2Point<C, E>,
        lambda: E,
        test_case: usize, // 0: add, 1: scalar_mul
    }

    impl<C: CurveAffine<Base = F, ScalarExt = E>, F: PrimeFieldBits, E: PrimeFieldBits>
        TestG2Circuit<C, F, E>
    {
        fn new(a: G2Point<C, E>, b: G2Point<C, E>, lambda: E, test_case: usize) -> Self {
            Self {
                a,
                b,
                lambda,
                test_case,
            }
        }
    }

    impl<C: CurveAffine<Base = F, ScalarExt = E>, F: PrimeFieldBits, E: PrimeFieldBits> Circuit<F>
        for TestG2Circuit<C, F, E>
    {
        type Config = TestG2CircuitConfig;
        type FloorPlanner = SimpleFloorPlanner;

        fn without_witnesses(&self) -> Self {
            Self::new(
                G2Point::default(),
                G2Point::default(),
                C::ScalarExt::ZERO,
                0,
            )
        }

        fn configure(meta: &mut ConstraintSystem<C::Base>) -> Self::Config {
            let instance = meta.instance_column();
            meta.enable_equality(instance);
            let config = MainGate::<C::Base, 5>::configure(meta);
            Self::Config { config, instance }
        }

        fn synthesize(
            &self,
            config: TestG2CircuitConfig,
            mut layouter: impl Layouter<C::Base>,
        ) -> Result<(), Error> {
            let g2_chip = G2EccChip::<C, C::Base, 5>::new(config.config);

            // Perform the G2 operations and return the result (output)
            let output = layouter.assign_region(
                || "G2 operations",
                |region| {
                    let ctx = &mut RegionCtx::new(region, 0);

                    let ax0 = ctx.assign_advice(
                        || "a.x.c0",
                        g2_chip.main_gate.config().state[0],
                        Value::known(self.a.x.c0),
                    )?;
                    let ax1 = ctx.assign_advice(
                        || "a.x.c1",
                        g2_chip.main_gate.config().state[1],
                        Value::known(self.a.x.c1),
                    )?;
                    let ay0 = ctx.assign_advice(
                        || "a.y.c0",
                        g2_chip.main_gate.config().state[2],
                        Value::known(self.a.y.c0),
                    )?;
                    let ay1 = ctx.assign_advice(
                        || "a.y.c1",
                        g2_chip.main_gate.config().state[3],
                        Value::known(self.a.y.c1),
                    )?;
                    let a = AssignedG2Point {
                        x: (ax0, ax1),
                        y: (ay0, ay1),
                    };
                    ctx.next();

                    // Perform G2 addition or scalar multiplication depending on `test_case`
                    let result = if self.test_case == 0 {
                        let bx0 = ctx.assign_advice(
                            || "b.x.c0",
                            g2_chip.main_gate.config().state[0],
                            Value::known(self.b.x.c0),
                        )?;
                        let bx1 = ctx.assign_advice(
                            || "b.x.c1",
                            g2_chip.main_gate.config().state[1],
                            Value::known(self.b.x.c1),
                        )?;
                        let by0 = ctx.assign_advice(
                            || "b.y.c0",
                            g2_chip.main_gate.config().state[2],
                            Value::known(self.b.y.c0),
                        )?;
                        let by1 = ctx.assign_advice(
                            || "b.y.c1",
                            g2_chip.main_gate.config().state[3],
                            Value::known(self.b.y.c1),
                        )?;
                        let b = AssignedG2Point {
                            x: (bx0, bx1),
                            y: (by0, by1),
                        };
                        ctx.next();
                        g2_chip.add_g2(ctx, &a, &b)?
                    } else {
                        let lambda: C::Base = fe_to_fe_safe(&self.lambda).unwrap();
                        let bit_len =
                            NonZeroUsize::new(lambda.to_le_bits().len()).expect("Non Zero");
                        let lambda_assigned = ctx.assign_advice(
                            || "lambda",
                            g2_chip.main_gate.config().state[4],
                            Value::known(lambda),
                        )?;
                        ctx.next();

                        let bits =
                            g2_chip
                                .main_gate
                                .le_num_to_bits(ctx, lambda_assigned, bit_len)?;
                        g2_chip.scalar_mul(ctx, &a, &bits)?
                    };

                    // Return the result from the region
                    Ok(result)
                },
            )?;

            // Use a separate call to constrain the result to the instance columns
            layouter.constrain_instance(output.x.0.cell(), config.instance, 0)?;
            layouter.constrain_instance(output.x.1.cell(), config.instance, 1)?;
            layouter.constrain_instance(output.y.0.cell(), config.instance, 2)?;
            layouter.constrain_instance(output.y.1.cell(), config.instance, 3)?;

            Ok(())
        }
    }

    #[test]
    fn test_g2_add() {
        let K: u32 = 14;
        let p: G2Point<G1Affine, Fr> = G2Point::random_vartime();
        let q: G2Point<G1Affine, Fr> = G2Point::random_vartime();
        let r = p.add(&q);

        let circuit = TestG2Circuit::new(p, q, Fr::zero(), 0);
        let public_inputs = vec![vec![r.x.c0, r.x.c1, r.y.c0, r.y.c1]];
        run_mock_prover_test!(K, circuit, public_inputs);
    }

    #[test]
    fn test_g2_scalar_mul() {
        let K: u32 = 17;
        let p: G2Point<G1Affine, Fr> = G2Point::random_vartime();
        let q: G2Point<G1Affine, Fr> = G2Point::default(); // Point at infinity
        let lambda = Fr::random(&mut OsRng);
        let r = p.scalar_mul(&lambda);
        let circuit = TestG2Circuit::new(p, q, lambda, 1);
        let public_inputs = vec![vec![r.x.c0, r.x.c1, r.y.c0, r.y.c1]];
        run_mock_prover_test!(K, circuit, public_inputs);
    }

    #[test]
    fn test_impl() {
        let p = g2_generator::<G1Affine, Fr>();
        let q = p.double();
        let r = p.add(&q);

        assert_eq!(p.add(&p), q, "P + P should equal 2P");
        assert_eq!(p.add(&q), r, "P + 2P should equal 3P");

        let two_p = p.scalar_mul(&Fr::from(2u64));
        assert_eq!(two_p, q, "2 * P should equal 2P (double of P)");

        let three_p = p.scalar_mul(&Fr::from(3u64));
        assert_eq!(three_p, r, "3 * P should equal P + 2P");

        let left = p.add(&q).add(&r);
        let right = p.add(&q.add(&r));
        assert_eq!(left, right, "Point addition should be associative");

        assert_eq!(p.add(&q), q.add(&p), "Point addition should be commutative");

        let infinity = G2Point::default();
        assert_eq!(p.add(&infinity), p, "P + infinity should equal P");

        let neg_p = p.scalar_mul(&(-Fr::one()));
        assert_eq!(p.add(&neg_p), infinity, "P + (-P) should equal infinity");

        let a = Fr::from(3u64);
        let b = Fr::from(4u64);
        let left = p.scalar_mul(&a).add(&p.scalar_mul(&b));
        let right = p.scalar_mul(&(a + b));
        assert_eq!(left, right, "a*P + b*P should equal (a+b)*P");

        println!("All G2 operation tests passed!");
    }
}
