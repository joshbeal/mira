use std::{fmt, num::NonZeroUsize};

use ff::{Field, FromUniformBytes, PrimeField, PrimeFieldBits};
use halo2_proofs::{
    circuit::{AssignedCell, Layouter, Value},
    plonk::{Column, ConstraintSystem, Instance},
};
use halo2curves::CurveAffine;
use itertools::Itertools;
use serde::Serialize;

use crate::{
    constants::NUM_CHALLENGE_BITS,
    ivc::{
        fold_relaxed_plonk_instance_chip::{
            AssignedRelaxedPlonkInstance, FoldRelaxedPlonkInstanceChip, FoldResult,
        },
        StepCircuit, SynthesisError,
    },
    main_gate::{AdviceCyclicAssignor, MainGate, MainGateConfig, RegionCtx, WrapValue},
    plonk::{PlonkInstance, RelaxedPlonkInstance},
    poseidon::ROCircuitTrait,
};

use super::SimpleFloorPlanner;

#[derive(Serialize)]
#[serde(bound(serialize = "RO::Args: Serialize"))]
pub(crate) struct StepParams<F, RO>
where
    F: PrimeFieldBits + FromUniformBytes<64>,
    RO: ROCircuitTrait<F>,
{
    pub limb_width: NonZeroUsize,
    pub limbs_count: NonZeroUsize,
    pub ro_constant: RO::Args,
}

impl<F, RO> fmt::Debug for StepParams<F, RO>
where
    F: PrimeFieldBits + FromUniformBytes<64>,
    RO: ROCircuitTrait<F>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SynthesizeStepParams")
            .field("limb_width", &self.limb_width)
            .field("n_limbs", &self.limbs_count)
            .field("ro_constant", &self.ro_constant)
            .finish()
    }
}

pub(crate) struct StepInputs<'link, const ARITY: usize, C, RO>
where
    C::Base: ff::PrimeFieldBits + ff::FromUniformBytes<64>,
    C: CurveAffine,
    RO: ROCircuitTrait<C::Base>,
{
    pub step: C::Base,
    pub step_pp: &'link StepParams<C::Base, RO>,
    pub public_params_hash: C,

    pub z_0: [C::Base; ARITY],
    /// Output of previous step & input of current one
    pub z_i: [C::Base; ARITY],

    // TODO docs
    pub U: RelaxedPlonkInstance<C>,

    // TODO docs
    pub u: PlonkInstance<C>,

    // TODO docs
    pub cross_term_commits: Vec<C>,
}

pub struct StepConfig<const ARITY: usize, F: PrimeField, SP: StepCircuit<ARITY, F>, const T: usize>
{
    pub instance: Column<Instance>,
    pub step_config: SP::Config,
    pub main_gate_config: MainGateConfig<T>,
}

impl<const ARITY: usize, F: PrimeField + Clone, SP: StepCircuit<ARITY, F>, const T: usize> Clone
    for StepConfig<ARITY, F, SP, T>
where
    SP::Config: Clone,
{
    fn clone(&self) -> Self {
        Self {
            instance: self.instance,
            step_config: self.step_config.clone(),
            main_gate_config: self.main_gate_config.clone(),
        }
    }
}

impl<const ARITY: usize, F: PrimeField + fmt::Debug, SP: StepCircuit<ARITY, F>, const T: usize>
    fmt::Debug for StepConfig<ARITY, F, SP, T>
where
    SP::Config: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("StepConfig")
            .field("step_config", &self.step_config)
            .field("main_gate_config", &self.main_gate_config)
            .finish()
    }
}

/// Circuit that fold [`StepCircuit`] to fold represented the augmented function `F'` in the IVC scheme.
///
/// TODO #32 I will add details abote actions after implementation of trait to directly link
/// methods
///
/// # References
/// - For a detailed understanding of IVC and the context in which a struct
///   [`StepFoldingCircuit`] might be used, refer to the 'Section 5' of
///   [Nova Whitepaper](https://eprint.iacr.org/2023/969.pdf).
///   This trait is representation of `F'` function at 'Figure 4'
///       - `F'` is an augmented version of `F` designed to produce the IVC proofs `P_i` at each step.
///         It takes additional parameters such as \( vk, Ui, ui, (i, z0, zi), \omega_i, T \) and outputs \( x \).
///
/// - For `F'` please look at [`StepFoldingCircuit`]
pub(crate) struct StepFoldingCircuit<'link, const ARITY: usize, C, SC, RO, const T: usize>
where
    C: CurveAffine,
    C::Base: ff::PrimeFieldBits + ff::FromUniformBytes<64>,
    C::Scalar: ff::PrimeFieldBits + ff::FromUniformBytes<64>,
    SC: StepCircuit<ARITY, C::Base> + Sized,
    RO: ROCircuitTrait<C::Base>,
{
    pub step_circuit: &'link SC,
    pub input: StepInputs<'link, ARITY, C, RO>,
}

impl<'link, const ARITY: usize, C, SC, RO, const T: usize> StepCircuit<ARITY, C::Base>
    for StepFoldingCircuit<'link, ARITY, C, SC, RO, T>
where
    C: CurveAffine,
    C::Base: ff::PrimeFieldBits + ff::FromUniformBytes<64>,
    C::Scalar: ff::PrimeFieldBits + ff::FromUniformBytes<64>,
    SC: StepCircuit<ARITY, C::Base> + Sized,
    RO: ROCircuitTrait<C::Base, Config = MainGateConfig<T>>,
{
    type Config = StepConfig<ARITY, C::Base, SC, T>;
    type FloorPlanner = SimpleFloorPlanner;

    /// Configure the step circuit. This method initializes necessary
    /// fixed columns and advice columns, but does not create any instance
    /// columns.
    ///
    /// This setup is crucial for the functioning of the IVC-based system.
    fn configure(cs: &mut ConstraintSystem<C::Base>) -> Self::Config {
        let before = cs.num_instance_columns();

        let main_gate_config = MainGate::configure(cs);
        let step_config = <SC as StepCircuit<ARITY, C::Base>>::configure(cs);

        if before != cs.num_instance_columns() {
            panic!("You can't use instance column");
        }

        let instance = cs.instance_column();
        cs.enable_equality(instance);
        StepConfig {
            instance,
            step_config,
            main_gate_config,
        }
    }

    /// Sythesize the circuit for a computation step and return variable
    /// that corresponds to the output of the step z_{i+1}
    /// this method will be called when we synthesize the IVC_Circuit
    ///
    /// Return `z_out` result
    fn synthesize_step(
        &self,
        config: Self::Config,
        layouter: &mut impl Layouter<C::Base>,
        assigned_z_i: &[AssignedCell<C::Base, C::Base>; ARITY],
    ) -> Result<[AssignedCell<C::Base, C::Base>; ARITY], SynthesisError> {
        let assigned_z_0: [_; ARITY] = layouter.assign_region(
            || "assigned_z0_primary",
            |region| {
                let mut region = RegionCtx::new(region, 0);

                config
                    .main_gate_config
                    .advice_cycle_assigner()
                    .assign_all_advice(&mut region, || "z0_primary", self.input.z_0.iter().copied())
                    .map(|inp| inp.try_into().unwrap())
            },
        )?;

        let chip = FoldRelaxedPlonkInstanceChip::new(
            self.input.U.clone(),
            self.input.step_pp.limb_width,
            self.input.step_pp.limbs_count,
            config.main_gate_config.clone(),
        );

        let (w, r) = layouter.assign_region(
            || "assign witness",
            |region| {
                Ok(chip.assign_witness_with_challenge(
                    &mut RegionCtx::new(region, 0),
                    &self.input.public_params_hash,
                    &self.input.u,
                    &self.input.cross_term_commits,
                    RO::new(
                        config.main_gate_config.clone(),
                        self.input.step_pp.ro_constant.clone(),
                    ),
                )?)
            },
        )?;

        // Synthesize the circuit for the base case and get the new running instance
        let U_new_base = w.assigned_relaxed.clone();

        let (assigned_step, assigned_next_step) = layouter.assign_region(
            || "generate input",
            |region| {
                let mut region = RegionCtx::new(region, 0);
                let gate = MainGate::new(config.main_gate_config.clone());

                let assigned_step = region.assign_advice(
                    || "step",
                    config.main_gate_config.input,
                    Value::known(self.input.step),
                )?;

                let assigned_next_step =
                    gate.add_with_const(&mut region, &assigned_step, C::Base::ONE)?;

                Ok((assigned_step, assigned_next_step))
            },
        )?;

        layouter.assign_region(
            || "generate input hash",
            |region| {
                let mut ctx = RegionCtx::new(region, 0);

                let bits = RO::new(
                    config.main_gate_config.clone(),
                    self.input.step_pp.ro_constant.clone(),
                )
                .absorb_point(WrapValue::from_assigned_point(&w.public_params_hash))
                .absorb_base(WrapValue::Assigned(assigned_step.clone()))
                .absorb_iter(assigned_z_0.iter())
                .absorb_iter(assigned_z_i.iter().cloned())
                .absorb_iter(w.assigned_relaxed.iter_wrap_values())
                .squeeze_n_bits(&mut ctx, NUM_CHALLENGE_BITS)?;

                let gate = MainGate::new(config.main_gate_config.clone());
                let expected_X0 = gate.le_bits_to_num(&mut ctx, &bits)?;
                ctx.constrain_equal(expected_X0.cell(), w.input_instance[0].0.cell())?;

                Ok(())
            },
        )?;

        // Synthesize the circuit for the non-base case and get the new running
        // instance along with a boolean indicating if all checks have passed
        let FoldResult {
            assigned_input: assigned_input_witness,
            assigned_result_of_fold: U_new_non_base,
        } = layouter.assign_region(
            || "synthesize_step_non_base_case",
            |region| Ok(chip.fold(&mut RegionCtx::new(region, 0), w.clone(), r.clone())?),
        )?;

        let (assigned_new_U, assigned_input) = layouter.assign_region(
            || "generate input",
            |region| {
                let mut region = RegionCtx::new(region, 0);
                let gate = MainGate::new(config.main_gate_config.clone());

                let assigned_is_zero_step =
                    gate.is_zero_term(&mut region, assigned_step.clone())?;

                let new_U = AssignedRelaxedPlonkInstance::<C>::conditional_select(
                    &mut region,
                    &config.main_gate_config,
                    &U_new_non_base,
                    &U_new_base,
                    assigned_is_zero_step.clone(),
                )?;

                let assigned_input: [_; ARITY] = assigned_z_0
                    .iter()
                    .zip_eq(assigned_z_i.iter())
                    .map(|(z_0, z_i)| {
                        gate.conditional_select(&mut region, z_0, z_i, &assigned_is_zero_step)
                    })
                    .collect::<Result<Vec<_>, _>>()?
                    .try_into()
                    .unwrap();

                Ok((new_U, assigned_input))
            },
        )?;

        let z_output =
            self.step_circuit
                .synthesize_step(config.step_config, layouter, &assigned_input)?;

        let output_hash = layouter.assign_region(
            || "generate output hash",
            |region| {
                let mut ctx = RegionCtx::new(region, 0);

                let bits = RO::new(
                    config.main_gate_config.clone(),
                    self.input.step_pp.ro_constant.clone(),
                )
                .absorb_point(WrapValue::from_assigned_point(
                    &assigned_input_witness.public_params_hash,
                ))
                .absorb_base(WrapValue::Assigned(assigned_next_step.clone()))
                .absorb_iter(assigned_z_0.iter())
                .absorb_iter(z_output.iter().cloned())
                .absorb_iter(assigned_new_U.iter_wrap_values())
                .squeeze_n_bits(&mut ctx, NUM_CHALLENGE_BITS)?;

                MainGate::new(config.main_gate_config.clone()).le_bits_to_num(&mut ctx, &bits)
            },
        )?;

        layouter.constrain_instance(
            assigned_input_witness.input_instance[1].0.cell(),
            config.instance,
            0,
        )?;
        layouter.constrain_instance(output_hash.cell(), config.instance, 1)?;

        Ok(z_output)
    }
}

impl<'link, const ARITY: usize, C, SC, RO, const T: usize>
    StepFoldingCircuit<'link, ARITY, C, SC, RO, T>
where
    C: CurveAffine,
    C::Base: ff::PrimeFieldBits + ff::FromUniformBytes<64>,
    C::Scalar: ff::PrimeFieldBits + ff::FromUniformBytes<64>,
    SC: StepCircuit<ARITY, C::Base> + Sized,
    RO: ROCircuitTrait<C::Base, Config = MainGateConfig<T>>,
{
    pub fn synthesize(
        &self,
        config: StepConfig<ARITY, C::Base, SC, T>,
        layouter: &mut impl Layouter<C::Base>,
    ) -> Result<[C::Base; ARITY], SynthesisError> {
        let assigned_z_i: [_; ARITY] = layouter.assign_region(
            || "assigned_zi",
            |region| {
                let mut region = RegionCtx::new(region, 0);

                config
                    .main_gate_config
                    .advice_cycle_assigner()
                    .assign_all_advice(&mut region, || "z0_primary", self.input.z_0.iter().copied())
                    .map(|inp| inp.try_into().unwrap())
            },
        )?;

        let assigned_z_next = <Self as StepCircuit<ARITY, C::Base>>::synthesize_step(
            self,
            config,
            layouter,
            &assigned_z_i,
        )?;

        Ok(assigned_z_next.map(|cell| cell.value().unwrap().cloned().unwrap()))
    }
}