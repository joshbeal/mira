use std::{collections::VecDeque, iter};

use rand::Rng;

use mira::{
    gadgets::merkle_tree_gadget::{
        chip::MerkleTreeUpdateChip,
        off_circuit::{NodeUpdate, Proof, Tree},
        *,
    },
    halo2_proofs::{circuit::*, plonk::*},
    halo2curves::ff::{FromUniformBytes, PrimeFieldBits},
    ivc::{StepCircuit, SynthesisError},
    main_gate::{MainGate, RegionCtx},
};

use tracing::info;

const INDEX_LIMIT: u32 = 1 << 31;
const ARITY: usize = 1;

type ProofBatch<F> = Box<[Proof<F>]>;

pub struct ProgramCounterUpdateCircuit<F>
where
    F: PrimeFieldBits + serde::Serialize + FromUniformBytes<64>,
{
    tree: Tree<F>,
    proofs_batches: VecDeque<ProofBatch<F>>,
    batch_size: usize,
}

impl<F> Default for ProgramCounterUpdateCircuit<F>
where
    F: PrimeFieldBits + serde::Serialize + FromUniformBytes<64>,
{
    fn default() -> Self {
        Self {
            tree: Default::default(),
            proofs_batches: VecDeque::new(),
            batch_size: 1,
        }
    }
}

impl<F> ProgramCounterUpdateCircuit<F>
where
    F: PrimeFieldBits + serde::Serialize + FromUniformBytes<64>,
{
    pub fn new(batch_size: usize) -> Self {
        Self {
            batch_size,
            ..Default::default()
        }
    }

    pub fn new_with_random_updates(
        rng: &mut impl Rng,
        batch_size: usize,
        batches_count: usize,
    ) -> Self {
        let mut self_ = Self::new(batch_size);

        for _ in 0..=batches_count {
            self_.random_update_leaves(rng);
        }

        self_
    }

    pub fn new_with_updates(
        updates: &Vec<Vec<(u32, F)>>,
        batch_size: usize,
        batches_count: usize,
    ) -> Self {
        let mut self_ = Self::new(batch_size);

        for i in 0..=batches_count {
            self_.update_leaves(updates[i].iter().cloned());
        }

        self_
    }

    pub fn pop_front_proof_batch(&mut self) -> bool {
        self.proofs_batches.pop_front().is_some()
    }

    pub fn front_proof_batch(&self) -> &ProofBatch<F> {
        self.proofs_batches.front().as_ref().unwrap()
    }

    pub fn random_update_leaves(&mut self, mut rng: &mut impl Rng) {
        // 'allow' is necessary, because otherwise the closure captures rnd and we have to copy it
        #[allow(clippy::needless_borrows_for_generic_args)]
        self.update_leaves(iter::repeat_with(move || {
            (rng.gen::<u32>() % INDEX_LIMIT, F::random(&mut rng))
        }));
    }

    pub fn update_leaves(&mut self, update: impl Iterator<Item = (u32, F)>) -> (F, F) {
        let proofs: Box<[_]> = update
            .take(self.batch_size)
            .map(|(index, data)| self.tree.update_leaf(index, data))
            .collect();

        if proofs.is_empty() {
            panic!("No updates provided");
        }

        let old = proofs.first().unwrap().root().old;
        let new = proofs.last().unwrap().root().new;

        self.proofs_batches.push_back(proofs);

        (old, new)
    }
}

impl<F> StepCircuit<1, F> for ProgramCounterUpdateCircuit<F>
where
    F: PrimeFieldBits + serde::Serialize + FromUniformBytes<64>,
{
    type Config = MainGateConfig;
    fn configure(cs: &mut ConstraintSystem<F>) -> Self::Config {
        MainGate::configure(cs)
    }

    fn process_step(
        &self,
        _z_i: &[F; ARITY],
        _k_table_size: u32,
    ) -> Result<[F; ARITY], SynthesisError> {
        Ok([self.front_proof_batch().last().as_ref().unwrap().root().new])
    }

    fn synthesize_step(
        &self,
        config: Self::Config,
        layouter: &mut impl Layouter<F>,
        z_i: &[AssignedCell<F, F>; 1],
    ) -> Result<[AssignedCell<F, F>; 1], SynthesisError> {
        layouter
            .assign_region(
                || "",
                |region| {
                    let mut region = RegionCtx::new(region, 0);

                    let mut prev = z_i[0].clone();
                    for proof in self.front_proof_batch().iter() {
                        let NodeUpdate { old, new, .. } = MerkleTreeUpdateChip::new(proof.clone())
                            .prove_next_update(&mut region, config.clone())?;

                        region.constrain_equal(prev.cell(), old.cell())?;
                        prev = new;
                    }
                    info!("offset = {}", region.offset);

                    Ok([prev])
                },
            )
            .map_err(SynthesisError::Halo2)
    }
}

impl<F> Circuit<F> for ProgramCounterUpdateCircuit<F>
where
    F: PrimeFieldBits + serde::Serialize + FromUniformBytes<64>,
{
    type Config = <Self as StepCircuit<1, F>>::Config;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        todo!()
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> Self::Config {
        <Self as StepCircuit<1, F>>::configure(meta)
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
        layouter.assign_region(
            || "",
            |region| {
                let mut region = RegionCtx::new(region, 0);

                region.next();

                let mut prev = Option::<AssignedCell<F, F>>::None;
                for proof in self.front_proof_batch().iter() {
                    let NodeUpdate { old, new, .. } = MerkleTreeUpdateChip::new(proof.clone())
                        .prove_next_update(&mut region, config.clone())?;

                    if let Some(prev) = prev {
                        region.constrain_equal(prev.cell(), old.cell())?;
                    }
                    prev = Some(new);
                }

                Ok([prev])
            },
        )?;

        Ok(())
    }
}
