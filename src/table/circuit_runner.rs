use halo2_proofs::plonk::{Circuit, ConstraintSystem, Error, FloorPlanner};
use tracing::*;

use super::{circuit_data::CircuitData, ConstraintSystemMetainfo, WitnessCollector};
use crate::{
    ff::PrimeField,
    plonk::{self, PlonkStructure},
    polynomial::sparse::SparseMatrix,
    util::batch_invert_assigned,
};

pub type Witness<F> = Vec<Vec<F>>;

// TODO(jbeal): Add correct group elements to circuit runner

#[derive(Debug, Clone)]
pub struct CircuitRunner<F: PrimeField, CT: Circuit<F>> {
    pub(crate) k: u32,
    pub(crate) circuit: CT,
    pub(crate) cs: ConstraintSystem<F>,
    pub(crate) config: CT::Config,
    pub(crate) instance: Vec<F>,
    pub(crate) num_g1: u32,
    pub(crate) num_g2: u32,
    pub(crate) gt_deg: u32,
    pub(crate) gt_cnt: u32,
}

impl<F: PrimeField, CT: Circuit<F>> CircuitRunner<F, CT> {
    pub fn new(
        k: u32,
        circuit: CT,
        instance: Vec<F>,
        num_g1: u32,
        num_g2: u32,
        gt_deg: u32,
        gt_cnt: u32,
    ) -> Self {
        let mut cs = ConstraintSystem::default();

        CircuitRunner {
            config: CT::configure(&mut cs),
            k,
            circuit,
            cs,
            instance,
            num_g1,
            num_g2,
            gt_deg,
            gt_cnt,
        }
    }

    #[instrument(name = "circuit_collect_plonk_struct", skip_all)]
    pub fn try_collect_plonk_structure(&self) -> Result<PlonkStructure<F>, Error> {
        debug!("start build metainfo");
        let ConstraintSystemMetainfo {
            num_challenges,
            round_sizes,
            gates,
            custom_gates_lookup_compressed,
            num_g1_elements,
            num_g2_elements,
            target_group_folding_degree,
            target_group_cross_terms,
            ..
        } = ConstraintSystemMetainfo::build(
            self.k as usize,
            self.num_g1 as usize,
            self.num_g2 as usize,
            self.gt_deg as usize,
            self.gt_cnt as usize,
            &self.cs,
        );
        debug!("meta info is ready");

        debug!("start preprocessing");
        let PreprocessingData {
            permutation_matrix,
            fixed_columns,
            selectors,
        } = self.try_collect_preprocessing()?;
        debug!("preprocessing is ready");

        Ok(PlonkStructure {
            k: self.k as usize,
            num_io: self.instance.len(),
            selectors,
            fixed_columns,
            num_advice_columns: self.cs.num_advice_columns(),
            num_challenges,
            round_sizes,
            custom_gates_lookup_compressed,
            gates,
            permutation_matrix,
            lookup_arguments: plonk::lookup::Arguments::compress_from(&self.cs),
            num_g1_elems: num_g1_elements,
            num_g2_elems: num_g2_elements,
            target_group_folding_degree,
            target_group_cross_terms,
        })
    }

    #[instrument(name = "circuit_collect_witness", skip_all)]
    pub fn try_collect_witness(&self) -> Result<Witness<F>, Error> {
        let mut witness = WitnessCollector {
            instance: self.instance.clone(),
            advice: vec![vec![F::ZERO.into(); 1 << self.k]; self.cs.num_advice_columns()],
        };

        CT::FloorPlanner::synthesize(&mut witness, &self.circuit, self.config.clone(), vec![])?;

        Ok(batch_invert_assigned(&witness.advice))
    }

    fn try_collect_preprocessing(&self) -> Result<PreprocessingData<F>, Error> {
        let nrow = 1 << self.k;

        let mut circuit_data = CircuitData {
            k: self.k,
            num_io: self.instance.len(),
            fixed: vec![vec![F::ZERO.into(); nrow]; self.cs.num_fixed_columns()],
            selector: vec![vec![false; nrow]; self.cs.num_selectors],
            permutation: plonk::permutation::Assembly::new(nrow, &self.cs.permutation),
        };

        CT::FloorPlanner::synthesize(
            &mut circuit_data,
            &self.circuit,
            self.config.clone(),
            vec![],
        )?;

        Ok(PreprocessingData {
            permutation_matrix: plonk::util::construct_permutation_matrix(
                self.k as usize,
                self.instance.len(),
                &self.cs,
                &circuit_data.permutation,
            ),
            fixed_columns: batch_invert_assigned(&circuit_data.fixed),
            selectors: circuit_data.selector,
        })
    }
}

struct PreprocessingData<F: PrimeField> {
    pub(crate) permutation_matrix: SparseMatrix<F>,
    pub(crate) fixed_columns: Vec<Vec<F>>,
    pub(crate) selectors: Vec<Vec<bool>>,
}
