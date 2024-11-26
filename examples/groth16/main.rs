#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::{
    fs::{self, File},
    path::{Path, PathBuf},
};

use clap::Parser;
use git2::Repository;
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt::format::FmtSpan, EnvFilter};

pub mod benchmark;
pub mod circuit;
pub mod conversion;
pub mod util;

mod mira_mod {
    use std::{io, num::NonZeroUsize, path::Path};

    use ark_bn254::{Bn254, Fr};
    use ark_ff::One;
    use ark_groth16::{prepare_verifying_key, Groth16};
    use ark_std::vec::Vec;
    use rand_core::SeedableRng;

    use crate::benchmark::Benchmark;
    use crate::conversion::{ark_to_ff_field, ark_to_halo2_g2_affine, ark_to_halo2_group};
    use crate::util::{Bn256Proof, Bn256VerifyingKey};

    use halo2_proofs::halo2curves::{bn256, grumpkin, CurveAffine};
    use mira::{
        commitment::CommitmentKey,
        ff::{Field, FromUniformBytes, PrimeField},
        gadgets::merkle_tree_gadget::{off_circuit::Tree, *},
        group::{prime::PrimeCurve, Group},
        ivc::{step_circuit::trivial, CircuitPublicParamsInput, PublicParams, IVC},
    };
    use tracing::info_span;

    use crate::circuit::Groth16AggregateCircuit;

    const ARITY: usize = 1;

    fn get_circuit_table_size(batch_size: usize) -> (usize, usize) {
        match batch_size {
            0 => (21, 21),
            1 => (19, 19),
            2 => (20, 20),
            4 => (21, 21),
            8 => (22, 22),
            16 => (23, 23),
            32 => (24, 24),
            _ => panic!(
                "Invalid batch size: {}. Supported batch sizes are 1, 2, 4, 8, 16, 32.",
                batch_size
            ),
        }
    }

    fn get_commitment_key_size(batch_size: usize) -> (usize, usize) {
        match batch_size {
            0 => (25, 24),
            1 => (23, 24),
            2 => (24, 24),
            4 => (25, 24),
            8 => (26, 25),
            16 => (27, 26),
            32 => (28, 27),
            _ => panic!(
                "Invalid batch size: {}. Supported batch sizes are 1, 2, 4, 8, 16, 32.",
                batch_size
            ),
        }
    }

    const LIMBS_COUNT_LIMIT: NonZeroUsize = unsafe { NonZeroUsize::new_unchecked(10) };
    const LIMB_WIDTH: NonZeroUsize = unsafe { NonZeroUsize::new_unchecked(32) };

    type C1 = bn256::G1;
    type C2 = grumpkin::G1;

    pub type C1Scalar = <C1 as Group>::Scalar;

    type C1Affine = <C1 as PrimeCurve>::Affine;
    type C2Affine = <C2 as PrimeCurve>::Affine;
    type C2Scalar = <C2 as Group>::Scalar;

    type RandomOracle = mira::poseidon::PoseidonRO<T, RATE>;
    type RandomOracleConstant<F> = <RandomOracle as ROPair<F>>::Args;

    fn get_or_create_commitment_key<C: CurveAffine>(
        k: usize,
        label: &'static str,
    ) -> io::Result<CommitmentKey<C>> {
        const FOLDER: &str = ".cache/examples";

        unsafe { CommitmentKey::load_or_setup_cache(Path::new(FOLDER), label, k) }
    }

    fn generate_groth16_aggregation_steps<F: PrimeField + FromUniformBytes<64>>(
        nproofs: usize,
        batch_size: usize,
    ) -> (
        Vec<Vec<Vec<F>>>,
        Vec<Vec<Bn256Proof>>,
        Vec<Vec<Bn256VerifyingKey>>,
    ) {
        let num_constraints = 1000;
        let mut rng = rand_chacha::ChaChaRng::seed_from_u64(1u64);
        let params = {
            let c = Benchmark::<Fr>::new(num_constraints);
            Groth16::<Bn254>::generate_random_parameters_with_reduction(c, &mut rng).unwrap()
        };
        // prepare the verification key
        let pvk = prepare_verifying_key(&params.vk);
        // create all the proofs
        let proofs = (0..nproofs)
            .map(|_| {
                let c = Benchmark::new(num_constraints);
                Groth16::<Bn254>::create_random_proof_with_reduction(c, &params, &mut rng)
                    .expect("Proof creation failed")
            })
            .collect::<Vec<_>>();
        // verify we can at least verify one
        let inputs: Vec<_> = [Fr::one(); 2].to_vec();
        let r = Groth16::<Bn254>::verify_proof(&pvk, &proofs[1], &inputs).unwrap();
        assert!(r, "Proof verification failed");

        let converted_inputs: Vec<_> = (0..nproofs)
            .map(|_| {
                inputs
                    .iter()
                    .map(|inp| ark_to_ff_field::<Fr, F, 64>(*inp).unwrap())
                    .collect()
            })
            .collect::<Vec<_>>();

        let converted_proofs: Vec<_> = proofs
            .iter()
            .map(|inp| Bn256Proof {
                a: ark_to_halo2_group::<ark_bn254::g1::Config, bn256::G1Affine, 64>(inp.a).unwrap(),
                b: ark_to_halo2_g2_affine(inp.b).unwrap(),
                c: ark_to_halo2_group::<ark_bn254::g1::Config, bn256::G1Affine, 64>(inp.c).unwrap(),
            })
            .collect();

        let converted_vks = (0..nproofs)
            .map(|_| Bn256VerifyingKey {
                alpha_g1: ark_to_halo2_group::<ark_bn254::g1::Config, bn256::G1Affine, 64>(
                    params.vk.alpha_g1,
                )
                .unwrap(),
                beta_g2: ark_to_halo2_g2_affine(params.vk.beta_g2).unwrap(),
                gamma_g2: ark_to_halo2_g2_affine(params.vk.gamma_g2).unwrap(),
                delta_g2: ark_to_halo2_g2_affine(params.vk.delta_g2).unwrap(),
                gamma_abc_g1: params
                    .vk
                    .gamma_abc_g1
                    .iter()
                    .map(|x| {
                        ark_to_halo2_group::<ark_bn254::g1::Config, bn256::G1Affine, 64>(*x)
                            .unwrap()
                    })
                    .collect(),
            })
            .collect::<Vec<_>>();

        let converted_inputs: Vec<Vec<Vec<F>>> = converted_inputs
            .chunks(batch_size)
            .map(|chunk| chunk.to_vec())
            .collect();

        let converted_proofs: Vec<Vec<Bn256Proof>> = converted_proofs
            .chunks(batch_size)
            .map(|chunk| chunk.to_vec())
            .collect();

        let converted_vks: Vec<Vec<Bn256VerifyingKey>> = converted_vks
            .chunks(batch_size)
            .map(|chunk| chunk.to_vec())
            .collect();

        (converted_inputs, converted_proofs, converted_vks)
    }

    // Assumes that all input vectors have length 2
    fn add_indices_to_elements<F: PrimeField>(v: Vec<Vec<Vec<F>>>) -> Vec<Vec<(u32, F)>> {
        v.into_iter()
            .enumerate()
            .flat_map(|(batch_idx, batch)| {
                let batch_len = batch.len();
                batch.into_iter().enumerate().map(move |(proof_idx, vec)| {
                    vec.into_iter()
                        .enumerate()
                        .map(|(j, element)| {
                            (
                                (batch_idx * batch_len * 2 + proof_idx * 2 + j) as u32,
                                element,
                            )
                        })
                        .collect::<Vec<(u32, F)>>()
                })
            })
            .collect()
    }

    pub fn run(repeat_count: usize, batch_size: usize, baseline: bool) {
        let _span = info_span!("groth16_example").entered();
        let prepare_span = info_span!("prepare").entered();

        let (inps, _, _) =
            generate_groth16_aggregation_steps((repeat_count + 1) * batch_size, batch_size);
        let updates = add_indices_to_elements(inps);

        let size_param = if baseline { 0 } else { batch_size };
        let (circuit_table_size1, circuit_table_size2) = get_circuit_table_size(size_param);
        let (commitment_key_size1, commitment_key_size2) = get_commitment_key_size(size_param);

        let mut sc1 = Groth16AggregateCircuit::new_with_updates(&updates, batch_size, repeat_count);

        let sc2 = trivial::Circuit::<ARITY, _>::default();

        let primary_spec = RandomOracleConstant::<C1Scalar>::new(10, 10);
        let secondary_spec = RandomOracleConstant::<C2Scalar>::new(10, 10);

        let primary_commitment_key =
            get_or_create_commitment_key::<bn256::G1Affine>(commitment_key_size1, "bn256")
                .expect("Failed to get secondary key");
        let secondary_commitment_key =
            get_or_create_commitment_key::<grumpkin::G1Affine>(commitment_key_size2, "grumpkin")
                .expect("Failed to get primary key");

        let pp = PublicParams::<
            '_,
            ARITY,
            ARITY,
            T,
            C1Affine,
            C2Affine,
            Groth16AggregateCircuit<_>,
            trivial::Circuit<ARITY, _>,
            RandomOracle,
            RandomOracle,
        >::new(
            CircuitPublicParamsInput::new(
                circuit_table_size1 as u32,
                0,
                0,
                0,
                0,
                &primary_commitment_key,
                primary_spec.clone(),
                &sc1,
            ),
            CircuitPublicParamsInput::new(
                circuit_table_size2 as u32,
                2 * batch_size as u32,
                1 * batch_size as u32,
                2,
                2 * batch_size as u32,
                &secondary_commitment_key,
                secondary_spec.clone(),
                &sc2,
            ),
            LIMB_WIDTH,
            LIMBS_COUNT_LIMIT,
        )
        .unwrap();

        prepare_span.exit();

        let mut ivc = IVC::new(
            &pp,
            &sc1,
            [*Tree::default().get_root()],
            &sc2,
            [C2Scalar::ZERO],
            false,
        )
        .unwrap();

        for _ in 0..repeat_count {
            sc1.pop_front_proof_batch();
            ivc.fold_step(&pp, &sc1, &sc2).unwrap();
        }

        ivc.verify(&pp).unwrap();
    }
}

#[derive(clap::Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(long, default_value_t = 1, global = true)]
    repeat_count: usize,
    #[arg(long, default_value_t = 1, global = true)]
    batch_size: usize,
    #[arg(long, default_value_t = false, global = true)]
    json_logs: bool,
    #[arg(long, default_value_t = false, global = true)]
    clean_cache: bool,
    /// Push all logs into file, with name builded from params
    #[arg(long, default_value_t = false, global = true)]
    file_logs: bool,
    #[arg(long, default_value_t = false, global = true)]
    baseline: bool,
}

pub fn build_log_folder() -> PathBuf {
    const LOGS_SUBFOLDER: &str = ".logs";

    let Ok(repo) = Repository::discover(".") else {
        return Path::new(LOGS_SUBFOLDER).to_path_buf();
    };

    // Get the current branch name
    let branch_name = repo
        .head()
        .ok()
        .and_then(|head| head.shorthand().map(String::from))
        .unwrap_or_else(|| "unknown".to_string());

    let branch_log_dir = repo
        .workdir()
        .unwrap()
        .join(LOGS_SUBFOLDER)
        .join(branch_name);

    fs::create_dir_all(&branch_log_dir)
        .unwrap_or_else(|err| panic!("Failed to create log directory {branch_log_dir:?}: {err:?}"));

    branch_log_dir
}

impl Args {
    fn build_log_filename(&self) -> Option<PathBuf> {
        if !self.file_logs {
            return None;
        }

        let Args {
            repeat_count,
            batch_size,
            ..
        } = &self;

        Some(build_log_folder().join(format!(
            "mira_groth16-1_trivial-1_{repeat_count}_{batch_size}.log"
        )))

        // Some(build_log_folder().join(format!(
        //     "mira_groth16-repeat-1_trivial-1_{repeat_count}_{batch_size}.log"
        // )))
    }

    fn init_logger(&self) {
        let mut builder = tracing_subscriber::fmt()
            // Adds events to track the entry and exit of the span, which are used to build
            // time-profiling
            .with_span_events(FmtSpan::ENTER | FmtSpan::CLOSE)
            // Changes the default level to INFO
            .with_env_filter(
                EnvFilter::builder()
                    .with_default_directive(LevelFilter::INFO.into())
                    .from_env_lossy(),
            );

        if let Some(log_filename) = self.build_log_filename() {
            let file = File::create(&log_filename).expect("Unable to create log file");

            builder = builder.with_ansi(false);

            if self.json_logs {
                builder.json().with_writer(file).init();
            } else {
                builder.with_writer(file).init();
            }
            println!(
                "logs will be written to: {}",
                log_filename.to_string_lossy()
            );
        } else if self.json_logs {
            builder.json().init();
        } else {
            builder.init();
        }

        info!("start with args: {self:?}");
    }
}

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let args = Args::parse();
    args.init_logger();

    mira_mod::run(args.repeat_count, args.batch_size, args.baseline);
}
