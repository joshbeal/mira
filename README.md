<h1 align="center">Mira: Efficient folding for pairing-based arguments</h1>

> [!NOTE]
> This repository is a fork of the Protostar implementation hosted at [https://github.com/snarkify/sirius](https://github.com/snarkify/sirius). 
> This is a research codebase that contributes an academic proof-of-concept implementation of the Mira proof system.

# Getting started

With Rust v1.80+ and Python v3.8+ installed:

```bash
python3 -m venv .env
source .env/bin/activate
pip install -r requirements.txt
```

# Run examples

For runnable examples, please check [examples](examples) folder.

```bash
# 're' is short for 'run example'

# Alias to run IVC for the SnarkStar protocol
cargo re-groth16-dev

# Alias to run IVC for the TensorStar protocol
cargo re-zkml-dev
```

# Time profiling 

Span lifetime tracking implemented, which allows you to understand in detail
how long a particular step takes to complete.

```bash
# 're' is short for 'run example'

# By default, it will output all spans with a lifetime greater than 1s
cargo re-groth16-dev | python3 .scripts/build_profiling.py

# It is possible to specify the bottom border of the output span
cargo re-groth16-dev | python3 .scripts/build_profiling.py --min-runtime 0.1s

# You can also store logs and process them at a later date
cargo re-groth16-dev > log; cat log | python3 .scripts/build_profiling.py 
```

# Memory profiling 
The [dhat](https://valgrind.org/docs/manual/dh-manual.html) utility and the [dhat-rs](https://github.com/nnethercote/dhat-rs) experimental crate are used

```bash
# Run dhat with default IVC protocol
cargo re-groth16-dhat
```
# Unit testing

```bash
cargo test --release
```

# Linting and formatting

```bash
cargo clippy -- -D warnings
ruff check --fix
ruff format
```
