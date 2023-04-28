--------------------------------------------------------------------------------
# VESPo: a C++ library for the Verified Evaluation of Secret Polynomials
--------------------------------------------------------------------------------

**Authors**:  Jean-Guillaume Dumas


**Requirements**:
- [RELIC](https://github.com/relic-toolkit/relic)
- [GMP](https://gmplib.org/), dev: headers & library



**Automatic linux install & run first benchmarks**:
- Fetch and run [auto-vm.run](https://raw.githubusercontent.com/jgdumas/vespo/main/auto-vm.run)

	- Requires a linux virtual machine or sudoer rights to install packages
	- Will install distribution packages: `wget`, `git`, `g++`, `cmake`, `libgmp-dev`.
	- Then clone and install `relic`, `libvespo`
	- Then run one small example and a large benchmark with increasing degrees



**Installation only**:
- Automatic (will install RELIC -- for now with x64-asm-4l arithmetic backend --)
	- Fetch and run [libvespo-auto-install.sh](https://raw.githubusercontent.com/jgdumas/vespo/main/libvespo-auto-install.sh)
		- Same requirements (git, cmake, c++ compiler, OpenMP, GMP --dev: headers & library--, ...) already installed.
		- Will clone and install `relic` (for now only with `x64-asm-4l` arithmetic backend)
		- Will clone and install `libvespo`

- By hand (supposing RELIC is already installed)
	1.  Clone the vespo directory and `cd vespo`
	2.  Set the RELIC/GMP path within the `Makefile`
	3.  Toggle compilation flags variants within the `Makefile`
	4.  Compile with `make`



**Benchmarking only**: `cd vespo`
- Automatic (two different range of benchmarks)
	1. run `./Linear-FDT.sh [numthreads (default=32)]` (Ref. Table 7).
	2. run `./PoR-FDT.sh [numthreads (default=32)]` (Ref. Table 10).
	3. Parse outputs via `./parse_bench.sh bench_*`

- By hand
	- Using executable `./vespo_bench`
	- Usage: `vespo_bench [Vector dimension] [Paillier bitsize] [Verification iterations] [tasks (default=threads)]`
	- Example: `./vespo_bench 8192 2048 3`

**Results**
- Parsed output
	1. ***Degree***: random polynomial degree
	2. ***Setup***: initialization time
	3. ***CStore***: client keys size
	4. ***CTime***: client verification time
	5. ***STime***: server computation time
	6. ***UTime***: client/server polynomial update time
	7. ***HTime***: reference Horner evaluation time

- Detailed Audit Benchmarks
	1. ***zeta***: encrypted evaluation by the server of P(r)
	2. ***xi***: proof by the server that the evaluation is correct
	3. ***C-gsum***: geometric sum checkpointing by the client
	4. ***C-powm***: pairings verification by the client
	5. ***C-H_dec***: deciphering the evaluation by the client
	6. ***horner***: reference Horner evaluation time
	

**References**:
- VESPo: Verified Evaluation of Secret Polynomials.  
  PETS 2023: 23rd Privacy Enhancing Technologies Symposium,  
  J-G. Dumas, A. Maignan, C. Pernet, D. S. Roche,
  [ArXiv 2110.02022](https://arxiv.org/abs/2110.02022)
