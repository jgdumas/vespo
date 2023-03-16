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
	- Then run one small example and a benchmark with increasing degrees



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
	1. run `./Linear-FDT.sh [numthreads (default=32)]`
	2. run `./PoR-FDT.sh [numthreads (default=32)]`
	3. Parse outputs via `./parse_bench.sh bench_*`

- By hand
	- Using executable `./vespo_bench`
	- Usage: `vespo_bench [Vector dimension] [Paillier bitsize] [Verification iterations] [tasks (default=threads)]`
	- Example: `./vespo_bench 8192 2048 3`  



**References**:
- VESPo: Verified Evaluation of Secret Polynomials.   
  J-G. Dumas, A. Maignan, C. Pernet, D. S. Roche,   
  [ArXiv 2110.02022](https://arxiv.org/abs/2110.02022)
