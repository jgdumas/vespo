--------------------------------------------------------------------------------
# VESPo: a C++ library for the Verified Evaluation of Secret Polynomials
--------------------------------------------------------------------------------

**Authors**:  Jean-Guillaume Dumas


**Requirements**:
- [RELIC](https://github.com/relic-toolkit/relic)
- [GMP](https://gmplib.org/)  


**Installation**:
- Automatic
	- Fetch and run [libvespo-auto-install.sh](https://raw.githubusercontent.com/jgdumas/vespo/main/libvespo-auto-install.sh)


- By hand (supposing RELIC is already installed)
	0.  Clone the vespo directory and `cd vespo`
	1.  Set the RELIC/GMP path within the Makefile
	2.  Toggle compilation flags variants within the Makefile
	3.  Compile with `make`


**Benchmarking**: `cd vespo`
- Automatic: two different range of benchmarks
	1. run `./Linear-FDT.sh [numthreads (default=32)]`
	2. run `./PoR-FDT.sh [numthreads (default=32)]`
	3. Parse outputs via `./parse_bench.sh bench_*`


- By hand
	- Using executable `./vespo_bench`
	- Usage: vespo_bench [Vector dimension] [Paillier bitsize] [Verification iterations] [tasks (default=threads)]
	- Example `./vespo_bench 8192 2048 3`  


**References**:
- VESPo: Verified Evaluation of Secret Polynomials.   
  J-G. Dumas, A. Maignan, C. Pernet, D. S. Roche,   
  [ArXiv 2110.02022](https://arxiv.org/abs/2110.02022)
