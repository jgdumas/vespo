--------------------------------------------------------------------------------
# VESPo: a C++ library for the Verified Evaluation of Secret Polynomials
--------------------------------------------------------------------------------


**Authors**:  Jean-Guillaume Dumas


*   Requirements:

    - [RELIC](https://github.com/relic-toolkit/relic)
    - [GMP](https://gmplib.org/)


*   Installation and running:

    1.  Set the RELIC/GMP path within the Makefile

    2.  Toggle compilation flags variants within the Makefile

    3.  Compile with
        `make`

*   Benchmarking

    - Using executable
	`vespo_bench`
    - Usage: vespo_bench <Vector dimension> <Paillier bitsize> <Verification iterations> [tasks (default=threads)]
    - Example
	`./vespo_bench 8192 2048 3`
