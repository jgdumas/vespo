# Location of the RELIC library
LIBS_DIR=/usr/local/soft/relic-0.6.0

# Used Libraries
LIBS_NAME=relic_s gmp rt

# Optimisation flags
OPTFLAGS += -Ofast -march=native

# Comment/Uncomment compilation variants
VARFLAGS += -DVESPO_SUB_TIMINGS 		# detailed timings
VARFLAGS += -DCLOCKTYPE=CLOCK_REALTIME -fopenmp	# using the parallel version
VARFLAGS += -DVESPO_RELIC_LIMIT_MAX_ALLOC=4096	# limiting RELIC allocation
VARFLAGS += -DVESPO_NOTSECURE=100		# only not benchmarking setup
#VARFLAGS += -DVESPO_CHECKERS			# Debug: adding checkers
#VARFLAGS += -DDEBUG					# Debug: general logs


# Benchmaring executable
EXE=vespo_bench
SRC=${EXE:%=%.cpp}
DEP=vespo_library.h vespo_library.inl

all: ${EXE}


# Generic compilation
LOADINGLIBS=${LIBS_NAME:%=-l%}
INCLUDES=${LIBS_DIR:%=-I%/include}
LIBFLAGS=${LIBS_DIR:%=-L%/lib}

CXXFLAGS += ${OPTFLAGS} ${VARFLAGS} ${INCLUDES}
LDFLAGS += ${LIBFLAGS}
LOADLIBES += ${LOADINGLIBS}

# Prevents the stack to be executable
LDFLAGS += -z noexecstack

${EXE}: ${SRC} ${DEP}
	$(LINK.cpp) $< $(LOADLIBES) $(LDLIBS) -o $@
