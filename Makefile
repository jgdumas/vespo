# Location of the RELIC library
LIBS_DIR=/usr/local/soft/relic-0.6.0

# Used Libraries
LIBS_NAME=relic_s gmp rt

# Optimisation flags
OPTFLAGS += -Ofast -march=native

# Comment/Uncomment compilation variants
OPTFLAGS += -DVESPO_SUB_TIMINGS 		# detailed timings
OPTFLAGS += -DCLOCKTYPE=CLOCK_REALTIME -fopenmp	# using the parallel version
OPTFLAGS += -DVESPO_RELIC_LIMIT_MAX_ALLOC=4096	# limiting RELIC allocation
OPTFLAGS += -DVESPO_NOTSECURE=100		# only not benchmarking setup
#OPTFLAGS += -DVESPO_CHECKERS			# Debug: adding checkers
#OPTFLAGS += -DDEBUG				# Debug: general logs


# Benchmaring executable
OBJ=vespo_bench

all: ${OBJ}


# Generic compilation
LOADINGLIBS=${LIBS_NAME:%=-l%}
INCLUDES=${LIBS_DIR:%=-I%/include}
LIBFLAGS=${LIBS_DIR:%=-L%/lib}

CFLAGS += ${OPTFLAGS} ${INCLUDES}
CXXFLAGS += ${OPTFLAGS} ${INCLUDES}
LDFLAGS += ${LIBFLAGS}
LOADLIBES += ${LOADINGLIBS}
