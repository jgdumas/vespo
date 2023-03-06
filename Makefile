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
OBJ=vespo_bench

all: ${OBJ}


# Generic compilation
LOADINGLIBS=${LIBS_NAME:%=-l%}
INCLUDES=${LIBS_DIR:%=-I%/include}
LIBFLAGS=${LIBS_DIR:%=-L%/lib}

CFLAGS += ${OPTFLAGS} ${VARFLAGS} ${INCLUDES}
CXXFLAGS += ${OPTFLAGS} ${VARFLAGS} ${INCLUDES}
LDFLAGS += ${LIBFLAGS}
LOADLIBES += ${LOADINGLIBS}
