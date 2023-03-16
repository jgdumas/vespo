#!/bin/bash -

# libvespo auto installer
# Copyright(c) 2023 Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>

### By default compile relic with 4 threads
MAKETHREADS=4
if [ ! -z "$1" ]
then
  MAKETHREADS=$1
fi


### VESPO git sources and commit ###
VESPO_GITDIR=https://github.com/jgdumas/vespo.git
VESPO_COMMIT=60791a2

### RELIC git sources and commit ###
RELIC_GITDIR=https://github.com/relic-toolkit/relic.git
RELIC_COMMIT=4e87ddc0

###  ###
LOGFILE=libvespo-auto-install.log


### die/cool ###
DONE="\033[0;36m done !\033[0m"
BEG="\033[1;32m * \033[0m"

die() {
    echo -ne "\n\033[1;31m * \033[0mfailed" ;
    if [[ -n $1 ]] ; then
	echo " ($1)"
    else
	echo "."
    fi
    exit -1 ;
}

cool() {
    echo -e $DONE
}


### Create local directory ###
echo "mkdir libvespo && cd libvespo"
mkdir libvespo && cd libvespo


### Extract RELIC sources ###
echo -en "${BEG}fetching RELIC..."| tee -a ${LOGFILE}
OK=0; git clone ${RELIC_GITDIR} 2>&1 >/dev/null && OK=1
[ "$OK" = "1" ] && cool | tee -a ${LOGFILE} || die
cd relic
OK=0; git reset --hard ${RELIC_COMMIT} 2>&1 >/dev/null && OK=1
[ "$OK" = "1" ] && cool | tee -a ../${LOGFILE} || die
cd ..


### Extract VESPO sources ###
echo -en "${BEG}fetching VESPO..."| tee -a ${LOGFILE}
OK=0; git clone ${VESPO_GITDIR} 2>&1 >/dev/null && OK=1
[ "$OK" = "1" ] && cool | tee -a ${LOGFILE} || die
cd vespo
OK=0; git reset --hard ${VESPO_COMMIT} 2>&1 >/dev/null && OK=1
[ "$OK" = "1" ] && cool | tee -a ../${LOGFILE} || die
cd ..


### Build RELIC ###
echo -e "${BEG}building RELIC ..."| tee -a ${LOGFILE}
cd relic
mkdir relic-target | tee -a ../${LOGFILE}
cd relic-target
cmake -DWITH=ALL -DALLOC=AUTO -DWSIZE=64 -DRAND=UDEV -DSHLIB=OFF -DSTBIN=ON -DTIMER=CYCLE -DCHECK=off -DVERBS=off -DARITH=x64-asm-4l -DFP_PRIME=254 -DFP_METHD="INTEG;INTEG;INTEG;MONTY;LOWER;LOWER;SLIDE" -DCFLAGS="-Ofast -funroll-loops -fomit-frame-pointer -finline-small-functions -march=native -mtune=native" -DFP_PMERS=off -DFP_QNRES=on -DFPX_METHD="INTEG;INTEG;LAZYR" -DPP_METHD="LAZYR;OATEP" -DBN_PRECI=4096 -DCMAKE_BUILD_TYPE=Release -DBENCH=10 -DTESTS=10 -DCMAKE_INSTALL_PREFIX="../" .. | tee -a ../../${LOGFILE}
make -j ${MAKETHREADS} | tee -a ../../${LOGFILE}
OK=0; make install | tee -a ../../${LOGFILE} || die && OK=1
[ "$OK" = "1" ] && cool | tee -a ../../${LOGFILE} || die
cd ../..

### Build VESPO ###
echo -e "${BEG}building VESPO ..."| tee -a ${LOGFILE}
cd vespo
OK=0; make LIBS_DIR=../relic | tee -a ../${LOGFILE} || die && OK=1
[ "$OK" = "1" ] && cool | tee -a ../${LOGFILE} || die
cd ..

### Done. ###
cd ..
