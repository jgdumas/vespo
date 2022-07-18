#!/bin/sh

# Example of RELIC presets

cmake -DWITH=ALL -DALLOC=AUTO -DWSIZE=64 -DRAND=UDEV -DSHLIB=OFF -DSTBIN=ON -DTIMER=CYCLE -DCHECK=off -DVERBS=off -DARITH=x64-asm-4l -DFP_PRIME=254 -DFP_METHD="INTEG;INTEG;INTEG;MONTY;LOWER;LOWER;SLIDE" -DCFLAGS="-Ofast -funroll-loops -fomit-frame-pointer -finline-small-functions -march=native -mtune=native" -DFP_PMERS=off -DFP_QNRES=on -DFPX_METHD="INTEG;INTEG;LAZYR" -DPP_METHD="LAZYR;OATEP" -DBN_PRECI=4096 -DCMAKE_BUILD_TYPE=Release -DBENCH=10 -DTESTS=10 -DCMAKE_INSTALL_PREFIX="/usr/local/soft/relic-0.6.0" $1
