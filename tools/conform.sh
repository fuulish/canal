#!/bin/bash

NVAL=`printf "%7.4f" $(grep 'CONDUCTIVITY IS' normal.out | awk '{print $3}')`
SVAL=`printf "%7.4f" $(grep 'CONDUCTIVITY IS' split.out | tail -n3 | awk '{sum += $3} END {print sum}')`
EVAL=`printf "%7.4f" $(grep 'actual conductivity' elmo.out | awk '{print $4}')`

echo "${NVAL} ${SVAL} ${EVAL}"

if [ "${NVAL}" == "${SVAL}" ]; then
    if [ "${EVAL}" == "${SVAL}" ]; then
	    echo "everything went well"
    else
    	echo "you suck"
        exit 1
    fi
else
	echo "you suck"
    exit 1
fi
