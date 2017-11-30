#!/bin/bash

NVAL=`printf "%7.4f" $(grep 'CONDUCTIVITY IS' normal.out | awk '{print $3}')`
SVAL=`printf "%7.4f" $(grep 'CONDUCTIVITY IS' split.out | awk '{sum += $3} END {print sum}')`
PVAL=`printf "%7.4f" $(grep 'CONDUCTIVITY IS' spatial.out | awk '{sum += $3} END {print sum}')`
EVAL=`printf "%7.4f" $(grep 'actual conductivity' elmo.out | awk '{print $4}')`

echo "${NVAL} ${SVAL} ${EVAL} ${PVAL}"

if [ "${NVAL}" == "${SVAL}" ]; then
    if [ "${EVAL}" == "${SVAL}" ]; then
        if [ "${EVAL}" == "${PVAL}" ]; then
	        echo "everything went well"
        else
            echo "you suck"
            exit 1
        fi
    else
    	echo "you suck"
        exit 1
    fi
else
	echo "you suck"
    exit 1
fi
