#!/bin/bash

NVAL=`printf "%7.4f" $(grep 'CONDUCTIVITY IS' normal.out | awk '{print $3}')`
SVAL=`printf "%7.4f" $(grep 'CONDUCTIVITY IS' split.out | tail -n3 | awk '{sum += $3} END {print sum}')`
echo "${NVAL} ${SVAL}"

if [ "${NVAL}" == "${SVAL}" ]; then
	echo "everything went well"
else
	echo "you suck"
    exit 1
fi
