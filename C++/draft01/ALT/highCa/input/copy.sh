#!/bin/bash

REDVOL=( 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 )
CONFIN=( 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 )
CAPNUM=( 00010 1 10 100 1000 )


#REDVOL=( 75 )
#CONFIN=( 90 )

# 1. create all folders
# 2. copy solution at high Ca from solving the 
#    full tenth-order system of ODEs
for i in "${REDVOL[@]}"
do
	mkdir v$i
	cd v$i
	for j in "${CONFIN[@]}"
	do
		mkdir conf$j
		cd conf$j
		for k in "${CAPNUM[@]}"
		do
			cp ../../../../allCa/output/v$i/conf$j/sln_v$i\_conf$j\_Ca$k\.dat ./
		#	sln_v$i\_conf$j\_CaInf.dat
		done
		cd ..
	done
	cd ..
done
