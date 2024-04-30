#!/usr/bin/env bash

## Shell script for extracting the ICC and SCALE_factor
## and exporting them into a list

set -xu

HERE=/ipl/ipl27/sfernandez/hvr_pet
DATA=${HERE}/data/data_2024
LIST=${HERE}/triad_2024.lst
VOLUMES=${HERE}/data/derivatives/icc_scale.csv

echo "PTID,VISIT,ICC,SCALEFACTOR" > $VOLUMES

mapfile -t IDS < $LIST

for id in ${IDS[@]}
do
	ptid=$(printf $id | cut -d, -f1)
	session=$(printf $id | cut -d, -f2)

	stx=${DATA}/${ptid}/${session}/stx_${ptid}_${session}_t1_n.mnc
	mask=${DATA}/${ptid}/${session}/stx_${ptid}_${session}_dmask.mnc
	xfm=${DATA}/${ptid}/${session}/stx_${ptid}_${session}_t1_lin.xfm

	# SCALEFACTOR from STX2 xfm
	scale=$(xfm2param $xfm |
		awk '/-scale/{print $2*$3*$4}')

	# ICC (native space)
	icc=$(print_all_labels $mask |
		awk -v scale=$scale '{printf "%.10f", $NF / scale}')

	printf "%s,%s,%f,%f\n" \
		$ptid $session $icc $scale >> $VOLUMES
done < $LIST
