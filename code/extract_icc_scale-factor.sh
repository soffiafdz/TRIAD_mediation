#!/usr/bin/env bash

## Shell script for extracting the ICC and SCALE_factor
## and exporting them into a list

set -xu

HERE=/ipl/ipl27/sfernandez/hvr_pet
LIST=${HERE}/triad.lst
VOLUMES=${HERE}/data/derivatives/icc_scale.csv

echo "PTID,SCANDATE,ICC,SCALEFACTOR" > $VOLUMES

mapfile -t IDS < $LIST

for id in ${IDS[@]}
do
	ptid=$(printf $id | cut -d, -f1)
	date=$(printf $id | cut -d, -f2)

	stx=${HERE}/data/t1/stx_${ptid}_${date}_t1_n.mnc
	mask=${HERE}/data/masks/stx_${ptid}_${date}_dmask.mnc
	xfm=${HERE}/data/xfms/stx_${ptid}_${date}_t1_lin.xfm

	# SCALEFACTOR from STX2 xfm
	scale=$(xfm2param $xfm |
		awk '/-scale/{print $2*$3*$4}')

	# ICC (native space)
	icc=$(print_all_labels $mask |
		awk -v scale=$scale '{printf "%.10f", $NF / scale}')

	printf "%s,%s,%f,%f\n" \
		$ptid $date $icc $scale >> $VOLUMES
done < $LIST
