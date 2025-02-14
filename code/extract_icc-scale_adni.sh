#!/usr/bin/env bash

## Shell script for extracting the ICC and SCALE_factor
## and exporting them into a list
## ADNI data (2024)

set -xu

HERE=/ipl/ipl27/sfernandez/hvr_pet
DATA=${HERE}/data/adni
LIST=${HERE}/lists/adni_2024.csv
OUTFILE=${HERE}/data/derivatives/icc_scale_adni_old.csv

echo "PTID,VISIT,ICC,SCALEFACTOR" > $OUTFILE

mapfile -t IDS < $LIST

for id in ${IDS[@]}
do
	dset=$(printf $id | cut -d, -f1)
	mrit=$(printf $id | cut -d, -f2)
	subj=$(printf $id | cut -d, -f3)
	sess=$(printf $id | cut -d, -f4)

	if [ $dset = adni3 ]
	then
		path=${DATA}/${dset}/${subj}/${sess}
	else
		path=${DATA}/${dset}/${mrit}/${subj}/${sess}
	fi

	stx=${path}/stx2_${subj}_${sess}_t1.mnc
	mask=${path}/stx2_${subj}_${sess}_mask.mnc
	xfm=${path}/stx2_${subj}_${sess}_t1.xfm

	# SCALEFACTOR from STX2 xfm
	scale=$(xfm2param $xfm |
		awk '/-scale/{print $2*$3*$4}')

	# ICC (native space)
	icc=$(print_all_labels $mask |
		awk -v scale=$scale '{printf "%.10f", $NF / scale}')

	printf "%s,%s,%f,%f\n" \
		$subj $sess $icc $scale >> $OUTFILE
done < $LIST
