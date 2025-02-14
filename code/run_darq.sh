#! /usr/bin/env bash

## Apply DARQ models to ADNI subjects from list

set -ux

BASE_DIR=/ipl/ipl27/sfernandez/hvr_pet
MRI_DIR=${BASE_DIR}/data/adni
QC_DIR=${BASE_DIR}/lists/darq
DARQ=${BASE_DIR}/code/DARQ/python/aqc_apply.py

LIST=${BASE_DIR}/lists/adni_2024.csv

# Main
[ -d $QC_DIR ] || mkdir $QC_DIR

mapfile -t SCANS < $LIST
SCANS=( $(shuf -e "${SCANS[@]}") )

for scan in "${SCANS[@]}"
do
	adni=$(printf $scan | cut -d, -f1)
	magn=$(printf $scan | cut -d, -f2)
	sub=$(printf $scan | cut -d, -f3)
	date=$(printf $scan | cut -d, -f4)

	if [ $adni = adni1_2_go ]
	then
		mri_dir=${MRI_DIR}/${adni}/${magn}/${sub}/${date}
	else
		mri_dir=${MRI_DIR}/${adni}/${sub}/${date}
	fi

	mri_scan=${mri_dir}/*_t1.mnc
	[ -f $mri_scan ] || continue

	out_qc=${QC_DIR}/${sub}_${date}.txt
	[ -f $out_qc ] && continue

	printf "%s,%s,%s\n" \
		$sub $date \
		$(python3 $DARQ --volume $mri_scan --net r152 --raw) > $out_qc
done
