#!/usr/bin/env bash

### Parse WMH mnc and extract their volumes

set -ex

ROOT=/ipl/ipl27/sfernandez/hvr_pet
DATA=${ROOT}/data/data_2024
LIST=${ROOT}/triad_2024.lst
OUTFILE=${ROOT}/data/derivatives/wmh_vols.csv

echo "PTID,SESSION,WMH" > $OUTFILE

mapfile -t IDS < $LIST

for id in ${IDS[@]}
do
	ptid=$(printf $id | cut -d, -f1)
	session=$(printf $id | cut -d, -f2)

	wmh_img=${DATA}/${ptid}/${session}/stx_${ptid}_${session}_WMH.mnc
	[ -f $wmh_img ] || continue

	wmh_vol=$(print_all_labels $wmh_img | awk '/Label: 1/ {print $3}')
	printf "%s,%s,%i\n" $ptid $session $wmh_vol >> $OUTFILE
done
