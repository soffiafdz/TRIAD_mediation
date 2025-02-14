#!/usr/bin/env bash

## Shell script for filtering the HCVC segmentation QC images
## of the ADNI subjects that are outliers (generally or by Dx)
## and copying them to a new folder to be parsed by Qrater

set -xu

HERE=/ipl/ipl27/sfernandez/hvr_pet
QCDATA=${HERE}/plots/qc_cnn/adni/hcvc/all
LIST=${HERE}/data/derivatives/adni_hcvc_outliers.csv

[ -f $LIST ] || Rscript ${HERE}/code/explore_segms_adni.R

if [ -f $LIST ]
then
	# Skip first line
	header=true
	while IFS= read line
	do
		$header && header=false && continue
		fname=$(printf $line | cut -d , -f 1)
		source_qc=$(printf "%s/%s.jpg" $QCDATA $fname)

		outlier=$(printf $line | cut -d , -f 2)
		if [ ! -z $outlier ]
		then
			outdir=${QCDATA/all/outlier}
			[ -d $outdir ] || mkdir $outdir
			outfile=$(printf "%s/%s_%s.jpg" $outdir $fname \
				$(echo $outlier | tr '[:upper:]' '[:lower:]' | sed "s/;/_/"))
			ln -f $source_qc $outfile
		fi

		outlier_dx=$(printf $line | cut -d , -f 3)
		if [ ! -z $outlier_dx ]
		then
			outdir=${QCDATA/all/outlier_dx}
			[ -d $outdir ] || mkdir $outdir
			outfile=$(printf "%s/%s_%s.jpg" $outdir $fname \
				$(echo $outlier | tr '[:upper:]' '[:lower:]' | sed "s/;/_/"))
			ln -f $source_qc $outfile
		fi
	done < $LIST
fi
set +xu
