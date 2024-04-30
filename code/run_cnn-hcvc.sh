#! /usr/bin/env bash

## Apply CNN ensemble models to TRIAD subjects from list
## Longitudinal list

## Need to load hvr_validation environment

TMPDIR=$(mktemp -d --tmpdir)
trap "rm -rf $TMPDIR" 0 1 2 15

set -ux

BASE_DIR=/ipl/ipl27/sfernandez/hvr_pet
MRI_DIR=${BASE_DIR}/data/data_2024
QC_DIR=${BASE_DIR}/plots/qc_cnn
OUT_DIR=${BASE_DIR}/data/derivatives
LIB_DIR=${BASE_DIR}/libraries/cnn-hcvc

function cnn_simple() {
	# Local variables
	local input=$1
	local output=$2

	local ref=${LIB_DIR}/ref_hcvc.mnc
	local model=${LIB_DIR}/ensemble_hcvc.pth

	local bname=$(basename $1 .mnc)

	local flip=$TMPDIR/flip.xfm
	[[ -e $flip ]] || param2xfm -scales -1 1 1 $flip

	# Left/Right
	local left=${TMPDIR}/${bname}_left.mnc
	local right=${TMPDIR}/${bname}_right.mnc

	[[ -e $right ]] || itk_resample $input --clobber --like $ref $right

	[[ -e $left ]] || itk_resample $input \
		--clobber --like $ref  \
		--transform $flip \
		$left

	# Apply on left
	date
	local left_seg=${TMPDIR}/${bname}_left_hcvc.mnc
	python ${BASE_DIR}/code/py_deep_seg/apply_multi_model.py \
		--stride 32 --patch 96 --cpu --crop 8 \
		$model $left $left_seg

	# Apply on right
	date
	local right_seg=${TMPDIR}/${bname}_right_hcvc.mnc
	python ${BASE_DIR}/code/py_deep_seg/apply_multi_model.py \
		--stride 32 --patch 96 --cpu --crop 8 \
		$model $right $right_seg

	date
	local left_rec=${TMPDIR}/${bname}_left_hcvc_recoded.mnc
	itk_resample --clobber --labels --byte --like $input \
		--transform $TMPDIR/flip.xfm \
		--lut-string "1 11; 2 12" \
		$left_seg $left_rec

	local right_rec=${TMPDIR}/${bname}_right_hcvc_recoded.mnc
	itk_resample --clobber --labels --byte --like $input \
		--lut-string "1 21; 2 22" \
		$right_seg $right_rec

	minccalc -q -clobber -labels -byte -express 'A[0]>0?A[0]:A[1]' \
		$left_rec $right_rec $output
}

function cnn_extra() {
	# Local variables
	local input=$1
	local output=$2

	local ref=${LIB_DIR}/ref_hcvc-ag.mnc
	local model=${LIB_DIR}/ensemble_hcvc-ag.pth

	local bname=$(basename $1 .mnc)

	local flip=$TMPDIR/flip.xfm
	[[ -e $flip ]] || param2xfm -scales -1 1 1 $flip

	# Left/Right
	local left=${TMPDIR}/${bname}_left.mnc
	local right=${TMPDIR}/${bname}_right.mnc

	[[ -e $right ]] || itk_resample $input --clobber --like $ref $right

	[[ -e $left ]] || itk_resample $input \
		--clobber --like $ref  \
		--transform $flip \
		$left

	# Apply on left
	date
	local left_seg=${TMPDIR}/${bname}_left_hcvc-ag.mnc
	python ${BASE_DIR}/code/py_deep_seg/apply_multi_model.py \
		--stride 32 --patch 96 --cpu --crop 8 \
		$model $left $left_seg

	# Apply on right
	date
	local right_seg=${TMPDIR}/${bname}_right_hcvc-ag.mnc
	python ${BASE_DIR}/code/py_deep_seg/apply_multi_model.py \
		--stride 32 --patch 96 --cpu --crop 8 \
		$model $right $right_seg

	date
	local left_rec=${TMPDIR}/${bname}_left_hcvc-ag_recoded.mnc
	itk_resample --clobber --labels --byte --like $input \
		--transform $TMPDIR/flip.xfm \
		--lut-string "1 111; 2 112; 3 113; 4 121; 5 122; 6 123; 7 130" \
		$left_seg $left_rec

	local right_rec=${TMPDIR}/${bname}_right_hcvc-ag_recoded.mnc
	itk_resample --clobber --labels --byte --like $input \
		--lut-string "1 211; 2 212; 3 213; 4 221; 5 222; 6 223; 7 230" \
		$right_seg $right_rec

	minccalc -q -clobber -labels -byte -express 'A[0]>0?A[0]:A[1]' \
		$left_rec $right_rec $output
}

function qc_plot() {
	# Local variables
	local img=$1

	local labels=$2
	local bname=$(basename $labels .mnc)

	local extra=$3
	if [[ $extra = true ]]
	then
		local plotter=${BASE_DIR}/code/qc/qc_plot.pl
		local outdir=${QC_DIR}/hcvc-ag/all
	else
		local plotter=${BASE_DIR}/code/qc/qc_plot_reduc.pl
		local outdir=${QC_DIR}/hcvc/all
	fi

	[ -d $outdir ] || mkdir -p $outdir
	local qcImg=${outdir}/${bname}.jpg
	[[ -e $qcImg ]] && rm $qcImg
	$plotter $img $labels -title $bname -clobber $qcImg
}

# Main
for in_mri in ${MRI_DIR}/*/*/*_t1_n.mnc
do
	bname=$(basename $in_mri .mnc)

	# HCVC CNN
	out_simple=${OUT_DIR}/hcvc/${bname}_hcvc.mnc
	[[ -d $(dirname $out_simple) ]] || mkd $(dirname $out_simple)
	if [[ ! -e $out_simple ]]
	then
		touch $out_simple
		cnn_simple $in_mri $out_simple
		qc_plot $in_mri $out_simple false
	fi

	# HCVC-AG CNN
	out_extra=${out_simple//hcvc/hcvc-ag}
	[[ -d $(dirname $out_extra) ]] || mkd $(dirname $out_extra)
	if [[ ! -e $out_extra ]]
	then
		touch $out_extra
		cnn_extra $in_mri $out_extra
		qc_plot $in_mri $out_extra true
	fi
done
