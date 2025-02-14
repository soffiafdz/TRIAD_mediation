#!/usr/bin/env bash

### Parse out_20240408 in ipl29 to link stx files for cube experiment
HERE=/ipl/ipl27/sfernandez/hvr_pet
ADNI3=${HERE}/data/adni/adni3
ADNI12GO=${HERE}/data/adni/adni1_2_go
SOURCEDIR=/ipl/ipl29/ADNI
LIST_NOTDIR=${HERE}/lists/missing_dir.txt
LIST_NOTXFM=${HERE}/lists/missing_xfm.txt

# Reinit lists
> $LIST_NOTDIR
> $LIST_NOTXFM

# Loop through local directories for ADNI12GO
fd "xfm" "$ADNI12GO" |
	# There are 9166 subdirectories
	pv -l -s 9166 |
	while read -r xfm_file
	do
		destdir=$(dirname $xfm_file)
		tesla=$(echo $destdir | cut -d / -f 9)
		sub=$(echo $destdir | cut -d / -f 10)
		date=$(echo $destdir | cut -d / -f 11)

		# Change to Uppercase
		[ $tesla = 3t ] && Tesla=3T || Tesla=15T

		# sourcedir
		printf -v sourcedir "%s/ADNI_1_2_GO/out_20240408/%s/%s/%s/stx" \
			$SOURCEDIR $Tesla $sub $date

		[ ! -d $sourcedir ] && echo $sourcedir >> $LIST_NOTDIR && continue

		# stx_xfm
		printf -v stx_xfm "%s/stx_%s_%s_t1.xfm" $sourcedir $sub $date

		[ ! -f $stx_xfm ] && echo $stx_xfm >> $LIST_NOTXFM && continue

		# Link
		ln -svf $stx_xfm $destdir
	done

# Loop through local directories for ADNI3
fd "xfm" "$ADNI3" |
	# There are 2148 subdirectories
	pv -l -s 2148 |
	while read -r xfm_file
	do
		destdir=$(dirname $xfm_file)
		sub=$(echo $destdir | cut -d / -f 9)
		date=$(echo $destdir | cut -d / -f 10)
		# sourcedir
		printf -v sourcedir "%s/ADNI3/out_20240408/%s/%s/stx" \
			$SOURCEDIR $sub $date
		[ ! -d $sourcedir ] && echo $sourcedir >> $LIST_NOTDIR && continue

		# stx_xfm
		printf -v stx_xfm "%s/stx_%s_%s_t1.xfm" $sourcedir $sub $date
		[ ! -f $stx_xfm ] && echo $stx_xfm >> $LIST_NOTXFM && continue

		# Link
		ln -svf $stx_xfm $destdir
	done
