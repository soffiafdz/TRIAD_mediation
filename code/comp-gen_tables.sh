#!/usr/bin/env bash

## Shell script to generate a single LaTeX file compiling all tables
## and compile it

# HERE
HERE=/ipl/ipl27/sfernandez/hvr_pet
TABLESDIR=${HERE}/tables
OUTPUT=${HERE}/tables/tables.tex

## Create the master file with the preamble
cat << 'EOF' > $OUTPUT
\documentclass{article}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{pdflscape}
\usepackage{longtable}
\usepackage{multirow}
\begin{document}
EOF

# Loop through all table files (adjust the pattern as needed)
for file in ${TABLESDIR}/table-*.tex
do
	# Remove stub line
	sed -i -E '/c{2,}/ s/\|//g' $file
	# Add title
	if [[ $file =~ table-([0-9]) ]]
	then
		tablenum="${BASH_REMATCH[1]}"
		case "$tablenum" in
			1) title="Table 1: Demographic data (HC integrity)" ;;
			2) title="Table 2: Demographic data (Cognition)" ;;
			*) title="Table $tablenum" ;;
		esac
	fi

	# Insert table
	if [[ $file == *"_lscape.tex" ]]
	then
		echo "\begin{landscape}" >> $OUTPUT
		echo "\begin{center}" >> $OUTPUT
		echo "\Large\textbf{$title}" >> $OUTPUT
		echo "\end{center}" >> $OUTPUT
		echo "\input{$file}" >> $OUTPUT
		echo "\end{landscape}" >> $OUTPUT
	else
		echo "\section*{$title}" >> $OUTPUT
		echo "\input{$file}" >> $OUTPUT
	fi

	# Flush
	echo "\clearpage" >> $OUTPUT

	done

	# Close the document
	echo "\end{document}" >> $OUTPUT

	# Optionally compile the document
	pdflatex -output-directory=$TABLESDIR $OUTPUT
