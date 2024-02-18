#!/bin/bash

set -euo pipefail
set -o errexit

# Store working directory so script can be run outside of dir.
wd=$(dirname $0)

if [[ $# -ne 4 ]]
then
    echo "usage: $0 input_folder output_folder path_and_name_of_the_profile number_of_threads" 1>&2
    exit 1
fi

if [ ! -d $1 ]
then
    echo "$1 is not a directory!" 1>&2
    exit 1
fi

if [ ! -d $2 ]
then
    echo "$2 is not a directory!" 1>&2
    exit 1
fi

if [ ! -s $3 ]
then
    echo "$3 file not found!" 1>&2
    exit 1
fi

input_dir=$1
output_dir=$2
hmm_profile_name=$3
threads=$4

# FEDOR: I usually use .fa extension
files=( $(find "$input_dir" -type f -name '*.fa' -print) )

p="$(basename "$hmm_profile_name")"
bp="${p%.*}"

for filename in "${files[@]}"
do
    n="$(basename "$filename")"
    bn="${n%.*}"

    echo "On $filename..."

#   HMMER analysis
    nhmmer --cpu $threads --notextw --noali --tblout $output_dir/nhmmer-$bp-vs-$bn-tbl.out -o /dev/null $hmm_profile_name $filename

#   Converting HMM table output to the BED format with filtering by threshold score to length
#   https://github.com/enigene/hmmertblout2bed
    awk -v th=0.7 -f $wd/hmmertblout2bed.awk $output_dir/nhmmer-$bp-vs-$bn-tbl.out > $output_dir/nhmmer-$bp-vs-$bn-tbl.bed

#   FEDOR: remove huge .out file as soon as it's used
    rm $output_dir/nhmmer-$bp-vs-$bn-tbl.out

#   Sorting by name and coordinates
#   recommended running sort with option --temporary-directory=/path/to/another/local/physical/disk
    sort -k 1.4,1 -k 2,2n $output_dir/nhmmer-$bp-vs-$bn-tbl.bed > $output_dir/_nhmmer-t0-$bn.bed

#   Filter by score in each region
    bedmap --max-element --fraction-either 0.1 $output_dir/_nhmmer-t0-$bn.bed > $output_dir/_nhmmer-t1-$bn.bed

#   Filter unique elements
    # FEDOR: don't skip SF monomers
    awk "{if(!(\$0 in a)){a[\$0]; print}}" $output_dir/_nhmmer-t1-$bn.bed > $output_dir/_nhmmer-t0-$bn.bed

#   FEDOR: addicional overlap filtering
    python3 $wd/overlap_filter.py $output_dir/_nhmmer-t0-$bn.bed > $output_dir/AS-SF-vs-$bn.bed

#   FEDOR: AS-strand annotation. "+" is blue, "-" is red
    awk -F $'\t' 'BEGIN {OFS = FS} {if ($6=="+") {$9="0,0,255"}; if ($6=="-") {$9="255,0,0"} print $0}' $output_dir/AS-SF-vs-$bn.bed > $output_dir/AS-strand-vs-$bn.bed

#   Delete temporary files
    rm $output_dir/_nhmmer-t1-$bn.bed
    rm $output_dir/_nhmmer-t0-$bn.bed
    rm $output_dir/nhmmer-$bp-vs-$bn-tbl.bed

done
