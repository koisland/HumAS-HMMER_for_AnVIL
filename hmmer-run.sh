#!/bin/bash

set -o pipefail
set -o errexit


# Store working directory so script can be run outside of dir.
wd=$(dirname $0)

if [[ $# -ne 3 ]]
then
    echo "usage: $0 input_folder path_and_name_of_the_profile number_of_threads" 1>&2
    exit 1
fi

if [ ! -d $1 ]
then
    echo "$1 is not a directory!" 1>&2
    exit 1
fi

if [ ! -s $2 ]
then
    echo "$2 file not found!" 1>&2
    exit 1
fi

dir=$1
hmm_profile_name=$2
threads=$3

# FEDOR: I usually use .fa extension
files=( $(find "$dir" -type f -name '*.fa' -print) )

p="$(basename "$hmm_profile_name")"
bp="${p%.*}"

for filename in "${files[@]}"
do
    n="$(basename "$filename")"
    bn="${n%.*}"

    echo "$filename"

#   HMMER analysis
    nhmmer --cpu $threads --notextw --noali --tblout $wd/nhmmer-$bp-vs-$bn-tbl.out -o /dev/null $hmm_profile_name $filename

#   Converting HMM table output to the BED format with filtering by threshold score to length
#   https://github.com/enigene/hmmertblout2bed
    awk -v th=0.7 -f $wd/hmmertblout2bed.awk $wd/nhmmer-$bp-vs-$bn-tbl.out > $wd/nhmmer-$bp-vs-$bn-tbl.bed

#   FEDOR: remove huge .out file as soon as it's used
    rm $wd/nhmmer-$bp-vs-$bn-tbl.out

#   Sorting by name and coordinates
#   recommended running sort with option --temporary-directory=/path/to/another/local/physical/disk
    sort -k 1.4,1 -k 2,2n $wd/nhmmer-$bp-vs-$bn-tbl.bed > $wd/_nhmmer-t0-$bn.bed

#   Filter by score in each region
    bedmap --max-element --fraction-either 0.1 $wd/_nhmmer-t0-$bn.bed > $wd/_nhmmer-t1-$bn.bed

#   Filter unique elements
    # FEDOR: don't skip SF monomers
    awk "{if(!(\$0 in a)){a[\$0]; print}}" $wd/_nhmmer-t1-$bn.bed > $wd/_nhmmer-t0-$bn.bed

#   FEDOR: addicional overlap filtering
    python3 ${wd}/overlap_filter.py $wd/_nhmmer-t0-$bn.bed > $wd/AS-HOR+SF-vs-$bn.bed

#   FEDOR: AS-HOR only (skip SF monomers)
    awk '{ if (length($4)==2) {next} print}' $wd/AS-HOR+SF-vs-$bn.bed > $wd/AS-HOR-vs-$bn.bed

#   Delete temporary files
    rm $wd/_nhmmer-t1-$bn.bed
    rm $wd/_nhmmer-t0-$bn.bed
    rm $wd/nhmmer-$bp-vs-$bn-tbl.bed

done
