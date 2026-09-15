#!/usr/bin/env bash
# Converts genome bins in fasta format (plain .fa or gzipped .fa.gz) to a
# DAS_Tool contig2bin TSV (contig_id<TAB>bin_id, no header).

function display_help() {
    echo " "
    echo "fasta_to_contig2bin: Converts genome bins in fasta format to contig-to-bin table."
    echo " "
    echo "Usage: fasta_to_contig2bin.sh -i <bins_folder> -e <extension> > my_contig2bin.tsv"
    echo " "
    echo "   -e, --extension            Extension of fasta files, without .gz (default: fasta)"
    echo "   -i, --input_folder         Folder with bins in fasta format. (default: ./)"
    echo "   -h, --help                 Show this message."
    echo " "
    echo "Matches both plain (<extension>) and gzipped (<extension>.gz) files."
    echo " "
    exit 1
}

extension="fasta"
folder="."

while [ "$1" != "" ]; do
    case $1 in
        -e | --extension )      shift
                                extension=$1
                                ;;
        -i | --input_folder )   shift
                                folder=$1
                                ;;
        -h | --help )           display_help
                                exit
                                ;;
        * )                     display_help
                                exit 1
    esac
    shift
done

shopt -s nullglob

for i in "$folder"/*."$extension" "$folder"/*."$extension".gz
do
    binname=$(basename "$i")
    binname="${binname%.gz}"
    binname="${binname%.$extension}"

    if [[ "$i" == *.gz ]]; then
        zcat "$i"
    else
        cat "$i"
    fi | grep ">" | perl -pe "s/\n/\t$binname\n/g" | perl -pe "s/>//g"
done
