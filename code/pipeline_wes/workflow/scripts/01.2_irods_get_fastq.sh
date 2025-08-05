#!/bin/bash

while getopts ":d:a:b:x:y:m:n:" opt; do
    case $opt in
	d) dir="$OPTARG"
	    ;;
	a) irods_1="$OPTARG"
	    ;;
	b) irods_2="$OPTARG"
	    ;;
	x) local_1="$OPTARG"
	    ;;
	y) local_2="$OPTARG"
	    ;;
	m) name_1="$OPTARG"
	    ;;
	n) name_2="$OPTARG"
	    ;;
	\?) echo "Invalid option -$OPTARG" >&2
	    exit 1
	    ;;
    esac

    case $OPTARG in
	-*) echo "Option $opt needs a valid argument"
	    exit 1
	    ;;
    esac
done

raw_1=${dir}/${name_1}
raw_2=${dir}/${name_2}


if [[ "${irods_1}" == *","* ]]; then
    # multiple fastqs - concatenate them
    IFS=',' read -r -a irods_1_all <<< "${irods_1}"
    IFS=',' read -r -a name_1_all <<< "${name_1}"
    for index in "${!irods_1_all[@]}"; do
        iget -vK ${irods_1_all[index]} ${dir}
        cat ${dir}/${name_1_all[index]} >> ${local_1}
        echo "concatenated ${dir}/${name_1_all[index]} to ${local_1}!"
        rm ${dir}/${name_1_all[index]}
    done
else
    # single fastq
    iget -vK ${irods_1} ${dir}

    if [[ "${raw_1}" != "${local_1}" ]]; then
        mv ${raw_1} ${local_1}
        echo "moved ${raw_1} to ${local_1}!"
    fi
fi


if [[ "${irods_2}" == *","*  ]]; then
    # multiple fastqs - concatenate them
    IFS=',' read -r -a irods_2_all <<< "${irods_2}"
    IFS=',' read -r -a name_2_all <<< "${name_2}"
    for index in "${!irods_2_all[@]}"; do
        iget -vK ${irods_2_all[index]} ${dir}
        cat ${dir}/${name_2_all[index]} >> ${local_2}
        echo "concatenated ${dir}/${name_2_all[index]} to ${local_2}!"
        rm ${dir}/${name_2_all[index]}
    done
else
    # single fastq
    iget -vK ${irods_2} ${dir}

    if [[ "${raw_2}" != "${local_2}" ]]; then
        mv ${raw_2} ${local_2}
        echo "moved ${raw_2} to ${local_2}!"
    fi
fi
