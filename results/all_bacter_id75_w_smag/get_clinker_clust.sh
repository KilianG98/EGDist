#!/bin/bash

# Define variables
superclusters="GCF_000012725_000000000001_gene_clusters.csv"
genomes=("GCF_000013885" "GCA_001897285" "GCF_002028545" "GCF_000152905" "GCF_000012725" "GCA_001896955" "GCF_006539545")
dir=$(pwd)
outputs=()

cur_acc=()
cur_starts=()
cur_names=()

read -r header < "$superclusters"
# Main loop to process each line in the superclusters file
while IFS=',' read -r line; do
    acc=$(echo "$line" | awk -F ',' '{print $NF}')
    start=$(echo "$line" | cut -d ',' -f 2)
    name=$(echo "$line" | cut -d ',' -f 5)
    echo $acc $cur_acc ${cur_acc[${#cur_acc[@]}-1]} 
    # Handle first entry or consecutive entries
    if [ -z "${cur_acc}" ] || [ "$acc" -le $(( ${cur_acc[${#cur_acc[@]}-1]} + 5 )) ]; then

        cur_acc+=("$acc")
        cur_starts+=("$start")
        cur_names+=("$name")
        continue
    fi

    # Check if there are at least two entries
    if [ ${#cur_acc[@]} -ge 2 ]; then
        start=$(( ${cur_starts[0]} - 2000 ))
        if [ "$start" -lt 0 ]; then
            start=0
        fi
        stop=$(( ${cur_starts[-1]} + 3000 ))

        # Define output array
        output=("GCF_000012725_000000000001" "$start" "$stop")

        # Iterate over genomes
        for genome in "${genomes[@]}"; do
            for file in "${genome}"*.csv; do
                c_starts=()
                contig=$(basename "${dir}/${file}" "_gene_clusters.csv")
                # Populate c_starts array with relevant data
                while read -r name; do
                    line=$(grep "$name" "$file")
		    echo $line $name $file
                    if [ -n "$line" ]; then
                        c_starts+=("$(echo "$line" | cut -d ',' -f 2)")
                    fi
                done <<< "${cur_names[@]}"

                # Check if c_starts array is not empty
                if [ ${#c_starts[@]} -ne 0 ]; then
                    start=$(printf "%s\n" "${c_starts[@]}" | sort -n | head -n 1)
                    start=$((start - 2000))
                    if [ "$start" -lt 0 ]; then
                        start=0
                    fi
                    end=$(printf "%s\n" "${c_starts[@]}" | sort -n | tail -n 1)
                    end=$((end + 3000))
                    output+=("$contig" "$start" "$end")
                fi
            done
        done

        # Append output to outputs array
        outputs+=("${output[@]}")
    fi
done < "$superclusters"

# Write outputs to file
for op in "${outputs[@]}"; do
    echo "$op" >> "extract_this.txt"
done
