#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ryalomair@gmail.com

conda deactivate

num=$1
dest=$2
ctype=$3
letter="a"

home="/home/users/ralomair1"
spec_dir="report"
base="${home}/${spec_dir}/output/${num}cells"
sim_base="${home}/reflow/sim_solns/${spec_dir}/${num}cells"

nm="c_1e"
f_nm="${num}${nm}_ct${ctype}"

for p in 1 ; do
    results=()
    total_k=0
    total_ties=0
    count=0
    possible_correct=$(( num - 2 ))
    total_ancs_missed=0
    real_total_ancs_missed=0
    output_tsv="${dest}_${num}cells.tsv"
    num_perf=0
    one_off=0
    good_rshfls=0
    diff_sum=0
    real_diff_sum=0
    
    # Tracking aggregates
    mdt_total=0
    mdt_inst_total=0
    total_real_tied_dirs=0
    total_num_tied_dirs=0

    for i in {0..49}; do
        a_dir="a${i}"
        
        min_dups=999999
        max_k=-1
        actual_max_k=-1
        
        # Track profiles across all 10 directories for this instance
        unset r_dups r_k valid_rs
	declare -A r_dups
        declare -A r_k
	declare -a valid_rs=()

        sim_dups=$(tail -n 1 "${sim_base}/${a_dir}/${f_nm}_sim_dups.txt" | awk '{print $NF}')
        
        # First Pass: Gather data and find the true minimum duplications
        for r in {0..9}; do
            r_dir="r${r}"
            target_dir="${base}/${a_dir}/${r_dir}"
            target_file="${target_dir}/sim_matches.txt"
            
            #the following will have to be adjusted if using k > 1 windows
	    latest_file=$(ls "${target_dir}"/*pass*"_${f_nm}_log.txt" 2>/dev/null | sort -V | tail -n 1)
            if [ -n "$latest_file" ]; then
                if ! grep -qE '^FINAL TOTAL DUPLICATIONS: [0-9]+(\.[0-9]+)?$' "$latest_file"; then
                    fallback_file="${target_dir}/${f_nm}_log.txt" # WILL ALSO HAVE TO ADJUST IF USING WINDOWS
                    if [ -f "$fallback_file" ]; then
                        line=$(tail -n 1 "$fallback_file")
                        echo "$line" >> "$latest_file"
                    fi
                fi
            fi

            if [ -f "$target_file" ] ; then
                dups=$(tail -n 1 "$latest_file" | awk '{print $NF}' | awk -F'.' '{print $1}')
                k=$(wc -l < "$target_file")
                
                r_dups["$r_dir"]=$dups
                r_k["$r_dir"]=$k
                valid_rs+=("$r_dir")

                if [ "$dups" -lt "$min_dups" ]; then
                    min_dups=$dups
                fi
            fi
        done

        # Second Pass: Evaluate ties strictly within the minimum duplication subset
        best_r_dirs=()
        min_dup_dirs_count=0
        real_dirs=()

        for r_dir in "${valid_rs[@]}"; do
            if [ "${r_dups[$r_dir]}" -eq "$min_dups" ]; then
                ((min_dup_dirs_count++))
                curr_k=${r_k[$r_dir]}
                real_dirs+=("$r_dir")
                
                if [ "$curr_k" -gt "$max_k" ]; then
                    max_k=$curr_k
                    best_r_dirs=("$r_dir")
                elif [ "$curr_k" -eq "$max_k" ]; then
                    best_r_dirs+=("$r_dir")
                fi
            fi
        done

        # Determine exact scenarios
        num_tied_dirs=${#best_r_dirs[@]}
        total_num_tied_dirs=$(( total_num_tied_dirs + num_tied_dirs ))
        real_tied_dirs=${#real_dirs[@]}
        total_real_tied_dirs=$(( total_real_tied_dirs + real_tied_dirs ))

        # Bad Ties: min_dup is shared, but we had to filter some out due to worse ancestor matches
        if [ "$min_dup_dirs_count" -gt 1 ] && [ "$min_dup_dirs_count" -ne "$num_tied_dirs" ]; then
            min_dup_ties=$(( min_dup_dirs_count - num_tied_dirs ))
            mdt_inst=1
        else
            min_dup_ties=0
            mdt_inst=0
        fi

        mdt_total=$(( mdt_total + min_dup_ties ))
        mdt_inst_total=$(( mdt_inst_total + mdt_inst ))

        # Process the selected directory
        if [ ${#best_r_dirs[@]} -gt 0 ]; then
            if [ $num_tied_dirs -eq 1 ]; then
                ((good_rshfls++))
            fi

            random_idx=$(( RANDOM % num_tied_dirs ))
            real_rand_idx=$(( RANDOM % real_tied_dirs ))
            
            selected_real_dir="${real_dirs[$real_rand_idx]}"
            real_max_k="${r_k[$selected_real_dir]}"

            selected_r_dir="${best_r_dirs[$random_idx]}"
            selected_dir_path="${base}/${a_dir}/${selected_r_dir}"
            selected_real_path="${base}/${a_dir}/${selected_real_dir}"

            final_dups=${r_dups[$selected_r_dir]}
            final_real_dups=${r_dups[$selected_real_dir]}
            final_diff=$(( final_dups - sim_dups ))
            final_real_diff=$(( final_real_dups - sim_dups ))
            diff_sum=$(( diff_sum + final_diff ))
            real_diff_sum=$(( real_diff_sum + final_real_diff ))

            max_k=$(( max_k - 1 )) # adjust for founder
            real_max_k=$(( real_max_k - 1 ))
            ancs_missed=$(( possible_correct - max_k ))
            real_ancs_missed=$(( possible_correct - real_max_k ))
            total_ancs_missed=$(( total_ancs_missed + ancs_missed ))
            real_total_ancs_missed=$(( real_total_ancs_missed + real_ancs_missed ))
            
            if [ $real_ancs_missed -eq 0 ]; then ((num_perf++)); fi
            if [ $real_ancs_missed -eq 1 ]; then ((one_off++)); fi

            # Format strictly maps columns: inst, real misses, lowest misses, dup_diff, real_ties, best_ties
            results+=("${a_dir}"$'\t'"${selected_real_dir}: ${real_ancs_missed}"$'\t'"${selected_r_dir}: ${ancs_missed}"$'\t'"${final_real_diff}"$'\t'"${real_tied_dirs}"$'\t'"${num_tied_dirs}")

            total_k=$(( total_k + max_k ))
            total_ties=$(( total_ties + num_tied_dirs ))
            count=$(( count + 1 ))
        else
            echo "Warning: No valid files found in ${a_dir} for ${num}cells" >&2
        fi
    done

    # Calculate averages cleanly
    if [ "$count" -gt 0 ]; then
        avg_misses=$(awk -v total="$total_ancs_missed" -v cnt="$count" 'BEGIN { printf "%.2f", total / cnt }')
        avg_ties=$(awk -v total="$total_ties" -v cnt="$count" 'BEGIN { printf "%.2f", total / cnt }')
        real_error=$(echo "scale=2; ( $real_total_ancs_missed * 100 ) / ( $count * $possible_correct )" | bc | sed 's/^\./0./; s/^-\./-0./')
        best_error=$(echo "scale=2; ( $total_ancs_missed * 100 ) / ( $count * $possible_correct )" | bc | sed 's/^\./0./; s/^-\./-0./')
        avg_diff=$(awk -v total="$real_diff_sum" -v cnt="$count" 'BEGIN { printf "%.2f", total / cnt }')
        avg_real_ties=$(awk -v total="$total_real_tied_dirs" -v cnt="$count" 'BEGIN { printf "%.2f", total / cnt }')
        avg_best_ties=$(awk -v total="$total_num_tied_dirs" -v cnt="$count" 'BEGIN { printf "%.2f", total / cnt }')
    else
        avg_misses=0; avg_ties=0; real_error=0; best_error=0; avg_diff=0; avg_real_ties=0; avg_best_ties=0
    fi
    total_possible=$(( possible_correct * count ))
    
    # Write output to TSV file
    echo "average misses: ${avg_misses}" > "$output_tsv"
    echo "average tied directories: ${avg_ties}" >> "$output_tsv"
    echo -e "\t${num} cells" >> "$output_tsv"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "inst" "real misses" "lowest misses" "dup_diff" "real ties" "best ties" >> "$output_tsv"
    
    for line in "${results[@]}"; do
        echo "$line" >> "$output_tsv"
    done
    
    # Aligned summary metrics line up underneath columns
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "" "${real_total_ancs_missed}" "${total_ancs_missed}" "${avg_diff} avg" "${avg_real_ties} avg" "${avg_best_ties} avg" >> "$output_tsv"
    printf "%s\t%s\t%s\t%s\t%s\n" "" "${num_perf} perf" "${one_off} 1_off" "${real_error}% err" "${best_error}% best err" >> "$output_tsv"
done
