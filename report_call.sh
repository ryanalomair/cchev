#!/bin/bash

# Define usage display functions for mode-specific feedback
usage_s() {
    echo "Usage (Simulated Annealing): $0 s <num_cells> <k> <window_size> <overlap> <temp> <cooling_rate> <iterations>"
    exit 1
}

usage_r() {
    echo "Usage (Recursive Rearrangement): $0 r <num_cells> <k> <window_size> <overlap>"
    exit 1
}

# 1. Ensure at least mode flag is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <mode: s|r> [mode-specific parameters...]"
    echo "  s : Simulated Annealing (requires 7 additional parameters)"
    echo "  r : Recursive Rearrangement (requires 4 additional parameters)"
    exit 1
fi

mode="$1"
ctype=0
# 2. Parse arguments dynamically based on selected mode (Option A)
if [ "$mode" == "s" ]; then
    if [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ] || [ -z "$8" ]; then
        usage_s
    fi
    num_cells=$2
    k=$3
    window_size=$4
    overlap=$5
    temp=$6
    cooling_rate=$7
    iterations=$8
    
    # Calculate integer crate representation for SA naming convention
    crate=$(echo "$cooling_rate" | awk -F'.' '{ if (NF>1) printf "%.0f\n", $0 * (10 ^ length($2)); else print $1 }')
    
    #ctype=0
    dest="busa_k${k}_ws${window_size}_o${overlap}_ct${ctype}_t${temp}_cr${crate}_i${iterations}"
    run_script="/home/users/ralomair1/report/sa_run.sh"

elif [ "$mode" == "r" ]; then
    if [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
        usage_r
    fi
    num_cells=$2
    k=$3
    window_size=$4
    overlap=$5
    
    #ctype=0
    dest="burr_k${k}_ws${window_size}_o${overlap}_ct${ctype}"
    run_script="/home/users/ralomair1/report/recur_run.sh"

else
    echo "Error: Invalid mode '$mode'. First parameter must be 's' (Simulated Annealing) or 'r' (Recursive Rearrangement)."
    exit 1
fi

# Configuration Constants
letter="a"
SEEDS=(42 101 555 777 1234 9999 2026 8888 711 13)
home="/home/users/ralomair1"
spec_dir="report"
base="${home}/${spec_dir}"
sim_base="${home}/reflow/sim_solns/${spec_dir}"
single_eval_script="${base}/eval_single.sh"
agg_eval_script="${base}/ag.sh"
nm="c_1e"
# Save execution root directory
start_dir=$(pwd)

# Array to collect all 1,000 single-evaluation Job IDs
eval_job_ids=()

for p in 1 ; do
    sim_cells="${sim_base}/${num_cells}cells"
    for i_num in {0..49} ; do
        a_dir="${letter}${i_num}"
        input_matrix="${sim_cells}/${a_dir}/${num_cells}${nm}.tsv"
        dir="/home/users/ralomair1/${spec_dir}/output/${num_cells}cells/${a_dir}"
        
        for i in {0..9} ; do
            r_dir="r${i}"
            mkdir -p "${dir}/${r_dir}"
            cd "${dir}/${r_dir}" || exit 1
            #rm -rf *
            mkdir -p logs
            seed=${SEEDS[$i]}

            # 1. Launch Optimization Job (Mode-Specific script & argument list)
            if [ "$mode" == "s" ]; then
                res_opt=$(sbatch "$run_script" "$num_cells" "$letter" "$i_num" "$input_matrix" "$seed" "$k" "$temp" "$cooling_rate" "$iterations" "$window_size" "$overlap" "$ctype")
            else
                res_opt=$(sbatch "$run_script" "$num_cells" "$letter" "$i_num" "$input_matrix" "$seed" "$k" "$window_size" "$overlap" "$ctype")
            fi
            opt_jid=$(echo "$res_opt" | awk '{print $4}')

            # 2. Launch Per-Directory Evaluation Job immediately dependent on optimization run completing
            res_eval=$(sbatch --dependency=afterany:"$opt_jid" "$single_eval_script" "$num_cells" "$dest" "$k" "$a_dir" "$r_dir" "$ctype")
            eval_jid=$(echo "$res_eval" | awk '{print $4}')

            # Store the evaluation Job ID for final aggregation dependency
            eval_job_ids+=("$eval_jid")
        done
    done
done

# Return to starting execution directory
cd "$start_dir" || exit 1

# 3. Format evaluation Job IDs into a colon-separated list
dep_list=$(IFS=:; echo "${eval_job_ids[*]}")

# 4. Submit ag.sh dependent on ALL 1,000 evaluation jobs completing
sbatch --dependency=afterany:"$dep_list" "$agg_eval_script" "$num_cells" "$dest" "$ctype"
