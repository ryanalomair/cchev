#!/bin/bash
#SBATCH --job-name=medicc2_batch
#SBATCH --output=medicc2_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ryalomair@gmail.com

# 1. Clear environment variables that conflict with Python
unset PYTHONPATH

# 2. Initialize Conda for subshell execution and activate environment
eval "$(conda shell.bash hook)"
conda activate medicc_env

letter="a"
spec_dir="report"
tr_script="/home/users/ralomair1/medicc2_stuff/fix_mdc_tree.py"
home="/home/users/ralomair1"
base="${home}/${spec_dir}"
nm="c_mdc_1e"

for num_cells in 5 10 25 ; do
    f_nm="${num_cells}${nm}"
    for i in {0..49} ; do
        a_dir="${letter}${i}"
        in_dir="${home}/reflow/sim_solns/${spec_dir}/${num_cells}cells/${a_dir}"
        out_dir="${base}/output/${num_cells}cells/${a_dir}"
        mkdir -p "$out_dir"
	
        profiles_file="${in_dir}/profiles.tsv"
        in_file="${out_dir}/${f_nm}.tsv"
        raw_tree="${out_dir}/${f_nm}_final_tree.new"
        adj_tree="${out_dir}/m_tree${num_cells}.nwk"

        # Check if the input directory exists
        if [ -d "$in_dir" ]; then

            # Auto-convert profiles.tsv if {num_cells}c_mdc_1e.tsv doesn't exist in out_dir
            if [ ! -f "$in_file" ] && [ -f "$profiles_file" ]; then
                awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print "sample_id","chrom","start","end","cn_a","cn_b"} NR>1 {split($5, cn, ","); print $1, $2, $3, $4, cn[1], cn[2]}' "$profiles_file" > "$in_file"
            fi

            # Execute MEDICC2 if the input file is ready
            if [ -f "$in_file" ]; then
                echo "Running MEDICC2 for ${num_cells} cells (instance ${letter}${i})..."
                medicc2 "$in_file" "$out_dir" -j "$SLURM_CPUS_PER_TASK"
            else
                echo "Warning: Neither ${f_nm}.tsv nor profiles.tsv found in ${in_dir}. Skipping."
            fi

        else
            echo "Warning: Directory ${in_dir} does not exist. Skipping."
        fi

        # Post-process tree if generated
        if [ -f "$raw_tree" ]; then
            python "$tr_script" "$raw_tree" > "$adj_tree"
        fi

    done
done

echo "Batch MEDICC2 processing complete."
