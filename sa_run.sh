#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ryalomair@gmail.com

conda deactivate

num_cells=$1
letter=$2
i_num=$3
input_matrix=$4
seed=$5
k=$6
temp=$7
cooling_rate=$8
iterations=$9
window_size=${10}
overlap=${11}
ctype=${12}

zero_one="1"
home="/home/users/ralomair1"
spec_dir="report"
base="${home}/${spec_dir}"
fp=$(pwd)
curr_dir=$(basename "$fp")
nm="c_1e"
f_nm="${num_cells}${nm}_ct${ctype}"
out="${f_nm}_ct${ctype}_out.txt"
in_dir="${home}/reflow/sim_solns/${spec_dir}/${num_cells}cells/${letter}${i_num}"

head -n 1 "${in_dir}/${num_cells}${nm}.tsv" > "anc_${f_nm}.tsv"
head -n 1 "${in_dir}/${num_cells}${nm}.tsv" > "sim_anc_${f_nm}.tsv"

all_output=$(python3 "${base}/bu.py" "$input_matrix" "${in_dir}/tree${num_cells}.nwk" "$zero_one" "$k" "$window_size" "$overlap" "$ctype")
sa_output=$(python3 "${base}/sa_search.py" "$input_matrix" "$seed" "$k" "$temp" "$cooling_rate" "$iterations" "$window_size" "$overlap" "$ctype")

echo "$all_output" | sed '/SIM_SCORE_START/,/SIM_SCORE_END/d' > "$out" ; echo "" >> "$out"
echo "$sa_output" >> "$out"

rm -f gmon.out
