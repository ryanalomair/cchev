#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/eval_%j.out
#SBATCH --error=logs/eval_%j.err

conda deactivate

num=$1
dest=$2
k=$3
a_dir=$4
r_dir=$5
ctype=$6

home="/home/users/ralomair1"
letter="a"
spec_dir="report"
base="${home}/${spec_dir}"
soln_dir="${home}/reflow/sim_solns/${spec_dir}/${num}cells/${a_dir}"
eval_script="/home/users/ralomair1/al/evaluate_vectors.py"
extract_script="/home/users/ralomair1/al/extract_vectors.py"

target_dir="${base}/output/${num}cells/${a_dir}/${r_dir}"

nm="c_1e"
f_nm="${num}${nm}_ct${ctype}"

if ! cd "$target_dir" 2>/dev/null; then
    echo "Error: Target directory not found: $target_dir" >&2
    exit 1
fi

pvals=""
fnames=""

# Extract vectors across all k windows
for ((j=1; j<=k; j++)) ; do
    if [[ $k -gt 1 ]] ; then
        pre="w${j}_"
    else
        pre=""
    fi

    nwk_file=$(ls "${pre}pass"*"_ct${ctype}.nwk" 2>/dev/null | sort -V | tail -n 1)

    if [ -z "$nwk_file" ]; then
        echo "Error: No tree file found matching pattern '${pre}pass*_ct${ctype}.nwk' in $target_dir" >&2
        exit 1
    fi

    pass_val=$(echo "$nwk_file" | sed -n 's/.*pass\([0-9]*\).*/\1/p')
    pass="pass${pass_val}"
    fn="${pre}${pass}_ancs_ct${ctype}.txt"

    if [ "$j" -eq 1 ] ; then
        pvals="$pass_val"
        fnames="$fn"
    else
        pvals="${pvals},${pass_val}"
        fnames="${fnames} ${fn}"
    fi

    python3 "$extract_script" "$nwk_file" "$fn"
done

# Evaluate extracted vectors against simulated ground truth
python3 "$eval_script" "$pvals" "${soln_dir}/sim_ancs.txt" "$ctype" $fnames
