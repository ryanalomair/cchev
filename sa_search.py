import pandas as pd
import numpy as np
import os
import sys
import re
import copy
import shutil
import random
import math
from Bio.Phylo.BaseTree import Clade, Tree
from Bio import Phylo

# --- Configuration & CLI Parsing ---
if len(sys.argv) < 4:
    print("Usage: python3 sa_search.py <input.tsv> <seed_val_or_-1> <K_WINDOWS> [TEMP] [COOLING_RATE] [ITERATIONS] [WINDOW_SIZE] [OVERLAP] [0/1 flag indicating duplication counting method. Default: 1]")
    print("  Use -1 for seed_val to run non-deterministically using system time.")
    sys.exit()

infile = sys.argv[1]
seed_val = int(sys.argv[2])
k_windows = int(sys.argv[3])

fn = os.path.basename(infile)
num = re.search(r'(\d+)', fn).group(1) if re.search(r'(\d+)', fn) else "out"

# Default SA Hyperparameters
init_temp = 10.0
cooling_rate = 0.98
max_iterations = 200

# Default arguments for window size extraction
window_size = 0
overlap = 0
ctype = 1

# Parse optional hyperparameters and trailing window properties
if len(sys.argv) > 4: init_temp = float(sys.argv[4])
if len(sys.argv) > 5: cooling_rate = float(sys.argv[5])
if len(sys.argv) > 6: max_iterations = int(sys.argv[6])
if len(sys.argv) > 7: window_size = int(sys.argv[7])
if len(sys.argv) > 8: overlap = int(sys.argv[8])
if len(sys.argv) > 9: ctype = int(sys.argv[9])

nm = fn.split('.')[0]
f_nm = f"{nm}_ct{ctype}"

# Set random seed if provided
if seed_val != -1:
    random.seed(seed_val)

# --- Window Slicing Logic ---
def get_windows(columns, k_windows, window_size, overlap):
    n_cols = len(columns)
    windows = []
    if window_size == 0 and overlap == 0:
        k_windows = max(1, k_windows)
        base_size = n_cols // k_windows
        remainder = n_cols % k_windows
        start = 0
        for i in range(k_windows):
            end = start + base_size + (1 if i < remainder else 0)
            windows.append(columns[start:end])
            start = end
    else:
        step = max(1, window_size - overlap)
        start = 0
        while start < n_cols:
            end = min(start + window_size, n_cols)
            windows.append(columns[start:end])
            if end == n_cols: break
            start += step
    return windows

def log_to_file(filename, msg):
    with open(filename, "a") as f:
        f.write(msg + "\n")

def append_ancestor_to_tsv(vector, label, n_cells, prefix="anc"):
    target_file = f"{prefix}_{f_nm}.tsv"
    write_header = not os.path.exists(target_file)
    row_data = vector.astype(int).to_frame().T
    row_data.insert(0, 'cell#', label)
    row_data.to_csv(target_file, sep='\t', mode='a', index=False, header=write_header)

def count_dups(vector_in, ctype):
    vec = np.array(vector_in, dtype=float).copy()
    total_dups = 0
    while np.any(vec > 0):
        is_pos = vec > 0
        runs = np.where(np.diff(np.concatenate(([False], is_pos, [False]))))[0]
        start, end = runs[0], runs[1]
        segment = vec[start:end]
        min_val = np.min(segment)
        vec[start:end] -= min_val
        if ctype == 0:
            total_dups += min_val
        else:
            total_dups += 1
    return total_dups

def get_dup_distance(vec_a, vec_b, ctype=1):
    ancestor = np.minimum(vec_a, vec_b)
    return count_dups(vec_a - ancestor, ctype) + count_dups(vec_b - ancestor, ctype)

def get_guided_simulated_merges(tree, initial_matrix, ctype=1):
    working_tree = copy.deepcopy(tree)
    curr_matrix = initial_matrix.copy()
    merges = []
    anc_idx = 1
    while True:
        ready_parents = [p for p in working_tree.get_nonterminals() if all(c.is_terminal() for c in p.clades)]
        if not ready_parents: break
        best_d = float('inf')
        best_parent, best_pair = None, None
        for parent in ready_parents:
            if len(parent.clades) >= 2:
                name_l, name_r = parent.clades[0].name, parent.clades[1].name
                if name_l in curr_matrix.index and name_r in curr_matrix.index:
                    d_val = get_dup_distance(curr_matrix.loc[name_l].values, curr_matrix.loc[name_r].values, ctype)
                    if d_val < best_d:
                        best_d = d_val
                        best_parent = parent
                        best_pair = (name_l, name_r)
        if best_parent and best_pair:
            name_l, name_r = best_pair
            anc_label = f"a{anc_idx}"
            merges.append((name_l, name_r, anc_label, anc_idx))
            anc_idx += 1
            new_vec = np.minimum(curr_matrix.loc[name_l].values, curr_matrix.loc[name_r].values)
            curr_matrix.loc[anc_label] = new_vec
            best_parent.clades = [] 
            best_parent.name = anc_label
        else:
            break
    return merges

def run_bottom_up_logic(matrix, merges, n_cells, log_file=None, prefix=None, write_disk=False, verbose=True, ctype=1):
    curr_matrix = matrix.copy()
    total_d = 0
    for name_l, name_r, anc_label, idx in merges:
        vec_l, vec_r = curr_matrix.loc[name_l].values, curr_matrix.loc[name_r].values
        new_vec = np.minimum(vec_l, vec_r)
        d_val = get_dup_distance(vec_l, vec_r, ctype)
        total_d += d_val
        if write_disk and log_file and verbose:
            log_to_file(log_file, "-" * 100)
            log_to_file(log_file, f"--- Merge Step {idx} ---")
            log_to_file(log_file, f"Pairing: ['{name_l}', '{name_r}'] -> {anc_label}")
            log_to_file(log_file, f"\nNewly Created Ancestor '{anc_label}' Profile:")
            anc_df = pd.Series(new_vec, index=matrix.columns).astype(int).to_frame().T
            log_to_file(log_file, anc_df.to_string(index=False))
            log_to_file(log_file, f"\nNormalized Profiles for Paired Cells (Cell - Ancestor '{anc_label}'):")
            diff_df = pd.DataFrame([vec_l - new_vec, vec_r - new_vec], index=[name_l, name_r], columns=matrix.columns).astype(int)
            log_to_file(log_file, diff_df.to_string())
            log_to_file(log_file, f"\nOptimization Score (Distance d({name_l},{name_r})): {d_val}\n")
        if write_disk and prefix and idx <= (n_cells - 1):
            append_ancestor_to_tsv(pd.Series(new_vec, index=matrix.columns), anc_label, n_cells, prefix=prefix)
        curr_matrix.loc[anc_label] = new_vec
        curr_matrix = curr_matrix.drop([name_l, name_r])
    return total_d

def get_parent_map(tree):
    parent_map = {}
    for clade in tree.find_clades():
        for child in clade.clades: parent_map[child] = clade
    return parent_map

def prune_leaf(tree, leaf_node):
    parent_map = get_parent_map(tree)
    parent = parent_map[leaf_node]
    sibling = parent.clades[0] if parent.clades[1] == leaf_node else parent.clades[1]
    if parent == tree.root: tree.root = sibling
    else:
        grandparent = parent_map[parent]
        idx = grandparent.clades.index(parent)
        grandparent.clades[idx] = sibling
    return Clade(name=leaf_node.name)

def regraft_leaf(tree, leaf_clade, target_node):
    parent_map = get_parent_map(tree)
    new_internal = Clade()
    if target_node == tree.root:
        new_internal.clades = [tree.root, leaf_clade]
        tree.root = new_internal
    else:
        parent = parent_map[target_node]
        idx = parent.clades.index(target_node)
        parent.clades[idx] = new_internal
        new_internal.clades = [target_node, leaf_clade]

def write_pass_files(idx, tree, merges, df_raw, global_min, norm_matrix, n_cells, score, is_final=False, w_prefix="", ctype=1):
    p_log = f"{w_prefix}pass{idx}_{f_nm}_log.txt"
    p_tsv = f"{w_prefix}pass{idx}_anc_{f_nm}.tsv"
    p_nwk = f"{w_prefix}pass{idx}_out_tree{num}_ct{ctype}.nwk"
    for f in [p_log, p_tsv, p_nwk]:
        if os.path.exists(f): os.remove(f)
    if is_final:
        log_to_file(p_log, "="*100)
        log_to_file(p_log, f"PASS {idx} VERBOSE LOG (FINAL SA OPTIMUM)")
        log_to_file(p_log, "="*100)
        log_to_file(p_log, "Initial Matrix (Raw Data):\n" + df_raw.astype(int).to_string() + "\n")
    append_ancestor_to_tsv(global_min, "f", n_cells, prefix=f"{w_prefix}pass{idx}_anc")
    run_bottom_up_logic(norm_matrix, merges, n_cells, log_file=p_log, prefix=f"{w_prefix}pass{idx}_anc", write_disk=True, verbose=is_final, ctype=ctype)
    log_to_file(p_log, f"FINAL TOTAL DUPLICATIONS: {score}")
    out_tree = copy.deepcopy(tree)
    for c in out_tree.find_clades():
        c.branch_length = None
        if not c.is_terminal(): c.name = None
    Phylo.write(out_tree, p_nwk, "newick")
    with open(p_nwk, "r") as f: content = f.read()
    with open(p_nwk, "w") as f: f.write(re.sub(r':None|:0\.0+|:\d+\.\d+', '', content).strip())

def process_sa_window(w_idx, win_cols, df_raw_full, n_cells, base_tree_file, ctype=1):
    win_num = w_idx + 1
    w_prefix = f"w{win_num}_" if k_windows > 1 else ""
    df_raw = df_raw_full[win_cols].copy()
    global_min = df_raw.min(axis=0)
    founder_d = count_dups(global_min.values, ctype)
    norm_matrix = df_raw.subtract(global_min, axis=1).astype(float)

    print(f"\nRunning Simulated Annealing on Window {win_num} using baseline: {base_tree_file}")
    
    current_tree = Phylo.read(base_tree_file, "newick")
    base_merges = get_guided_simulated_merges(current_tree, norm_matrix, ctype=ctype)
    base_score = founder_d + run_bottom_up_logic(norm_matrix, base_merges, n_cells, write_disk=False, ctype=ctype)
    
    current_best_score = base_score
    current_best_tree = copy.deepcopy(current_tree)
    current_best_merges = base_merges
    
    active_tree = copy.deepcopy(current_tree)
    active_score = base_score
    
    temp = init_temp
    pass_index = 0

    for iteration in range(max_iterations):
        if temp <= 0.0001:
            break
            
        # 1. Propose an SPR move completely at random
        test_tree = copy.deepcopy(active_tree)
        terminals = test_tree.get_terminals()
        if len(terminals) < 3:
            break
            
        leaf_to_move = random.choice(terminals)
        
        # Prune from configuration
        pruned_clade = prune_leaf(test_tree, leaf_to_move)
        
        # Pick target node uniformly at random from remaining target structures
        all_targets = list(test_tree.find_clades())
        chosen_target = random.choice(all_targets)
        
        # Regraft to create variant topology
        regraft_leaf(test_tree, pruned_clade, chosen_target)
        
        # Strip structural labeling markers to ensure deterministic merge logic paths
        for c in test_tree.find_clades():
            if not c.is_terminal(): c.name = None
            
        test_merges = get_guided_simulated_merges(test_tree, norm_matrix, ctype=ctype)
        test_score = founder_d + run_bottom_up_logic(norm_matrix, test_merges, n_cells, write_disk=False, ctype=ctype)
        
        delta_e = test_score - active_score
        
        # 2. Evaluation using the Boltzmann distribution criteria 
        accept_move = False
        if delta_e < 0:
            accept_move = True
        else:
            prob = math.exp(-delta_e / temp)
            if random.random() < prob:
                accept_move = True
                
        if accept_move:
            active_tree = test_tree
            active_score = test_score
            
            # If we uncovered a new global optimum, log its milestone pass
            if active_score < current_best_score:
                current_best_score = active_score
                current_best_tree = copy.deepcopy(active_tree)
                current_best_merges = test_merges
                
                write_pass_files(pass_index, current_best_tree, current_best_merges, df_raw, global_min, norm_matrix, n_cells, current_best_score, is_final=False, w_prefix=w_prefix, ctype=ctype)
                pass_index += 1
                
        # 3. Apply cooling scheduler properties
        temp *= cooling_rate

    # Final outputs rendering pass matching original metrics execution layout cleanly
    if pass_index == 0:
        p0_log = f"{w_prefix}pass0_{f_nm}_log.txt"
        p0_tsv = f"{w_prefix}pass0_anc_{f_nm}.tsv"
        p0_nwk = f"{w_prefix}pass0_out_tree{num}_ct{ctype}.nwk"
        for f in [p0_log, p0_tsv, p0_nwk]:
            if os.path.exists(f): os.remove(f)
        append_ancestor_to_tsv(global_min, "f", n_cells, prefix=f"{w_prefix}pass0_anc")
        run_bottom_up_logic(norm_matrix, current_best_merges, n_cells, log_file=p0_log, prefix=f"{w_prefix}pass0_anc", write_disk=True, verbose=True, ctype=ctype)
        shutil.copy(base_tree_file, p0_nwk)
    else:
        final_idx = pass_index - 1
        write_pass_files(final_idx, current_best_tree, current_best_merges, df_raw, global_min, norm_matrix, n_cells, current_best_score, is_final=True, w_prefix=w_prefix, ctype=ctype)

    print(f"SA Window Complete -> Base Agglomerative: {base_score} | Annealed Minimum: {current_best_score}")

def main():
    df_raw_full = pd.read_csv(infile, sep='\t', index_col=0)
    df_raw_full.index = df_raw_full.index.astype(str)
    n_cells = len(df_raw_full.index)

    windows = get_windows(df_raw_full.columns.tolist(), k_windows, window_size, overlap)

    if k_windows == 1:
        expected_nwk = f"out_tree{num}_ct{ctype}.nwk"
        if os.path.exists(expected_nwk):
            process_sa_window(0, df_raw_full.columns.tolist(), df_raw_full, n_cells, expected_nwk, ctype=ctype)
        else:
            print(f"Error: Baseline tree file '{expected_nwk}' not found. Ensure bu.py has run successfully.")
    else:
        for w_idx, win_cols in enumerate(windows):
            win_num = w_idx + 1
            w_prefix = f"w{win_num}_"
            expected_nwk = f"{w_prefix}out_tree{num}_ct{ctype}.nwk"
            if os.path.exists(expected_nwk):
                process_sa_window(w_idx, win_cols, df_raw_full, n_cells, expected_nwk, ctype=ctype)
            else:
                print(f"Warning: Expected baseline tree file '{expected_nwk}' not found. Skipping window {win_num}.")

if __name__ == "__main__":
    main()
