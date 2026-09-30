import pandas as pd
import numpy as np
import os
import sys
import re
import copy
import shutil
from Bio.Phylo.BaseTree import Clade, Tree
from Bio import Phylo

# --- Configuration & CLI Parsing ---
if len(sys.argv) < 2:
    print("Usage: python3 bu.py <input.tsv> [simulated_tree.nwk] [0/1 flag] [K_WINDOWS] [WINDOW_SIZE] [OVERLAP] [0/1 flag for duplication counting method. Default 1 -> count 1]")
    sys.exit()

infile = sys.argv[1]
fn = os.path.basename(infile)
num = re.search(r'(\d+)', fn).group(1) if re.search(r'(\d+)', fn) else "out"

# Default arguments
SIM_TREE_FILE = None
first_improvement_flag = 1
k_windows = 1
window_size = 0
overlap = 0
ctype = 1

# Parse optional positional arguments
curr_idx = 2
if curr_idx < len(sys.argv) and not sys.argv[curr_idx].isdigit():
    SIM_TREE_FILE = sys.argv[curr_idx]
    curr_idx += 1

if curr_idx < len(sys.argv):
    first_improvement_flag = int(sys.argv[curr_idx])
    curr_idx += 1
if curr_idx < len(sys.argv):
    k_windows = int(sys.argv[curr_idx])
    curr_idx += 1
if curr_idx < len(sys.argv):
    window_size = int(sys.argv[curr_idx])
    curr_idx += 1
if curr_idx < len(sys.argv):
    overlap = int(sys.argv[curr_idx])
    curr_idx += 1
if curr_idx < len(sys.argv):
    ctype = int(sys.argv[curr_idx])

nm = fn.split('.')[0]
f_nm = f"{nm}_ct{ctype}"

# --- Window Logic ---
def get_windows(columns, k_windows, window_size, overlap):
    n_cols = len(columns)
    windows = []
    
    if window_size == 0 and overlap == 0:
        # Case A: Divide into non-overlapping contiguous slices
        k_windows = max(1, k_windows)
        base_size = n_cols // k_windows
        remainder = n_cols % k_windows
        start = 0
        for i in range(k_windows):
            end = start + base_size + (1 if i < remainder else 0)
            windows.append(columns[start:end])
            start = end
    else:
        # Case B: Sliding window across the sequence
        step = max(1, window_size - overlap)
        start = 0
        while start < n_cols:
            end = min(start + window_size, n_cols)
            windows.append(columns[start:end])
            if end == n_cols:
                break
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

def count_dups(vector_in, ctype=1):
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

def get_dup_distance(vec_a, vec_b, ctype):
    ancestor = np.minimum(vec_a, vec_b)
    return count_dups(vec_a - ancestor, ctype) + count_dups(vec_b - ancestor, ctype)

def get_guided_simulated_merges(tree, initial_matrix, ctype):
    working_tree = copy.deepcopy(tree)
    curr_matrix = initial_matrix.copy()
    merges = []
    anc_idx = 1
    
    while True:
        ready_parents = [p for p in working_tree.get_nonterminals() if all(c.is_terminal() for c in p.clades)]
        if not ready_parents:
            break
            
        best_d = float('inf')
        best_parent = None
        best_pair = None
        
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
        vec_l = curr_matrix.loc[name_l].values
        vec_r = curr_matrix.loc[name_r].values
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
            diff_df = pd.DataFrame([vec_l - new_vec, vec_r - new_vec], 
                                   index=[name_l, name_r], columns=matrix.columns).astype(int)
            log_to_file(log_file, diff_df.to_string())
            log_to_file(log_file, f"\nOptimization Score (Distance d({name_l},{name_r})): {d_val}\n")

        if write_disk and prefix and idx <= (n_cells - 1):
            append_ancestor_to_tsv(pd.Series(new_vec, index=matrix.columns), anc_label, n_cells, prefix=prefix)
            
        curr_matrix.loc[anc_label] = new_vec
        curr_matrix = curr_matrix.drop([name_l, name_r])
        
        if write_disk and log_file and verbose and len(curr_matrix) > 1:
            log_to_file(log_file, "Current Matrix State (Active Nodes for next merge):")
            log_to_file(log_file, curr_matrix.astype(int).to_string() + "\n")
            
    return total_d

# --- Tree Manipulation ---
def get_parent_map(tree):
    parent_map = {}
    for clade in tree.find_clades():
        for child in clade.clades:
            parent_map[child] = clade
    return parent_map

def prune_leaf(tree, leaf_node):
    parent_map = get_parent_map(tree)
    parent = parent_map[leaf_node]
    sibling = parent.clades[0] if parent.clades[1] == leaf_node else parent.clades[1]
    
    if parent == tree.root:
        tree.root = sibling
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

def main(infile):
    # Clear out any old passX_ files or windowed files to avoid conflicts
    #for f in os.listdir('.'):
    #    if ('pass' in f or re.match(r'^w\d+_', f)) and (f.endswith('.txt') or f.endswith('.tsv') or f.endswith('.nwk')):
    #        os.remove(f)
            
    df_raw_full = pd.read_csv(infile, sep='\t', index_col=0)
    df_raw_full.index = df_raw_full.index.astype(str)
    n_cells = len(df_raw_full.index)

    windows = get_windows(df_raw_full.columns.tolist(), k_windows, window_size, overlap)

    for w_idx, win_cols in enumerate(windows):
        win_num = w_idx + 1
        w_prefix = f"w{win_num}_" if len(windows) > 1 else ""
        print(f"Processing Window {win_num} of {k_windows}...")
        LOG_SIM = f"{w_prefix}sim_{f_nm}_log.txt"
        LOG_BU = f"{w_prefix}{f_nm}_log.txt"

        for f in [LOG_SIM, LOG_BU]:
            if os.path.exists(f): os.remove(f)

        df_raw = df_raw_full[win_cols].copy()

        global_min = df_raw.min(axis=0)
        founder_d = count_dups(global_min.values, ctype)
        norm_matrix = df_raw.subtract(global_min, axis=1).astype(float)
        
        print(f"\n{'='*50}\nProcessing Window {w_idx + 1}/{len(windows)}: {len(win_cols)} bins\n{'='*50}")

        # --- 1. SIMULATED LOGGING ---
        sim_total_dups = 0
        if SIM_TREE_FILE:
            log_to_file(LOG_SIM, "="*100)
            log_to_file(LOG_SIM, "SIMULATED SOLUTION SCORING (BOTTOM-UP TRAVERSAL)")
            log_to_file(LOG_SIM, "="*100)
            append_ancestor_to_tsv(global_min, "f", n_cells, prefix=f"{w_prefix}sim_anc")
            sim_tree = Phylo.read(SIM_TREE_FILE, "newick")
            for leaf in sim_tree.get_terminals():
                leaf.name = re.search(r'(\d+)', leaf.name).group(1)

            sim_merges = get_guided_simulated_merges(sim_tree, norm_matrix, ctype=ctype)
            sim_total_dups = founder_d + run_bottom_up_logic(norm_matrix, sim_merges, n_cells, LOG_SIM, prefix=f"{w_prefix}sim_anc", write_disk=True, ctype=ctype)
            log_to_file(LOG_SIM, f"FINAL SIMULATED TOTAL DUPLICATIONS: {sim_total_dups}")

        # --- 2. BOTTOM-UP LOGGING (Base Agglomerative) ---
        log_to_file(LOG_BU, "="*100)
        log_to_file(LOG_BU, "AGGLOMERATIVE DUPLICATION CLUSTERING (VERBOSE LOG)")
        log_to_file(LOG_BU, "="*100)
        
        TSV_BU = f"{w_prefix}anc_{f_nm}.tsv"
        if os.path.exists(TSV_BU): os.remove(TSV_BU)
        append_ancestor_to_tsv(global_min, "f", n_cells, prefix=f"{w_prefix}anc")
        
        current_matrix = norm_matrix.copy()
        nodes = {idx: Clade(name=idx) for idx in current_matrix.index}
        dynamic_merges = []
        ancestor_idx = 1

        while len(current_matrix) > 1:
            best_d, pair = float('inf'), (None, None)
            indices = current_matrix.index.tolist()
            for i in range(len(indices)):
                for j in range(i + 1, len(indices)):
                    d_ij = get_dup_distance(current_matrix.loc[indices[i]].values, current_matrix.loc[indices[j]].values, ctype)
                    if d_ij < best_d:
                        best_d, pair = d_ij, (indices[i], indices[j])
            
            idx_i, idx_j = pair
            anc_label = f"a{ancestor_idx}"
            dynamic_merges.append((idx_i, idx_j, anc_label, ancestor_idx))
            
            new_clade = Clade(name=anc_label)
            new_clade.clades = [nodes[idx_i], nodes[idx_j]]
            nodes[anc_label] = new_clade
            
            current_matrix.loc[anc_label] = np.minimum(current_matrix.loc[idx_i].values, current_matrix.loc[idx_j].values)
            current_matrix = current_matrix.drop([idx_i, idx_j])
            del nodes[idx_i], nodes[idx_j]
            ancestor_idx += 1

        total_inferred_d = founder_d + run_bottom_up_logic(norm_matrix, dynamic_merges, n_cells, LOG_BU, prefix=f"{w_prefix}anc", write_disk=True, ctype=ctype)
        log_to_file(LOG_BU, f"FINAL TOTAL DUPLICATIONS: {total_inferred_d}")
        

        final_tree = Tree(root=nodes[list(nodes.keys())[0]])
        for c in final_tree.find_clades():
            c.branch_length = None
            if not c.is_terminal(): c.name = None
        out_nwk = f"{w_prefix}out_tree{num}_ct{ctype}.nwk"
        Phylo.write(final_tree, out_nwk, "newick")
        with open(out_nwk, "r") as f: content = f.read()
        with open(out_nwk, "w") as f: f.write(re.sub(r':None|:0\.0+|:\d+\.\d+', '', content).strip())

        print(f"-----------------------------------")
        print(f"Agglomerative Total: {total_inferred_d}")
        if SIM_TREE_FILE: print(f"Simulated Total: {sim_total_dups}")
        print(f"SIM_SCORE_START\n{sim_total_dups}\nSIM_SCORE_END")
        print(f"-----------------------------------\n")

if __name__ == "__main__":
    main(sys.argv[1])
