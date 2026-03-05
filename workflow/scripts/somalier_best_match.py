
import sys
import pandas as pd
import os

def get_sample_type(sample_name):
    """
    Identifies sample type based on suffix:
    _T or _N -> DNA
    _R       -> RNA
    """
    if sample_name.endswith('_T') or sample_name.endswith('_N'):
        return 'DNA'
    elif sample_name.endswith('_R'):
        return 'RNA'
    return 'Unknown'

def main(input_file, output_file):
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found.")
        sys.exit(1)
        
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        sys.exit(1)

    if df.empty:
        print(f"Warning: Input file {input_file} is empty.")
        # Create empty output with headers
        pd.DataFrame(columns=['Sample', 'Best_Match', 'Relatedness', 'ibs0', 'ibs2', 'n']).to_csv(output_file, sep='\t', index=False)
        return

    # Column names might have leading #
    rename_dict = {col: col.lstrip('#') for col in df.columns}
    df = df.rename(columns=rename_dict)

    # Required columns
    required = ['sample_a', 'sample_b', 'relatedness', 'ibs0', 'ibs2', 'n']
    for col in required:
        if col not in df.columns:
            print(f"Error: Missing required column {col}")
            sys.exit(1)

    # Get unique samples
    samples = pd.concat([df['sample_a'], df['sample_b']]).unique()
    
    results = []
    for sample in samples:
        sample_type = get_sample_type(sample)
        
        # Matches involving this sample
        matches = df[(df['sample_a'] == sample) | (df['sample_b'] == sample)].copy()
        
        if matches.empty:
            continue
            
        # Identify the other sample in the pair
        matches['other_sample'] = matches.apply(
            lambda row: row['sample_b'] if row['sample_a'] == sample else row['sample_a'], 
            axis=1
        )
        
        # Filter for cross-type matches (DNA matches RNA, RNA matches DNA)
        if sample_type == 'DNA':
            filtered_matches = matches[matches['other_sample'].apply(lambda s: get_sample_type(s) == 'RNA')]
        elif sample_type == 'RNA':
            filtered_matches = matches[matches['other_sample'].apply(lambda s: get_sample_type(s) == 'DNA')]
        else:
            # If type is unknown, we don't apply the cross-match filter to allow some fallback?
            # Or should we exclude it? The user specifically asked for DNA/RNA matches.
            filtered_matches = pd.DataFrame() 
        
        if filtered_matches.empty:
            continue
            
        # Sort by relatedness descending to find the best match
        best_match_row = filtered_matches.sort_values(by='relatedness', ascending=False).iloc[0]
        
        results.append({
            'Sample': sample,
            'Best_Match': best_match_row['other_sample'],
            'Relatedness': best_match_row['relatedness'],
            'ibs0': best_match_row['ibs0'],
            'ibs2': best_match_row['ibs2'],
            'n': best_match_row['n']
        })
        
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        # Sort results by sample name for consistency
        res_df = res_df.sort_values(by='Sample')
    res_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python somalier_best_match.py <input_tsv> <output_tsv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
