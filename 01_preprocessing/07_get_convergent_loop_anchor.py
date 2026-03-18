import pandas as pd

# 1. Read data; use sep=r'\s+' to match any whitespace (space or tab)
file_path = 'K562_CTCF_merge_distance_Convergent_Loops_all_7col.bedpe'
columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'type']

# Suggest reading one row to verify column count, or read directly
df = pd.read_csv(file_path, sep=r'\s+', header=None, names=columns, engine='python')

# Verify the data load
print(f"Total rows read: {len(df)}")
if len(df) > 0:
    print("First row 'type' value:", repr(df['type'].iloc[0]))

# 2. Preprocess: strip hidden spaces in 'type' and lowercase
df['type'] = df['type'].str.strip()

# simple comment
convergent_df = df[df['type'] == 'convergent']
print(f"Rows after filtering 'convergent': {len(convergent_df)}")

if not convergent_df.empty:
    # Extract anchor1 and anchor2
    anchor1 = convergent_df[['chrom1', 'start1', 'end1']].rename(columns={'chrom1': 'chrom', 'start1': 'start', 'end1': 'end'})
    anchor2 = convergent_df[['chrom2', 'start2', 'end2']].rename(columns={'chrom2': 'chrom', 'start2': 'start', 'end2': 'end'})

    # Merge and deduplicate
    anchors_bed = pd.concat([anchor1, anchor2]).drop_duplicates().sort_values(['chrom', 'start'])

    # Save as BED
    anchors_bed.to_csv('K562.CTCF-loop-anchors.bed', sep='\t', index=False, header=False)
    print("Success: File saved.")
else:
    print("Warning: No 'convergent' rows found. Output file not created.")