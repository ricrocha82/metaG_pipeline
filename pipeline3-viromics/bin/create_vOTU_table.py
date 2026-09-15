#!/usr/bin/env python3

# modified from MVP (https://gitlab.com/ccoclet/mvp)

import pandas as pd
import glob
import os
import argparse
from functools import reduce

def create_votu_abundance_table(coverm_dir, output_dir, covered_fraction_thresholds=None, 
                               normalization='Trimmed Mean', quality_file=None):
    """
    Create vOTU abundance tables from CoverM outputs
    
    Parameters:
    coverm_dir: Directory containing CoverM TSV files
    output_dir: Directory to save output tables
    covered_fraction_thresholds: List of coverage thresholds to apply
    normalization: Normalization method (Trimmed Mean, RPKM or TPM)
    quality_file: Optional CheckV quality summary file
    """
    
    print("Step 1: Reading and merging CoverM files...")
    
    # Find all CoverM files
    coverm_files = glob.glob(os.path.join(coverm_dir, "*.tsv"))
    if not coverm_files:
        coverm_files = glob.glob(os.path.join(coverm_dir, "*", "*.tsv"))
    
    print(f"Found {len(coverm_files)} CoverM files")
    
    dfs = []
    for filename in coverm_files:
        print(f"Processing: {os.path.basename(filename)}")
        
        # Read CoverM file
        df = pd.read_csv(filename, sep='\t')
        
        # Select relevant columns
        required_cols = ['Sample', 'Contig', 'Covered Fraction', 'Length', normalization]
        df = df[required_cols]
        
        # Remove duplicates (keep last)
        df = df.drop_duplicates(subset=['Sample', 'Contig', 'Length'], keep='last')
        
        # Pivot to have samples as columns
        df_pivot = df.pivot(index=['Contig', 'Length'], columns='Sample')
        
        # Flatten column names
        df_pivot.columns = ['_'.join(str(col) for col in cols).strip() for cols in df_pivot.columns.values]
        
        dfs.append(df_pivot)
    
    # Merge all dataframes
    print("Step 2: Merging all CoverM tables...")
    if len(dfs) > 1:
        df_merged = reduce(lambda left, right: pd.merge(left, right, on=['Contig', 'Length'], how='outer'), dfs)
    else:
        df_merged = dfs[0]
    
    df_merged = df_merged.reset_index()
    
    # Fill NaN values with 0
    df_merged = df_merged.fillna(0)
    
    # Merge with quality data if provided
    if quality_file and os.path.exists(quality_file):
        print("Step 3: Merging with quality data...")
        quality_df = pd.read_csv(quality_file, sep='\t')
        df_merged = pd.merge(quality_df, df_merged, left_on='virus_id', right_on='Contig', how='right')
        df_merged = df_merged.loc[:, ~df_merged.columns.duplicated(keep='last')]
        # Remove redundant columns
        if 'Contig' in df_merged.columns and 'virus_id' in df_merged.columns:
            df_merged = df_merged.drop(columns=['Contig'])
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Save unfiltered table
    unfiltered_output = os.path.join(output_dir, 'Unfiltered_vOTU_table.tsv')
    df_unfiltered = df_merged.loc[:, ~df_merged.columns.str.contains('Covered Fraction')]
    df_unfiltered = df_unfiltered.rename(columns = lambda x: x.strip(f'{normalization}_'))
    df_unfiltered.to_csv(unfiltered_output, sep='\t', index=False)
    print(f"Unfiltered table saved: {unfiltered_output}")
    
    print("Step 4: Creating coverage-filtered tables...")
    
    # Create filtered tables for each coverage threshold
    if not covered_fraction_thresholds:
        print("Skipping coverage fraction filtering (no thresholds provided).")
        print("\nvOTU abundance tables created successfully!")
        return df_unfiltered
    else:
        for cf in covered_fraction_thresholds:
            print(f"Processing coverage threshold: {cf}")
            
            # Copy the dataframe
            df_filtered = df_merged.copy()
            
            # Get coverage and abundance columns
            covered_fraction_cols = [col for col in df_filtered.columns if 'Covered Fraction' in col]
            abundance_cols = [col for col in df_filtered.columns if normalization in col]
            
            # Apply coverage filtering
            for i, cf_col in enumerate(covered_fraction_cols):
                # Find corresponding abundance column
                sample_name = cf_col.replace('Covered Fraction_', '')
                abundance_col = f'{normalization}_{sample_name}'
                
                if abundance_col in df_filtered.columns:
                    # Set abundance to 0 where coverage is below threshold
                    mask = df_filtered[cf_col] < cf
                    df_filtered.loc[mask, abundance_col] = 0
            
            # Remove rows where all abundance values are 0
            if abundance_cols:
                non_zero_mask = (df_filtered[abundance_cols] > 0).any(axis=1)
                df_filtered = df_filtered[non_zero_mask]
            
            # Remove coverage fraction columns from final output
            df_filtered = df_filtered.loc[:, ~df_filtered.columns.str.contains('Covered Fraction')]
            df_filtered_renamed = df_filtered.rename(columns = lambda x: x.strip(f'{normalization}_'))
            
            # Save filtered table
            filtered_output = os.path.join(output_dir, f'Filtered_CF_{cf}_vOTU_table.tsv')
            df_filtered_renamed.to_csv(filtered_output, sep='\t', index=False)
            print(f"Filtered table (CF={cf}) saved: {filtered_output}")
            
            # Print summary
            print(f"  - vOTUs retained: {len(df_filtered)}")
            if abundance_cols:
                total_abundance = df_filtered[abundance_cols].sum().sum()
                print(f"  - Total abundance: {total_abundance:.2f}")
        
            print("\nvOTU abundance tables created successfully!")
            return df_unfiltered

def apply_quality_filtering(df, filtration='conservative', viral_min_genes=1, host_viral_genes_ratio=1):
    """
    Apply quality-based filtering to the vOTU table
    """
    print(f"Applying {filtration} quality filtering...")
    
    # Check if required columns exist
    required_cols = ['viral_genes']
    if not all(col in df.columns for col in required_cols):
        print("Warning: Quality filtering columns not found. Skipping quality filtering.")
        return df
    
    # Apply basic filtering
    df_filtered = df.copy()
    
    # Filter by viral genes
    if 'viral_genes' in df_filtered.columns:
        df_filtered = df_filtered[df_filtered['viral_genes'] >= viral_min_genes]
    
    # Filter by host/viral genes ratio
    if 'host_genes' in df_filtered.columns and 'viral_genes' in df_filtered.columns:
        mask = df_filtered['host_genes'] / df_filtered['viral_genes'] <= host_viral_genes_ratio
        df_filtered = df_filtered[mask]
    
    # Apply conservative filtering if specified
    if filtration == 'conservative':
        if 'checkv_quality' in df_filtered.columns and 'virus_length' in df_filtered.columns:
            # Keep high quality sequences OR sequences >= 5kb
            quality_mask = df_filtered['checkv_quality'].isin(['High-quality', 'Medium-quality', 'Complete'])
            length_mask = df_filtered['virus_length'] >= 5000
            
            if 'completeness_method' in df_filtered.columns:
                completeness_mask = (df_filtered['completeness_method'].str.contains('AAI-based', na=False) | 
                                   df_filtered['completeness_method'].str.contains('DTR', na=False))
                conservative_mask = (completeness_mask & quality_mask) | length_mask
            else:
                conservative_mask = quality_mask | length_mask
            
            df_filtered = df_filtered[conservative_mask]
    
    print(f"Quality filtering: {len(df)} -> {len(df_filtered)} vOTUs retained")
    return df_filtered

def main():
    parser = argparse.ArgumentParser(description='Create vOTU abundance tables from CoverM outputs')
    parser.add_argument('-i', '--input_dir', required=True, 
                        help='Directory containing CoverM TSV files')
    parser.add_argument('-o', '--output_dir', required=True,
                    help='Output directory for vOTU tables')
    parser.add_argument('-q', '--quality_file', 
                    help='CheckV quality summary file (optional)')
    parser.add_argument('--covered_fraction', nargs='+', type=float, 
                    default=[],
                    help='Coverage fraction threshold (default: no filtering)')
    parser.add_argument('--normalization', choices=['RPKM', 'TPM', "Trimmed Mean"], 
                    default='Trimmed Mean',
                    help='Normalization method (default: Trimmed Mean)')
    parser.add_argument('--filtration', choices=['relaxed', 'conservative'], 
                    default='conservative',
                    help='Quality filtration level (default: conservative)')
    
    args = parser.parse_args()
    
    # Create vOTU abundance tables
    df_merged = create_votu_abundance_table(
        coverm_dir=args.input_dir,
        output_dir=args.output_dir,
        covered_fraction_thresholds=args.covered_fraction,
        normalization=args.normalization,
        quality_file=args.quality_file
    )

    # Normalize
    
    
    # Apply quality filtering if quality file is provided
    if args.quality_file and os.path.exists(args.quality_file):
        df_filtered = apply_quality_filtering(df_merged, filtration=args.filtration)
        
        # Save quality-filtered table
        quality_output = os.path.join(args.output_dir, f'Quality_Filtered_{args.filtration}_{args.normalization}_vOTU_table.tsv')
        df_filtered.to_csv(quality_output, sep='\t', index=False)
        print(f"Quality-filtered table saved: {quality_output}")


if __name__ == "__main__":
    main()