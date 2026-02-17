#!/usr/bin/env python3
"""
Script to parse BEAST tip-dating log files and summarize age estimates.
Parsing results from log files:
- joint (Mean, Stdev, HPD_95_Lower, HPD_95_Upper, ESS)
- prior (Mean, Stdev, HPD_95_Lower, HPD_95_Upper, ESS)
- likelihood (Mean, Stdev, HPD_95_Lower, HPD_95_Upper, ESS)
- skygrid.logPopSize1 (Mean, Stdev, HPD_95_Lower, HPD_95_Upper, ESS)
- age(root) (Mean, Stdev, HPD_95_Lower, HPD_95_Upper, ESS)
- clock.rate (Mean, Stdev, HPD_95_Lower, HPD_95_Upper, ESS)
- age(sample_id_ND) (Mean, Stdev, HPD_95_Lower, HPD_95_Upper, ESS)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
import argparse


class BeastLogParser:
    def __init__(self, log_file, burnin_percent=10):
        """
        Initialize the parser for a single BEAST log file.
        
        Args:
            log_file: Path to the log file
            burnin_percent: Percentage of chain to discard as burnin (default: 10%)
        """
        self.log_file = Path(log_file)
        self.burnin_percent = burnin_percent
        self.data = None
        self.data_post_burnin = None
        
    def calculate_hpd_95(self, values):
        """
        Calculate 95% Highest Posterior Density interval.
        
        Args:
            values: Array of MCMC samples
            
        Returns:
            Tuple of (lower, upper) bounds
        """
        sorted_values = np.sort(values)
        n = len(sorted_values)
        interval_size = int(n * 0.95)
        
        # Find the narrowest interval containing 95% of samples
        min_width = float('inf')
        best_interval = (sorted_values[0], sorted_values[-1])
        
        for i in range(n - interval_size + 1):
            width = sorted_values[i + interval_size - 1] - sorted_values[i]
            if width < min_width:
                min_width = width
                best_interval = (sorted_values[i], sorted_values[i + interval_size - 1])
        
        return best_interval
    
    def calculate_ess(self, values):
        """
        Calculate Effective Sample Size using integrated autocorrelation time.
        This method matches Tracer's implementation.
        
        Args:
            values: Array of MCMC samples
            
        Returns:
            ESS value
        """
        n = len(values)
        if n < 2:
            return 0.0
        
        mean_val = np.mean(values)
        variance = np.var(values, ddof=1)
        
        if variance == 0:
            return 0.0
        
        # Calculate autocorrelation at multiple lags
        max_lag = min(n - 1, int(n / 3))  # Don't go beyond n/3 lags
        autocorr_sum = 0.0
        
        for lag in range(1, max_lag):
            # Calculate autocorrelation at this lag
            if n - lag < 1:
                break
                
            autocovariance = np.sum((values[:-lag] - mean_val) * (values[lag:] - mean_val)) / n
            autocorr = autocovariance / variance
            
            # Stop if autocorrelation becomes negative (standard practice)
            if autocorr < 0:
                break
            
            autocorr_sum += autocorr
            
            # Early stopping if autocorrelation becomes very small
            if autocorr < 0.05:
                break
        
        # Calculate ESS using integrated autocorrelation time
        # ESS = n / (1 + 2 * sum(rho_k))
        act = 1.0 + 2.0 * autocorr_sum  # Autocorrelation time
        ess = n / act
        
        return max(ess, 1.0)
    
    def calculate_stats(self, column_name):
        """
        Calculate statistics for a given parameter.
        
        Args:
            column_name: Name of the column/parameter
            
        Returns:
            Dictionary with statistics
        """
        values = self.data_post_burnin[column_name].values
        
        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1)
        hpd_lower, hpd_upper = self.calculate_hpd_95(values)
        ess_val = self.calculate_ess(values)
        
        return {
            'Parameter': column_name,
            'Mean': mean_val,
            'Stdev': std_val,
            'HPD_95_Lower': hpd_lower,
            'HPD_95_Upper': hpd_upper,
            'ESS': ess_val
        }
    
    def parse_log(self):
        """
        Parse the log file and extract all relevant parameters.
        
        Returns:
            DataFrame with statistics for all parameters
        """
        print(f"Reading log file: {self.log_file}")
        
        # Read log file (skip comment lines starting with #)
        self.data = pd.read_csv(self.log_file, sep='\t', comment='#')
        
        # Apply burnin
        burnin_index = int(len(self.data) * (self.burnin_percent / 100))
        self.data_post_burnin = self.data.iloc[burnin_index:].reset_index(drop=True)
        
        print(f"Total samples: {len(self.data)}, After {self.burnin_percent}% burnin: {len(self.data_post_burnin)}")
        
        results = []
        
        # Parse key parameters
        key_params = ['joint', 'prior', 'likelihood', 'age(root)', 'clock.rate', 'skygrid.logPopSize1']
        for param in key_params:
            if param in self.data_post_burnin.columns:
                results.append(self.calculate_stats(param))
                print(f"  Processed: {param}")
            else:
                print(f"  Warning: {param} not found in log file")
        
        # Find all age(..._ND) columns
        age_nd_pattern = re.compile(r'age\((.+_ND)\)')
        age_nd_columns = []
        
        for col in self.data_post_burnin.columns:
            match = age_nd_pattern.match(col)
            if match:
                sample_id = match.group(1)  # Extract full ID with _ND suffix
                age_nd_columns.append((col, sample_id))
        
        print(f"\nFound {len(age_nd_columns)} sample age parameters with _ND suffix")
        
        # Process all age(_ND) parameters
        for col, sample_id in age_nd_columns:
            stats = self.calculate_stats(col)
            stats['Sample_ID'] = sample_id
            results.append(stats)
            print(f"  Processed: {sample_id}")
        
        df = pd.DataFrame(results)
        print(f"\nTotal parameters processed: {len(df)}")
        
        return df
    
    def save_results(self, df, output_file=None):
        """
        Save results to CSV and Excel files.
        
        Args:
            df: DataFrame with results
            output_file: Output filename (without extension). 
                        If None, saves to same directory as log file.
        """
        if output_file is None:
            output_file = self.log_file.parent / f"{self.log_file.stem}_summary"
        else:
            output_file = Path(output_file)
        
        # Save as CSV
        csv_path = output_file.with_suffix('.csv')
        df.to_csv(csv_path, index=False)
        print(f"\nResults saved to: {csv_path}")
        
        
        # Print summary
        summary_lines = [
            "",
            "=== Log Result Parsing Summary ===",
            f"Total parameters: {len(df)}",
        ]

        sample_df = df[df['Parameter'].str.contains(r'age\(.*_ND\)', regex=True, na=False)]
        if len(sample_df) > 0:
            summary_lines.extend([
                f"Sample ages with _ND suffix: {len(sample_df)}",
                "",
                "ESS Statistics for samples:",
                f"  Mean ESS: {sample_df['ESS'].mean():.2f}",
                f"  Min ESS: {sample_df['ESS'].min():.2f}",
                f"  Max ESS: {sample_df['ESS'].max():.2f}",
                f"  Samples with ESS < 200: {(sample_df['ESS'] < 200).sum()}",
            ])

        # Print summary to stdout
        for line in summary_lines:
            print(line)

        # Also save summary to a separate log file (text)
        summary_log_path = output_file.with_suffix('.log')
        with open(summary_log_path, 'w', encoding='utf-8') as fh:
            fh.write("\n".join(summary_lines) + "\n")
        print(f"Summary log saved to: {summary_log_path}")

    


def main():
    """
    Main function to run the parser.
    """
    parser = argparse.ArgumentParser(
        description='Parse BEAST log file and extract age estimates for samples with _ND suffix'
    )
    parser.add_argument('log_file', type=str, help='Path to BEAST log file')
    parser.add_argument('-b', '--burnin', type=int, default=10, 
                       help='Burnin percentage (default: 10)')
    parser.add_argument('-o', '--output', type=str, default=None,
                       help='Output file path (without extension). Default: log_file_summary')
    
    args = parser.parse_args()
    
    # Parse log file
    log_parser = BeastLogParser(args.log_file, burnin_percent=args.burnin)
    df_results = log_parser.parse_log()
    
    # Save results
    log_parser.save_results(df_results, output_file=args.output)
    
    # Display results
    #print("\n=== Results Preview ===")
    #print(df_results.to_string(index=False))


if __name__ == "__main__":
    main()