import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')


class DifferentialExpressionAnalyzer:
    """Performs differential gene expression analysis."""

    def __init__(self):
        self.results = None
        self.significant_genes = None
        self.summary = {}

    def compute_fold_change(self, df, condition_col, case_group, control_group):
        """Compute fold change between case and control groups."""
        case_samples = df.columns[df.columns.str.contains(case_group, case=False)]
        control_samples = df.columns[df.columns.str.contains(control_group, case=False)]

        case_mean = df[case_samples].mean(axis=1)
        control_mean = df[control_samples].mean(axis=1)

        fold_change = case_mean / control_mean.replace(0, 0.001)
        log2_fold_change = np.log2(fold_change.replace(0, 0.001))

        return case_mean, control_mean, fold_change, log2_fold_change

    def compute_p_values(self, df, condition_col, case_group, control_group):
        """Compute p-values using independent t-test."""
        case_samples = df.columns[df.columns.str.contains(case_group, case=False)]
        control_samples = df.columns[df.columns.str.contains(control_group, case=False)]

        p_values = []
        t_stats = []

        for idx in df.index:
            case_vals = df.loc[idx, case_samples].values.astype(float)
            control_vals = df.loc[idx, control_samples].values.astype(float)

            if len(case_vals) > 1 and len(control_vals) > 1:
                t_stat, p_val = stats.ttest_ind(case_vals, control_vals, equal_var=False)
            else:
                t_stat, p_val = 0, 1.0

            p_values.append(p_val)
            t_stats.append(t_stat)

        return np.array(p_values), np.array(t_stats)

    def adjust_p_values(self, p_values, method='fdr_bh'):
        """Adjust p-values for multiple testing."""
        reject, p_adjusted, _, _ = multipletests(p_values, method=method)
        return p_adjusted, reject

    def identify_significant_genes(self, log2_fc, p_adjusted, reject,
                                    lfc_threshold=1.0, p_threshold=0.05):
        """Identify significantly differentially expressed genes."""
        results_df = pd.DataFrame({
            'log2_fold_change': log2_fc,
            'p_adjusted': p_adjusted,
            'reject_null': reject
        })
        results_df['significant'] = (
            (np.abs(log2_fc) >= lfc_threshold) &
            (p_adjusted <= p_threshold) &
            reject
        )
        results_df['regulation'] = 'Not Significant'
        results_df.loc[
            (results_df['significant']) & (results_df['log2_fold_change'] > 0),
            'regulation'
        ] = 'Upregulated'
        results_df.loc[
            (results_df['significant']) & (results_df['log2_fold_change'] < 0),
            'regulation'
        ] = 'Downregulated'

        self.results = results_df
        self.significant_genes = results_df[results_df['significant']].copy()

        n_up = int((results_df['regulation'] == 'Upregulated').sum())
        n_down = int((results_df['regulation'] == 'Downregulated').sum())
        n_total = int(results_df['significant'].sum())

        self.summary = {
            'total_genes': int(len(results_df)),
            'significant_genes': n_total,
            'upregulated': n_up,
            'downregulated': n_down,
            'not_significant': int(len(results_df) - n_total),
            'lfc_threshold': lfc_threshold,
            'p_threshold': p_threshold
        }

        return results_df

    def run_analysis(self, df, condition_col='condition',
                     case_group='AD', control_group='Control',
                     lfc_threshold=1.0, p_threshold=0.05):
        """Run complete differential expression analysis."""
        case_mean, control_mean, fc, lfc = self.compute_fold_change(
            df, condition_col, case_group, control_group
        )
        p_values, t_stats = self.compute_p_values(
            df, condition_col, case_group, control_group
        )
        p_adjusted, reject = self.adjust_p_values(p_values)

        results = self.identify_significant_genes(
            lfc, p_adjusted, reject, lfc_threshold, p_threshold
        )

        results['case_mean'] = case_mean
        results['control_mean'] = control_mean
        results['fold_change'] = fc
        results['p_value'] = p_values
        results['t_statistic'] = t_stats

        self.results = results
        return results

    def get_top_genes(self, n=50, by='log2_fold_change'):
        """Get top N genes sorted by specified metric."""
        if self.results is None:
            return None
        sorted_results = self.results.sort_values(by, ascending=False)
        return sorted_results.head(n)

    def get_summary(self):
        """Get analysis summary."""
        return self.summary

    def export_results(self, filepath, significant_only=False):
        """Export results to CSV."""
        if significant_only and self.significant_genes is not None:
            self.significant_genes.to_csv(filepath)
        elif self.results is not None:
            self.results.to_csv(filepath)
        return filepath
