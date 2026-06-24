import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import os
import json


class DataPreprocessor:
    """Handles all data preprocessing operations for gene expression data."""

    def __init__(self):
        self.summary = {}
        self.scaler = None

    def load_data(self, filepath):
        """Load gene expression data from CSV file."""
        df = pd.read_csv(filepath, index_col=0) if 'index_col' in open(filepath).readline() else pd.read_csv(filepath)
        if df.columns[0] == '' or 'Unnamed' in df.columns[0]:
            df = pd.read_csv(filepath, index_col=0)
        self.raw_data = df.copy()
        self.summary['original_shape'] = list(df.shape)
        self.summary['original_genes'] = int(df.shape[0])
        self.summary['original_samples'] = int(df.shape[1])
        return df

    def handle_missing_values(self, df, strategy='mean'):
        """Handle missing values using specified strategy."""
        missing_before = int(df.isnull().sum().sum())
        if strategy == 'mean':
            df = df.fillna(df.mean())
        elif strategy == 'median':
            df = df.fillna(df.median())
        elif strategy == 'zero':
            df = df.fillna(0)
        elif strategy == 'drop':
            df = df.dropna()
        missing_after = int(df.isnull().sum().sum())
        self.summary['missing_before'] = missing_before
        self.summary['missing_after'] = missing_after
        self.summary['missing_handled'] = missing_before - missing_after
        return df

    def filter_genes(self, df, min_expression=1.0, min_samples=0.5):
        """Filter lowly expressed genes."""
        threshold_samples = int(df.shape[1] * min_samples)
        mask = (df > min_expression).sum(axis=1) >= threshold_samples
        genes_before = int(df.shape[0])
        df = df.loc[mask]
        genes_after = int(df.shape[0])
        self.summary['genes_before_filter'] = genes_before
        self.summary['genes_after_filter'] = genes_after
        self.summary['genes_filtered_out'] = genes_before - genes_after
        return df

    def normalize_data(self, df, method='log2'):
        """Normalize gene expression data."""
        if method == 'log2':
            df = df.copy()
            df = df.replace(0, 0.001)
            df = np.log2(df)
        elif method == 'quantile':
            rank_mean = df.stack().groupby(
                df.rank(method='first').stack().astype(int)
            ).mean()
            df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
        self.summary['normalization_method'] = method
        return df

    def scale_data(self, df, method='standard'):
        """Scale gene expression data."""
        if method == 'standard':
            self.scaler = StandardScaler()
            scaled_values = self.scaler.fit_transform(df.T).T
            df_scaled = pd.DataFrame(scaled_values, index=df.index, columns=df.columns)
        elif method == 'minmax':
            self.scaler = MinMaxScaler()
            scaled_values = self.scaler.fit_transform(df.T).T
            df_scaled = pd.DataFrame(scaled_values, index=df.index, columns=df.columns)
        self.summary['scaling_method'] = method
        return df_scaled

    def run_pipeline(self, filepath, missing_strategy='mean',
                     min_expression=1.0, min_samples=0.5,
                     normalization='log2', scaling='standard'):
        """Run complete preprocessing pipeline."""
        df = self.load_data(filepath)
        df = self.handle_missing_values(df, missing_strategy)
        df = self.filter_genes(df, min_expression, min_samples)
        df = self.normalize_data(df, normalization)
        df = self.scale_data(df, scaling)
        self.processed_data = df
        self.summary['final_shape'] = list(df.shape)
        self.summary['final_genes'] = int(df.shape[0])
        self.summary['final_samples'] = int(df.shape[1])
        return df

    def get_summary(self):
        """Get preprocessing summary as JSON-serializable dict."""
        return self.summary

    def save_processed_data(self, filepath):
        """Save processed data to CSV."""
        self.processed_data.to_csv(filepath)
        return filepath
