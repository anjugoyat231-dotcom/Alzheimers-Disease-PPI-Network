import pandas as pd
import numpy as np

np.random.seed(42)

n_genes = 1000
n_ad = 6
n_control = 6

ad_samples = [f'AD_{i+1}' for i in range(n_ad)]
control_samples = [f'Control_{i+1}' for i in range(n_control)]
all_samples = ad_samples + control_samples

genes = [f'GENE_{i+1}' for i in range(n_genes)]

data = np.random.lognormal(mean=2, sigma=0.8, size=(n_genes, len(all_samples)))
df = pd.DataFrame(data, index=genes, columns=all_samples)

ad_related = [
    'APP', 'PSEN1', 'PSEN2', 'BACE1', 'MAPT', 'APOE', 'GSK3B', 'CDK5',
    'BDNF', 'SNCA', 'CLU', 'BIN1', 'TREM2', 'ABCA1', 'ADAM10',
    'IL6', 'TNF', 'IL1B', 'TLR4', 'MYD88', 'NFKB1', 'RELA',
    'BCL2', 'BAX', 'CASP3', 'MTOR', 'SOD1', 'SOD2', 'CAT',
    'GRIN1', 'GRIN2B', 'DLG4', 'SYNGAP1', 'SHANK3', 'SNAP25', 'VAMP2',
    'CREB1', 'CAMK2A', 'HOMER1', 'ARC', 'EGR1', 'FOS', 'JUN',
    'HSPA1A', 'DNAJB1', 'SQSTM1', 'ATG5', 'ATG7', 'BECN1', 'PINK1',
    'NGF', 'NTRK1', 'NTRK2', 'MAPK1', 'MAPK3', 'AKT1', 'TP53',
    'PARK2', 'UCHL1', 'LRRK2', 'DYRK1A', 'PIN1', 'PP2A', 'FKBP4'
]

for i, gene in enumerate(ad_related):
    if gene in df.index:
        continue
    df.loc[gene] = np.random.lognormal(mean=2, sigma=0.8, size=len(all_samples))

for gene in ad_related[:20]:
    if gene in df.index:
        df.loc[gene, ad_samples] *= np.random.uniform(1.5, 3.0)

for gene in ad_related[20:35]:
    if gene in df.index:
        df.loc[gene, ad_samples] *= np.random.uniform(0.3, 0.7)

for gene in df.index:
    missing_cols = np.random.choice(len(all_samples), size=np.random.randint(0, 3), replace=False)
    df.loc[gene, all_samples[missing_cols]] = np.nan

df.to_csv('data/raw/alzheimer_sample_data.csv')
print(f"Sample dataset created: {df.shape[0]} genes, {df.shape[1]} samples")
print(f"AD-related genes included: {len(ad_related)}")
