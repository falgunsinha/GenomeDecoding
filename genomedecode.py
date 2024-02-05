!pip install bio

## Download & Convert Fasta files to CSV files

import os
import requests
import gzip
from io import BytesIO
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import csv
import math

# Function to calculate GC1, GC2, GC3
def calculate_codon_gc(sequence):
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    gc1 = gc_fraction("".join(codons[::3])) * 100
    gc2 = gc_fraction("".join(codons[1::3])) * 100
    gc3 = gc_fraction("".join(codons[2::3])) * 100
    return gc1, gc2, gc3

# Function to calculate ENc
def calculate_enc(sequence):
    codon_count = {}
    total_codons = len(sequence) // 3

    for i in range(0, len(sequence) - 3, 3):
        codon = sequence[i:i + 3]
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            codon_count[codon] = 1

    f_values = list(codon_count.values())
    f_sum = sum(f_values)
    f_sum_squared = sum([x * x for x in f_values])

    if f_sum > 1:
        en_c = 2 / ((total_codons - 1) * (1 / f_sum_squared))
    else:
        en_c = math.nan

    return en_c

# Download and unzip function
def download_and_unzip(url, output_directory):
    response = requests.get(url)
    if response.status_code == 200:
        file_name = url.split('/')[-1]
        file_path = os.path.join(output_directory, file_name.replace('.gz', ''))
        with gzip.open(BytesIO(response.content), 'rb') as f_in, open(file_path, 'wb') as f_out:
            f_out.write(f_in.read())

        return file_path
    else:
        print(f"Failed to download {url}")
        return None

# Function to parse fasta file
def parse_fasta_file(input_file_path, output_file_path, organism_group, organism):
    with open(output_file_path, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['cds ID', 'cds length', 'Organism', 'Organism Group', 'GC%', 'ENc', 'GC1', 'GC2', 'GC3'])

        with open(input_file_path, "r") as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                cds_id = record.id.split()[0]
                cds_length = len(record.seq)
                gc_percent = gc_fraction(str(record.seq)) * 100
                enc_score = calculate_enc(str(record.seq))
                gc1, gc2, gc3 = calculate_codon_gc(str(record.seq))

                writer.writerow([cds_id, cds_length, organism, organism_group, gc_percent, enc_score, gc1, gc2, gc3])

# Organism data - Organism Group, Organism, and URL
organism_data = [
    ('Mammalia', 'Pongo_abelii', 'https://ftp.ensembl.org/pub/release-108/fasta/pongo_abelii/cds/Pongo_abelii.Susie_PABv2.cds.all.fa.gz'),
    ('Mammalia', 'Sarcophilus_harrisii', 'https://ftp.ensembl.org/pub/release-108/fasta/sarcophilus_harrisii/cds/Sarcophilus_harrisii.mSarHar1.11.cds.all.fa.gz'),
    ('Mammalia', 'Ornithorthynchus_anatinus', 'https://ftp.ensembl.org/pub/release-108/fasta/ornithorhynchus_anatinus/cds/Ornithorhynchus_anatinus.mOrnAna1.p.v1.cds.all.fa.gz'),
    ('Mammalia', 'Ursus_maritimus', 'https://ftp.ensembl.org/pub/release-108/fasta/ursus_maritimus/cds/Ursus_maritimus.UrsMar_1.0.cds.all.fa.gz'),
    ('Mammalia', 'Rattus_norvegicus', 'https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/cds/Rattus_norvegicus.mRatBN7.2.cds.all.fa.gz'),
    ('Aves', 'Apteryx_haastii', 'https://ftp.ensembl.org/pub/release-108/fasta/apteryx_haastii/cds/Apteryx_haastii.aptHaa1.cds.all.fa.gz'),
    ('Aves', 'Parus_major', 'https://ftp.ensembl.org/pub/release-108/fasta/parus_major/cds/Parus_major.Parus_major1.1.cds.all.fa.gz'),
    ('Aves', 'Otus_sunia', 'https://ftp.ensembl.org/pub/release-108/fasta/otus_sunia/cds/Otus_sunia.OtuSun1.0.cds.all.fa.gz'),
    ('Aves', 'Anas_platyrhynchos', 'https://ftp.ensembl.org/pub/release-108/fasta/anas_platyrhynchos/cds/Anas_platyrhynchos.ASM874695v1.cds.all.fa.gz'),
    ('Aves', 'Numida_meleagris', 'https://ftp.ensembl.org/pub/release-108/fasta/numida_meleagris/cds/Numida_meleagris.NumMel1.0.cds.all.fa.gz'),
    ('Molusca', 'Octupus_bimaculoides', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/octopus_bimaculoides_gca001194135v1/cds/Octopus_bimaculoides_gca001194135v1.Octopus_bimaculoides_v2_0.cds.all.fa.gz'),
    ('Molusca', 'Lottia_gigantea', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/lottia_gigantea/cds/Lottia_gigantea.Lotgi1.cds.all.fa.gz'),
    ('Molusca', 'Haliotis_rubra', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/haliotis_rubra_gca003918875v1rs/cds/Haliotis_rubra_gca003918875v1rs.ASM391887v1.cds.all.fa.gz'),
    ('Molusca', 'Crassostrea_gigas', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/crassostrea_gigas/cds/Crassostrea_gigas.cgigas_uk_roslin_v1.cds.all.fa.gz'),
    ('Molusca', 'Mizuhopecten_yessoensis', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/mizuhopecten_yessoensis_gca002113885v2/cds/Mizuhopecten_yessoensis_gca002113885v2.ASM211388v2.cds.all.fa.gz'),
    ('Osteichthyes', 'Anabas_testudineous', 'https://ftp.ensembl.org/pub/release-108/fasta/anabas_testudineus/cds/Anabas_testudineus.fAnaTes1.2.cds.all.fa.gz'),
    ('Osteichthyes', 'Cyprinus_carpio_carpio', 'https://ftp.ensembl.org/pub/release-108/fasta/cyprinus_carpio_carpio/cds/Cyprinus_carpio_carpio.Cypcar_WagV4.0.cds.all.fa.gz'),
    ('Osteichthyes', 'Cottoperca_gobio', 'https://ftp.ensembl.org/pub/release-108/fasta/cottoperca_gobio/cds/Cottoperca_gobio.fCotGob3.1.cds.all.fa.gz'),
    ('Osteichthyes', 'Esox_lucius', 'https://ftp.ensembl.org/pub/release-108/fasta/esox_lucius/cds/Esox_lucius.Eluc_v4.cds.all.fa.gz'),
    ('Osteichthyes', 'Denticeps_clupeoides', 'https://ftp.ensembl.org/pub/release-108/fasta/denticeps_clupeoides/cds/Denticeps_clupeoides.fDenClu1.1.cds.all.fa.gz'),
    ('Crustacea', 'Penaeus_monodon', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/penaeus_monodon_gca015228065v1/cds/Penaeus_monodon_gca015228065v1.NSTDA_Pmon_1.cds.all.fa.gz'),
    ('Crustacea', 'Daphnia_magna', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/daphnia_magna_gca020631705v2/cds/Daphnia_magna_gca020631705v2.ASM2063170v1.1.cds.all.fa.gz'),
    ('Crustacea', 'Pollicipes_pollicipes', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/pollicipes_pollicipes_gca011947565v2/cds/Pollicipes_pollicipes_gca011947565v2.Ppol_2.cds.all.fa.gz'),
    ('Crustacea', 'Hyalella_azteca', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/hyalella_azteca_gca000764305v2/cds/Hyalella_azteca_gca000764305v2.Hazt_2.0.cds.all.fa.gz'),
    ('Crustacea', 'Lepeophtheirus_salmonis', 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/fasta/lepeophtheirus_salmonis/cds/Lepeophtheirus_salmonis.LSalAtl2s.cds.all.fa.gz')
]

# Function to process organism data
def process_organism_data(organism_data, output_directory):
    for organism_group, organism, url in organism_data:
        fasta_file_name = url.split('/')[-1]
        print(f"Processing Organism Group: {organism_group}, Organism: {organism}, Fasta File: {fasta_file_name}")

        file_path = download_and_unzip(url, output_directory)
        if file_path:
            output_file_name = os.path.basename(file_path).replace('.fa', '.csv')
            output_csv_path = os.path.join(output_directory, output_file_name)
            parse_fasta_file(file_path, output_csv_path, organism_group, organism)
            csv_file_name = os.path.basename(output_csv_path)
            print(f"Generating Organism Group: {organism_group}, Organism: {organism}, CSV File: {csv_file_name}")

output_directory = '/content/drive/MyDrive/FE_Data'
process_organism_data(organism_data, output_directory)

## Combine all CSV files

import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Folder path containing all 25 CSV files
folder_path = '/content/drive/MyDrive/FE_Data'

# Combined data frame
all_combined_data = pd.DataFrame()

# Combine CSVs
for file_name in os.listdir(folder_path):
    if file_name.endswith('.csv'):
        file_path = os.path.join(folder_path, file_name)
        data = pd.read_csv(file_path)
        all_combined_data = pd.concat([all_combined_data, data])

## K-Means & Scatter Plots with ENc, GC1, GC2 & GC3

# Normalize data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(all_combined_data[['ENc', 'GC1', 'GC2', 'GC3']])

# Perform K-means clustering based on 5 clusters (Organism Groups)
kmeans = KMeans(n_clusters=5, random_state=0, n_init=1)
cluster_labels = kmeans.fit_predict(scaled_data)

# Add the cluster labels
all_combined_data['Cluster'] = cluster_labels

# 2D Scatter Plot for Organism Groups
plt.figure(figsize=(12, 9))
for group_name in all_combined_data['Organism Group'].unique():
    group_data = all_combined_data[all_combined_data['Organism Group'] == group_name]
    plt.scatter(group_data['ENc'], group_data['GC3'], label=f'{group_name}')
plt.xlabel('ENc')
plt.ylabel('GC3')
plt.title('2D Scatter Plot: ENc vs GC3 for Organism Groups')
plt.legend(title='Organism Groups')
plt.grid(True)
plt.show()

# 3D Scatter Plot for Organism Groups
fig = plt.figure(figsize=(18, 12))
ax = fig.add_subplot(111, projection='3d')
for group_name in all_combined_data['Organism Group'].unique():
    group_data = all_combined_data[all_combined_data['Organism Group'] == group_name]
    ax.scatter(group_data['GC1'], group_data['GC2'], group_data['GC3'], label=f'{group_name}')
ax.set_xlabel('GC1')
ax.set_ylabel('GC2')
ax.set_zlabel('GC3', labelpad=-3)
plt.title('3D Scatter Plot: with GC1, GC2, GC3 for Organism Groups')
ax.legend(title='Organism Groups')
plt.show()

## K-Means & Scatter Plots with GC%, ENc & GC1

# Normalize data using StandardScaler for all combined data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(all_combined_data[['GC%', 'ENc', 'GC1', 'GC2', 'GC3']])

# Perform K-means clustering based on 5 clusters (Organism Groups)
kmeans = KMeans(n_clusters=5, random_state=0, n_init=1)
cluster_labels = kmeans.fit_predict(scaled_data)  # Clustering

# Adding the cluster labels to the combined data
all_combined_data['Cluster'] = cluster_labels

# Plotting 2D Scatter Plot for all Organism Groups together
plt.figure(figsize=(12, 9))
for group_name in all_combined_data['Organism Group'].unique():
    group_data = all_combined_data[all_combined_data['Organism Group'] == group_name]
    plt.scatter(group_data['GC%'], group_data['ENc'], label=f'{group_name}')
plt.xlabel('GC%')
plt.ylabel('ENc')
plt.title('2D Scatter Plot: GC% vs ENc for Organism Groups')
plt.legend(title='Organism Groups')
plt.grid(True)
plt.show()

# Plotting 3D Scatter Plot with adjusted axes ranges
fig = plt.figure(figsize=(18, 12))
ax = fig.add_subplot(111, projection='3d')
for group_name in all_combined_data['Organism Group'].unique():
    group_data = all_combined_data[all_combined_data['Organism Group'] == group_name]
    ax.scatter(group_data['GC%'], group_data['ENc'], group_data['GC1'], label=f'{group_name}')
ax.set_xlabel('GC%')
ax.set_ylabel('ENc')
ax.set_zlabel('GC1', labelpad=-3)
plt.title('3D Scatter Plot: with GC%, ENc, GC1 for Organism Groups')
ax.legend(title='Organism Groups')
plt.show()

## PCA & Scatter Plots

from sklearn.decomposition import PCA

# Normalize data using StandardScaler for all combined data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(all_combined_data[['ENc', 'GC1', 'GC2', 'GC3']])

# Perform PCA
pca = PCA(n_components=3)
principal_components = pca.fit_transform(scaled_data)

# Adding the PCA components to the combined data
all_combined_data['PC1'] = principal_components[:, 0]
all_combined_data['PC2'] = principal_components[:, 1]
all_combined_data['PC3'] = principal_components[:, 2]

# Plotting 2D Scatter Plot for all Organism Groups together based on PCA components
plt.figure(figsize=(12, 9))
for group_name in all_combined_data['Organism Group'].unique():
    group_data = all_combined_data[all_combined_data['Organism Group'] == group_name]
    plt.scatter(group_data['PC1'], group_data['PC2'], label=f'{group_name}')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('2D Scatter Plot: PC1 vs PC2 for Organism Groups')
plt.legend(title='Organism Groups')
plt.grid(True)
plt.show()

# Plotting 3D Scatter Plot with adjusted axes ranges based on PCA components
fig = plt.figure(figsize=(18, 12))
ax = fig.add_subplot(111, projection='3d')
for group_name in all_combined_data['Organism Group'].unique():
    group_data = all_combined_data[all_combined_data['Organism Group'] == group_name]
    ax.scatter(group_data['PC1'], group_data['PC2'], group_data['PC3'], label=f'{group_name}')
ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')
ax.set_zlabel('Principal Component 3', labelpad=-3)
plt.title('3D Scatter Plot: with PC1, PC2, PC3 for Organism Groups')
ax.legend(title='Organism Groups')
plt.show()

## Drawing Correlation Matrix for GC%, ENc, GC1, GC2, GC3

import seaborn as sns

# Normalize data using StandardScaler for all combined data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(all_combined_data[['GC%', 'ENc', 'GC1', 'GC2', 'GC3']])

# Calculate the correlation matrix
correlation_matrix = all_combined_data[['GC%', 'ENc', 'GC1', 'GC2', 'GC3']].corr()

# Plotting the correlation matrix as a heatmap with adjusted range
plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1,
            fmt='.2f', annot_kws={"size": 10})
plt.title('Correlation Matrix')
plt.show()

