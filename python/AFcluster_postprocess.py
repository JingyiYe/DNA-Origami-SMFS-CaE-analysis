# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 22:02:55 2024

@author: ccana
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from matplotlib.colors import LinearSegmentedColormap

from polyleven import levenshtein

from pandarallel import pandarallel
from numpy import genfromtxt

pandarallel.initialize(use_memory_fs=False)

sns.set_style('ticks')
sns.set_context('paper')

def get_type(pdb):
    end = pdb.split('_')[-1]
    if end=='REF.pdb':
        return 'REF'
    elif end.startswith('U'):
        return end.split('-')[0]
    else:
        return 'Tree'

def read_b_factor(pdb_file):
    vals=[]

    with open(pdb_file,'r') as f:
        for lin in f.readlines()[1:-3]:
            fields = lin.split()
            if fields[2] =='CA':
                vals.append(float(fields[10]))
    return vals

def calc_dssp(decoy):
    if isinstance(decoy,str):
        decoy = md.load_pdb(decoy)

    return ''.join([x for x in md.compute_dssp(decoy).tolist()[0]])

def calc_rmsd(decoy, ref_obj, debug=False):
    if isinstance(decoy,str):
        decoy = md.load_pdb(decoy)

    ref_atoms = [x.index for x in ref_obj.top.residues]
    ref_seq = ref_obj.top.to_fasta()[0]
    decoy_seq = decoy.top.to_fasta()[0]

    string='name CA and ('+ ' or '.join("resid %d" % x for x in ref_atoms)+')'

    if debug:
        
        print(ref_seq)
        print(decoy_seq)
        if len(decoy.top.select(string)) != len(ref_obj.top.select(string)):
            print(len(decoy.top.select(string)), len(ref_obj.top.select(string)), string)

        print([x.index for x in ref_obj.top.residues])
        print([x.index for x in decoy.top.residues])

        print(ref_obj.top.select(string))
        print(decoy.top.select(string))
    return md.rmsd(decoy, ref_obj,  atom_indices=decoy.top.select(string), ref_atom_indices=ref_obj.top.select(string))[0]

    
#%%

#filename='C:/Users/ccana/Downloads/Rfah_3recycles.json.zip'
filename='C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/NCS1_234/rmsds.json.zip'
#filename='C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/NCS1_234/df_with_PCs.json.zip'
filename='C:\Users\ccana\Desktop\DNA Origami\AF-Cluster\CaE\output.zip'
df = pd.read_json(filename)


    
df['Type'] = df.apply(lambda row: get_type(row['pdb']), axis=1)

df = df.loc[df.Type=='Tree']

#%%
# Path to your PDB files
#CaE
#CaE 234
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaE234/CaE234_output/content/output'
#ref_struct= 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaE234/CaE234_1e9a1.result/CaE234_1e9a1/CaE234_1e9a1_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
#CaE 34
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaE34/CaE34_output/content/output'
#ref_struct= 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaE34/CaE34_792af.result/CaE34_792af/CaE34_792af_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
#CaE
folder_path ='C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaE/output/content/output'
ref_struct= 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaE/CaE_d1f8c.result/CaE_d1f8c/CaE_d1f8c_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'

#CaM
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaM/CaM_output/content/output'
#ref_struct = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaM/CaM_b5361.result/CaM_b5361/CaM_b5361_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb'
#CaM34
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaM34/CaM34_output/content/output'
#ref_struct = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaM34/CaM34_9e85b.result/CaM34_9e85b/CaM34_9e85b_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb'
#CaM123
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaM123/CaM123_output/content/output'
#ref_struct = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaM123/CaM123_a8a54.result/CaM123_a8a54/CaM123_a8a54_unrelaxed_rank_002_alphafold2_ptm_model_5_seed_000.pdb'

#NCS1
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/NCS1/NCS1_output/content/output'
#ref_struct = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/NCS1/NCS1_e4f4b.result/NCS1_e4f4b/NCS1_e4f4b_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
#NCS1_234
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/NCS1_234/NCS1_234_output/content/output'
#ref_struct= 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/NCS1_234/NCS1_234_1d0ea.result/NCS1_234_1d0ea/NCS1_234_1d0ea_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
#NCS1_34
#folder_path = 'C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/NCS1_34/NCS1_34_output/content/output'
#ref_struct= 'C:/Users/ccana\Desktop/DNA Origami/AF-Cluster/NCS1_34/NCS1_34_71bdf.result/NCS1_34_71bdf/NCS1_34_71bdf_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
#ref_struct = ref_struct.replace('\\', '/')


ref_obj = md.load_pdb(ref_struct)

# List all files in the folder
all_files = os.listdir(folder_path)
# Filter out .pdb files that do not start with 'EX-U'
selected_files = [file for file in all_files if file.endswith('.pdb') and not file.startswith('EX_U')]


pdb_files = sorted([os.path.join(folder_path, f) for f in selected_files if f.endswith('.pdb')])

print("Found %d pdbs" % len(pdb_files))

df = pd.DataFrame({'pdb': pdb_files})
#['pLDDT_vector'] = df.apply(lambda row: read_b_factor(row['pdb']), axis=1)

df['pLDDT_vector'] = df.apply(lambda row: read_b_factor(row['pdb']), axis=1)
df['mean_pLDDT'] = df.apply(lambda row: np.mean(row['pLDDT_vector']), axis=1)
#df['dssp_string'] = df.apply(lambda row: calc_dssp(row['pdb']), axis=1)
df['rmsd_ref'] = df.apply(lambda row: calc_rmsd(row['pdb'], ref_obj), axis=1)

        

# Load all PDB files as a trajectory
trajectory = md.load(pdb_files)
contacts = md.compute_contacts(trajectory, contacts='all', scheme='ca')

contact_distances, contact_indices = contacts

#%%
# Define the colors for the custom colormap (red to blue)
colors = ["red", "orange", "yellow", "cornflowerblue", "blue"]
n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins

# Create the colormap
cmap_name = 'confidence'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=100)

#%%
# Each row represents a frame from the trajectory, and each column represents a contact distance
contact_matrix = np.stack(contact_distances)

# Initialize PCA
pca = PCA(n_components=3)  # For example, let's reduce it to 3 components

# Fit PCA on the contact distance matrix
pca_result = pca.fit_transform(contact_matrix)

# pca_result now contains your principal components with shape (n_samples, n_components)
# You can access the explained variance ratio to see how much variance each principal component accounts for
explained_variance_ratio = pca.explained_variance_ratio_

# with plt.style.context(['science']): 
fig = plt.figure(figsize=(2,1.5))

ax = fig.add_subplot(111)

# Scatter plot using the first three PCA components
# X-axis: PC1, Y-axis: PC2, Z-axis: PC3
#ax.scatter(pca_result[:, 0], pca_result[:, 1], pca_result[:, 2])

plt.scatter(pca_result[:, 0], pca_result[:, 1], c= df['mean_pLDDT'], cmap=cm, vmin=0, vmax=100)
plot_xlabel='PC1 ('+ str(round(100*explained_variance_ratio[0],1))+ '%)'
plot_ylabel='PC2 ('+ str(round(100*explained_variance_ratio[1],1))+ '%)'


plt.xlabel(plot_xlabel)
plt.ylabel(plot_ylabel)
plt.colorbar(label='pLDDT Confidence Score')
plt.xlim(-100,100)
plt.ylim(-100,100)
#ax.set_zlabel('PC3')

plt.show()



#%%plt.figure(figsize=(4,4))



#%%
plt.figure(figsize=(2.5,1.25))
df = df.sort_values('mean_pLDDT')
# Using seaborn to create a scatter plot
#sns.scatterplot(x='rmsd_ref', y='mean_pLDDT', data=df, c=df['mean_pLDDT'], cmap=cm, vmin=0, vmax=100)
plt.scatter(df['rmsd_ref'], df['mean_pLDDT'], c= df['mean_pLDDT'], cmap=cm, vmin=0, vmax=100)

# Setting the x and y labels
plt.xlabel('r.m.s.d')
plt.ylabel('pLDDT')
plt.xlim(0,2)
plt.ylim(40,100)

#%%
plt.figure(figsize=(0.5,1.5))
data_zeros=np.zeros(len(df['rmsd_ref']))
#sns.violinplot(data=df, y="rmsd_ref",fill=False,color='#bababa',inner=None, linewidth=1.5, alpha=0)
# Create the violin plot
violin_parts = sns.violinplot(data=df, y="rmsd_ref", inner=None, color='#bababa', linewidth=1.5)

# Remove the fill
for violin in violin_parts.collections:
    violin.set_facecolor('none')
sns.stripplot(data=df, y="rmsd_ref", c= df['mean_pLDDT'], cmap=cm, alpha=0.7,s=2, vmin=0, vmax=100)
#    hue_order=[True, False], split=True)
#sns.kdeplot(data=df, y="rmsd_ref")
#plt.scatter(df['rmsd_ref'],data_zeros, c= df['mean_pLDDT'],s=40, cmap=cm, vmin=0, vmax=100)
plt.ylabel('r.m.s.d')
plt.ylim(0,2)
#%%

df.to_pickle('C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaE34/CaE34_data.pkl')
#savefig('rfah_landscape.pdf',bbox_inches='tight')

#%%
df = pd.read_pickle('C:/Users/ccana/Desktop/DNA Origami/AF-Cluster/CaM34/CaM34_data.pkl')

#%%
filename='C:/Users/ccana/Desktop/DNA Origami/CaE_contact_maps_table.csv'
#filename='C:/Users/ccana/Desktop/DNA Origami/CaM_contact_maps_table.csv'
#filename='C:/Users/ccana/Desktop/DNA Origami/NCS1_contact_maps_table.csv'
df = pd.read_csv(filename, sep=',')
plt.figure(figsize=(3,3))
#plt.imshow(df,cmap='RdYlBu_r',vmin=0, vmax=20)
plt.imshow(df,cmap='YlGnBu',vmin=0, vmax=20)
plt.colorbar(label='Distance (A)')