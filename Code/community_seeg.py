# %%
import numpy as np
import scipy
import mat73
import os
import bct
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
# %%
DATA_DIR = '../Data/'
INP_F = 'Epat08_10_FAS_imcoh.mat'
BANDS = ['delta', 'theta', 'alpha','beta','gamma_l','gamma_h']
band_ind = dict(zip(BANDS,[i for i in range(6)]))
# %%

print(os.path.join(DATA_DIR,INP_F))
print("Hello world, how many modules are there today?")
conndata = mat73.loadmat(os.path.join(DATA_DIR,INP_F))

# %%
network_type = 'ImCoh'
conn_matrix = conndata["seizure"]["imcoh_all_windowed"]
np.shape(conn_matrix)
# results will show (num_freq_bin, num_windows, n_contact, n_contact)
# meaning that I have over 1000 connectivity matrices. 
# first interesting step would be to examine the connectivity of 
# one time slice and one freq. Let's pick gamma for a random
# middle time slice
# %%
BAND = 'alpha'
b_ind = band_ind[BAND]
slice_num = 100
bandSlice = conn_matrix[b_ind, slice_num,:,:]
sns.histplot(bandSlice.flatten())
plt.title(f"Histogram of {network_type} Values for {slice_num}th slice in {BAND}")

# %%
#Interestingly this is a skewed distribution with values mostly 
# concentrated around 0.1 with a fat tail an da bump at 0.5 and 0.4
# rich get richer?
# let's see a heatmap
bandSlice[np.isnan(bandSlice)] = .0000000001
sns.heatmap(bandSlice)
plt.title(f"Connectivity PDC  of {slice_num}  slice in {BAND} band")
# %%
# that is really weird, looks like all contacs are heavily linked to 
# 32 but 32 does not project so widely out. Wonder what the labels on here are
# Is 32 part of the soz? 
# let's get communities
from bct import BCTParamError
C = np.zeros((1000,42))
Q = np.zeros((1000,1))
ci = np.random.randint(1,42,(42,1))
param_errors = np.zeros((1000,1))
print(ci.shape)
gammas = np.linspace(0.01,3,1000)
#setup

# %% 
def consensus_modularity(W, ci, g,n=10):
    d = min(W.shape) # assume n > d here
    
    C = np.zeros((n,d))
    Q = np.zeros((n, 1))
    
    for i in range(n):
        ci,q = bct.community_louvain(W,gamma=g,ci=ci)
        C[i,:] = ci
        Q[i,:] = q
    W1 = np.zeros(W.shape)

    for i in range(n):
        ci = np.expand_dims(C[i,:],1)
        KKi = ci == ci.T
        W1 += KKi.astype(float)
    W1 = W1/n
    c_consennsus, q_consensus = bct.community_louvain(W1, gamma=g)
    return c_consennsus, q_consensus

for i,g in enumerate(gammas):
    try:
        
        ci, q = consensus_modularity(bandSlice, ci, g)
        C[i,:] = ci
        Q[i,:] = q
    except BCTParamError:
        param_errors[i] = g
        continue

num_errors = np.sum(param_errors > 0)
print(f"number of param errors: {num_errors}")
print(np.sum(~np.isnan(Q)))
# %%

gammas = np.linspace(0.01,3,1000)
valid_run_inds = np.where(param_errors==0)[0]

data = np.zeros((max(valid_run_inds.shape),2))
data[:,0] =  [len(set(C[i,:])) for i in valid_run_inds]

data[:,1] = [gammas[i] for i in valid_run_inds]
modularity_df = pd.DataFrame(data,\
                             columns=['N_modules', 'gamma_vals'])

# 
sns.lineplot(modularity_df, x='gamma_vals', y='N_modules')

# %%

modules = np.reshape(C[valid_run_inds[700], : ], (42,1))+1
mod_matrix = np.outer(modules, modules)
mod_heatmap = np.zeros(mod_matrix.shape)
for i in np.unique(modules):
    mod_inds = np.where(mod_matrix == i**2)
    mod_heatmap[mod_inds] = i -1


modules = modules == modules.T

# %%
df = pd.read_csv("../Data/Epat08_bip_labels.csv")
labels = df['labels']

sns.heatmap(mod_heatmap,cmap="crest", xticklabels=labels, yticklabels=labels)

# %%
