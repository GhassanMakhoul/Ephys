# %%
import numpy as np
import scipy
import mat73
import os
import bct
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score as NMI
from bct import partition_distance as bct_nmi
# %%
DATA_DIR = '/home/ghassan/Documents/Ephys/Data/'
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
N_GAM = 1000
gammas = np.linspace(0.01,2,N_GAM)
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

    #let's get consensus of all runs
    mut_infos = np.zeros(int(n**2/2 -n/2))
    ind = 0
    for i in range(n):
        c_i = C[i,:]
        for j in range(i+1,n):
            c_j = C[j,:]
            [_, norm_mi] = bct_nmi(c_i, c_j)
            mut_infos[ind] =norm_mi
            ind +=1

    for i in range(n):
        ci = np.expand_dims(C[i,:],1)
        KKi = ci == ci.T
        W1 += KKi.astype(float)
    W1 = W1/n
    c_consennsus, q_consensus = bct.community_louvain(W1)
    return c_consennsus, q_consensus, mut_infos
NMI = []
N_ITER= 10 
for i,g in enumerate(gammas):
    try:
        
        ci, q, norm_mut_infos = consensus_modularity(bandSlice, ci, g,n=N_ITER)
        C[i,:] = ci
        Q[i,:] = q
        df = pd.DataFrame(data=norm_mut_infos.flatten())
        df['gamma'] = g
        NMI.append(df)
    except BCTParamError:
        param_errors[i] = g
        continue

num_errors = np.sum(param_errors > 0)
print(f"number of param errors: {num_errors}")
print(np.sum(~np.isnan(Q)))
# %%
#TODO refactor NMI calcs
N_GAM =1000
partition_sim = np.zeros((N_GAM, N_GAM))
valid_run_inds = np.where(param_errors==0)[0]
pdist =0
for i in range(len(gammas)):
    c_i = C[i,:]
    for j in range(i+1, len(gammas)):
        c_j = C[j,:]
        [_, pdist] = bct_nmi(c_i, c_j)
        partition_sim[i,j] = pdist

    

data = np.zeros((max(valid_run_inds.shape),3))
data[:,0] =  [len(set(C[i,:])) for i in valid_run_inds]
data[:,1] = [gammas[i] for i in valid_run_inds]
# data[:,2] = partition_sim.ravel()
modularity_df = pd.DataFrame(data,\
                              columns=['N_modules', 'gamma_vals','NMI'])

# 
# %%
sns.lineplot(modularity_df, x='gamma_vals', y='N_modules')
plt.show()
# sns.lineplot(modularity_df, x='gamma_vals', y='NMI', ax=ax2)
gticks = gammas[np.linspace(0,N_GAM-1,25 ,dtype=int)]
# gticks = [str(g) for g in gticks]
df = pd.DataFrame(data=partition_sim, columns = gammas, index=gammas)
ax2 = sns.heatmap(df)

plt.show()

# %%

modules = np.reshape(C[valid_run_inds[700], : ], (42,1))+1
mod_matrix = np.outer(modules, modules)
mod_heatmap = np.zeros(mod_matrix.shape)
for i in np.unique(modules):
    mod_inds = np.where(mod_matrix == i**2)
    mod_heatmap[mod_inds] = i -1


modules = modules == modules.T


# %%
# looks like the best gammas are around 1.12
best_gamma = 1.12
best_g_ind = np.where(np.logical_and(gammas < 1.12 , gammas >1.11))[0]
print(best_g_ind)
modules = np.reshape(C[best_g_ind[0], : ], (42,1))+1
mod_matrix = np.outer(modules, modules)
mod_heatmap = np.zeros(mod_matrix.shape)
for i in np.unique(modules):
    mod_inds = np.where(mod_matrix == i**2)
    mod_heatmap[mod_inds] = i -1


# %%
df = pd.read_csv("../Data/Epat08_bip_labels.csv")
labels = df['labels']

sns.heatmap(mod_heatmap,cmap="crest", xticklabels=labels, yticklabels=labels)


# %%
nmi_df = pd.concat(NMI, axis=0)
nmi_df.columns = ['NMI', 'gamma']

#TODO fix df
sns.lineplot(data=nmi_df.fillna(-1), y='NMI', x='gamma')
# %%
