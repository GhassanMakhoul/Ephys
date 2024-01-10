# %%
import pandas as pd
import seaborn as sns
# %%
label_f = '~/Ephys/Data/all_pats_bipole_soz_labels.csv'
label_df = pd.read_csv(label_f, header=None)
label_df.columns = ['subj','bipole','label']
label_df.head(5)
# %%

def split_bipole(bip_df: pd.DataFrame):
    """splits the bipole column of a bipole df
    duplicates rows

    """
    assert 'bipole' in bip_df.columns, "Need bipole column!"

    contact1 = bip_df.bipole.apply(lambda x: x.split("-")[0])
    contact2 = bip_df.bipole.apply(lambda x: x.split("-")[0])

    df2 = bip_df.copy(deep=True)

    bip_df['contact'] = contact1
    df2['contact'] = contact2

    return pd.concat([bip_df,df2])

label_df = split_bipole(label_df)
label_df.head()

# %%
res_f  = '/mnt/ernie_main/Ghassan/ephys/test/Epat26_stim.csv'
res_df = pd.read_csv(res_f)
subj = res_df.subj.values[0]
label_df = label_df[label_df.subj == subj]

res_df.head(5)
# %%
label_df.head(5)
# %%
stim_res_df = res_df.merge(label_df[['label','contact']],\
              left_on='resp_reg', right_on='contact')


sns.boxplot(stim_res_df,x='label', y='alphas_prime')

#%%


