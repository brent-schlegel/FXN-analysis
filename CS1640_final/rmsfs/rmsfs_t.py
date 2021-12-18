import pandas as pd
import numpy  as np
import scipy.stats as stats
import math

wt_data = pd.read_csv('all_wt_rmsfs.csv')
m_data  = pd.read_csv('rmsfs.csv')
rmsfs_scored = m_data.filter('resnum', axis=1)
rmsfs_scored["wt_mean"]     = ""
rmsfs_scored["mutant_mean"] = ""
rmsfs_scored["p_score"]     = ""


m_ar  = m_data.iloc[: , 1:].values
wt_ar = wt_data.iloc[: , 1:].values

for index, row in rmsfs_scored.iterrows():
        t_test = stats.ttest_ind(a=wt_ar[index], b=m_ar[index], equal_var=True)
	pval = t_test.pvalue
	if math.isnan(pval):
		pval = 1
        rmsfs_scored.iloc[index] = {'wt_mean':np.mean(wt_ar[index]),'mutant_mean':np.mean(m_ar[index]),'p_score':pval}

rmsfs_scored.to_csv(r'rmsfs_scored.csv')
