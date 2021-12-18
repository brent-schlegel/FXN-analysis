import pandas as pd
import numpy  as np
import scipy.stats as stats
import math

wt_data = pd.read_csv('all_wt_waters.csv')
m_data  = pd.read_csv('waters.csv')
waters_scored = m_data.filter('resnum', axis=1)
waters_scored["wt_mean"]     = ""
waters_scored["mutant_mean"] = ""
waters_scored["p_score"]     = ""


m_ar  = m_data.iloc[: , 1:].values
wt_ar = wt_data.iloc[: , 1:].values

for index, row in waters_scored.iterrows():
        t_test = stats.ttest_ind(a=wt_ar[index], b=m_ar[index], equal_var=True)
	pval = t_test.pvalue
	if math.isnan(pval):
		pval = 1
        waters_scored.iloc[index] = {'wt_mean':np.mean(wt_ar[index]),'mutant_mean':np.mean(m_ar[index]),'p_score':pval}

waters_scored.to_csv(r'waters_scored.csv')
