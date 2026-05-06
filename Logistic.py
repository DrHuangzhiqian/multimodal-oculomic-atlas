# 和既往疾病的关联分析——logistic
import pandas as pd
import numpy as np
import os
import statsmodels.api as sm
from joblib import Parallel, delayed
from tqdm import tqdm

def process(tmp, tg_df, cov_lst):
    import warnings
    warnings.filterwarnings('error')

    tmp_x_df = OCT_df[['eid', tmp]]
    tmp_df = pd.merge(tg_df, tmp_x_df, how='inner', on=['eid'])
    tmp_df.rename(columns={tmp: 'x'}, inplace=True)
    rm_eid_idx = tmp_df.index[tmp_df.x.isnull() == True]
    tmp_df.drop(rm_eid_idx, axis=0, inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    nb_all, nb_case = len(tmp_df), tmp_df.target_y.sum()
    prop_case = np.round(nb_case / nb_all * 100, 3)
    Y = tmp_df.target_y
    X = tmp_df[cov_lst + ['x']]
    try:
        try:
            log_mod = sm.Logit(Y, sm.add_constant(X)).fit(disp=False,maxiter=1000)
            oratio = np.round(np.exp(log_mod.params).loc['x'], 5)
            pval = log_mod.pvalues.loc['x']
            ci_mod = log_mod.conf_int(alpha=0.05)
            lbd, ubd = np.round(np.exp(ci_mod.loc['x'][0]), 5), np.round(np.exp(ci_mod.loc['x'][1]), 5)
            tmpout = [tmp, nb_all, nb_case, prop_case, oratio, lbd, ubd, pval,'success']
        except :
            log_mod = sm.Logit(Y, sm.add_constant(X)).fit(method='lbfgs',disp=False,maxiter=1000)
            oratio = np.round(np.exp(log_mod.params).loc['x'], 5)
            pval = log_mod.pvalues.loc['x']
            ci_mod = log_mod.conf_int(alpha=0.05)
            lbd, ubd = np.round(np.exp(ci_mod.loc['x'][0]), 5), np.round(np.exp(ci_mod.loc['x'][1]), 5)
            tmpout = [tmp, nb_all, nb_case, prop_case, oratio, lbd, ubd, pval,'success_lbfgs']
    except:
        tmpout = [tmp, nb_all, nb_case, prop_case, np.nan, np.nan, np.nan, np.nan,'fail']
    return tmpout

disease_df = pd.read_excel('')
disease_names = disease_df['NAME'].tolist()
target_dir = ""
disease_files = [os.path.join(target_dir, f"{disease}.csv") for disease in disease_names]

OCT_df = pd.read_csv('')
eid_df = pd.read_csv('')
eid_discovery_eid = eid_df.loc[eid_df['group']==0,'eid']
OCT_df = OCT_df[OCT_df['eid'].isin(eid_discovery_eid)]
feature_lst = OCT_df.columns[1:].tolist()
#print(len(OCT_df))

cov_df = pd.read_csv('')
disease_sex = pd.read_csv('')
disease_sex= disease_sex['NAME'].tolist()

all_results = []
for file in tqdm(disease_files, desc="Processing diseases"):
    disease_name = os.path.basename(file).replace('.csv', '')
    print(f'正在处理{disease_name}')
    disease_df = pd.read_csv(file,usecols = ['eid','target_y','BL2Target_yrs'])
    disease_df = disease_df[ (disease_df['target_y'] == 0) | ((disease_df['target_y'] == 1) & (disease_df['BL2Target_yrs'] < 0))]
    tg_df = pd.merge(disease_df, cov_df, how='inner', on=['eid'])
    if disease_name in disease_sex:
        cov_lst = ['Age','TDI','BMI','smoke','drink','avMSE','Race']
    else:
        cov_lst = ['Age','Sex','TDI','BMI','smoke','drink','avMSE','Race']
    tgt_out_df = Parallel(n_jobs=5)(delayed(process)(tmp,tg_df,cov_lst) for tmp in feature_lst)
    tgt_out_df = pd.DataFrame(tgt_out_df)
    tgt_out_df.columns = ['feature', 'nb_individuals', 'nb_case', 'prop_case(%)', 'oratio', 'or_lbd', 'or_ubd', 'pval_raw','status']
    tgt_out_df['NAME'] = disease_name
    all_results.append(tgt_out_df)

final_df = pd.concat(all_results, ignore_index=True)
columns_order = ['NAME'] + [col for col in final_df.columns if col != 'NAME']
final_df = final_df[columns_order]
final_df.to_csv(f'', index=False)