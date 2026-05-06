# 和新发疾病的关联分析——COX
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
import os
from joblib import Parallel, delayed
from tqdm import tqdm

def process(tmp):
    import warnings
    warnings.filterwarnings('error')
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    tmp_x_df = OCT_df[['eid',tmp]]
    tmp_df = pd.merge(tg_df, tmp_x_df, how='inner', on=['eid'])
    tmp_df.rename(columns={tmp: 'x'}, inplace=True)
    rm_eid_idx = tmp_df.index[tmp_df.x.isnull() == True]
    tmp_df.drop(rm_eid_idx, axis=0, inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    nb_all, nb_case = len(tmp_df), tmp_df.target_y.sum()
    prop_case = np.round(nb_case / nb_all * 100, 3)

    # 按顺序尝试不同的惩罚值：从无惩罚开始，逐步增加
    penalties = [0, 0.001, 0.01, 0.1]
    for penalizer in penalties:
        try:
            cph = CoxPHFitter(penalizer=penalizer)
            cph.fit(tmp_df, duration_col='BL2Target_yrs', event_col='target_y', formula=my_formula)
            hr = np.round(cph.hazard_ratios_.x, 5)
            ci = np.exp(cph.confidence_intervals_)
            lbd = np.round(ci.loc['x'].iloc[0], 5)
            ubd = np.round(ci.loc['x'].iloc[1], 5)
            pval = cph.summary.p.x
            tmpout = [tmp, nb_all, nb_case, prop_case, hr, lbd, ubd, pval, f'success_{penalizer}']
            return tmpout
        except:
            continue
    tmpout = [tmp, nb_all, nb_case, prop_case, np.nan, np.nan, np.nan, np.nan, 'fail']
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

cov_df = pd.read_csv('')
disease_sex = pd.read_csv('')
disease_sex= disease_sex['NAME'].tolist()

all_results = []
for file in tqdm(disease_files, desc="Processing diseases"):
    disease_name = os.path.basename(file).replace('.csv', '')
    print(f'正在处理{disease_name}')
    disease_df = pd.read_csv(file,usecols = ['eid','target_y','BL2Target_yrs'])
    disease_df = disease_df[ (disease_df['target_y'] == 0) | ((disease_df['target_y'] == 1) & (disease_df['BL2Target_yrs'] > 0))]
    tg_df = pd.merge(disease_df, cov_df, how='inner', on=['eid'])
    if disease_name in disease_sex:
        my_formula = "Age + TDI + BMI + smoke + drink + avMSE + x + Race"
    else:
        my_formula = "Age + Sex + TDI + BMI + smoke + drink + x + avMSE +  Race"
    tgt_out_df = Parallel(n_jobs=5)(delayed(process)(tmp) for tmp in feature_lst)
    tgt_out_df = pd.DataFrame(tgt_out_df)
    tgt_out_df.columns = ['feature', 'nb_individuals', 'nb_case', 'prop_case(%)', 'hr', 'hr_lbd', 'hr_ubd', 'pval_raw','status']
    tgt_out_df['NAME'] = disease_name
    all_results.append(tgt_out_df)

final_df = pd.concat(all_results, ignore_index=True)
columns_order = ['NAME'] + [col for col in final_df.columns if col != 'NAME']
final_df = final_df[columns_order]
final_df.to_csv(f'', index=False)