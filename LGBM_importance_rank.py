# 利用LGBM五折交叉验证计算重要性排序
import pandas as pd
import os
from collections import Counter
from lightgbm import LGBMClassifier
import numpy as np
import shap
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='shap')
warnings.filterwarnings('ignore', category=FutureWarning)

def normal_imp(mydict):
    mysum = sum(mydict.values())
    mykeys = mydict.keys()
    for key in mykeys:
        mydict[key] = mydict[key] / mysum
    return mydict

type_name = 'prevalent'
dpath = ''
disease_df = pd.read_csv(dpath + f'{type_name}_disease.csv')
disease_names = disease_df['NAME'].tolist()
target_dir = ""
disease_files = [os.path.join(target_dir, f"{disease}.csv") for disease in disease_names]

data = pd.read_csv(dpath + 'instance0_after_QC_preprocess.csv')
overlap_df = pd.read_csv(dpath + 'overlap_eid.csv')
data = data[data['eid'].isin(overlap_df['eid'])]
print(len(data))
feature_df = pd.read_excel(dpath + f'feature_dictionary.xlsx')
feature_lst = feature_df['feature'].tolist()
print(len(feature_lst))

for file in disease_files:
    disease_name = os.path.basename(file).replace('.csv', '')
    disease_df = pd.read_csv(file,usecols = ['eid','target_y','BL2Target_yrs'])
    if type_name == 'incident':
        disease_df = disease_df[(disease_df['target_y'] == 0) | ((disease_df['target_y'] == 1) & (disease_df['BL2Target_yrs'] > 0))]
    else:
        disease_df = disease_df[(disease_df['target_y'] == 0) | ((disease_df['target_y'] == 1) & (disease_df['BL2Target_yrs'] < 0))]
    merged_df = pd.merge(disease_df,data,on=['eid'],how = 'inner')
    rng = np.random.default_rng(42)
    merged_df['fold_id'] = rng.integers(0, 5, size=len(merged_df))

    all_X_test = []
    total_importance_counter = Counter()

    os.makedirs(dpath + 'importance/each_fold', exist_ok=True)
    file_path = os.path.join(dpath, f'importance/each_fold/{disease_name}.xlsx')
    with pd.ExcelWriter(file_path) as writer:
        for fold_id in range(5):
            train_df = merged_df[merged_df['fold_id'] != fold_id].copy()
            test_df = merged_df[merged_df['fold_id'] == fold_id].copy()
            X_train, y_train = train_df[feature_lst], train_df['target_y']
            X_test, y_test = test_df[feature_lst], test_df['target_y']

            my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, verbosity=1, seed=2023)
            my_lgb.fit(X_train, y_train)

            gain_imp = dict(zip(my_lgb.booster_.feature_name(),
                                my_lgb.booster_.feature_importance(importance_type='gain').tolist()))
            fold_imp_norm = normal_imp(gain_imp)
            total_importance_counter += Counter(fold_imp_norm)
            explainer = shap.TreeExplainer(my_lgb)
            shap_v = explainer.shap_values(X_test)

    # 5. 计算 5 折的平均重要性并保存
    avg_imp = {k: v / 5 for k, v in total_importance_counter.items()}
    avg_imp_df = pd.DataFrame({
        f'feature': list(avg_imp.keys()),
        'Average_TotalGain': list(avg_imp.values())
    })
    avg_imp_df.sort_values(by='Average_TotalGain', ascending=False, inplace=True)
    os.makedirs(dpath + 'importance/average', exist_ok=True)
    avg_imp_df.to_csv(dpath + f"/importance/average/{disease_name}.csv", index=False)

    print(f'====================================={disease_name}_Finished=======================================')