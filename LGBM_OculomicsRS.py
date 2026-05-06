# 在测试集上用LGBM计算Oculomics_RS
import pandas as pd
import os
from lightgbm import LGBMClassifier
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV
from sklearn.metrics import roc_auc_score, average_precision_score
import numpy as np

dpath = ''
disease_file = ""
disease_df = pd.read_csv(disease_file,usecols = ['eid','target_y','BL2Target_yrs'])
data = pd.read_csv(dpath + f'')
data = pd.merge(disease_df,data,how='inner', on=['eid'])

split_df = pd.read_csv(dpath + '/data/eid_split.csv')
overlap_df = pd.read_csv(dpath + f'/data/overlap_eid.csv')
data = data[data['eid'].isin(overlap_df['eid'])]
train_eid = split_df.loc[split_df['group']==0,'eid']
train_df = data[data['eid'].isin(train_eid)]
test_eid = split_df.loc[split_df['group']==1,'eid']
test_df = data[data['eid'].isin(test_eid)]

feature_df = pd.read_excel(dpath + f'/data/dictionary.xlsx',sheet_name='Sheet1')
feature_lst = feature_df['feature'].tolist()
X_train, X_test = train_df[feature_lst], test_df[feature_lst]
y_train, y_test = train_df['target_y'], test_df['target_y']


base_lgb = LGBMClassifier(
    objective='binary',
    metric='auc',
    is_unbalance=True,
    verbosity=-1,
    seed=2023
)

param_dist = {
    'n_estimators': [100, 200, 300, 500],
    'learning_rate': [0.01, 0.03, 0.05, 0.1],
    'num_leaves': [15, 31, 63],
    'max_depth': [-1, 3, 5, 7],
    'min_child_samples': [20, 50, 100],
    'subsample': [0.7, 0.8, 1.0],
    'colsample_bytree': [0.7, 0.8, 1.0],
    'reg_alpha': [0, 0.1, 1],
    'reg_lambda': [0, 0.1, 1, 5]
}

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=2023)

search = RandomizedSearchCV(
    estimator=base_lgb,
    param_distributions=param_dist,
    n_iter=50,
    scoring='roc_auc',
    cv=cv,
    n_jobs=-1,
    random_state=2023,
    verbose=1
)

search.fit(X_train, y_train)

print("Best CV AUC:", search.best_score_)
print("Best params:", search.best_params_)

best_lgb = search.best_estimator_
feature_pred = best_lgb.predict_proba(X_test)[:, 1]

print("Test AUC:", roc_auc_score(y_test, feature_pred))
print("Test PR-AUC:", average_precision_score(y_test, feature_pred))

final_pred_df = pd.DataFrame({
    'eid': test_df.eid,
    'target_y': test_df.target_y,
    'BL2Target_yrs': test_df.BL2Target_yrs,
    'y_pred_feature': feature_pred
})
final_pred_df.to_csv('',index=False)