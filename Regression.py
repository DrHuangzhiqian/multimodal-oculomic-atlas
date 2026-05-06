# 蛋白质和眼部feature之间的关联
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection,multipletests
#from mne.stats import bonferroni_correction
from tqdm import tqdm
# ============================================
#   线性回归函数：feature ~ protein + covs
# ============================================
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel

def run_regression(df, feature, protein, covariates,
                   binary_features, ordinal_features):

    cols = [feature, protein] + covariates
    tmp = df[cols].dropna()

    y = tmp[feature]
    X = tmp[[protein] + covariates]
    X = sm.add_constant(X)

    res = {
        "OCT_feature": feature,
        "Protein": protein,
        "N": len(tmp)
    }

    try:
        # =========================
        # 1️⃣ 二分类 → Logistic
        # =========================
        if feature in binary_features:
            model = sm.Logit(y, X).fit(disp=0)
            res["Coef_protein"] = model.params[protein]
            res["P_protein"] = model.pvalues[protein]

        # =========================
        # 2️⃣ 有序多分类 → Ordinal Logistic
        # =========================
        elif feature in ordinal_features:
            model = OrderedModel(
                y,
                tmp[[protein] + covariates],
                distr='logit'
            )
            fit = model.fit(method='bfgs', disp=False)

            res["Coef_protein"] = fit.params[protein]
            res["P_protein"] = fit.pvalues[protein]

        # =========================
        # 3️⃣ 连续 → 线性回归（原来的）
        # =========================
        else:
            model = sm.OLS(y, X).fit()
            res["Coef_protein"] = model.params[protein]
            res["P_protein"] = model.pvalues[protein]
            res["R2"] = model.rsquared

    except Exception as e:
        # 防止某些模型报错（比如完全分离）
        res["Coef_protein"] = None
        res["P_protein"] = None

    return res


# === 数据读取 ===
pro_df = pd.read_csv('')
pro_lst = pro_df.columns[1:].tolist()

dpath = ''
OCT_cov = pd.read_csv(dpath + 'data/cov_0.csv')
egfr_df = pd.read_csv('')
egfr_df = egfr_df.dropna()
OCT_cov = pd.merge(OCT_cov, egfr_df, on='eid', how='inner')
covariates = OCT_cov.columns[1:].tolist()

# 合并蛋白质 + 协变量
merged_df = pd.merge(pro_df, OCT_cov, on='eid', how='inner')

# 读取 OCT 特征
OCT_df = pd.read_csv(dpath + '').rename(columns={'VIT':'VIT_OCT'})
feature_df = pd.read_csv('')
feature_lst = feature_df['feature'].tolist()
#feature_lst = ["Drusen","HIS","VMA","NM","AMD","RH","ERM"]
binary_features = ["Drusen","HIS","VMA","NM","AMD","RH"]
ordinal_features = ["ERM"]

# 合并 OCT 特征
merged_df = pd.merge(merged_df, OCT_df, on='eid', how='inner')
# ============================================
#      遍历所有 (OCT_feature × protein)
# ============================================
for feature in tqdm(feature_lst, desc="Processing diseases"):
    print(f'正在处理{feature}')
    results = []
    for pro in pro_lst:
        res = run_regression(
            merged_df,
            feature,
            pro,
            covariates,
            binary_features,
            ordinal_features
        )
        results.append(res)
    results_df = pd.DataFrame(results)
    _, p_f_bfi, _, _ = multipletests(
    results_df['P_protein'].fillna(1),
    alpha=0.05,
    method='bonferroni')
    results_df['pval_bfi'] = p_f_bfi
    results_df.loc[results_df['pval_bfi'] >= 1, 'pval_bfi'] = 1

    _, p_fdr = fdrcorrection(results_df['P_protein'].fillna(1),alpha=0.05)
    results_df['pval_fdr'] = p_fdr
    results_df.to_csv(f'', index=False)

print("完成！")