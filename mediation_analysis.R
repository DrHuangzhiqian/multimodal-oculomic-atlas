# 加载包
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(mediation)
  library(foreach)
  library(doParallel)
})
cat("--- 所有依赖包处理完成 ---\n")

# --- 1. 定义参数和路径 ---

# 分析参数
MIN_N <- 50      # 最小样本量
SIMS <- 1000     # Bootstrap次数
set.seed(2025)   # 保证可重复性

# 并行计算设置
num_cores <- max(1, min(3, parallel::detectCores(logical = TRUE) - 1))

if (exists("cl")) try(stopCluster(cl), silent = TRUE)
if (.Platform$OS.type == "windows") {
  cl <- tryCatch({
    parallel::makePSOCKcluster(
      num_cores,
      rscript = file.path(R.home("bin"), "Rscript")
    )
  }, error = function(e) {
    NULL
  })
  if (!is.null(cl)) {
    doParallel::registerDoParallel(cl)
    cat(sprintf("已注册 %d 个核心用于并行计算 (PSOCK)。进度将输出到: %s\n", num_cores, log_file))
  } else {
    doParallel::registerDoSEQ()
    cat("并行集群创建失败，切换为顺序执行。\n")
  }
} else {
  # Linux/Mac 优先使用 Fork 模式
  doParallel::registerDoParallel(cores = num_cores)
  cat(sprintf("已注册 %d 个核心用于并行计算 (Fork)。进度将输出到: %s\n", num_cores, log_file))
}

# --- 2. 读取数据 ---
cat("正在读取数据...\n")
data_x <- fread(path_x)
data_m <- fread(path_m)
data_y <- fread(path_y)
data_cov <- fread(path_cov)
pair_list <- fread(path_pairs)

cat(sprintf("读取完成:\n X: %d 行, %d 列\n M: %d 行, %d 列\n Y: %d 行, %d 列\n Cov: %d 行, %d 列\n Pairs: %d 对\n",
            nrow(data_x), ncol(data_x), nrow(data_m), ncol(data_m), 
            nrow(data_y), ncol(data_y), nrow(data_cov), ncol(data_cov), nrow(pair_list)))

# --- 3. 数据合并与预处理 ---
cat("正在合并数据...\n")

# 使用 Left Join 合并所有数据 (基于 eid)
# 以 X 数据集为主 (保留所有 X 样本)，依次左连接 Cov, M, Y
# 注意: mediation 包在分析时会自动剔除缺失值 (na.omit)，所以这里我们先保留所有样本
full_data <- merge(data_x, data_cov, by = "eid", all.x = TRUE)
full_data <- merge(full_data, data_m, by = "eid", all.x = TRUE)
full_data <- merge(full_data, data_y, by = "eid", all.x = TRUE)

cat(sprintf("合并后数据维度: %d 行, %d 列\n", nrow(full_data), ncol(full_data)))

# 提取 X 变量列表 (直接从 data_x 中获取，排除 eid)
x_vars <- setdiff(names(data_x), "eid")
cat(sprintf("识别到 %d 个 X 变量 (直接从 X 文件提取)。\n", length(x_vars)))

# 清理不再需要的单独数据集以释放内存
rm(data_x, data_m, data_y, data_cov)
gc()

# --- 4. 变量因子化 ---
cat("正在处理变量因子化...\n")

# 指定需要因子化的 X 变量 (列名可能是数字字符串)
x_vars_to_factor <- c("1180", "1418", "1428", "1448", "1468", "1508")
# 指定需要因子化的 Cov 变量
cov_vars_to_factor <- c("drink", "smoke", "Race", "Sex")

all_vars_to_factor <- unique(c(x_vars_to_factor, cov_vars_to_factor))

# 检查并转换
existing_factor_vars <- intersect(all_vars_to_factor, names(full_data))
if(length(existing_factor_vars) > 0){
  cat("  将以下变量转换为因子: ", paste(existing_factor_vars, collapse = ", "), "\n")
  full_data[, (existing_factor_vars) := lapply(.SD, as.factor), .SDcols = existing_factor_vars]
} else {
  cat("  未找到指定的因子变量，跳过因子化步骤。\n")
}

# --- 5. 构建任务列表 ---
cat("正在构建分析任务列表...\n")

# 构建任务清单 (Grid Expand)
# 每个 M-Y Pair 都要跑所有的 X
tasks <- list()
counter <- 1

for (i in 1:nrow(pair_list)) {
  m_var <- pair_list$M[i]
  y_var <- pair_list$Y[i]
  
  for (x_var in x_vars) {
    tasks[[counter]] <- list(X = x_var, M = m_var, Y = y_var)
    counter <- counter + 1
  }
}

# 转为 data.frame 方便 foreach 使用
task_df <- rbindlist(tasks)
cat(sprintf("总计生成 %d 个中介分析任务。\n", nrow(task_df)))

# --- 6. 并行执行中介分析 ---
cat("--- 开始并行计算 ---\n")

# 显式清理内存
gc()

# 定义协变量公式部分
# 为了保险，我们假设 cov.csv 的 header 我们已经看过: drink, smoke, Race, Age, TDI, BMI, Sex, avMSE
base_covars <- c("drink", "smoke", "Race", "Age", "TDI", "BMI", "Sex", "avMSE")

results <- foreach(i = 1:nrow(task_df),
                   .combine = "rbind",
                   .packages = c("mediation", "dplyr", "data.table"),
                   .errorhandling = "remove",
                   .inorder = FALSE) %dopar% {
                     
                     # 日志记录 (每 100 个任务记录一次)
                     if (i %% 100 == 0) {
                       try({
                         sink(log_file, append = TRUE)
                         cat(sprintf("[%s] 正在处理第 %d / %d 个任务...\n", format(Sys.time(), "%H:%M:%S"), i, nrow(task_df)))
                         sink()
                       }, silent = TRUE)
                     }
                     
                     # 获取任务变量
                     X_var <- task_df$X[i]
                     M_var <- task_df$M[i]
                     Y_var <- task_df$Y[i]
                     
                     # 检查变量是否存在
                     if (!all(c(X_var, M_var, Y_var, base_covars) %in% names(full_data))) {
                       return(NULL)
                     }
                     
                     # 准备子数据集 (转为 data.frame 避免指针问题，保持与原脚本一致)
                     analysis_vars <- c(X_var, M_var, Y_var, base_covars)
                     dt_subset <- as.data.frame(full_data)[, analysis_vars]
                     dt_subset <- na.omit(dt_subset)
                     
                     # 检查样本量
                     if (nrow(dt_subset) < MIN_N) return(NULL)
                     
                     # 检查 Y 变量单一性 (逻辑回归要求至少有两个水平)
                     if (length(unique(dt_subset[[Y_var]])) < 2) return(NULL)
                     
                     # 动态筛选协变量: 移除常数协变量
                     covars_subset <- dt_subset[, base_covars, drop = FALSE]
                     covars_lenuniq <- sapply(covars_subset, function(x) length(unique(x)))
                     valid_covars <- names(covars_lenuniq[covars_lenuniq > 1])
                     
                     # 构建公式
                     covars_str <- ""
                     if (length(valid_covars) > 0) {
                       covars_str <- paste(" +", paste(valid_covars, collapse = " + "))
                     }
                     
                     tryCatch({
                       # 1. Mediator Model (M ~ X + Covs)
                       # 线性回归 (假设 M 是连续变量)
                       f_med <- as.formula(paste0("`", M_var, "` ~ `", X_var, "`", covars_str))
                       med.fit <- lm(f_med, data = dt_subset)
                       
                       # 2. Outcome Model (Y ~ X + M + Covs)
                       # 逻辑回归 (假设 Y 是二分类变量 0/1)
                       f_out <- as.formula(paste0("`", Y_var, "` ~ `", X_var, "` + `", M_var, "`", covars_str))
                       out.fit <- glm(f_out, data = dt_subset, family = binomial())
                       
                       # 检查收敛性
                       if (!out.fit$converged) return(NULL)
                       
                       # 3. Mediation Analysis
                       # treat: X, mediator: M
                       mo <- mediate(med.fit, out.fit, treat = X_var, mediator = M_var, robustSE = TRUE, sims = SIMS)
                       
                       # 整理结果
                       res_row <- data.frame(
                         X = X_var,
                         M = M_var,
                         Y = Y_var,
                         n_obs = nrow(dt_subset),
                         
                         # ACME (间接效应)
                         acme_est = mo$d.avg,
                         acme_p = mo$d.avg.p,
                         acme_ci_low = mo$d.avg.ci[1],
                         acme_ci_high = mo$d.avg.ci[2],
                         
                         # ADE (直接效应)
                         ade_est = mo$z.avg,
                         ade_p = mo$z.avg.p,
                         ade_ci_low = mo$z.avg.ci[1],
                         ade_ci_high = mo$z.avg.ci[2],
                         
                         # Total Effect (总效应)
                         total_est = mo$tau.coef,
                         total_p = mo$tau.p,
                         total_ci_low = mo$tau.ci[1],
                         total_ci_high = mo$tau.ci[2],
                         
                         # Prop. Mediated (中介占比)
                         prop_med_est = mo$n.avg,
                         prop_med_ci_low = mo$n.avg.ci[1],
                         prop_med_ci_high = mo$n.avg.ci[2]
                       )
                       
                       return(res_row)
                       
                     }, error = function(e) {
                       return(NULL)
                     })
                   }

# --- 7. 保存结果 ---
if (!is.null(results) && nrow(results) > 0) {
  fwrite(results, output_file)
  
  cat(sprintf("\n分析全部完成！\n成功计算路径数: %d / %d\n结果已保存至: %s\n", 
              nrow(results), nrow(task_df), output_file))
  
  # 简单预览显著结果
  sig_count <- sum(results$acme_p < 0.05)
  cat(sprintf("其中 ACME 显著 (P < 0.05) 的路径数: %d\n", sig_count))
  
} else {
  cat("\n分析完成，但未生成任何有效结果。\n")
}

# --- 8. 清理 ---
if (exists("cl")) try(stopCluster(cl), silent = TRUE)
