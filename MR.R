library(TwoSampleMR)
library(dplyr)
library(stringr)
library(vroom)
library(data.table)
library(ieugwasr)

# ==================== 1. 路径设置 ====================

# ==================== 2. 主循环分析 ====================
pairs <- fread(sign_list_path) 
unique_outcomes <- unique(pairs$NAME)

for (current_out in unique_outcomes) {
  message("\n>>> 正在处理 Outcome: ", current_out)
  
  log_list <- list()
  res_list <- list()
  
  # A. 加载 Outcome 与 映射表 (保持不变)
  out_file <- paste0(out_gwas_dir, current_out, ".txt")
  map_file <- paste0(anno_map_dir, current_out, "_rsid_map.csv")
  
  if(!file.exists(out_file) | !file.exists(map_file)) {
    warning("文件或映射表缺失，跳过 Outcome: ", current_out)
    next
  }
  
  out_raw <- vroom(out_file, show_col_types = FALSE)
  mapping <- fread(map_file)
  out_with_rsid <- out_raw %>% 
    left_join(mapping, by = "id") %>% 
    filter(!is.na(rsid))
  
  outcome_dat <- format_data(
    dat = out_with_rsid, 
    type = "outcome", 
    snp_col = "rsid",
    beta_col = "beta", 
    se_col = "se", 
    effect_allele_col = "EA",
    other_allele_col = "NEA", 
    pval_col = "p",
    eaf_col = "maf"  # 明确指定 eaf 列
  )
  
  # B. 循环处理该 Outcome 对应的 Exposures
  current_exps <- pairs$feature[pairs$NAME == current_out]
  
  for (current_exp in current_exps) {
    message("  -> Exposure: ", current_exp)
    
    log_entry <- data.table(
      outcome = current_out,
      exposure = current_exp,
      snps_post_clump = 0, # 这里改为 .clumps 文件里的数量
      snps_harmonized = 0,
      status = "Pending",
      note = ""
    )
    
    # 1. 设置文件路径
    exp_file <- paste0(exp_gwas_dir, current_exp, ".glm.linear")
    clump_file <- paste0(exp_clump_dir, current_exp, ".clumps") # 你的已clump文件
    
    if(!file.exists(exp_file) | !file.exists(clump_file)) {
      log_entry$status <- "Clump_or_GWAS_File_Missing"
      log_list[[current_exp]] <- log_entry
      next
    }
    
    # 2. 读取已 Clump 的结果 (获取 SNP 列表)
    # 注意：请根据你 .clumps 文件的实际列名修改（通常是 SNP 或 ID）
    clumped_snps <- fread(clump_file)
    if(nrow(clumped_snps) == 0) {
      log_entry$status <- "No_SNPs_In_Clump_File"
      log_list[[current_exp]] <- log_entry
      next
    }
    log_entry$snps_post_clump <- nrow(clumped_snps)
    
    # 3. 读取原始 GWAS 并提取这些 SNP 的统计量
    exp_raw <- vroom(exp_file, show_col_types = FALSE) %>%
      # 先改名现有的列
      rename(SNP=ID, beta=BETA, se=SE, P=P, EA=A1, EAF=A1_FREQ, N=OBS_CT) %>%
      # 动态生成 NEA 列：如果 EA 是 REF，那么 NEA 就是 ALT；反之亦然
      mutate(NEA = ifelse(EA == REF, ALT, REF)) %>%
      # 此时可以选出你需要的列，或者直接进入 format_data
      select(SNP, beta, se, P, EA, NEA, EAF, N)
    
    # 核心步骤：只保留在 clump 文件中出现的 SNP
    exp_filtered <- exp_raw %>% filter(SNP %in% clumped_snps$ID) 
    
    # 4. 格式化并 Harmonise
    exp_dat <- format_data(
      dat = exp_filtered,
      type = "exposure", snp_col = "SNP", beta_col = "beta",
      se_col = "se", effect_allele_col = "EA", other_allele_col = "NEA",
      pval_col = "P", eaf_col = "EAF", samplesize_col = "N",phenotype_col = "trait"
    )
    
    dat <- harmonise_data(exp_dat, outcome_dat)
    dat <- dat[dat$mr_keep == TRUE, ] 
    if(nrow(dat) > 0){
      dat$id.exposure <- current_exp  # 强制设为 "27800_avg"
      dat$exposure    <- current_exp
      dat$id.outcome  <- current_out  # 同时也建议把结局 ID 统一
      dat$outcome     <- current_out
    }
    
    log_entry$snps_harmonized <- nrow(dat)
    
    # 5. 运行 MR
    if(nrow(dat) < 1) {
      log_entry$status <- "No_SNPs_After_Harmonise"
      log_list[[current_exp]] <- log_entry
      next
    }
    
    mr_res <- tryCatch({ mr(dat) }, error = function(e) return(NULL))
    
    if(!is.null(mr_res) && nrow(mr_res) > 0) {
      log_entry$status <- "Success"
      res_list[[current_exp]] <- mr_res
    } else {
      log_entry$status <- "MR_Failed"
    }
    
    log_list[[current_exp]] <- log_entry
  }
  
  # C. 保存结果与日志 (保持不变)
  if(length(res_list) > 0) {
    fwrite(bind_rows(res_list), paste0(result_dir, current_out, "_mr_results.csv"))
  }
  
  if(length(log_list) > 0) {
    fwrite(bind_rows(log_list), paste0(result_dir_log, current_out, "_summary_log.csv"))
  }
  
  message("✅ Outcome: ", current_out, " 处理完毕。")
}
