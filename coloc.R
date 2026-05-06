library(data.table)
library(dplyr)
library(coloc)
library(vroom)

# ================== 路径 ==================

log_file <- file.path(result_dir, "coloc_log_new.csv")
log_df <- data.frame(Disease=character(), Feature=character(), Status=character(), Message=character(), stringsAsFactors=FALSE)

# ================== 参数 ==================
p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5

# 读取关联对
pairs <- fread(sign_list_path) 
unique_outcomes <- unique(pairs$NAME)


total_outcomes <- length(unique_outcomes)
outcome_count <- 0
for (current_out in unique_outcomes) {
  outcome_count <- outcome_count + 1
  start_time <- Sys.time()
  message(sprintf("\n[%s] >>> 进度: %d/%d | 正在处理 Outcome: %s", 
                  format(start_time, "%H:%M:%S"), outcome_count, total_outcomes, current_out))
  
  out_file <- paste0(out_gwas_dir, current_out, ".txt")
  map_file <- paste0(anno_map_dir, current_out, "_rsid_map.csv")
  
  if(!file.exists(out_file) | !file.exists(map_file)) {
    log_df <- rbind(log_df, data.frame(Disease=current_out, Feature=NA, Status="Failed", Message="Outcome GWAS 或映射文件缺失"))
    next
  }
  
  # 读取 outcome 数据并加上 rsid
  out_raw <- fread(out_file)
  mapping <- fread(map_file)
  out_gwas <- out_raw %>% 
    left_join(mapping, by="id")
  out_gwas <- out_gwas[rsid != "" & !is.na(rsid)]
  out_gwas <- out_gwas %>% 
    rename(SNP=rsid, chr=chr, POS=pos, BETA=beta, SE=se, MAF=maf, P=p)
  
  current_exps <- pairs$feature[pairs$NAME == current_out]
  exp_count <- 0
  total_exps <- length(current_exps)
  
  for(current_exp in current_exps) {
    exp_count <- exp_count + 1
    # 实时显示 Exposure 进度
    message(sprintf("  -> [%d/%d] Exposure: %s", exp_count, total_exps, current_exp))
    exp_file <- paste0(exp_gwas_dir, current_exp, ".glm.linear")
    clump_file <- paste0(exp_clump_dir, current_exp, ".clumps")
    
    if(!file.exists(exp_file)) {
      log_df <- rbind(log_df, data.frame(Disease=current_out, Feature=current_exp, Status="Failed", Message="Exposure GWAS文件缺失"))
      next
    }
    
    if(!file.exists(clump_file)) {
      log_df <- rbind(log_df, data.frame(Disease=current_out, Feature=current_exp, Status="Failed", Message="Exposure clump 文件缺失"))
      next
    }
    
    exp_gwas <- fread(exp_file)
    colnames(exp_gwas)[1] <- "chr"  # 确保列名一致
    
    clumped_snps <- fread(clump_file)
    snv_list <- clumped_snps$ID  # clump文件的SNP ID列
    
    coloc_results_final_table <- data.frame()
    for(snv in snv_list) {
      lead_snv_data <- exp_gwas[exp_gwas$ID == snv,]
      if(nrow(lead_snv_data) == 0) {
        log_df <- rbind(log_df, data.frame(Disease=current_out, Feature=current_exp, Status="Failed", Message=paste0("Lead SNP ", snv, " not found in exposure GWAS")))
        next
      }
      
      lead_chr <- lead_snv_data$chr[1]
      lead_pos <- lead_snv_data$POS[1]
      
      # 提取 ±500kb 区域
      exp_sub <- exp_gwas %>% filter(chr == lead_chr & POS >= (lead_pos - 500000) & POS <= (lead_pos + 500000)) %>%
        distinct(ID, .keep_all = TRUE)
      
      
      out_sub <- out_gwas %>% filter(chr == lead_chr & POS >= (lead_pos - 500000) & POS <= (lead_pos + 500000)) %>%
        distinct(SNP, .keep_all = TRUE)
      
      common_snps <- intersect(exp_sub$ID, out_sub$SNP)
      if(length(common_snps) == 0) {
        log_df <- rbind(
          log_df,
          data.frame(
            Disease = current_out,
            Feature = current_exp,
            Status = "Failed",
            Message = paste0("No common SNPs for lead SNP ", snv)
          )
        )
        next
      }
      
      exp_sub_flitered <- exp_sub %>%
        filter(ID %in% common_snps) %>%
        distinct(ID, .keep_all = TRUE) %>%
        arrange(match(ID, common_snps))
      
      out_sub_flitered <- out_sub %>%
        filter(SNP %in% common_snps) %>%
        distinct(SNP, .keep_all = TRUE) %>%
        arrange(match(SNP, common_snps))
      
      dataset1 <- list(
          snp = exp_sub_flitered$ID,
          beta = exp_sub_flitered$BETA,
          varbeta = exp_sub_flitered$SE^2,
          #position = exp_sub_flitered$POS,
          type = "quant",
          sdY = 1  # 因为使用 --variance-standardize
        )
      
      dataset2 <- list(
        snp = out_sub_flitered$SNP,
        beta = out_sub_flitered$BETA,
        varbeta = out_sub_flitered$SE^2,
        #position = out_sub_flitered$POS,
        type = "cc"
      )
      
      # 运行 coloc
      coloc_res <- NULL
      tryCatch({
        coloc_res <- coloc.abf(dataset1=dataset1, dataset2=dataset2, p1=p1, p2=p2, p12=p12)
        
        # 【全量提取结果】
        if(!is.null(coloc_res$summary)) {
          summary_stats <- as.data.frame(t(coloc_res$summary))
          results_df <- coloc_res$results
          
          # 找到该区域内最强的 SNP
          top_snp_row <- results_df[which.max(results_df$SNP.PP.H4), ]
          
          res_entry <- data.frame(
            Lead_SNV      = snv,
            NSNPS         = summary_stats$nsnps,
            PP_H0         = summary_stats$PP.H0.abf,
            PP_H1         = summary_stats$PP.H1.abf,
            PP_H2         = summary_stats$PP.H2.abf,
            PP_H3         = summary_stats$PP.H3.abf,
            PP_H4         = summary_stats$PP.H4.abf,
            Best_SNP      = top_snp_row$snp,
            Best_SNP_PPH4 = top_snp_row$SNP.PP.H4
          )
          coloc_results_final_table <- rbind(coloc_results_final_table, res_entry)
        }
      }, error=function(e) {
        message("    !! SNP ", snv, " Error: ", e$message)
        log_df <<- rbind(log_df, data.frame(Disease=current_out, Feature=current_exp, 
                                            Status="Failed", Message=paste0(snv, ": ", e$message)))
      })
    } # 结束 snv 循环
    
    if(nrow(coloc_results_final_table) > 0) {
      coloc_results_final_table$Exposure <- current_exp
      coloc_results_final_table$Outcome  <- current_out
      out_path <- file.path(result_dir, paste0(current_out, "_", current_exp, "_coloc.csv"))
      fwrite(coloc_results_final_table, out_path)
      log_df <- rbind(log_df, data.frame(Disease=current_out, Feature=current_exp, Status="Success", Message="Done"))
    }
    
    if(exists("exp_gwas")) rm(exp_gwas)
    if(exists("clumped_snps")) rm(clumped_snps)
    # 强制回收内存
    gc(verbose = FALSE)
    
  } # 结束 current_exp 循环
  rm(out_gwas, out_raw, mapping)
  gc(verbose = FALSE)
  
  end_time <- Sys.time()
  message(sprintf(">>> Outcome [%s] 处理完毕，耗时: %.2f 分钟", 
                  current_out, as.numeric(difftime(end_time, start_time, units="mins"))))
} # 结束 current_out 循环

# 保存日志
fwrite(log_df, log_file)