library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(openxlsx)

# --- 1. 环境准备 ---
if (!dir.exists(res_base_folder)) {
  dir.create(res_base_folder, recursive = TRUE)
}

# 设置日志文件路径
log_file <- file.path(res_base_folder, paste0("analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

# 【关键】启动日志记录
# split = TRUE 表示同时在控制台显示和写入文件
sink(log_file, append = FALSE, split = TRUE)

# 打印日志头部信息
cat("================================================\n")
cat("任务启动时间:", as.character(Sys.time()), "\n")
cat("输入文件:", input_file, "\n")
cat("================================================\n\n")

# 使用 tryCatch 确保即使出错也能正常关闭 sink
tryCatch({
  
  # --- 2. 读取数据 ---
  all_data <- read.csv(input_file, stringsAsFactors = FALSE)
  root_list <- unique(all_data$Root)
  cat("检测到共有", length(root_list), "个疾病需要进行富集分析。\n")
  
  # --- 3. 循环遍历每个疾病 ---
  for (root in root_list) {
    cat("\n[", as.character(Sys.time()), "] 正在处理:", root, "...\n")
    
    # ... (这里保持你原来的分析代码不变) ...
    protein_names <- all_data %>% 
      filter(Root == root) %>% 
      pull(Protein) %>% 
      na.omit() %>% 
      unique()
    
    protein_names <- protein_names[protein_names != ""]
    
    if (length(protein_names) < 10) {
      cat(">>> 跳过:", dis, "(蛋白质数量 ", length(protein_names), " < 10)\n")
      next
    }
    
    # ID 转换
    gene_df <- tryCatch({
      bitr(protein_names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat(">>> 错误: ID转换失败 (", e$message, ")\n")
      return(NULL)
    })
    
    if (is.null(gene_df) || nrow(gene_df) == 0) next
    
    entrez_ids <- unique(gene_df$ENTREZID)
    
    # 富集分析 (GO/KEGG/Reactome)
    # 注意：为了防止报错中断日志，建议给每个分析也加一层 tryCatch
    GO_all <- tryCatch({ enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE) }, error = function(e) NULL)
    KEGG <- tryCatch({ enrichKEGG(gene = entrez_ids, organism = "hsa") }, error = function(e) NULL)
    REACTOME <- tryCatch({ enrichPathway(gene = entrez_ids, organism = "human", readable = TRUE) }, error = function(e) NULL)
    
    # 保存 Excel
    safe_dis_name <- gsub("[^[:alnum:]]", "_", root)
    res_path <- file.path(res_base_folder, paste0(safe_dis_name, ".xlsx"))
    
    wb <- createWorkbook()
    results_list <- list(GO_all = GO_all, KEGG = KEGG, REACTOME = REACTOME)
    
    has_data <- FALSE
    for (name in names(results_list)) {
      if (!is.null(results_list[[name]])) {
        df <- as.data.frame(results_list[[name]])
        if (nrow(df) > 0) {
          addWorksheet(wb, name)
          writeData(wb, name, df)
          has_data <- TRUE
        }
      }
    }
    
    if (has_data) {
      saveWorkbook(wb, res_path, overwrite = TRUE)
      cat(">>> 成功: 结果已保存至 Excel。\n")
    } else {
      cat(">>> 提示: 无显著富集项。\n")
    }
  }
  
}, error = function(e) {
  cat("\n程序发生严重错误:\n", e$message, "\n")
}, finally = {
  cat("\n================================================\n")
  cat("任务结束时间:", as.character(Sys.time()), "\n")
  cat("================================================\n")
  
  # 【关键】关闭日志记录
  sink() 
})

cat("\n所有任务已执行完毕！Log文件保存在:", log_file, "\n")

