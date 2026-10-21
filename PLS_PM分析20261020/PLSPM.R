#按照 Treat (处理) → Water (水分) → Stress (胁迫) → Micro (微生物)→ Grow (生长) 的顺序排列
#=================== 路径顺序：处理→水分→胁迫→微生物→生长 ===================
# ========== 0. 加载包 ==========
library(plspm)
library(car)

# ========== 1. 读取与预处理 ==========
raw_data <- read.csv("PLS_SM.csv", stringsAsFactors = FALSE)
numeric_vars <- c("Richness", "Shannon", "PCA1_micro",
                  "RWC", "Chl", "Biomass",
                  "Sugar", "Protein", "Proline", "MDA",
                  "SOD", "POD", "CAT", "FvFm",
                  "Gene", "MGEs", "Microbe",
                  "Total_edges", "Total_nodes", "Keystone_gene_links")
data_std <- as.data.frame(scale(raw_data[, numeric_vars]))
data_std$Treatment <- raw_data$Treatment
rownames(data_std) <- raw_data$Sample

# ========== 2. 5潜变量级联模型（新顺序） ==========
# 统一潜变量顺序：1.Treat  2.Water  3.Stress  4.Micro  5.Grow
# 所有模块严格按此顺序对应，杜绝错位

# 2.1 结构模型：纯链式级联，无越级路径
# 通路：Treat(处理) → Water(水分) → Stress(胁迫) → Micro(微生物) → Grow(生长)
path_mat <- matrix(
  c(0, 0, 0, 0, 0,   # 第1行 Treat：外生，无输入
    1, 0, 0, 0, 0,   # 第2行 Water ← Treat  处理提升水分
    0, 1, 0, 0, 0,   # 第3行 Stress ← Water  水分越充足，胁迫越低
    0, 1, 0, 0, 0,   # 第4行 Micro ← Stress  胁迫程度调控微生物组
    0, 0, 1, 1, 0),  # 第5行 Grow ← Micro  微生物组影响生长状态
  nrow = 5, byrow = TRUE,
  dimnames = list(
    c("Treat", "Water", "Stress", "Micro", "Grow"),
    c("Treat", "Water", "Stress", "Micro", "Grow")
  )
)

# 可视化验证路径方向
innerplot(path_mat)

# 2.2 测量模型（顺序与路径矩阵完全一一对应！）
blocks <- list(
  # 1. 处理分组：单指标构成型
  Treat  = "Treatment",
  # 2. 水分状态：单指标构成型
  Water  = "RWC",
  # 3. 氧化渗透胁迫：6指标反映型
  Stress = c("SOD", "POD", "CAT", "MDA", "Proline", "Sugar"),
  # 4. 微生物组响应：4指标构成型
  Micro  = c("PCA1_micro", "Total_nodes", "Keystone_gene_links", "Gene"),
  # 5. 光合生长：3指标反映型
  Grow   = c("Chl", "FvFm", "Biomass")
)

# 测量模式（顺序与blocks严格对应）
modes <- c("B",   # Treat：构成型
           "B",   # Water：单指标，构成型
           "A",   # Stress：多指标反映型
           "B",   # Micro：多指标构成型
           "A")   # Grow：多指标反映型

# ========== 必做自检：确认顺序完全匹配 ==========
cat("blocks潜变量顺序：", names(blocks), "\n")
cat("路径矩阵行名顺序：", rownames(path_mat), "\n")
# 两行输出的名称、顺序必须完全一致，再运行后续代码

# 2.3 运行基础模型
set.seed(123)
pls_final <- plspm(
  Data = data_std,
  path_matrix = path_mat,
  blocks = blocks,
  modes = modes,
  scaled = FALSE
)

# 核心结果
summary(pls_final)
pls_final$unidim              # 信效度：Stress、Grow有α/CR/AVE
pls_final$path_coefs          # 标准化路径系数
pls_final$inner_summary$R2    # 内生潜变量R²

# 2.4 Bootstrap 1000次置信区间
set.seed(123)
pls_boot <- plspm(
  Data = data_std,
  path_matrix = path_mat,
  blocks = blocks,
  modes = modes,
  scaled = FALSE,
  boot.val = TRUE,
  br = 1000
)
cat("\n===== 路径系数Bootstrap 95%置信区间 =====\n")
print(pls_boot$boot$paths)

# 2.5 构成型块VIF检验
calc_vif <- function(var_names, df) {
  if (length(var_names) < 3) return(NA)
  form <- as.formula(paste(var_names[1], "~", paste(var_names[-1], collapse = " + ")))
  return(max(vif(lm(form, data = df))))
}
cat("\n===== 构成型块最大VIF =====\n")
cat("Micro块最大VIF：", calc_vif(blocks$Micro, data_std), "\n")


# ========== 中介分析（匹配分叉结构：Water→Stress+Micro，两者共同→Grow） ==========
paths <- pls_final$path_coefs

# 提取所有实际存在的直接路径
a1 <- paths["Water", "Treat"]    # Treat → Water
a2 <- paths["Stress", "Water"]   # Water → Stress
a3 <- paths["Micro", "Water"]    # Water → Micro
b1 <- paths["Grow", "Stress"]    # Stress → Grow
b2 <- paths["Grow", "Micro"]     # Micro → Grow

# 计算各条间接效应
ind_treat_stress <- a1 * a2    # Treat → Water → Stress
ind_treat_micro  <- a1 * a3    # Treat → Water → Micro
ind_stress_grow  <- a2 * b1    # Water → Stress → Grow
ind_micro_grow   <- a3 * b2    # Water → Micro → Grow
ind_total_grow   <- a1*(a2*b1 + a3*b2)  # Treat对Grow的总间接效应（两条通路之和）

# Bootstrap置信区间
boot_paths <- pls_boot$boot$paths
boot_a1 <- boot_paths["Treat -> Water", ]
boot_a2 <- boot_paths["Water -> Stress", ]
boot_a3 <- boot_paths["Water -> Micro", ]
boot_b1 <- boot_paths["Stress -> Grow", ]
boot_b2 <- boot_paths["Micro -> Grow", ]

boot_ind_treat_stress <- boot_a1 * boot_a2
boot_ind_treat_micro  <- boot_a1 * boot_a3
boot_ind_stress_grow  <- boot_a2 * boot_b1
boot_ind_micro_grow   <- boot_a3 * boot_b2
boot_ind_total_grow   <- boot_a1 * (boot_a2*boot_b1 + boot_a3*boot_b2)

calc_ci <- function(boot_vec) {
  quantile(boot_vec, c(0.025, 0.975), na.rm = TRUE)
}

# 结果表
mediation_result <- data.frame(
  中介路径 = c(
    "Treat → Water → Stress",
    "Treat → Water → Micro",
    "Water → Stress → Grow",
    "Water → Micro → Grow",
    "Treat → Water → (Stress+Micro) → Grow（总间接效应）"
  ),
  效应值 = round(c(
    ind_treat_stress, ind_treat_micro,
    ind_stress_grow, ind_micro_grow,
    ind_total_grow
  ), 4),
  CI_2.5 = round(c(
    calc_ci(boot_ind_treat_stress)[1], calc_ci(boot_ind_treat_micro)[1],
    calc_ci(boot_ind_stress_grow)[1], calc_ci(boot_ind_micro_grow)[1],
    calc_ci(boot_ind_total_grow)[1]
  ), 4),
  CI_97.5 = round(c(
    calc_ci(boot_ind_treat_stress)[2], calc_ci(boot_ind_treat_micro)[2],
    calc_ci(boot_ind_stress_grow)[2], calc_ci(boot_ind_micro_grow)[2],
    calc_ci(boot_ind_total_grow)[2]
  ), 4),
  显著性 = c(
    ifelse(prod(calc_ci(boot_ind_treat_stress)) > 0, "显著", "不显著"),
    ifelse(prod(calc_ci(boot_ind_treat_micro)) > 0, "显著", "不显著"),
    ifelse(prod(calc_ci(boot_ind_stress_grow)) > 0, "显著", "不显著"),
    ifelse(prod(calc_ci(boot_ind_micro_grow)) > 0, "显著", "不显著"),
    ifelse(prod(calc_ci(boot_ind_total_grow)) > 0, "显著", "不显著")
  )
)

print(mediation_result, row.names = FALSE)









