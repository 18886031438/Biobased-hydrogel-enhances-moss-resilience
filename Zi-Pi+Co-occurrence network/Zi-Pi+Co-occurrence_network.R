#=================Nodes==================
# 加载必要的包
library(dplyr)
library(tidyr)
library(readxl)

# 设置工作目录（修改为您的文件路径）
setwd("D:\\科研项目\\屋顶绿化\\2025-11-19\\模块三\\3.ZI-PI筛选\\绘图数据")

# 1. 读取文件
key_nodes <- read_excel("RJG_1.xlsx")  # 关键网络节点文件

setwd("D:\\科研项目\\屋顶绿化\\2025-11-19\\模块三\\3.ZI-PI筛选\\总数据-拆分")
microbes <- read_excel("RJ_核心微生物-G.xlsx")         # 微生物文件
genes <- read_excel("RJG-关键基因.xlsx")              # 基因文件
mges <- read_excel("RJ_核心-G.xlsx")               # MGEs文件
network <- read_excel("RJG.xlsx")      # 网络图文件

# 2. 确保所有节点名称都是字符类型
key_nodes[] <- lapply(key_nodes, as.character)
microbes[] <- lapply(microbes, as.character)
genes[] <- lapply(genes, as.character)
mges[] <- lapply(mges, as.character)
network[] <- lapply(network, as.character)

# 3. 根据关键网络节点筛选网络图中的行
# 使用关键节点文件的第一列(nodes_id)筛选网络文件
filtered_network <- network %>%
  filter(Data1 %in% key_nodes[[1]] & Data2 %in% key_nodes[[1]])

# 4. 添加Etype列（考虑方向）
# 首先判断每个节点的类型
get_node_type <- function(node_name) {
  if (node_name %in% microbes[[1]]) {
    return("F")  # 微生物
  } else if (node_name %in% genes[[1]]) {
    return("T")  # 基因
  } else if (node_name %in% mges[[1]]) {
    return("S")  # MGEs
  } else {
    return(NA)
  }
}

# 应用函数获取节点类型
filtered_network_with_type <- filtered_network %>%
  mutate(
    Type1 = sapply(Data1, get_node_type),
    Type2 = sapply(Data2, get_node_type)
  )

# 根据方向组合确定Etype
filtered_network_with_type <- filtered_network_with_type %>%
  mutate(
    Etype = case_when(
      Type1 == "F" & Type2 == "T" ~ "FT",  # 微生物→基因
      Type1 == "T" & Type2 == "F" ~ "TF",  # 基因→微生物
      Type1 == "F" & Type2 == "S" ~ "FS",  # 微生物→MGEs
      Type1 == "S" & Type2 == "F" ~ "SF",  # MGEs→微生物
      Type1 == "T" & Type2 == "S" ~ "TS",  # 基因→MGEs
      Type1 == "S" & Type2 == "T" ~ "ST",  # MGEs→基因
      Type1 == Type2 ~ paste0(Type1, Type2),  # 同类型连接，如FF、TT、SS
      TRUE ~ "Unknown"  # 未知类型
    )
  ) %>%
  select(-Type1, -Type2)  # 移除中间列

# 5. 创建节点文件
# 获取所有唯一的节点
all_nodes <- unique(c(filtered_network_with_type$Data1, filtered_network_with_type$Data2))

# 创建节点数据框
node_list <- data.frame(
  Id = all_nodes,
  stringsAsFactors = FALSE
)

# 添加shape和Label列
node_list <- node_list %>%
  mutate(
    shape = sapply(Id, get_node_type),
    Label = Id
  ) %>%
  arrange(Id)  # 按节点名称排序

# 6. 保存结果
setwd("D:\\科研项目\\屋顶绿化\\2025-11-19\\模块三\\3.ZI-PI筛选\\绘图数据")
write.csv(filtered_network_with_type, "RJG_带Etype.csv", row.names = FALSE)
write.csv(node_list, "RJG-分组.csv", row.names = FALSE)

# 7. 打印结果预览
cat("========== 提取的网络图（前10行） ==========\n")
print(head(filtered_network_with_type, 10))

cat("\n========== 节点列表（前10行） ==========\n")
print(head(node_list, 10))

# 8. 统计信息
cat("\n========== 各类型边的数量 ==========\n")
edge_counts <- table(filtered_network_with_type$Etype, useNA = "ifany")
print(edge_counts)








#=====================Zi-Pi===============
##定义函数
zi.pi<-function(nodes_bulk, z.bulk, modularity_class, degree){
  
  z.bulk[abs(z.bulk)>0]<-1
  module<-which(colnames(nodes_bulk)==modularity_class)
  module.max<-max(nodes_bulk[,module])
  degree<-which(colnames(nodes_bulk)==degree)
  
  #按照模块将相关矩阵分割
  bulk.module<-list(NA)
  length(bulk.module)<-module.max
  
  for(i in 1:max(nodes_bulk[,module])){
    bulk.module[[i]]<-z.bulk[which(nodes_bulk[,module]==i),which(nodes_bulk[,module]==i)]
    bulk.module[[i]]<-as.data.frame(bulk.module[[i]])
    rownames(bulk.module[[i]])<-rownames(z.bulk)[which(nodes_bulk[,module]==i)]
    colnames(bulk.module[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
  }
  
  # within-module degree z
  z_bulk<-list(NA)
  length(z_bulk)<-module.max
  
  for(i in 1:length(z_bulk)){
    z_bulk[[i]]<-bulk.module[[i]][,1]
    z_bulk[[i]]<-as.data.frame(z_bulk[[i]])
    colnames(z_bulk[[i]])<-"z"
    rownames(z_bulk[[i]])<-rownames(bulk.module[[i]])
  }
  
  #计算z值
  for(i in 1:max(nodes_bulk[,module])){
    if(length(bulk.module[[i]])==1){
      z_bulk[[i]][,1]<-0
    }else if(sum(bulk.module[[i]])==0){
      z_bulk[[i]][,1]<-0
    }else{
      k<-rowSums(bulk.module[[i]]) -1
      mean<-mean(k)
      sd<-sd(k)
      if (sd==0){
        z_bulk[[i]][,1]<-0
      }else{
        z_bulk[[i]][,1]<-(k-mean)/sd
      }
    }
  }
  
  #z值合并
  for(i in 2:max(nodes_bulk[,module])) {
    z_bulk[[i]]<-rbind(z_bulk[[i-1]],z_bulk[[i]])
  }
  z_bulk<-z_bulk[[module.max]]
  
  #按照模块将相关矩阵列分割
  bulk.module1<-list(NA)
  length(bulk.module1)<-module.max
  
  for(i in 1:max(nodes_bulk[,module])){
    bulk.module1[[i]]<-z.bulk[,which(nodes_bulk[,module]==i)]
    bulk.module1[[i]]<-as.data.frame(bulk.module1[[i]])
    rownames(bulk.module1[[i]])<-rownames(z.bulk)
    colnames(bulk.module1[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
  }
  
  #among-module connectivity c
  c_bulk<-list(NA)
  length(c_bulk)<-module.max
  
  for(i in 1:length(c_bulk)){
    c_bulk[[i]]<-z.bulk[,1]
    c_bulk[[i]]<-as.matrix(c_bulk[[i]])
    colnames(c_bulk[[i]])<-"c"
    rownames(c_bulk[[i]])<-rownames(z.bulk)
    c_bulk[[i]][,1]<-NA
  }
  
  #每个节点各模块连接数平方
  for(i in 1:max(nodes_bulk[,module])){
    c_bulk[[i]]<-rowSums(bulk.module1[[i]])
    c_bulk[[i]]<-as.matrix(c_bulk[[i]])
    c_bulk[[i]][which(nodes_bulk$modularity == i),] = c_bulk[[i]][which(nodes_bulk$modularity == i),] -1
    c_bulk[[i]]<-c_bulk[[i]]*c_bulk[[i]]
    colnames(c_bulk[[i]])<-"c"
    rownames(c_bulk[[i]])<-rownames(z.bulk)
  }
  
  #平方和
  for(i in 2:max(nodes_bulk[,module])){
    c_bulk[[i]]<-c_bulk[[i]]+c_bulk[[i-1]]
  }
  c_bulk<-c_bulk[[module.max]]
  
  for(i in 1: length(c_bulk)){
    if(nodes_bulk$degree[i]==1){
      c_bulk[i] <- 0
    }else{
      c_bulk[i] <- 1-(c_bulk[i]/(nodes_bulk$degree[i]*nodes_bulk$degree[i]))
    }
  }
  colnames(c_bulk)<-"c"
  
  #z,c整合
  z_c_bulk<-c_bulk
  z_c_bulk<-as.data.frame(z_c_bulk)
  z_c_bulk$z<-z_bulk[match(rownames(c_bulk),rownames(z_bulk)),]
  z_c_bulk<-z_c_bulk[,c(2,1)]
  names(z_c_bulk)[1:2]<-c('within_module_connectivities','among_module_connectivities')
  
  z_c_bulk$nodes_id<-rownames(z_c_bulk)
  nodes_bulk$nodes_id<-rownames(nodes_bulk)
  z_c_bulk<-merge(z_c_bulk,nodes_bulk,by='nodes_id')
  z_c_bulk
  
}


#===============igraph 包计算网络模块
##igraph 包计算网络模块
library(igraph)

#输入数据示例，邻接矩阵
#这是一个微生物互作网络，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
setwd("D:\\科研项目\\屋顶绿化\\2025-11-19\\模块三\\3.ZI-PI筛选\\ZI-PI")
# 1. 读取原始数据（使用check.names = FALSE避免R自动修改）
raw_data <- read.csv('RJG互作表.csv', row.names = 1, check.names = FALSE)

# 2. 检查并修复行列名中的纯数字
# 函数：检查字符串是否为纯数字
is_pure_number <- function(x) {
  grepl("^[0-9]+$", as.character(x))
}

# 函数：为纯数字添加前缀
add_prefix_to_numbers <- function(names_vector, prefix = "A_") {
  sapply(names_vector, function(name) {
    if (is_pure_number(name)) {
      return(paste0(prefix, name))
    } else {
      return(as.character(name))
    }
  })
}

# 获取当前的行名和列名
current_row_names <- rownames(raw_data)
current_col_names <- colnames(raw_data)

# 检查是否有纯数字的行列名
num_rows_pure <- sum(is_pure_number(current_row_names))
num_cols_pure <- sum(is_pure_number(current_col_names))

cat("检查结果:\n")
cat("- 行名中纯数字的数量:", num_rows_pure, "\n")
cat("- 列名中纯数字的数量:", num_cols_pure, "\n")

# 如果有纯数字的行列名，添加前缀
if (num_rows_pure > 0 || num_cols_pure > 0) {
  cat("正在为纯数字的行列名添加'A_'前缀...\n")
  
  # 修改行名
  new_row_names <- add_prefix_to_numbers(current_row_names, "A_")
  rownames(raw_data) <- new_row_names
  
  # 修改列名
  new_col_names <- add_prefix_to_numbers(current_col_names, "A_")
  colnames(raw_data) <- new_col_names
  
  # 显示修改前后的对比
  cat("\n修改前后的行列名对比示例:\n")
  cat("修改前（前5个）:", head(current_row_names, 5), "\n")
  cat("修改后（前5个）:", head(new_row_names, 5), "\n")
  
  # 如果有纯数字的行名，显示更多示例
  if (num_rows_pure > 0) {
    pure_number_rows <- current_row_names[is_pure_number(current_row_names)]
    cat("\n纯数字行名示例（修改前）:", head(pure_number_rows, 10), "\n")
    cat("对应修改后:", head(new_row_names[is_pure_number(current_row_names)], 10), "\n")
  }
  
  cat("\n行列名修改完成!\n")
} else {
  cat("没有发现纯数字的行列名，无需修改。\n")
}

# 3. 将处理后的数据赋值给adjacency_unweight
adjacency_unweight <- raw_data

# 4. 显示处理后的数据结构
cat("\n处理后的数据维度:", dim(adjacency_unweight), "\n")
cat("处理后的行名示例（前10个）:", head(rownames(adjacency_unweight), 10), "\n")
cat("处理后的列名示例（前10个）:", head(colnames(adjacency_unweight), 10), "\n")

# 5. 继续后续分析（原代码保持不变）
cat("\n--- 开始网络分析 ---\n")
head(adjacency_unweight)[1:6]    #邻接矩阵类型的网络文件

#邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph    #igraph 的邻接列表

#计算节点度
V(igraph)$degree <- degree(igraph)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块

write.table(nodes_list, 'RJG-nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')

#上述的邻接矩阵类型的网络文件
adjacency_unweight <- raw_data

#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('RJG-nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

write.table(zi_pi, 'RJG-zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
write.table(zi_pi, 'RJG-zi_pi_result1.txt', sep = '\t', row.names = FALSE, quote = FALSE)

p<- ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 15) +
  scale_color_manual(values = c('#67D5B5','#EE7785','#C89EC4','#84B1ED'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)
p
ggsave(filename = "RJG.pdf", width = 25, height = 20, units = "cm",dpi = 300)


cat("\n========== 各类型节点的数量 ==========\n")
node_counts <- table(node_list$shape)
print(node_counts)

# 9. 汇总信息
cat("\n========== 汇总信息 ==========\n")
cat(paste("原始网络图行数：", nrow(network), "\n"))
cat(paste("提取后网络图行数：", nrow(filtered_network_with_type), "\n"))
cat(paste("提取比例：", round(nrow(filtered_network_with_type)/nrow(network)*100, 2), "%\n"))
cat(paste("提取的唯一节点数：", nrow(node_list), "\n"))

cat("\n边类型分布：\n")
edge_summary <- as.data.frame(edge_counts)
colnames(edge_summary) <- c("Etype", "Count")
print(edge_summary)

cat("\n节点类型分布：\n")
node_summary <- as.data.frame(node_counts)
colnames(node_summary) <- c("Shape", "Count")
print(node_summary)

