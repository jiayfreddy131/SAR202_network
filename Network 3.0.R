####### Create T0-30 order MFs networks Mar 22 #########

library(vegan)
library(reshape)
library(igraph)
library(fdrtool)
library(ggplot2)
library(Hmisc)
library(data.table)
library(dplyr)
library(rgexf)
library(psych)

## Import data from Network analysi 2.0
library(readxl)
network_raw_3 <- read_excel("network_raw_cleaned.xlsx")

# select data for T0-30
network_3 <- network_raw_3[1:5,]


## Create correlation matrix
matrix_dist_3 <- rcorr(as.matrix(network_3), type = "spearman")
matrix_corr_3 <- matrix_dist_3$r
matrix_p_3 <- matrix_dist_3$P
matrix_p_3 <- p.adjust(matrix_p_3, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_3[which(matrix_corr_3>-0.7 & matrix_corr_3<0.7)] = 0
matrix_corr_3[which(matrix_p_3>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_3[is.nan(matrix_corr_3)] = 0
# delete variables with no correlations
# matrix_corr_3 <- matrix_corr_3[which(rowSums(matrix_corr_3)!=0),]
# matrix_corr_3 <- matrix_corr_3[, which(colSums(matrix_corr_3)!=0)]

## Need to filter out the OTU/OTU and MF/MF corrs, or the edge table result in 600000 rows
# transform matrix corr to data frame to preserve OTU/MF corrs
matrix_test <- as.data.frame(matrix_corr_3)
matrix_test_OTU_row <- matrix_test[1:177,] # 177 rows of OTU
matrix_test_OTU_row_MF_col <- matrix_test_OTU_row[,178:6190]
matrix_test_matrix <- as.matrix(matrix_test_OTU_row_MF_col)
# see if a igraph can be generated
g_3_test <- graph.adjacency(matrix_test_matrix, weighted = T, mode = "upper", diag = F)
# NOPE! must be a square matrix, can't work around the issue here

## Create network attributes
g_network_3 <- graph.adjacency(matrix_corr_3, weighted = T, mode = "upper", diag = F)
e_3 <- get.edgelist(g_network_3)
edges_3 <- as.data.frame(cbind(e_3,E(g_network_3)$weight)) 
colnames(edges_3) <- c("source", "target", "weight")

nodes1 <-as.data.frame(unique(edges_3[,1]));colnames(nodes1)<-c("id")
nodes2 <-as.data.frame(unique(edges_3[,2]));colnames(nodes2)<-c("id")
nodes3=rbind(nodes1,nodes2)
nodes <-as.data.frame(unique(nodes3[,1]));colnames(nodes)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_test <- edges_3
edges_test_filter <- dplyr::filter(edges_test, grepl(';', source))
edges_test_filter_2 <- dplyr::filter(edges_test_filter, !grepl(';', target))

edges_3_filter <- edges_test_filter_2

# rebuild igraph with filtered edges
g_3 <- graph.data.frame(d=edges_3_filter, vertices=nodes, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_3, "Network 3.0")

## Using only first 30 day data showed more separated networks
## filtered low-intensity data (<10^7)


####### Create T30-364 order MFs networks Mar 31 #########

# select data for T30-364
network_30 <- network_raw_3[5:8,]
## !!! rcorr needs at least 5 data points, so have to select T10-364
## (this data will probably be not convincing)
network_30 <- network_raw_3[4:8,]

## Create correlation matrix
matrix_dist_30 <- rcorr(as.matrix(network_30), type = "spearman")
matrix_corr_30 <- matrix_dist_30$r
matrix_p_30 <- matrix_dist_30$P
matrix_p_30 <- p.adjust(matrix_p_30, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_30[which(matrix_corr_30>-0.7 & matrix_corr_30<0.7)] = 0
matrix_corr_30[which(matrix_p_30>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_30[is.nan(matrix_corr_30)] = 0
# delete variables with no correlations
# matrix_corr_3 <- matrix_corr_3[which(rowSums(matrix_corr_3)!=0),]
# matrix_corr_3 <- matrix_corr_3[, which(colSums(matrix_corr_3)!=0)]

## Create network attributes
g_network_30 <- graph.adjacency(matrix_corr_30, weighted = T, mode = "upper", diag = F)
e_30 <- get.edgelist(g_network_30)
edges_30 <- as.data.frame(cbind(e_30,E(g_network_30)$weight)) 
colnames(edges_30) <- c("source", "target", "weight")

nodes1_30 <-as.data.frame(unique(edges_30[,1]));colnames(nodes1_30)<-c("id")
nodes2_30 <-as.data.frame(unique(edges_30[,2]));colnames(nodes2_30)<-c("id")
nodes3_30=rbind(nodes1_30,nodes2_30)
nodes_30 <-as.data.frame(unique(nodes3_30[,1]));colnames(nodes_30)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_30_filter <- dplyr::filter(edges_30, grepl(';', source))
edges_30_filter_2 <- dplyr::filter(edges_30_filter, !grepl(';', target))


# rebuild igraph with filtered edges
g_30 <- graph.data.frame(d=edges_30_filter_2, vertices=nodes_30, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_30, "Network 3.0 After 30 Days")


######## Build non-AOA network April 11 ########
## Import non_AOA rel abun network table
library(readxl)
network_raw_non_AOA <- read_excel("network_non_AOA_cleaned.xlsx")

## Create correlation matrix
matrix_dist_non_AOA <- rcorr(as.matrix(network_raw_non_AOA), type = "spearman")
matrix_corr_non_AOA <- matrix_dist_non_AOA$r
matrix_p_non_AOA <- matrix_dist_non_AOA$P
matrix_p_non_AOA <- p.adjust(matrix_p_non_AOA, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_non_AOA[which(matrix_corr_non_AOA>-0.7 & matrix_corr_non_AOA<0.7)] = 0
matrix_corr_non_AOA[which(matrix_p_non_AOA>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_non_AOA[is.nan(matrix_corr_non_AOA)] = 0

## Create network attributes
g_network_non_AOA <- graph.adjacency(matrix_corr_non_AOA, weighted = T, mode = "upper", diag = F)
e_non_AOA <- get.edgelist(g_network_non_AOA)
edges_non_AOA <- as.data.frame(cbind(e_non_AOA,E(g_network_non_AOA)$weight)) 
colnames(edges_non_AOA) <- c("source", "target", "weight")

nodes1_non_AOA <-as.data.frame(unique(edges_non_AOA[,1]));colnames(nodes1_non_AOA)<-c("id")
nodes2_non_AOA <-as.data.frame(unique(edges_non_AOA[,2]));colnames(nodes2_non_AOA)<-c("id")
nodes3_non_AOA=rbind(nodes1_non_AOA,nodes2_non_AOA)
nodes_non_AOA <-as.data.frame(unique(nodes3_non_AOA[,1]));colnames(nodes_non_AOA)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_non_AOA_filter <- dplyr::filter(edges_non_AOA, grepl(';', source))
edges_non_AOA_filter_2 <- dplyr::filter(edges_non_AOA_filter, !grepl(';', target))


# rebuild igraph with filtered edges
g_non_AOA <- graph.data.frame(d=edges_non_AOA_filter_2, vertices=nodes_non_AOA, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_non_AOA, "Network 3.0 AOA filtered")


######## Build Control Network Non-AOA April 14 ########

## Import Control rel abun network table
library(readxl)
network_control <- read_excel("network_control.xlsx")

## Create correlation matrix
matrix_dist_control <- rcorr(as.matrix(network_control), type = "spearman")
matrix_corr_control <- matrix_dist_control$r
matrix_p_control <- matrix_dist_control$P
matrix_p_control <- p.adjust(matrix_p_control, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_control[which(matrix_corr_control>-0.7 & matrix_corr_control<0.7)] = 0
matrix_corr_control[which(matrix_p_control>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_control[is.nan(matrix_corr_control)] = 0

## Create network attributes
g_network_control <- graph.adjacency(matrix_corr_control, weighted = T, mode = "upper", diag = F)
e_control <- get.edgelist(g_network_control)
edges_control <- as.data.frame(cbind(e_control,E(g_network_control)$weight)) 
colnames(edges_control) <- c("source", "target", "weight")

nodes1_control <-as.data.frame(unique(edges_control[,1]));colnames(nodes1_control)<-c("id")
nodes2_control <-as.data.frame(unique(edges_control[,2]));colnames(nodes2_control)<-c("id")
nodes3_control=rbind(nodes1_control,nodes2_control)
nodes_control <-as.data.frame(unique(nodes3_control[,1]));colnames(nodes_control)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_control_filter <- dplyr::filter(edges_control, grepl(';', source))
edges_control_filter_2 <- dplyr::filter(edges_control_filter, !grepl(';', target))


# rebuild igraph with filtered edges
g_control <- graph.data.frame(d=edges_control_filter_2, vertices=nodes_control, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_control, "Network 3.0 Control Network fixed")


### Just realized the OTU rel abun was taken from AOA-filtered table
# assume this won't make much difference, but will need to build one with AOA

######## Build Treatment Phase I (0-30) & II (30-364) Networks without AOA April 22 ########

### Treatment non-AOA T0-30
## Import rel abun network table
library(readxl)
network_non_AOA_phase_1 <- read_excel("network_non_AOA_0-30.xlsx")

## Create correlation matrix
matrix_dist_non_AOA_phase_1 <- rcorr(as.matrix(network_non_AOA_phase_1), type = "spearman")
matrix_corr_non_AOA_phase_1 <- matrix_dist_non_AOA_phase_1$r
matrix_p_non_AOA_phase_1 <- matrix_dist_non_AOA_phase_1$P
matrix_p_non_AOA_phase_1 <- p.adjust(matrix_p_non_AOA_phase_1, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_non_AOA_phase_1[which(matrix_corr_non_AOA_phase_1>-0.7 & matrix_corr_non_AOA_phase_1<0.7)] = 0
matrix_corr_non_AOA_phase_1[which(matrix_p_non_AOA_phase_1>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_non_AOA_phase_1[is.nan(matrix_corr_non_AOA_phase_1)] = 0

## Create network attributes
g_network_non_AOA_phase_1 <- graph.adjacency(matrix_corr_non_AOA_phase_1, weighted = T, mode = "upper", diag = F)
e_non_AOA_phase_1 <- get.edgelist(g_network_non_AOA_phase_1)
edges_non_AOA_phase_1 <- as.data.frame(cbind(e_non_AOA_phase_1,E(g_network_non_AOA_phase_1)$weight)) 
colnames(edges_non_AOA_phase_1) <- c("source", "target", "weight")

nodes1_non_AOA_phase_1 <-as.data.frame(unique(edges_non_AOA_phase_1[,1]));colnames(nodes1_non_AOA_phase_1)<-c("id")
nodes2_non_AOA_phase_1 <-as.data.frame(unique(edges_non_AOA_phase_1[,2]));colnames(nodes2_non_AOA_phase_1)<-c("id")
nodes3_non_AOA_phase_1=rbind(nodes1_non_AOA_phase_1,nodes2_non_AOA_phase_1)
nodes_non_AOA_phase_1 <-as.data.frame(unique(nodes3_non_AOA_phase_1[,1]));colnames(nodes_non_AOA_phase_1)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_non_AOA_phase_1_filter <- dplyr::filter(edges_non_AOA_phase_1, grepl(';', source))
edges_non_AOA_phase_1_filter_2 <- dplyr::filter(edges_non_AOA_phase_1_filter, !grepl(';', target))


# rebuild igraph with filtered edges
g_non_AOA_phase_1 <- graph.data.frame(d=edges_non_AOA_phase_1_filter_2, vertices=nodes_non_AOA_phase_1, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_non_AOA_phase_1, "Network 3.0 non-AOA T0-30")


### Treatment non-AOA T10-364

## Import Control rel abun network table
library(readxl)
network_non_AOA_phase_2 <- read_excel("network_non_AOA_10-364.xlsx")

## Create correlation matrix
matrix_dist_non_AOA_phase_2 <- rcorr(as.matrix(network_non_AOA_phase_2), type = "spearman")
matrix_corr_non_AOA_phase_2 <- matrix_dist_non_AOA_phase_2$r
matrix_p_non_AOA_phase_2 <- matrix_dist_non_AOA_phase_2$P
matrix_p_non_AOA_phase_2 <- p.adjust(matrix_p_non_AOA_phase_2, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_non_AOA_phase_2[which(matrix_corr_non_AOA_phase_2>=-0.7 & matrix_corr_non_AOA_phase_2<=0.7)] = 0
matrix_corr_non_AOA_phase_2[which(matrix_p_non_AOA_phase_2>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_non_AOA_phase_2[is.nan(matrix_corr_non_AOA_phase_2)] = 0

## Create network attributes
g_network_non_AOA_phase_2 <- graph.adjacency(matrix_corr_non_AOA_phase_2, weighted = T, mode = "upper", diag = F)
e_non_AOA_phase_2 <- get.edgelist(g_network_non_AOA_phase_2)
edges_non_AOA_phase_2 <- as.data.frame(cbind(e_non_AOA_phase_2,E(g_network_non_AOA_phase_2)$weight)) 
colnames(edges_non_AOA_phase_2) <- c("source", "target", "weight")

nodes1_non_AOA_phase_2 <-as.data.frame(unique(edges_non_AOA_phase_2[,1]));colnames(nodes1_non_AOA_phase_2)<-c("id")
nodes2_non_AOA_phase_2 <-as.data.frame(unique(edges_non_AOA_phase_2[,2]));colnames(nodes2_non_AOA_phase_2)<-c("id")
nodes3_non_AOA_phase_2=rbind(nodes1_non_AOA_phase_2,nodes2_non_AOA_phase_2)
nodes_non_AOA_phase_2 <-as.data.frame(unique(nodes3_non_AOA_phase_2[,1]));colnames(nodes_non_AOA_phase_2)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_non_AOA_phase_2_filter <- dplyr::filter(edges_non_AOA_phase_2, grepl(';', source))
edges_non_AOA_phase_2_filter_2 <- dplyr::filter(edges_non_AOA_phase_2_filter, !grepl(';', target))


# rebuild igraph with filtered edges
g_non_AOA_phase_2 <- graph.data.frame(d=edges_non_AOA_phase_2_filter_2, vertices=nodes_non_AOA_phase_2, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_non_AOA_phase_2, "Network 3.0 non-AOA T10-364")


######## Build Control Phase I (0-30) & II (30-364) Networks without AOA May 9 ########

### Control non-AOA T0-30
## Import rel abun network table
library(readxl)
network_non_AOA_con_phase_1 <- read_excel("network_control_non_AOA_0-30.xlsx")

## Create correlation matrix
matrix_dist_non_AOA_con_phase_1 <- rcorr(as.matrix(network_non_AOA_con_phase_1), type = "spearman")
matrix_corr_non_AOA_con_phase_1 <- matrix_dist_non_AOA_con_phase_1$r
matrix_p_non_AOA_con_phase_1 <- matrix_dist_non_AOA_con_phase_1$P
matrix_p_non_AOA_con_phase_1 <- p.adjust(matrix_p_non_AOA_con_phase_1, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_non_AOA_con_phase_1[which(matrix_corr_non_AOA_con_phase_1>-0.7 & matrix_corr_non_AOA_con_phase_1<0.7)] = 0
matrix_corr_non_AOA_con_phase_1[which(matrix_p_non_AOA_con_phase_1>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_non_AOA_con_phase_1[is.nan(matrix_corr_non_AOA_con_phase_1)] = 0

## Create network attributes
g_network_non_AOA_con_phase_1 <- graph.adjacency(matrix_corr_non_AOA_con_phase_1, weighted = T, mode = "upper", diag = F)
e_non_AOA_con_phase_1 <- get.edgelist(g_network_non_AOA_con_phase_1)
edges_non_AOA_con_phase_1 <- as.data.frame(cbind(e_non_AOA_con_phase_1,E(g_network_non_AOA_con_phase_1)$weight)) 
colnames(edges_non_AOA_con_phase_1) <- c("source", "target", "weight")

nodes1_non_AOA_con_phase_1 <-as.data.frame(unique(edges_non_AOA_con_phase_1[,1]));colnames(nodes1_non_AOA_con_phase_1)<-c("id")
nodes2_non_AOA_con_phase_1 <-as.data.frame(unique(edges_non_AOA_con_phase_1[,2]));colnames(nodes2_non_AOA_con_phase_1)<-c("id")
nodes3_non_AOA_con_phase_1=rbind(nodes1_non_AOA_con_phase_1,nodes2_non_AOA_con_phase_1)
nodes_non_AOA_con_phase_1 <-as.data.frame(unique(nodes3_non_AOA_con_phase_1[,1]));colnames(nodes_non_AOA_con_phase_1)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_non_AOA_con_phase_1_filter <- dplyr::filter(edges_non_AOA_con_phase_1, grepl(';', source))
edges_non_AOA_con_phase_1_filter_2 <- dplyr::filter(edges_non_AOA_con_phase_1_filter, !grepl(';', target))


# rebuild igraph with filtered edges
g_non_AOA_con_phase_1 <- graph.data.frame(d=edges_non_AOA_con_phase_1_filter_2, vertices=nodes_non_AOA_con_phase_1, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_non_AOA_con_phase_1, "Network 3.0 non-AOA Control T0-30")


### Control non-AOA T10-364

## Import Control rel abun network table
library(readxl)
network_non_AOA_con_phase_2 <- read_excel("network_control_non_AOA_10-364.xlsx")

## Create correlation matrix
matrix_dist_non_AOA_con_phase_2 <- rcorr(as.matrix(network_non_AOA_con_phase_2), type = "spearman")
matrix_corr_non_AOA_con_phase_2 <- matrix_dist_non_AOA_con_phase_2$r
matrix_p_non_AOA_con_phase_2 <- matrix_dist_non_AOA_con_phase_2$P
matrix_p_non_AOA_con_phase_2 <- p.adjust(matrix_p_non_AOA_con_phase_2, method = "BH")

# select strong correlation (|r| >= 0.7) and significant (p < 0.05)
matrix_corr_non_AOA_con_phase_2[which(matrix_corr_non_AOA_con_phase_2>=-0.7 & matrix_corr_non_AOA_con_phase_2<=0.7)] = 0
matrix_corr_non_AOA_con_phase_2[which(matrix_p_non_AOA_con_phase_2>=0.05)] = 0

# delete the NaN (no need to filter diag 1 because it will be removed by the arg diag = F)
matrix_corr_non_AOA_con_phase_2[is.nan(matrix_corr_non_AOA_con_phase_2)] = 0

## Create network attributes
g_network_non_AOA_con_phase_2 <- graph.adjacency(matrix_corr_non_AOA_con_phase_2, weighted = T, mode = "upper", diag = F)
e_non_AOA_con_phase_2 <- get.edgelist(g_network_non_AOA_con_phase_2)
edges_non_AOA_con_phase_2 <- as.data.frame(cbind(e_non_AOA_con_phase_2,E(g_network_non_AOA_con_phase_2)$weight)) 
colnames(edges_non_AOA_con_phase_2) <- c("source", "target", "weight")

nodes1_non_AOA_con_phase_2 <-as.data.frame(unique(edges_non_AOA_con_phase_2[,1]));colnames(nodes1_non_AOA_con_phase_2)<-c("id")
nodes2_non_AOA_con_phase_2 <-as.data.frame(unique(edges_non_AOA_con_phase_2[,2]));colnames(nodes2_non_AOA_con_phase_2)<-c("id")
nodes3_non_AOA_con_phase_2=rbind(nodes1_non_AOA_con_phase_2,nodes2_non_AOA_con_phase_2)
nodes_non_AOA_con_phase_2 <-as.data.frame(unique(nodes3_non_AOA_con_phase_2[,1]));colnames(nodes_non_AOA_con_phase_2)<-c("id")

## delete MF interactions (takes up so much edges)
# testing dply filter to remove MF/MF and OTU/OTU corr
edges_non_AOA_con_phase_2_filter <- dplyr::filter(edges_non_AOA_con_phase_2, grepl(';', source))
edges_non_AOA_con_phase_2_filter_2 <- dplyr::filter(edges_non_AOA_con_phase_2_filter, !grepl(';', target))


# rebuild igraph with filtered edges
g_non_AOA_con_phase_2 <- graph.data.frame(d=edges_non_AOA_con_phase_2_filter_2, vertices=nodes_non_AOA_con_phase_2, directed=T)

library(RCy3)
cytoscapePing()
createNetworkFromIgraph(g_non_AOA_con_phase_2, "Network 3.0 non-AOA Control T10-364")

