#Title: JPB1 2mo Cytokine Outlier Detection and Removal Script
#File: JPB1 2mo Cytokine Outlier Detection and Removal Script 2023-01-04
#Authors: Liam North

#REFERRENCES###############################################################
  #https://rstudio-pubs-static.s3.amazonaws.com/30836_e99fbbe452ba481aab475f0dbcf5459b.html
  #https://rpubs.com/michaelmallari/anomaly-detection-r
  #https://www.nature.com/articles/s41598-019-46097-9

#LOAD PACKAGES#############################################################
pacman::p_load(pacman, tidyverse, FNN)

#DATA IMPORT###############################################################
#DF1 for import of cytokine csv named "01_Compiled 2MO CHLA+BCH Cytokine Data LN 2023-01-04.CSV"
df1 <- data.frame(read.csv(file.choose()))
head(df1)
tail(df1)
as_tibble(df1)

#OUTLIER DETECTION USING KNN###############################################
#IDENTIFY IFN-GAMMA OUTLIERS KNN
nn_ifngamma <- get.knn(df1$X2mo.Obs.Conc.IFNgamma..25., k=20)
head(nn_ifngamma)
# Create score by averaging distances
nn_ifngamma_d <- rowMeans(nn_ifngamma$nn.dist)
# Print row index of the most anomalous point
which.max(nn_ifngamma_d)
# Print the 5-nearest neighbor distance score
nn_ifngamma_d[1:5]
# Append the score as a new column 
df1$ifngamma_knn_score <- nn_ifngamma_d
# Visualize knn
plot(nn_ifngamma_d ~ 1, data = df1, cex = sqrt(ifngamma_knn_score), pch = 20)
ggplot(df1, aes(nn_ifngamma_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("IFN-gamma KNN Distances")+
  xlab("IFN-gamma KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 50
df1 <- within(df1, {   
  ifngamma.outlier <- NA
  ifngamma.outlier[ifngamma_knn_score == 50 | ifngamma_knn_score > 50] <- 1
  ifngamma.outlier[ifngamma_knn_score < 50] <- 0
} )

df1$ifngamma.outlier


#IDENTIFY IL-6 OUTLIERS KNN
nn_il6 <- get.knn(df1$X2mo.Obs.Conc.IL.6..57., k=20)
head(nn_il6)
# Create score by averaging distances
nn_il6_d <- rowMeans(nn_il6$nn.dist)
# Print row index of the most anomalous point
which.max(nn_il6_d)
# Print the 5-nearest neighbor distance score
nn_il6_d[1:5]
# Append the score as a new column 
df1$il6_knn_score <- nn_il6_d
# Visualize knn
plot(nn_il6_d ~ 1, data = df1, cex = sqrt(il6_knn_score), pch = 20)
ggplot(df1, aes(nn_il6_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("IL-6 KNN Distances")+
  xlab("IL-6 KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 50
df1 <- within(df1, {   
  il6.outlier <- NA
  il6.outlier[il6_knn_score == 50 | il6_knn_score > 50] <- 1
  il6.outlier[il6_knn_score < 50] <- 0
} )

df1$il6.outlier

#IDENTIFY IL-1beta OUTLIERS KNN
nn_il1beta <- get.knn(df1$X2mo.Obs.Conc.IL.1beta..46., k=20)
head(nn_il1beta)
# Create score by averaging distances
nn_il1beta_d <- rowMeans(nn_il1beta$nn.dist)
# Print row index of the most anomalous point
which.max(nn_il1beta_d)
# Print the 5-nearest neighbor distance score
nn_il1beta_d[1:5]
# Append the score as a new column 
df1$il1beta_knn_score <- nn_il1beta_d
# Visualize knn
plot(nn_il1beta_d ~ 1, data = df1, cex = sqrt(il1beta_knn_score), pch = 20)
ggplot(df1, aes(nn_il1beta_d))+
  geom_histogram(binwidth = 0.1) +
  ggtitle("IL-1beta KNN Distances")+
  xlab("IL-1beta KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 6
df1 <- within(df1, {   
  il1beta.outlier <- NA
  il1beta.outlier[il1beta_knn_score == 6 | il1beta_knn_score > 6] <- 1
  il1beta.outlier[il1beta_knn_score < 6] <- 0
} )

df1$il1beta.outlier

#IDENTIFY TNF-alpha OUTLIERS KNN
nn_tnfalpha <- get.knn(df1$X2mo.Obs.Conc.TNFalpha..75., k=20)
head(nn_tnfalpha)
# Create score by averaging distances
nn_tnfalpha_d <- rowMeans(nn_tnfalpha$nn.dist)
# Print row index of the most anomalous point
which.max(nn_tnfalpha_d)
# Print the 5-nearest neighbor distance score
nn_tnfalpha_d[1:5]
# Append the score as a new column 
df1$tnfalpha_knn_score <- nn_tnfalpha_d
# Visualize knn
plot(nn_tnfalpha_d ~ 1, data = df1, cex = sqrt(tnfalpha_knn_score), pch = 20)
ggplot(df1, aes(nn_tnfalpha_d))+
  geom_histogram(binwidth = 0.1) +
  ggtitle("TNF-alpha KNN Distances")+
  xlab("TNF-alpha KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 10
df1 <- within(df1, {   
  tnfalpha.outlier <- NA
  tnfalpha.outlier[tnfalpha_knn_score == 10 | tnfalpha_knn_score > 10] <- 1
  tnfalpha.outlier[tnfalpha_knn_score <10] <- 0
} )

df1$tnfalpha.outlier

#IDENTIFY MCP-1 OUTLIERS KNN
nn_mcp1 <- get.knn(df1$X2mo.Obs.Conc.MCP.1..67., k=20)
head(nn_mcp1)
# Create score by averaging distances
nn_mcp1_d <- rowMeans(nn_mcp1$nn.dist)
# Print row index of the most anomalous point
which.max(nn_mcp1_d)
# Print the 5-nearest neighbor distance score
nn_mcp1_d[1:5]
# Append the score as a new column 
df1$mcp1_knn_score <- nn_mcp1_d
# Visualize knn
plot(nn_mcp1_d ~ 1, data = df1, cex = sqrt(mcp1_knn_score), pch = 20)
ggplot(df1, aes(nn_mcp1_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("MCP-1 KNN Distances")+
  xlab("MCP-1 KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 100
df1 <- within(df1, {   
  mcp1.outlier <- NA
  mcp1.outlier[mcp1_knn_score == 100 | mcp1_knn_score > 100] <- 1
  mcp1.outlier[mcp1_knn_score < 100] <- 0
} )

df1$mcp1.outlier

#IDENTIFY IL-10 OUTLIERS KNN
nn_il10 <- get.knn(df1$X2mo.Obs.Conc.IL.10..27., k=20)
head(nn_il10)
# Create score by averaging distances
nn_il10_d <- rowMeans(nn_il10$nn.dist)
# Print row index of the most anomalous point
which.max(nn_il10_d)
# Print the 5-nearest neighbor distance score
nn_il10_d[1:5]
# Append the score as a new column 
df1$il10_knn_score <- nn_il10_d
# Visualize knn
plot(nn_il10_d ~ 1, data = df1, cex = sqrt(il10_knn_score), pch = 20)
ggplot(df1, aes(nn_il10_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("IL-10 KNN Distances")+
  xlab("IL-10 KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 10
df1 <- within(df1, {   
  il10.outlier <- NA
  il10.outlier[il10_knn_score == 10 | il10_knn_score > 10] <- 1
  il10.outlier[il10_knn_score < 10] <- 0
} )

df1$il10.outlier

#IDENTIFY IL-4 OUTLIERS KNN
nn_il4 <- get.knn(df1$X2mo.Obs.Conc.IL.4..53., k=20)
head(nn_il10)
# Create score by averaging distances
nn_il4_d <- rowMeans(nn_il4$nn.dist)
# Print row index of the most anomalous point
which.max(nn_il4_d)
# Print the 5-nearest neighbor distance score
nn_il4_d[1:5]
# Append the score as a new column 
df1$il4_knn_score <- nn_il4_d
# Visualize knn
plot(nn_il4_d ~ 1, data = df1, cex = sqrt(il4_knn_score), pch = 20)
ggplot(df1, aes(nn_il4_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("IL-4 KNN Distances")+
  xlab("IL-4 KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 50
df1 <- within(df1, {   
  il4.outlier <- NA
  il4.outlier[il4_knn_score == 50 | il4_knn_score > 50] <- 1
  il4.outlier[il4_knn_score < 50] <- 0
} )

df1$il4.outlier

#IDENTIFY IL-8 OUTLIERS KNN
nn_il8 <- get.knn(df1$X2mo.Obs.Conc.IL.8..63., k=20)
head(nn_il8)
# Create score by averaging distances
nn_il8_d <- rowMeans(nn_il8$nn.dist)
# Print row index of the most anomalous point
which.max(nn_il8_d)
# Print the 5-nearest neighbor distance score
nn_il8_d[1:5]
# Append the score as a new column 
df1$il8_knn_score <- nn_il8_d
# Visualize knn
plot(nn_il8_d ~ 1, data = df1, cex = sqrt(il8_knn_score), pch = 20)
ggplot(df1, aes(nn_il8_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("IL-8 KNN Distances")+
  xlab("IL-8 KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 30
df1 <- within(df1, {   
  il8.outlier <- NA
  il8.outlier[il8_knn_score == 30 | il8_knn_score > 30] <- 1
  il8.outlier[il8_knn_score < 30] <- 0
} )

df1$il8.outlier

#IDENTIFY IFN-alpha OUTLIERS KNN
nn_ifnalpha <- get.knn(df1$X2mo.Obs.Conc.IFNalpha2..22., k=20)
head(nn_ifnalpha)
# Create score by averaging distances
nn_ifnalpha_d <- rowMeans(nn_ifnalpha$nn.dist)
# Print row index of the most anomalous point
which.max(nn_ifnalpha_d)
# Print the 5-nearest neighbor distance score
nn_ifnalpha_d[1:5]
# Append the score as a new column 
df1$ifnalpha_knn_score <- nn_ifnalpha_d
# Visualize knn
plot(nn_ifnalpha_d ~ 1, data = df1, cex = sqrt(ifnalpha_knn_score), pch = 20)
ggplot(df1, aes(nn_ifnalpha_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("IFN-alpha KNN Distances")+
  xlab("IFN-alpha KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 100
df1 <- within(df1, {   
  ifnalpha.outlier <- NA
  ifnalpha.outlier[ifnalpha_knn_score == 100 | ifnalpha_knn_score > 100] <- 1
  ifnalpha.outlier[ifnalpha_knn_score < 100] <- 0
} )

df1$ifnalpha.outlier

#Make new variable IFNg:IL-10 ratio
df1 <- within(df1, {
  X2mo.Obs.Conc.ratio.IFNgIL10ratio <- X2mo.Obs.Conc.IFNgamma..25./df1$X2mo.Obs.Conc.IL.10..27.
} )

df1$X2mo.Obs.Conc.ratio.IFNgIL10ratio

df1 <- within(df1, {
  X2mo.Obs.Conc.ratio.IFNgIL10ratio[X2mo.Obs.Conc.ratio.IFNgIL10ratio == Inf | X2mo.Obs.Conc.ratio.IFNgIL10ratio == -Inf] <- 0
} )

df1$X2mo.Obs.Conc.ratio.IFNgIL10ratio

#IDENTIFY IFN-gamma/IL10 ratio OUTLIERS KNN
nn_ifngamma_il10_ratio <- get.knn(df1$X2mo.Obs.Conc.ratio.IFNgIL10ratio, k=20)
head(nn_ifngamma_il10_ratio)
# Create score by averaging distances
nn_ifngamma_il10_ratio_d <- rowMeans(nn_ifngamma_il10_ratio$nn.dist)
# Print row index of the most anomalous point
which.max(nn_ifngamma_il10_ratio_d)
# Print the 5-nearest neighbor distance score
nn_ifngamma_il10_ratio_d[1:5]
# Append the score as a new column 
df1$ifngamma_il10_ratio_knn_score <- nn_ifngamma_il10_ratio_d
# Visualize knn
plot(nn_ifngamma_il10_ratio_d ~ 1, data = df1, cex = sqrt(ifngamma_il10_ratio_knn_score), pch = 20)
ggplot(df1, aes(nn_ifngamma_il10_ratio_d))+
  geom_histogram(binwidth = 0.9) +
  ggtitle("IFN-gamma:IL-10 ratio KNN Distances")+
  xlab("IFN-gamma:IL-10 ratio KNN Distances")+ 
  ylab("Frequency")+
  theme(plot.title = element_text(size=26), axis.text = element_text(size = 18),axis.title = element_text(size = 20), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Create variable for outlier(s), cut off of 15
df1 <- within(df1, {   
  ifngamma_il10_ratio.outlier <- NA
  ifngamma_il10_ratio.outlier[ifngamma_il10_ratio_knn_score == 15 | ifngamma_il10_ratio_knn_score > 15] <- 1
  ifngamma_il10_ratio.outlier[ifngamma_il10_ratio_knn_score < 15] <- 0
} )

df1$ifngamma_il10_ratio.outlier

#REMOVE OUTLIERS###########################################################
df1_ifngamma_sansoutliers <- subset(df1, ifngamma.outlier == 0)
df1_ifngamma_sansoutliers

df1_il6_sansoutliers <- subset(df1, il6.outlier == 0)
df1_il6_sansoutliers 

df1_il1beta_sansoutliers <- subset(df1, il1beta.outlier == 0)
df1_il1beta_sansoutliers

df1_tnfalpha_sansoutliers <- subset(df1, tnfalpha.outlier == 0)
df1_tnfalpha_sansoutliers

df1_mcp1_sansoutliers <- subset(df1, mcp1.outlier == 0)
df1_mcp1_sansoutliers

df1_il10_sansoutliers <- subset(df1, il10.outlier == 0)
df1_il10_sansoutliers

df1_il4_sansoutliers <- subset(df1, il4.outlier == 0)
df1_il4_sansoutliers

df1_il8_sansoutliers <- subset(df1, il8.outlier == 0)
df1_il8_sansoutliers

df1_ifnalpha_sansoutliers <- subset(df1, ifnalpha.outlier == 0)
df1_ifnalpha_sansoutliers

df1_il10_ifngamma_ratio_sansoutliers <- subset(df1, ifngamma_il10_ratio.outlier == 0 )
df1_il10_ifngamma_ratio_sansoutliers


