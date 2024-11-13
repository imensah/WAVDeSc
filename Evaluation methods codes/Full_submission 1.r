#Load Libraries

library(ggplot2)
library(reshape2)
library(matrixStats)
library(vioplot)
library(plotrix)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(aplpack)
library(coop)
library(Rfast)



############# Simulation using Symsim ################################








########################## Visualization of a signal #######################################################
Real1<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Noisy_1000_UMI.csv")
View(Real1)
Real1 = Real1[1,]
Real1 = Real1[,-1]
Real1 = t(Real1)
Real1 = as.matrix(Real1)
plot(Real, type = "l")
plot(Real, type ='l', ylab = "Gene Expression Values", xlab = "Cell", lwd = 4)
plot(Real1, type ='l', ylab = "Gene Expression Values", xlab = "Cells", cex.lab = 1.5, lwd = 2, font.axis = 2, cex.axis = 1.7,
     main = "Signal Representationof Gene 1 across the cells")




#Plot the gene signal and its density plot for 500genes by 300cells 1
reads <- read.delim("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Real_Data/Tung/reads.csv")
A = reads[1,]
A=A[,-1]
A=t(A)
colnames(A) = 'G'
A=as.data.frame(A)
signal <-A$G


# Set the margins for individual plots within the mfrow layout
#par(mfrow = c(2, 1), mar = c(4, 4, 4, 4))
t = plot(Real1, type ='l', ylab = "Gene Expression Values", xlab = "Cell")
corner.label("A",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)
d <- density(signal) # returns the density A
plot(d, main = NA)
corner.label("B",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)







##############################Plotting the distributions of the various datasets#########################
#Violin plot for 500 by 300

#Truth data
D1 <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Truth_umi.csv")
d1=reshape::melt(D1,id.vars = "Genes")
Truth=d1$value

#UMI data
D2<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Noisy_umi.csv")
d2 <- reshape::melt(D2,id.vars = "Genes")
UMI=d2$value

#nonUMI data
D3 = read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Noisy_nonUMI.csv")
d3 <- reshape::melt(D3,id.vars = "Genes")
nonUMI=d3$value

#Plot using the violin to check distributions
par(mfrow = c(3, 1))
par(cex.lab = 1.35)
vioplot(UMI,nonUMI, names = c("UMI","nonUMI"),
        col="gold",xlab = 'Datasets', ylab = 'Gene Expression levels', cex.axis=1.4)
grid()
corner.label("A",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)


#Violin plot for 5000 by 1000
D4<-read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Truth_1000.csv")
d4=reshape::melt(D4,id.vars = "Genes")
Truth=d4$value

D5<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Noisy_1000_UMI.csv")
d5=reshape::melt(D5,id.vars = "Genes")
UMI=d5$value

D6<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Noisy_1000_nonUMI.csv")
d6=reshape::melt(D6,id.vars = "Genes")
nonUMI=d6$value

#Plot using the violin to check distributions
par(cex.lab = 1.35)
vioplot(UMI,nonUMI, names = c("UMI","nonUMI"),
        col="#00BFC4",xlab = 'Datasets', ylab = "Gene Expression levels", cex.axis=1.7, main = "Distribution of Data with Dimension 5000 genes by 1000 cells")
grid()
corner.label("A",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)



#Violin plot for Tung data
D7<- read.delim("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Real_Data/Tung/reads.csv")
d7=reshape::melt(D7,id.vars = "Transcript")
Tung=d7$value
#summary(Tung)
#sd(Tung)
par(cex.lab = 1.35)
vioplot(Tung,  names= c("Pluripotent Stem Cells Data"), col="#7CAE00",xlab = 'Datasets', ylab = "Gene Expression levels", cex.axis=1.7,main = "Distribution of Data on PSC")
corner.label("B",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)


D8 = ("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Real_Data/Autism/Combined_data.csv")
d8=reshape::melt(D8,id.vars = "gene_ID")
Autism = d8$value
summary(Autism)
sd(Autism)
par(cex.lab = 1.35)
vioplot(AUtism, names= c("Autism Data"),col="skyblue",xlab = 'Datasets', ylab = NA, cex.axis=1.7)
corner.label("E",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)





############################# Variability Analysis (Mean-Variance Plot) ######################################
par(mfrow = c(3, 2), mar = c(4, 4, 4, 4))
## 500 by 300
# Calculate the mean expression and variance for each gene
D2<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Noisy_umi.csv")
mean_expression <- apply(D2[,2:301], 1, mean)
variance <- apply(D2[,2:301], 1, var)

# Create the mean-variance plot
plot(log10(variance) ~ log10(mean_expression), pch = 16, cex = 0.5, xlab = " log10(mean_expression)", ylab = "log10(variance)")#, main = "Mean-Variance Plot")

grid()

# Add a trend line
#abline(lm(log10(variance) ~ log10(mean_expression)), col = "red")

#lines(lowess(log10(mean_expression), log10(variance)), col = "red")

# Add circles around lower, middle, and upper regions
lower_genes <- log10(variance) < quantile(log10(variance), 0.25)
middle_genes <- log10(variance) >= quantile(log10(variance), 0.25) & log10(variance) <= quantile(log10(variance), 0.66)
upper_genes <- log10(variance) > quantile(log10(variance), 0.50)


points( log10(mean_expression)[lower_genes], log10(variance)[lower_genes],  col = "deepskyblue3",pch = 1)
points( log10(mean_expression)[middle_genes], log10(variance)[middle_genes],  pch = 16)
points( log10(mean_expression)[upper_genes], log10(variance)[upper_genes], col = "azure4", pch = 9)
corner.label("A",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)


D3 = read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Noisy_nonUMI.csv")
mean_expression <- apply(D3[,2:301], 1, mean)
variance <- apply(D3[,2:301], 1, var)

# Create the mean-variance plot
plot(log10(variance) ~ log10(mean_expression), pch = 16, cex = 0.5, xlab = " log10(mean_expression)", ylab = "log10(variance)")#, main = "Mean-Variance Plot")

grid()

# Add circles around lower, middle, and upper regions
lower_genes <- log10(variance) < quantile(log10(variance), 0.25)
middle_genes <- log10(variance) >= quantile(log10(variance), 0.25) & log10(variance) <= quantile(log10(variance), 0.66)
upper_genes <- log10(variance) > quantile(log10(variance), 0.50)


points( log10(mean_expression)[lower_genes], log10(variance)[lower_genes], col = "deepskyblue3",pch = 1)
points( log10(mean_expression)[middle_genes], log10(variance)[middle_genes],  pch = 16)
points( log10(mean_expression)[upper_genes], log10(variance)[upper_genes], col = "azure4", pch = 9)
corner.label("B",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)



##5000 by 1000
D5<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Noisy_1000_UMI.csv")
mean_expression <- apply(D5[,2:1001], 1, mean)
variance <- apply(D5[,2:1001], 1, var)

# Create the mean-variance plot
plot(log10(variance) ~ log10(mean_expression), pch = 16, cex = 0.5, xlab = " log10(mean_expression)", ylab = "log10(variance)")#, main = "Mean-Variance Plot")

grid()

# Add circles around lower, middle, and upper regions
lower_genes <- log10(variance) < quantile(log10(variance), 0.25)
middle_genes <- log10(variance) >= quantile(log10(variance), 0.25) & log10(variance) <= quantile(log10(variance), 0.66)
upper_genes <- log10(variance) > quantile(log10(variance), 0.50)


points( log10(mean_expression)[lower_genes], log10(variance)[lower_genes], col = "deepskyblue3",pch = 1)
points( log10(mean_expression)[middle_genes], log10(variance)[middle_genes],  pch = 16)
points( log10(mean_expression)[upper_genes], log10(variance)[upper_genes],col = "azure4", pch = 9)
corner.label("C",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)


D6<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Noisy_1000_nonUMI.csv")
mean_expression <- apply(D6[,2:1001], 1, mean)
variance <- apply(D6[,2:1001], 1, var)
# Create the mean-variance plot
plot(log10(variance) ~ log10(mean_expression), pch = 16, cex = 0.5, xlab = " log10(mean_expression)", ylab = "log10(variance)", main = "Mean-Variance Plot for D6")
grid()

# Add circles around lower, middle, and upper regions
lower_genes <- log10(variance) < quantile(log10(variance), 0.25)
middle_genes <- log10(variance) >= quantile(log10(variance), 0.25) & log10(variance) <= quantile(log10(variance), 0.66)
upper_genes <- log10(variance) > quantile(log10(variance), 0.50)
points( log10(mean_expression)[lower_genes], log10(variance)[lower_genes], col =  "#F8766D", pch = 20)
points( log10(mean_expression)[middle_genes], log10(variance)[middle_genes],  pch = 20)
points( log10(mean_expression)[upper_genes], log10(variance)[upper_genes], col = "#00BFC4", pch = 20)
corner.label("A ",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)


D7<- read.delim("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Real_Data/Tung/reads.csv")
mean_expression <- apply(D7[,2:865], 1, mean)
variance <- apply(D7[,2:865], 1, var)

# Create the mean-variance plot
plot(log10(variance) ~ log10(mean_expression), pch = 16, cex = 0.5, xlab = " log10(mean_expression)", ylab = "log10(variance)", main = "Mean-Variance Plot for PSC Data")
grid()

 # Add circles around lower, middle, and upper regions
lower_genes <- log10(variance) < quantile(log10(variance), 0.25)
middle_genes <- log10(variance) >= quantile(log10(variance), 0.25) & log10(variance) <= quantile(log10(variance), 0.66)
upper_genes <- log10(variance) > quantile(log10(variance), 0.50)

points( log10(mean_expression)[lower_genes], log10(variance)[lower_genes], col = "#F8766D",pch = 16)
points( log10(mean_expression)[middle_genes], log10(variance)[middle_genes],  pch = 20)
points( log10(mean_expression)[upper_genes], log10(variance)[upper_genes], col = "#00BFC4", pch = 20)
corner.label("B",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)

# #D8 = read.csv("~/Downloads/Data_submission/Real_Data/Autism/Combined_data.csv")
# mean_expression <- apply(D8[,2:865], 1, mean)
# variance <- apply(D8[,2:865], 1, var)
# 
# # Create the mean-variance plot
# plot(log10(variance) ~ log10(mean_expression), pch = 16, cex = 0.5, xlab = " log10(mean_expression)", ylab = "log10(variance)")#, main = "Mean-Variance Plot")
# grid()
# 
# # Add circles around lower, middle, and upper regions
# lower_genes <- log10(variance) < quantile(log10(variance), 0.25)
# middle_genes <- log10(variance) >= quantile(log10(variance), 0.25) & log10(variance) <= quantile(log10(variance), 0.66)
# upper_genes <- log10(variance) > quantile(log10(variance), 0.50)
# 
# points( log10(mean_expression)[lower_genes], log10(variance)[lower_genes], col = "deepskyblue3",pch = 16)
# points( log10(mean_expression)[middle_genes], log10(variance)[middle_genes],  pch = 16)
# points( log10(mean_expression)[upper_genes], log10(variance)[upper_genes], col = "azure4", pch = 16)
# corner.label("F",y=1,x=-1,figcorner=TRUE, font = 2, cex = 1.8)








#################################### cell and gene filtering #################################################
#filtering genes
D3 <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Truth_umi.csv")
row_cv <- apply(D3[,2:301], 1, function(x) sd(x) / mean(x))
sum(row_cv==0)
filtering$CV=as.numeric(filtering$CV)
filtering$SNR=as.numeric(filtering$SNR)
filtering_na=na.omit(filtering)
View(filtering_na)


##filtering cells
cells_filter<- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Noisy_umi.csv")
View(cells_filter)
cells_filter$CV=as.numeric(cells_filter$CV)
cells_filter$SNR=as.numeric(cells_filter$SNR)
cells_filter_na=na.omit(cells_filter)
View(cells_filter_na)

##filtering_nonUMI
filternonUMI <- read.csv("~/Desktop/Final_Glory/Final_simulation/Data/nonUMI/nonUMI_Noisy_300.csv")
View(filternonUMI)
filternonUMI$CV=as.numeric(filternonUMI$CV)
filternonUMI$SNR=as.numeric(filternonUMI$SNR)
filternonUMI_na=na.omit(filternonUMI)
dim(filternonUMI_na)

nonUMI_transpose <- read.csv("~/Desktop/Final_Glory/Final_simulation/Data/nonUMI/nonUMI_transpose.csv")
View(nonUMI_transpose)
nonUMI_transpose$cv=as.numeric(nonUMI_transpose$cv)
nonUMI_transpose$SNR=as.numeric(nonUMI_transpose$SNR)
nonUMI_transpose_na=na.omit(nonUMI_transpose)
dim(nonUMI_transpose_na)





D1 <- matrix(rnorm(1000), nrow = 100, ncol = 10)

# Calculate the coefficient of variation for each row
row_cv <- apply(D1, 1, function(x) sd(x) / mean(x))

# Print the coefficient of variation for each row
print(row_cv)



################## Check for unique cell populations for PSC data################################
D7<- read.delim("~/Downloads/Data_submission/Real_Data/Tung/reads.csv")
A = colnames(D7) #extract column names
extracted_string <- gsub("^(NA.{0,5}).*", "\\1", A) #Extract the string up to "NA........9"
unique(extracted_string) #Sort unique cells









##################################### Results for Objective 2############################################
#### Comparison of truth and noisy data to denoised datasets using the CV


#### Simulated Data 5000 genes by 1000 cells
Truth_500 <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Truth_umi.csv")
Truth_500 = Truth_500[,-1]
CV_Truth_500 = rowcvs(as.matrix(Truth_500), ln = FALSE, unbiased = TRUE) 
View(CV_Truth_500)

Noisy_500_UMI <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Noisy_umi.csv")
Noisy_500_UMI = Noisy_500_UMI[,-1]
CV_Noisy_500_UMI = rowcvs(as.matrix(Noisy_500_UMI), ln = FALSE, unbiased = TRUE) 
#CV_Noisy = Truth$CV_R
#plot(CV_Noisy_1000_UMI,type = 'l')

Noisy_500_nonUMI <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simlated_500by300/Noisy_nonUMI.csv")
Noisy_500_nonUMI = Noisy_500_nonUMI[,-1]
CV_Noisy_500_nonUMI = rowcvs(as.matrix(Noisy_500_nonUMI), ln = FALSE, unbiased = TRUE) 
#CV_Noisy = Truth$CV_R
#plot(CV_Noisy_1000_UMI,type = 'l')


Gene = paste(1:500, sep = " ")

df = cbind(Gene,CV_Truth_500,CV_Noisy_500_nonUMI)
View(df)
df5= as.data.frame(df)
View(df5)

#melting Data for ggplot
df_dat <- melt(df5,id.vars = "Gene")
View(df_dat)

n1 = t.test(CV_Truth_500,CV_Noisy_500_nonUMI)

















Gene = paste(1:5000, sep = " ")

df = cbind(Gene,CV_Truth_5000,CV_Noisy_5000_UMI)
df5= as.data.frame(df)
View(df5)

#melting Data for ggplot
df_dat <- melt(df5,id.vars = "Gene")
#View(df_dat)

class(df_dat$Gene)
class(df_dat$variable)
df_dat$Gene=as.integer(df_dat$Gene)
df_dat$value = as.numeric(df_dat$value)
class(df_dat$value)




A=ggplot(df_dat, aes(x = Gene, y = value))+
  geom_line(aes(color = variable), show.legend = FALSE)+
  geom_point(data = df_dat[df_dat$Gene == max(df_dat$Gene),],
             aes(x =Gene, y = value, color = variable), show.legend = TRUE)+ theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20))+theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+theme(legend.title = element_blank())+ theme(legend.text = element_text(size=15))+
  theme(legend.position = c(.92, .9))+theme(legend.key.size = unit(1, 'cm'),legend.key.width= unit(1, 'cm'))
A+ ylab("Coefficient of Variation")
















UMI_magic500 <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/Magic/denoised_UMI_magic1000.csv")
UMI_magic1000 = UMI_magic1000[,-1]
UMI_magic1000=t(UMI_magic1000)
UMI_magic1000 = as.data.frame(UMI_magic1000)
CV_UMI_magic1000 = rowcvs(as.matrix(UMI_magic1000), ln = FALSE, unbiased = TRUE) 
#plot(CV_UMI_magic1000,type = 'l')



saver_UMI <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/saver/denoised_saver_UMI.csv")
CV_saver_1000_UMI = rowcvs(as.matrix(saver_UMI), ln = FALSE, unbiased = TRUE) 
#plot(CV_saver_1000_UMI,type = 'l')


db61000_UMI <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/db6/denoised_db61000_UMI.csv", header=FALSE)
CV_db61000_UMI  = rowcvs(as.matrix(db61000_UMI ), ln = FALSE, unbiased = TRUE) 
#plot(CV_db61000_UMI,type = 'l')


bior1000_UMI <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/bior2.6/denoised_bior1000_UMI.csv", header=FALSE)
CV_bior1000_UMI = rowcvs(as.matrix(bior1000_UMI), ln = FALSE, unbiased = TRUE) 
#plot(CV_bior1000_UMI,type = 'l')

enhance1000_UMI <- read.delim("~/Downloads/5000 1000/enhance/enhance_UMI/denoised_enhance1000_UMI.tsv")
enhance1000_UMI=enhance1000_UMI[,-1]
CV_enhance1000_UMI = rowcvs(as.matrix(enhance1000_UMI), ln = FALSE, unbiased = TRUE) 
#plot(CV_enhance1000_UMI,type = 'l')





#plotting
Gene = paste(1:500, sep = " ")

df = cbind(Gene,CV_Truth_500,CV_Noisy_500_UMI)
df3 = cbind(Genes,CV_Truth_1000_UMI,CV_enhance1000_UMI)
df5= as.data.frame(df)
#View(df5)


#melting Data for ggplot
df_dat <- melt(df5,id.vars = "Gene")
#View(df_dat)

class(df_dat$Gene)
class(df_dat$variable)
df_dat$Gene=as.integer(df_dat$Gene)
df_dat$value = as.numeric(df_dat$value)
class(df_dat$value)



A=ggplot(df_dat, aes(x = Gene, y = value))+
  geom_line(aes(color = variable), show.legend = FALSE)+
  geom_point(data = df_dat[df_dat$Gene == max(df_dat$Gene),],
             aes(x =Gene, y = value, color = variable), show.legend = TRUE)+ theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20))+theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+theme(legend.title = element_blank())+ theme(legend.text = element_text(size=15))+
  theme(legend.position = c(.92, .9))+theme(legend.key.size = unit(1, 'cm'),legend.key.width= unit(1, 'cm'))
A+ ylab("Coefficient of Variation")

B=ggplot(df_dat, aes(x = Genes, y = value))+
  geom_line(aes(color = variable), show.legend = FALSE)+
  geom_point(data = df_dat[df_dat$Genes == max(df_dat$Genes),],
             aes(x =Genes, y = value, color = variable), show.legend = TRUE)+ theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20))+theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+theme(legend.title = element_blank())+ theme(legend.text = element_text(size=15))+
  theme(legend.position = c(.92, .9))+theme(legend.key.size = unit(1, 'cm'),legend.key.width= unit(1, 'cm'))
A+ ylab("Coefficient of Variation")


### Finding significance between CV of models with noisy data

n1 = t.test(CV_Truth_1000_UMI,CV_Noisy_1000_UMI)

n2 = t.test(CV_Noisy_1000_UMI,CV_saver_1000_UMI)

n3 = t.test(CV_Noisy_1000_UMI,CV_enhance1000_UMI)

n4 = t.test(CV_Noisy_1000_UMI,CV_UMI_magic1000)

n5 = t.test(CV_Noisy_1000_UMI,CV_db61000_UMI)

n6 = t.test(CV_Noisy_1000_UMI,CV_bior1000_UMI)


### Finding significance between CV of models with noisy data

t1 = t.test(CV_Truth,CV_Noisy)

t2 = t.test(CV_Truth,CV_saver)

t3 = t.test(CV_Truth,CV_enhance)

t4 = t.test(CV_Truth,CV_Magic)

t5 = t.test(CV_Truth,CV_db6)

t6 = t.test(CV_Truth,CV_bior2.6)









#### Simulated Data 5000 by 1000
Truth_1000 <-read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Truth_1000.csv")
Truth_1000 = Truth_1000[,-1]
Truth_1000 =t(Truth_1000)

CV_Truth_1000 = rowcvs(as.matrix(Truth_1000), ln = FALSE, unbiased = TRUE) 


Noisy_5000_UMI <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Data_submission/Simulated_5000by1000/Noisy_1000_UMI.csv")
Noisy_1000_UMI = Noisy_5000_UMI[,-1]
Noisy_1000_UMI = t(Noisy_1000_UMI)

CV_Noisy_1000_UMI = rowcvs(as.matrix(Noisy_1000_UMI), ln = FALSE, unbiased = TRUE) 
#CV_Noisy = Truth$CV_R
#plot(CV_Noisy_1000_UMI,type = 'l')


Gene = paste(1:1000, sep = " ")
df11 = cbind(Gene,CV_Truth_5000,CV_Noisy_5000_UMI)


class(df_dat$Gene)
class(df_dat$variable)
df_dat$Gene=as.integer(df_dat$Gene)
df_dat$value = as.numeric(df_dat$value)
class(df_dat$value)



A=ggplot(df_dat, aes(x = Gene, y = value))+
  geom_line(aes(color = variable), show.legend = FALSE)+
  geom_point(data = df_dat[df_dat$Gene == max(df_dat$Gene),],
             aes(x =Gene, y = value, color = variable), show.legend = TRUE)+ theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20))+theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+theme(legend.title = element_blank())+ theme(legend.text = element_text(size=15))+
  theme(legend.position = c(.92, .9))+theme(legend.key.size = unit(1, 'cm'),legend.key.width= unit(1, 'cm'))
A+ ylab("Coefficient of Variation")
















UMI_magic500 <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/Magic/denoised_UMI_magic1000.csv")
UMI_magic1000 = UMI_magic1000[,-1]
UMI_magic1000=t(UMI_magic1000)
UMI_magic1000 = as.data.frame(UMI_magic1000)
CV_UMI_magic1000 = rowcvs(as.matrix(UMI_magic1000), ln = FALSE, unbiased = TRUE) 
#plot(CV_UMI_magic1000,type = 'l')



saver_UMI <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/saver/denoised_saver_UMI.csv")
CV_saver_1000_UMI = rowcvs(as.matrix(saver_UMI), ln = FALSE, unbiased = TRUE) 
#plot(CV_saver_1000_UMI,type = 'l')


db61000_UMI <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/db6/denoised_db61000_UMI.csv", header=FALSE)
CV_db61000_UMI  = rowcvs(as.matrix(db61000_UMI ), ln = FALSE, unbiased = TRUE) 
#plot(CV_db61000_UMI,type = 'l')


bior1000_UMI <- read.csv("~/Downloads/5000 1000/5000*1000/UMI/bior2.6/denoised_bior1000_UMI.csv", header=FALSE)
CV_bior1000_UMI = rowcvs(as.matrix(bior1000_UMI), ln = FALSE, unbiased = TRUE) 
#plot(CV_bior1000_UMI,type = 'l')

enhance1000_UMI <- read.delim("~/Downloads/5000 1000/enhance/enhance_UMI/denoised_enhance1000_UMI.tsv")
enhance1000_UMI=enhance1000_UMI[,-1]
CV_enhance1000_UMI = rowcvs(as.matrix(enhance1000_UMI), ln = FALSE, unbiased = TRUE) 
#plot(CV_enhance1000_UMI,type = 'l')





#plotting
Gene = paste(1:500, sep = " ")

df = cbind(Gene,CV_Truth_500,CV_Noisy_500_UMI)
df3 = cbind(Genes,CV_Truth_1000_UMI,CV_enhance1000_UMI)
df5= as.data.frame(df)
#View(df5)


#melting Data for ggplot
df_dat <- melt(df5,id.vars = "Gene")
#View(df_dat)

class(df_dat$Gene)
class(df_dat$variable)
df_dat$Gene=as.integer(df_dat$Gene)
df_dat$value = as.numeric(df_dat$value)
class(df_dat$value)



A=ggplot(df_dat, aes(x = Gene, y = value))+
  geom_line(aes(color = variable), show.legend = FALSE)+
  geom_point(data = df_dat[df_dat$Gene == max(df_dat$Gene),],
             aes(x =Gene, y = value, color = variable), show.legend = TRUE)+ theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20))+theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+theme(legend.title = element_blank())+ theme(legend.text = element_text(size=15))+
  theme(legend.position = c(.92, .9))+theme(legend.key.size = unit(1, 'cm'),legend.key.width= unit(1, 'cm'))
A+ ylab("Coefficient of Variation")

B=ggplot(df_dat, aes(x = Genes, y = value))+
  geom_line(aes(color = variable), show.legend = FALSE)+
  geom_point(data = df_dat[df_dat$Genes == max(df_dat$Genes),],
             aes(x =Genes, y = value, color = variable), show.legend = TRUE)+ theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20))+theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+theme(legend.title = element_blank())+ theme(legend.text = element_text(size=15))+
  theme(legend.position = c(.92, .9))+theme(legend.key.size = unit(1, 'cm'),legend.key.width= unit(1, 'cm'))
A+ ylab("Coefficient of Variation")


### Finding significance between CV of models with noisy data

n1 = t.test(CV_Truth_1000_UMI,CV_Noisy_1000_UMI)

n2 = t.test(CV_Noisy_1000_UMI,CV_saver_1000_UMI)

n3 = t.test(CV_Noisy_1000_UMI,CV_enhance1000_UMI)

n4 = t.test(CV_Noisy_1000_UMI,CV_UMI_magic1000)

n5 = t.test(CV_Noisy_1000_UMI,CV_db61000_UMI)

n6 = t.test(CV_Noisy_1000_UMI,CV_bior1000_UMI)


### Finding significance between CV of models with noisy data

t1 = t.test(CV_Truth,CV_Noisy)

t2 = t.test(CV_Truth,CV_saver)

t3 = t.test(CV_Truth,CV_enhance)

t4 = t.test(CV_Truth,CV_Magic)

t5 = t.test(CV_Truth,CV_db6)

t6 = t.test(CV_Truth,CV_bior2.6)

























############# DEGs Identification on ASD data ####
Methods = c("Groundtruth","Raw Data","MAGIC","ENHANCE", "SAVER", "WAVDeSC")
Significant_Genes = c( 1882, 995, 3214, 3277, 6359, 2935)
SigGe = data.frame(Methods , Significant_Genes)
View(SigGe)

# Barplot
SigGe$Methods <- factor(SigGe$Methods, levels = c("Groundtruth", "Raw Data", "MAGIC", "ENHANCE", "SAVER", "WAVDeSC"))

ggplot(SigGe, aes(x = Methods, y = Significant_Genes, fill = Significant_Genes)) + 
  geom_bar(stat = "identity", width = 0.4) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Define the color gradient
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20, angle = 45, hjust = 1))+
  # theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))+
  theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+
  theme(axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))

# Barplot
SigGe$Methods <- factor(SigGe$Methods, levels = c("Groundtruth","Raw Data","MAGIC","ENHANCE", "SAVER"))#, "WAVDeSC"))

ggplot(SigGe, aes(x=Methods, y=Significant_Genes)) + 
  geom_bar(stat = "identity",width=0.4)+
  scale_fill_grey(start = 0.25, end = 0.75)+
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20, angle = 45, hjust = 1))+
 # theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))+
  theme(axis.title.y = element_text(size=25))+
  theme(axis.title.x = element_text(size=25))+
  theme(axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))









