## GitHub deposits

## set your path to the data
setwd("C:/Users/kls/Desktop/kls Daten/kls research/R/data/Human/HIV_flu_UNC/GitHub")
dir()

## expression matrices are provided as compressed files, need to be unziped - see below

## data tables
######### QC table for HIV RV reads and maps
dd1 <- read.csv2("QC_HIV_110820.csv", header=T)
names(dd1)
head(dd1)
dim(dd1)

########## DEG tables ####################
######### limma HIV-RV table
dd2 <- read.table("Limma_HIV_vs_HIV_RVI_180820.txt", header=T)
names(dd2)
head(dd2)
dim(dd2)
######### limma Tsalik table
dd3 <- read.table("Limma_Tsalik_110820.txt", header=T)
names(dd3)
head(dd3)
dim(dd3)


############### scripts ##############################################

########################### DEGs & Volcano plots ###########################

############# HIV RV data ###########################
setwd("C:/Users/kls/Desktop/kls Daten/kls research/R/data/Human/HIV_flu_UNC/GitHub")
dir()

library(limma);

## data is provided as compressed 7zip file, unzip before using for this script
## normalized data set
dat  <- read.table("norm_HIV_UNC_181119_2.txt", header=TRUE)
names(dat)
dim(dat)
dat[1:5,c(3,16:21)]
## data description
ta <- read.table("sample_descr_github_180820.txt", header=T)
names(ta)
head(ta)
dim(ta)
table(ta[,c(2,3)])
identical(names(dat)[8:69], as.character(ta[,1]))

## select two groups
ta$status
sel1 <- as.character(c("HIV" ,"HIV_vir")); sel1
ta1 <- ta[ta$status %in% sel1,]
dim(ta); dim(ta1)
ta1 <- droplevels(ta1)
ta1$status
sel2 <- ta1$sample_ID; sel2

test1 <- dat[,colnames(dat) %in% sel2]; names(test1)
identical(colnames(test1),as.character(ta1$sample_ID))
n <- nrow(test1); n

#######  Differential Expression for RV infected versus non-infected HIV patients
names(ta)
ta1$status
xx1 <- ta1$status; xx1
xx2 <- gsub("HIV_vir", 1,xx1); xx2
xx3 <- as.numeric(gsub("HIV", 0,xx2)); xx3
infct <- xx3
xx4 <- gsub("HIV_vir", 0,xx1); xx4
xx5 <- as.numeric(gsub("HIV", 1,xx4)); xx5
non_inf <- xx5
design <-as.data.frame(cbind(infct,non_inf)); design

fit<-lmFit(test1,design);
contrast.matrix<- makeContrasts(infct-non_inf, levels=design);
fit2 <- contrasts.fit(fit, contrast.matrix);
options(warn=-1)# to suppress the warning message if any;
fit2 <- eBayes(fit2);
y1 <- topTable (fit2, coef=1, number=n, adjust="BH");
y1 [1:10,]; dim(y1)
names(y1)

## retrieve data from dat table that match the DEG results and merge results
names(dat)
dat[1:5,c(3,10:12)]
y2 <- merge(dat, y1, by.x = "row.names", by.y = "row.names", all = FALSE)
y2[1:3,]
dim(y2)
names(y2)

############  selection of DEGs with certain logFC and p-value:
lim.FC <- log2(1.5); lim.FC
lim.P <- 0.05; lim.P
#### select significance and logFC adjusted p-value: no hits
# y3 <- subset(y2,abs(y2$logFC)>lim.FC & y2$adj.P.Val<lim.P)
## use non-adjusted p-value
y3 <- subset(y2,abs(y2$logFC)>lim.FC & y2$P.Value<lim.P)
names(y3); dim(y3)

y5 <- y3[order(y3$logFC, decreasing=T),]; 
names(y5);  #sort by logFC
y5[1:20,c(2,4,71,74)]
dim(y5)

##### list with unique genes (top probeset only)
dup1 <- duplicated(y5[,5]); dup1
y66 <- y5[dup1,c(2,4,71,74)]; head(y66) ##some genes are duplicates!!
y6 <- y5[!dup1,]; head(y6)
names(y6)
dim(y5);dim(y6);
y6[1:20,c(2,4,71,74)]

#calculate up and down-regulated
up.reg <- which(y5$logFC>0) ## genes with lgfc > 1
up1 <- y5[up.reg,]
down.reg <- which(y5$logFC<0) ## genes with lgfc < 1
down1 <- y5[down.reg,]
dim(y5)
dim(up1)
dim(down1)

# Make a basic volcano plot
lim.FC; lim.P
y2[1:3,]; dim(y2);
names(y2)
par(mfrow=c(1,1), font=2,font.axis=2,font.lab=3, mar=c(4,4,2,2)); #bottom,left, top, right ; xpd=TRUE
with(y2, plot(logFC, -log10(P.Value), pch=20, 
              main=paste("infected", "vs heathy controls"))) #, xlim=c(-2.5,2)

# Add colored points:
with(subset(y2, abs(logFC)>lim.FC & P.Value<lim.P), 
     points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(y2, adj.P.Val<0.001 & abs(logFC)>2), 
     points(logFC, -log10(adj.P.Val), pch=20, col="green"))



########## Tsalik data #############################
setwd("C:/Users/kls/Desktop/kls Daten/kls research/R/data/Human/HIV_flu_UNC/GitHub")
dir()
library(limma)

## data is provided as compressed 7zip file, unzip before using for this script
################ log2 transformed, normalized data, controls and viral infections only
dat1  <- read.table("sbst1_GSE63990_Tsalik_2016_140820.txt", header=TRUE)
names(dat1)
dim(dat1)
dat1[1:5,c(3,16:19)]

#select colums with the data
names(dat1)
dat2 <- as.matrix(dat1[,6:ncol(dat1)])
colnames(dat2)
dim(dat2)
head(dat2)

########## Tsalik data description ######################
ta2 <- read.table("sbst1_target_GSE63990_Tsalik_2016_220118.txt", header=TRUE);
names(ta2)
head(ta2)
table(ta2[,c(3)])
identical(colnames(dat2), as.character(ta2[,2]))

test1 <- dat2; colnames(test1)
n <- nrow(test1); n

#######  Differential Expression for infected versus healthy
names(ta2)
xx1 <- ta2$status; xx1
xx2 <- gsub("viral", 1,xx1); xx2
xx3 <- as.numeric(gsub("non-inf_ill", 0,xx2)); xx3
viral <- xx3
xx4 <- gsub("viral", 0,xx1); xx4
xx5 <- as.numeric(gsub("non-inf_ill", 1,xx4)); xx5
non_inf_ill <- xx5
design <-as.data.frame(cbind(viral,non_inf_ill)); design
## test if ok
test <- as.data.frame(cbind(design,as.character(ta2$status))); test
library(limma);
fit<-lmFit(test1,design);
contrast.matrix<- makeContrasts(viral-non_inf_ill, levels=design);
fit2 <- contrasts.fit(fit, contrast.matrix);
options(warn=-1)# to suppress the warning message if any;
fit2 <- eBayes(fit2);
y1 <- topTable (fit2, coef=1, number=n, adjust="BH");
y1 [1:10,]; dim(y1)
names(y1)

## retrieve data from dat table that match the DEG results and merge results
names(dat1)
dat1[1:5,c(1,3,10:12)]
y2 <- merge(dat1, y1, by.x = "row.names", by.y = "row.names", all = FALSE)
y2[1:3,]
dim(y2)
names(y2)

############  selection of DEGs with certain logFC and p-value:
lim.FC <- log2(1.5); lim.FC
lim.P <- 0.05; lim.P
reslog<-which(abs(y2$logFC)>lim.FC); reslog #log2(1.1)); reslog   ## log2(1.2) = 0.2630344 ## log2FoldChange: >1 means 2-fold regulated
resp<-which(y2$adj.P.Val<lim.P); resp    #### adj.p is FDR with BH,
RowID<-intersect(reslog,resp)####intersektion, um die richtige RowID zu finden

y3 <- y2[RowID,] ## view data
names(y3); dim(y3)

y5 <- y3[order(y3$logFC, decreasing=T),]; 
names(y5);  #sort by logFC
y5[1:20,c(4,210,214)]
dim(y5)

##### list with unique genes (top probeset only)
dup1 <- duplicated(y5[,3]); dup1
y66 <- y5[dup1,c(4,210,214)]; head(y66) ##some genes are duplicates!!
y6 <- y5[!dup1,]; head(y6)
names(y6)
dim(y5);dim(y6);
y6[1:20,c(4,210,214)]

#calculate up and down-regulated
up.reg <- which(y5$logFC>0) ## genes with lgfc > 1
up1 <- y5[up.reg,]
down.reg <- which(y5$logFC<0) ## genes with lgfc < 1
down1 <- y5[down.reg,]
dim(y5)
dim(up1)
dim(down1)

# Make a basic volcano plot
lim.FC; lim.P
y2[1:3,]; dim(y2);
names(y2)
par(mfrow=c(1,1), font=2,font.axis=2,font.lab=3, mar=c(4,4,2,2)); #bottom,left, top, right ; xpd=TRUE
with(y2, plot(logFC, -log10(adj.P.Val), pch=20, 
              main=paste("infected", "vs heathy controls"))) #, xlim=c(-2.5,2)

# Add colored points:
with(subset(y2, abs(logFC)>lim.FC & adj.P.Val<lim.P), 
     points(logFC, -log10(adj.P.Val), pch=20, col="orange"))
with(subset(y2, adj.P.Val<0.001 & abs(logFC)>2), 
     points(logFC, -log10(adj.P.Val), pch=20, col="green"))


