# R script for 10-fold CV (randomForest and XGBoost)
# setwd("C:/Users/Abraham Moller/Documents/ReadLab")
# keep all data in the same directory 
# stuff for phylogenetic signal and CC association
require(picante)
require(ape)
require(adephylo)
require(ade4)
require(phylobase)
require(geiger)
require(phytools)
# stuff for random forest models
require(randomForest)
# stuff for extreme gradient boosted decision trees
require(data.table)
require(xgboost)
require(dplyr)
require(caret)
# stuff to write xlsx
require(xlsx)
phagetree <- read.tree(file="core_gene_alignment.filtered.iqtree")
phagetree$node.label <- NULL
# remove NRS268 - bad strain and no phage resistance data
phagetree <- drop.tip(phagetree,"NRS268")
# get all NARSA/CF/VISA phage resistance data
phagedata <- read.delim(file="phagedata.txt", header=TRUE, sep="\t", strip.white=TRUE)
# get qualitative (ternary) phage resistance
phagedata_qual <- read.delim(file="phagedata_qual.txt", header=TRUE, sep="\t", strip.white=TRUE)
# remove strain missing in phylogeny but present in the data
phagedata <- phagedata[ !grepl("NRS408", phagedata$Strain) , ]
phagedata_qual <- phagedata_qual[ !grepl("NRS408", phagedata_qual$Strain) , ]
# set row names of phage data to strain names
row.names(phagedata) <- phagedata[,1]
phagedata <- phagedata[,-1]
row.names(phagedata_qual) <- phagedata_qual[,1]
phagedata_qual <- phagedata_qual[,-1]
phagedata_qual <- phagedata_qual[ !grepl("NRS268", row.names(phagedata_qual)) , ]
# remove strains with no match between tree and phenotype datasets
phagedata_qual <- phagedata_qual[ !grepl("NA.*", row.names(phagedata_qual)), ]
# sort phage data row names based on 259 labels for phagetree tips
phagetree <- as(phagetree, "phylo4")
# phagetreelabels <- data.frame(labels=labels(phagetree)[1:124])
phagetreelabels <- data.frame(labels=labels(phagetree)[1:261])
# wasn't able to sort by phage tree labels in R
orderedphagedata <- phagedata[match(row.names(phagetreelabels),row.names(phagedata)),]
orderedphagedata <- orderedphagedata[ !grepl("NRS268", row.names(orderedphagedata)) , ]
# remove strains with no match between tree and phenotype datasets
orderedphagedata <- orderedphagedata[ !grepl("NA.*", row.names(orderedphagedata)), ]
# make random forest models
# first load feature presence/absence tables
p11presabs <- read.csv(file="p11pyseer_tr_presabs.txt",header=TRUE,sep="\t")
p0006presabs <- read.csv(file="p0006pyseer_tr_presabs.txt",header=TRUE,sep="\t")
p0017presabs <- read.csv(file="p0017pyseer_tr_presabs.txt",header=TRUE,sep="\t")
p0017Spresabs <- read.csv(file="p0017Spyseer_tr_presabs.txt",header=TRUE,sep="\t")
p002ypresabs <- read.csv(file="p002ypyseer_tr_presabs.txt",header=TRUE,sep="\t")
p003ppresabs <- read.csv(file="p003ppyseer_tr_presabs.txt",header=TRUE,sep="\t")
p0040presabs <- read.csv(file="p0040pyseer_tr_presabs.txt",header=TRUE,sep="\t")
pyopresabs <- read.csv(file="pyopyseer_tr_presabs.txt",header=TRUE,sep="\t")
# then load MLST conversion and ST to CC conversion
mlstold <- read.csv(file="mlst.tsv.txt", header=TRUE, sep="\t")
mlst <- read.csv(file="phageGWAS.tsv", header=TRUE, sep="\t")
STtoCC <- read.csv(file="STtoCC.tsv.txt", header=TRUE, sep="\t")
row.names(mlst) <- mlst[,1]
row.names(STtoCC) <- STtoCC[,1]
#orderedmlst <- mlst[match(row.names(mlst),row.names(orderedphagedata)),]
#orderedmlst <- orderedmlst[ !grepl("NA.*", row.names(orderedmlst)), ]
mlstST <- as.character(mlst$ST)
# add column with CC using hash
hash = new.env(hash=TRUE, parent=emptyenv(), size=100L)
# assign(key, value, hash)
# get(key, hash)
STconv <- as.character(STtoCC$ST)
CCconv <- as.list(STtoCC$CC)
STtoCCh <- hash( keys=STconv, values=CCconv )
# use hash for each ST in orderedmlst
for (i in 1:length(mlst$ST)) {
  mlst$CC[i] <- STtoCCh[[ mlstST[i] ]]
}
# orderedmlst$CC <- STtoCCh[ orderedmlstST ]
# set row names of presence/absence tables
row.names(p11presabs) <- p11presabs[,1]
row.names(p0006presabs) <- p0006presabs[,1]
row.names(p0017presabs) <- p0017presabs[,1]
row.names(p0017Spresabs) <- p0017Spresabs[,1]
row.names(p002ypresabs) <- p002ypresabs[,1]
row.names(p003ppresabs) <- p003ppresabs[,1]
row.names(p0040presabs) <- p0040presabs[,1]
row.names(pyopresabs) <- pyopresabs[,1]
# just load phagedata strains in presence/absence tables
p11presabs <- p11presabs[match(row.names(mlst),row.names(p11presabs)),]
p0006presabs <- p0006presabs[match(row.names(mlst),row.names(p0006presabs)),]
p0017presabs <- p0017presabs[match(row.names(mlst),row.names(p0017presabs)),]
p0017Spresabs <- p0017Spresabs[match(row.names(mlst),row.names(p0017Spresabs)),]
p002ypresabs <- p002ypresabs[match(row.names(mlst),row.names(p002ypresabs)),]
p003ppresabs <- p003ppresabs[match(row.names(mlst),row.names(p003ppresabs)),]
p0040presabs <- p0040presabs[match(row.names(mlst),row.names(p0040presabs)),]
pyopresabs <- pyopresabs[match(row.names(mlst),row.names(pyopresabs)),]
# remove NA strains from list
# example from earlier: orderedphagedata <- orderedphagedata[ !grepl("NA.*", row.names(orderedphagedata)), ]
p11presabs <- p11presabs[ !grepl("NA.*", row.names(p11presabs)), ]
p0006presabs <- p0006presabs[ !grepl("NA.*", row.names(p0006presabs)), ]
p0017presabs <- p0017presabs[ !grepl("NA.*", row.names(p0017presabs)), ]
p0017Spresabs <- p0017Spresabs[ !grepl("NA.*", row.names(p0017Spresabs)), ]
p002ypresabs <- p002ypresabs[ !grepl("NA.*", row.names(p002ypresabs)), ]
p003ppresabs <- p003ppresabs[ !grepl("NA.*", row.names(p003ppresabs)), ]
p0040presabs <- p0040presabs[ !grepl("NA.*", row.names(p0040presabs)), ]
pyopresabs <- pyopresabs[ !grepl("NA.*", row.names(pyopresabs)), ]
# clean up phagedata likewise
phagedata_qual_presabs <- phagedata_qual[match(row.names(mlst),row.names(phagedata_qual)),]
phagedata_qual_presabs <- phagedata_qual_presabs[ !grepl("NA.*", row.names(phagedata_qual_presabs)), ]
phagedata_presabs <- phagedata[match(row.names(mlst),row.names(phagedata)),]
phagedata_presabs <- phagedata_presabs[ !grepl("NA.*", row.names(phagedata_presabs)), ]
# make presabs data frames with ST and CC for separate models
p11presabsSTCC <- p11presabs
p0006presabsSTCC <- p0006presabs
p0017presabsSTCC <- p0017presabs
p0017SpresabsSTCC <- p0017Spresabs
p002ypresabsSTCC <- p002ypresabs
p003ppresabsSTCC <- p003ppresabs
p0040presabsSTCC <- p0040presabs
pyopresabsSTCC <- pyopresabs
# add ST and CC to datasets
p11presabsSTCC$ST <- mlst$ST
p0006presabsSTCC$ST <- mlst$ST
p0017presabsSTCC$ST <- mlst$ST
p0017SpresabsSTCC$ST <- mlst$ST
p002ypresabsSTCC$ST <- mlst$ST
p003ppresabsSTCC$ST <- mlst$ST
p0040presabsSTCC$ST <- mlst$ST
pyopresabsSTCC$ST <- mlst$ST
p11presabsSTCC$CC <- mlst$CC
p0006presabsSTCC$CC <- mlst$CC
p0017presabsSTCC$CC <- mlst$CC
p0017SpresabsSTCC$CC <- mlst$CC
p002ypresabsSTCC$CC <- mlst$CC
p003ppresabsSTCC$CC <- mlst$CC
p0040presabsSTCC$CC <- mlst$CC
pyopresabsSTCC$CC <- mlst$CC
# match pheno and presabs row names
p11presabs <- p11presabs[match(row.names(phagedata_qual_presabs),row.names(p11presabs)),]
p0006presabs <- p0006presabs[match(row.names(phagedata_qual_presabs),row.names(p0006presabs)),]
p0017presabs <- p0017presabs[match(row.names(phagedata_qual_presabs),row.names(p0017presabs)),]
p0017Spresabs <- p0017Spresabs[match(row.names(phagedata_qual_presabs),row.names(p0017Spresabs)),]
p002ypresabs <- p002ypresabs[match(row.names(phagedata_qual_presabs),row.names(p002ypresabs)),]
p003ppresabs <- p003ppresabs[match(row.names(phagedata_qual_presabs),row.names(p003ppresabs)),]
p0040presabs <- p0040presabs[match(row.names(phagedata_qual_presabs),row.names(p0040presabs)),]
pyopresabs <- pyopresabs[match(row.names(phagedata_qual_presabs),row.names(pyopresabs)),]
# add data column as final one
p11presabs$pheno <- as.factor(phagedata_qual_presabs$p11)
p0006presabs$pheno <- as.factor(phagedata_qual_presabs$p0006)
p0017presabs$pheno <- as.factor(phagedata_qual_presabs$p0017)
p0017Spresabs$pheno <- as.factor(phagedata_qual_presabs$p0017S)
p002ypresabs$pheno <- as.factor(phagedata_qual_presabs$p002y)
p003ppresabs$pheno <- as.factor(phagedata_qual_presabs$p003p)
p0040presabs$pheno <- as.factor(phagedata_qual_presabs$p0040)
pyopresabs$pheno <- as.factor(phagedata_qual_presabs$pyo)
# match pheno and presabs row names
p11presabsSTCC <- p11presabsSTCC[match(row.names(phagedata_qual_presabs),row.names(p11presabsSTCC)),]
p0006presabsSTCC <- p0006presabsSTCC[match(row.names(phagedata_qual_presabs),row.names(p0006presabsSTCC)),]
p0017presabsSTCC <- p0017presabsSTCC[match(row.names(phagedata_qual_presabs),row.names(p0017presabsSTCC)),]
p0017SpresabsSTCC <- p0017SpresabsSTCC[match(row.names(phagedata_qual_presabs),row.names(p0017SpresabsSTCC)),]
p002ypresabsSTCC <- p002ypresabsSTCC[match(row.names(phagedata_qual_presabs),row.names(p002ypresabsSTCC)),]
p003ppresabsSTCC <- p003ppresabsSTCC[match(row.names(phagedata_qual_presabs),row.names(p003ppresabsSTCC)),]
p0040presabsSTCC <- p0040presabsSTCC[match(row.names(phagedata_qual_presabs),row.names(p0040presabsSTCC)),]
pyopresabsSTCC <- pyopresabsSTCC[match(row.names(phagedata_qual_presabs),row.names(pyopresabsSTCC)),]
# also for STCC
p11presabsSTCC$pheno <- as.factor(phagedata_qual_presabs$p11)
p0006presabsSTCC$pheno <- as.factor(phagedata_qual_presabs$p0006)
p0017presabsSTCC$pheno <- as.factor(phagedata_qual_presabs$p0017)
p0017SpresabsSTCC$pheno <- as.factor(phagedata_qual_presabs$p0017S)
p002ypresabsSTCC$pheno <- as.factor(phagedata_qual_presabs$p002y)
p003ppresabsSTCC$pheno <- as.factor(phagedata_qual_presabs$p003p)
p0040presabsSTCC$pheno <- as.factor(phagedata_qual_presabs$p0040)
pyopresabsSTCC$pheno <- as.factor(phagedata_qual_presabs$pyo)
# remove Feature column - don't want to include this in random forest prediction
p11presabs <- p11presabs[,-1]
p0006presabs <- p0006presabs[,-1]
p0017presabs <- p0017presabs[,-1]
p0017Spresabs <- p0017Spresabs[,-1]
p002ypresabs <- p002ypresabs[,-1]
p003ppresabs <- p003ppresabs[,-1]
p0040presabs <- p0040presabs[,-1]
pyopresabs <- pyopresabs[,-1]
# also for STCC
p11presabsSTCC <- p11presabsSTCC[,-1]
p0006presabsSTCC <- p0006presabsSTCC[,-1]
p0017presabsSTCC <- p0017presabsSTCC[,-1]
p0017SpresabsSTCC <- p0017SpresabsSTCC[,-1]
p002ypresabsSTCC <- p002ypresabsSTCC[,-1]
p003ppresabsSTCC <- p003ppresabsSTCC[,-1]
p0040presabsSTCC <- p0040presabsSTCC[,-1]
pyopresabsSTCC <- pyopresabsSTCC[,-1]
# now just do k-mers
# first load feature presence/absence tables
p11kpresabs <- read.csv(file="p11pyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
p0006kpresabs <- read.csv(file="p0006pyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
p0017kpresabs <- read.csv(file="p0017pyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
p0017Skpresabs <- read.csv(file="p0017Spyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
p002ykpresabs <- read.csv(file="p002ypyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
p003pkpresabs <- read.csv(file="p003ppyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
p0040kpresabs <- read.csv(file="p0040pyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
pyokpresabs <- read.csv(file="pyopyseer_tr_kmer_presabs.txt",header=TRUE,sep="\t")
# set row names of presence/absence tables
row.names(p11kpresabs) <- p11kpresabs[,1]
row.names(p0006kpresabs) <- p0006kpresabs[,1]
row.names(p0017kpresabs) <- p0017kpresabs[,1]
row.names(p0017Skpresabs) <- p0017Skpresabs[,1]
row.names(p002ykpresabs) <- p002ykpresabs[,1]
row.names(p003pkpresabs) <- p003pkpresabs[,1]
row.names(p0040kpresabs) <- p0040kpresabs[,1]
row.names(pyokpresabs) <- pyokpresabs[,1]
# just load phagedata strains in presence/absence tables
p11kpresabs <- p11kpresabs[match(row.names(mlst),row.names(p11kpresabs)),]
p0006kpresabs <- p0006kpresabs[match(row.names(mlst),row.names(p0006kpresabs)),]
p0017kpresabs <- p0017kpresabs[match(row.names(mlst),row.names(p0017kpresabs)),]
p0017Skpresabs <- p0017Skpresabs[match(row.names(mlst),row.names(p0017Skpresabs)),]
p002ykpresabs <- p002ykpresabs[match(row.names(mlst),row.names(p002ykpresabs)),]
p003pkpresabs <- p003pkpresabs[match(row.names(mlst),row.names(p003pkpresabs)),]
p0040kpresabs <- p0040kpresabs[match(row.names(mlst),row.names(p0040kpresabs)),]
pyokpresabs <- pyokpresabs[match(row.names(mlst),row.names(pyokpresabs)),]
# remove NA strains from list
# example from earlier: orderedphagedata <- orderedphagedata[ !grepl("NA.*", row.names(orderedphagedata)), ]
p11kpresabs <- p11kpresabs[ !grepl("NA.*", row.names(p11kpresabs)), ]
p0006kpresabs <- p0006kpresabs[ !grepl("NA.*", row.names(p0006kpresabs)), ]
p0017kpresabs <- p0017kpresabs[ !grepl("NA.*", row.names(p0017kpresabs)), ]
p0017Skpresabs <- p0017Skpresabs[ !grepl("NA.*", row.names(p0017Skpresabs)), ]
p002ykpresabs <- p002ykpresabs[ !grepl("NA.*", row.names(p002ykpresabs)), ]
p003pkpresabs <- p003pkpresabs[ !grepl("NA.*", row.names(p003pkpresabs)), ]
p0040kpresabs <- p0040kpresabs[ !grepl("NA.*", row.names(p0040kpresabs)), ]
pyokpresabs <- pyokpresabs[ !grepl("NA.*", row.names(pyokpresabs)), ]
# clean up phagedata likewise
phagedata_qual_kpresabs <- phagedata_qual[match(row.names(mlst),row.names(phagedata_qual)),]
phagedata_qual_kpresabs <- phagedata_qual_kpresabs[ !grepl("NA.*", row.names(phagedata_qual_kpresabs)), ]
phagedata_kpresabs <- phagedata_kpresabs[ !grepl("NA.*", row.names(phagedata_kpresabs)), ]
# make presabs data frames with ST and CC for separate models
p11kpresabsSTCC <- p11kpresabs
p0006kpresabsSTCC <- p0006kpresabs
p0017kpresabsSTCC <- p0017kpresabs
p0017SkpresabsSTCC <- p0017Skpresabs
p002ykpresabsSTCC <- p002ykpresabs
p003pkpresabsSTCC <- p003pkpresabs
p0040kpresabsSTCC <- p0040kpresabs
pyokpresabsSTCC <- pyokpresabs
# add ST and CC to datasets
p11kpresabsSTCC$ST <- mlst$ST
p0006kpresabsSTCC$ST <- mlst$ST
p0017kpresabsSTCC$ST <- mlst$ST
p0017SkpresabsSTCC$ST <- mlst$ST
p002ykpresabsSTCC$ST <- mlst$ST
p003pkpresabsSTCC$ST <- mlst$ST
p0040kpresabsSTCC$ST <- mlst$ST
pyokpresabsSTCC$ST <- mlst$ST
p11kpresabsSTCC$CC <- mlst$CC
p0006kpresabsSTCC$CC <- mlst$CC
p0017kpresabsSTCC$CC <- mlst$CC
p0017SkpresabsSTCC$CC <- mlst$CC
p002ykpresabsSTCC$CC <- mlst$CC
p003pkpresabsSTCC$CC <- mlst$CC
p0040kpresabsSTCC$CC <- mlst$CC
pyokpresabsSTCC$CC <- mlst$CC
# match pheno and presabs row names
p11kpresabs <- p11kpresabs[match(row.names(phagedata_qual_kpresabs),row.names(p11kpresabs)),]
p0006kpresabs <- p0006kpresabs[match(row.names(phagedata_qual_kpresabs),row.names(p0006kpresabs)),]
p0017kpresabs <- p0017kpresabs[match(row.names(phagedata_qual_kpresabs),row.names(p0017kpresabs)),]
p0017Skpresabs <- p0017Skpresabs[match(row.names(phagedata_qual_kpresabs),row.names(p0017Skpresabs)),]
p002ykpresabs <- p002ykpresabs[match(row.names(phagedata_qual_kpresabs),row.names(p002ykpresabs)),]
p003pkpresabs <- p003pkpresabs[match(row.names(phagedata_qual_kpresabs),row.names(p003pkpresabs)),]
p0040kpresabs <- p0040kpresabs[match(row.names(phagedata_qual_kpresabs),row.names(p0040kpresabs)),]
pyokpresabs <- pyokpresabs[match(row.names(phagedata_qual_kpresabs),row.names(pyokpresabs)),]
# add data column as final one
p11kpresabs$pheno <- as.factor(phagedata_qual_kpresabs$p11)
p0006kpresabs$pheno <- as.factor(phagedata_qual_kpresabs$p0006)
p0017kpresabs$pheno <- as.factor(phagedata_qual_kpresabs$p0017)
p0017Skpresabs$pheno <- as.factor(phagedata_qual_kpresabs$p0017S)
p002ykpresabs$pheno <- as.factor(phagedata_qual_kpresabs$p002y)
p003pkpresabs$pheno <- as.factor(phagedata_qual_kpresabs$p003p)
p0040kpresabs$pheno <- as.factor(phagedata_qual_kpresabs$p0040)
pyokpresabs$pheno <- as.factor(phagedata_qual_kpresabs$pyo)
# match pheno and presabs row names
p11kpresabsSTCC <- p11kpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(p11kpresabsSTCC)),]
p0006kpresabsSTCC <- p0006kpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(p0006kpresabsSTCC)),]
p0017kpresabsSTCC <- p0017kpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(p0017kpresabsSTCC)),]
p0017SkpresabsSTCC <- p0017SkpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(p0017SkpresabsSTCC)),]
p002ykpresabsSTCC <- p002ykpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(p002ykpresabsSTCC)),]
p003pkpresabsSTCC <- p003pkpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(p003pkpresabsSTCC)),]
p0040kpresabsSTCC <- p0040kpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(p0040kpresabsSTCC)),]
pyokpresabsSTCC <- pyokpresabsSTCC[match(row.names(phagedata_qual_kpresabs),row.names(pyokpresabsSTCC)),]
# also for STCC
p11kpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$p11)
p0006kpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$p0006)
p0017kpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$p0017)
p0017SkpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$p0017S)
p002ykpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$p002y)
p003pkpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$p003p)
p0040kpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$p0040)
pyokpresabsSTCC$pheno <- as.factor(phagedata_qual_kpresabs$pyo)
# remove Feature column - don't want to include this in random forest prediction
p11kpresabs <- p11kpresabs[,-1]
p0006kpresabs <- p0006kpresabs[,-1]
p0017kpresabs <- p0017kpresabs[,-1]
p0017Skpresabs <- p0017Skpresabs[,-1]
p002ykpresabs <- p002ykpresabs[,-1]
p003pkpresabs <- p003pkpresabs[,-1]
p0040kpresabs <- p0040kpresabs[,-1]
pyokpresabs <- pyokpresabs[,-1]
# also STCC
p11kpresabsSTCC <- p11kpresabsSTCC[,-1]
p0006kpresabsSTCC <- p0006kpresabsSTCC[,-1]
p0017kpresabsSTCC <- p0017kpresabsSTCC[,-1]
p0017SkpresabsSTCC <- p0017SkpresabsSTCC[,-1]
p002ykpresabsSTCC <- p002ykpresabsSTCC[,-1]
p003pkpresabsSTCC <- p003pkpresabsSTCC[,-1]
p0040kpresabsSTCC <- p0040kpresabsSTCC[,-1]
pyokpresabsSTCC <- pyokpresabsSTCC[,-1]
# make a function to build each kind of predictive model
# do 10-fold cross-validation
# get validation accuracy from randomForest model
randomForestValidAcc <- function(train, valid) {
  model <- randomForest(pheno ~ ., data = train, importance = TRUE)
  predvalid <- predict(model, valid, type="class")
  validacc <- mean(predvalid == valid$pheno)
  return(validacc)
}
# get validation accuracy from XGBoost model
# use softmax objective, multi-class classification with three classes
xgboostValidAcc <- function(train, valid) {
  params <- list(booster = "gbtree", objective = "reg:squarederror", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)
  trainlabel <- as.numeric(train$pheno)-1
  #trainlabel <- as.factor(train$pheno)
  #levels(trainlabel) <- c("0","1","2")
  setDT(train)
  dtrain <-  xgb.DMatrix(data = model.matrix(~.+0,data = train[,-c("pheno"),with=F]), label = trainlabel) 
  validlabel <- as.numeric(valid$pheno)-1
  #validlabel <- as.factor(valid$pheno)
  #levels(validlabel) <- c("0","1","2")
  setDT(valid)
  dvalid <-  xgb.DMatrix(data = model.matrix(~.+0,data = valid[,-c("pheno"),with=F]), label = validlabel)  
  xgbcv <- xgb.cv( objective = "multi:softmax", num_class = 3, data = dtrain, nrounds = 100, nfold = 5, showsd = T, stratified = T, print_every_n = 10, early_stopping_round = 20, maximize = F, eval_metric = "merror")
  xgb1 <- xgb.train( objective = "multi:softmax", num_class = 3, data = dtrain, nrounds = xgbcv$best_iteration, watchlist = list(val=dvalid,train=dtrain), print_every_n = 10, early_stopping_round = 10, maximize = F , eval_metric = "merror")
  xgbpred <- predict (xgb1,dvalid)
  # fix label levels here
  validtest <- as.factor(valid$pheno)
  levels(validtest) <- c("0","1","2")
  xgbpred <- as.factor(xgbpred)
  levels(xgbpred) <- c("0","1","2")
  xgbCM <- confusionMatrix (xgbpred, validtest)
  xgbacc <- xgbCM$overall["Accuracy"]
}
# now do CV on randomForest
kFoldrf <- function(dataset,k) {
  flds <- createFolds(dataset$pheno, k, list = TRUE, returnTrain = FALSE)
  vec_valid <- vector(length=k)
  for (i in 1:k) {
    names(flds)[i] <- "valid"
    valid <- dataset[flds$valid,]
    train <- dataset[-flds$valid,]
    vec_valid[i] <- randomForestValidAcc(train,valid)
  }
  return(vec_valid)
}
# now do CV on xgboost
kFoldxgb <- function(dataset,k) {
  dataset$pheno <- as.factor(dataset$pheno)
  #levels(dataset$pheno) <- list("0","1","2")
  flds <- createFolds(dataset$pheno, k, list = TRUE, returnTrain = FALSE)
  vec_valid <- vector(length=k)
  for (i in 1:k) {
    names(flds)[i] <- "valid"
    valid <- dataset[flds$valid,]
    #levels(valid$pheno) <- list("0","1","2")
    train <- dataset[-flds$valid,]
    #levels(train$pheno) <- list("0","1","2")
    vec_valid[i] <- xgboostValidAcc(train,valid)
  }
  return(vec_valid)
}
# loop through each dataset this time
phages <- list(p11presabs,p0006presabs,p0017presabs,p0017Spresabs,p002ypresabs,p003ppresabs,p0040presabs,pyopresabs)
phagesSTCC <- list(p11presabsSTCC,p0006presabsSTCC,p0017presabsSTCC,p0017SpresabsSTCC,p002ypresabsSTCC,p003ppresabsSTCC,p0040presabsSTCC,pyopresabsSTCC)
phagesK <- list(p11kpresabs,p0006kpresabs,p0017kpresabs,p0017Skpresabs,p002ykpresabs,p003pkpresabs,p0040kpresabs,pyokpresabs)
phagesKSTCC <- list(p11kpresabsSTCC,p0006kpresabsSTCC,p0017kpresabsSTCC,p0017SkpresabsSTCC,p002ykpresabsSTCC,p003pkpresabsSTCC,p0040kpresabsSTCC,pyokpresabsSTCC)
# now make arrays and go through them for cross validation
phagesrfAccCV <- matrix(nrow=length(phages),ncol=10)
phagesrfSTCCAccCV <- matrix(nrow=length(phagesSTCC),ncol=10)
phagesrfKAccCV <- matrix(nrow=length(phagesK),ncol=10)
phagesrfKSTCCAccCV <- matrix(nrow=length(phagesKSTCC),ncol=10)
phagesxgbAccCV <- matrix(nrow=length(phages),ncol=10)
phagesxgbSTCCAccCV <- matrix(nrow=length(phagesSTCC),ncol=10)
phagesxgbKAccCV <- matrix(nrow=length(phagesK),ncol=10)
phagesxgbKSTCCAccCV <- matrix(nrow=length(phagesKSTCC),ncol=10)
# make array list for cross validation replicates
rfCVreplicates <- array(data=list(list(1,1,1,1),list(1,1,1,1),list(1,1,1,1),list(1,1,1,1)),dim=4)
xgbCVreplicates <- array(data=list(list(1,1,1,1),list(1,1,1,1),list(1,1,1,1),list(1,1,1,1)),dim=4)
for (j in 1:4) {
  for (i in 1:length(phages)) {
    phagesrfAccCV[i,] <- kFoldrf(phages[[i]],10)
    phagesrfSTCCAccCV[i,] <- kFoldrf(phages[[i]],10)
    phagesrfKAccCV[i,] <- kFoldrf(phages[[i]],10)
    phagesrfKSTCCAccCV[i,] <- kFoldrf(phages[[i]],10)
    phagesxgbAccCV[i,] <- kFoldxgb(phages[[i]],10)
    phagesxgbSTCCAccCV[i,] <- kFoldxgb(phages[[i]],10)
    phagesxgbKAccCV[i,] <- kFoldxgb(phages[[i]],10)
    phagesxgbKSTCCAccCV[i,] <- kFoldxgb(phages[[i]],10)
  } 
  rfCVreplicate <- list(phagesrfAccCV,phagesrfSTCCAccCV,phagesrfKAccCV,phagesrfKSTCCAccCV)
  rfCVreplicates[[j]] <- rfCVreplicate
  xgbCVreplicate <- list(phagesxgbAccCV,phagesxgbSTCCAccCV,phagesxgbKAccCV,phagesxgbKSTCCAccCV)
  xgbCVreplicates[[j]] <- xgbCVreplicate
}
fwrite(rfCVreplicates,file="rfCVreplicates.csv")
fwrite(xgbCVreplicates,file="xgbCVreplicates.csv")
# write to xlsx
write.xlsx(rfCVreplicates[[1]], file = "rfCVreplicates.xlsx",
           sheetName = "replicate1", append = FALSE)
# Add a second data set in a new worksheet
write.xlsx(rfCVreplicates[[2]], file = "rfCVreplicates.xlsx",
           sheetName = "replicate2", append = TRUE)
# Add a third data set
write.xlsx(rfCVreplicates[[3]], file = "rfCVreplicates.xlsx",
           sheetName = "replicate3", append = TRUE)
# Add a fourth data set
write.xlsx(rfCVreplicates[[4]], file = "rfCVreplicates.xlsx",
           sheetName = "replicate4", append = TRUE)
# now do xgb
write.xlsx(xgbCVreplicates[[1]], file = "xgbCVreplicates.xlsx",
           sheetName = "replicate1", append = FALSE)
# Add a second data set in a new worksheet
write.xlsx(xgbCVreplicates[[2]], file = "xgbCVreplicates.xlsx",
           sheetName = "replicate2", append = TRUE)
# Add a third data set
write.xlsx(xgbCVreplicates[[3]], file = "xgbCVreplicates.xlsx",
           sheetName = "replicate3", append = TRUE)
# Add a fourth data set
write.xlsx(xgbCVreplicates[[4]], file = "xgbCVreplicates.xlsx",
           sheetName = "replicate4", append = TRUE)
# Write down strains, ST, and CC
# write.table(strain_list, file = "strain_list.txt", sep = "\t")
