#This file runs differential expression on the data
# File also outputs key files for python classification
#BiocManager::install("edgeR")

library("edgeR")
library("limma")
require(xlsx)
library("tidyverse")
library(dplyr)
library(e1071)


library(mlbench)
library(mlr3)
library(mlr3viz)
library(precrec)
library(mlr3learners)
require(xlsx)
library("caret")
library("mlbench")
library("pROC")
library(ROCR)
library("rpart")
library("caretEnsemble")
library("mlbench")
library("randomForest")
library("nnet")
library("caTools")
library("gbm")
library(Boruta)
library(blkbox)

# Required Packages

library('GEOquery')
library('Biobase')

#for PSPBP data import paths


nih_data<-"/Users/michaelallwright/Dropbox (Sydney Uni)/michael_PhD/Projects/Parkinson's Longitudinal Study/Data/Raw/PDBP/"
nih_files <- dir(path=nih_data, pattern="*\\.txt",recursive = TRUE)

p_dat<-"/Users/michaelallwright/Documents/python/PDBP/Project Final/data/"

nih_lkup<-read.xlsx(paste0(p_dat,"Mapping Sample to Diagnosis 20200122.xlsx"),sheetName = "Calcs")

View(nih_lkup)

#setwd("~/Dropbox (Sydney Uni)/PDBP Hack/Data/Raw PDBP")
x_nih <- read.maimages(nih_files, source="agilent", green.only=TRUE)

x_nih$genes

############### background correct and normalise
#Both on biofind files and nih data

y_nih <- backgroundCorrect(x_nih,method="normexp")
y_nih <- normalizeBetweenArrays(y_nih,method="quantile")

#We will filter out control probes as indicated by the ControlType column

Control <- y_nih$genes$ControlType==1L

#Finally, we will filter probes that don’t appear to be expressed. 
#We keep probes that are above background on at least four arrays (because there are four replicates of each treatment):

IsExpr <- rowSums(y_nih$other$gIsWellAboveBG > 0) >= 4

nih_full=as.data.frame(y_nih$E)
nih_full$ProbeName=y_nih$genes$ProbeName
nih_full = nih_full[!duplicated(nih_full$ProbeName)&!Control,]

nih_lkup<-read.csv("/Users/admin/Dropbox (Sydney Uni)/PDBP Hack/Data/Labels/LabelsPDBP.csv")


nih_lkupA<-nih_lkup[(nih_lkup$Sample.Set=="A"),]
#nih_lkupA<-nih_lkup
nih_lkup2<-select(nih_lkupA,RNAFileName, Diagnosis)
nih_lkup3<-select(nih_lkupA,RNAFileName, Diagnosis,Sample.Set)

saveRDS(nih_lkup2,file="label.rds")
saveRDS(nih_full,file="data.rds")

############### Log2 transform
#Both on biofind files and nih data
y_nih_lt<-y_nih
ex <- y_nih_lt$E
qx <-as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  y_nih_lt$E <- log2(ex)
}

nih_full=as.data.frame(y_nih_lt$E)
nih_full$ProbeName=y_nih_lt$genes$ProbeName
nih_full = nih_full[!duplicated(nih_full$ProbeName),]
rownames(nih_full)=nih_full$ProbeName
nih_full_int=nih_full[,names(nih_full)!="ProbeName"]
nih_full_int_t=as.data.frame(t(nih_full))
nih_full_int_t$sample=rownames(nih_full_int_t)




#####################################################################
mod_data_full<-merge(nih_full_int_t, nih_lkup3, by.x="sample", by.y="RNAFileName")
write_csv(mod_data_full,paste0(p_dat,"mod_dataPDBPv2.csv"))


mod_data<-merge(nih_full_int_t, nih_lkup2, by.x="sample", by.y="RNAFileName")
mod_data <- mod_data[((mod_data$Diagnosis == "PD")|(mod_data$Diagnosis == "HC")),]
mod_data$Diagnosis<-as.factor(mod_data$Diagnosis)

#dim(mod_data)

#mod_data2=mod_data[((mod_data$Diagnosis != "PSP")|(mod_data$Diagnosis != "MSA"))]
levels(mod_data$Diagnosis)
mod_data <- mod_data %>%
  droplevels
levels(mod_data$Diagnosis)
#View(mod_data$Diagnosis)

dim(mod_data)

mod_datafin=mod_data
mod_data<-mod_data[,names(mod_data)!="sample"]
write_csv(mod_datafin,paste0(p_dat,"mod_dataPDBP.csv"))

saveRDS(mod_data,"/Users/michaelallwright/Dropbox (Sydney Uni)/michael_PhD/PDBP Hack/Data/Processed Data/mod_data.RDS")

write.csv()
###### Below is to perform differential expression and run models on result

### function to perform differential expression

diffexfunc <- function(data,num) {
  design <- model.matrix( ~ Diagnosis+0, data)
  fit <- lmFit(t(data[,names(data)!="Diagnosis"]), design)
  cont.matrix <- makeContrasts( DiagnosisPD-DiagnosisHC, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  
  #100 most differentially expressed genes
  tT <- topTable(fit2,
                 adjust = "fdr",
                 sort.by = "B",
                 number = num)
  out <- list()
  out$data_sel <- data[,c(rownames(tT),"Diagnosis")]
  out$feats=rownames(tT)
  return(out)
}

#change variables for models


#DE on entire data
moddata_de=diffexfunc(mod_data,2000)
mod_data_defirst<-moddata_de$data_sel

#split into training and testing after
data_set_size <- floor(2*nrow(mod_data_defirst)/3)
indexes <- sample(1:nrow(mod_data_defirst), size = data_set_size)
training_defirst <- mod_data_defirst[indexes,]
validation1_defirst <- mod_data_defirst[-indexes,]

dim(mod_data)
#split to training and testing first
data_set_size <- floor(2*nrow(mod_data)/3)
indexes <- sample(1:nrow(mod_data), size = data_set_size)
training <- mod_data[indexes,]
validation1 <- mod_data[-indexes,]

#DE on training after split
training_de=diffexfunc(training,2000)
feats<-training_de$feats
training_de_after<-training_de$data_sel
validation_de_after<-validation1[,c(feats,"Diagnosis")]

#Output

write_csv(training_defirst,paste0(p_dat,"NIH_trainingDE.csv"))
write_csv(validation1_defirst,paste0(p_dat,"NIH_validation1DE.csv"))
write_csv(training_de_after,paste0(p_dat,"NIH_training_de_after.csv"))
write_csv(validation_de_after,paste0(p_dat,"NIH_validation_de_after.csv"))

require(xgboost)

#### Run model

svmmodel <- function(traindata,valdata) {
  valdata <- valdata %>%
    mutate(Diagnosis = ifelse(Diagnosis == "PD",0,1))
  
  traindata <- traindata %>%
    mutate(Diagnosis = ifelse(Diagnosis == "PD",0,1))
  
  model <- svm(Diagnosis ~ ., data = traindata)
  valdata$pred <- predict(model, valdata[,names(valdata)!="Diagnosis"])
  table(valdata$pred,valdata$Diagnosis)
  roc_obj <- roc(valdata$Diagnosis,valdata$pred)
  auc=roc_obj$auc[1]
  return(auc)
}


svmmodel(training_de_after,validation_de_after)
svmmodel(training_defirst,validation1_defirst)


#Loop all 100 times outputting AUC
n=100
auc_bef=1
auc_aft=1
for(i in 1:n)
{
  #DE on entire data
  moddata_de=diffexfunc(mod_data,2000)
  mod_data_defirst<-moddata_de$data_sel
  
  #split into training and testing after
  data_set_size <- floor(2*nrow(mod_data_defirst)/3)
  indexes <- sample(1:nrow(mod_data_defirst), size = data_set_size)
  training_defirst <- mod_data_defirst[indexes,]
  validation1_defirst <- mod_data_defirst[-indexes,]
  
  #split to training and testing first
  data_set_size <- floor(2*nrow(mod_data)/3)
  indexes <- sample(1:nrow(mod_data), size = data_set_size)
  training <- mod_data[indexes,]
  validation1 <- mod_data[-indexes,]
  
  #DE on training after split
  training_de=diffexfunc(training,2000)
  feats<-training_de$feats
  training_de_after<-training_de$data_sel
  validation_de_after<-validation1[,c(feats,"Diagnosis")]
  
  auc_bef[i] = svmmodel(training_defirst,validation1_defirst)
  auc_aft[i] = svmmodel(training_de_after,validation_de_after)
}

boxplot(auc_bef,auc_aft)

boxplot(auc_bef, auc_aft,
        main = "Comparing Differential Expression Methods SVM",
        #at = c(1,2,4,5),
        names = c("Old V", "New V"),
        las = 2,
        col = c("orange","red"),
        border = "brown",
        horizontal = TRUE,
        notch = TRUE
)

#######################
#3####################
require(dplyr)
validation1_defirst <- validation1_defirst %>%
  mutate(Diagnosis = ifelse(Diagnosis == "PD",0,1))

training_defirst <- training_defirst %>%
  mutate(Diagnosis = ifelse(Diagnosis == "PD",0,1))

validation_de_after <- validation_de_after %>%
  mutate(Diagnosis = ifelse(Diagnosis == "PD",0,1))

training_de_after <- training_de_after %>%
  mutate(Diagnosis = ifelse(Diagnosis == "PD",0,1))



model <- svm(Diagnosis ~ ., data = training_defirst)
validation1_defirst$pred <- predict(model, X_test)
table(validation1_defirst$pred,validation1_defirst$Diagnosis)
roc_obj <- roc(validation1_defirst$Diagnosis,validation1_defirst$pred)
auc=roc_obj$auc[1]


#bstSparse <- xgboost(data = data.matrix(training_defirst), label = factor(training_defirst$Diagnosis, levels = c("yes", "no")), max.depth = 2,
#                     eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic")

library(randomForest)
require(caTools)
rf <- randomForest(Diagnosis ~ .,data=training_defirst)

View(training_defirst)


View(validation1_sel)

validation1_sel <- validation1[,c(rownames(tT),"Diagnosis")]

View(training_sel)
library(randomForest)
# Perform training:
rf_classifier = randomForest(Diagnosis ~ ., data=training, ntree=20, mtry=2, importance=TRUE)

prediction_for_table <- predict(rf_classifier,validation1[,c(features)])
table(observed=validation1[,c("Diagnosis")],predicted=prediction_for_table)



#############run simple cv on this



ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     number = 5,
                     classProbs = TRUE,
                     savePredictions = "final",
                     allowParallel=TRUE,
                     returnData=FALSE)

tl <- 30
#tl for xgbtree shorter because it has 7 tuning parameters compared to 1-4 for others
tlx <- 10

rf_model <- train(x = training[, names(training) != "Diagnosis"],
                  y = training$Diagnosis,
                  method = "rf",
                  tuneLength = tl,
                  trControl = ctrl,
                  metric = "Kappa")

rf_model$results

#predict on validation data




task = TaskClassif$new("mod_data_sel", backend = mod_data_sel, target = "Diagnosis", positive = "PD")
#task = TaskClassif$new("BF_mod", backend = BF_mod, target = "diagnosis", positive = "PD")

#train and test

train_set = sample(task$nrow, 0.6 * task$nrow)
test_set = setdiff(seq_len(task$nrow), train_set)
learner = lrn("classif.rpart", predict_type = "prob")
#learner = mlr_learners$get("classif.xgboost")
#learner = mlr_learners$get("classif.rpart")

#key to set row_ids as trian_set 0- if you don't you're training on all data and putting bias in
learner$train(task, row_ids = train_set)

#learner$train(task)

prediction = learner$predict(task, row_ids = test_set)
prediction$confusion
autoplot(prediction)
autoplot(prediction, type = "roc", row_ids = test_set)

View(mod_data_sel)



######################################################################

# Load Data

## load series and platform data from GEO
gset = GEOquery::getGEO("GSE99039", GSEMatrix =TRUE, AnnotGPL=TRUE)

if (length(gset) > 1) idx = grep("GPL570", attr(gset, "names")) else idx = 1
gset = gset[[idx]]

View(gset)

Diagnosis<-gset$`disease label:ch1`
Sample<-gset$geo_accession

labelsfile<-data.frame(cbind(Sample,Diagnosis))


View(labelsfile)
## Make coloumn names match toptable
fvarLabels(gset) = make.names(fvarLabels(gset))

## group names for all samples
gsms = base::paste0(
  "011011111X1010010XX00001X000X00001X111XXX001101111",
  "XX0X01X1010XX1111111110XXXX1X1110110XXXXXXXXX01010",
  "XXX01000010010000X00000000001X0X1111X1001011111000",
  "00010011001XX10000X000010XXXXX11111X11X11111XXXX11",
  "0110XX10111110110111110111011X11011011X110X0111011",
  "1100011000X101000000001000X11000000000001011110000",
  "0010XX00110X0100000010000000XXX110XX1111X100001011",
  "10000X1101011000100111000010000X0X001X1X1011001111",
  "0010001010111X00X000100011101111101X00000X1XXXXXXX",
  "XXXXXXXXXXXX10XX11X0XX000X0000000000111101XX011X11",
  "111111111111X10000000000X00XX0X011X1X0110X00XX000X",
  "XX0X0XXX"
)
sml = c()
for (i in 1:nchar(gsms)) {
  sml[i] = base::substr(gsms, i, i)
}

## Eliminate samples marked as "X"
sel = which(sml != "X")
sml = sml[sel]
gset = gset[, sel]

## log2 transform
ex = Biobase::exprs(gset)
qx = as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}


View(ex)

transcripts<-data.frame(t(ex))
View(transcripts)
dim(transcripts)

transcripts$Sample=rownames(transcripts)

#merge to create training data

mod_data_geo<-merge(transcripts, labelsfile, by.x="Sample", by.y="Sample")
mod_data_geo$Diagnosis

mod_data_geo <- mod_data_geo[((mod_data_geo$Diagnosis == "IPD")|(mod_data_geo$Diagnosis == "CONTROL"))]
mod_data_geo$Diagnosis<-as.factor(mod_data_geo$Diagnosis)

#mod_data2=mod_data[((mod_data$Diagnosis != "PSP")|(mod_data$Diagnosis != "MSA"))]
levels(mod_data_geo$Diagnosis)
mod_data_geo <- mod_data_geo %>%
  droplevels
levels(mod_data_geo$Diagnosis)


View(mod_data_geo)

write.csv(mod_data_geo,paste0(p_dat,"mod_data_geo.csv"))

#split to training and testing first
data_set_size <- floor(2*nrow(mod_data_geo)/3)
indexes <- sample(1:nrow(mod_data_geo), size = data_set_size)
training_geo <- mod_data_geo[indexes,]
validation1_geo <- mod_data_geo[-indexes,]

svmmodel(training_geo,validation1_geo)

