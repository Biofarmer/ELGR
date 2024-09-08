library(plyr)
library(dplyr)
library(tibble)
library(caret)
library(randomForest)
library(ROCR)

#subtype
window_db_subtype <- read.csv("data/window_db_subtype.csv",check.names=FALSE,stringsAsFactors = FALSE)
#species
metaphlan4_preterm <- read.csv("data/metaphlan4_preterm.csv")
metaphlan4_preterm[is.na(metaphlan4_preterm)] <- 0
metaphlan4_preterm_2 <- metaphlan4_preterm%>%dplyr::select(clade_name,unique(window_db_subtype$SequenceID))
metaphlan4_preterm_2_kingdom <- metaphlan4_preterm_2%>%filter(grepl("k__",clade_name))%>%filter(!grepl("p__",clade_name))
metaphlan4_species <- metaphlan4_preterm_2%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
metaphlan4_species$species <- sapply(strsplit(metaphlan4_species$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
metaphlan4_species <- metaphlan4_species%>%dplyr::select(-clade_name)%>%column_to_rownames("species")

metadata_prenec <- read.csv("data/metadata_prenec_closest.csv")

window_db <- read.csv("data/window_db.csv")
window_db_047 <- window_db%>%filter(window%in%c(0,4,7))

#subtype
window_db_subtype_rel <- window_db_subtype%>%column_to_rownames("SequenceID")
window_db_subtype_rel <-  as.data.frame(prop.table(as.matrix(window_db_subtype_rel), 1))
window_db_subtype_rel <- window_db_subtype_rel%>%rownames_to_column("SequenceID")
window_db_subtype_047 <- window_db_subtype_rel%>%filter(SequenceID%in%window_db_047$SequenceID)%>%column_to_rownames("SequenceID")
window_db_subtype_047 <- window_db_subtype_047[,which(colSums(window_db_subtype_047)!=0)]
window_db_subtype_047 <- window_db_subtype_047%>%rownames_to_column("SequenceID")
window_db_sel_subtype_meta <- right_join(metadata_prenec%>%dplyr::select(SequenceID,Study,nec),window_db_subtype_047,by="SequenceID")
window_db_sel_subtype_meta <- window_db_sel_subtype_meta%>%filter(Study%in%c("BrooksB_2017","MasiAC_2021","RahmanSF_2018","RavehSadkaT_2016"))%>%
   dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)

#species
metaphlan4_species_t <- metaphlan4_species%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
metaphlan4_species_t[-1] <- metaphlan4_species_t[-1]/100
window_db_species_047 <- metaphlan4_species_t%>%filter(SequenceID%in%window_db_047$SequenceID)%>%column_to_rownames("SequenceID")
window_db_species_047 <- window_db_species_047[,which(colSums(window_db_species_047)!=0)]
window_db_species_047 <- window_db_species_047%>%rownames_to_column("SequenceID")
window_db_sel_species_meta <- right_join(metadata_prenec%>%dplyr::select(SequenceID,Study,nec),window_db_species_047,by="SequenceID")
window_db_sel_species_meta <- window_db_sel_species_meta%>%filter(Study%in%c("BrooksB_2017","MasiAC_2021","RahmanSF_2018","RavehSadkaT_2016"))%>%
   dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)

#subtype + species
window_db_sel_subtype_species_meta <- left_join(window_db_sel_subtype_meta,window_db_species_047%>%dplyr::select(SequenceID,names(window_db_sel_species_meta)[-c(1:3)]),by="SequenceID")

########################################################################################################################
#predictive model training based on four public cohorts: 
#features: subtype, species, subtype+species 
#approaches: intra cohort, combined cohort, LOSO
########################################################################################################################
#subtype
#intra cohort prediction
rf_auc_intra_subtype_047 <- data.frame("study"="x","time"="x","auc"=0,"fold"=0,"mtry"=0,stringsAsFactors = F)[-1,]
for (s in unique(window_db_sel_subtype_meta$Study)) {
   for (c in seq(1,10,1)) {
      window_db_sel_subtype_meta_study <- window_db_sel_subtype_meta%>%filter(Study==s)%>%dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)
      set.seed(c*100)
      dataset_index <- createFolds(window_db_sel_subtype_meta_study$nec, k = 5, list = TRUE, returnTrain = FALSE)
      for (f in paste0("Fold",seq(1,5,1))) {
         test_set <- window_db_sel_subtype_meta_study[dataset_index[[f]],]%>%dplyr::select(-Study)
         rownames(test_set) <- NULL
         train_set <- window_db_sel_subtype_meta_study%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
         x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
         y_train <- as.factor(train_set$nec)
         for (m in c(seq(1,ncol(x_train),1))) {
            set.seed(1000)
            rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m)
            x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
            y_test <- test_set$nec
            rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
            rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
            names(rf_pred) <- "pred"
            #to calucate AUC
            comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
            names(comb)[1] <- "obs"
            comb$obs <- as.factor(comb$obs)
            comb$pred <- as.factor(comb$pred)
            comb$no <- as.numeric(comb$no)
            comb$yes <- as.numeric(comb$yes)
            pred <- ROCR::prediction(comb[,3], comb[,1])
            auc <- ROCR::performance(pred,"auc")@y.values[[1]]
            rf_auc_intra_subtype_047[nrow(rf_auc_intra_subtype_047)+1,] <- c(s,c,auc,f,m)
         }
      }
   }
}

#combined cohort
rf_auc_cb_subtype_047 <- data.frame("time"="x","auc"=0,"fold"=0,"mtry"=0,stringsAsFactors = F)[-1,]
for (c in seq(1,10,1)) {
   # c=3
   print(c)
   set.seed(c*100)
   dataset_index <- createFolds(window_db_sel_subtype_meta$nec, k = 5, list = TRUE, returnTrain = FALSE)
   for (f in paste0("Fold",seq(1,5,1))) {
      # f="Fold1"
      print(f)
      test_set <- window_db_sel_subtype_meta[dataset_index[[f]],]%>%dplyr::select(-Study)
      rownames(test_set) <- NULL
      train_set <- window_db_sel_subtype_meta%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
      x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
      y_train <- as.factor(train_set$nec)
      for (m in c(seq(1,ncol(x_train),1))) {
         set.seed(1000)
         rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m)
         x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
         y_test <- test_set$nec
         rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
         rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
         names(rf_pred) <- "pred"
         #to calucate AUC
         comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
         names(comb)[1] <- "obs"
         comb$obs <- as.factor(comb$obs)
         comb$pred <- as.factor(comb$pred)
         comb$no <- as.numeric(comb$no)
         comb$yes <- as.numeric(comb$yes)
         pred <- ROCR::prediction(comb[,3], comb[,1])
         auc <- ROCR::performance(pred,"auc")@y.values[[1]]
         rf_auc_cb_subtype_047[nrow(rf_auc_cb_subtype_047)+1,] <- c(c,auc,f,m)
      }
   }
}
#species
#intra cohort prediction
rf_auc_intra_species_047 <- data.frame("study"="x","time"="x","auc"=0,"fold"=0,"mtry"=0,stringsAsFactors = F)[-1,]
for (s in unique(window_db_sel_species_meta$Study)) {
   for (c in seq(1,10,1)) {
      window_db_sel_species_meta_study <- window_db_sel_species_meta%>%filter(Study==s)%>%dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)
      set.seed(c*100)
      dataset_index <- createFolds(window_db_sel_species_meta_study$nec, k = 5, list = TRUE, returnTrain = FALSE)
      for (f in paste0("Fold",seq(1,5,1))) {
         test_set <- window_db_sel_species_meta_study[dataset_index[[f]],]%>%dplyr::select(-Study)
         rownames(test_set) <- NULL
         train_set <- window_db_sel_species_meta_study%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
         x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
         y_train <- as.factor(train_set$nec)
         for (m in c(seq(1,ncol(x_train),1))) {
            set.seed(1000)
            rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m)
            x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
            y_test <- test_set$nec
            rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
            rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
            names(rf_pred) <- "pred"
            #to calucate AUC
            comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
            names(comb)[1] <- "obs"
            comb$obs <- as.factor(comb$obs)
            comb$pred <- as.factor(comb$pred)
            comb$no <- as.numeric(comb$no)
            comb$yes <- as.numeric(comb$yes)
            pred <- ROCR::prediction(comb[,3], comb[,1])
            auc <- ROCR::performance(pred,"auc")@y.values[[1]]
            rf_auc_intra_species_047[nrow(rf_auc_intra_species_047)+1,] <- c(s,c,auc,f,m)
         }
      }
   }
}
#combined cohort
rf_auc_cb_species_047 <- data.frame("time"="x","auc"=0,"fold"=0,"mtry"=0,stringsAsFactors = F)[-1,]
for (c in seq(1,10,1)) {
   # c=3
   print(c)
   set.seed(c*100)
   dataset_index <- createFolds(window_db_sel_species_meta$nec, k = 5, list = TRUE, returnTrain = FALSE)
   for (f in paste0("Fold",seq(1,5,1))) {
      # f="Fold1"
      print(f)
      test_set <- window_db_sel_species_meta[dataset_index[[f]],]%>%dplyr::select(-Study)
      rownames(test_set) <- NULL
      train_set <- window_db_sel_species_meta%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
      x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
      y_train <- as.factor(train_set$nec)
      for (m in c(seq(1,ncol(x_train),1))) {
         set.seed(1000)
         rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m)
         x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
         y_test <- test_set$nec
         rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
         rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
         names(rf_pred) <- "pred"
         #to calucate AUC
         comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
         names(comb)[1] <- "obs"
         comb$obs <- as.factor(comb$obs)
         comb$pred <- as.factor(comb$pred)
         comb$no <- as.numeric(comb$no)
         comb$yes <- as.numeric(comb$yes)
         pred <- ROCR::prediction(comb[,3], comb[,1])
         auc <- ROCR::performance(pred,"auc")@y.values[[1]]
         rf_auc_cb_species_047[nrow(rf_auc_cb_species_047)+1,] <- c(c,auc,f,m)
      }
   }
}
#subtype+species
#intra cohort prediction
rf_auc_intra_subtype_species_047 <- data.frame("study"="x","time"="x","auc"=0,"fold"=0,"mtry"=0,stringsAsFactors = F)[-1,]
for (s in unique(window_db_sel_subtype_species_meta$Study)) {
   for (c in seq(1,10,1)) {
      window_db_sel_subtype_species_meta_study <- window_db_sel_subtype_species_meta%>%filter(Study==s)%>%dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)
      set.seed(c*100)
      dataset_index <- createFolds(window_db_sel_subtype_species_meta_study$nec, k = 5, list = TRUE, returnTrain = FALSE)
      for (f in paste0("Fold",seq(1,5,1))) {
         test_set <- window_db_sel_subtype_species_meta_study[dataset_index[[f]],]%>%dplyr::select(-Study)
         rownames(test_set) <- NULL
         train_set <- window_db_sel_subtype_species_meta_study%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
         x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
         y_train <- as.factor(train_set$nec)
         for (m in c(seq(1,ncol(x_train),1))) {
            set.seed(1000)
            rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m)
            x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
            y_test <- test_set$nec
            rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
            rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
            names(rf_pred) <- "pred"
            #to calucate AUC
            comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
            names(comb)[1] <- "obs"
            comb$obs <- as.factor(comb$obs)
            comb$pred <- as.factor(comb$pred)
            comb$no <- as.numeric(comb$no)
            comb$yes <- as.numeric(comb$yes)
            pred <- ROCR::prediction(comb[,3], comb[,1])
            auc <- ROCR::performance(pred,"auc")@y.values[[1]]
            rf_auc_intra_subtype_species_047[nrow(rf_auc_intra_subtype_species_047)+1,] <- c(s,c,auc,f,m)
         }
      }
   }
}
#combined cohort
rf_auc_cb_subtype_species_047 <- data.frame("time"="x","auc"=0,"fold"=0,"mtry"=0,stringsAsFactors = F)[-1,]
for (c in seq(1,10,1)) {
   # c=3
   print(c)
   set.seed(c*100)
   dataset_index <- createFolds(window_db_sel_subtype_species_meta$nec, k = 5, list = TRUE, returnTrain = FALSE)
   for (f in paste0("Fold",seq(1,5,1))) {
      # f="Fold1"
      print(f)
      test_set <- window_db_sel_subtype_species_meta[dataset_index[[f]],]%>%dplyr::select(-Study)
      rownames(test_set) <- NULL
      train_set <- window_db_sel_subtype_species_meta%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
      x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
      y_train <- as.factor(train_set$nec)
      for (m in c(seq(1,ncol(x_train),1))) {
         set.seed(1000)
         rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m)
         x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
         y_test <- test_set$nec
         rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
         rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
         names(rf_pred) <- "pred"
         #to calucate AUC
         comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
         names(comb)[1] <- "obs"
         comb$obs <- as.factor(comb$obs)
         comb$pred <- as.factor(comb$pred)
         comb$no <- as.numeric(comb$no)
         comb$yes <- as.numeric(comb$yes)
         pred <- ROCR::prediction(comb[,3], comb[,1])
         auc <- ROCR::performance(pred,"auc")@y.values[[1]]
         rf_auc_cb_subtype_species_047[nrow(rf_auc_cb_subtype_species_047)+1,] <- c(c,auc,f,m)
      }
   }
}

####################################################################################################################################
#validation with independent cohort
#featues: subtype, species, subtype+species 
####################################################################################################################################
#subtype
window_db_subtype <- read.csv("data/window_db_subtype.csv",check.names=FALSE,stringsAsFactors = FALSE)
window_db_subtype_validation <- read.csv("data/window_db_sel_subtype_meta_validation.csv",check.names=FALSE,stringsAsFactors = FALSE)
#species
metaphlan4_preterm <- read.csv("data/metaphlan4_preterm.csv")
metaphlan4_preterm[is.na(metaphlan4_preterm)] <- 0
metaphlan4_preterm_2 <- metaphlan4_preterm%>%dplyr::select(clade_name,unique(window_db_subtype$SequenceID))
metaphlan4_species <- metaphlan4_preterm_2%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
metaphlan4_species$species <- sapply(strsplit(metaphlan4_species$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
metaphlan4_species <- metaphlan4_species%>%dplyr::select(-clade_name)%>%column_to_rownames("species")
#validation
metaphlan4_validation <- read.csv("data/metaphlan4_validation.csv")
metaphlan4_validation[is.na(metaphlan4_validation)] <- 0
metaphlan4_validation_2 <- metaphlan4_validation%>%dplyr::select(clade_name,unique(window_db_subtype_validation$SequenceID))
metaphlan4_species_validation <- metaphlan4_validation_2%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
metaphlan4_species_validation$species <- sapply(strsplit(metaphlan4_species_validation$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
metaphlan4_species_validation <- metaphlan4_species_validation%>%dplyr::select(-clade_name)%>%column_to_rownames("species")
metaphlan4_species_validation_rel <- metaphlan4_species_validation/100
metaphlan4_species_validation_rel_t <- metaphlan4_species_validation_rel%>%t()%>%as.data.frame()
metaphlan4_species_validation_rel_t <- metaphlan4_species_validation_rel_t[,which(colSums(metaphlan4_species_validation_rel_t)!=0)]
metaphlan4_species_validation_rel_t <- metaphlan4_species_validation_rel_t%>%rownames_to_column("SequenceID")
metaphlan4_species_validation_rel_t_meta <- right_join(window_db_subtype_validation%>%dplyr::select(1:3),
                                                       metaphlan4_species_validation_rel_t,by="SequenceID")
window_db_047 <- window_db%>%filter(window%in%c(0,4,7))

#subtype
window_db_subtype_rel <- window_db_subtype%>%column_to_rownames("SequenceID")
window_db_subtype_rel <-  as.data.frame(prop.table(as.matrix(window_db_subtype_rel), 1))
window_db_subtype_rel <- window_db_subtype_rel%>%rownames_to_column("SequenceID")
window_db_subtype_047 <- window_db_subtype_rel%>%filter(SequenceID%in%window_db_047$SequenceID)%>%column_to_rownames("SequenceID")
window_db_subtype_047 <- window_db_subtype_047[,which(colSums(window_db_subtype_047)!=0)]
window_db_subtype_047 <- window_db_subtype_047%>%rownames_to_column("SequenceID")
window_db_sel_subtype_meta <- right_join(metadata_prenec%>%dplyr::select(SequenceID,Study,nec),window_db_subtype_047,by="SequenceID")
window_db_sel_subtype_meta <- window_db_sel_subtype_meta%>%filter(Study%in%c("BrooksB_2017","MasiAC_2021","RahmanSF_2018","RavehSadkaT_2016"))%>%
   dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)
window_db_sel_subtype_meta <- rbind.fill(window_db_sel_subtype_meta,window_db_subtype_validation)
window_db_sel_subtype_meta[is.na(window_db_sel_subtype_meta)] <- 0

#species
metaphlan4_species_t <- metaphlan4_species%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
metaphlan4_species_t[-1] <- metaphlan4_species_t[-1]/100
window_db_species_047 <- metaphlan4_species_t%>%filter(SequenceID%in%window_db_047$SequenceID)%>%column_to_rownames("SequenceID")
window_db_species_047 <- window_db_species_047[,which(colSums(window_db_species_047)!=0)]
window_db_species_047 <- window_db_species_047%>%rownames_to_column("SequenceID")
window_db_sel_species_meta <- right_join(metadata_prenec%>%dplyr::select(SequenceID,Study,nec),window_db_species_047,by="SequenceID")
window_db_sel_species_meta <- window_db_sel_species_meta%>%filter(Study%in%c("BrooksB_2017","MasiAC_2021","RahmanSF_2018","RavehSadkaT_2016"))%>%
   dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)
window_db_sel_species_meta <- rbind.fill(window_db_sel_species_meta,metaphlan4_species_validation_rel_t_meta)
window_db_sel_species_meta[is.na(window_db_sel_species_meta)] <- 0

#subtype + species
window_db_sel_subtype_species_meta <- left_join(window_db_sel_subtype_meta,window_db_sel_species_meta%>%dplyr::select(-c(2:3)),by="SequenceID")

#subtype
rf_auc_subtype_validation_m <- data.frame("study"="x","auc"=0,"mtry"=0,stringsAsFactors = F)[-1,]
train_set <- window_db_sel_subtype_meta %>% filter(Study!="validation")%>%dplyr::select(-Study)
x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train <- as.factor(train_set$nec)
for (m in c(seq(1,ncol(x_train),1))) {
   set.seed(1000)
   rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE,mtry=m)
   test_set <- window_db_sel_subtype_meta %>% filter(Study=="validation")%>%dplyr::select(-Study)
   x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
   y_test <- test_set$nec
   rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
   rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
   names(rf_pred) <- "pred"
   #to calucate AUC
   comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
   names(comb)[1] <- "obs"
   comb$obs <- as.factor(comb$obs)
   comb$pred <- as.factor(comb$pred)
   comb$no <- as.numeric(comb$no)
   comb$yes <- as.numeric(comb$yes)
   pred <- ROCR::prediction(comb[,3], comb[,1])
   auc <- ROCR::performance(pred,"auc")@y.values[[1]]
   rf_auc_subtype_validation_m[nrow(rf_auc_subtype_validation_m)+1,] <- c("validation",auc,m)
}
#species
rf_auc_species_validation_m <- data.frame("study"="x","auc"=0,"mtry"=0,stringsAsFactors = F)[-1,]
train_set <- window_db_sel_species_meta %>% filter(Study!="validation")%>%dplyr::select(-Study)
x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train <- as.factor(train_set$nec)
for (m in c(seq(1,ncol(x_train),1))) {
   set.seed(1000)
   rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE,mtry=m)
   test_set <- window_db_sel_species_meta %>%filter(Study=="validation")%>%dplyr::select(-Study)
   x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
   y_test <- test_set$nec
   rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
   rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
   names(rf_pred) <- "pred"
   #to calucate AUC
   comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
   names(comb)[1] <- "obs"
   comb$obs <- as.factor(comb$obs)
   comb$pred <- as.factor(comb$pred)
   comb$no <- as.numeric(comb$no)
   comb$yes <- as.numeric(comb$yes)
   pred <- ROCR::prediction(comb[,3], comb[,1])
   auc <- ROCR::performance(pred,"auc")@y.values[[1]]
   rf_auc_species_validation_m[nrow(rf_auc_species_validation_m)+1,] <- c("validation",auc,m)
}
#subtype+species
rf_auc_subtype_species_validation_m <- data.frame("study"="x","auc"=0,"mtry"=0,stringsAsFactors = F)[-1,]
train_set <- window_db_sel_subtype_species_meta %>% filter(Study!="validation")%>%dplyr::select(-Study)
x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train <- as.factor(train_set$nec)
for (m in c(seq(1,ncol(x_train),1))) {
   set.seed(1000)
   rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE,mtry=m)
   test_set <- window_db_sel_subtype_species_meta %>% filter(Study=="validation")%>%dplyr::select(-Study)
   x_test <- test_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
   y_test <- test_set$nec
   rf_pred_prob <- predict(rf_train, newdata = x_test, type = "prob")
   rf_pred <- predict(rf_train, newdata = x_test)%>%as.data.frame()
   names(rf_pred) <- "pred"
   #to calucate AUC
   comb <- cbind(y_test,rf_pred_prob,rf_pred)%>%as.data.frame()
   names(comb)[1] <- "obs"
   comb$obs <- as.factor(comb$obs)
   comb$pred <- as.factor(comb$pred)
   comb$no <- as.numeric(comb$no)
   comb$yes <- as.numeric(comb$yes)
   pred <- ROCR::prediction(comb[,3], comb[,1])
   auc <- ROCR::performance(pred,"auc")@y.values[[1]]
   rf_auc_subtype_species_validation_m[nrow(rf_auc_subtype_species_validation_m)+1,] <- c("validation",auc,m)
}

