library(dplyr)
library(Rmisc)
library(coin)
library(ggplot2)
library(tibble)
library(caret)
library(randomForest)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(pROC)
library(readxl)
library(vegan)

std.error <- function(x) sd(x)/sqrt(length(x))
window_db_sel_subtype_pval_1 <- read.csv("data/window_db_sel_subtype_pval_1.csv")
window_db_sel_subtype_pval_1$fdr <- p.adjust(window_db_sel_subtype_pval_1$pval,method = "fdr")
window_db_sel_subtype_pval_1_0.01 <- window_db_sel_subtype_pval_1%>%filter(fdr<0.01)

auc_cb_intra <- read.csv("data/auc_cb_intra.csv")
auc_cb_intra%>%group_by(approach,study,feature)%>%dplyr::summarise(mean=mean(auc),median=quantile(auc)[3],
                                                                   q25=quantile(auc)[2],q75=quantile(auc)[4],
                                                                   CI_lower=CI(auc)[3],CI_upper=CI(auc)[1])
auc_cb_intra%>%group_by(approach,feature)%>%dplyr::summarise(mean=mean(auc),median=quantile(auc)[3],
                                                                   q25=quantile(auc)[2],q75=quantile(auc)[4],
                                                                   CI_lower=CI(auc)[3],CI_upper=CI(auc)[1])
#comparison of auc between subtype, species, and subtype+species
auc_cb_intra$feature <- as.factor(auc_cb_intra$feature)
wilcox_test(auc ~ feature, data=auc_cb_intra%>%filter(approach=="intra")%>%filter(feature%in%c("subtype","species")))
wilcox_test(auc ~ feature, data=auc_cb_intra%>%filter(approach=="intra")%>%filter(feature%in%c("subtype","subtype_species")))
wilcox_test(auc ~ feature, data=auc_cb_intra%>%filter(approach=="intra")%>%filter(feature%in%c("species","subtype_species")))
wilcox_test(auc ~ feature, data=auc_cb_intra%>%filter(approach=="cb")%>%filter(feature%in%c("subtype","species")))
wilcox_test(auc ~ feature, data=auc_cb_intra%>%filter(approach=="cb")%>%filter(feature%in%c("subtype","subtype_species")))
wilcox_test(auc ~ feature, data=auc_cb_intra%>%filter(approach=="cb")%>%filter(feature%in%c("species","subtype_species")))
auc_cb_intra$approach <- factor(auc_cb_intra$approach,levels = c("intra","cb"))
auc_cb_intra$feature <- factor(auc_cb_intra$feature,levels = c("subtype","species","subtype_species"))
#Supplementary Figure 12a
auc_compa <- ggplot(auc_cb_intra, aes(x=feature,y=auc,color=feature)) + 
   geom_boxplot(linewidth=1)+
   facet_wrap(~approach,
              labeller = labeller(approach=c("intra"="Intra-cohort","cb"="Combined cohort")))+
   stat_summary(fun.y=mean, geom="point",size=2, color="red")+
   scale_color_manual(values=c("#D95F02", "#1B9E77", "#7570B3"))+
   scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
   theme_bw()+
   labs(x = "", y="AU-ROC")+
   theme(panel.border = element_rect(linewidth = 1)) +
   theme(panel.grid.minor = element_blank())+
   theme(legend.position="none",
         strip.text = element_text(color="black",size=15),
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_text(size=15, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_text(size = 16),
         axis.title.x = element_blank())

#importance feature
#subtype
window_db_subtype <- read.csv("data/window_db_subtype.csv",check.names=FALSE,stringsAsFactors = FALSE)
metadata_prenec <- read.csv("data/metadata_prenec_closest.csv")
window_db <- read.csv("data/window_db.csv")
window_db_047 <- window_db%>%filter(window%in%c(0,4,7))
window_db_subtype_rel <- window_db_subtype%>%column_to_rownames("SequenceID")
window_db_subtype_rel <-  as.data.frame(prop.table(as.matrix(window_db_subtype_rel), 1))
window_db_subtype_rel <- window_db_subtype_rel%>%rownames_to_column("SequenceID")
window_db_subtype_047 <- window_db_subtype_rel%>%filter(SequenceID%in%window_db_047$SequenceID)%>%column_to_rownames("SequenceID")
window_db_subtype_047 <- window_db_subtype_047[,which(colSums(window_db_subtype_047)!=0)]
window_db_subtype_047 <- window_db_subtype_047%>%rownames_to_column("SequenceID")
window_db_sel_subtype_meta <- right_join(metadata_prenec%>%dplyr::select(SequenceID,Study,nec),window_db_subtype_047,by="SequenceID")
window_db_sel_subtype_meta <- window_db_sel_subtype_meta%>%filter(Study%in%c("BrooksB_2017","MasiAC_2021","RahmanSF_2018","RavehSadkaT_2016"))%>%
   dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)

#intra cohort prediction
rf_auc_intra_subtype_047 <- data.frame("study"="x","time"="x","auc"=0,"fold"=0,"mtry"=0,"test"="x",stringsAsFactors = F)[-1,]
rf_imp_intra_subtype_047 <- data.frame()
for (s in unique(window_db_sel_subtype_meta$Study)) {
   # s="BrooksB_2017"
   print(s)
   for (c in seq(1,10,1)) {
      # c=1
      print(c)
      window_db_sel_subtype_meta_study <- window_db_sel_subtype_meta%>%filter(Study==s)%>%dplyr::select(1:3,which(colSums(.[-c(1:3)])!=0)+3)
      set.seed(c*100)
      dataset_index <- createFolds(window_db_sel_subtype_meta_study$nec, k = 5, list = TRUE, returnTrain = FALSE)
      for (f in paste0("Fold",seq(1,5,1))) {
         # f="Fold3"
         print(f)
         m <- auc_intra_pub%>%filter(feature=="subtype")%>%filter(study==s,time==c,fold==f)
         test_set <- window_db_sel_subtype_meta_study[dataset_index[[f]],]%>%dplyr::select(-Study)
         rownames(test_set) <- NULL
         train_set <- window_db_sel_subtype_meta_study%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
         x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
         y_train <- as.factor(train_set$nec)
         set.seed(1000)
         rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m$mtry)
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
         print(round(auc,3)==round(m$auc,3))
         if (round(auc,3)==round(m$auc,3)) {
            rf_auc_intra_subtype_047[nrow(rf_auc_intra_subtype_047)+1,] <- c(s,c,auc,f,m$mtry,"true")
         }else{
            rf_auc_intra_subtype_047[nrow(rf_auc_intra_subtype_047)+1,] <- c(s,c,auc,f,m$mtry,"false")
         }
         imp_fea <- as.data.frame(rf_train$importance)
         imp_fea <- imp_fea%>%rownames_to_column("subtype")%>%
            dplyr::select(subtype,MeanDecreaseAccuracy, MeanDecreaseGini)%>%
            dplyr::arrange(desc(MeanDecreaseGini))%>%
            mutate("study"=s)%>%mutate("time"=c)%>%mutate("fold"=f)%>%mutate("mtry"=m$mtry)
         rf_imp_intra_subtype_047 <- rbind(rf_imp_intra_subtype_047,imp_fea)
      }
   }
}
rf_imp_intra_subtype_047_sum_n <- rf_imp_intra_subtype_047%>%group_by(study,subtype)%>%dplyr::summarise(n=n())
rf_imp_intra_subtype_047_sum_prop <- rf_imp_intra_subtype_047%>%group_by(study,time,fold)%>%
   dplyr::mutate(prop=MeanDecreaseGini/sum(MeanDecreaseGini))%>%filter(subtype%in%window_db_sel_subtype_pval_1_0.01$subtype)%>%
   group_by(study,time,fold)%>%dplyr::summarise(total=sum(prop))%>%
   group_by(study)%>%dplyr::summarise(mean=mean(total),se=std.error(total))

rf_imp_intra_subtype_047_sum <- rf_imp_intra_subtype_047%>%group_by(study,subtype)%>%
   dplyr::summarise(mean=mean(MeanDecreaseGini))%>%
   dplyr::group_by(study)%>%dplyr::arrange(desc(mean),.by_group = TRUE)%>%dplyr::mutate(seq=rev(dense_rank(mean)))
which(!window_db_sel_subtype_pval_1_0.01$subtype%in%rf_imp_intra_subtype_047_sum$subtype)
rf_imp_intra_subtype_047_sum_sel <- rf_imp_intra_subtype_047_sum%>%filter(subtype%in%window_db_sel_subtype_pval_1_0.01$subtype)

#combined cohort
rf_auc_cb_subtype_047 <- data.frame("time"="x","auc"=0,"fold"=0,"mtry"=0,"test"="x",stringsAsFactors = F)[-1,]
rf_imp_cb_subtype_047 <- data.frame()
for (c in seq(1,10,1)) {
   # c=3
   print(c)
   set.seed(c*100)
   dataset_index <- createFolds(window_db_sel_subtype_meta$nec, k = 5, list = TRUE, returnTrain = FALSE)
   for (f in paste0("Fold",seq(1,5,1))) {
      # f="Fold1"
      print(f)
      m <- auc_cb_pub%>%filter(feature=="subtype")%>%filter(time==c,fold==f)
      test_set <- window_db_sel_subtype_meta[dataset_index[[f]],]%>%dplyr::select(-Study)
      rownames(test_set) <- NULL
      train_set <- window_db_sel_subtype_meta%>%filter(!SequenceID%in%test_set$SequenceID)%>%dplyr::select(-Study)
      x_train <- train_set%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
      y_train <- as.factor(train_set$nec)
      set.seed(1000)
      rf_train <- randomForest(x=x_train,y=y_train, importance = TRUE, ntree = 500, proximity = TRUE, mtry = m$mtry)
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
      print(round(auc,3)==round(m$auc,3))
      if (round(auc,3)==round(m$auc,3)) {
         rf_auc_cb_subtype_047[nrow(rf_auc_cb_subtype_047)+1,] <- c(c,auc,f,m$mtry,"true")
      }else{
         rf_auc_cb_subtype_047[nrow(rf_auc_cb_subtype_047)+1,] <- c(c,auc,f,m$mtry,"false")
      }
      imp_fea <- as.data.frame(rf_train$importance)
      imp_fea <- imp_fea%>%rownames_to_column("subtype")%>%
         dplyr::select(subtype,MeanDecreaseAccuracy, MeanDecreaseGini)%>%
         dplyr::arrange(desc(MeanDecreaseGini))%>%
         mutate("time"=c)%>%mutate("fold"=f)%>%mutate("mtry"=m$mtry)
      rf_imp_cb_subtype_047 <- rbind(rf_imp_cb_subtype_047,imp_fea)
   }
}
rf_imp_cb_subtype_047_sum_n <- rf_imp_cb_subtype_047%>%group_by(subtype)%>%dplyr::summarise(n=n())
rf_imp_cb_subtype_047_sum_prop <- rf_imp_cb_subtype_047%>%group_by(time,fold)%>%
   dplyr::mutate(prop=MeanDecreaseGini/sum(MeanDecreaseGini))%>%filter(subtype%in%window_db_sel_subtype_pval_1_0.01$subtype)%>%
   group_by(time,fold)%>%dplyr::summarise(total=sum(prop))%>%ungroup()%>%
   dplyr::summarise(mean=mean(total),se=std.error(total))%>%mutate(study="all")

rf_imp_cb_subtype_047_sum <- rf_imp_cb_subtype_047%>%group_by(subtype)%>%
   dplyr::summarise(mean=mean(MeanDecreaseGini))%>%
   dplyr::arrange(desc(mean),.by_group = TRUE)%>%dplyr::mutate(seq=rev(dense_rank(mean)))
rf_imp_cb_subtype_047_sum_sel <- rf_imp_cb_subtype_047_sum%>%filter(subtype%in%window_db_sel_subtype_pval_1_0.01$subtype)%>%mutate(study="all")

subtype_imp <- rbind(rf_imp_intra_subtype_047_sum_sel,rf_imp_cb_subtype_047_sum_sel)
subtype_imp$study[subtype_imp$study=="all"] <- "Combined cohorts"
subtype_imp_hp <- subtype_imp%>%dplyr::select(-seq)
names(subtype_imp_hp)[3] <- "MeanDecreaseGini"
subtype_imp_hp$MeanDecreaseGini[subtype_imp_hp$MeanDecreaseGini>1] <- 1
subtype_imp_hp_s <- subtype_imp_hp%>%spread(.,subtype,MeanDecreaseGini)%>%column_to_rownames("study")
rownames(subtype_imp_hp_s)
subtype_imp_hp_s_2 <- subtype_imp_hp_s[c(1,3,4,5,2),]
range(subtype_imp_hp_s_2,na.rm = T)#0.0004612298 2.9231589015
names(subtype_imp_hp_s_2) <- gsub("subtype_","",names(subtype_imp_hp_s_2),fixed = T)
names(subtype_imp_hp_s_2) <- gsub("other_peptide_antibiotics","OPA",names(subtype_imp_hp_s_2),fixed = T)
col_fun_type <- colorRamp2(c(0,0.0625,0.125,0.25,0.5,0.75,1),c("#0E6197","#5199C8","#8CC0E1","#EFBBB4","#D88B80","#CD7467","#B6594C"))
#Figure 5e
subtype_rf_imp <- Heatmap(as.matrix(subtype_imp_hp_s_2),cluster_rows=F,cluster_columns=T,col = col_fun_type,
                          na_col = "white",
                          show_row_names = T,show_column_names = T,row_names_side = "left",
                          row_split = c(rep(1,4),2),
                          border = T,row_names_gp = gpar(fontsize =10),column_names_gp = gpar(fontsize =10),
                          rect_gp = gpar(col = "grey", lwd = 0.7),
                          border_gp=gpar(col = "black", lwd = 1),column_gap = unit(2, "mm"),
                          width = unit(20, "cm"), height = unit(5, "cm"),
                          show_heatmap_legend = T,
                          heatmap_legend_param = list(title = "MeanDecreaseGini",direction = "horizontal", 
                                                      legend_hight = unit(1, "cm"),legend_width = unit(3, "cm"),
                                                      at = c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.5","0.75",">1"))
)

#auc and proportion of MeanDecreaseGini for intra and combined cohorts
auc_cb_intra_sum <- auc_cb_intra%>%group_by(approach,study,feature)%>%
   dplyr::summarise(mean=mean(auc),median=quantile(auc)[3],q25=quantile(auc)[2],q75=quantile(auc)[4],se=std.error(auc))
auc_cb_intra_sum_subtype <- auc_cb_intra_sum%>%filter(feature=="subtype")%>%ungroup()%>%dplyr::select(study,mean,se)%>%mutate(cata="auc")
prop_cb_intra <- rbind(rf_imp_intra_subtype_047_sum_prop%>%mutate(cata="prop"),rf_imp_cb_subtype_047_sum_prop%>%mutate(cata="prop",study="all"))
auc_cb_intra_sum_subtype_prop <- rbind(auc_cb_intra_sum_subtype,prop_cb_intra)
auc_cb_intra_sum_subtype_prop$study <- factor(auc_cb_intra_sum_subtype_prop$study,levels = rev(c("BrooksB_2017","MasiAC_2021","RahmanSF_2018","RavehSadkaT_2016","all")))
#Figure 5e
subtype_auc <- ggplot(auc_cb_intra_sum_subtype_prop, aes(x=mean, y=study, fill=cata)) + 
   geom_bar(stat="identity", width = 0.7, color="black", position=position_dodge(0.7)) +
   geom_errorbar(aes(xmin=mean, xmax=mean+se), width=0.5,
                 position=position_dodge(0.7))+
   geom_vline(xintercept = 0.9,linetype="dashed",size=0.5,color="#C1C4C3")+
   scale_fill_manual(values = c("#4994C4","#E04B37"),labels=c("auc"="AU-ROC","prop"="Relative weight"),guide = guide_legend(reverse = TRUE))+
   scale_x_continuous(expand = c(0,0),breaks = c(0,0.25,0.5,0.75,0.9),limits = c(0,1.1))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(color = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=12, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=10, colour = "black",angle = 90,vjust = 0.5),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x = element_blank(),
         plot.margin = unit(c(1,5,1,5),"mm"))

########################################################################################################################################
#validation
########################################################################################################################################
metadata_prenec <- read_excel("data/metadata_NEC_validation.xlsx")
metadata_prenec_pre <- metadata_prenec%>%filter(PreNEC!="postNEC")
metadata_prenec$NEC_onset <- as.numeric(metadata_prenec$NEC_onset)
metadata_prenec_yes <- metadata_prenec%>%filter(nec=="yes")%>%mutate(preday=NEC_onset-DOL)
metadata_prenec_yes_sum <- metadata_prenec_yes%>%group_by(preday,SubjectID)%>%dplyr::summarise(n=n())

window_db <- data.frame()
window_list <-list(c(0,6),c(7,9),c(10,12),c(13,15),c(16,21),c(22,30),c(31,60))
for (w in 1:length(window_list)) {
   # w=1
   print(window_list[[w]])
   pre_yes <- metadata_prenec_yes%>%filter(preday>=window_list[[w]][1] & preday<=window_list[[w]][2])
   if (nrow(pre_yes)==0) {
      next
   }else{
      pre_yes_1 <- pre_yes%>%
         mutate(window=window_list[[w]][1])
      window_db_no <- data.frame()
      for (s in unique(pre_yes$SubjectID)) {
         print(s)
         a <- pre_yes%>%filter(SubjectID==s)%>%top_n(-1,preday)%>%mutate(window=window_list[[w]][1])%>%distinct(DOL,.keep_all = T)
         length(unique(a$SubjectID))
         b <- metadata_prenec_no
         b_1 <- b%>%filter(Gestational_age>=a$Gestational_age-0.5 & Gestational_age<=a$Gestational_age+0.5)%>%
            filter(DOL>=a$DOL-1 & DOL<=a$DOL+1)%>%
            group_by(SubjectID)%>%dplyr::top_n(1,DOL)%>%
            mutate(preday=NA)%>%mutate(window=window_list[[w]][1])
         print(nrow(b_1))
         window_db_no <- rbind(window_db_no,b_1)
      }
      window_db_no_unique <- window_db_no%>%distinct(SampleID,.keep_all = T)#some samples are selected by >1 times
      pre_yes_2 <- rbind(pre_yes_1,window_db_no_unique)
      window_db <- rbind(window_db,pre_yes_2)
   }
}
window_db_yes <- window_db%>%filter(!is.na(preday))
window_db_yes_sum <- window_db_yes%>%group_by(window)%>%dplyr::summarise(n=n())
window_db_no <- window_db%>%filter(is.na(preday))
window_db_no_sum <- window_db_no%>%group_by(window)%>%dplyr::summarise(n=n())

#subtype
preterm_resistome_prenec_subtype <- read.delim("data/normalized_cell.subtype.txt")
names(preterm_resistome_prenec_subtype) <- gsub(".kneaddata_paired","",names(preterm_resistome_prenec_subtype),fixed = T)
preterm_resistome_prenec_subtype <- preterm_resistome_prenec_subtype%>%dplyr::select(subtype,metadata_prenec$SampleID)%>%column_to_rownames("subtype")
preterm_resistome_prenec_subtype_t <- preterm_resistome_prenec_subtype[which(rowSums(preterm_resistome_prenec_subtype)!=0),]%>%
   t()%>%as.data.frame()
names(preterm_resistome_prenec_subtype_t) <- paste0("subtype_",names(preterm_resistome_prenec_subtype_t))
subtype_beta <- vegdist(preterm_resistome_prenec_subtype_t, method = "bray")
subtype_beta_2 <- subtype_beta%>%as.matrix()%>%as.data.frame()
subtype_beta_2[upper.tri(subtype_beta_2, diag=T)] <- NA
subtype_beta_2_g <- subtype_beta_2%>%rownames_to_column("run1")%>%gather(.,run2,distance,-run1)%>%filter(!is.na(distance))
subtype_beta_2_g <- left_join(subtype_beta_2_g,metadata_prenec%>%dplyr::select(SampleID,SubjectID,DOL),by=c("run1"="SampleID"))
subtype_beta_2_g <- left_join(subtype_beta_2_g,metadata_prenec%>%dplyr::select(SampleID,SubjectID,DOL),by=c("run2"="SampleID"))

table(window_db_yes$window)
nec_beta_yes_no <- data.frame()
nec_beta_yes_no_pval <- data.frame("preday"=0,"median_yes"=0,"mean_yes"=0,"median_no"=0,"mean_no"=0,"pval"=0,stringsAsFactors = F)[-1,]
for (i in c(0,10,13,16,22,31)) {
   # i=0
   print(i)
   window_db_yes_a <- window_db%>%filter(nec=="yes")%>%filter(window==i)
   dist_yes <- subtype_beta_2_g%>%filter(run1%in%window_db_yes_a$SampleID)%>%filter(run2%in%window_db_yes_a$SampleID)%>%
      mutate(nec="yes")%>%mutate(window=i)
   print(which(dist_yes$SubjectID.x==dist_yes$SubjectID.y))
   window_db_no_a <- window_db%>%filter(nec=="no")%>%filter(window==i)
   dist_no <- subtype_beta_2_g%>%filter(run1%in%window_db_no_a$SampleID)%>%filter(run2%in%window_db_no_a$SampleID)%>%
      mutate(nec="no")%>%mutate(window=i)
   print(which(dist_no$SubjectID.x==dist_no$SubjectID.y))
   dist_yes_no <- rbind(dist_no,dist_yes)
   ptest <- wilcox.test(dist_yes_no$distance ~ dist_yes_no$nec)
   pval <- ptest$p.value
   nec_beta_yes_no_pval[nrow(nec_beta_yes_no_pval)+1,] <- c(i,median(dist_yes$distance),mean(dist_yes$distance),
                                                            median(dist_no$distance),mean(dist_no$distance),
                                                            pval)
   nec_beta_yes_no <- rbind(nec_beta_yes_no,dist_yes_no)
}
nec_beta_yes_no_pval
nec_beta_yes_no$nec <- factor(nec_beta_yes_no$nec,levels = c("yes","no"))
nec_beta_yes_no <- nec_beta_yes_no%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Figure 6a
nec_preday <- ggplot(data=nec_beta_yes_no, aes(x=window_2, y=distance, color=nec)) +
   geom_jitter(alpha=0.2,size=1,width = 0.2)+
   geom_smooth(method=loess,formula = 'y ~ x', level=0.95,alpha = 0.5,aes(fill = nec))+
   scale_color_manual(values = c("#E04B37","#4994C4"),labels=c("yes"="NEC","no"="No NEC"))+
   labs(x ="Days before NEC onset", y="Bray-Curtis dissimilarity\n(based on ARG subtypes)")+
   theme_bw()+
   scale_y_continuous(limits = c(0,1.2),breaks = c(0,0.3,0.6,0.9,1.2))+
   scale_x_reverse(breaks = unique(nec_beta_yes_no$window_2),
                   labels=c("6-0","12-10","15-13","21-16","30-22","60-31"))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(legend.key = element_rect(fill = "transparent",colour = NA))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90, vjust = 0.5,hjust=1),
         axis.title.x = element_text(size=12, colour = "black"))

##############################################################################################################################
#random forest of validation
##############################################################################################################################
vali_auc <- read.csv("data/vali_auc.csv")
#subtype
window_db_subtype <- read.csv("data/window_db_subtype.csv",check.names=FALSE,stringsAsFactors = FALSE)
window_db_subtype_validation <- read.csv("data/window_db_sel_subtype_meta_validation.csv",check.names=FALSE,stringsAsFactors = FALSE)
window_db_subtype_validation <- window_db_subtype_validation[c(1:3,which(colSums(window_db_subtype_validation[-c(1:3)])!=0)+3)]

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

#random forest
train_set_subtype <- window_db_sel_subtype_meta %>% filter(Study!="validation")%>%dplyr::select(-Study)
x_train_subtype <- train_set_subtype%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train_subtype <- as.factor(train_set_subtype$nec)
set.seed(1000)
rf_train_subtype <- randomForest(x=x_train_subtype,y=y_train_subtype, importance = TRUE, ntree = 500, proximity = TRUE,mtry=vali_auc$mtry[vali_auc$feature=="subtype"])
test_set_subtype <- window_db_sel_subtype_meta %>% filter(Study=="validation")%>%dplyr::select(-Study)
x_test_subtype <- test_set_subtype%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_test_subtype <- test_set_subtype$nec
rf_pred_prob_subtype <- predict(rf_train_subtype, newdata = x_test_subtype, type = "prob")
imp_fea_subtype <- as.data.frame(rf_train_subtype$importance)
imp_fea_subtype <- imp_fea_subtype%>%rownames_to_column("subtype")%>%
   dplyr::select(subtype,MeanDecreaseAccuracy, MeanDecreaseGini)%>%dplyr::arrange(desc(MeanDecreaseGini))

#species
train_set_species <- window_db_sel_species_meta %>% filter(Study!="validation")%>%dplyr::select(-Study)
x_train_species <- train_set_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train_species <- as.factor(train_set_species$nec)
set.seed(1000)
rf_train_species <- randomForest(x=x_train_species,y=y_train_species, importance = TRUE, ntree = 500, proximity = TRUE,mtry=vali_auc$mtry[vali_auc$feature=="species"])
test_set_species <- window_db_sel_species_meta %>% filter(Study=="validation")%>%dplyr::select(-Study)
x_test_species <- test_set_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_test_species <- test_set_species$nec
rf_pred_prob_species <- predict(rf_train_species, newdata = x_test_species, type = "prob")
imp_fea_species <- as.data.frame(rf_train_species$importance)
imp_fea_species <- imp_fea_species%>%rownames_to_column("species")%>%
   dplyr::select(species,MeanDecreaseAccuracy, MeanDecreaseGini)%>%dplyr::arrange(desc(MeanDecreaseGini))

#subtype+species
train_set_subtype_species <- window_db_sel_subtype_species_meta %>% filter(Study!="validation")%>%dplyr::select(-Study)
x_train_subtype_species <- train_set_subtype_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train_subtype_species <- as.factor(train_set_subtype_species$nec)
set.seed(1000)
rf_train_subtype_species <- randomForest(x=x_train_subtype_species,y=y_train_subtype_species, importance = TRUE, ntree = 500, proximity = TRUE,mtry=vali_auc$mtry[vali_auc$feature=="subtype_species"])
test_set_subtype_species <- window_db_sel_subtype_species_meta %>% filter(Study=="validation")%>%dplyr::select(-Study)
x_test_subtype_species <- test_set_subtype_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_test_subtype_species <- test_set_subtype_species$nec
rf_pred_prob_subtype_species <- predict(rf_train_subtype_species, newdata = x_test_subtype_species, type = "prob")
imp_fea_subtype_species <- as.data.frame(rf_train_subtype_species$importance)
imp_fea_subtype_species <- imp_fea_subtype_species%>%rownames_to_column("subtype_species")%>%
   dplyr::select(subtype_species,MeanDecreaseAccuracy, MeanDecreaseGini)%>%dplyr::arrange(desc(MeanDecreaseGini))

#combined three type of featurs
rf_pred_prob_subtype_df <- rf_pred_prob_subtype%>%as.data.frame()%>%dplyr::select(yes)%>%setNames("subtype")%>%rownames_to_column("SequenceID")
rf_pred_prob_species_df <- rf_pred_prob_species%>%as.data.frame()%>%dplyr::select(yes)%>%setNames("species")%>%rownames_to_column("SequenceID")
rf_pred_prob_subtype_species_df <- rf_pred_prob_subtype_species%>%as.data.frame()%>%dplyr::select(yes)%>%setNames("subtype_species")%>%rownames_to_column("SequenceID")
rf_pred_prob_cb <- left_join(window_db_sel_subtype_meta%>%filter(Study=="validation")%>%dplyr::select(SequenceID,nec),rf_pred_prob_subtype_df,by="SequenceID")
rf_pred_prob_cb <- left_join(rf_pred_prob_cb,rf_pred_prob_species_df,by="SequenceID")
rf_pred_prob_cb <- left_join(rf_pred_prob_cb,rf_pred_prob_subtype_species_df,by="SequenceID")
rocobj_list <- roc(nec ~ subtype+species+subtype_species,ci=T,smooth=F,data=rf_pred_prob_cb)
rocobj_list$subtype$auc
rocobj_list$subtype$ci
rocobj_list$species$auc
rocobj_list$species$ci
rocobj_list$subtype_species$auc
rocobj_list$subtype_species$ci

ci_list <- lapply(rocobj_list, ci.se, specificities = seq(0, 1, 0.01))
ci_list_df <- lapply(ci_list, function(ciobj)
   data.frame(prop = as.numeric(rownames(ciobj)),
              lower = ciobj[, 1],
              upper = ciobj[, 3]))
#Figure 6b
plot_curve <- ggroc(rocobj_list,size=0.5,legacy.axes =T)+
   scale_colour_manual(values = c("#0B74BA", "#42A90C", "#CC2B0A"))+
   geom_segment(aes(x=0, y=0, xend=1, yend=1),colour="grey",linetype = "dashed")+
   # geom_ribbon(data=data_ci_cb,aes(x=1-prop,ymin=`2.5%`,ymax=`97.5%`,fill=type),alpha=0.5)+
   theme_bw()+
   labs(x ="False positive rate", y="True positive rate")+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.border = element_rect(colour = "black",size = 1),
         panel.grid = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.x = element_text(size=12, colour = "black"))
for(i in 1:3) {
   plot_curve <- plot_curve + geom_ribbon(
      data = ci_list_df[[i]],
      aes(x = 1-prop, ymin = lower, ymax = upper),
      fill = i + 1,
      alpha = 0.1,
      inherit.aes = F) 
} 

#feature importance: subtype
imp_fea_subtype_sel <- imp_fea_subtype%>%mutate(rank=rownames(.))%>%filter(subtype%in%window_db_sel_subtype_pval_1_0.01$subtype)
imp_fea_subtype_sel$subtype <- gsub("subtype_","",imp_fea_subtype_sel$subtype,fixed = T)
imp_fea_subtype_sel$subtype <- gsub("other_peptide_antibiotics","OPA",imp_fea_subtype_sel$subtype,fixed = T)
imp_fea_subtype_sel$subtype <- factor(imp_fea_subtype_sel$subtype,levels = imp_fea_subtype_sel$subtype)
#Figure 6c
subtype_imp_plot <- ggplot(imp_fea_subtype_sel, aes(x=subtype, y=MeanDecreaseGini)) + 
   geom_bar(stat = "identity",fill="#E09137",color="black",width = 0.8)+
   geom_text(aes(label=rank,y=MeanDecreaseGini+0.06),size=3)+
   theme_bw() + # remove the backgroud
   labs(x = element_blank(), y="MeanDecreaseGini")+
   scale_y_continuous(expand = c(0,0),limits = c(0,1.5),breaks = c(0,0.5,1.0,1.5),labels = c("0","0.5","1.0","1.5"))+
   theme(panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth = 0.5)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   theme(aspect.ratio=1/3.8)+
   theme(legend.position="none",
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.title.y =  element_text(size=12, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour = "black",angle = 90,hjust = 1,vjust=0.5),
         plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

