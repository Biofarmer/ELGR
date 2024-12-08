library(readxl)
library(dplyr)
library(tibble)
library(vegan)
library(tidyr)
library(ggplot2)
library(lmerTest)
library(plyr)
########################################################################################################################################
#convergence of validation cohort
########################################################################################################################################
metadata_prenec <- read_excel("data/metadata_validation.xlsx")
metadata_prenec_pre <- metadata_prenec%>%filter(PreNEC!="postNEC")
metadata_prenec$NEC_onset <- as.numeric(metadata_prenec$NEC_onset)
metadata_prenec_yes <- metadata_prenec%>%filter(nec=="yes")%>%mutate(preday=NEC_onset-DOL)
metadata_prenec_no <- metadata_prenec%>%filter(nec=="no")
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
         a <- pre_yes%>%filter(SubjectID==s)%>%top_n(-1,preday)%>%mutate(window=window_list[[w]][1])%>%distinct(DOL,.keep_all = T)
         b <- metadata_prenec_no
         b_1 <- b%>%filter(Gestational_age>=a$Gestational_age-0.5 & Gestational_age<=a$Gestational_age+0.5)%>%
            filter(DOL>=a$DOL-1 & DOL<=a$DOL+1)%>%
            group_by(SubjectID)%>%dplyr::top_n(1,DOL)%>%
            mutate(preday=NA)%>%mutate(window=window_list[[w]][1])%>%
            mutate(case=s)
         window_db_no <- rbind(window_db_no,b_1)
      }
      window_db_no_dup <- window_db_no%>%distinct(SampleID,.keep_all = T)
      print(paste0("case number: ", length(unique(window_db_no_dup$case))))
      print(paste0("control number: ", length(unique(window_db_no_dup$SubjectID))))
      #to make sure each case has only one matched sample from each control infant.
      window_db_no_dup_1 <- window_db_no_dup
      window_db_no_unique_cb <- data.frame()
      window_db_no_sum_1 <- window_db_no_dup[1,]
      window_db_no_sum_2 <- window_db_no_dup[1,]
      while (nrow(window_db_no_sum_1)!=0 & nrow(window_db_no_sum_2)!=0) {
         window_db_no_sum <- window_db_no_dup_1%>%group_by(SubjectID)%>%dplyr::summarise(n=n())
         window_db_no_sum_1 <- window_db_no_sum%>%filter(n==1)
         window_db_no_sum_2 <- window_db_no_sum%>%filter(n>1)
         window_db_no_unique_1 <- window_db_no_dup_1%>%filter(SubjectID%in%window_db_no_sum_1$SubjectID)
         window_db_no_unique_cb <- rbind(window_db_no_unique_cb,window_db_no_unique_1)
         window_db_no_dup_1 <- window_db_no_dup_1%>%filter(SubjectID%in%window_db_no_sum_2$SubjectID)%>%filter(!case%in%window_db_no_unique_1$case)
      }
      print(paste0("window_db_no_sum_2: ",nrow(window_db_no_sum_2)))
      print(paste0("case number: ", length(unique(window_db_no_unique_cb$case))))
      #after run code above, some control infants left and then pick only one sample from each left control infants.
      window_db_no_unique_left <- window_db_no_dup%>%filter(!SubjectID%in%window_db_no_unique_cb$SubjectID)%>%distinct(SubjectID,.keep_all = T)
      window_db_no_unique <- rbind(window_db_no_unique_cb,window_db_no_unique_left)
      print(paste0("control number: ", length(unique(window_db_no_unique$SubjectID))))
      # window_db_no_unique2 <- window_db_no%>%distinct(SubjectID,.keep_all = T)#some samples are selected by >1 times
      pre_yes_2 <- rbind.fill(pre_yes_1,window_db_no_unique)
      window_db <- rbind.fill(window_db,pre_yes_2)
   }
}
window_db_yes <- window_db%>%filter(!is.na(preday))
window_db_yes_sum <- window_db_yes%>%group_by(window)%>%dplyr::summarise(n=n())
window_db_no <- window_db%>%filter(is.na(preday))
window_db_no_sum <- window_db_no%>%group_by(window)%>%dplyr::summarise(n=n())

#subtype
preterm_resistome_prenec_subtype <- read.csv("data/normalized_cell.subtype_sel.csv",check.names = F)
preterm_resistome_prenec_subtype_t <- preterm_resistome_prenec_subtype%>%column_to_rownames("sample")
names(preterm_resistome_prenec_subtype_t) <- paste0("subtype_",names(preterm_resistome_prenec_subtype_t))
subtype_beta <- vegdist(preterm_resistome_prenec_subtype_t, method = "bray")
subtype_beta_2 <- subtype_beta%>%as.matrix()%>%as.data.frame()
subtype_beta_2[upper.tri(subtype_beta_2, diag=T)] <- NA
subtype_beta_2_g <- subtype_beta_2%>%rownames_to_column("run1")%>%gather(.,run2,distance,-run1)%>%filter(!is.na(distance))
subtype_beta_2_g <- left_join(subtype_beta_2_g,metadata_prenec%>%dplyr::select(SampleID,SubjectID,DOL),by=c("run1"="SampleID"))
subtype_beta_2_g <- left_join(subtype_beta_2_g,metadata_prenec%>%dplyr::select(SampleID,SubjectID,DOL),by=c("run2"="SampleID"))
#beta diversity
nec_beta_yes_no <- data.frame()
nec_beta_yes_no_pval <- data.frame("preday"=0,"median_yes"=0,"mean_yes"=0,"median_no"=0,"mean_no"=0,"pval"=0,stringsAsFactors = F)[-1,]
for (i in c(0,10,13,16,22,31)) {
   # i=16
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
nec_beta_yes_no_pval$fdr <- p.adjust(nec_beta_yes_no_pval$pval,method = "fdr")
nec_beta_yes_no$nec <- factor(nec_beta_yes_no$nec,levels = c("yes","no"))
nec_beta_yes_no <- nec_beta_yes_no%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Figure 5d
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
########################################################################################################################################
#prediction
########################################################################################################################################
vali_subtype <- read.csv("data/rf_auc_subtype_validation_m_noadd_six_newrule.csv")
vali_subtype_sel <- vali_subtype%>%filter(auc==max(auc))%>%filter(mtry==min(mtry))
vali_species <- read.csv("data/rf_auc_species_validation_m_noadd_six_newrule.csv")
vali_species_sel <- vali_species%>%filter(auc==max(auc))%>%filter(mtry==min(mtry))
vali_subtype_species <- read.csv("data/rf_auc_subtype_species_validation_m_noadd_six_newrule.csv")
vali_subtype_species_sel <- vali_subtype_species%>%filter(auc==max(auc))%>%filter(mtry==min(mtry))

#subtype
window_db_subtype <- read.csv("data/window_db_subtype_new.csv",check.names=FALSE,stringsAsFactors = FALSE)#1382 runs
#species
metaphlan4_preterm <- read.csv("data/metaphlan4_preterm_new.csv")
metaphlan4_preterm[is.na(metaphlan4_preterm)] <- 0
metaphlan4_preterm_2 <- metaphlan4_preterm%>%dplyr::select(clade_name,unique(window_db_subtype$SequenceID))
metaphlan4_species <- metaphlan4_preterm_2%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
metaphlan4_species$species <- sapply(strsplit(metaphlan4_species$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
metaphlan4_species <- metaphlan4_species%>%dplyr::select(-clade_name)%>%column_to_rownames("species")
metaphlan4_species <- metaphlan4_species[which(colSums(metaphlan4_species)!=0)]
#validation
window_db_subtype_validation <- read.csv("data/window_db_sel_subtype_meta_validation_newrule.csv",check.names=FALSE,stringsAsFactors = FALSE)#1382 runs
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
#metadata
metadata_prenec <- read.csv("data/metadata_prenec_closest_new.csv")
#windows
window_db <- read.csv("data/window_db_new.csv")#1434 runs
length(unique(window_db$SequenceID))#1229
#build model for windows 0,4,7 only with subtype
window_db_047 <- window_db%>%filter(window%in%c(0,4,7))%>%
   filter(Study%in%c("BrooksB_2017","LouYC_2024","MasiAC_2021","RahmanSF_2018","RavehSadkaT_2016","ThanertR_2024"))
#subtype
window_db_subtype_rel <- window_db_subtype%>%column_to_rownames("SequenceID")
window_db_subtype_rel <-  as.data.frame(prop.table(as.matrix(window_db_subtype_rel), 1))
window_db_subtype_rel <- window_db_subtype_rel%>%rownames_to_column("SequenceID")
window_db_subtype_047 <- window_db_subtype_rel%>%filter(SequenceID%in%window_db_047$SequenceID)%>%column_to_rownames("SequenceID")
window_db_subtype_047 <- window_db_subtype_047[,which(colSums(window_db_subtype_047)!=0)]
metadata_prenec_047_yes <- metadata_prenec%>%filter(SequenceID%in%window_db_047$SequenceID)%>%filter(nec=="yes")
window_db_subtype_047_yes <- window_db_subtype_047[metadata_prenec_047_yes$SequenceID,]
window_db_subtype_047_add <- rbind(window_db_subtype_047, 
                                   subtype_median = apply(window_db_subtype_047_yes,MARGIN=2,median,na.rm = TRUE),
                                   subtype_mean = apply(window_db_subtype_047_yes,MARGIN=2,mean,na.rm = TRUE))
window_db_sel_subtype_meta <- right_join(metadata_prenec%>%dplyr::select(SequenceID,Study,nec),
                                         window_db_subtype_047_add%>%rownames_to_column("SequenceID"),by="SequenceID")
window_db_sel_subtype_meta$Study[is.na(window_db_sel_subtype_meta$Study)] <- "cb" 
window_db_sel_subtype_meta <- rbind.fill(window_db_sel_subtype_meta,window_db_subtype_validation)
window_db_sel_subtype_meta[is.na(window_db_sel_subtype_meta)] <- 0
window_db_sel_subtype_meta_dist <- window_db_sel_subtype_meta%>%dplyr::select(-Study,-nec)%>%column_to_rownames("SequenceID")
window_db_sel_subtype_meta_dist_richness <- rowSums(window_db_sel_subtype_meta_dist != 0)%>%as.data.frame()%>%setNames("subtype_Richness")%>%rownames_to_column("SequenceID")
window_db_sel_subtype_meta_dist_shannon <- diversity(window_db_sel_subtype_meta_dist, index = "shannon", MARGIN=1)%>%as.data.frame()%>%setNames("subtype_Shannon")%>%rownames_to_column("SequenceID")
subtype_alpha <- left_join(window_db_sel_subtype_meta_dist_richness,window_db_sel_subtype_meta_dist_shannon,by="SequenceID")
window_db_sel_subtype_meta_dist_beta <- vegdist(window_db_sel_subtype_meta_dist, method = "bray")
window_db_sel_subtype_meta_dist_beta_2 <- window_db_sel_subtype_meta_dist_beta%>%as.matrix()%>%as.data.frame()
window_db_sel_subtype_meta_dist_beta_3 <- window_db_sel_subtype_meta_dist_beta_2[c("subtype_median","subtype_mean")]
window_db_sel_subtype_meta_feature <- left_join(window_db_sel_subtype_meta,
                                                window_db_sel_subtype_meta_dist_beta_3%>%rownames_to_column("SequenceID"),by="SequenceID")
window_db_sel_subtype_meta_feature <- left_join(window_db_sel_subtype_meta_feature,subtype_alpha,by="SequenceID")
window_db_sel_subtype_meta_feature <- window_db_sel_subtype_meta_feature%>%filter(!SequenceID%in%c("subtype_median","subtype_mean"))
window_db_sel_subtype_meta_feature <- window_db_sel_subtype_meta_feature%>%dplyr::select(-subtype_mean)

#species
metaphlan4_species_t <- metaphlan4_species%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
metaphlan4_species_t[-1] <- metaphlan4_species_t[-1]/100
window_db_species_047 <- metaphlan4_species_t%>%filter(SequenceID%in%window_db_047$SequenceID)%>%column_to_rownames("SequenceID")
window_db_species_047 <- window_db_species_047[,which(colSums(window_db_species_047)!=0)]
window_db_species_047_yes <- window_db_species_047[metadata_prenec_047_yes$SequenceID,]
window_db_species_047_add <- rbind(window_db_species_047, 
                                   species_median = apply(window_db_species_047_yes,MARGIN=2,median,na.rm = TRUE),
                                   species_mean = apply(window_db_species_047_yes,MARGIN=2,mean,na.rm = TRUE))
window_db_sel_species_meta <- right_join(metadata_prenec%>%dplyr::select(SequenceID,Study,nec),
                                         window_db_species_047_add%>%rownames_to_column("SequenceID"),by="SequenceID")
window_db_sel_species_meta <- rbind.fill(window_db_sel_species_meta,metaphlan4_species_validation_rel_t_meta)
window_db_sel_species_meta[is.na(window_db_sel_species_meta)] <- 0
window_db_sel_species_meta_dist <- window_db_sel_species_meta%>%dplyr::select(-Study,-nec)%>%column_to_rownames("SequenceID")
window_db_sel_species_meta_dist_richness <- rowSums(window_db_sel_species_meta_dist != 0)%>%as.data.frame()%>%setNames("species_Richness")%>%rownames_to_column("SequenceID")
window_db_sel_species_meta_dist_shannon <- diversity(window_db_sel_species_meta_dist, index = "shannon", MARGIN=1)%>%as.data.frame()%>%setNames("species_Shannon")%>%rownames_to_column("SequenceID")
species_alpha <- left_join(window_db_sel_species_meta_dist_richness,window_db_sel_species_meta_dist_shannon,by="SequenceID")
window_db_sel_species_meta_dist_beta <- vegdist(window_db_sel_species_meta_dist, method = "bray")
window_db_sel_species_meta_dist_beta_2 <- window_db_sel_species_meta_dist_beta%>%as.matrix()%>%as.data.frame()
window_db_sel_species_meta_dist_beta_3 <- window_db_sel_species_meta_dist_beta_2[c("species_median","species_mean")]
window_db_sel_species_meta_feature <- left_join(window_db_sel_species_meta,
                                                window_db_sel_species_meta_dist_beta_3%>%rownames_to_column("SequenceID"),by="SequenceID")
window_db_sel_species_meta_feature <- left_join(window_db_sel_species_meta_feature,species_alpha,by="SequenceID")
window_db_sel_species_meta_feature <- window_db_sel_species_meta_feature%>%filter(!SequenceID%in%c("species_median","species_mean"))
window_db_sel_species_meta_feature <- window_db_sel_species_meta_feature%>%dplyr::select(-species_mean)

#subtype + species
window_db_sel_subtype_species_meta <- left_join(window_db_sel_subtype_meta_feature,window_db_sel_species_meta_feature%>%dplyr::select(-c(2:3)),by="SequenceID")
window_db_sel_subtype_species_meta[is.na(window_db_sel_subtype_species_meta)] <- 0

#random forest
train_set_subtype <- window_db_sel_subtype_meta_feature %>% filter(Study!="validation")%>%dplyr::select(-Study,-subtype_median,-subtype_Richness,-subtype_Shannon)
x_train_subtype <- train_set_subtype%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train_subtype <- as.factor(train_set_subtype$nec)
set.seed(1000)
rf_train_subtype <- randomForest(x=x_train_subtype,y=y_train_subtype, importance = TRUE, ntree = 500, proximity = TRUE,mtry = vali_subtype_sel$mtry)
test_set_subtype <- window_db_sel_subtype_meta_feature%>%filter(Study=="validation")%>%
   dplyr::select(-Study,-subtype_median,-subtype_Richness,-subtype_Shannon)
x_test_subtype <- test_set_subtype%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_test_subtype <- test_set_subtype$nec
rf_pred_prob_subtype <- predict(rf_train_subtype, newdata = x_test_subtype, type = "prob")
rf_pred_subtype <- predict(rf_train_subtype, newdata = x_test_subtype)%>%as.data.frame()
names(rf_pred_subtype) <- "pred"
#to calucate AUC
comb <- cbind(y_test_subtype,rf_pred_prob_subtype,rf_pred_subtype)%>%as.data.frame()
names(comb)[1] <- "obs"
comb$obs <- as.factor(comb$obs)
comb$pred <- as.factor(comb$pred)
comb$no <- as.numeric(comb$no)
comb$yes <- as.numeric(comb$yes)
pred <- ROCR::prediction(comb[,3], comb[,1])
auc <- ROCR::performance(pred,"auc")@y.values[[1]]
print(auc)
imp_fea_subtype <- as.data.frame(rf_train_subtype$importance)
imp_fea_subtype <- imp_fea_subtype%>%rownames_to_column("subtype")%>%
   dplyr::select(subtype,MeanDecreaseAccuracy, MeanDecreaseGini)%>%dplyr::arrange(desc(MeanDecreaseGini))

#species
train_set_species <- window_db_sel_species_meta_feature %>% filter(Study!="validation")%>%dplyr::select(-Study,-species_median,-species_Richness,-species_Shannon)
x_train_species <- train_set_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train_species <- as.factor(train_set_species$nec)
set.seed(1000)
rf_train_species <- randomForest(x=x_train_species,y=y_train_species, importance = TRUE, ntree = 500, proximity = TRUE,mtry = vali_species_sel$mtry)
test_set_species <- window_db_sel_species_meta_feature %>% filter(Study=="validation")%>%
   dplyr::select(-Study,-species_median,-species_Richness,-species_Shannon)
x_test_species <- test_set_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_test_species <- test_set_species$nec
rf_pred_prob_species <- predict(rf_train_species, newdata = x_test_species, type = "prob")
rf_pred_species <- predict(rf_train_species, newdata = x_test_species)%>%as.data.frame()
names(rf_pred_species) <- "pred"
#to calucate AUC
comb <- cbind(y_test_species,rf_pred_prob_species,rf_pred_species)%>%as.data.frame()
names(comb)[1] <- "obs"
comb$obs <- as.factor(comb$obs)
comb$pred <- as.factor(comb$pred)
comb$no <- as.numeric(comb$no)
comb$yes <- as.numeric(comb$yes)
pred <- ROCR::prediction(comb[,3], comb[,1])
auc <- ROCR::performance(pred,"auc")@y.values[[1]]
print(auc)
imp_fea_species <- as.data.frame(rf_train_species$importance)
imp_fea_species <- imp_fea_species%>%rownames_to_column("species")%>%
   dplyr::select(species,MeanDecreaseAccuracy, MeanDecreaseGini)%>%dplyr::arrange(desc(MeanDecreaseGini))

#subtype+species
train_set_subtype_species <- window_db_sel_subtype_species_meta %>% filter(Study!="validation")%>%
   dplyr::select(-Study,-subtype_median,-subtype_Richness,-subtype_Shannon,-species_median,-species_Richness,-species_Shannon)
x_train_subtype_species <- train_set_subtype_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_train_subtype_species <- as.factor(train_set_subtype_species$nec)
set.seed(1000)
rf_train_subtype_species <- randomForest(x=x_train_subtype_species,y=y_train_subtype_species, importance = TRUE, ntree = 500, proximity = TRUE,mtry=vali_subtype_species_sel$mtry)
test_set_subtype_species <- window_db_sel_subtype_species_meta %>% filter(Study=="validation")%>%
   dplyr::select(-Study,-subtype_median,-subtype_Richness,-subtype_Shannon,-species_median,-species_Richness,-species_Shannon)
x_test_subtype_species <- test_set_subtype_species%>%dplyr::select(-nec)%>%column_to_rownames("SequenceID")%>%dplyr::select(sort(names(.)))
y_test_subtype_species <- test_set_subtype_species$nec
rf_pred_prob_subtype_species <- predict(rf_train_subtype_species, newdata = x_test_subtype_species, type = "prob")
rf_pred_subtype_species <- predict(rf_train_species, newdata = x_test_species)%>%as.data.frame()
names(rf_pred_subtype_species) <- "pred"
#to calucate AUC
comb <- cbind(y_test_subtype_species,rf_pred_prob_subtype_species,rf_pred_subtype_species)%>%as.data.frame()
names(comb)[1] <- "obs"
comb$obs <- as.factor(comb$obs)
comb$pred <- as.factor(comb$pred)
comb$no <- as.numeric(comb$no)
comb$yes <- as.numeric(comb$yes)
pred <- ROCR::prediction(comb[,3], comb[,1])
auc <- ROCR::performance(pred,"auc")@y.values[[1]]
print(auc)
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
# Area under the curve: 0.8182
rocobj_list$subtype$ci
# 95% CI: 0.6276-1 (DeLong)
rocobj_list$species$auc
# Area under the curve: 0.7424
rocobj_list$species$ci
# 95% CI: 0.5173-0.9675 (DeLong)
rocobj_list$subtype_species$auc
# Area under the curve: 0.8384
rocobj_list$subtype_species$ci
# 95% CI: 0.6559-1 (DeLong)

ci_list <- lapply(rocobj_list, ci.se, specificities = seq(0, 1, 0.01))
ci_list_df <- lapply(ci_list, function(ciobj)
   data.frame(prop = as.numeric(rownames(ciobj)),
              lower = ciobj[, 1],
              upper = ciobj[, 3]))
#Figure 5e
plot_curve <- ggroc(rocobj_list,size=0.5,legacy.axes =T)+
   scale_colour_manual(values = c("#0B74BA", "#42A90C", "#CC2B0A"))+
   geom_segment(aes(x=0, y=0, xend=1, yend=1),colour="grey",linetype = "dashed")+
   # geom_ribbon(data=data_ci_cb,aes(x=1-prop,ymin=`2.5%`,ymax=`97.5%`,fill=type),alpha=0.5)+
   theme_bw()+
   labs(x ="False positive rate", y="True positive rate")+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.border = element_rect(colour = "black",linewidth = 1),
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
window_db_sel_subtype_pval_1 <- read.csv("data/window_db_sel_subtype_pval_1.csv")
window_db_sel_subtype_pval_1$fdr <- p.adjust(window_db_sel_subtype_pval_1$pval_lm,method = "fdr")
window_db_sel_subtype_pval_1_0.05 <- window_db_sel_subtype_pval_1%>%filter(fdr<0.05)
which(!window_db_sel_subtype_pval_1_0.05$subtype%in%imp_fea_subtype$subtype)
imp_fea_subtype_sel <- imp_fea_subtype%>%mutate(rank=rownames(.))%>%filter(subtype%in%window_db_sel_subtype_pval_1_0.05$subtype)
imp_fea_subtype_sel$subtype <- gsub("subtype_","",imp_fea_subtype_sel$subtype,fixed = T)
imp_fea_subtype_sel$subtype <- gsub("other_peptide_antibiotics","OPA",imp_fea_subtype_sel$subtype,fixed = T)
imp_fea_subtype_sel$subtype <- gsub("macrolide-lincosamide-streptogramin","MLS",imp_fea_subtype_sel$subtype,fixed = T)
imp_fea_subtype_sel$subtype <- gsub("Klebsiella pneumoniae","KP",imp_fea_subtype_sel$subtype,fixed = T)

imp_fea_subtype_sel$subtype <- factor(imp_fea_subtype_sel$subtype,levels = imp_fea_subtype_sel$subtype)
sum(imp_fea_subtype_sel$MeanDecreaseGini)/sum(imp_fea_subtype$MeanDecreaseGini)#0.1426551
sum(imp_fea_subtype_sel$MeanDecreaseGini[1:3])/sum(imp_fea_subtype$MeanDecreaseGini)
#Figure 5f
subtype_imp_plot <- ggplot(imp_fea_subtype_sel, aes(x=subtype, y=MeanDecreaseGini)) + 
   geom_bar(stat = "identity",fill="#E09137",color="black",width = 0.6)+
   geom_text(aes(label=rank,y=MeanDecreaseGini+0.1),size=3)+
   theme_bw() + # remove the backgroud
   labs(x = element_blank(), y="MeanDecreaseGini")+
   scale_y_continuous(expand = c(0,0),limits = c(0,2.2),breaks = c(0,0.5,1.0,1.5,2.0),labels = c("0","0.5","1.0","1.5","2.0"))+
   theme(panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth = 0.5)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   # theme(aspect.ratio=1/3.8)+
   theme(legend.position="none",
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.title.y =  element_text(size=12, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour = "black",angle = 90,hjust = 1,vjust=0.5),
         plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
