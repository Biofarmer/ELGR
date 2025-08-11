library(RColorBrewer)
library(dplyr)
library(tibble)
library(Rmisc)
library(vegan)
library(ggplot2)
library(coin)
library(tidyr)
library(stringr)


##############################################################################################################################
#ARG dynamics
##############################################################################################################################
preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")
preterm_metadata <- preterm_metadata%>%
   mutate(month=case_when(DOL>=0 & DOL<=7 ~ "0.25",
                          DOL>7 & DOL<=14 ~ "0.5",
                          DOL>14 & DOL<=21 ~ "0.75",
                          DOL>21 & DOL<=30 ~ "1",
                          DOL>30 & DOL<=60 ~ "2",
                          DOL>60 & DOL<=90 ~ "3",
                          DOL>90 & DOL<=180 ~ "6",
                          DOL>180 & DOL<=360 ~ "12",
                          DOL>360 & DOL<=540 ~ "18",
                          DOL>540 ~ "36"))%>%
   mutate(st_subjectID = paste0(Study,"_",SubjectID))
preterm_metadata$month <- as.numeric(preterm_metadata$month)
#subtypes
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_type <- preterm_resistome%>%filter(group=="type")%>%dplyr::select(-group)
type_sum_abund <- preterm_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("SequenceID"=sample,"type_total_abund"=rowSums(.[-1]))
preterm_subtype <- preterm_resistome%>%filter(group=="subtype")%>%dplyr::select(-group)
preterm_subtype_t <- preterm_subtype%>%column_to_rownames("names")%>%t()%>%as.data.frame()
subtype_abund <- preterm_subtype_t%>%rownames_to_column("sample")%>%dplyr::summarise("SequenceID"=sample,"subtype_total_abund"=rowSums(.[-1]))
subtype_richness <- rowSums(preterm_subtype_t != 0)
subtype_shannon <- diversity(preterm_subtype_t, index = "shannon", MARGIN=1) 
subtype_richness_df <- subtype_richness%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(subtype_richness_df)[2] <- "Richness"
subtype_shannon_df <- subtype_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(subtype_shannon_df)[2] <- "Shannon"
subtype_alpha <- left_join(subtype_richness_df,subtype_shannon_df,by="SequenceID")
subtype_alpha <- left_join(subtype_alpha,subtype_abund,by="SequenceID")
subtype_alpha <- left_join(subtype_alpha,type_sum_abund,by="SequenceID")
subtype_alpha_meta <- left_join(preterm_metadata,subtype_alpha,by="SequenceID")
subtype_alpha_meta <- subtype_alpha_meta%>%filter(!is.na(month))
pairwise.wilcox.test(subtype_alpha_meta$Shannon, subtype_alpha_meta$month,p.adjust.method="fdr")
pairwise.wilcox.test(subtype_alpha_meta$subtype_total_abund, subtype_alpha_meta$month,p.adjust.method="fdr")
pairwise.wilcox.test(subtype_alpha_meta$type_total_abund, subtype_alpha_meta$month,p.adjust.method="fdr")

subtype_alpha_meta$month <- as.factor(subtype_alpha_meta$month)
#Figure 4b
subtype_shannon_time_p <- ggplot(subtype_alpha_meta, aes(x=month, y=Shannon)) +
   geom_violin(aes(fill=month),alpha=0.8,trim=T)+
   geom_boxplot(width=0.1,outlier.size = 1)+
   scale_fill_manual(values = c("#556997","#8286A3","#ACA4AE","#D5C3B9","#FFE4C4","#FBC39D","#F4A277","#E98152","#DC5F2D","#CD3700"))+
   labs(x ="Timepoints (month)", y="Shannon index (based on ARG subtypes)")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"))

subtype_alpha_meta <- subtype_alpha_meta%>%group_by(month)%>%dplyr::mutate(month_2=cur_group_id())
#Figure 4a
subtype_abund_time_p <- ggplot(subtype_alpha_meta, aes(x=month_2, y=subtype_total_abund)) +
   geom_smooth(method = "loess",color="#4994C4")+
   labs(x ="Timepoints (month)", y="ARG abundance (ARGs-OAP,CAPC)")+
   theme_bw()+
   scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10),
                      labels = c("1"="0.25","2"="0.5","3"="0.75","4"="1","5"="2",
                                 "6"="3","7"="6","8"="12","9"="18","10"="36"))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

#beta diversity
preterm_subtype_t_sel <- preterm_subtype_t%>%filter(rownames(.)%in%subtype_alpha_meta$SequenceID)
preterm_subtype_beta_dist <- vegdist(preterm_subtype_t_sel, method = "bray")
preterm_metadata_subtype <- preterm_metadata%>%filter(!is.na(DOL))%>%arrange(match(SequenceID,rownames(preterm_subtype_t_sel)))
preterm_metadata_subtype$month <- as.factor(preterm_metadata_subtype$month)
set.seed(1000)
vegan::adonis2(preterm_subtype_beta_dist ~ month,data=preterm_metadata_subtype,permutations = 1000)
subtype_beta_pco <- cmdscale(preterm_subtype_beta_dist, k=2, eig = T)
subtype_beta_pco_axis.1.title <- paste('PCoA1 [', 
                      round((subtype_beta_pco$eig[1]/sum(subtype_beta_pco$eig))*100,1),
                      '%]', sep='')
subtype_beta_pco_axis.2.title <- paste('PCoA2 [', 
                      round((subtype_beta_pco$eig[2]/sum(subtype_beta_pco$eig))*100,1),
                      '%]', sep='')

subtype_beta_pco.point <- as.data.frame(subtype_beta_pco$points)%>%rownames_to_column('SequenceID')
subtype_beta_pco.tabble <- right_join(preterm_metadata, subtype_beta_pco.point, by="SequenceID")

subtype_beta_df.plot <- tibble(Axis1 = subtype_beta_pco.tabble$V1,
                               Axis2 = subtype_beta_pco.tabble$V2,
                               Sample_ID = subtype_beta_pco.tabble$SubjectID,
                               Time=subtype_beta_pco.tabble$month)
subtype_beta_df.plot <- subtype_beta_df.plot%>%group_by(Time)%>%dplyr::mutate(Time2=cur_group_id())

#get the color used by beta diversity and then for alpha diversity as above.
dd <- as.data.frame(ggplot_build(subtype_beta_time_p)$data[1])
dd$cb <- paste0(dd$x,"_",dd$y) 
subtype_beta_df.plot_2 <- subtype_beta_df.plot
subtype_beta_df.plot_2$cb <- paste0(subtype_beta_df.plot_2$Axis1,"_",subtype_beta_df.plot_2$Axis2)
subtype_beta_df.plot_2 <- left_join(subtype_beta_df.plot_2,dd%>%dplyr::select(cb,colour),by="cb")
subtype_beta_df.plot_2_sum <- subtype_beta_df.plot_2%>%distinct(colour)%>%arrange(Time)
#Figure 4c
subtype_beta_time_p <- ggplot(data=subtype_beta_df.plot,aes(x=Axis1, y=Axis2, col=Time2)) +
   geom_point(size=1.5, alpha=1) + 
   scale_color_gradient2(midpoint=5, low="dodgerblue4", mid="bisque",high="orangered3", space ="Lab",
                         breaks=c(1,2,3,4,5,6,7,8,9,10),
                         labels=c("0.25","0.5","0.75","1","2","3","6","12","18","36"))+
   theme_bw()+ 
   xlab(subtype_beta_pco_axis.1.title) + ylab(subtype_beta_pco_axis.2.title) +
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "top",
         legend.text = element_text(size=10, colour = "black",angle = 90,hjust = 1,vjust = 1),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,1,1,1),"mm"))

#######################################################################################################################################
#amr
#######################################################################################################################################
amr_profile_coverm <- read.csv("data/amr_preterm_resistome_coverm_bins_gtdb_stat_plasmid.csv")
amr_profile_coverm_gene_sum <- amr_profile_coverm%>%
   group_by(sample,Element_symbol)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%spread(.,sample,rpkm_total)
amr_profile_coverm_gene_sum[is.na(amr_profile_coverm_gene_sum)] <- 0
amr_profile_coverm_gene_sum_t <- amr_profile_coverm_gene_sum%>%column_to_rownames("Element_symbol")%>%t()%>%as.data.frame()
#beta diversity
amr_profile_coverm_gene_sum_t_2 <- amr_profile_coverm_gene_sum_t%>%filter(rownames(.)%in%amr_alpha_meta$SequenceID)
amr_beta_dist <- vegdist(amr_profile_coverm_gene_sum_t_2, method = "bray")
preterm_metadata_amr <- preterm_metadata%>%
   filter(SequenceID%in%rownames(amr_profile_coverm_gene_sum_t_2))%>%
   arrange(match(SequenceID,rownames(amr_profile_coverm_gene_sum_t_2)))
preterm_metadata_amr$month <- as.factor(preterm_metadata_amr$month)
set.seed(1000)
vegan::adonis2(amr_beta_dist ~ month,data=preterm_metadata_amr,permutations = 1000)
amr_beta_pco <- cmdscale(amr_beta_dist, k=2, eig = T)
amr_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                   round((amr_beta_pco$eig[1]/sum(amr_beta_pco$eig))*100,1),
                                   '%]', sep='')
amr_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                   round((amr_beta_pco$eig[2]/sum(amr_beta_pco$eig))*100,1),
                                   '%]', sep='')

amr_beta_pco.point <- as.data.frame(amr_beta_pco$points)%>%rownames_to_column('SequenceID')
amr_beta_pco.tabble <- merge(preterm_metadata_amr, amr_beta_pco.point, by="SequenceID", all = TRUE)

amr_beta_df.plot <- tibble(Axis1 = amr_beta_pco.tabble$V1,
                           Axis2 = amr_beta_pco.tabble$V2,
                           Sample_ID = amr_beta_pco.tabble$SubjectID,
                           Time=amr_beta_pco.tabble$month)
amr_beta_df.plot <- amr_beta_df.plot%>%group_by(Time)%>%dplyr::mutate(Time2=cur_group_id())
#Supplementary Figure 9b 
amr_beta_time_p <- ggplot(data=amr_beta_df.plot,aes(x=Axis1, y=Axis2, col=Time2)) +
   geom_point(size=1.5, alpha=1) + 
   scale_color_gradient2(midpoint=5, low="dodgerblue4", mid="bisque",high="orangered3", space ="Lab",
                         breaks=c(1,2,3,4,5,6,7,8,9,10),
                         labels=c("0.25","0.5","0.75","1","2","3","6","12","18","36"))+
   theme_bw()+ 
   xlab(amr_beta_pco_axis.1.title) + ylab(amr_beta_pco_axis.2.title) +
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "top",
         legend.text = element_text(size=10, colour = "black",angle = 90,hjust = 1,vjust = 1),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"))
##############################################################################################################################
#covariate analysis
##############################################################################################################################
preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")
preterm_metadata <- preterm_metadata%>%
   mutate(month=case_when(DOL>=0 & DOL<=7 ~ "0.25",
                          DOL>7 & DOL<=14 ~ "0.5",
                          DOL>14 & DOL<=21 ~ "0.75",
                          DOL>21 & DOL<=30 ~ "1",
                          DOL>30 & DOL<=60 ~ "2",
                          DOL>60 & DOL<=90 ~ "3",
                          DOL>90 & DOL<=180 ~ "6",
                          DOL>180 & DOL<=360 ~ "12",
                          DOL>360 & DOL<=540 ~ "18",
                          DOL>540 ~ "36"))%>%
   mutate(Feed_2=case_when(Feed=="breast"~"Exclusive",
                           Feed=="combined"~"Mix",
                           Feed=="nofeed"~"fasting",
                           Feed=="formula" | Feed =="solid"~"Non-breastmilk"))%>%
   mutate(Birth_weight_2=case_when(Birth_weight < 1000 ~"ELBW",
                                   Birth_weight >= 1000 & Birth_weight < 1500 ~"VLBW",
                                   Birth_weight >= 1500 & Birth_weight < 2500 ~"LBW",
                                   Birth_weight >= 2500 ~"NBW"))%>%
   mutate(Gestational_age_2=case_when(Gestational_age < 28~"extremely preterm",
                                      Gestational_age >= 28 & Gestational_age < 32 ~"very preterm",
                                      Gestational_age >= 32 & Gestational_age < 34 ~"moderate preterm",
                                      Gestational_age >= 34 ~"late preterm"))

preterm_metadata$month <- as.numeric(preterm_metadata$month)

subtype_adonis_strata <- read.csv("data/figure/covariate/subtype_covar_contri_adonis_strata.csv")
amr_adonis_strata <- read.csv("data/figure/covariate/amr_covar_contri_adonis_strata.csv")

subtype_adonis_strata_sel <- subtype_adonis_strata%>%
   filter(covariates%in%c("cb","month","Delivery","Country","Gender","Feed_2","Infant_AB","Birth_weight_2","Gestational_age_2"))%>%
   mutate(fdr=p.adjust(.$p_val,method = "fdr"))%>%
   arrange(desc(r2))%>%dplyr::select(covariates,r2,fdr)
subtype_adonis_strata_sel$r2 <- subtype_adonis_strata_sel$r2*100
subtype_adonis_strata_sel$covariates[subtype_adonis_strata_sel$covariates=="cb"] <- "Subjects"
subtype_adonis_strata_sel$covariates[subtype_adonis_strata_sel$covariates=="month"] <- "Infant age (month)"
subtype_adonis_strata_sel$covariates[subtype_adonis_strata_sel$covariates=="Gestational_age_2"] <- "Gestational age"
subtype_adonis_strata_sel$covariates[subtype_adonis_strata_sel$covariates=="Infant_AB"] <- "Antibiotics"
subtype_adonis_strata_sel$covariates[subtype_adonis_strata_sel$covariates=="Birth_weight_2"] <- "Birth weight"
subtype_adonis_strata_sel$covariates[subtype_adonis_strata_sel$covariates=="Feed_2"] <- "Feeding pattern"
subtype_adonis_strata_sel$covariates[subtype_adonis_strata_sel$covariates=="Gender"] <- "Sex"
subtype_adonis_strata_sel_2 <- subtype_adonis_strata_sel%>%dplyr::select(-fdr)%>%column_to_rownames("covariates")%>%t()%>%as.data.frame()
subtype_adonis_strata_sel_2 <- rbind(rep(60,9),rep(0,9),subtype_adonis_strata_sel_2)

#Supplementary Figure 9c
radarchart(subtype_adonis_strata_sel_2, axistype=0, seg=3,
           #custom polygon
           pcol="#C62916", pfcol=NA, plwd=2, 
           #custom the grid
           cglcol="#959494", cglty=2, axislabcol="grey", caxislabels=seq(0,60,20), cglwd=1,calcex = c(0),
           #custom labels
           vlcex=0.8)

amr_adonis_strata_sel <- amr_adonis_strata%>%
   filter(covariates%in%c("cb","month","Delivery","Country","Gender","Feed_2","Infant_AB","Birth_weight_2","Gestational_age_2"))%>%
   mutate(fdr=p.adjust(.$p_val,method = "fdr"))%>%
   arrange(desc(r2))%>%dplyr::select(covariates,r2,fdr)
amr_adonis_strata_sel$r2 <- amr_adonis_strata_sel$r2*100
amr_adonis_strata_sel$covariates[amr_adonis_strata_sel$covariates=="cb"] <- "Subjects"
amr_adonis_strata_sel$covariates[amr_adonis_strata_sel$covariates=="month"] <- "Infant age (month)"
amr_adonis_strata_sel$covariates[amr_adonis_strata_sel$covariates=="Gestational_age_2"] <- "Gestational age"
amr_adonis_strata_sel$covariates[amr_adonis_strata_sel$covariates=="Infant_AB"] <- "Antibiotics"
amr_adonis_strata_sel$covariates[amr_adonis_strata_sel$covariates=="Birth_weight_2"] <- "Birth weight"
amr_adonis_strata_sel$covariates[amr_adonis_strata_sel$covariates=="Feed_2"] <- "Feeding pattern"
amr_adonis_strata_sel$covariates[amr_adonis_strata_sel$covariates=="Gender"] <- "Sex"
amr_adonis_strata_sel_2 <- amr_adonis_strata_sel%>%dplyr::select(-fdr)%>%column_to_rownames("covariates")%>%t()%>%as.data.frame()
amr_adonis_strata_sel_2 <- rbind(rep(60,9),rep(0,9),amr_adonis_strata_sel_2)

#Supplementary Figure 10a
radarchart(amr_adonis_strata_sel_2, axistype=0, seg=3,
           #custom polygon
           pcol="#C62916", pfcol=NA, plwd=2, 
           #custom the grid
           cglcol="#959494", cglty=2, axislabcol="grey", caxislabels=seq(0,60,20), cglwd=1,calcex = c(0),
           #custom labels
           vlcex=0.8)

cor.test(preterm_metadata$Gestational_age,preterm_metadata$Birth_weight)#r = 0.87401, p-value < 2.2e-16
#Supplementary Figure 10b
corr_gestation_birthweight <- ggplot(preterm_metadata, aes(x=Gestational_age, y=Birth_weight)) +
   geom_point(shape=21,color=alpha("black",0.5),fill=alpha("black",0.5),size=2)+
   geom_smooth(method=lm, formula =  y ~ poly(x,1), level=0,colour="#DA732D",linewidth=2) +
   labs(x ="Gestational age (week)", y="Birth weight (gram)")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"))


#subtypes
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_subtype <- preterm_resistome%>%filter(group=="subtype")%>%dplyr::select(-group)
preterm_subtype_t <- preterm_subtype%>%column_to_rownames("names")%>%t()%>%as.data.frame()
subtype_dist_beta <- vegdist(preterm_subtype_t, method = "bray")
subtype_dist_beta_2 <- subtype_dist_beta%>%as.matrix()%>%as.data.frame()%>%rownames_to_column("run_1")
#adonis2
subtype_covar_contri <- list()
month <- c("0.25","0.5","0.75","1","2","3","6","12","18","36")
variables <- c("Delivery","Country","Gender","Feed_2","Infant_AB","Birth_weight_2","Gestational_age_2")
for (p in month) {
   # p="1"
   print(p)
   subtype_covar_contri[[p]] <- data.frame("covariates"="x","level"=0,"level_min"=0,"Df"=0,"n"=0, "r2"=0,"P_val"=0,"permdisp_p_val"=0,stringsAsFactors=F)[-1,]
   for (v in variables) {
      # v="Birth_weight_2"
      print(v)
      subtype_db_1 <- preterm_metadata%>%filter(month==p)%>%dplyr::select(SequenceID,Study,v)%>%setNames(c("run","Study","covariate"))%>%filter(!is.na(covariate))
      subtype_dist_beta_3 <- subtype_dist_beta_2%>%filter(run_1%in%subtype_db_1$run)
      subtype_dist_beta_4 <- subtype_dist_beta_3%>%dplyr::select(run_1,subtype_dist_beta_3$run_1)%>%column_to_rownames("run_1")
      subtype_dist_beta_5 <- subtype_dist_beta_4%>%as.matrix()%>%as.dist()
      subtype_db_2 <- subtype_db_1%>%filter(run%in%subtype_dist_beta_3$run_1)%>%arrange(match(run,subtype_dist_beta_3$run_1))
      print(identical(rownames(subtype_dist_beta_4),subtype_db_2$run))
      subtype_num_per_group <- subtype_db_2%>%group_by(covariate)%>%dplyr::summarise(n_per_cat=n())
      if (min(subtype_num_per_group$n_per_cat)<1 | length(unique(subtype_db_2$covariate))<2) {
         next
      }else{
         print(min(subtype_num_per_group$n_per_cat))
         if (!v%in%c("Birth_weight","Gestational_age")) {
            subtype_db_2$covariate <- as.factor(subtype_db_2$covariate)
            set.seed(1000)
            subtype_a <- try(with(subtype_db_2, adonis2(subtype_dist_beta_5 ~ covariate, data = subtype_db_2, permutations = 1000, strata = Study)))
            if (unique(class(subtype_a)=="try-error")) {
               print("error")
               next
            }else{
               subtype_d <- with(subtype_db_2,betadisper(subtype_dist_beta_5, covariate))
               set.seed(1000)
               subtype_d_1 <- permutest(subtype_d,permutations = 1000)
               subtype_covar_contri[[p]][nrow(subtype_covar_contri[[p]])+1,] <- c(v,length(unique(subtype_db_2$covariate)),min(table(subtype_db_2$covariate)),subtype_a$Df[1],nrow(subtype_db_2),
                                                                                  subtype_a$R2[1],subtype_a$`Pr(>F)`[1],subtype_d_1$tab$`Pr(>F)`[1])
            }
         }else{
            set.seed(1000)
            subtype_a <- try(with(subtype_db_2, adonis2(subtype_dist_beta_5 ~ covariate, data = subtype_db_2, permutations = 1000, strata = Study)))
            if (unique(class(subtype_a)=="try-error")) {
               print("error")
               next
            }else{
               subtype_d <- with(subtype_db_2,betadisper(subtype_dist_beta_5, covariate))
               set.seed(1000)
               subtype_d_1 <- permutest(subtype_d,permutations = 1000)
               subtype_covar_contri[[p]][nrow(subtype_covar_contri[[p]])+1,] <- c(v,NA,NA,subtype_a$Df[1],nrow(subtype_db_2),
                                                                                  subtype_a$R2[1],subtype_a$`Pr(>F)`[1],subtype_d_1$tab$`Pr(>F)`[1])
            }
         }
      }
   }
}

for (i in 1:length(subtype_covar_contri)) {
   subtype_covar_contri[[i]]$P_val <- as.numeric(subtype_covar_contri[[i]]$P_val)
   subtype_covar_contri[[i]]$permdisp_p_val <- as.numeric(subtype_covar_contri[[i]]$permdisp_p_val)
   subtype_covar_contri[[i]]$p_adj <- p.adjust(subtype_covar_contri[[i]]$P_val, method = "fdr") 
   subtype_covar_contri[[i]]$p_adj_disp <- p.adjust(subtype_covar_contri[[i]]$permdisp_p_val, method = "fdr") 
   name <- names(subtype_covar_contri[i])
   subtype_covar_contri[[i]]$month <- name
}

subtype_covar_contri_cb <- rbind(subtype_covar_contri[[1]],subtype_covar_contri[[2]],subtype_covar_contri[[3]],subtype_covar_contri[[4]],
                                 subtype_covar_contri[[5]],subtype_covar_contri[[6]],subtype_covar_contri[[7]],
                                 subtype_covar_contri[[8]],subtype_covar_contri[[9]],subtype_covar_contri[[10]])
subtype_covar_contri_cb$r2 <- as.numeric(subtype_covar_contri_cb$r2)
subtype_covar_contri_cb$level_min <- as.numeric(subtype_covar_contri_cb$level_min)
subtype_covar_contri_cb$month <- factor(subtype_covar_contri_cb$month,levels = c("0.25","0.5","0.75","1","2","3","6","12","18","36"))
subtype_covar_contri_cb$padj_per <- p.adjust(subtype_covar_contri_cb$permdisp_p_val,method = "fdr")
subtype_covar_contri_cb_sel <- subtype_covar_contri_cb%>%filter(!covariates%in%c("Birth_weight","Gestational_age","Instrument_2"))
subtype_covar_contri_cb_sel <- subtype_covar_contri_cb_sel%>%group_by(month)%>%dplyr::mutate(month_2=cur_group_id())
subtype_covar_contri_cb_sel_sig <- subtype_covar_contri_cb_sel%>%filter(p_adj<0.05)
subtype_covar_contri_cb_sel$covariates <- factor(subtype_covar_contri_cb_sel$covariates,
                                                 levels = rev(c("Birth_weight_2","Gestational_age_2","Infant_AB","Feed_2","Gender",
                                                                "Delivery","Country")))
#Figure 4e
subtype_covariates <- ggplot()+
   geom_tile(data=subtype_covar_contri_cb_sel,aes(x=month_2, y=covariates),fill="white",color="white",linewidth = 0.5)+
   geom_point(data=subtype_covar_contri_cb_sel,aes(x=month_2, y=covariates,size=r2, fill=covariates),shape=21)+
   geom_tile(data=subtype_covar_contri_cb_sel_sig,aes(x=month_2, y=covariates),fill=NA,color="#C62916",linewidth = 0.5)+
   theme_bw()+
   theme(panel.grid = element_blank())+
   scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10),
                      labels = c("1"="0.25","2"="0.5","3"="0.75","4"="1","5"="2","6"="3","7"="6","8"="12","9"="18","10"="36"))+
   scale_y_discrete(labels=c("Birth_weight_2"="Birth weight","Gestational_age_2"="Gestational age","Infant_AB"="Antibiotics",
                             "Feed_2"="Feeding pattern","Gender"="Sex","Delivery"="Delivery mode","Country"="Country"))+
   scale_fill_manual(values = c("#9E4546","#984EA3","#657AD7","#30CAB3","#DEB141","#4DAF4A","#377EB8"))+
   labs(x ="Timepoints (month)")+
   theme(panel.border = element_rect(colour = "black",linewidth=1),axis.line = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=14, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=14, colour = "black"))

#amr
amr_preterm_profile_coverm <- read.csv("data/amr_preterm_resistome_coverm_bins_gtdb_stat_plasmid.csv")
amr_preterm_profile_coverm_gene_sum <- amr_preterm_profile_coverm%>%
   group_by(sample,Element_symbol)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%spread(.,sample,rpkm_total)
amr_preterm_profile_coverm_gene_sum[is.na(amr_preterm_profile_coverm_gene_sum)] <- 0
amr_preterm_profile_coverm_gene_sum_t <- amr_preterm_profile_coverm_gene_sum%>%
   dplyr::select(Element_symbol,which(colSums(.[-1])!=0)+1)%>%column_to_rownames("Element_symbol")%>%t()%>%as.data.frame()
amr_preterm_profile_coverm_gene_sum_t <- amr_preterm_profile_coverm_gene_sum_t[which(colSums(amr_preterm_profile_coverm_gene_sum_t)!=0)]

amr_dist_beta <- vegdist(amr_preterm_profile_coverm_gene_sum_t, method = "bray")
amr_dist_beta_2 <- amr_dist_beta%>%as.matrix()%>%as.data.frame()%>%rownames_to_column("run_1")

#adonis2
amr_covar_contri <- list()
month <- c("0.25","0.5","0.75","1","2","3","6","12","18","36")
variables <- c("Delivery","Country","Gender","Feed_2","Infant_AB","Birth_weight_2","Gestational_age_2")
for (p in month) {
   # p="1"
   print(p)
   amr_covar_contri[[p]] <- data.frame("covariates"="x","level"=0,"level_min"=0,"Df"=0,"n"=0, "r2"=0,"P_val"=0,"permdisp_p_val"=0,stringsAsFactors=F)[-1,]
   for (v in variables) {
      # v="Birth_weight_2"
      print(v)
      amr_db_1 <- preterm_metadata%>%filter(month==p)%>%dplyr::select(SequenceID,Study,v)%>%setNames(c("run","Study","covariate"))%>%filter(!is.na(covariate))
      amr_dist_beta_3 <- amr_dist_beta_2%>%filter(run_1%in%amr_db_1$run)
      amr_dist_beta_4 <- amr_dist_beta_3%>%dplyr::select(run_1,amr_dist_beta_3$run_1)%>%column_to_rownames("run_1")
      amr_dist_beta_5 <- amr_dist_beta_4%>%as.matrix()%>%as.dist()
      amr_db_2 <- amr_db_1%>%filter(run%in%amr_dist_beta_3$run_1)%>%arrange(match(run,amr_dist_beta_3$run_1))
      print(identical(rownames(amr_dist_beta_4),amr_db_2$run))
      amr_num_per_group <- amr_db_2%>%group_by(covariate)%>%dplyr::summarise(n_per_cat=n())
      if (min(amr_num_per_group$n_per_cat)<1 | length(unique(amr_db_2$covariate))<2) {
         next
      }else{
         print(min(amr_num_per_group$n_per_cat))
         if (!v%in%c("Birth_weight","Gestational_age")) {
            amr_db_2$covariate <- as.factor(amr_db_2$covariate)
            set.seed(1000)
            amr_a <- try(with(amr_db_2, adonis2(amr_dist_beta_5 ~ covariate, data = amr_db_2, permutations = 1000, strata = Study)))
            if (unique(class(amr_a)=="try-error")) {
               print("error")
               next
            }else{
               amr_d <- with(amr_db_2,betadisper(amr_dist_beta_5, covariate))
               set.seed(1000)
               amr_d_1 <- permutest(amr_d,permutations = 1000)
               amr_covar_contri[[p]][nrow(amr_covar_contri[[p]])+1,] <- c(v,length(unique(amr_db_2$covariate)),min(table(amr_db_2$covariate)),amr_a$Df[1],nrow(amr_db_2),
                                                                          amr_a$R2[1],amr_a$`Pr(>F)`[1],amr_d_1$tab$`Pr(>F)`[1])
            }
         }else{
            set.seed(1000)
            amr_a <- try(with(amr_db_2, adonis2(amr_dist_beta_5 ~ covariate, data = amr_db_2, permutations = 1000, strata = Study)))
            if (unique(class(amr_a)=="try-error")) {
               print("error")
               next
            }else{
               amr_d <- with(amr_db_2,betadisper(amr_dist_beta_5, covariate))
               set.seed(1000)
               amr_d_1 <- permutest(amr_d,permutations = 1000)
               amr_covar_contri[[p]][nrow(amr_covar_contri[[p]])+1,] <- c(v,NA,NA,amr_a$Df[1],nrow(amr_db_2),
                                                                          amr_a$R2[1],amr_a$`Pr(>F)`[1],amr_d_1$tab$`Pr(>F)`[1])
            }
         }
      }
   }
}

for (i in 1:length(amr_covar_contri)) {
   amr_covar_contri[[i]]$P_val <- as.numeric(amr_covar_contri[[i]]$P_val)
   amr_covar_contri[[i]]$permdisp_p_val <- as.numeric(amr_covar_contri[[i]]$permdisp_p_val)
   amr_covar_contri[[i]]$p_adj <- p.adjust(amr_covar_contri[[i]]$P_val, method = "fdr") 
   amr_covar_contri[[i]]$p_adj_disp <- p.adjust(amr_covar_contri[[i]]$permdisp_p_val, method = "fdr") 
   name <- names(amr_covar_contri[i])
   amr_covar_contri[[i]]$month <- name
}

amr_covar_contri_cb <- rbind(amr_covar_contri[[1]],amr_covar_contri[[2]],amr_covar_contri[[3]],amr_covar_contri[[4]],
                             amr_covar_contri[[5]],amr_covar_contri[[6]],amr_covar_contri[[7]],amr_covar_contri[[8]],
                             amr_covar_contri[[9]],amr_covar_contri[[10]])
amr_covar_contri_cb$r2 <- as.numeric(amr_covar_contri_cb$r2)
amr_covar_contri_cb$level_min <- as.numeric(amr_covar_contri_cb$level_min)
amr_covar_contri_cb$month <- factor(amr_covar_contri_cb$month,levels = c("0.25","0.5","0.75","1","2","3","6","12","18","36"))
amr_covar_contri_cb$padj_per <- p.adjust(amr_covar_contri_cb$permdisp_p_val,method = "fdr")
amr_covar_contri_cb_sel <- amr_covar_contri_cb%>%filter(!covariates%in%c("Birth_weight","Gestational_age","Instrument_2"))
amr_covar_contri_cb_sel <- amr_covar_contri_cb_sel%>%group_by(month)%>%dplyr::mutate(month_2=cur_group_id())
amr_covar_contri_cb_sel_sig <- amr_covar_contri_cb_sel%>%filter(p_adj<0.05)
amr_covar_contri_cb_sel$covariates <- factor(amr_covar_contri_cb_sel$covariates,
                                             levels = rev(c("Birth_weight_2","Gestational_age_2","Infant_AB","Feed_2","Gender",
                                                            "Delivery","Country")))
#Supplementary Figure 9d
amr_covariates <- ggplot()+
   geom_tile(data=amr_covar_contri_cb_sel,aes(x=month_2, y=covariates),fill="white",color="white",linewidth = 0.5)+
   geom_point(data=amr_covar_contri_cb_sel,aes(x=month_2, y=covariates,size=r2, fill=covariates),shape=21)+
   geom_tile(data=amr_covar_contri_cb_sel_sig,aes(x=month_2, y=covariates),fill=NA,color="#C62916",linewidth = 0.5)+
   theme_bw()+
   theme(panel.grid = element_blank())+
   scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10),
                      labels = c("1"="0.25","2"="0.5","3"="0.75","4"="1","5"="2","6"="3","7"="6","8"="12","9"="18","10"="36"))+
   scale_y_discrete(labels=c("Birth_weight_2"="Birth weight","Gestational_age_2"="Gestational age","Infant_AB"="Antibiotics",
                             "Feed_2"="Feeding pattern","Gender"="Sex","Delivery"="Delivery mode","Country"="Country"))+
   scale_fill_manual(values = c("#9E4546","#984EA3","#657AD7","#30CAB3","#DEB141","#4DAF4A","#377EB8"))+
   labs(x ="Timepoints (month)")+
   theme(panel.border = element_rect(colour = "black",linewidth=1),axis.line = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=14, colour="black"),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=14, colour = "black"))
########################################################################################################################################
#Maaslin2
########################################################################################################################################
preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")
preterm_metadata <- preterm_metadata%>%
   mutate(month=case_when(DOL>=0 & DOL<=7 ~ "0.25",
                          DOL>7 & DOL<=14 ~ "0.5",
                          DOL>14 & DOL<=21 ~ "0.75",
                          DOL>21 & DOL<=30 ~ "1",
                          DOL>30 & DOL<=60 ~ "2",
                          DOL>60 & DOL<=90 ~ "3",
                          DOL>90 & DOL<=180 ~ "6",
                          DOL>180 & DOL<=360 ~ "12",
                          DOL>360 & DOL<=540 ~ "18",
                          DOL>540 ~ "36"))%>%
   mutate(Feed_2=case_when(Feed=="breast"~"Exclusive",
                           Feed=="combined"~"Mix",
                           Feed=="nofeed"~"nofeed",
                           Feed=="formula" | Feed =="solid"~"Non_breastmilk"))%>%
   mutate(Birth_weight_2=case_when(Birth_weight < 1000 ~"ELBW",
                                   Birth_weight >= 1000 & Birth_weight < 1500 ~"VLBW",
                                   Birth_weight >= 1500 & Birth_weight < 2500 ~"LBW",
                                   Birth_weight >= 2500 ~"NBW"))%>%
   mutate(Gestational_age_2=case_when(Gestational_age < 28~"extremely_preterm",
                                      Gestational_age >= 28 & Gestational_age < 32 ~"very_preterm",
                                      Gestational_age >= 32 & Gestational_age < 34 ~"moderate_preterm",
                                      Gestational_age >= 34 ~"late_preterm"))%>%
   mutate(st_subjectID = paste0(Study,"_",SubjectID))
preterm_metadata$month <- as.numeric(preterm_metadata$month)
preterm_metadata$st_subjectID <- gsub("_Repeat","",preterm_metadata$st_subjectID,fixed = T)
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
meta_input <- preterm_metadata%>%
   dplyr::select(SequenceID,st_subjectID,Delivery,Country,Gestational_age_2,Feed_2,month,Gender,Birth_weight_2,Infant_AB)%>%
   column_to_rownames("SequenceID")
meta_input$month <- as.factor(meta_input$month)
meta_input$Delivery <- gsub("vaginal","a_vaginal",meta_input$Delivery,fixed = T)
meta_input_sum_1 <- meta_input%>%
   dplyr::summarise_at(vars(-st_subjectID),n_distinct,na.rm = T)%>%t()%>%as.data.frame()%>%rownames_to_column("covariate")
meta_input_sum_2 <- meta_input%>%
   dplyr::summarise_at(vars(-st_subjectID),~sum(is.na(.)))%>%
   t()%>%as.data.frame()%>%rownames_to_column("covariate")%>%
   mutate(total=nrow(meta_input))%>%
   mutate(perc=V1/total)
names(meta_input_sum_1)[2] <- "levels"
names(meta_input_sum_2)[2] <- "NA_number"   
meta_input_sum <- left_join(meta_input_sum_2,meta_input_sum_1,by="covariate")
meta_input_sum <- meta_input_sum%>%
   mutate(reference=case_when(covariate=="Country" ~ "Country,USA",
                              covariate=="Feed_2" ~ "Feed_2,Non_breastmilk",
                              covariate=="Birth_weight_2" ~ "Birth_weight_2,NBW",
                              covariate=="Gestational_age_2" ~ "Gestational_age_2,late_preterm",
                              covariate=="Delivery" ~ "Delivery,a_vaginal",
                              covariate=="Gender" ~ "Gender,female",
                              covariate=="Infant_AB" ~ "Infant_AB,no",
                              covariate=="month" ~ "month,0"))
meta_input_sum_sel <- meta_input_sum%>%filter(perc<0.5)%>%filter(levels>=2)
#type
preterm_type <- preterm_resistome%>%filter(group=="type")%>%dplyr::select(-group)
preterm_type_t <- preterm_type%>%column_to_rownames("names")%>%t()%>%as.data.frame()
meta_input_type <- meta_input%>%rownames_to_column("run")%>%arrange(match(run,rownames(preterm_type_t)))%>%column_to_rownames("run")
identical(rownames(meta_input_type),rownames(preterm_type_t))
table(meta_input_type$Feed_2)
fit_data_type = Maaslin2(
   input_data = preterm_type_t, 
   input_metadata = meta_input_type,
   output = 'maaslin2/type_prev5/all/',
   normalization = 'TSS',
   transform = 'AST',
   min_abundance=0,
   min_prevalence=0.05,
   plot_scatter = F,
   fixed_effects = c('Delivery','Country','Gestational_age_2','Feed_2','month','Gender','Birth_weight_2','Infant_AB'),
   random_effects = c('st_subjectID'),
   reference = c('Country,USA','Feed_2,Non_breastmilk','month,0','Birth_weight_2,NBW','Gestational_age_2,late_preterm',
                 'Gender,female','Infant_AB,no','Delivery,a_vaginal'))
month <- c("0.25","0.5","0.75","1","2","3","6","12","18","36")
for (m in month) {
   # m="18"
   print(m)
   meta_input_type_sel <- meta_input_type%>%filter(month==m)%>%dplyr::select(-month)
   meta_input_type_sel_sum_1 <- meta_input_type_sel%>%
      dplyr::summarise_at(vars(-st_subjectID),n_distinct,na.rm = T)%>%t()%>%as.data.frame()%>%rownames_to_column("covariate")
   meta_input_type_sel_sum_2 <- meta_input_type_sel%>%
      dplyr::summarise_at(vars(-st_subjectID),~sum(is.na(.)))%>%
      t()%>%as.data.frame()%>%rownames_to_column("covariate")%>%
      mutate(total=nrow(meta_input_type_sel))%>%
      mutate(perc=V1/total)
   names(meta_input_type_sel_sum_1)[2] <- "levels"
   names(meta_input_type_sel_sum_2)[2] <- "NA_number"   
   meta_input_type_sel_sum <- left_join(meta_input_type_sel_sum_2,meta_input_type_sel_sum_1,by="covariate")
   meta_input_type_sel_sum <- meta_input_type_sel_sum%>%
      mutate(reference=case_when(covariate=="Country" ~ "Country,USA",
                                 covariate=="Feed_2" ~ "Feed_2,Non_breastmilk",
                                 covariate=="Birth_weight_2" ~ "Birth_weight_2,NBW",
                                 covariate=="Gestational_age_2" ~ "Gestational_age_2,late_preterm",
                                 covariate=="Delivery" ~ "Delivery,a_vaginal",
                                 covariate=="Gender" ~ "Gender,female",
                                 covariate=="Infant_AB" ~ "Infant_AB,no"))
   meta_input_type_sel_sum_sel <- meta_input_type_sel_sum%>%filter(perc<0.5)%>%filter(levels>=2)
   if (m=="36") {
      meta_input_type_sel_sum_sel <- meta_input_type_sel_sum_sel%>%filter(!covariate%in%c("Country"))#Country is removed as only USA left after filtering the samples with NA for Birth_weight_2.
   }
   preterm_type_t_sel <- preterm_type_t%>%filter(rownames(.)%in%rownames(meta_input_type_sel))%>%dplyr::select(which(colSums(.)!=0))
   identical(rownames(meta_input_type_sel),rownames(preterm_type_t_sel))
   fit_data_type_month = Maaslin2(
      input_data = preterm_type_t_sel, 
      input_metadata = meta_input_type_sel,
      output = paste0('maaslin2/type_prev5/period/',m),
      normalization = 'TSS',
      transform = 'AST',
      min_abundance=0,
      min_prevalence=0.05,
      plot_scatter = F,
      fixed_effects = meta_input_type_sel_sum_sel$covariate,
      random_effects = c('st_subjectID'),
      reference = meta_input_type_sel_sum_sel$reference)
}
sig_month_type_cb <- data.frame()
all_month_type_cb <- data.frame()
for (m in month) {
   # m=1
   print(m)
   sig_month_type <- read.delim(paste0("maaslin2/type_prev5/period/",m,"/significant_results.tsv"))
   sig_month_type$month <- m
   sig_month_type_cb <- rbind(sig_month_type_cb,sig_month_type)
   all_month_type <- read.delim(paste0("maaslin2/type_prev5/period/",m,"/all_results.tsv"))
   all_month_type$month <- m
   all_month_type_cb <- rbind(all_month_type_cb,all_month_type)
}
all_results_type <- read.delim("maaslin2/type_prev5/all/all_results.tsv")
sig_results_type <- read.delim("maaslin2/type_prev5/all/significant_results.tsv")
length(unique(sig_results_type$feature))#24

type_all_period_all <- rbind(all_results_type%>%mutate(month="all"),all_month_type_cb)
type_all_period_sig <- rbind(sig_results_type%>%mutate(month="all"),sig_month_type_cb)
type_all_period_all$feature <- gsub("macrolide.lincosamide.streptogramin","macrolide-lincosamide-streptogramin",type_all_period_all$feature,fixed = T)
type_all_period_sig$feature <- gsub("macrolide.lincosamide.streptogramin","macrolide-lincosamide-streptogramin",type_all_period_sig$feature,fixed = T)

type_all_period_sig$hp <- -log10(type_all_period_sig$qval)*sign(type_all_period_sig$coef)
type_all_period_sig$meta <- paste0(type_all_period_sig$metadata,"_",type_all_period_sig$value,"_",type_all_period_sig$month)
type_all_period_sig$feature <- gsub("macrolide-lincosamide-streptogramin","MLS",type_all_period_sig$feature,fixed = T)
type_all_period_sig$feature <- gsub("other_peptide_antibiotics","OPA",type_all_period_sig$feature,fixed = T)
type_all_period_sig$feature <- gsub("antibacterial_fatty_acid","AFA",type_all_period_sig$feature,fixed = T)

type_all_period_sig_s <- type_all_period_sig%>%
   filter(metadata%in%c("Gestational_age_2","Birth_weight_2","Infant_AB","Feed_2"))%>%
   dplyr::select(feature,meta,hp)%>%spread(.,meta,hp)
type_all_period_sig_s <- type_all_period_sig_s%>%
   dplyr::select(feature,
                 Birth_weight_2_ELBW_all,Birth_weight_2_ELBW_0.25,Birth_weight_2_ELBW_0.5,Birth_weight_2_ELBW_0.75,Birth_weight_2_ELBW_1,
                 Birth_weight_2_ELBW_2,Birth_weight_2_ELBW_12,Birth_weight_2_ELBW_18,Birth_weight_2_ELBW_36,
                 Birth_weight_2_VLBW_all,Birth_weight_2_VLBW_0.25,Birth_weight_2_VLBW_0.5,Birth_weight_2_VLBW_0.75,Birth_weight_2_VLBW_1,
                 Birth_weight_2_VLBW_2,Birth_weight_2_VLBW_12,Birth_weight_2_VLBW_18,Birth_weight_2_VLBW_36,
                 Birth_weight_2_LBW_all,Birth_weight_2_LBW_0.25,Birth_weight_2_LBW_0.5,Birth_weight_2_LBW_2,Birth_weight_2_LBW_3,Birth_weight_2_LBW_36,
                 Gestational_age_2_extremely_preterm_all,Gestational_age_2_extremely_preterm_0.25,Gestational_age_2_extremely_preterm_0.75,
                 Gestational_age_2_extremely_preterm_2,Gestational_age_2_extremely_preterm_12,Gestational_age_2_extremely_preterm_18,Gestational_age_2_extremely_preterm_36,
                 Gestational_age_2_very_preterm_all,Gestational_age_2_very_preterm_0.25,Gestational_age_2_very_preterm_0.75,Gestational_age_2_very_preterm_2,
                 Gestational_age_2_very_preterm_12,Gestational_age_2_very_preterm_18,Gestational_age_2_very_preterm_36,
                 Gestational_age_2_moderate_preterm_all,Gestational_age_2_moderate_preterm_0.25,Gestational_age_2_moderate_preterm_0.75,
                 Gestational_age_2_moderate_preterm_2,Gestational_age_2_moderate_preterm_12,Gestational_age_2_moderate_preterm_18,
                 Infant_AB_yes_all,Infant_AB_yes_0.25,Infant_AB_yes_6,Infant_AB_yes_12,Infant_AB_yes_18,
                 Feed_2_Exclusive_all,Feed_2_Exclusive_6,
                 Feed_2_Mix_all,
                 Feed_2_nofeed_all)%>%
   column_to_rownames("feature")%>%as.matrix()
type_all_period_sig_s[is.na(type_all_period_sig_s)] <- 0

type_all_period_sig <- type_all_period_sig%>%
   mutate(sig=case_when(qval<0.25 & coef >0 ~ 1,qval<0.25 & coef <0 ~ -1,TRUE ~ 0))
type_all_period_sig_s_plus <- type_all_period_sig%>%
   filter(metadata%in%c("Gestational_age_2","Birth_weight_2","Infant_AB","Feed_2"))%>%
   dplyr::select(feature,meta,sig)%>%spread(.,meta,sig)
type_all_period_sig_s_plus <- type_all_period_sig_s_plus%>%
   dplyr::select(feature,
                 Birth_weight_2_ELBW_all,Birth_weight_2_ELBW_0.25,Birth_weight_2_ELBW_0.5,Birth_weight_2_ELBW_0.75,Birth_weight_2_ELBW_1,
                 Birth_weight_2_ELBW_2,Birth_weight_2_ELBW_12,Birth_weight_2_ELBW_18,Birth_weight_2_ELBW_36,
                 Birth_weight_2_VLBW_all,Birth_weight_2_VLBW_0.25,Birth_weight_2_VLBW_0.5,Birth_weight_2_VLBW_0.75,Birth_weight_2_VLBW_1,
                 Birth_weight_2_VLBW_2,Birth_weight_2_VLBW_12,Birth_weight_2_VLBW_18,Birth_weight_2_VLBW_36,
                 Birth_weight_2_LBW_all,Birth_weight_2_LBW_0.25,Birth_weight_2_LBW_0.5,Birth_weight_2_LBW_2,Birth_weight_2_LBW_3,Birth_weight_2_LBW_36,
                 Gestational_age_2_extremely_preterm_all,Gestational_age_2_extremely_preterm_0.25,Gestational_age_2_extremely_preterm_0.75,
                 Gestational_age_2_extremely_preterm_2,Gestational_age_2_extremely_preterm_12,Gestational_age_2_extremely_preterm_18,Gestational_age_2_extremely_preterm_36,
                 Gestational_age_2_very_preterm_all,Gestational_age_2_very_preterm_0.25,Gestational_age_2_very_preterm_0.75,Gestational_age_2_very_preterm_2,
                 Gestational_age_2_very_preterm_12,Gestational_age_2_very_preterm_18,Gestational_age_2_very_preterm_36,
                 Gestational_age_2_moderate_preterm_all,Gestational_age_2_moderate_preterm_0.25,Gestational_age_2_moderate_preterm_0.75,
                 Gestational_age_2_moderate_preterm_2,Gestational_age_2_moderate_preterm_12,Gestational_age_2_moderate_preterm_18,
                 Infant_AB_yes_all,Infant_AB_yes_0.25,Infant_AB_yes_6,Infant_AB_yes_12,Infant_AB_yes_18,
                 Feed_2_Exclusive_all,Feed_2_Exclusive_6,
                 Feed_2_Mix_all,
                 Feed_2_nofeed_all)%>%
   column_to_rownames("feature")%>%as.matrix()
type_all_period_sig_s_plus[is.na(type_all_period_sig_s_plus)] <- 0

hp_split_type <- data.frame(a=c("ELBW_all",rep("ELBW_month",8),"VLBW_all",rep("VLBW_month",8),"LBW_all",rep("LBW_month",5),
                                "extre_all",rep("extre_month",6),"very_all",rep("very_month",6),"moderate_all",rep("moderate_month",5),
                                "ab",rep("ab_month",4),
                                "exclu_all",rep("exclu_month",1),"mix_all","nofeed_all"),
                            b=names(as.data.frame(type_all_period_sig_s_plus)))
hp_split_type$`Timepoints (month)` <- sapply(strsplit(hp_split_type$b, split='_', fixed=TRUE), tail,1)
hp_split_type$`Timepoints (month)` <- factor(hp_split_type$`Timepoints (month)`,levels = c("all","0.25","0.5","0.75","1","2","3","6","12","18","36"))
hp_split_type$Covariates <- sapply(strsplit(hp_split_type$a, split='_', fixed=TRUE), function(x)(x[1]))
hp_split_type$Covariates <- factor(hp_split_type$Covariates,levels = c("ELBW","VLBW","LBW","extre","very","moderate","ab","exclu","mix","nofeed"))
hp_split_type$a <- factor(hp_split_type$a,levels = c("ELBW_all","ELBW_month","VLBW_all","VLBW_month","LBW_all","LBW_month",
                                                     "extre_all","extre_month","very_all","very_month","moderate_all","moderate_month",
                                                     "ab","ab_month",
                                                     "exclu_all","exclu_month","mix_all","nofeed_all"))

palette_type = brewer.pal(11,'RdBu')
col_fun_type <- colorRamp2(c(-10,-5,-2,0,2,5,10),c("#066A9B","#1E9FDE","#5AA4C8","white","#D6604D","#E6391D","#A71B05"))

hp_type_top_anno <- 
   HeatmapAnnotation(
      Covariates=hp_split_type%>%pull(Covariates),
      `Timepoints (month)`=hp_split_type%>%pull(`Timepoints (month)`),
      annotation_legend_param =list(nrow = 1),
      col=list("Covariates"=c("ELBW"="#066A9B","VLBW"="#54A5CD","LBW"="#9ED3EE",
                              "extre"="#098D3B","very"="#4CC279","moderate"="#A3F5C2",
                              "ab"="#B95912",
                              "exclu"="#F0E2A1","mix"="#D2BB4B","nofeed"="#8C7505"),
               `Timepoints (month)`=c("all"="#252525","0.25"="#556997","0.5"="#8286A3","0.75"="#ACA4AE","1"="#D5C3B9",
                                      "2"="#FFE4C4","3"="#FBC39D","6"="#F4A277","12"="#E98152","18"="#DC5F2D","36"="#CD3700")),
      simple_anno_size = unit(0.3, "cm"), height = unit(1.6, "cm"),
      annotation_name_gp=gpar(fontsize = 10))
#Figure 4f
ht.main_type_sig <- Heatmap(type_all_period_sig_s,cluster_rows=T,cluster_columns=F,col = col_fun_type,
                            show_row_names = T,show_column_names = F,row_names_side = "right",row_dend_side="left",
                            top_annotation = hp_type_top_anno,
                            column_split =hp_split_type$a, 
                            border = T,row_names_gp = gpar(fontsize =12),
                            border_gp=gpar(col = "black", lwd = 0.7),column_gap = unit(2, "mm"),
                            width = unit(24, "cm"), height = unit(11, "cm"),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = list(title = "coef",direction = "horizontal", legend_hight = unit(1, "cm"),legend_width = unit(3, "cm")),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                               if(type_all_period_sig_s_plus[i,j]>0) {
                                  grid.text("+", x = x, y = y,gp = gpar(fontsize = 8))
                               } else if(type_all_period_sig_s_plus[i,j]<0) {
                                  grid.text("-", x = x, y = y,gp = gpar(fontsize = 8))
                               } else {
                                  grid.text("", x = x, y = y,gp = gpar(fontsize = 8))
                               }
                            }
)




