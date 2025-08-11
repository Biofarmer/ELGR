library(readxl)
library(dplyr)
library(tibble)
library(vegan)
library(coin)
library(ggplot2)
library(tidyr)
library(stringr)


##########################################################################################################################################
#alpha and beta diversity between nec and control
##########################################################################################################################################
metadata_nec_sel <- read.csv("data/metadata_nec_sel_new.csv")#4180 runs
metadata_nec_sel <- metadata_nec_sel%>%
   mutate(Feed_2=case_when(Feed=="breast"~"Exclusive",
                           Feed=="combined"~"Mix",
                           Feed=="nofeed"~"fasting",
                           Feed=="formula" | Feed =="solid"~"Non-breastmilk"))%>%
   mutate(nec=case_when(PreNEC!="noNEC"~"yes",
                        PreNEC=="noNEC"~"no"))%>%
   mutate(nec2=case_when(PreNEC=="preNEC"~"preNEC",
                         PreNEC=="onsetNEC"~"postNEC",
                         PreNEC=="postNEC"~"postNEC",
                         PreNEC=="noNEC"~"no"))%>%
   mutate(st_subject=paste0(Study,"_",SubjectID))

########################################################################################################################################
#amr
########################################################################################################################################
amr_preterm_profile <- read.csv("data/amr_preterm_resistome_coverm_bins_gtdb_stat_plasmid.csv")
amr_preterm_profile_coverm <- amr_preterm_profile%>%filter(sample%in%metadata_nec_sel$SequenceID)
amr_preterm_profile_coverm_gene_sum <- amr_preterm_profile_coverm%>%
   group_by(sample,Element_symbol)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%spread(.,sample,rpkm_total)
amr_preterm_profile_coverm_gene_sum[is.na(amr_preterm_profile_coverm_gene_sum)] <- 0
amr_preterm_profile_coverm_gene_sum_t <- amr_preterm_profile_coverm_gene_sum%>%column_to_rownames("Element_symbol")%>%t()%>%as.data.frame()
#richness is based amr regardless of the abundance from coverm
amr_preterm_profile_coverm_gene_sum2 <- amr_preterm_profile_coverm%>%
   group_by(sample,Element_symbol)%>%dplyr::summarise(rpkm_total=sum(rpkm))
amr_richness <- amr_preterm_profile_coverm_gene_sum2%>%group_by(sample)%>%dplyr::summarise(n=n())
names(amr_richness) <- c("SequenceID", "Richness")
amr_shannon <- diversity(amr_preterm_profile_coverm_gene_sum_t, index = "shannon", MARGIN=1) 
amr_shannon_df <- amr_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(amr_shannon_df) <- c("SequenceID", "Shannon")
amr_alpha <- left_join(amr_richness,amr_shannon_df,by="SequenceID")
amr_alpha_meta <- left_join(metadata_nec_sel,amr_alpha,by="SequenceID")
amr_alpha_meta <- amr_alpha_meta%>%filter(!is.na(Richness))
amr_alpha_meta%>%group_by(nec2)%>%dplyr::summarise(mean=mean(Shannon),median=median(Shannon))
alpha_amr_pval <- data.frame("index"="x","comp1"="x","comp2"="x","pval_wil"=0,"pval_lm"=0,stringsAsFactors = F)[-1,]
com <- list(c("no","preNEC"),c("no","postNEC"),c("preNEC","postNEC"))
for (c in seq(1:3)) {
   a <- amr_alpha_meta%>%filter(nec2%in%com[[c]])
   a$nec2 <- as.character(a$nec2)
   a$Study <- as.factor(a$Study)
   a$nec2 <- as.factor(a$nec2)
   #richness
   a_sum <- a%>%group_by(nec2)%>%dplyr::summarise(mean=mean(Richness),median=median(Richness))
   pval_wil <- pvalue(wilcox_test(Richness ~ nec2 | Study, data = a))
   lm_fit <- lmerTest::lmer(Richness ~ nec2+(1|Study)+(1|st_subject),data = a)
   pval_lm <- anova(lm_fit)
   alpha_amr_pval[nrow(alpha_amr_pval)+1,] <- c("Richness",com[[c]],pval_wil,pval_lm$`Pr(>F)`)
   #shannon
   a_sum <- a%>%group_by(nec2)%>%dplyr::summarise(mean=mean(Shannon),median=median(Shannon))
   pval_wil <- pvalue(wilcox_test(Shannon ~ nec2 | Study, data = a))
   lm_fit <- lmerTest::lmer(Shannon ~ nec2+(1|Study)+(1|st_subject),data = a)
   pval_lm <- anova(lm_fit)
   alpha_amr_pval[nrow(alpha_amr_pval)+1,] <- c("Shannon",com[[c]],pval_wil,pval_lm$`Pr(>F)`)
}
alpha_amr_pval <- alpha_amr_pval%>%mutate_at(vars(-index,-comp1,-comp2),as.numeric)
amr_alpha_meta_g <- amr_alpha_meta%>%dplyr::select(Richness,Shannon,SequenceID,nec2)%>%gather(.,group,value,-SequenceID,-nec2)
amr_alpha_meta_g$group[amr_alpha_meta_g$group=="Shannon"] <- "Shannon index"
amr_alpha_meta_g$nec2 <- factor(amr_alpha_meta_g$nec2,levels = c("no","preNEC","postNEC"))
#Supplementary Figure 11b
amr_gene_alpha_p <- ggplot(data=amr_alpha_meta_g,aes(x=nec2, y=value, col=nec2)) +
   geom_boxplot(outlier.size = 0.1) + 
   stat_summary(fun.y=mean, geom="point",size=2, color="black")+
   facet_wrap(~group,scales = "free")+
   scale_color_manual(values = c("#4994C4","#E04B37","#2ABF2A"),labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   scale_x_discrete(labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   theme_bw()+ 
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   labs(x="",y="Alpha diversity")+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         strip.text = element_text(size=12, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 45,hjust = 1),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))

#beta diversity
#no NEC vs. pre-NEC
metadata_nec_sel_no_pre <- metadata_nec_sel%>%filter(nec2%in%c("no","preNEC"))
amr_t_no_pre <- amr_preterm_profile_coverm_gene_sum_t%>%filter(rownames(.)%in%metadata_nec_sel_no_pre$SequenceID)
amr_t_no_pre <- amr_t_no_pre[which(rowSums(amr_t_no_pre)!=0),which(colSums(amr_t_no_pre)!=0)]
amr_t_no_pre_dist_beta <- vegdist(amr_t_no_pre, method = "bray")
amr_t_no_pre_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(amr_t_no_pre))%>%
   arrange(match(SequenceID,rownames(amr_t_no_pre)))
identical(rownames(amr_t_no_pre),amr_t_no_pre_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(amr_t_no_pre_dist_beta_meta, vegan::adonis2(amr_t_no_pre_dist_beta ~ nec2,data=amr_t_no_pre_dist_beta_meta,permutations = 1000, strata = Study))

#no NEC vs. post-NEC
metadata_nec_sel_no_post <- metadata_nec_sel%>%filter(nec2%in%c("no","postNEC"))
amr_t_no_post <- amr_preterm_profile_coverm_gene_sum_t%>%filter(rownames(.)%in%metadata_nec_sel_no_post$SequenceID)
amr_t_no_post <- amr_t_no_post[which(rowSums(amr_t_no_post)!=0),which(colSums(amr_t_no_post)!=0)]
amr_t_no_post_dist_beta <- vegdist(amr_t_no_post, method = "bray")
amr_t_no_post_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(amr_t_no_post))%>%
   arrange(match(SequenceID,rownames(amr_t_no_post)))
identical(rownames(amr_t_no_post),amr_t_no_post_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(amr_t_no_post_dist_beta_meta, vegan::adonis2(amr_t_no_post_dist_beta ~ nec2,data=amr_t_no_post_dist_beta_meta,permutations = 1000, strata = Study))
#pre-NEC vs. post-NEC
metadata_nec_sel_pre_post <- metadata_nec_sel%>%filter(nec2%in%c("preNEC","postNEC"))
amr_t_pre_post <- amr_preterm_profile_coverm_gene_sum_t%>%filter(rownames(.)%in%metadata_nec_sel_pre_post$SequenceID)
amr_t_pre_post <- amr_t_pre_post[which(rowSums(amr_t_pre_post)!=0),which(colSums(amr_t_pre_post)!=0)]
amr_t_pre_post_dist_beta <- vegdist(amr_t_pre_post, method = "bray")
amr_t_pre_post_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(amr_t_pre_post))%>%
   arrange(match(SequenceID,rownames(amr_t_pre_post)))
identical(rownames(amr_t_pre_post),amr_t_pre_post_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(amr_t_pre_post_dist_beta_meta, vegan::adonis2(amr_t_pre_post_dist_beta ~ nec2,data=amr_t_pre_post_dist_beta_meta,permutations = 1000, strata = Study))
amr_t <- amr_preterm_profile_coverm_gene_sum_t[which(rowSums(amr_preterm_profile_coverm_gene_sum_t)!=0),which(colSums(amr_preterm_profile_coverm_gene_sum_t)!=0)]
amr_t_dist_beta <- vegdist(amr_t, method = "bray")
amr_t_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(amr_t))%>%
   arrange(match(SequenceID,rownames(amr_t)))
identical(rownames(amr_t),amr_t_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(amr_t_dist_beta_meta, vegan::adonis2(amr_t_dist_beta ~ nec2,data=amr_t_dist_beta_meta,permutations = 1000, strata = Study))
amr_beta_pco <- cmdscale(amr_t_dist_beta, k=2, eig = T)
amr_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                   round((amr_beta_pco$eig[1]/sum(amr_beta_pco$eig))*100,1),
                                   '%]', sep='')
amr_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                   round((amr_beta_pco$eig[2]/sum(amr_beta_pco$eig))*100,1),
                                   '%]', sep='')
amr_beta_pco.point <- as.data.frame(amr_beta_pco$points)%>%rownames_to_column('SequenceID')
amr_beta_pco.tabble <- right_join(amr_t_dist_beta_meta, amr_beta_pco.point, by="SequenceID")
amr_beta_df.plot <- tibble(Axis1 = amr_beta_pco.tabble$V1,
                           Axis2 = amr_beta_pco.tabble$V2,
                           Sample_ID = amr_beta_pco.tabble$SubjectID,
                           nec2=amr_beta_pco.tabble$nec2)
amr_beta_df.plot$nec2 <- factor(amr_beta_df.plot$nec2,levels = c("no","preNEC","postNEC"))
#Supplementary Figure 11b
amr_gene_beta_time_p <- ggplot(data=amr_beta_df.plot,aes(x=Axis1, y=Axis2, col=nec2)) +
   geom_point(size=1.5, alpha=1) + 
   stat_ellipse(aes(fill=nec2),geom = "polygon",level = 0.95, alpha=0.05)+
   scale_color_manual(values = c("#4994C4","#E04B37","#2ABF2A"),labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   theme_bw()+ 
   xlab(amr_beta_pco_axis.1.title) + ylab(amr_beta_pco_axis.2.title) +
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
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))

#subtype
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_resistome_nec <- preterm_resistome%>%select(names,group,metadata_nec_sel$SequenceID)
preterm_type <- preterm_resistome_nec%>%filter(group=="type")%>%dplyr::select(-group)
type_sum_abund <- preterm_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("SequenceID"=sample,"type_total_abund"=rowSums(.[-1]))
preterm_subtype <- preterm_resistome_nec%>%filter(group=="subtype")%>%dplyr::select(-group)
preterm_subtype_t <- preterm_subtype%>%column_to_rownames("names")%>%t()%>%
   as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
subtype_richness <- rowSums(preterm_subtype_t != 0)
subtype_shannon <- diversity(preterm_subtype_t, index = "shannon", MARGIN=1) 
subtype_richness_df <- subtype_richness%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(subtype_richness_df)[2] <- "Richness"
subtype_shannon_df <- subtype_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(subtype_shannon_df)[2] <- "Shannon"
subtype_alpha <- left_join(subtype_richness_df,subtype_shannon_df,by="SequenceID")
subtype_alpha_meta <- left_join(metadata_nec_sel,subtype_alpha,by="SequenceID")
alpha_subtype_pval <- data.frame("index"="x","comp1"="x","comp2"="x","pval_wil"=0,"pval_lm"=0,stringsAsFactors = F)[-1,]
com <- list(c("no","preNEC"),c("no","postNEC"),c("preNEC","postNEC"))
for (c in seq(1:3)) {
   a <- subtype_alpha_meta%>%filter(nec2%in%com[[c]])
   a$nec2 <- as.character(a$nec2)
   a$Study <- as.factor(a$Study)
   a$nec2 <- as.factor(a$nec2)
   #richness
   a_sum <- a%>%group_by(nec2)%>%dplyr::summarise(mean=mean(Richness),median=median(Richness))
   pval_wil <- pvalue(wilcox_test(Richness ~ nec2 | Study, data = a))
   lm_fit <- lmerTest::lmer(Richness ~ nec2+(1|Study)+(1|st_subject),data = a)
   pval_lm <- anova(lm_fit)
   alpha_subtype_pval[nrow(alpha_subtype_pval)+1,] <- c("Richness",com[[c]],pval_wil,pval_lm$`Pr(>F)`)
   #shannon
   a_sum <- a%>%group_by(nec2)%>%dplyr::summarise(mean=mean(Shannon),median=median(Shannon))
   pval_wil <- pvalue(wilcox_test(Shannon ~ nec2 | Study, data = a))
   lm_fit <- lmerTest::lmer(Shannon ~ nec2+(1|Study)+(1|st_subject),data = a)
   pval_lm <- anova(lm_fit)
   alpha_subtype_pval[nrow(alpha_subtype_pval)+1,] <- c("Shannon",com[[c]],pval_wil,pval_lm$`Pr(>F)`)
}
alpha_subtype_pval <- alpha_subtype_pval%>%mutate_at(vars(-index,-comp1,-comp2),as.numeric)
subtype_alpha_meta_g <- subtype_alpha_meta%>%dplyr::select(Richness,Shannon,SequenceID,nec2)%>%gather(.,group,value,-SequenceID,-nec2)
subtype_alpha_meta_g$group[subtype_alpha_meta_g$group=="Shannon"] <- "Shannon index"
subtype_alpha_meta_g$nec2 <- factor(subtype_alpha_meta_g$nec2,levels = c("no","preNEC","postNEC"))
#Supplementary Figure 11a
subtype_alpha_p <- ggplot(data=subtype_alpha_meta_g,aes(x=nec2, y=value, col=nec2)) +
   geom_boxplot(outlier.size = 0.1) + 
   stat_summary(fun.y=mean, geom="point",size=2, color="black")+
   facet_wrap(~group,scales = "free")+
   scale_color_manual(values = c("#4994C4","#E04B37","#2ABF2A"),labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   scale_x_discrete(labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   theme_bw()+ 
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   labs(x="",y="Alpha diversity")+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         strip.text = element_text(size=12, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 45,hjust = 1),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))

#beta diversity
metadata_nec_sel_no_pre <- metadata_nec_sel%>%filter(nec2%in%c("no","preNEC"))
subtype_t_no_pre <- preterm_subtype_t%>%filter(rownames(.)%in%metadata_nec_sel_no_pre$SequenceID)
subtype_t_no_pre <- subtype_t_no_pre[which(rowSums(subtype_t_no_pre)!=0),which(colSums(subtype_t_no_pre)!=0)]
subtype_t_no_pre_dist_beta <- vegdist(subtype_t_no_pre, method = "bray")
subtype_t_no_pre_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(subtype_t_no_pre))%>%
   arrange(match(SequenceID,rownames(subtype_t_no_pre)))
identical(rownames(subtype_t_no_pre),subtype_t_no_pre_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(subtype_t_no_pre_dist_beta_meta, vegan::adonis2(subtype_t_no_pre_dist_beta ~ nec2,data=subtype_t_no_pre_dist_beta_meta,permutations = 1000, strata = Study))
metadata_nec_sel_no_post <- metadata_nec_sel%>%filter(nec2%in%c("no","postNEC"))
subtype_t_no_post <- preterm_subtype_t%>%filter(rownames(.)%in%metadata_nec_sel_no_post$SequenceID)
subtype_t_no_post <- subtype_t_no_post[which(rowSums(subtype_t_no_post)!=0),which(colSums(subtype_t_no_post)!=0)]
subtype_t_no_post_dist_beta <- vegdist(subtype_t_no_post, method = "bray")
subtype_t_no_post_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(subtype_t_no_post))%>%
   arrange(match(SequenceID,rownames(subtype_t_no_post)))
identical(rownames(subtype_t_no_post),subtype_t_no_post_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(subtype_t_no_post_dist_beta_meta, vegan::adonis2(subtype_t_no_post_dist_beta ~ nec2,data=subtype_t_no_post_dist_beta_meta,permutations = 1000, strata = Study))
#pre-NEC vs. post-NEC
metadata_nec_sel_pre_post <- metadata_nec_sel%>%filter(nec2%in%c("preNEC","postNEC"))
subtype_t_pre_post <- preterm_subtype_t%>%filter(rownames(.)%in%metadata_nec_sel_pre_post$SequenceID)
subtype_t_pre_post <- subtype_t_pre_post[which(rowSums(subtype_t_pre_post)!=0),which(colSums(subtype_t_pre_post)!=0)]
subtype_t_pre_post_dist_beta <- vegdist(subtype_t_pre_post, method = "bray")
subtype_t_pre_post_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(subtype_t_pre_post))%>%
   arrange(match(SequenceID,rownames(subtype_t_pre_post)))
identical(rownames(subtype_t_pre_post),subtype_t_pre_post_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(subtype_t_pre_post_dist_beta_meta, vegan::adonis2(subtype_t_pre_post_dist_beta ~ nec2,data=subtype_t_pre_post_dist_beta_meta,permutations = 1000, strata = Study))
preterm_subtype_beta_dist <- vegdist(preterm_subtype_t, method = "bray")
subtype_alpha_meta_order <- subtype_alpha_meta%>%arrange(match(SequenceID,rownames(preterm_subtype_t)))
identical(rownames(preterm_subtype_t),subtype_alpha_meta_order$SequenceID)#TRUE
set.seed(1000)
with(subtype_alpha_meta_order, vegan::adonis2(preterm_subtype_beta_dist ~ nec2,data=subtype_alpha_meta_order,permutations = 1000, strata = Study))
subtype_beta_pco <- cmdscale(preterm_subtype_beta_dist, k=2, eig = T)
subtype_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                       round((subtype_beta_pco$eig[1]/sum(subtype_beta_pco$eig))*100,1),
                                       '%]', sep='')
subtype_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                       round((subtype_beta_pco$eig[2]/sum(subtype_beta_pco$eig))*100,1),
                                       '%]', sep='')
subtype_beta_pco.point <- as.data.frame(subtype_beta_pco$points)%>%rownames_to_column('SequenceID')
subtype_beta_pco.tabble <- right_join(subtype_alpha_meta_order, subtype_beta_pco.point, by="SequenceID")
subtype_beta_df.plot <- tibble(Axis1 = subtype_beta_pco.tabble$V1,
                               Axis2 = subtype_beta_pco.tabble$V2,
                               Sample_ID = subtype_beta_pco.tabble$SubjectID,
                               nec2=subtype_beta_pco.tabble$nec2)
subtype_beta_df.plot$nec2 <- factor(subtype_beta_df.plot$nec2,levels = c("no","preNEC","postNEC"))
#Supplementary Figure 11a
subtype_beta_time_p <- ggplot(data=subtype_beta_df.plot,aes(x=Axis1, y=Axis2, col=nec2)) +
   geom_point(size=1.5, alpha=1) + 
   stat_ellipse(aes(fill=nec2),geom = "polygon",level = 0.95, alpha=0.05)+
   scale_color_manual(values = c("#4994C4","#E04B37","#2ABF2A"),labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   theme_bw()+ 
   xlab(subtype_beta_pco_axis.1.title) + ylab(subtype_beta_pco_axis.2.title) +
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
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))

#species from metaphlan4
mpl_preterm <- read.csv("data/metaphlan4_preterm_new.csv")
mpl_term <- read.csv("data/metaphlan4_term_new.csv")
mpl_cb <- left_join(mpl_preterm,mpl_term,by="clade_name")
mpl_cb[is.na(mpl_cb)] <- 0
mpl_nec <- mpl_cb%>%dplyr::select(clade_name,metadata_nec_sel$SequenceID)
species_level <- mpl_nec%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
species_level$species <- sapply(strsplit(species_level$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
species_level_1 <- species_level%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
species_level_richness <- rowSums(species_level_1!=0)%>%as.data.frame()%>%setNames("Richness")%>%rownames_to_column("run")
species_level_shannon <- vegan::diversity(species_level_1, index = "shannon", MARGIN=1)%>%as.data.frame()%>%setNames("Shannon")%>%rownames_to_column("run")
species_alpha <- left_join(species_level_richness,species_level_shannon,by="run")
species_alpha_meta <- left_join(metadata_nec_sel,species_alpha,by=c("SequenceID"="run"))
alpha_species_pval <- data.frame("index"="x","comp1"="x","comp2"="x","pval_wil"=0,"pval_lm"=0,stringsAsFactors = F)[-1,]
com <- list(c("no","preNEC"),c("no","postNEC"),c("preNEC","postNEC"))
for (c in seq(1:3)) {
   a <- species_alpha_meta%>%filter(nec2%in%com[[c]])
   a$nec2 <- as.character(a$nec2)
   a$Study <- as.factor(a$Study)
   a$nec2 <- as.factor(a$nec2)
   #richness
   a_sum <- a%>%group_by(nec2)%>%dplyr::summarise(mean=mean(Richness),median=median(Richness))
   pval_wil <- pvalue(wilcox_test(Richness ~ nec2 | Study, data = a))
   lm_fit <- lmerTest::lmer(Richness ~ nec2+(1|Study)+(1|st_subject),data = a)
   pval_lm <- anova(lm_fit)
   alpha_species_pval[nrow(alpha_species_pval)+1,] <- c("Richness",com[[c]],pval_wil,pval_lm$`Pr(>F)`)
   #shannon
   a_sum <- a%>%group_by(nec2)%>%dplyr::summarise(mean=mean(Shannon),median=median(Shannon))
   pval_wil <- pvalue(wilcox_test(Shannon ~ nec2 | Study, data = a))
   lm_fit <- lmerTest::lmer(Shannon ~ nec2+(1|Study)+(1|st_subject),data = a)
   pval_lm <- anova(lm_fit)
   alpha_species_pval[nrow(alpha_species_pval)+1,] <- c("Shannon",com[[c]],pval_wil,pval_lm$`Pr(>F)`)
}
alpha_species_pval <- alpha_species_pval%>%mutate_at(vars(-index,-comp1,-comp2),as.numeric)
species_alpha_meta_g <- species_alpha_meta%>%dplyr::select(Richness,Shannon,SequenceID,nec2)%>%gather(.,group,value,-SequenceID,-nec2)
species_alpha_meta_g$group[species_alpha_meta_g$group=="Shannon"] <- "Shannon index"
species_alpha_meta_g$nec2 <- factor(species_alpha_meta_g$nec2,levels = c("no","preNEC","postNEC"))
#Supplementary Figure 11c
species_alpha_p <- ggplot(data=species_alpha_meta_g,aes(x=nec2, y=value, col=nec2)) +
   geom_boxplot(outlier.size = 0.1) + 
   stat_summary(fun.y=mean, geom="point",size=2, color="black")+
   facet_wrap(~group,scales = "free")+
   scale_color_manual(values = c("#4994C4","#E04B37","#2ABF2A"),labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   scale_x_discrete(labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   theme_bw()+ 
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   labs(x="",y="Alpha diversity")+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         strip.text = element_text(size=12, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 45,hjust = 1),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))
#beta diversity
metadata_nec_sel_no_pre <- metadata_nec_sel%>%filter(nec2%in%c("no","preNEC"))
species_t_no_pre <- species_level_1%>%filter(rownames(.)%in%metadata_nec_sel_no_pre$SequenceID)
species_t_no_pre <- species_t_no_pre[which(rowSums(species_t_no_pre)!=0),which(colSums(species_t_no_pre)!=0)]
species_t_no_pre_dist_beta <- vegdist(species_t_no_pre, method = "bray")
species_t_no_pre_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(species_t_no_pre))%>%
   arrange(match(SequenceID,rownames(species_t_no_pre)))
identical(rownames(species_t_no_pre),species_t_no_pre_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(species_t_no_pre_dist_beta_meta, vegan::adonis2(species_t_no_pre_dist_beta ~ nec2,data=species_t_no_pre_dist_beta_meta,permutations = 1000, strata = Study))
metadata_nec_sel_no_post <- metadata_nec_sel%>%filter(nec2%in%c("no","postNEC"))
species_t_no_post <- species_level_1%>%filter(rownames(.)%in%metadata_nec_sel_no_post$SequenceID)
species_t_no_post <- species_t_no_post[which(rowSums(species_t_no_post)!=0),which(colSums(species_t_no_post)!=0)]
species_t_no_post_dist_beta <- vegdist(species_t_no_post, method = "bray")
species_t_no_post_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(species_t_no_post))%>%
   arrange(match(SequenceID,rownames(species_t_no_post)))
identical(rownames(species_t_no_post),species_t_no_post_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(species_t_no_post_dist_beta_meta, vegan::adonis2(species_t_no_post_dist_beta ~ nec2,data=species_t_no_post_dist_beta_meta,permutations = 1000, strata = Study))
metadata_nec_sel_pre_post <- metadata_nec_sel%>%filter(nec2%in%c("preNEC","postNEC"))
species_t_pre_post <- species_level_1%>%filter(rownames(.)%in%metadata_nec_sel_pre_post$SequenceID)
species_t_pre_post <- species_t_pre_post[which(rowSums(species_t_pre_post)!=0),which(colSums(species_t_pre_post)!=0)]
species_t_pre_post_dist_beta <- vegdist(species_t_pre_post, method = "bray")
species_t_pre_post_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(species_t_pre_post))%>%
   arrange(match(SequenceID,rownames(species_t_pre_post)))
identical(rownames(species_t_pre_post),species_t_pre_post_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(species_t_pre_post_dist_beta_meta, vegan::adonis2(species_t_pre_post_dist_beta ~ nec2,data=species_t_pre_post_dist_beta_meta,permutations = 1000, strata = Study))
species_level_2 <- species_level_1[which(rowSums(species_level_1)!=0),which(colSums(species_level_1)!=0)]
species_alpha_meta_order <- species_alpha_meta%>%filter(SequenceID%in%rownames(species_level_2))%>%arrange(match(SequenceID,rownames(species_level_1)))
identical(rownames(species_level_2),species_alpha_meta_order$SequenceID)#TRUE
species_dist_beta <- vegdist(species_level_2, method = "bray")
set.seed(1000)
with(species_alpha_meta_order, vegan::adonis2(species_dist_beta ~ nec2,data=species_alpha_meta_order,permutations = 1000, strata = Study))
species_beta_pco <- cmdscale(species_dist_beta, k=2, eig = T)
species_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                       round((species_beta_pco$eig[1]/sum(species_beta_pco$eig))*100,1),
                                       '%]', sep='')
species_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                       round((species_beta_pco$eig[2]/sum(species_beta_pco$eig))*100,1),
                                       '%]', sep='')
species_beta_pco.point <- as.data.frame(species_beta_pco$points)%>%rownames_to_column('SequenceID')
species_beta_pco.tabble <- right_join(species_alpha_meta_order, species_beta_pco.point, by="SequenceID")
species_beta_df.plot <- tibble(Axis1 = species_beta_pco.tabble$V1,
                               Axis2 = species_beta_pco.tabble$V2,
                               Sample_ID = species_beta_pco.tabble$SubjectID,
                               nec2=species_beta_pco.tabble$nec2)
species_beta_df.plot$nec2 <- factor(species_beta_df.plot$nec2,levels = c("no","preNEC","postNEC"))
#Supplementary Figure 11c
species_beta_time_p <- ggplot(data=species_beta_df.plot,aes(x=Axis1, y=Axis2, col=nec2)) +
   geom_point(size=1.5, alpha=1) + 
   stat_ellipse(aes(fill=nec2),geom = "polygon",level = 0.95, alpha=0.05)+
   scale_color_manual(values = c("#4994C4","#E04B37","#2ABF2A"),labels=c("no"="Control","preNEC"="pre-NEC","postNEC"="post-NEC"))+
   theme_bw()+ 
   xlab(species_beta_pco_axis.1.title) + ylab(species_beta_pco_axis.2.title) +
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))
#########################################################################################################################################
#Accelerated convergence of gut resistome preceding NEC onset
#########################################################################################################################################
metadata_nec_sel <- read.csv("data/metadata_nec_sel_new.csv")#4180 runs
metadata_nec_sel <- metadata_nec_sel%>%
   mutate(Feed_2=case_when(Feed=="breast"~"Exclusive",
                           Feed=="combined"~"Mix",
                           Feed=="formula" | Feed =="solid"~"Non-breastmilk"))%>%
   mutate(nec=case_when(PreNEC!="noNEC"~"yes",
                        PreNEC=="noNEC"~"no"))%>%
   mutate(st_subject=paste0(Study,"_",SubjectID))

closest_preNEC <- read_excel("data/closest preNEC.xlsx")
#for infants with two times of NEC from lou 2024, only the first time is analyzed
closest_preNEC <- closest_preNEC%>%filter(!preNEC_closest_DOL%in%c("NA","nopre"))%>%filter(start!=53)%>%
   mutate(st_subject=paste0(Study,"_",infant))%>%filter(st_subject%in%metadata_nec_sel$st_subject)
closest_preNEC$diff <- closest_preNEC$start- as.numeric(closest_preNEC$preNEC_closest_DOL)
metadata_prenec <- left_join(metadata_nec_sel,closest_preNEC%>%dplyr::select(st_subject,start),by="st_subject")
metadata_prenec$Gestational_age <- round(metadata_prenec$Gestational_age,2)

metadata_prenec_yes <- metadata_prenec%>%filter(nec=="yes")%>%mutate(preday=start-DOL)
metadata_prenec_no <- metadata_prenec%>%filter(nec=="no")
metadata_prenec_yes_sum <- metadata_prenec_yes%>%group_by(preday,st_subject)%>%dplyr::summarise(n=n())
table(metadata_prenec_yes$preday)#0 97

#in each window, we only keep one sample from each case, and one sample from each control that matched with each case.
window_list <-list(c(0,3),c(4,6),c(7,9),c(10,12),c(13,15),c(16,18),c(19,21),c(22,24),c(25,27),c(28,30),c(31,40),c(41,60))
#part 1
window_db_1 <- data.frame()
for (w in c(1,7,8,10,11,12)) {
   # w=1
   print(window_list[[w]])
   pre_yes <- metadata_prenec_yes%>%filter(preday>=window_list[[w]][1] & preday<=window_list[[w]][2])
   pre_yes_1 <- pre_yes%>%
      group_by(st_subject)%>%top_n(-1,preday)%>%distinct(st_subject,.keep_all = T)%>%
      mutate(window=window_list[[w]][1])
   window_db_no <- data.frame()
   for (s in unique(pre_yes$st_subject)) {
      # s="This_study_73"
      a <- pre_yes%>%filter(st_subject==s)%>%top_n(-1,preday)%>%mutate(window=window_list[[w]][1])%>%distinct(DOL,.keep_all = T)
      b <- metadata_prenec_no%>%filter(Study==a$Study)
      b_1 <- b%>%filter(Gestational_age>=a$Gestational_age-0.5 & Gestational_age<=a$Gestational_age+0.5)%>%
         filter(DOL>=a$DOL-1 & DOL<=a$DOL+1)%>%
         group_by(SubjectID)%>%dplyr::top_n(1,DOL)%>%mutate(preday=NA)%>%mutate(window=window_list[[w]][1])%>%
         mutate(case=s)
      window_db_no <- rbind.fill(window_db_no,b_1)
   }
   window_db_no_dup <- window_db_no%>%distinct(SequenceID,.keep_all = T)
   print(paste0("case number: ", length(unique(window_db_no_dup$case))))
   print(paste0("control number: ", length(unique(window_db_no_dup$st_subject))))
   #to make sure the
   window_db_no_dup_1 <- window_db_no_dup
   window_db_no_unique_cb <- data.frame()
   window_db_no_sum_1 <- window_db_no_dup[1,]
   window_db_no_sum_2 <- window_db_no_dup[1,]
   while (nrow(window_db_no_sum_1)!=0 & nrow(window_db_no_sum_2)!=0) {
      window_db_no_sum <- window_db_no_dup_1%>%group_by(st_subject)%>%dplyr::summarise(n=n())
      window_db_no_sum_1 <- window_db_no_sum%>%filter(n==1)
      window_db_no_sum_2 <- window_db_no_sum%>%filter(n>1)
      window_db_no_unique_1 <- window_db_no_dup_1%>%filter(st_subject%in%window_db_no_sum_1$st_subject)
      window_db_no_unique_cb <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_1)
      window_db_no_dup_1 <- window_db_no_dup_1%>%filter(st_subject%in%window_db_no_sum_2$st_subject)%>%filter(!case%in%window_db_no_unique_1$case)
   }
   print(paste0("window_db_no_sum_2: ",nrow(window_db_no_sum_2)))
   print(paste0("case number: ", length(unique(window_db_no_unique_cb$case))))
   window_db_no_unique_left <- window_db_no_dup%>%filter(!st_subject%in%window_db_no_unique_cb$st_subject)%>%distinct(st_subject,.keep_all = T)
   window_db_no_unique <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_left)
   print(paste0("control number: ", length(unique(window_db_no_unique$st_subject))))
   # window_db_no_unique2 <- window_db_no%>%distinct(st_subject,.keep_all = T)#some samples are selected by >1 times
   pre_yes_2 <- rbind.fill(pre_yes_1,window_db_no_unique)
   window_db_1 <- rbind.fill(window_db_1,pre_yes_2)
}
#part 2
window_db_2 <- data.frame()
for (w in c(2,3,4,5,6,9)) {
   # w=6
   print(window_list[[w]])
   pre_yes <- metadata_prenec_yes%>%filter(preday>=window_list[[w]][1] & preday<=window_list[[w]][2])
   pre_yes_1 <- pre_yes%>%
      group_by(st_subject)%>%top_n(-1,preday)%>%distinct(st_subject,.keep_all = T)%>%
      mutate(window=window_list[[w]][1])
   window_db_no <- data.frame()
   for (s in unique(pre_yes$st_subject)) {
      # s="This_study_73"
      a <- pre_yes%>%filter(st_subject==s)%>%top_n(-1,preday)%>%mutate(window=window_list[[w]][1])%>%distinct(DOL,.keep_all = T)
      b <- metadata_prenec_no%>%filter(Study==a$Study)
      b_1 <- b%>%filter(Gestational_age>=a$Gestational_age-0.5 & Gestational_age<=a$Gestational_age+0.5)%>%
         filter(DOL>=a$DOL-1 & DOL<=a$DOL+1)%>%
         group_by(SubjectID)%>%dplyr::top_n(1,DOL)%>%mutate(preday=NA)%>%mutate(window=window_list[[w]][1])%>%
         mutate(case=s)
      window_db_no <- rbind.fill(window_db_no,b_1)
   }
   window_db_no_dup <- window_db_no%>%distinct(SequenceID,.keep_all = T)
   print(paste0("case number: ", length(unique(window_db_no_dup$case))))
   print(paste0("control number: ", length(unique(window_db_no_dup$st_subject))))
   #to make sure the
   window_db_no_dup_1 <- window_db_no_dup
   window_db_no_unique_cb <- data.frame()
   window_db_no_sum_1 <- window_db_no_dup[1,]
   window_db_no_sum_2 <- window_db_no_dup[1,]
   while (nrow(window_db_no_sum_1)!=0 & nrow(window_db_no_sum_2)!=0) {
      window_db_no_sum <- window_db_no_dup_1%>%group_by(st_subject)%>%dplyr::summarise(n=n())
      window_db_no_sum_1 <- window_db_no_sum%>%filter(n==1)
      window_db_no_sum_2 <- window_db_no_sum%>%filter(n>1)
      window_db_no_unique_1 <- window_db_no_dup_1%>%filter(st_subject%in%window_db_no_sum_1$st_subject)
      window_db_no_unique_cb <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_1)
      window_db_no_dup_1 <- window_db_no_dup_1%>%filter(st_subject%in%window_db_no_sum_2$st_subject)%>%filter(!case%in%window_db_no_unique_1$case)
   }
   print(paste0("window_db_no_sum_2: ",nrow(window_db_no_sum_2)))
   if (nrow(window_db_no_sum_2)==0) {
      print("window_db_no_sum_2: 0")
      window_db_no_unique_left <- window_db_no_dup%>%filter(!st_subject%in%window_db_no_unique_cb$st_subject)%>%distinct(st_subject,.keep_all = T)
      window_db_no_unique <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_left)
   }else{
      window_db_no_dup_2 <- window_db_no_dup_1
      window_db_no_unique_cb_2 <- data.frame()
      window_db_no_sum_sum_1 <- window_db_no_dup_1[1,]
      window_db_no_sum_sum_2 <- window_db_no_dup_1[2,]
      while(nrow(window_db_no_sum_sum_1)!=0 & nrow(window_db_no_sum_sum_2)!=0){
         window_db_no_sum_sum <- window_db_no_dup_2%>%group_by(case)%>%dplyr::summarise(n=n())
         window_db_no_sum_sum_1 <- window_db_no_sum_sum%>%filter(n==1)
         window_db_no_sum_sum_2 <- window_db_no_sum_sum%>%filter(n>1)
         window_db_no_unique_2 <- window_db_no_dup_2%>%filter(case%in%window_db_no_sum_sum_1$case)
         print(paste0("some cases shared same st_subjects: ",length(unique(window_db_no_unique_2$st_subject))!=nrow(window_db_no_unique_2)))
         window_db_no_unique_2 <- window_db_no_unique_2%>%distinct(st_subject,.keep_all = T)#some cases shared the same st_subjects: ThanertR_2024_415-01 and ThanertR_2024_2202-01 from window 2
         window_db_no_unique_cb_2 <- rbind.fill(window_db_no_unique_cb_2,window_db_no_unique_2)
         window_db_no_dup_2 <- window_db_no_dup_2%>%filter(case%in%window_db_no_sum_sum_2$case)%>%filter(!st_subject%in%window_db_no_unique_2$st_subject)
      }
   }
   print(paste0("window_db_no_sum_sum_2: ",nrow(window_db_no_sum_2)))
   if (nrow(window_db_no_sum_sum_2)==0) {
      window_db_no_unique_cb <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_cb_2)
      window_db_no_unique_left <- window_db_no_dup%>%filter(!st_subject%in%window_db_no_unique_cb$st_subject)%>%distinct(st_subject,.keep_all = T)
      window_db_no_unique <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_left)
   }else{
      window_db_no_unique_cb_3 <- data.frame()
      window_db_no_dup_2_sum <- window_db_no_dup_2%>%group_by(st_subject)%>%dplyr::summarise(n=n())%>%arrange(n)
      a <- window_db_no_dup_2%>%filter(st_subject==window_db_no_dup_2_sum$st_subject[1])
      window_db_no_unique_cb_3 <- rbind.fill(window_db_no_unique_cb_3,a[1,])
      for (sub in window_db_no_dup_2_sum$st_subject[-1]) {
         b <- window_db_no_dup_2%>%filter(st_subject==sub)%>%filter(!case%in%window_db_no_unique_cb_3$case)
         window_db_no_unique_cb_3 <- rbind.fill(window_db_no_unique_cb_3,b[1,])
      }
      window_db_no_unique_cb_3 <- window_db_no_unique_cb_3%>%filter(!is.na(case))
   }
   window_db_no_unique_cb <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_cb_2,window_db_no_unique_cb_3)
   print(paste0("case number: ", length(unique(window_db_no_unique_cb$case))))
   window_db_no_unique_left <- window_db_no_dup%>%filter(!st_subject%in%window_db_no_unique_cb$st_subject)%>%distinct(st_subject,.keep_all = T)
   window_db_no_unique <- rbind.fill(window_db_no_unique_cb,window_db_no_unique_left)
   print(paste0("control number: ", length(unique(window_db_no_unique$st_subject))))
   # window_db_no_unique2 <- window_db_no%>%distinct(st_subject,.keep_all = T)#some samples are selected by >1 times
   pre_yes_2 <- rbind.fill(pre_yes_1,window_db_no_unique)
   window_db_2 <- rbind.fill(window_db_2,pre_yes_2)
}
window_db <- rbind(window_db_1,window_db_2)
window_db_yes <- window_db%>%dplyr::select(-case)%>%filter(!is.na(preday))
length(unique(window_db_yes$st_subject))#94
which(!closest_preNEC$infant%in%window_db_yes$SubjectID)
window_db_yes_sum <- window_db_yes%>%group_by(window)%>%dplyr::summarise(n=n())
window_db_no <- window_db%>%filter(is.na(preday))
length(unique(window_db_no$st_subject))#195

preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")
window_db <- left_join(window_db,preterm_metadata%>%dplyr::select(SequenceID,raw_reads,clean_reads),by="SequenceID")

#subtype
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_resistome_prenec <- preterm_resistome%>%select(names,group,window_db$SequenceID)

preterm_resistome_prenec_subtype <- preterm_resistome_prenec%>%filter(group=="subtype")%>%dplyr::select(-group)%>%column_to_rownames("names")
preterm_resistome_prenec_subtype_t <- preterm_resistome_prenec_subtype[which(rowSums(preterm_resistome_prenec_subtype)!=0),]%>%
   t()%>%as.data.frame()
names(preterm_resistome_prenec_subtype_t) <- paste0("subtype_",names(preterm_resistome_prenec_subtype_t))
subtype_beta <- vegdist(preterm_resistome_prenec_subtype_t, method = "bray")
subtype_beta_2 <- subtype_beta%>%as.matrix()%>%as.data.frame()
subtype_beta_2[upper.tri(subtype_beta_2, diag=T)] <- NA
subtype_beta_2_g <- subtype_beta_2%>%rownames_to_column("run1")%>%gather(.,run2,distance,-run1)%>%filter(!is.na(distance))
subtype_beta_2_g <- left_join(subtype_beta_2_g,metadata_prenec%>%dplyr::select(SequenceID,st_subject),by=c("run1"="SequenceID"))
subtype_beta_2_g <- left_join(subtype_beta_2_g,metadata_prenec%>%dplyr::select(SequenceID,st_subject),by=c("run2"="SequenceID"))

nec_beta_yes_no <- data.frame()
nec_beta_yes_no_pval <- data.frame("preday"=0,"median_yes"=0,"mean_yes"=0,"median_no"=0,"mean_no"=0,"pval"=0,stringsAsFactors = F)[-1,]
for (i in c(0,seq(1,28,3)[-1],31,41)) {
   # i=0
   print(i)
   window_db_yes_a <- window_db%>%filter(nec=="yes")%>%filter(window==i)
   dist_yes <- subtype_beta_2_g%>%filter(run1%in%window_db_yes_a$SequenceID)%>%filter(run2%in%window_db_yes_a$SequenceID)%>%
      mutate(nec="yes")%>%mutate(window=i)
   print(which(dist_yes$st_subject.x==dist_yes$st_subject.y))
   window_db_no_a <- window_db%>%filter(nec=="no")%>%filter(window==i)
   dist_no <- subtype_beta_2_g%>%filter(run1%in%window_db_no_a$SequenceID)%>%filter(run2%in%window_db_no_a$SequenceID)%>%
      mutate(nec="no")%>%mutate(window=i)
   print(which(dist_no$st_subject.x==dist_no$st_subject.y))
   dist_yes_no <- rbind(dist_no,dist_yes)
   ptest <- wilcox.test(dist_yes_no$distance ~ dist_yes_no$nec)
   pval <- ptest$p.value
   nec_beta_yes_no_pval[nrow(nec_beta_yes_no_pval)+1,] <- c(i,median(dist_yes$distance),mean(dist_yes$distance),
                                                            median(dist_no$distance),mean(dist_no$distance),
                                                            pval)
   nec_beta_yes_no <- rbind(nec_beta_yes_no,dist_yes_no)
}
nec_beta_yes_no_pval$fdr <- p.adjust(nec_beta_yes_no_pval$pval,method = "fdr")
nec_beta_yes_no_pval$fc_mean <- nec_beta_yes_no_pval$mean_no/nec_beta_yes_no_pval$mean_yes
nec_beta_yes_no_pval$fc_median <- nec_beta_yes_no_pval$median_no/nec_beta_yes_no_pval$median_yes
nec_beta_yes_no_sum <- nec_beta_yes_no%>%group_by(window,nec)%>%dplyr::summarise(median=median(distance))
nec_beta_yes_no_yes <- nec_beta_yes_no%>%filter(nec=="yes")
set.seed(10)
nec_beta_yes_no_no <- nec_beta_yes_no%>%filter(nec=="no")%>%dplyr::sample_n(nrow(nec_beta_yes_no_yes))
nec_beta_yes_no_2 <- rbind(nec_beta_yes_no_yes,nec_beta_yes_no_no)
nec_beta_yes_no_2$nec <- factor(nec_beta_yes_no_2$nec,levels = c("yes","no"))
nec_beta_yes_no_2 <- nec_beta_yes_no_2%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Figure 5a
nec_preday <- ggplot(data=nec_beta_yes_no_2, aes(x=window_2, y=distance, color=nec)) +
   geom_jitter(alpha=0.2,size=1,width = 0.2)+
   geom_smooth(method=loess,formula = 'y ~ x', level=0.95,alpha = 0.5,aes(fill = nec))+
   scale_color_manual(values = c("#E04B37","#4994C4"),labels=c("yes"="NEC","no"="No NEC"))+
   labs(x ="Days before NEC onset", y="Bray-Curtis dissimilarity\n(based on ARG subtypes)")+
   theme_bw()+
   scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.25,0.5,0.75,1))+
   scale_x_reverse(breaks = unique(nec_beta_yes_no_2$window_2),
                   labels=c("3-0","6-4","9-7","12-10","15-13","18-16","21-19","24-22","27-25","30-28","40-31","60-41"))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(colour = "black"))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90, vjust = 0.5,hjust=1),
         axis.title.x = element_text(size=12, colour = "black"))

#differential subtypes from window 0, 4, 7
window_db_sel <- window_db%>%filter(window%in%c(0,4,7))
table(window_db_sel$Study,window_db_sel$nec)
window_db_sel_sum <- window_db_sel%>%group_by(SequenceID)%>%dplyr::summarise(n=length(unique(window)))
window_db_sel_unique <- window_db_sel%>%distinct(SequenceID,.keep_all = T)
window_db_sel_subtype <- preterm_resistome_prenec_subtype_t%>%rownames_to_column("SequenceID")%>%filter(SequenceID%in%window_db_sel$SequenceID)
window_db_sel_subtype <- window_db_sel_subtype%>%dplyr::select(SequenceID,which(colSums(.[-1])!=0)+1)
window_db_subtype <- preterm_resistome_prenec_subtype_t%>%rownames_to_column("SequenceID")%>%filter(SequenceID%in%window_db$SequenceID)
window_db_subtype <- window_db_subtype%>%dplyr::select(SequenceID,which(colSums(.[-1])!=0)+1)
std.error <- function(x) sd(x)/sqrt(length(x))
window_db_sel_subtype_pval <- data.frame("subtype"="x","mean_nec"=0,"median_nec"=0,"q25_nec"=0,"q75_nec"=0,"se_nec"=0,"freq_nec"=0,
                                         "mean_no"=0,"median_no"=0,"q25_no"=0,"q75_no"=0,"se_no"=0,"freq_no"=0,"pval_wil"=0,"pval_lm"=0,"mean"=0,"max"=0,"prevalence"=0,
                                         stringsAsFactors = F)[-1,]
count=0
for (s in names(window_db_sel_subtype)[-1]) {
   # s="subtype_beta_lactam__ACT-17"
   count=count+1
   print(count)
   db_1 <- window_db_sel_subtype%>%dplyr::select(SequenceID,s)%>%
      left_join(.,metadata_prenec%>%dplyr::select(SequenceID,Study,st_subject,nec),by="SequenceID")
   names(db_1)[2] <- "subtype"
   db_1_sum <- db_1%>%group_by(nec)%>%dplyr::summarise(mean=mean(subtype),median=quantile(subtype)[3],q25=quantile(subtype)[2],
                                                       q75=quantile(subtype)[4],se=std.error(subtype),freq=sum(subtype>0))
   db_1$Study <- as.factor(db_1$Study)
   db_1$nec <- as.factor(db_1$nec)
   pval <- pvalue(wilcox_test(subtype ~ nec | Study, data = db_1))
   # db_1_sum <- db_1%>%group_by(st_subject)%>%dplyr::summarise(n=n())
   lm_fit <- try(lmerTest::lmer(subtype ~ nec+(1|Study)+(1|st_subject),data = db_1))
   if (class(lm_fit) == "try-error") {
      window_db_sel_subtype_pval[nrow(window_db_sel_subtype_pval)+1,] <- c(s,
                                                                           db_1_sum$mean[db_1_sum$nec=="yes"],db_1_sum$median[db_1_sum$nec=="yes"],db_1_sum$q25[db_1_sum$nec=="yes"],db_1_sum$q75[db_1_sum$nec=="yes"],db_1_sum$se[db_1_sum$nec=="yes"],db_1_sum$freq[db_1_sum$nec=="yes"],
                                                                           db_1_sum$mean[db_1_sum$nec=="no"],db_1_sum$median[db_1_sum$nec=="no"],db_1_sum$q25[db_1_sum$nec=="no"],db_1_sum$q75[db_1_sum$nec=="no"],db_1_sum$se[db_1_sum$nec=="no"],db_1_sum$freq[db_1_sum$nec=="no"],
                                                                           pval,NA,
                                                                           mean(db_1$subtype),max(db_1$subtype),sum(db_1$subtype>0)/nrow(db_1))
   }else{
      pval_lm <- anova(lm_fit)
      window_db_sel_subtype_pval[nrow(window_db_sel_subtype_pval)+1,] <- c(s,
                                                                           db_1_sum$mean[db_1_sum$nec=="yes"],db_1_sum$median[db_1_sum$nec=="yes"],db_1_sum$q25[db_1_sum$nec=="yes"],db_1_sum$q75[db_1_sum$nec=="yes"],db_1_sum$se[db_1_sum$nec=="yes"],db_1_sum$freq[db_1_sum$nec=="yes"],
                                                                           db_1_sum$mean[db_1_sum$nec=="no"],db_1_sum$median[db_1_sum$nec=="no"],db_1_sum$q25[db_1_sum$nec=="no"],db_1_sum$q75[db_1_sum$nec=="no"],db_1_sum$se[db_1_sum$nec=="no"],db_1_sum$freq[db_1_sum$nec=="no"],
                                                                           pval,pval_lm$`Pr(>F)`,
                                                                           mean(db_1$subtype),max(db_1$subtype),sum(db_1$subtype>0)/nrow(db_1))
   }
}
window_db_sel_subtype_pval_1 <- window_db_sel_subtype_pval%>%mutate_at(vars(-subtype),as.numeric)
window_db_sel_subtype_pval_1$subtype <- gsub("subtype_","",window_db_sel_subtype_pval_1$subtype,fixed = T)
window_db_sel_subtype_pval_1$subtype <- gsub("other_peptide_antibiotics","OPA",window_db_sel_subtype_pval_1$subtype,fixed = T)
window_db_sel_subtype_pval_1$subtype <- gsub("macrolide-lincosamide-streptogramin","MLS",window_db_sel_subtype_pval_1$subtype,fixed = T)
window_db_sel_subtype_pval_1$fdr <- p.adjust(window_db_sel_subtype_pval_1$pval_lm,method = "fdr")
window_db_sel_subtype_pval_1$mean_diff <- window_db_sel_subtype_pval_1$mean_nec-window_db_sel_subtype_pval_1$mean_no
window_db_sel_subtype_pval_1$median_diff <- window_db_sel_subtype_pval_1$median_nec-window_db_sel_subtype_pval_1$median_no
window_db_sel_subtype_pval_1$rate_nec <- window_db_sel_subtype_pval_1$freq_nec/199*100
window_db_sel_subtype_pval_1$rate_no <- window_db_sel_subtype_pval_1$freq_no/367*100
window_db_sel_subtype_pval_1$rate_diff <- window_db_sel_subtype_pval_1$rate_nec-window_db_sel_subtype_pval_1$rate_no
window_db_sel_subtype_pval_1_sig_0.05 <- window_db_sel_subtype_pval_1%>%filter(fdr<0.05)
window_db_sel_subtype_pval_1_sig <- window_db_sel_subtype_pval_1%>%filter(fdr<0.05)
window_db_sel_subtype_pval_1_sig$type <- sapply(strsplit(window_db_sel_subtype_pval_1_sig$subtype, split='__', fixed=TRUE), function(x)(x[1]))
subtype_nec_1_sig_1 <- window_db_sel_subtype_pval_1_sig%>%filter(mean_diff>0)%>%arrange(mean_nec)
subtype_nec_1_sig_2 <- window_db_sel_subtype_pval_1_sig%>%filter(mean_diff<0)%>%arrange(desc(mean_no))
subtype_nec_1_sig_3 <- window_db_sel_subtype_pval_1_sig%>%filter(mean_diff==0)%>%arrange(desc(mean_no))
subtype_nec_1_sig_sel <- rbind(subtype_nec_1_sig_2,subtype_nec_1_sig_3,subtype_nec_1_sig_1)
subtype_nec_1_sig_sel$subtype <- gsub("Klebsiella pneumoniae","KP",subtype_nec_1_sig_sel$subtype,fixed = T)
subtype_nec_1_sig_yes <- subtype_nec_1_sig_sel%>%dplyr::select(subtype,mean_nec,se_nec)%>%mutate("group"="NEC")
names(subtype_nec_1_sig_yes)[2] <- "mean"
names(subtype_nec_1_sig_yes)[3] <- "se"
subtype_nec_1_sig_no <- subtype_nec_1_sig_sel%>%dplyr::select(subtype,mean_no,se_no)%>%mutate("group"="No NEC")
names(subtype_nec_1_sig_no)[2] <- "mean"
names(subtype_nec_1_sig_no)[3] <- "se"
subtype_nec_1_sig_p <- rbind(subtype_nec_1_sig_yes,subtype_nec_1_sig_no)%>%arrange(desc(mean))
subtype_nec_1_sig_p$subtype <- factor(subtype_nec_1_sig_p$subtype,levels = rev(subtype_nec_1_sig_sel$subtype))
subtype_nec_1_sig_p$group <- factor(subtype_nec_1_sig_p$group,levels = c("NEC","No NEC"))
#Figure 5c
subtype_window_nec <- ggplot(subtype_nec_1_sig_p, aes(x=subtype, y=mean, fill=group)) + 
   # geom_boxplot(outlier.size = 0)+
   geom_bar(stat="identity", width = 0.7, color="black", position=position_dodge()) +
   geom_errorbar(aes(ymin=mean, ymax=mean+se), width=0.5,
                 position=position_dodge())+
   scale_y_continuous(expand = c(0,0),limits = c(0,1))+
   scale_fill_manual(values = c("#E04B37","#4994C4"),labels=c("no"="No NEC","yes"="NEC"),guide = guide_legend(reverse = TRUE))+
   labs(x ="Mean abundance of ARG subtypes", y="")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour = "black",angle = 90,hjust = 1,vjust=0.5),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_blank(),
         plot.margin = unit(c(5,1,1,1),"mm"))
###########################################################################################################################################
#species from metaphlan4
###########################################################################################################################################
metadata_prenec <- read.csv("data/metadata_prenec_closest_new.csv")
window_db <- read.csv("data/window_db_new.csv")
window_db_047 <- window_db%>%filter(window%in%c(0,4,7))
metaphlan4_preterm <- read.csv("data/metaphlan4_preterm_new.csv")
metaphlan4_preterm[is.na(metaphlan4_preterm)] <- 0
mpl_prenec <- metaphlan4_preterm%>%dplyr::select(clade_name,metadata_prenec$SequenceID)
#species with kingdom: bacteria, not archaea, Eukaryotes
species_level <- mpl_prenec%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
species_level$species <- sapply(strsplit(species_level$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
species_level_1 <- species_level%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
species_level_1 <- species_level_1[which(rowSums(species_level_1)!=0),]
names(species_level_1) <- paste0("species_",names(species_level_1))
species_beta <- vegdist(species_level_1, method = "bray")
species_beta_2 <- species_beta%>%as.matrix()%>%as.data.frame()
species_beta_2[upper.tri(species_beta_2, diag=T)] <- NA
species_beta_2_g <- species_beta_2%>%rownames_to_column("run1")%>%gather(.,run2,distance,-run1)%>%filter(!is.na(distance))

species_nec_beta_yes_no <- data.frame()
species_nec_beta_yes_no_pval <- data.frame("preday"=0,"median_yes"=0,"mean_yes"=0,"median_no"=0,"mean_no"=0,"pval"=0,stringsAsFactors = F)[-1,]
for (i in c(0,seq(1,28,3)[-1],31,41)) {
   # i=1
   print(i)
   window_db_yes_a <- window_db%>%filter(nec=="yes")%>%filter(window==i)
   dist_yes <- species_beta_2_g%>%filter(run1%in%window_db_yes_a$SequenceID)%>%filter(run2%in%window_db_yes_a$SequenceID)%>%
      mutate(nec="yes")%>%mutate(window=i)
   window_db_no_a <- window_db%>%filter(nec=="no")%>%filter(window==i)
   dist_no <- species_beta_2_g%>%filter(run1%in%window_db_no_a$SequenceID)%>%filter(run2%in%window_db_no_a$SequenceID)%>%
      mutate(nec="no")%>%mutate(window=i)
   dist_yes_no <- rbind(dist_no,dist_yes)
   ptest <- wilcox.test(dist_yes_no$distance ~ dist_yes_no$nec)
   pval <- ptest$p.value
   species_nec_beta_yes_no_pval[nrow(species_nec_beta_yes_no_pval)+1,] <- c(i,median(dist_yes$distance),mean(dist_yes$distance),
                                                                            median(dist_no$distance),mean(dist_no$distance),
                                                                            pval)
   species_nec_beta_yes_no <- rbind(species_nec_beta_yes_no,dist_yes_no)
}
species_nec_beta_yes_no_pval$fdr <- p.adjust(species_nec_beta_yes_no_pval$pval,method = "fdr")
species_nec_beta_yes_no_sum <- species_nec_beta_yes_no%>%group_by(window,nec)%>%dplyr::summarise(median=median(distance))
species_nec_beta_yes_no_yes <- species_nec_beta_yes_no%>%filter(nec=="yes")
species_nec_beta_yes_no_no <- species_nec_beta_yes_no%>%filter(nec=="no")
set.seed(10)
species_nec_beta_yes_no_no <- species_nec_beta_yes_no_no%>%dplyr::sample_frac(nrow(species_nec_beta_yes_no_yes)/nrow(species_nec_beta_yes_no_no))
species_nec_beta_yes_no_2 <- rbind(species_nec_beta_yes_no_yes,species_nec_beta_yes_no_no)
species_nec_beta_yes_no_2$nec <- factor(species_nec_beta_yes_no_2$nec,levels = c("yes","no"))
species_nec_beta_yes_no_2 <- species_nec_beta_yes_no_2%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Figure 5b
species_nec_preday <- ggplot(data=species_nec_beta_yes_no_2, aes(x=window_2, y=distance, color=nec)) +
   geom_jitter(alpha=0.2,size=1,width = 0.2)+
   geom_smooth(method=loess,formula = 'y ~ x', level=0.95,alpha = 0.5,aes(fill = nec))+
   scale_color_manual(values = c("#E04B37","#4994C4"),labels=c("yes"="NEC","no"="No NEC"))+
   labs(x ="Days before NEC onset", y="Bray-Curtis dissimilarity\n(based on species)")+
   theme_bw()+
   scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.25,0.5,0.75,1))+
   scale_x_reverse(breaks = unique(species_nec_beta_yes_no_2$window_2),
                   labels=c("3-0","6-4","9-7","12-10","15-13","18-16","21-19","24-22","27-25","30-28","40-31","60-41"))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(colour = "black"))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90, vjust = 0.5,hjust=1),
         axis.title.x = element_text(size=12, colour = "black"))


