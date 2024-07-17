library(dplyr)
library(coin)
library(tidyr)
library(stringr)
library(tibble)
library(vegan)
library(ggplot2)

metadata_nec_sel <- read.csv("data/metadata_nec_sel.csv")
metadata_nec_sel <- metadata_nec_sel%>%
   mutate(Feed_2=case_when(Feed=="breast"~"Exclusive",
                           Feed=="combined"~"Mix",
                           Feed=="nofeed"~"fasting",
                           Feed=="formula" | Feed =="solid"~"Non-breastmilk"))%>%
   mutate(nec=case_when(PreNEC!="noNEC"~"yes",
                        PreNEC=="noNEC"~"no"))%>%
   mutate(st_subject=paste0(Study,"_",SubjectID))
metadata_nec_sel_yes <- metadata_nec_sel%>%filter(nec=="yes")
length(unique(metadata_nec_sel_yes$st_subject))#45
metadata_nec_sel_no <- metadata_nec_sel%>%filter(nec=="no")
length(unique(metadata_nec_sel_no$st_subject))#179
metadata_nec_sel_unique <- metadata_nec_sel%>%dplyr::select(Study,SubjectID,st_subject,Delivery,Gender,Gestational_age,Birth_weight,Feed_2,Infant_AB,nec)%>%
   distinct(st_subject,.keep_all = T)
metadata_nec_sel_unique$Study <- as.factor(metadata_nec_sel_unique$Study)
metadata_nec_sel_unique$nec <- as.factor(metadata_nec_sel_unique$nec)
wilcox_test(Gestational_age ~ nec|Study,data=metadata_nec_sel_unique)
metadata_nec_sel_unique$Delivery <- as.factor(metadata_nec_sel_unique$Delivery)
cmh_test(Delivery ~ nec|Study,data=metadata_nec_sel_unique)
metadata_nec_sel_unique$Feed_2 <- as.factor(metadata_nec_sel_unique$Feed_2)
cmh_test(Feed_2 ~ nec|Study,data=metadata_nec_sel_unique)
metadata_nec_sel_unique$Infant_AB <- as.factor(metadata_nec_sel_unique$Infant_AB)
cmh_test(Infant_AB ~ nec|Study,data=metadata_nec_sel_unique)
########################################################################################################################################
#ARGs comparisons between cases and controls
########################################################################################################################################
#rgi
rgi_preterm_profile <- read.csv("~/Desktop/work/preterm/analysis/revision/rgi_preterm_resistome_coverm_bins_gtdb_stats_plasmid.csv")
rgi_preterm_profile_coverm <- rgi_preterm_profile%>%filter(sample%in%metadata_nec_sel$SequenceID)
length(unique(rgi_preterm_profile_coverm$sample))#1792
#organize the function
split_into_multiple <- function(column, pattern, into_prefix){
   cols <- str_split_fixed(column, pattern, n = Inf)
   cols[which(cols == "")] <- NA
   cols <- tibble::as.tibble(cols)
   m <- dim(cols)[2]
   names(cols) <- paste(into_prefix, 1:m, sep = "_")
   return(cols)
}
rgi_preterm_profile_coverm_gene <- rgi_preterm_profile_coverm %>%
   dplyr::select(cb,sample,rpkm,AMR.Gene.Family)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(cb,sample,rpkm,AMR.Gene.Family))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)
rgi_preterm_profile_coverm_gene_sum <- rgi_preterm_profile_coverm_gene%>%
   group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%spread(.,sample,rpkm_total)
rgi_preterm_profile_coverm_gene_sum[is.na(rgi_preterm_profile_coverm_gene_sum)] <- 0
rgi_preterm_profile_coverm_gene_sum_t <- rgi_preterm_profile_coverm_gene_sum%>%column_to_rownames("gene")%>%t()%>%as.data.frame()
rgi_preterm_profile_coverm_gene_sum2 <- rgi_preterm_profile_coverm_gene%>%
   group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))
rgi_richness <- rgi_preterm_profile_coverm_gene_sum2%>%group_by(sample)%>%dplyr::summarise(n=n())
names(rgi_richness) <- c("SequenceID", "Richness")
rgi_shannon <- diversity(rgi_preterm_profile_coverm_gene_sum_t, index = "shannon", MARGIN=1) 
rgi_shannon_df <- rgi_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(rgi_shannon_df) <- c("SequenceID", "Shannon")
rgi_alpha <- left_join(rgi_richness,rgi_shannon_df,by="SequenceID")
rgi_alpha_meta <- left_join(metadata_nec_sel,rgi_alpha,by="SequenceID")
rgi_alpha_meta <- rgi_alpha_meta%>%filter(!is.na(Richness))
rgi_alpha_meta$nec <- as.factor(rgi_alpha_meta$nec)
rgi_alpha_meta$Study <- as.factor(rgi_alpha_meta$Study)
wilcox_test(Richness ~ nec | Study, data = rgi_alpha_meta)
wilcox_test(Shannon ~ nec | Study, data = rgi_alpha_meta)
rgi_alpha_meta_g <- rgi_alpha_meta%>%dplyr::select(Richness,Shannon,SequenceID,nec)%>%gather(.,group,value,-SequenceID,-nec)
rgi_alpha_meta_g$group[rgi_alpha_meta_g$group=="Shannon"] <- "Shannon index"
#Supplementary Figure 11a
rgi_gene_alpha_p <- ggplot(data=rgi_alpha_meta_g,aes(x=nec, y=value, col=nec)) +
   geom_boxplot() + 
   facet_wrap(~group,scales = "free")+
   scale_color_manual(values = c("#4994C4","#E04B37"),labels=c("no"="Control","yes"="NEC"))+
   scale_x_discrete(labels=c("no"="Control","yes"="NEC"))+
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
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))
#beta diversity
rgi_t <- rgi_preterm_profile_coverm_gene_sum_t[which(rowSums(rgi_preterm_profile_coverm_gene_sum_t)!=0),which(colSums(rgi_preterm_profile_coverm_gene_sum_t)!=0)]
rgi_dist_beta <- vegdist(rgi_t, method = "bray")
rgi_dist_beta_meta <- metadata_nec_sel%>%filter(SequenceID%in%rownames(rgi_t))%>%
   arrange(match(SequenceID,rownames(rgi_t)))
identical(rownames(rgi_t),rgi_dist_beta_meta$SequenceID)#TRUE
set.seed(1000)
with(rgi_dist_beta_meta, vegan::adonis2(rgi_dist_beta ~ nec,data=rgi_dist_beta_meta,permutations = 1000, strata = Study))
rgi_beta_pco <- cmdscale(rgi_dist_beta, k=2, eig = T)
rgi_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                       round((rgi_beta_pco$eig[1]/sum(rgi_beta_pco$eig))*100,1),
                                       '%]', sep='')
rgi_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                       round((rgi_beta_pco$eig[2]/sum(rgi_beta_pco$eig))*100,1),
                                       '%]', sep='')
rgi_beta_pco.point <- as.data.frame(rgi_beta_pco$points)%>%rownames_to_column('SequenceID')
rgi_beta_pco.tabble <- right_join(rgi_dist_beta_meta, rgi_beta_pco.point, by="SequenceID")
rgi_beta_df.plot <- tibble(Axis1 = rgi_beta_pco.tabble$V1,
                               Axis2 = rgi_beta_pco.tabble$V2,
                               Sample_ID = rgi_beta_pco.tabble$SubjectID,
                               nec=rgi_beta_pco.tabble$nec)
#Supplementary Figure 11a
rgi_gene_beta_time_p <- ggplot(data=rgi_beta_df.plot,aes(x=Axis1, y=Axis2, col=nec)) +
   geom_point(size=1.5, alpha=1) + 
   stat_ellipse(aes(fill=nec),geom = "polygon",level = 0.95, alpha=0.05)+
   scale_color_manual(values = c("#4994C4","#E04B37"),labels=c("no"="Xontrol","yes"="NEC"))+
   theme_bw()+ 
   xlab(rgi_beta_pco_axis.1.title) + ylab(rgi_beta_pco_axis.2.title) +
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
subtype_alpha_meta$nec <- as.factor(subtype_alpha_meta$nec)
subtype_alpha_meta$Study <- as.factor(subtype_alpha_meta$Study)
wilcox_test(Richness ~ nec | Study, data = subtype_alpha_meta)
wilcox_test(Shannon ~ nec | Study, data = subtype_alpha_meta)
subtype_alpha_meta_g <- subtype_alpha_meta%>%dplyr::select(Richness,Shannon,SequenceID,nec)%>%gather(.,group,value,-SequenceID,-nec)
subtype_alpha_meta_g$group[subtype_alpha_meta_g$group=="Shannon"] <- "Shannon index"
#Supplementary Figure 11b
subtype_alpha_p <- ggplot(data=subtype_alpha_meta_g,aes(x=nec, y=value, col=nec)) +
   geom_boxplot() + 
   facet_wrap(~group,scales = "free")+
   scale_color_manual(values = c("#4994C4","#E04B37"),labels=c("no"="Control","yes"="NEC"))+
   scale_x_discrete(labels=c("no"="Control","yes"="NEC"))+
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
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))

#beta diversity
preterm_subtype_beta_dist <- vegdist(preterm_subtype_t, method = "bray")
subtype_alpha_meta_order <- subtype_alpha_meta%>%arrange(match(SequenceID,rownames(preterm_subtype_t)))
identical(rownames(preterm_subtype_t),subtype_alpha_meta_order$SequenceID)#TRUE
set.seed(1000)
with(subtype_alpha_meta_order, vegan::adonis2(preterm_subtype_beta_dist ~ nec,data=subtype_alpha_meta_order,permutations = 1000, strata = Study))
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
                               nec=subtype_beta_pco.tabble$nec)
#Supplementary Figure 11b
subtype_beta_time_p <- ggplot(data=subtype_beta_df.plot,aes(x=Axis1, y=Axis2, col=nec)) +
   geom_point(size=1.5, alpha=1) + 
   stat_ellipse(aes(fill=nec),geom = "polygon",level = 0.95, alpha=0.05)+
   scale_color_manual(values = c("#4994C4","#E04B37"),labels=c("no"="Control","yes"="NEC"))+
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
mpl_preterm <- read.csv("data/metaphlan4_preterm.csv")
mpl_term <- read.csv("data/metaphlan4_term.csv")
mpl_cb <- left_join(mpl_preterm,mpl_term,by="clade_name")
mpl_cb[is.na(mpl_cb)] <- 0
mpl_nec <- mpl_cb%>%dplyr::select(clade_name,metadata_nec_sel$SequenceID)
#species with kingdom: bacteria, no archaea, Eukaryotes
species_level <- mpl_nec%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
species_level$species <- sapply(strsplit(species_level$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
species_level_1 <- species_level%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
species_level_richness <- rowSums(species_level_1!=0)%>%as.data.frame()%>%setNames("Richness")%>%rownames_to_column("run")
species_level_shannon <- vegan::diversity(species_level_1, index = "shannon", MARGIN=1)%>%as.data.frame()%>%setNames("Shannon")%>%rownames_to_column("run")
species_alpha <- left_join(species_level_richness,species_level_shannon,by="run")
species_alpha_meta <- left_join(metadata_nec_sel,species_alpha,by=c("SequenceID"="run"))
species_alpha_meta$nec <- as.factor(species_alpha_meta$nec)
species_alpha_meta$Study <- as.factor(species_alpha_meta$Study)
wilcox_test(Richness ~ nec | Study, data = species_alpha_meta)
wilcox_test(Shannon ~ nec | Study, data = species_alpha_meta)
species_alpha_meta_g <- species_alpha_meta%>%dplyr::select(Richness,Shannon,SequenceID,nec)%>%gather(.,group,value,-SequenceID,-nec)
species_alpha_meta_g$group[species_alpha_meta_g$group=="Shannon"] <- "Shannon index"
#Supplementary Figure 11c
species_alpha_p <- ggplot(data=species_alpha_meta_g,aes(x=nec, y=value, col=nec)) +
   geom_boxplot() + 
   facet_wrap(~group,scales = "free")+
   scale_color_manual(values = c("#4994C4","#E04B37"),labels=c("no"="Control","yes"="NEC"))+
   scale_x_discrete(labels=c("no"="Control","yes"="NEC"))+
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
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))

#beta diversity
species_dist_beta <- vegdist(species_level_1, method = "bray")
species_alpha_meta_order <- species_alpha_meta%>%arrange(match(SequenceID,rownames(species_level_1)))
identical(rownames(species_level_1),species_alpha_meta_order$SequenceID)#TRUE
set.seed(1000)
with(species_alpha_meta_order, vegan::adonis2(species_dist_beta ~ nec,data=species_alpha_meta_order,permutations = 1000, strata = Study))
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
                               nec=species_beta_pco.tabble$nec)
#Supplementary Figure 11c
species_beta_time_p <- ggplot(data=species_beta_df.plot,aes(x=Axis1, y=Axis2, col=nec)) +
   geom_point(size=1.5, alpha=1) + 
   stat_ellipse(aes(fill=nec),geom = "polygon",level = 0.95, alpha=0.05)+
   scale_color_manual(values = c("#4994C4","#E04B37"),labels=c("no"="Control","yes"="NEC"))+
   theme_bw()+ 
   xlab(species_beta_pco_axis.1.title) + ylab(species_beta_pco_axis.2.title) +
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))
#########################################################################################################################################
#Accelerated convergence of gut resistome preceding NEC onset
#########################################################################################################################################
metadata_prenec <- read.csv("data/metadata_prenec_closest.csv")
metadata_prenec_yes <- metadata_prenec%>%filter(nec=="yes")%>%mutate(preday=start-DOL)
metadata_prenec_no <- metadata_prenec%>%filter(nec=="no")
window_db <- data.frame()
window_list <-list(c(0,3),c(4,6),c(7,9),c(10,12),c(13,15),c(16,18),c(19,21),c(22,24),c(25,27),c(28,30),
                   c(31,40),c(41,60))
for (w in 1:length(window_list)) {
   print(window_list[[w]])
   pre_yes <- metadata_prenec_yes%>%filter(preday>=window_list[[w]][1] & preday<=window_list[[w]][2])
   pre_yes_1 <- pre_yes%>%
      mutate(window=window_list[[w]][1])
   window_db_no <- data.frame()
   for (s in unique(pre_yes$st_subject)) {
      print(s)
      a <- pre_yes%>%filter(st_subject==s)%>%top_n(-1,preday)%>%mutate(window=window_list[[w]][1])%>%distinct(DOL,.keep_all = T)
      length(unique(a$st_subject))
      b <- metadata_prenec_no%>%filter(Study==a$Study)
      b_1 <- b%>%filter(Gestational_age>=a$Gestational_age-0.5 & Gestational_age<=a$Gestational_age+0.5)%>%
         filter(DOL>=a$DOL-1 & DOL<=a$DOL+1)%>%
         group_by(SubjectID)%>%dplyr::top_n(1,DOL)%>%mutate(preday=NA)%>%mutate(window=window_list[[w]][1])
      print(nrow(b_1))
      window_db_no <- rbind(window_db_no,b_1)
   }
   window_db_no_unique <- window_db_no%>%distinct(SequenceID,.keep_all = T)#some samples are selected by >1 times
   pre_yes_2 <- rbind(pre_yes_1,window_db_no_unique)
   window_db <- rbind(window_db,pre_yes_2)
}
window_db_yes <- window_db%>%filter(!is.na(preday))
length(unique(window_db_yes$st_subject))#41
window_db_yes_sum <- window_db_yes%>%group_by(window)%>%dplyr::summarise(n=n())
window_db_no <- window_db%>%filter(is.na(preday))
length(unique(window_db_no$st_subject))#97
window_db_no_sum <- window_db_no%>%group_by(window)%>%dplyr::summarise(n=n())

preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")
window_db <- left_join(window_db,preterm_metadata%>%dplyr::select(SequenceID,raw_reads,clean_reads),by="SequenceID")

#subtype
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_resistome_prenec <- preterm_resistome%>%select(names,group,metadata_prenec$SequenceID)
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
   # i=4
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
nec_beta_yes_no_pval
nec_beta_yes_no_sum <- nec_beta_yes_no%>%group_by(window,nec)%>%dplyr::summarise(median=median(distance))
nec_beta_yes_no$nec <- factor(nec_beta_yes_no$nec,levels = c("yes","no"))
nec_beta_yes_no <- nec_beta_yes_no%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Figure 5a
nec_preday <- ggplot(data=nec_beta_yes_no, aes(x=window_2, y=distance, color=nec)) +
   geom_jitter(alpha=0.2,size=1,width = 0.2)+
   geom_smooth(method=loess,formula = 'y ~ x', level=0.95,alpha = 0.5,aes(fill = nec))+
   scale_color_manual(values = c("#E04B37","#4994C4"),labels=c("yes"="NEC","no"="No NEC"))+
   labs(x ="Days before NEC onset", y="Bray-Curtis dissimilarity\n(based on ARG subtypes)")+
   theme_bw()+
   scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.25,0.5,0.75,1))+
   scale_x_reverse(breaks = unique(nec_beta_yes_no$window_2),
                   labels=c("3-0","6-4","9-7","12-10","15-13","18-16","21-19","24-22","27-25","30-28","40-31","60-41"))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90, vjust = 0.5,hjust=1),
         axis.title.x = element_text(size=12, colour = "black"))
#differential subtypes from window 0, 4, 7
window_db_sel <- window_db%>%filter(window%in%c(0,4,7))
window_db_sel_unique <- window_db_sel%>%distinct(SequenceID,.keep_all = T)
window_db_sel_subtype <- preterm_resistome_prenec_subtype_t%>%rownames_to_column("SequenceID")%>%filter(SequenceID%in%window_db_sel$SequenceID)
window_db_sel_subtype <- window_db_sel_subtype%>%dplyr::select(SequenceID,which(colSums(.[-1])!=0)+1)
std.error <- function(x) sd(x)/sqrt(length(x))
window_db_sel_subtype_pval <- data.frame("subtype"="x","mean_nec"=0,"median_nec"=0,"q25_nec"=0,"q75_nec"=0,"se_nec"=0,"freq_nec"=0,
                                         "mean_no"=0,"median_no"=0,"q25_no"=0,"q75_no"=0,"se_no"=0,"freq_no"=0,"pval"=0,"mean"=0,"max"=0,"prevalence"=0,
                                         stringsAsFactors = F)[-1,]
count=0
for (s in names(window_db_sel_subtype)[-1]) {
   # s="subtype_beta_lactam__ACT-17"
   count=count+1
   print(count)
   db_1 <- window_db_sel_subtype%>%dplyr::select(SequenceID,s)%>%
      left_join(.,metadata_prenec%>%dplyr::select(SequenceID,Study,nec),by="SequenceID")
   names(db_1)[2] <- "subtype"
   db_1_sum <- db_1%>%group_by(nec)%>%dplyr::summarise(mean=mean(subtype),median=quantile(subtype)[3],q25=quantile(subtype)[2],
                                                       q75=quantile(subtype)[4],se=std.error(subtype),freq=sum(subtype>0))
   db_1$Study <- as.factor(db_1$Study)
   db_1$nec <- as.factor(db_1$nec)
   pval <- pvalue(wilcox_test(subtype ~ nec | Study, data = db_1))
   window_db_sel_subtype_pval[nrow(window_db_sel_subtype_pval)+1,] <- c(s,
                                                                        db_1_sum$mean[db_1_sum$nec=="yes"],db_1_sum$median[db_1_sum$nec=="yes"],db_1_sum$q25[db_1_sum$nec=="yes"],db_1_sum$q75[db_1_sum$nec=="yes"],db_1_sum$se[db_1_sum$nec=="yes"],db_1_sum$freq[db_1_sum$nec=="yes"],
                                                                        db_1_sum$mean[db_1_sum$nec=="no"],db_1_sum$median[db_1_sum$nec=="no"],db_1_sum$q25[db_1_sum$nec=="no"],db_1_sum$q75[db_1_sum$nec=="no"],db_1_sum$se[db_1_sum$nec=="no"],db_1_sum$freq[db_1_sum$nec=="no"],
                                                                        pval,
                                                                        mean(db_1$subtype),max(db_1$subtype),sum(db_1$subtype>0)/nrow(db_1))
}
window_db_sel_subtype_pval_1 <- window_db_sel_subtype_pval%>%mutate_at(vars(-subtype),as.numeric)
window_db_sel_subtype_pval_1$subtype <- gsub("subtype_","",window_db_sel_subtype_pval_1$subtype,fixed = T)
window_db_sel_subtype_pval_1$subtype <- gsub("other_peptide_antibiotics","OPA",window_db_sel_subtype_pval_1$subtype,fixed = T)
window_db_sel_subtype_pval_1$subtype <- gsub("macrolide-lincosamide-streptogramin","MLS",window_db_sel_subtype_pval_1$subtype,fixed = T)
window_db_sel_subtype_pval_1$fdr <- p.adjust(window_db_sel_subtype_pval_1$pval,method = "fdr")
window_db_sel_subtype_pval_1$mean_diff <- window_db_sel_subtype_pval_1$mean_nec-window_db_sel_subtype_pval_1$mean_no
window_db_sel_subtype_pval_1_sig_0.05 <- window_db_sel_subtype_pval_1%>%filter(fdr<0.05)
nrow(window_db_sel_subtype_pval_1_sig_0.05)#92
window_db_sel_subtype_pval_1_sig <- window_db_sel_subtype_pval_1%>%filter(fdr<0.01)
nrow(window_db_sel_subtype_pval_1_sig)#47
window_db_sel_subtype_pval_1_sig$type <- sapply(strsplit(window_db_sel_subtype_pval_1_sig$subtype, split='__', fixed=TRUE), function(x)(x[1]))
subtype_nec_1_sig_1 <- window_db_sel_subtype_pval_1_sig%>%filter(mean_diff>0)%>%arrange(mean_nec)
subtype_nec_1_sig_2 <- window_db_sel_subtype_pval_1_sig%>%filter(mean_diff<0)%>%arrange(desc(mean_no))
subtype_nec_1_sig_sel <- rbind(subtype_nec_1_sig_2,subtype_nec_1_sig_1)
subtype_nec_1_sig_yes <- subtype_nec_1_sig_sel%>%dplyr::select(subtype,mean_nec,se_nec)%>%mutate("group"="NEC")
names(subtype_nec_1_sig_yes)[2] <- "mean"
names(subtype_nec_1_sig_yes)[3] <- "se"
subtype_nec_1_sig_no <- subtype_nec_1_sig_sel%>%dplyr::select(subtype,mean_no,se_no)%>%mutate("group"="No NEC")
names(subtype_nec_1_sig_no)[2] <- "mean"
names(subtype_nec_1_sig_no)[3] <- "se"
subtype_nec_1_sig_p <- rbind(subtype_nec_1_sig_yes,subtype_nec_1_sig_no)%>%arrange(desc(mean))
subtype_nec_1_sig_p$subtype <- factor(subtype_nec_1_sig_p$subtype,levels = subtype_nec_1_sig_sel$subtype)
subtype_nec_1_sig_p$group <- factor(subtype_nec_1_sig_p$group,levels = c("No NEC","NEC"))
subtype_nec_1_sig_p_sum <- subtype_nec_1_sig_p%>%group_by(subtype)%>%dplyr::summarise(total=sum(mean))
subtype_nec_1_sig_p_sum_2 <- subtype_nec_1_sig_p_sum%>%filter(total<0.01)
#Figure 5c
subtype_window_nec <- ggplot(subtype_nec_1_sig_p, aes(x=mean, y=subtype, fill=group)) + 
   geom_bar(stat="identity", width = 0.7, color="black", position=position_dodge()) +
   geom_errorbar(aes(xmin=mean, xmax=mean+se), width=0.5,
                 position=position_dodge())+
   scale_x_continuous(expand = c(0,0),limits = c(0,0.5),position = "top")+
   scale_fill_manual(values = c("#4994C4","#E04B37"),labels=c("no"="No NEC","yes"="NEC"),guide = guide_legend(reverse = TRUE))+
   labs(x ="Mean abundance of ARG subtypes", y="")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=10, colour = "black"),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         plot.margin = unit(c(1,5,1,1),"mm"))

subtype_nec_1_sig_p_2 <- subtype_nec_1_sig_p%>%filter(subtype%in%subtype_nec_1_sig_p_sum_2$subtype)
#Figure 5c
subtype_window_nec_2 <- ggplot(subtype_nec_1_sig_p_2, aes(x=mean, y=subtype, fill=group)) + 
   geom_bar(stat="identity", width = 0.7, color="black", position=position_dodge()) +
   geom_errorbar(aes(xmin=mean, xmax=mean+se), width=0.5,
                 position=position_dodge())+
   scale_x_continuous(expand = c(0,0),limits = c(0,0.012),position = "top")+
   scale_fill_manual(values = c("#4994C4","#E04B37"),labels=c("no"="No NEC","yes"="NEC"),guide = guide_legend(reverse = TRUE))+
   labs(x ="Mean abundance of ARG subtypes", y="")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=10, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         plot.margin = unit(c(1,5,1,5),"mm"))

#rgi
window_db <- read.csv("data/window_db.csv")
window_db_047 <- window_db%>%filter(window%in%c(0,4,7))
#rgi
rgi_profile <- read.csv("data/rgi_preterm_resistome_coverm_bins_gtdb_stats_plasmid.csv")
rgi_profile_sel <- rgi_profile%>%filter(sample%in%metadata_prenec$SequenceID)
rgi_profile_sel$family <- sapply(strsplit(rgi_profile_sel$GTDB_classification, split=';', fixed=TRUE), function(x)(x[5]))
rgi_profile_sel$genus <- sapply(strsplit(rgi_profile_sel$GTDB_classification, split=';', fixed=TRUE), function(x)(x[6]))
rgi_profile_sel$species <- sapply(strsplit(rgi_profile_sel$GTDB_classification, split=';', fixed=TRUE), function(x)(x[7]))
rgi_profile_sel$family[is.na(rgi_profile_sel$family)] <- "unknown"
rgi_profile_sel$genus[is.na(rgi_profile_sel$genus)] <- "unknown"
rgi_profile_sel$species[rgi_profile_sel$species=="s__"] <- "unknown"
rgi_profile_sel$species[is.na(rgi_profile_sel$species)] <- "unknown"

rgi_profile_sel_sum <- rgi_profile_sel%>%group_by(sample,Best_Hit_ARO)%>%
   group_by(sample,Best_Hit_ARO)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%
   spread(.,sample,rpkm_total)
rgi_profile_sel_sum[is.na(rgi_profile_sel_sum)] <- 0
rgi_profile_sel_sum_t <- rgi_profile_sel_sum%>%
   dplyr::select(Best_Hit_ARO,which(colSums(.[-1])!=0)+1)%>%column_to_rownames("Best_Hit_ARO")%>%t()%>%as.data.frame()%>%
   dplyr::select(which(colSums(.)!=0))

#organize the function
split_into_multiple <- function(column, pattern, into_prefix){
   cols <- str_split_fixed(column, pattern, n = Inf)
   cols[which(cols == "")] <- NA
   cols <- tibble::as.tibble(cols)
   m <- dim(cols)[2]
   names(cols) <- paste(into_prefix, 1:m, sep = "_")
   return(cols)
}
rgi_preterm_profile_coverm_gene <- rgi_profile_sel %>%
   dplyr::select(cb,sample,rpkm,AMR.Gene.Family,family,genus,species)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(cb,sample,rpkm,AMR.Gene.Family,family,genus,species))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)
rgi_preterm_profile_coverm_gene_sum <- rgi_preterm_profile_coverm_gene%>%
   group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%
   spread(.,sample,rpkm_total)
rgi_preterm_profile_coverm_gene_sum[is.na(rgi_preterm_profile_coverm_gene_sum)] <- 0
rgi_preterm_profile_coverm_gene_sum_t <- rgi_preterm_profile_coverm_gene_sum%>%
   dplyr::select(gene,which(colSums(.[-1])!=0)+1)%>%column_to_rownames("gene")%>%t()%>%as.data.frame()%>%
   dplyr::select(which(colSums(.)!=0))

#based on the aro abundance
rgi_gene_beta <- vegdist(rgi_profile_sel_sum_t, method = "bray")
rgi_gene_beta_2 <- rgi_gene_beta%>%as.matrix()%>%as.data.frame()
rgi_gene_beta_2[upper.tri(rgi_gene_beta_2, diag=T)] <- NA
rgi_gene_beta_2_g <- rgi_gene_beta_2%>%rownames_to_column("run1")%>%gather(.,run2,distance,-run1)%>%filter(!is.na(distance))
rgi_nec_beta_yes_no <- data.frame()
rgi_nec_beta_yes_no_pval <- data.frame("preday"=0,"median_yes"=0,"mean_yes"=0,"median_no"=0,"mean_no"=0,"pval"=0,stringsAsFactors = F)[-1,]
for (i in c(0,seq(1,28,3)[-1],31,41)) {
   # i=1
   print(i)
   window_db_yes_a <- window_db%>%filter(nec=="yes")%>%filter(window==i)
   dist_yes <- rgi_gene_beta_2_g%>%filter(run1%in%window_db_yes_a$SequenceID)%>%filter(run2%in%window_db_yes_a$SequenceID)%>%
      mutate(nec="yes")%>%mutate(window=i)
   window_db_no_a <- window_db%>%filter(nec=="no")%>%filter(window==i)
   dist_no <- rgi_gene_beta_2_g%>%filter(run1%in%window_db_no_a$SequenceID)%>%filter(run2%in%window_db_no_a$SequenceID)%>%
      mutate(nec="no")%>%mutate(window=i)
   dist_yes_no <- rbind(dist_no,dist_yes)
   ptest <- wilcox.test(dist_yes_no$distance ~ dist_yes_no$nec)
   pval <- ptest$p.value
   rgi_nec_beta_yes_no_pval[nrow(rgi_nec_beta_yes_no_pval)+1,] <- c(i,median(dist_yes$distance),mean(dist_yes$distance),
                                                                    median(dist_no$distance),mean(dist_no$distance),
                                                                    pval)
   rgi_nec_beta_yes_no <- rbind(rgi_nec_beta_yes_no,dist_yes_no)
}
rgi_nec_beta_yes_no_pval$fdr <- p.adjust(rgi_nec_beta_yes_no_pval$pval,method = "fdr")
rgi_nec_beta_yes_no_pval
rgi_nec_beta_yes_no$nec <- factor(rgi_nec_beta_yes_no$nec,levels = c("yes","no"))
rgi_nec_beta_yes_no <- rgi_nec_beta_yes_no%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Figure 5b
rgi_nec_preday <- ggplot(data=rgi_nec_beta_yes_no, aes(x=window_2, y=distance, color=nec)) +
   geom_jitter(alpha=0.2,size=1,width = 0.2)+
   geom_smooth(method=loess,formula = 'y ~ x', level=0.95,alpha = 0.5,aes(fill = nec))+
   scale_color_manual(values = c("#E04B37","#4994C4"),labels=c("yes"="NEC","no"="No NEC"))+
   labs(x ="Days before NEC onset", y="Bray-Curtis dissimilarity\n(based on assembled ARGs)")+
   theme_bw()+
   scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.25,0.5,0.75,1))+
   scale_x_reverse(breaks = unique(rgi_nec_beta_yes_no$window_2),
                   labels=c("3-0","6-4","9-7","12-10","15-13","18-16","21-19","24-22","27-25","30-28","40-31","60-41"))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90, vjust = 0.5,hjust=1),
         axis.title.x = element_text(size=12, colour = "black"))

#based on the gene families abundance
rgi_gene_beta_family <- vegdist(rgi_preterm_profile_coverm_gene_sum_t, method = "bray")
rgi_gene_beta_family_2 <- rgi_gene_beta_family%>%as.matrix()%>%as.data.frame()
rgi_gene_beta_family_2[upper.tri(rgi_gene_beta_family_2, diag=T)] <- NA
rgi_gene_beta_family_2_g <- rgi_gene_beta_family_2%>%rownames_to_column("run1")%>%gather(.,run2,distance,-run1)%>%filter(!is.na(distance))
rgi_gene_beta_family_yes_no <- data.frame()
rgi_gene_beta_family_yes_no_pval <- data.frame("preday"=0,"median_yes"=0,"mean_yes"=0,"median_no"=0,"mean_no"=0,"pval"=0,stringsAsFactors = F)[-1,]
for (i in c(0,seq(1,28,3)[-1],31,41)) {
   # i=1
   print(i)
   window_db_yes_a <- window_db%>%filter(nec=="yes")%>%filter(window==i)
   dist_yes <- rgi_gene_beta_family_2_g%>%filter(run1%in%window_db_yes_a$SequenceID)%>%filter(run2%in%window_db_yes_a$SequenceID)%>%
      mutate(nec="yes")%>%mutate(window=i)
   window_db_no_a <- window_db%>%filter(nec=="no")%>%filter(window==i)
   dist_no <- rgi_gene_beta_family_2_g%>%filter(run1%in%window_db_no_a$SequenceID)%>%filter(run2%in%window_db_no_a$SequenceID)%>%
      mutate(nec="no")%>%mutate(window=i)
   dist_yes_no <- rbind(dist_no,dist_yes)
   ptest <- wilcox.test(dist_yes_no$distance ~ dist_yes_no$nec)
   pval <- ptest$p.value
   rgi_gene_beta_family_yes_no_pval[nrow(rgi_gene_beta_family_yes_no_pval)+1,] <- c(i,median(dist_yes$distance),mean(dist_yes$distance),
                                                                                    median(dist_no$distance),mean(dist_no$distance),
                                                                                    pval)
   rgi_gene_beta_family_yes_no <- rbind(rgi_gene_beta_family_yes_no,dist_yes_no)
}
rgi_gene_beta_family_yes_no_pval$fdr <- p.adjust(rgi_gene_beta_family_yes_no_pval$pval,method = "fdr")
rgi_gene_beta_family_yes_no_pval

rgi_gene_beta_family_yes_no$nec <- factor(rgi_gene_beta_family_yes_no$nec,levels = c("yes","no"))
rgi_gene_beta_family_yes_no <- rgi_gene_beta_family_yes_no%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Supplementary Figure 11d
rgi_nec_preday <- ggplot(data=rgi_gene_beta_family_yes_no, aes(x=window_2, y=distance, color=nec)) +
   geom_jitter(alpha=0.2,size=1,width = 0.2)+
   geom_smooth(method=loess,formula = 'y ~ x', level=0.95,alpha = 0.5,aes(fill = nec))+
   scale_color_manual(values = c("#E04B37","#4994C4"),labels=c("yes"="NEC","no"="No NEC"))+
   labs(x ="Days before NEC onset", y="Bray-Curtis dissimilarity\n(based on ARG families)")+
   theme_bw()+
   scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.25,0.5,0.75,1))+
   scale_x_reverse(breaks = unique(rgi_nec_beta_yes_no$window_2),
                   labels=c("3-0","6-4","9-7","12-10","15-13","18-16","21-19","24-22","27-25","30-28","40-31","60-41"))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90, vjust = 0.5,hjust=1),
         axis.title.x = element_text(size=12, colour = "black"))

#differential ARG families
window_db_sel_rgi <- rgi_preterm_profile_coverm_gene_sum_t%>%rownames_to_column("SequenceID")%>%filter(SequenceID%in%window_db_047$SequenceID)
window_db_sel_rgi <- window_db_sel_rgi%>%dplyr::select(SequenceID,which(colSums(.[-1])!=0)+1)
window_db_sel_unique_rgi <- window_db_047%>%filter(SequenceID%in%window_db_sel_rgi$SequenceID)%>%distinct(SequenceID,.keep_all = T)
std.error <- function(x) sd(x)/sqrt(length(x))
window_db_sel_rgi_pval <- data.frame("rgi"="x","mean_nec"=0,"median_nec"=0,"q25_nec"=0,"q75_nec"=0,"se_nec"=0,"freq_nec"=0,
                                     "mean_no"=0,"median_no"=0,"q25_no"=0,"q75_no"=0,"se_no"=0,"freq_no"=0,"pval"=0,
                                     stringsAsFactors = F)[-1,]
count=0
for (s in names(window_db_sel_rgi)[-1]) {
   # s="rgi_vanY"
   count=count+1
   print(count)
   db_1 <- window_db_sel_rgi%>%dplyr::select(SequenceID,s)%>%
      left_join(.,metadata_prenec%>%dplyr::select(SequenceID,Study,nec),by="SequenceID")
   names(db_1)[2] <- "rgi"
   db_1_sum <- db_1%>%group_by(nec)%>%dplyr::summarise(mean=mean(rgi),median=quantile(rgi)[3],q25=quantile(rgi)[2],
                                                       q75=quantile(rgi)[4],se=std.error(rgi),freq=sum(rgi>0))
   db_1$Study <- as.factor(db_1$Study)
   db_1$nec <- as.factor(db_1$nec)
   pval <- pvalue(wilcox_test(rgi ~ nec | Study, data = db_1))
   window_db_sel_rgi_pval[nrow(window_db_sel_rgi_pval)+1,] <- c(s,
                                                                db_1_sum$mean[db_1_sum$nec=="yes"],db_1_sum$median[db_1_sum$nec=="yes"],db_1_sum$q25[db_1_sum$nec=="yes"],db_1_sum$q75[db_1_sum$nec=="yes"],db_1_sum$se[db_1_sum$nec=="yes"],db_1_sum$freq[db_1_sum$nec=="yes"],
                                                                db_1_sum$mean[db_1_sum$nec=="no"],db_1_sum$median[db_1_sum$nec=="no"],db_1_sum$q25[db_1_sum$nec=="no"],db_1_sum$q75[db_1_sum$nec=="no"],db_1_sum$se[db_1_sum$nec=="no"],db_1_sum$freq[db_1_sum$nec=="no"],
                                                                pval)
}
window_db_sel_rgi_pval_1 <- window_db_sel_rgi_pval%>%mutate_at(vars(-rgi),as.numeric)
window_db_sel_rgi_pval_1$rgi <- gsub("rgi_","",window_db_sel_rgi_pval_1$rgi,fixed = T)
window_db_sel_rgi_pval_1$fdr <- p.adjust(window_db_sel_rgi_pval_1$pval,method = "fdr")
window_db_sel_rgi_pval_1$mean_diff <- window_db_sel_rgi_pval_1$mean_nec-window_db_sel_rgi_pval_1$mean_no
window_db_sel_rgi_pval_1$rate_nec <- window_db_sel_rgi_pval_1$freq_nec/117*100
window_db_sel_rgi_pval_1$rate_no <- window_db_sel_rgi_pval_1$freq_no/162*100
window_db_sel_rgi_pval_1$rate_diff <- window_db_sel_rgi_pval_1$rate_nec-window_db_sel_rgi_pval_1$rate_no
window_db_sel_rgi_pval_1_sig <- window_db_sel_rgi_pval_1%>%filter(fdr<0.05)
rgi_nec_1_sig_1 <- window_db_sel_rgi_pval_1_sig%>%filter(mean_diff>0)%>%arrange(mean_nec)
rgi_nec_1_sig_2 <- window_db_sel_rgi_pval_1_sig%>%filter(mean_diff<0)%>%arrange(desc(mean_no))
rgi_nec_1_sig_sel <- rbind(rgi_nec_1_sig_2,rgi_nec_1_sig_1)

rgi_nec_1_sig_yes <- rgi_nec_1_sig_sel%>%dplyr::select(rgi,mean_nec,se_nec)%>%mutate("group"="NEC")
names(rgi_nec_1_sig_yes)[2] <- "mean"
names(rgi_nec_1_sig_yes)[3] <- "se"
rgi_nec_1_sig_no <- rgi_nec_1_sig_sel%>%dplyr::select(rgi,mean_no,se_no)%>%mutate("group"="No NEC")
names(rgi_nec_1_sig_no)[2] <- "mean"
names(rgi_nec_1_sig_no)[3] <- "se"
rgi_nec_1_sig_p <- rbind(rgi_nec_1_sig_yes,rgi_nec_1_sig_no)%>%arrange(desc(mean))
rgi_nec_1_sig_p$rgi <- factor(rgi_nec_1_sig_p$rgi,levels = rgi_nec_1_sig_sel$rgi)
rgi_nec_1_sig_p$group <- factor(rgi_nec_1_sig_p$group,levels = c("No NEC","NEC"))
rgi_nec_1_sig_p_sum <- rgi_nec_1_sig_p%>%group_by(rgi)%>%dplyr::summarise(total=sum(mean))
rgi_nec_1_sig_p_sum_2 <- rgi_nec_1_sig_p_sum%>%filter(total<120)
#Figure 5d
rgi_window_nec <- ggplot(rgi_nec_1_sig_p, aes(x=mean, y=rgi, fill=group)) + 
   geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
   geom_errorbar(aes(xmin=mean, xmax=mean+se), width=0.5,
                 position=position_dodge(0.9))+
   scale_x_continuous(expand = c(0,0),limits = c(0,44000),breaks = c(0,10000,20000,30000,40000),
                      labels = c("0","10,000","20,000","30,000","40,000"))+
   scale_fill_manual(values = c("#4994C4","#E04B37"),labels=c("no"="No NEC","yes"="NEC"),guide = guide_legend(reverse = TRUE))+
   labs(x ="Mean abundance of ARG gene families", y="")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         plot.margin = unit(c(1,5,1,1),"mm"))

rgi_nec_1_sig_p_2 <- rgi_nec_1_sig_p%>%filter(rgi%in%rgi_nec_1_sig_p_sum_2$rgi)
#Figure 5d
rgi_window_nec_2 <- ggplot(rgi_nec_1_sig_p_2, aes(x=mean, y=rgi, fill=group)) + 
   geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
   geom_errorbar(aes(xmin=mean, xmax=mean+se), width=0.5,
                 position=position_dodge(0.9))+
   scale_x_continuous(expand = c(0,0),limits = c(0,170),position = "top")+
   scale_fill_manual(values = c("#4994C4","#E04B37"),labels=c("no"="No NEC","yes"="NEC"),guide = guide_legend(reverse = TRUE))+
   labs(x ="Mean abundance of ARG gene families", y="")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         plot.margin = unit(c(1,5,1,1),"mm"))

#drug for differential ARG gene familes
rgi_preterm_profile_coverm_gene_2 <- rgi_profile_sel %>%
   dplyr::select(cb,sample,rpkm,Drug.Class,AMR.Gene.Family,family,genus,species)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(cb,sample,rpkm,Drug.Class,AMR.Gene.Family,family,genus,species))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)
rgi_preterm_profile_coverm_gene_drug <- rgi_preterm_profile_coverm_gene_2 %>%
   dplyr::bind_cols(split_into_multiple(.$Drug.Class, "; ", "cat"))%>%
   gather(.,cat,drug,-c(cb,sample,rpkm,Drug.Class,AMR.Gene.Family,family,genus,species,gene))%>%
   filter(!is.na(drug))%>%
   dplyr::select(-cat)
rgi_profile_sel_2 <- rgi_preterm_profile_coverm_gene_drug%>%filter(sample%in%window_db_047$SequenceID)%>%
   filter(gene%in%window_db_sel_rgi_pval_1_sig$rgi)
rgi_profile_sel_2_1 <- rgi_profile_sel_2%>%group_by(sample,gene,drug)%>%dplyr::summarise(total=sum(rpkm))%>%
   left_join(.,metadata_prenec%>%dplyr::select(SequenceID,Study,nec),by=c("sample"="SequenceID"))
rgi_profile_sel_2_2 <- rgi_profile_sel_2_1%>%group_by(nec,gene,drug)%>%dplyr::summarise(total=sum(total))
rgi_profile_sel_2_2_yes <- rgi_profile_sel_2_2%>%filter(nec=="yes")%>%mutate(mean=total/117)
rgi_profile_sel_2_2_no <- rgi_profile_sel_2_2%>%filter(nec=="no")%>%mutate(mean=total/162)
rgi_profile_sel_2_3 <- rbind(rgi_profile_sel_2_2_yes,rgi_profile_sel_2_2_no)
rgi_profile_sel_2_3$gene <- factor(rgi_profile_sel_2_3$gene,levels = rgi_nec_1_sig_sel$rgi)

rgi_profile_sel_2_3 <- rgi_profile_sel_2_3%>%group_by(gene)%>%dplyr::mutate(gene_2=cur_group_id())
rgi_profile_sel_2_3$drug[rgi_profile_sel_2_3$drug=="disinfecting agents and antiseptics"] <- "DAA"
rgi_profile_sel_2_3$drug <- gsub("antibiotic","",rgi_profile_sel_2_3$drug,fixed = T)

barwidth=0.45
#Figure 5d
rgi_window_nec_drug <- ggplot() + 
   geom_bar(data = rgi_profile_sel_2_3%>%filter(nec=="no"), color="black",
            mapping = aes(x = gene_2-barwidth/2, y = mean, fill = as.factor(drug)), 
            stat="identity", 
            position='stack', 
            width = barwidth) + 
   geom_bar(data = rgi_profile_sel_2_3%>%filter(nec=="yes"), color="black",
            mapping = aes(x = gene_2 + barwidth/2, y = mean, fill = as.factor(drug)), 
            stat="identity", 
            position='stack' , 
            width = barwidth)+
   scale_fill_manual(values = c("#8DD3C7", "#BEBADA", "#F1F1A8", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
                                "#FCCDE5", "#D9D9D9", "#BC80BD",  "#CCEBC5", "#E8D863","#FFA833","#BE643A","#9C76EC"))+
   scale_y_continuous(expand = c(0,0),limits = c(0,450000),breaks = c(0,100000,200000,300000,400000),
                      labels = c("0","100,000","200,000","300,000","400,000"))+
   scale_x_continuous(breaks = c(1,2,3,4,5,6),labels = c("1","2","3","4","5","6"))+
   theme_bw()+
   guides(fill=guide_legend(ncol=1))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   guides(fill=guide_legend(ncol=2))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=12, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         plot.margin = unit(c(1,5,1,5),"mm"))+
   coord_flip()
#Figure 5d
rgi_profile_sel_2_3_2 <- rgi_profile_sel_2_3%>%filter(gene_2==2)
rgi_window_nec_drug_2 <- ggplot(rgi_profile_sel_2_3_2, aes(x=gene, y=mean, fill=nec)) + 
   geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
   scale_fill_manual(values = c("#BC80BD","#BC80BD"))+
   scale_y_continuous(expand = c(0,0),limits = c(0,120),position = "right")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         plot.margin = unit(c(1,5,1,5),"mm"))+
   coord_flip()

rgi_profile_sel_4 <- rgi_preterm_profile_coverm_gene_2%>%filter(sample%in%window_db_047$SequenceID)%>%filter(gene%in%window_db_sel_rgi_pval_1_sig$rgi)
rgi_profile_sel_4%>%group_by(family)%>%dplyr::summarise(n=n())%>%arrange(desc(n))
rgi_profile_sel_4%>%filter(family=="f__Enterobacteriaceae")%>%group_by(genus)%>%dplyr::summarise(n=n())%>%arrange(desc(n))
rgi_profile_sel_4%>%group_by(genus)%>%dplyr::summarise(n=n())%>%filter(genus!="unknown")%>%arrange(desc(n))%>%top_n(5)
rgi_profile_sel_4%>%group_by(species)%>%dplyr::summarise(n=n())%>%arrange(desc(n))%>%top_n(10)
nrow(rgi_profile_sel_4)-306#1009

###########################################################################################################################################
#species from metaphlan4
###########################################################################################################################################
metaphlan4_preterm <- read.csv("data/metaphlan4_preterm.csv")
metaphlan4_preterm[is.na(metaphlan4_preterm)] <- 0
mpl_prenec <- metaphlan4_preterm%>%dplyr::select(clade_name,metadata_prenec$SequenceID)
#species with kingdom: bacteria, no archaea, Eukaryotes
species_level <- mpl_prenec%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
species_level$species <- sapply(strsplit(species_level$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
species_level_1 <- species_level%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
names(species_level_1) <- paste0("species_",names(species_level_1))
species_level_richness <- rowSums(species_level_1!=0)%>%as.data.frame()%>%setNames("species_Richness")
species_level_shannon <- vegan::diversity(species_level_1, index = "shannon", MARGIN=1)%>%as.data.frame()%>%setNames("species_Shannon")

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
species_nec_beta_yes_no_pval
species_nec_beta_yes_no_sum <- species_nec_beta_yes_no%>%group_by(window,nec)%>%dplyr::summarise(median=median(distance))
species_nec_beta_yes_no$nec <- factor(species_nec_beta_yes_no$nec,levels = c("yes","no"))
species_nec_beta_yes_no <- species_nec_beta_yes_no%>%group_by(window)%>%dplyr::mutate(window_2=cur_group_id())
#Supplementary Figure 11e
species_nec_preday <- ggplot(data=species_nec_beta_yes_no, aes(x=window_2, y=distance, color=nec)) +
   geom_jitter(alpha=0.2,size=1,width = 0.2)+
   geom_smooth(method=loess,formula = 'y ~ x', level=0.95,alpha = 0.5,aes(fill = nec))+
   scale_color_manual(values = c("#E04B37","#4994C4"),labels=c("yes"="NEC","no"="No NEC"))+
   labs(x ="Days before NEC onset", y="Bray-Curtis dissimilarity\n(based on species)")+
   theme_bw()+
   scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.25,0.5,0.75,1))+
   scale_x_reverse(breaks = unique(species_nec_beta_yes_no$window_2),
                   labels=c("3-0","6-4","9-7","12-10","15-13","18-16","21-19","24-22","27-25","30-28","40-31","60-41"))+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90, vjust = 0.5,hjust=1),
         axis.title.x = element_text(size=12, colour = "black"))

