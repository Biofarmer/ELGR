library(dplyr)
library(tibble)
library(cowplot)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(plyr)
library(coin)
library(ggtrendline)
library(vegan)
library(lmerTest)

std.error <- function(x) sd(x)/sqrt(length(x))
term_metadata <- read.csv("data/metadata_cb_infant_term.csv")
term_metadata <- term_metadata%>%
   mutate(month=case_when(DOL>=0 & DOL<=7 ~ "0.25",
                          DOL>7 & DOL<=14 ~ "0.5",
                          DOL>14 & DOL<=21 ~ "0.75",
                          DOL>21 & DOL<=30 ~ "1",
                          DOL>30 & DOL<=60 ~ "2",
                          DOL>60 & DOL<=90 ~ "3",
                          DOL>90 & DOL<=180 ~ "6",
                          DOL>180 & DOL<=360 ~ "12",
                          DOL>360 ~ "18"))
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
                          DOL>360 ~ "18"))
metadata_cb <- rbind(preterm_metadata,term_metadata)
##############################################################################################################################
#ARGs profile from full-term infants
##############################################################################################################################
term_resistome <- read.csv("data/args_oap_term_resistome.csv")
term_type <- term_resistome%>%filter(group=="type")%>%dplyr::select(-group)
term_type[is.na(term_type)] <- 0
term_type$names[term_type$names=="macrolide-lincosamide-streptogramin"] <- "MLS"
term_type$names[term_type$names=="other_peptide_antibiotics"] <- "OPA"
term_type$names[term_type$names=="antibacterial_fatty_acid"] <- "AFA"
term_type$names[term_type$names=="pleuromutilin_tiamulin"] <- "PT"
term_type_sum <- term_type%>%
   dplyr::summarise("type"=names,
                    "n"=rowSums(.[-1]>0),
                    "prevalence"=rowSums(.[-1]>0)/(ncol(.)-1),
                    "Mean_copies_per_cell"=rowMeans(.[-1]))%>%
   arrange(desc(prevalence))

term_type_se <- data.frame("type"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in term_type$names) {
   a <- term_type%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%setNames("type")
   term_type_se[nrow(term_type_se)+1,] <- c(t,std.error(a$type))
}
term_type_sum_se <- left_join(term_type_sum,term_type_se,by="type")
term_type_sum_se$se <- as.numeric(term_type_sum_se$se)
term_type_sum_abund <- term_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"type_total_abund"=rowSums(.[-1]))
term_type_sum_abund <- left_join(term_type_sum_abund,metadata_cb%>%dplyr::select(SequenceID,Study),by=c("sample"="SequenceID"))
term_type_sum_abund$group <- "term"
term_type_sum_se$type <- factor(term_type_sum_se$type,levels = c("multidrug","tetracycline","beta_lactam","polymyxin","MLS","bacitracin",
                                                                 term_type_sum_se$type[7:26]))
#Supplementary Figure 5a
term_type_p <- ggplot() +
   geom_bar(data=term_type_sum_se, aes(x=type, y=Mean_copies_per_cell),stat="identity",fill="#E09137",color="black",width = 0.8)+
   geom_errorbar(data=term_type_sum_se,aes(x=type,ymin=Mean_copies_per_cell, ymax=Mean_copies_per_cell+se),width=0.5)+
   geom_line(data=term_type_sum_se, aes(x=type, y=prevalence/1.5,group=1),linetype = "dashed",color="#E09137",linewidth=0.5)+
   geom_point(data=term_type_sum_se, aes(x=type, y=prevalence/1.5),size=1.5)+
   geom_hline(yintercept=0.9/1.5,linetype="dashed",size=0.5,color="black")+
   scale_y_continuous(expand = c(0,0),name= "Mean ARG abundance (capc)",limits = c(0,0.9),
                      sec.axis = sec_axis(~ . *1.5, breaks = c(0,0.25,0.5,0.75,0.9,1)))+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         axis.text.y = element_text(size=10, colour = "black"),
         axis.text.x = element_text(size=10, colour="black",angle = 90,vjust = 0.4,hjust = 1),
         axis.title.y = element_text(size=10, colour = "black"),
         axis.title.x = element_blank())

term_subtype <- term_resistome%>%filter(group=="subtype")%>%dplyr::select(-group)
term_subtype_sum <- term_subtype%>%
   dplyr::summarise("subtype"=names,
                    "n"=rowSums(.[-1]>0),
                    "prevalence"=rowSums(.[-1]>0)/(ncol(.)-1),
                    "Mean_copies_per_cell"=rowMeans(.[-1]))%>%
   arrange(desc(prevalence))

#comparisons between preterm and term infants
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_type <- preterm_resistome%>%filter(group=="type")%>%dplyr::select(-group)
preterm_type[is.na(preterm_type)] <- 0
preterm_type$names[preterm_type$names=="macrolide-lincosamide-streptogramin"] <- "MLS"
preterm_type$names[preterm_type$names=="other_peptide_antibiotics"] <- "OPA"
preterm_type$names[preterm_type$names=="antibacterial_fatty_acid"] <- "AFA"
preterm_type$names[preterm_type$names=="pleuromutilin_tiamulin"] <- "PT"
preterm_type_sum <- preterm_type%>%
   dplyr::summarise("type"=names,
                    "n"=rowSums(.[-1]>0),
                    "prevalence"=rowSums(.[-1]>0)/(ncol(.)-1),
                    "Mean_copies_per_cell"=rowMeans(.[-1]))%>%
   arrange(desc(Mean_copies_per_cell))

preterm_type_se <- data.frame("type"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in preterm_type$names) {
   a <- preterm_type%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%setNames("type")
   preterm_type_se[nrow(preterm_type_se)+1,] <- c(t,std.error(a$type))
}
preterm_type_sum_se <- left_join(preterm_type_sum,preterm_type_se,by="type")
preterm_type_sum_se$se <- as.numeric(preterm_type_sum_se$se)
preterm_type_sum_richness <- preterm_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"richness"=rowSums(.[-1]>0))
preterm_type_sum_abund <- preterm_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"type_total_abund"=rowSums(.[-1]))
preterm_type_sum_abund <- left_join(preterm_type_sum_abund,metadata_cb%>%dplyr::select(SequenceID,Study),by=c("sample"="SequenceID"))
preterm_type_sum_abund$group <- "preterm"
preterm_type_sum_abund_term <- rbind(preterm_type_sum_abund,term_type_sum_abund)
preterm_type_sum_abund_term$group <- as.factor(preterm_type_sum_abund_term$group)
preterm_type_sum_abund_term$Study <- as.factor(preterm_type_sum_abund_term$Study)
wilcox_test(type_total_abund ~ group|Study, data=preterm_type_sum_abund_term)
term_preterm_type_pval <- data.frame("type"="x","term_mean"=0,"term_se"=0,"term_prevalence"=0,"preterm_mean"=0,"preterm_se"=0,"preterm_prevalence"=0,"pval"=0,stringsAsFactors = F)[-1,]
type_list <- unique(preterm_type$names)
count=0
for (t in type_list) {
   count=count+1
   print(count)
   db_term <- term_type%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
   if (ncol(db_term)==1) {
      db_term$type <- 0
      print(t)
   }
   db_term$group <- "term"
   names(db_term)[2] <- "type"
   db_preterm <- preterm_type%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
   db_preterm$group <- "preterm"
   names(db_preterm)[2] <- "type"
   db_cb <- rbind.fill(db_term,db_preterm)
   db_cb <- left_join(db_cb,metadata_cb%>%dplyr::select(SequenceID,Study),by="SequenceID")
   db_cb$Study <- as.factor(db_cb$Study)
   db_cb$group <- as.factor(db_cb$group)
   pval <- pvalue(wilcox_test(type ~ group|Study, data=db_cb))
   term_preterm_type_pval[nrow(term_preterm_type_pval)+1,] <- c(t,mean(db_term$type),std.error(db_term$type),sum(db_term$type>0)/1356,
                                                                mean(db_preterm$type),std.error(db_preterm$type),sum(db_preterm$type>0)/3302,pval)
}
term_preterm_type_pval <- term_preterm_type_pval%>%mutate_at(vars(-type),as.numeric)
term_preterm_type_pval$pval[term_preterm_type_pval$pval=="NaN"] <- NA
term_preterm_type_pval$fdr <- p.adjust(term_preterm_type_pval$pval,method = "fdr")
term_preterm_type_pval_sig <- term_preterm_type_pval%>%filter(fdr<0.05)
term_preterm_type_pval_sig$diff <- term_preterm_type_pval_sig$preterm_mean-term_preterm_type_pval_sig$term_mean

term_preterm_type_pval_sel <- term_preterm_type_pval%>%filter(fdr<0.05)%>%arrange(desc(preterm_mean))
term_preterm_type_sum <- rbind(term_type_sum_se%>%mutate(group="full_term")%>%filter(type%in%term_preterm_type_pval_sel$type),
                               preterm_type_sum_se%>%mutate(group="preterm")%>%filter(type%in%term_preterm_type_pval_sel$type))
term_preterm_type_sum$type <- factor(term_preterm_type_sum$type,levels = term_preterm_type_pval_sel$type)
term_preterm_type_sum$group <- factor(term_preterm_type_sum$group,levels = c("preterm","full_term"))
term_preterm_type_sum$group2 <- paste0(term_preterm_type_sum$group,"_",term_preterm_type_sum$type)
term_preterm_type_sum$Mean_copies_per_cell[term_preterm_type_sum$group2=="preterm_multidrug"]#3.549179
term_preterm_type_sum$Mean_copies_per_cell[term_preterm_type_sum$group2=="preterm_multidrug"] <- 1.2
#Figure 2a
term_preterm_type_p <- ggplot(term_preterm_type_sum, aes(x=type, y=Mean_copies_per_cell,fill=group)) +
   geom_bar(stat="identity",position=position_dodge(),color="black",width = 0.8)+
   geom_errorbar(aes(ymin=Mean_copies_per_cell, ymax=Mean_copies_per_cell+se),position=position_dodge(0.8),width=0.5)+
   scale_fill_manual(values = c("#4994C4","#E09137"),labels=c("preterm"="Preterm","full_term"="Full term"))+
   scale_y_continuous(expand = c(0,0),name= "Mean abundance (capc)",limits = c(0,1.3),
                      breaks = c(0,0.25,0.5,0.75,1,1.2),labels = c("0","0.25","0.50","0.75","1.00","3.55"))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90,vjust = 0.4,hjust = 1),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_blank())

term_preterm_type_sum_2 <- term_preterm_type_sum%>%filter(type%in%term_preterm_type_pval_sel$type[13:19])
#Figure 2a
term_preterm_type_p_2 <- ggplot(term_preterm_type_sum_2, aes(x=type, y=Mean_copies_per_cell,fill=group)) +
   geom_bar(stat="identity",position=position_dodge(),color="black",width = 0.8)+
   geom_errorbar(aes(ymin=Mean_copies_per_cell, ymax=Mean_copies_per_cell+se),position=position_dodge(0.8),width=0.5)+
   scale_fill_manual(values = c("#4994C4","#E09137"),labels=c("preterm"="Preterm","full_term"="Full term"))+
   scale_y_continuous(expand = c(0,0),name= "Mean abundance (capc)",limits = c(0,0.023),
                      breaks = c(0,0.005,0.010,0.015,0.020),labels = c("0","0.005","0.010","0.015","0.020"))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

#alpha diversity based on ARG subtypes
term_subtype_t <- term_subtype%>%column_to_rownames("names")%>%t()%>%as.data.frame()
term_subtype_richness <- rowSums(term_subtype_t != 0)
term_subtype_shannon <- diversity(term_subtype_t, index = "shannon", MARGIN=1) 
term_subtype_richness_df <- term_subtype_richness%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(term_subtype_richness_df)[2] <- "Richness"
term_subtype_shannon_df <- term_subtype_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(term_subtype_shannon_df)[2] <- "Shannon"
term_subtype_alpha <- left_join(term_subtype_richness_df,term_subtype_shannon_df,by="SequenceID")

preterm_subtype <- preterm_resistome%>%filter(group=="subtype")%>%dplyr::select(-group)
preterm_subtype_sum <- preterm_subtype%>%
   dplyr::summarise("subtype"=names,
                    "n"=rowSums(.[-1]>0),
                    "prevalence"=rowSums(.[-1]>0)/(ncol(.)-1),
                    "Mean_copies_per_cell"=rowMeans(.[-1]))%>%
   arrange(desc(prevalence))
length(which(!term_subtype_sum$subtype%in%preterm_subtype_sum$subtype))
#48 ARG subtypes are exclusively in full-term infants.
length(which(!preterm_subtype_sum$subtype%in%term_subtype_sum$subtype))
#567 ARG subtypes are exclusively in preterm infants.
preterm_subtype_t <- preterm_subtype%>%column_to_rownames("names")%>%t()%>%as.data.frame()
preterm_subtype_richness <- rowSums(preterm_subtype_t != 0)
preterm_subtype_shannon <- diversity(preterm_subtype_t, index = "shannon", MARGIN=1) 
preterm_subtype_richness_df <- preterm_subtype_richness%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(preterm_subtype_richness_df)[2] <- "Richness"
preterm_subtype_shannon_df <- preterm_subtype_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(preterm_subtype_shannon_df)[2] <- "Shannon"
preterm_subtype_alpha <- left_join(preterm_subtype_richness_df,preterm_subtype_shannon_df,by="SequenceID")

subtype_alpha <- rbind(preterm_subtype_alpha,term_subtype_alpha)
subtype_alpha <- left_join(subtype_alpha,metadata_cb,by="SequenceID")
subtype_alpha$Study <- as.factor(subtype_alpha$Study)
subtype_alpha$month <- as.factor(subtype_alpha$month)
subtype_alpha$Term <- factor(subtype_alpha$Term,levels = c("preterm","full_term"))
wilcox_test(Richness ~ Term | Study ,data=subtype_alpha)
wilcox_test(Shannon ~ Term | Study ,data=subtype_alpha)
wilcox_test(Richness ~ Term | month ,data=subtype_alpha)
wilcox_test(Shannon ~ Term | month ,data=subtype_alpha)
#Figure 2b
preterm_full_richness <- ggplot(subtype_alpha, aes(x=Term, y=Richness)) +
   geom_violin(aes(fill=Term),alpha=1,trim=T)+
   geom_boxplot(width=0.1,outlier.size = 1)+
   scale_fill_manual(values = c("#4994C4","#E09137"))+
   scale_x_discrete(labels=c("preterm"="Preterm","full_term"="Full term"))+
   labs(x ="", y="Richness (number of ARG subtypes)")+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_blank())
#Figure 2b
preterm_full_richness_day <- ggplot(subtype_alpha, aes(x=DOL, y=Richness,color=Term)) + 
   geom_point(alpha=0.5,size=0.8)+
   geom_smooth(method = loess,size=1.5,level=0.95)+
   scale_color_manual(values = c("#4994C4","#E09137"),labels=c("preterm"="Preterm","full_term"="Full term"))+
   guides(color = guide_legend(reverse=T))+
   labs(x ="Age (days)", y="Richness (number of ARG subtypes)")+
   scale_x_continuous(labels=c("0","300","600","900","1,200"))+
   theme_bw()+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid.minor = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour="black"),
         axis.ticks.y = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=12, colour = "black"))
#Supplementary Figure 5b
preterm_full_shannon <- ggplot(subtype_alpha, aes(x=Term, y=Shannon)) +
   geom_violin(aes(fill=Term),alpha=1,trim=T)+
   geom_boxplot(width=0.1,outlier.size = 1)+
   scale_fill_manual(values = c("#4994C4","#E09137"))+
   scale_x_discrete(labels=c("preterm"="Preterm","full_term"="Full term"))+
   labs(x ="", y="Shannon index (based on ARG subtypes)")+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_blank())
#Supplementary Figure 5b
preterm_full_shannon_day <- ggplot(subtype_alpha, aes(x=DOL, y=Shannon,color=Term)) + 
   geom_point(alpha=0.5,size=0.8)+
   geom_smooth(method = loess,size=1.5,level=0.95)+
   scale_color_manual(values = c("#4994C4","#E09137"),labels=c("preterm"="Preterm","full_term"="Full term"))+
   guides(color = guide_legend(reverse=T))+
   labs(x ="Age (days)", y="Richness (number of ARG subtypes)")+
   scale_x_continuous(labels=c("0","300","600","900","1,200"))+
   theme_bw()+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid.minor = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour="black"),
         axis.ticks.y = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=12, colour = "black"))
#differential subtypes
term_preterm_subtype_pval <- data.frame("subtype"="x","term_mean"=0,"term_se"=0,"term_prevalence"=0,"preterm_mean"=0,"preterm_se"=0,"preterm_prevalence"=0,"pval"=0,stringsAsFactors = F)[-1,]
subtype_list <- c(term_subtype$names,preterm_subtype$names[-which(preterm_subtype$names%in%term_subtype$names)])#2041
count=0
for (t in subtype_list) {
   # t="beta_lactam__CfxA"
   count=count+1
   print(count)
   db_term <- term_subtype%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
   if (ncol(db_term)==1) {
      db_term$subtype <- 0
      print(t)
   }
   db_term$group <- "term"
   names(db_term)[2] <- "subtype"
   db_preterm <- preterm_subtype%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
   if (ncol(db_preterm)==1) {
      db_preterm$subtype <- 0
      print(t)
   }
   db_preterm$group <- "preterm"
   names(db_preterm)[2] <- "subtype"
   db_cb <- rbind.fill(db_term,db_preterm)
   db_cb <- left_join(db_cb,metadata_cb%>%dplyr::select(SequenceID,Study),by="SequenceID")
   db_cb$Study <- as.factor(db_cb$Study)
   db_cb$group <- as.factor(db_cb$group)
   pval <- pvalue(wilcox_test(subtype ~ group|Study, data=db_cb))
   term_preterm_subtype_pval[nrow(term_preterm_subtype_pval)+1,] <- c(t,mean(db_term$subtype),std.error(db_term$subtype),sum(db_term$subtype>0)/1356,
                                                                      mean(db_preterm$subtype),std.error(db_preterm$subtype),sum(db_preterm$subtype>0)/3302,
                                                                      pval)
}
term_preterm_subtype_pval_new <- term_preterm_subtype_pval%>%mutate_at(vars(-subtype),as.numeric)
term_preterm_subtype_pval_new$pval[term_preterm_subtype_pval_new$pval=="NaN"] <- NA
term_preterm_subtype_pval_new_sel <- term_preterm_subtype_pval_new%>%filter(term_prevalence>0.05 | preterm_prevalence>0.05)
term_preterm_subtype_pval_new_sel$fdr <- p.adjust(term_preterm_subtype_pval_new_sel$pval,method = "fdr")
term_preterm_subtype_pval_new_sel$logfdr <- -log10(term_preterm_subtype_pval_new_sel$fdr)
term_preterm_subtype_pval_new_sel$logfdr[term_preterm_subtype_pval_new_sel$logfdr=="Inf"] <- 20
term_preterm_subtype_pval_new_sel$logfold <- log2(term_preterm_subtype_pval_new_sel$preterm_mean/term_preterm_subtype_pval_new_sel$term_mean)
term_preterm_subtype_pval_new_sig <- term_preterm_subtype_pval_new_sel%>%filter(fdr<0.05)
term_preterm_subtype_pval_new_sig$diff <- term_preterm_subtype_pval_new_sig$preterm_mean-term_preterm_subtype_pval_new_sig$term_mean
term_preterm_subtype_pval_new2 <- term_preterm_subtype_pval_new%>%filter(term_mean==0)
term_preterm_subtype_pval_new3 <- term_preterm_subtype_pval_new%>%filter(preterm_mean==0)
term_preterm_subtype_pval_new_sel$type <- sapply(strsplit(term_preterm_subtype_pval_new_sel$subtype, split='__', fixed=TRUE), function(x)(x[1]))
term_preterm_subtype_pval_new_sel$type[term_preterm_subtype_pval_new_sel$type=="macrolide-lincosamide-streptogramin"] <- "MLS"
term_preterm_subtype_pval_new_sel$type[term_preterm_subtype_pval_new_sel$type=="other_peptide_antibiotics"] <- "OPA"
term_preterm_subtype_pval_new_sel$type[term_preterm_subtype_pval_new_sel$type=="antibacterial_fatty_acid"] <- "AFA"
term_preterm_subtype_pval_new_sel$gene <- sapply(strsplit(term_preterm_subtype_pval_new_sel$subtype, split='__', fixed=TRUE), function(x)(x[2]))
term_preterm_subtype_pval_new_sel <- term_preterm_subtype_pval_new_sel%>%
   mutate("sig"=case_when(fdr<0.05 & logfold>0 ~ "Up",
                          fdr<0.05 & logfold<0 ~ "Down",
                          fdr>=0.05 ~ "Nosig"))%>%
   mutate("gene2"=case_when(fdr<0.05 ~ gene,
                            fdr>=0.05 ~ NA))%>%
   mutate("type2"=case_when(fdr<0.05 ~ type,
                            fdr>=0.05 ~ NA))%>%
   mutate("prevalence"=case_when(sig=="Up" ~ preterm_prevalence,
                                 sig=="Down" ~ term_prevalence,
                                 sig=="Nosig" ~ 0.2))%>%
   mutate("shape"=case_when(sig=="Up" ~ "Up",
                            sig=="Down" ~ "Down",
                            sig=="Nosig" ~ NA))
term_preterm_subtype_pval_new_sel$sig <- factor(term_preterm_subtype_pval_new_sel$sig,levels = c("Nosig","Up","Down"))
term_preterm_subtype_pval_new_sel$type2 <- as.factor(term_preterm_subtype_pval_new_sel$type2)
#Figure 2c
term_preterm_subtype_sig_p <- ggplot() +
   geom_point(data=term_preterm_subtype_pval_new_sel, aes(x=logfold,y=logfdr,size=prevalence,shape=sig,color=type2),stroke=2,alpha=1)+
   scale_shape_manual(values = c(16,16,10))+
   scale_color_manual(values = c("#684795","#386CB0","#7FC97F","#43760A","#FCCDE5","#E6AB02","#D9D9D9","#DA732D","#BC80BD",
                                 "#934B07","#BEAED4","#1B9E77","#E7298A","#F3A2CD","#D95F02","#80B1D3",
                                 "#FDB462","#C4C462","#800404","#F54646"))+
   geom_hline(yintercept=-log10(0.05),linetype="dashed",size=1,color="grey") +
   geom_vline(xintercept=0,linetype="dashed",size=1,color="grey") +
   labs(x ="log2FC", y="log10FDR")+
   theme_bw()+
   guides(color = guide_legend(override.aes = list(size=3)),
          shape = guide_legend(override.aes = list(size=5)))+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1.5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=12, colour = "black"),
         axis.text.y = element_text(size=20, colour="black"),
         axis.text.x = element_text(size=20, colour="black"),
         axis.title = element_text(size=20, colour = "black"))

#beta diversity based on ARG subtypes
subtype_t <- rbind.fill(preterm_subtype_t%>%rownames_to_column("run"),term_subtype_t%>%rownames_to_column("run"))
subtype_t[is.na(subtype_t)] <- 0
subtype_t <- subtype_t%>%column_to_rownames("run")
subtype_t_beta_dist <- vegdist(subtype_t, method = "bray")
subtype_t_beta_dist_db <- subtype_t_beta_dist%>%as.matrix()%>%as.data.frame()
subtype_t_beta_dist_db[upper.tri(subtype_t_beta_dist_db, diag=T)] <- NA
subtype_t_beta_dist_db <- subtype_t_beta_dist_db%>%rownames_to_column("run1")
subtype_t_beta_dist_db_g <- gather(subtype_t_beta_dist_db,run2,distance,-run1)
subtype_t_beta_dist_db_g_2 <- subtype_t_beta_dist_db_g%>%filter(!is.na(distance))
subtype_t_beta_dist_db_g_2 <- left_join(subtype_t_beta_dist_db_g_2,metadata_cb%>%dplyr::select(SequenceID,Term),by=c("run1"="SequenceID"))
subtype_t_beta_dist_db_g_2 <- left_join(subtype_t_beta_dist_db_g_2,metadata_cb%>%dplyr::select(SequenceID,Term),by=c("run2"="SequenceID"))
subtype_t_beta_dist_db_g_2_sel <- subtype_t_beta_dist_db_g_2%>%filter(Term.x==Term.y)
subtype_t_beta_dist_db_g_2_sel$Term.x <- factor(subtype_t_beta_dist_db_g_2_sel$Term.x,levels = c("preterm","full_term"))
wilcox_test(distance ~ Term.x,data = subtype_t_beta_dist_db_g_2_sel)
#Figure 2e
subytpe_distance <- ggplot(subtype_t_beta_dist_db_g_2_sel, aes(x=Term.x, y=distance, fill=Term.x)) +
   geom_boxplot(outlier.size = 1)+
   scale_fill_manual(values = c("#4994C4","#E09137"))+
   scale_x_discrete(labels=c("preterm"="Preterm","full_term"="Full term"))+
   labs(x ="", y="Bray-Curtis distance based on ARG subtypes")+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_blank())

identical(metadata_cb$SequenceID,rownames(subtype_t))
metadata_cb$Term <- as.factor(metadata_cb$Term)
set.seed(1000)
with(metadata_cb,vegan::adonis2(subtype_t_beta_dist ~ Term,data=metadata_cb,permutations = 1000,strata = Study))
subtype_beta_pco <- cmdscale(subtype_t_beta_dist, k=2, eig = T)
subtype_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                       round((subtype_beta_pco$eig[1]/sum(subtype_beta_pco$eig))*100,1),
                                       '%]', sep='')
subtype_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                       round((subtype_beta_pco$eig[2]/sum(subtype_beta_pco$eig))*100,1),
                                       '%]', sep='')
subtype_beta_pco.point <- as.data.frame(subtype_beta_pco$points)%>%rownames_to_column('SequenceID')
subtype_beta_pco.tabble <- right_join(metadata_cb, subtype_beta_pco.point, by="SequenceID")
subtype_beta_df.plot <- tibble(Axis1 = subtype_beta_pco.tabble$V1,
                               Axis2 = subtype_beta_pco.tabble$V2,
                               Sample_ID = subtype_beta_pco.tabble$SubjectID,
                               Term=subtype_beta_pco.tabble$Term)
#Figure 2d
subtype_beta_time_p <- ggplot(data=subtype_beta_df.plot,aes(x=Axis1, y=Axis2, col=Term)) +
   geom_point(size=1.5, alpha=1) + 
   scale_color_manual(values = c("#E09137","#4994C4"),labels=c("Full term","Preterm"),guide = guide_legend(reverse = TRUE))+
   theme_bw()+ 
   xlab(subtype_beta_pco_axis.1.title) + ylab(subtype_beta_pco_axis.2.title) +
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1),"mm"))

#rgi
rgi_term_profile_coverm <- read.csv("data/rgi_term_resistome_coverm_bins_gtdb_stats_plasmid.csv")
length(unique(rgi_term_profile_coverm$Best_Hit_ARO))#1352
#aro
rgi_term_profile_coverm_aro <- rgi_term_profile_coverm%>%group_by(sample,Best_Hit_ARO)%>%dplyr::summarise(rpkm_total=sum(rpkm))
rgi_term_profile_coverm_aro_richness <- rgi_term_profile_coverm_aro%>%group_by(sample)%>%dplyr::summarise(n=n())
rgi_term_profile_coverm_aro_sum <- rgi_term_profile_coverm_aro%>%
   group_by(Best_Hit_ARO)%>%
   dplyr::summarise("n"=n(),
                    "prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(rpkm_total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))

#organize the function
split_into_multiple <- function(column, pattern, into_prefix){
   cols <- str_split_fixed(column, pattern, n = Inf)
   cols[which(cols == "")] <- NA
   cols <- tibble::as.tibble(cols)
   m <- dim(cols)[2]
   names(cols) <- paste(into_prefix, 1:m, sep = "_")
   return(cols)
}
rgi_term_profile_coverm_gene <- rgi_term_profile_coverm %>%
   dplyr::select(cb,sample,rpkm,AMR.Gene.Family)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(cb,sample,rpkm,AMR.Gene.Family))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)
length(unique(rgi_term_profile_coverm_gene$gene))#276 gene families
rgi_term_profile_coverm_gene_sum <- rgi_term_profile_coverm_gene%>%
   group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%spread(.,sample,rpkm_total)
rgi_term_profile_coverm_gene_sum[is.na(rgi_term_profile_coverm_gene_sum)] <- 0
rgi_term_profile_coverm_gene_sum_t <- rgi_term_profile_coverm_gene_sum%>%column_to_rownames("gene")%>%t()%>%as.data.frame()
rgi_term_profile_coverm_gene_sum2 <- rgi_term_profile_coverm_gene%>%group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))
rgi_term_richness <- rgi_term_profile_coverm_gene_sum2%>%group_by(sample)%>%dplyr::summarise(n=n())
rgi_term_shannon <- diversity(rgi_term_profile_coverm_gene_sum_t, index = "shannon", MARGIN=1) 
names(rgi_term_richness) <- c("SequenceID","Richness")
rgi_term_shannon_df <- rgi_term_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(rgi_term_shannon_df)[2] <- "Shannon"
rgi_term_alpha <- left_join(rgi_term_richness,rgi_term_shannon_df,by="SequenceID")

rgi_preterm_profile_coverm <- read.csv("data/rgi_preterm_resistome_coverm_bins_gtdb_stats_plasmid.csv")
rgi_preterm_profile_coverm_gene <- rgi_preterm_profile_coverm %>%
   dplyr::select(cb,sample,rpkm,AMR.Gene.Family)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(cb,sample,rpkm,AMR.Gene.Family))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)
length(unique(rgi_preterm_profile_coverm_gene$gene))#282 gene families
rgi_preterm_profile_coverm_gene_sum <- rgi_preterm_profile_coverm_gene%>%
   group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%spread(.,sample,rpkm_total)
rgi_preterm_profile_coverm_gene_sum[is.na(rgi_preterm_profile_coverm_gene_sum)] <- 0
rgi_preterm_profile_coverm_gene_sum_t <- rgi_preterm_profile_coverm_gene_sum%>%column_to_rownames("gene")%>%t()%>%as.data.frame()
rgi_preterm_profile_coverm_gene_sum2 <- rgi_preterm_profile_coverm_gene%>%
   group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))
rgi_preterm_richness <- rgi_preterm_profile_coverm_gene_sum2%>%group_by(sample)%>%dplyr::summarise(n=n())
rgi_preterm_shannon <- diversity(rgi_preterm_profile_coverm_gene_sum_t, index = "shannon", MARGIN=1) 
names(rgi_preterm_richness) <- c("SequenceID","Richness")
rgi_preterm_shannon_df <- rgi_preterm_shannon%>%as.data.frame()%>%rownames_to_column("SequenceID")
names(rgi_preterm_shannon_df)[2] <- "Shannon"
rgi_preterm_alpha <- left_join(rgi_preterm_richness,rgi_preterm_shannon_df,by="SequenceID")

rgi_alpha <- rbind(rgi_preterm_alpha,rgi_term_alpha)
rgi_alpha <- left_join(rgi_alpha,metadata_cb,by="SequenceID")
rgi_alpha$Study <- as.factor(rgi_alpha$Study)
rgi_alpha$month <- as.factor(rgi_alpha$month)
rgi_alpha$Term <- factor(rgi_alpha$Term,levels = c("preterm","full_term"))
wilcox_test(Richness ~ Term | Study ,data=rgi_alpha)
wilcox_test(Shannon ~ Term | Study ,data=rgi_alpha)
wilcox_test(Richness ~ Term |month,data=rgi_alpha)
wilcox_test(Shannon ~ Term |month,data=rgi_alpha)
#Supplementary Figure 6a
preterm_full_richness_rgi <- ggplot(rgi_alpha, aes(x=Term, y=Richness)) +
   geom_violin(aes(fill=Term),alpha=1,trim=T)+
   geom_boxplot(width=0.1,outlier.size = 1)+
   scale_fill_manual(values = c("#4994C4","#E09137"))+
   scale_x_discrete(labels=c("preterm"="Preterm","full_term"="Full term"))+
   labs(x ="", y="Richness (number of ARG familes)")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_blank())
#Supplementary Figure 6a
preterm_full_richness_day_rgi <- ggplot(rgi_alpha, aes(x=DOL, y=Richness,color=Term)) + 
   geom_point(alpha=0.5,size=0.8)+
   geom_smooth(method = loess,size=1.5,level=0.95)+
   scale_color_manual(values = c("#4994C4","#E09137"),labels=c("preterm"="Preterm","full_term"="Full term"))+
   guides(color = guide_legend(reverse=T))+
   labs(x ="Age (days)", y="Richness (number of ARG families)")+
   scale_x_continuous(labels=c("0","300","600","900","1,200"))+
   theme_bw()+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour="black"),
         axis.ticks.y = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=12, colour = "black"))
#Supplementary Figure 6b
preterm_full_shannon_rgi <- ggplot(rgi_alpha, aes(x=Term, y=Shannon)) +
   geom_violin(aes(fill=Term),alpha=1,trim=T)+
   geom_boxplot(width=0.1,outlier.size = 1)+
   scale_fill_manual(values = c("#4994C4","#E09137"))+
   scale_x_discrete(labels=c("preterm"="Preterm","full_term"="Full term"))+
   labs(x ="", y="Shannon index (based on ARG families)")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_blank())
#Supplementary Figure 6b
preterm_full_shannon_day <- ggplot(rgi_alpha, aes(x=DOL, y=Shannon,color=Term)) + 
   geom_point(alpha=0.5,size=0.8)+
   geom_smooth(method = loess,size=1.5,level=0.95)+
   scale_color_manual(values = c("#4994C4","#E09137"),labels=c("preterm"="Preterm","full_term"="Full term"))+
   guides(color = guide_legend(reverse=T))+
   labs(x ="Age (days)", y="Richness (number of ARG families)")+
   scale_x_continuous(labels=c("0","300","600","900","1,200"))+
   theme_bw()+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=12, colour="black"),
         axis.ticks.y = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=12, colour = "black"))

#beta diversity based on rgi families
rgi_t <- rbind.fill(rgi_preterm_profile_coverm_gene_sum_t%>%rownames_to_column("run"),rgi_term_profile_coverm_gene_sum_t%>%rownames_to_column("run"))
rgi_t[is.na(rgi_t)] <- 0
rgi_t <- rgi_t%>%column_to_rownames("run")
rgi_t <- rgi_t[which(rowSums(rgi_t)!=0),which(colSums(rgi_t)!=0)]
rgi_t_beta_dist <- vegdist(rgi_t, method = "bray")
metadata_cb_sel <- metadata_cb%>%filter(SequenceID%in%rownames(rgi_t))%>%arrange(match(SequenceID,rownames(rgi_t)))
identical(rownames(rgi_t),metadata_cb_sel$SequenceID)
metadata_cb_sel$Term <- as.factor(metadata_cb_sel$Term)
set.seed(1000)
with(metadata_cb_sel,vegan::adonis2(rgi_t_beta_dist ~ Term,data=metadata_cb_sel,permutations = 1000,strata = Study))
rgi_beta_pco <- cmdscale(rgi_t_beta_dist, k=2, eig = T)
rgi_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                   round((rgi_beta_pco$eig[1]/sum(rgi_beta_pco$eig))*100,1),
                                   '%]', sep='')
rgi_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                   round((rgi_beta_pco$eig[2]/sum(rgi_beta_pco$eig))*100,1),
                                   '%]', sep='')
rgi_beta_pco.point <- as.data.frame(rgi_beta_pco$points)%>%rownames_to_column('SequenceID')
rgi_beta_pco.tabble <- right_join(metadata_cb, rgi_beta_pco.point, by="SequenceID")
rgi_beta_df.plot <- tibble(Axis1 = rgi_beta_pco.tabble$V1,
                           Axis2 = rgi_beta_pco.tabble$V2,
                           Sample_ID = rgi_beta_pco.tabble$SubjectID,
                           Term=rgi_beta_pco.tabble$Term)
#Supplementary Figure 6c
rgi_beta_time_p <- ggplot(data=rgi_beta_df.plot,aes(x=Axis1, y=Axis2, col=Term)) +
   geom_point(size=1.5, alpha=1) + 
   scale_color_manual(values = c("#E09137","#4994C4"),labels=c("Full term","Preterm"),guide = guide_legend(reverse = TRUE))+
   theme_bw()+ 
   xlab(rgi_beta_pco_axis.1.title) + ylab(rgi_beta_pco_axis.2.title) +
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = unit(c(1,5,1,1), "mm"))
#rgi drug class
rgi_term_profile_coverm_2 <- rgi_term_profile_coverm %>%
   dplyr::select(cb,sample,rpkm,Drug.Class)%>%
   dplyr::bind_cols(split_into_multiple(.$Drug.Class, "; ", "cat"))%>%
   gather(.,cat,drug,-c(cb,sample,rpkm,Drug.Class))%>%
   filter(!is.na(drug))%>%
   dplyr::select(-cat)
rgi_term_profile_coverm_2_sum <- rgi_term_profile_coverm_2%>%group_by(sample,drug)%>%dplyr::summarise(rpkm_total=sum(rpkm))
length(unique(rgi_term_profile_coverm_2_sum$drug))#36 drug class
rgi_term_profile_coverm_3 <- rgi_term_profile_coverm_2_sum%>%
   group_by(drug)%>%
   dplyr::summarise("prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(rpkm_total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))

rgi_term_profile_coverm_2_sum_s <- spread(rgi_term_profile_coverm_2_sum,sample,rpkm_total)
rgi_term_profile_coverm_2_sum_s[is.na(rgi_term_profile_coverm_2_sum_s)] <- 0
term_rgi_se <- data.frame("drug"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in rgi_term_profile_coverm_2_sum_s$drug) {
   a <- rgi_term_profile_coverm_2_sum_s%>%filter(drug==t)%>%column_to_rownames("drug")%>%t()%>%as.data.frame()%>%setNames("drug")
   term_rgi_se[nrow(term_rgi_se)+1,] <- c(t,std.error(a$drug))
}
rgi_term_profile_coverm_3_se <- left_join(rgi_term_profile_coverm_3,term_rgi_se,by="drug")
rgi_term_profile_coverm_3_se$se <- as.numeric(rgi_term_profile_coverm_3_se$se)
rgi_term_profile_coverm_3_se$drug <- gsub(" antibiotic","",rgi_term_profile_coverm_3_se$drug,fixed = T)
rgi_term_profile_coverm_3_se$drug[rgi_term_profile_coverm_3_se$drug=="disinfecting agents and antiseptics"] <- "DAA"
rgi_term_profile_coverm_3_se$drug[rgi_term_profile_coverm_3_se$drug=="antibacterial free fatty acids"] <- "AFFA"
rgi_term_profile_coverm_3_se$drug <- factor(rgi_term_profile_coverm_3_se$drug,levels = rgi_term_profile_coverm_3_se$drug)
rgi_term_profile_coverm_3_se$Mean_rpkm[rgi_term_profile_coverm_3_se$drug=="tetracycline"]#328940
rgi_term_profile_coverm_3_se$Mean_rpkm[rgi_term_profile_coverm_3_se$drug=="tetracycline"] <- 220000
#Supplementary Figure 5c
rgi_term_drug_class_p <- ggplot() +
   geom_bar(data=rgi_term_profile_coverm_3_se, aes(x=drug, y=Mean_rpkm),stat="identity",fill="#E09137",color="black",width = 0.8)+
   geom_errorbar(data=rgi_term_profile_coverm_3_se,aes(x=drug,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=rgi_term_profile_coverm_3_se, aes(x=drug, y=prevalence/0.000006,group=1),linetype = "dashed",color="#E09137",size=0.5)+
   geom_point(data=rgi_term_profile_coverm_3_se, aes(x=drug, y=prevalence/0.000006),size=1.5)+
   geom_hline(yintercept=0.9/0.000006,linetype="dashed",size=0.5,color="black")+
   scale_y_continuous(expand = c(0,0),name= "Mean ARG family abundance (rpkm)",limits = c(0,240000),
                      breaks = c(0,50000,100000,150000,220000),labels = c("0","50,000","100,000","150,000","328,940"),
                      sec.axis = sec_axis(~ . *0.000006, breaks = c(0,0.25,0.5,0.75,0.9,1)))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         plot.background = element_blank(),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90,vjust = 0.4,hjust = 1),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

####################################################################################################################################
#ARGs host from full-term infants
####################################################################################################################################
std.error <- function(x) sd(x)/sqrt(length(x))
rgi_term_profile <- read.csv("data/rgi_term_resistome_coverm_bins_gtdb_stats_plasmid.csv")
rgi_term_profile$bin[is.na(rgi_term_profile$bin)] <- "nobin" 
rgi_term_profile_nobin <- rgi_term_profile%>%filter(bin=="nobin")
rgi_preterm_profile <- read.csv("data/rgi_preterm_resistome_coverm_bins_gtdb_stats_plasmid.csv")
plasmid_term_preterm <- rbind(rgi_preterm_profile%>%dplyr::select(plasmid)%>%mutate(group="preterm"),
                              rgi_term_profile%>%dplyr::select(plasmid)%>%mutate(group="term"))
chisq.test(plasmid_term_preterm$plasmid,plasmid_term_preterm$group,correct = T)
rgi_term_profile$domain <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[1]))
rgi_term_profile$phylum <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[2]))
rgi_term_profile$class <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[3]))
rgi_term_profile$order <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[4]))
rgi_term_profile$family <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[5]))
rgi_term_profile$genus <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[6]))
rgi_term_profile$species <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[7]))
rgi_term_profile$domain[is.na(rgi_term_profile$domain)] <- "unknown"
rgi_term_profile$phylum[is.na(rgi_term_profile$phylum)] <- "unknown"
rgi_term_profile$class[is.na(rgi_term_profile$class)] <- "unknown"
rgi_term_profile$order[is.na(rgi_term_profile$order)] <- "unknown"
rgi_term_profile$family[is.na(rgi_term_profile$family)] <- "unknown"
rgi_term_profile$genus[is.na(rgi_term_profile$genus)] <- "unknown"
rgi_term_profile$species[rgi_term_profile$species=="s__"] <- "unknown"
rgi_term_profile$species[is.na(rgi_term_profile$species)] <- "unknown"
rgi_term_profile_aa <- rgi_term_profile%>%filter(species=="s__")#109
rgi_term_profile_domain <- rgi_term_profile%>%group_by(domain)%>%dplyr::summarise(n=n())
tax_number <- data.frame("tax"="x","n_tax"=0,"n_args"=0,stringsAsFactors = F)[-1,]
tax_number[1,] <- c("Phylum",length(table(rgi_term_profile$phylum))-1,sum(table(rgi_term_profile$phylum))-nrow(rgi_term_profile%>%filter(phylum=="unknown")))
tax_number[2,] <- c("Class",length(table(rgi_term_profile$class))-1,sum(table(rgi_term_profile$class))-nrow(rgi_term_profile%>%filter(class=="unknown")))
tax_number[3,] <- c("Order",length(table(rgi_term_profile$order))-1,sum(table(rgi_term_profile$order))-nrow(rgi_term_profile%>%filter(order=="unknown")))
tax_number[4,] <- c("Family",length(table(rgi_term_profile$family))-1,sum(table(rgi_term_profile$family))-nrow(rgi_term_profile%>%filter(family=="unknown")))
tax_number[5,] <- c("Genus",length(table(rgi_term_profile$genus))-1,sum(table(rgi_term_profile$genus))-nrow(rgi_term_profile%>%filter(genus=="unknown")))
tax_number[6,] <- c("Species",length(table(rgi_term_profile$species))-1,sum(table(rgi_term_profile$species))-nrow(rgi_term_profile%>%filter(species=="unknown")))
tax_number$n_tax <- as.numeric(tax_number$n_tax)
tax_number$n_args <- as.numeric(tax_number$n_args)
tax_number$tax <- factor(tax_number$tax,levels = c("Phylum", "Class","Order","Family","Genus","Species"))
#Supplementary Figure 5d
tax_number_plot <- ggplot(tax_number, aes(x=tax, y=n_tax)) + 
   geom_bar(stat = "identity",fill="#E09137",color="black",width = 0.8)+
   geom_text(aes(label=paste0("(n=",formattable::comma(n_tax,digits=0),"; n=",formattable::comma(n_args,digits=0),")"),hjust=-0.1,size=3))+
   theme_bw() + # remove the backgroud
   labs(x = element_blank(), y=element_blank())+
   theme(axis.title.y = element_text(size = 16))+
   theme(panel.border = element_blank(), axis.line = element_blank()) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the grid
   coord_flip()+
   theme(legend.position="none",
         #aspect.ratio = 1.5/1,
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_text(size=14, colour = "black"),
         axis.text.x = element_blank(),
         axis.ticks = element_blank())

tax=c("phylum","class","order","family","genus")
tax_diver <- data.frame()
for (t in tax) {
   #t="domain"
   a <- rgi_term_profile%>%dplyr::group_by(get(t))%>%dplyr::summarise(n=n())%>%arrange(desc(n))
   print(sum(a$n))
   b <- a[c(1:11),] 
   b$propor <- b$n/sum(a$n)
   b[nrow(b)+1,] <- list("Others",nrow(rgi_term_profile)-sum(b$n,na.rm = T),1-sum(b$propor,na.rm = T))
   b$taxonomy <- t
   tax_diver <- rbind(tax_diver,b)
}
names(tax_diver)[1] <- "name" 
tax_diver$name[tax_diver$name=="unknown"] <- "Unknown"
tax_diver$propor <- as.numeric(tax_diver$propor)*100
tax_diver$propor <- round(tax_diver$propor,2)
tax_diver_genus <- tax_diver%>%filter(taxonomy=="genus")%>%mutate(group="term")%>%filter(!name%in%c("Unknown","Others"))
tax_diver_genus$propor <- paste0("-",tax_diver_genus$propor)
tax_diver_genus$propor <- as.numeric(tax_diver_genus$propor)
tax_diver_genus <- tax_diver_genus%>%arrange(desc(propor))
sum(tax_diver_genus$n)/34103#0.8724746 ARGs from the 10 genera in full-term infants
tax_diver_preterm <- read.csv("data/tax_diversity_preterm.csv")
tax_diver_preterm_genus <- tax_diver_preterm%>%filter(taxonomy=="genus")%>%mutate(group="preterm")%>%filter(!name%in%c("Unknown","Others"))
sum(tax_diver_preterm_genus$n)/156355 #0.9028813 ARGs from the 10 genera in preterm infants
tax_diver_cb <- rbind(tax_diver_preterm_genus,tax_diver_genus)
tax_diver_cb$name <- gsub("g__","",tax_diver_cb$name,fixed = T)
tax_diver_cb$name2 <- paste0(tax_diver_cb$group,"_",tax_diver_cb$name)
tax_diver_cb$name2 <- factor(tax_diver_cb$name2,levels = rev(tax_diver_cb$name2))
#Figure 2f
tax_diver_cb_plot <- ggplot(tax_diver_cb, aes(x=propor, y=name2, fill=group)) +
   geom_bar(stat = "identity")+
   geom_vline(xintercept = 0)+
   scale_fill_manual(values=c("#4994C4","#E09137"))+
   labs(x ="Proportion (%)", y="")+
   scale_x_continuous(expand = c(0,0),limits = c(-40,30),position = "top",
                      breaks = c(-40,-20,-10,0,10,20,30),labels = c("40","20","10","0","10","20","30"))+
   theme_bw()+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line.x = element_line(color="black"),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=12, colour = "black"),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=13, colour="black"),
         axis.ticks.y = element_blank(),
         axis.title.x = element_text(size=13, colour = "black"))

tax_diver_2 <- tax_diver%>%group_by(taxonomy)%>%dplyr::mutate(new_t=c(1:12))
tax_diver_2$new_t[tax_diver_2$new_t==1] <- 2
tax_diver_2$new_t[tax_diver_2$name=="Unknown"] <- 1
tax_diver_2$new_t <- as.character(tax_diver_2$new_t)
tax_diver_2$name <- gsub(".__","",tax_diver_2$name)
tax_diver_2$taxonomy <- paste(toupper(substr(tax_diver_2$taxonomy, 1, 1)), substr(tax_diver_2$taxonomy, 2, nchar(tax_diver_2$taxonomy)), sep="")
tax_diver_2$taxonomy <- factor(tax_diver_2$taxonomy,levels = c("Phylum", "Class","Order","Family","Genus"))
tax_diver_2$new_t <- factor(tax_diver_2$new_t, levels = c("1","12","11","10","9","8","7","6","5","4","3","2"))
#Supplementary Figure 6d
rgi_taxa_all <- ggplot(tax_diver_2, aes(x=taxonomy, y=n, fill=new_t)) + 
   geom_bar(stat = "identity", width=0.8)+
   scale_fill_manual(values = c("#CACACA","#8F8F8F","#6AABF5","#BEBADA","#FB8072","#80B1D3","#FDB462",
                                "#B3DE69","#FCCDE5","#CABB57","#BC80BD","#8DD3C7"))+
   theme_bw() +
   labs(x = element_blank(), y="Number of ARG ORFs")+
   theme(axis.title.y = element_text(size = 16))+
   theme(panel.border = element_blank(), axis.line = element_line()) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the grid
   scale_y_continuous(expand = c(0,0),breaks=c(0,10000,20000,30000,40000,50000,60000),
                      labels=c("0","10,000", "20,000", "30,000","40,000","50,000","60,000"))+
   theme(legend.position="none",
         aspect.ratio = 1.3/1,
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 45,vjust = 0.7))

#only species
rgi_term_profile_species <- rgi_term_profile%>%dplyr::select(sample,species,rpkm)
rgi_term_profile_species$species <- gsub("s__","",rgi_term_profile_species$species,fixed = T)
rgi_term_profile_species_sum <- rgi_term_profile_species%>%group_by(sample,species)%>%dplyr::summarise(total=sum(rpkm))
length(unique(rgi_term_profile_species_sum$sample))#1355
rgi_term_profile_species_sum_1 <- rgi_term_profile_species_sum%>%
   group_by(species)%>%
   dplyr::summarise("prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))
rgi_term_profile_species_sum_s <- spread(rgi_term_profile_species_sum,sample,total)
rgi_term_profile_species_sum_s[is.na(rgi_term_profile_species_sum_s)] <- 0
rgi_term_species_se <- data.frame("species"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in rgi_term_profile_species_sum_s$species) {
   a <- rgi_term_profile_species_sum_s%>%filter(species==t)%>%column_to_rownames("species")%>%t()%>%as.data.frame()%>%setNames("species")
   rgi_term_species_se[nrow(rgi_term_species_se)+1,] <- c(t,std.error(a$species))
}
rgi_term_profile_species_sum_2 <- left_join(rgi_term_profile_species_sum_1,rgi_term_species_se,by="species")
rgi_term_profile_species_sum_2$se <- as.numeric(rgi_term_profile_species_sum_2$se)
rgi_term_profile_species_sum_2 <- rgi_term_profile_species_sum_2%>%filter(species!="unknown")
rgi_term_profile_species_sum_3 <- rgi_term_profile_species_sum_2%>%dplyr::filter(prevalence>0.01)
#get the number of ARG ORFs for each species
rgi_term_profile_species_2 <- rgi_term_profile_species%>%group_by(species)%>%dplyr::summarise(n_args=n())%>%
   filter(species%in%rgi_term_profile_species_sum_3$species)%>%arrange(desc(n_args))
rgi_term_profile_species_sum_3 <- left_join(rgi_term_profile_species_sum_3,rgi_term_profile_species_2)
rgi_term_profile_species_sum_3$species <- factor(rgi_term_profile_species_sum_3$species,levels = rgi_term_profile_species_sum_3$species)
rgi_term_profile_species_sum_3$label <- paste0("(n=",formattable::comma(rgi_term_profile_species_sum_3$n_args,digits=0),")")
rgi_term_profile_species_sum_3$label[rgi_term_profile_species_sum_3$species=="Escherichia coli"]#"(n=19,725)"
rgi_term_profile_species_sum_3$label[rgi_term_profile_species_sum_3$species=="Escherichia coli"] <- NA
rgi_term_profile_species_sum_3$Mean_rpkm[rgi_term_profile_species_sum_3$species=="Escherichia coli"]#129159.2
rgi_term_profile_species_sum_3$Mean_rpkm[rgi_term_profile_species_sum_3$species=="Escherichia coli"] <- 33000
#Supplementary Figure 5e
rgi_taxa_species <- ggplot() +
   geom_bar(data=rgi_term_profile_species_sum_3, aes(x=species, y=Mean_rpkm),stat="identity",fill="#E09137",color="black",width = 0.8)+
   geom_errorbar(data=rgi_term_profile_species_sum_3,aes(x=species,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=rgi_term_profile_species_sum_3, aes(x=species, y=prevalence*80000,group=1),linetype = "dashed",color="#C43A0A",size=1)+
   geom_point(data=rgi_term_profile_species_sum_3, aes(x=species, y=prevalence*80000),size=1)+
   geom_text(data=rgi_term_profile_species_sum_3,aes(x=species, y=Mean_rpkm+se,label=label),angle=90,hjust=-0.1,vjust=0.5,size=3)+
   scale_y_continuous(expand = c(0,0),name= "Mean abundance (rpkm)",limits = c(0,40000),
                      breaks = c(0,10000,20000,30000,33000),labels = c("0","10,000","20,000","30,000","129,159"),
                      sec.axis = sec_axis(~ . /80000, breaks = c(0,0.1,0.2,0.3)))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(legend.position = "none",
         panel.background = element_blank(),
         plot.background = element_blank(),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=11, colour = "black",angle = 90,hjust = 1,vjust = 0.5,face="italic"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

#plasmid ARG ORFs cluster
rgi_term_profile_yes_orf <- rgi_term_profile%>%filter(bin!="nobin")%>%filter(plasmid=="Yes")
rgi_term_profile_yes_orf$sample_orf <- paste0(rgi_term_profile_yes_orf$sample,"_rgi_card_",rgi_term_profile_yes_orf$ORF_ID)
rgi_cluster_plas <- read.delim("data/rgi_term_plasmid_nr_cluster.txt")
length(unique(rgi_cluster_plas$clstr))#204 clusters formed
length(unique(rgi_cluster_plas$id))#423 args
rgi_term_profile_yes_orf_cluster <- left_join(rgi_term_profile_yes_orf,rgi_cluster_plas%>%dplyr::select(id,clstr),by=c("sample_orf"="id"))

for (v in c("phylum","class","order","family","genus","species")) {
   print(v)
   cluster_sum_1 <- rgi_term_profile_yes_orf_cluster%>%group_by(clstr)%>%dplyr::summarise(n_taxa_all=length(unique(get(v))),n_args=n())
   cluster_sum_2 <- rgi_term_profile_yes_orf_cluster%>%filter(get(v)!="unknown")%>%group_by(clstr)%>%dplyr::summarise(n_taxa_sel=length(unique(get(v))),n_args_2=n())
   print(sum(cluster_sum_1$n_args))
   print(sum(cluster_sum_2$n_args_2))
   cluster_sum <- left_join(cluster_sum_1,cluster_sum_2,by="clstr")
   cluster_sum$n_taxa_sel[is.na(cluster_sum$n_taxa_sel)] <- 0
   cluster_sum$n_args_2[is.na(cluster_sum$n_args_2)] <- 0
   cluster_sum <- cluster_sum%>%mutate(status=case_when(n_taxa_sel>1 ~ "multiple_taxa",
                                                        n_taxa_sel==1 ~ "single_taxa",
                                                        n_taxa_sel==0 ~ "unassigned"))%>%
      filter(n_taxa_all>1)
   print(nrow(cluster_sum))
   cluster_sum_sel <- cluster_sum%>%filter(n_taxa_sel>1)
   print(nrow(cluster_sum_sel))
}   
#32 clusters with at least 2 species including unknown species; 32 clusters with at least 2 known species
cluster_sum_sel_species <- rgi_term_profile_yes_orf_cluster%>%filter(clstr%in%cluster_sum_sel$clstr)
length(unique(cluster_sum_sel_species$species))-1#34
cluster_sum_sel_species_cl <- cluster_sum_sel_species%>%group_by(clstr)%>%dplyr::summarise(n_species=length(unique(species)))%>%arrange(desc(n_species))
cluster_sum_sel_species_cl_orf <- cluster_sum_sel_species%>%group_by(clstr,Best_Hit_ARO)%>%dplyr::summarise(n_ARGs=n())
cluster_sum_sel_species_cl_cb <- left_join(cluster_sum_sel_species_cl,cluster_sum_sel_species_cl_orf,by="clstr")

cluster_sum_sel_species_sum <- cluster_sum_sel_species%>%filter(species!="unknown")%>%group_by(clstr,species)%>%dplyr::summarise(n_args=n())
cluster_sum_sel_species_sum_1 <- cluster_sum_sel_species_sum%>%group_by(species)%>%dplyr::summarise(n=n())%>%arrange(desc(n))
cluster_sum_sel_species_sum_2 <- cluster_sum_sel_species_sum%>%group_by(clstr)%>%dplyr::summarise(n=n())%>%filter(n>=4)
cluster_sum_sel_species_sum_sel <- cluster_sum_sel_species_sum%>%
   filter(clstr%in%cluster_sum_sel_species_sum_2$clstr)%>%
   mutate(n_args_2=case_when(n_args>20 ~ "10",
                             TRUE ~ as.character(n_args)))%>%
   dplyr::select(-n_args)
cluster_sum_sel_species_sum_sel$species <- gsub("s__","",cluster_sum_sel_species_sum_sel$species,fixed = T)
cluster_sum_sel_species_sum_sel$n_args_2 <- as.numeric(cluster_sum_sel_species_sum_sel$n_args_2)
quantile(cluster_sum_sel_species_sum_sel$n_args_2)
cluster_sum_sel_species_s <- spread(cluster_sum_sel_species_sum_sel,species,n_args_2)%>%column_to_rownames("clstr")
cluster_sum_sel_species_s[is.na(cluster_sum_sel_species_s)] <- 0
names(cluster_sum_sel_species_s)
col_fun_type <- colorRamp2(c(0,1,2,4,6,8,10),c("white","#A8C6EC","#7FA4D3","#5A87BF","#3768A4","#184B8A","#073876"))
#Figure 2g
plasmid_cluster <- Heatmap(as.matrix(cluster_sum_sel_species_s),cluster_rows=T,cluster_columns=T,col = col_fun_type,
                           show_row_names = T,show_column_names = T,
                           border = T,row_names_gp = gpar(fontsize =10),
                           rect_gp = gpar(col = "grey", lwd = 0.7),
                           column_names_gp = gpar(fontface="italic",fontsize=12),
                           border_gp=gpar(col = "black", lwd = 1),column_gap = unit(2, "mm"),
                           width = unit(12, "cm"), height = unit(4, "cm"),
                           show_heatmap_legend = TRUE,
                           heatmap_legend_param = list(title = "Number of ARG ORFs",direction = "vertical", 
                                                       legend_height = unit(3, "cm"),legend_width = unit(0.3, "cm"),
                                                       at = c(0,2,4,6,8,10), labels = c("0","2","4","6","8",">10"))
)
cluster_shared <- list()
length(unique(cluster_sum_sel_species_sum$species))#34
for (s in unique(cluster_sum_sel_species_sum$species)) {
   a <- cluster_sum_sel_species_sum%>%filter(species==s)
   b <- cluster_sum_sel_species_sum%>%filter(species!=s)%>%filter(clstr%in%a$clstr)
   b_sum <- b%>%group_by(species)%>%dplyr::summarise(n=n())%>%arrange(desc(n))
   cluster_shared[[s]] <- b_sum
}

####################################################################################################################################
#Differential carriage of ARG potentials from bacterial species in preterm infants
####################################################################################################################################
term_metadata <- read.csv("data/metadata_cb_infant_term.csv")
preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")
metadata_cb <- rbind(preterm_metadata,term_metadata)

rgi_preterm_profile <- read.csv("data/rgi_preterm_resistome_coverm_bins_gtdb_stats_plasmid.csv")
rgi_preterm_profile$species <- sapply(strsplit(rgi_preterm_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[7]))
rgi_preterm_profile$species[is.na(rgi_preterm_profile$species)] <- "unknown"
rgi_preterm_profile$species[rgi_preterm_profile$species=="s__"] <- "unknown"
rgi_preterm_profile$species <- gsub("s__","",rgi_preterm_profile$species,fixed = T)
rgi_preterm_profile$ncbi_species <- sapply(strsplit(rgi_preterm_profile$NCBI_classification, split=';', fixed=TRUE), function(x)(x[7]))
rgi_preterm_profile$ncbi_species[is.na(rgi_preterm_profile$ncbi_species)] <- "unknown"
rgi_preterm_profile$ncbi_species[rgi_preterm_profile$ncbi_species=="s__"] <- "unknown"
rgi_preterm_profile$ncbi_species <- gsub("s__","",rgi_preterm_profile$ncbi_species,fixed = T)
rgi_preterm_profile_species <- rgi_preterm_profile%>%dplyr::select(sample,species,rpkm)
rgi_preterm_profile_species_sum <- rgi_preterm_profile_species%>%group_by(sample,species)%>%dplyr::summarise(total=sum(rpkm))
length(unique(rgi_preterm_profile_species_sum$species))#601 species with unknown
rgi_preterm_profile_species_sum_1 <- rgi_preterm_profile_species_sum%>%
   group_by(species)%>%
   dplyr::summarise("prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))%>%filter(species!="unknown")%>%
   filter(prevalence>0.01)

rgi_term_profile <- read.csv("data/rgi_term_resistome_coverm_bins_gtdb_stats_plasmid.csv")
rgi_term_profile$species <- sapply(strsplit(rgi_term_profile$GTDB_classification, split=';', fixed=TRUE), function(x)(x[7]))
rgi_term_profile$species[is.na(rgi_term_profile$species)] <- "unknown"
rgi_term_profile$species[rgi_term_profile$species=="s__"] <- "unknown"
rgi_term_profile$species <- gsub("s__","",rgi_term_profile$species,fixed = T)
rgi_term_profile$ncbi_species <- sapply(strsplit(rgi_term_profile$NCBI_classification, split=';', fixed=TRUE), function(x)(x[7]))
rgi_term_profile$ncbi_species[is.na(rgi_term_profile$ncbi_species)] <- "unknown"
rgi_term_profile$ncbi_species[rgi_term_profile$ncbi_species=="s__"] <- "unknown"
rgi_term_profile$ncbi_species <- gsub("s__","",rgi_term_profile$ncbi_species,fixed = T)
rgi_term_profile_species <- rgi_term_profile%>%dplyr::select(sample,species,rpkm)
rgi_term_profile_species_sum <- rgi_term_profile_species%>%group_by(sample,species)%>%dplyr::summarise(total=sum(rpkm))
length(unique(rgi_term_profile_species_sum$species))#472 species with unknown
rgi_term_profile_species_sum_1 <- rgi_term_profile_species_sum%>%
   group_by(species)%>%
   dplyr::summarise("prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))%>%filter(species!="unknown")%>%
   filter(prevalence>0.01)
shared_list <- data.frame("species"=intersect(rgi_term_profile_species_sum_1$species,rgi_preterm_profile_species_sum_1$species))
shared_list$genus <- sapply(strsplit(shared_list$species, split=' ', fixed=TRUE), function(x)(x[1]))
shared_list%>%group_by(genus)%>%dplyr::summarise(n=n())%>%arrange(desc(n))
#rgi class
rgi_preterm_profile_drug <- rgi_preterm_profile %>%
   dplyr::select(sample,species,cb,rpkm,Drug.Class)%>%
   dplyr::bind_cols(split_into_multiple(.$Drug.Class, "; ", "cat"))%>%
   gather(.,cat,drug,-c(sample,species,cb,rpkm,Drug.Class))%>%
   filter(!is.na(drug))%>%
   dplyr::select(-cat)%>%filter(species!="unknown")

rgi_term_profile_drug <- rgi_term_profile %>%
   dplyr::select(sample,species,cb,rpkm,Drug.Class)%>%
   dplyr::bind_cols(split_into_multiple(.$Drug.Class, "; ", "cat"))%>%
   gather(.,cat,drug,-c(sample,species,cb,rpkm,Drug.Class))%>%
   filter(!is.na(drug))%>%
   dplyr::select(-cat)%>%filter(species!="unknown")

preterm_term_drug_pval <- data.frame("species"="x","preterm_mean"=0,"term_mean"=0,"pval"=0,stringsAsFactors = F)[-1,]
count=0
for (s in shared_list$species) {
   count=count+1
   print(count)
   a_preterm_1 <- rgi_preterm_profile_drug%>%filter(species==s)%>%group_by(sample)%>%
      dplyr::summarise(total=sum(rpkm))%>%mutate(group="preterm")
   a_preterm_2 <- data.frame("sample"=unique(rgi_preterm_profile$sample),"total"=0,"group"="preterm")%>%filter(!sample%in%a_preterm_1$sample)
   a_preterm <- rbind(a_preterm_1,a_preterm_2)
   a_term_1 <- rgi_term_profile_drug%>%filter(species==s)%>%group_by(sample)%>%
      dplyr::summarise(total=sum(rpkm))%>%mutate(group="term")
   a_term_2 <- data.frame("sample"=unique(rgi_term_profile$sample),"total"=0,"group"="term")%>%filter(!sample%in%a_term_1$sample)
   a_term <- rbind(a_term_1,a_term_2)
   a <- rbind(a_preterm,a_term)
   a <- left_join(a,metadata_cb%>%dplyr::select(Study,SequenceID),by=c("sample"="SequenceID"))
   a_sum <- a%>%group_by(group)%>%dplyr::summarise(mean=mean(total))
   a_fit <- lmerTest::lmer(total ~ group+(1|Study),a)
   a_pval <- anova(a_fit)
   preterm_term_drug_pval[nrow(preterm_term_drug_pval)+1,] <- c(s,a_sum$mean[1],a_sum$mean[2],a_pval$`Pr(>F)`)
}
preterm_term_drug_pval <- preterm_term_drug_pval%>%mutate_at(vars(-species),as.numeric)
preterm_term_drug_pval$fdr <- p.adjust(preterm_term_drug_pval$pval,method = "fdr")
preterm_term_drug_pval$diff <- preterm_term_drug_pval$preterm_mean-preterm_term_drug_pval$term_mean
preterm_term_drug_pval_sig <- preterm_term_drug_pval%>%filter(pval<0.05)%>%arrange(desc(diff))

preterm_term_drug_species_pval <- data.frame("species"="x","drug"="x","preterm_total"=0,"term_total"=0,"pval"=0,stringsAsFactors = F)[-1,]
count=0
for (s in shared_list$species) {
   count=count+1
   print(count)
   a_preterm_1 <- rgi_preterm_profile_drug%>%filter(species==s)%>%group_by(sample,drug)%>%
      dplyr::summarise(total=sum(rpkm))%>%mutate(group="preterm")
   a_term_1 <- rgi_term_profile_drug%>%filter(species==s)%>%group_by(sample,drug)%>%
      dplyr::summarise(total=sum(rpkm))%>%mutate(group="term")
   a_1 <- rbind(a_preterm_1,a_term_1)
   for (d in unique(a_1$drug)) {
      # d="monobactam"
      b <- a_1%>%filter(drug==d)
      b_sum <- b%>%group_by(group)%>%dplyr::summarise(total=sum(total))
      a_preterm_2 <- data.frame("sample"=unique(rgi_preterm_profile$sample),"drug"=NA,"total"=0,"group"="preterm")%>%filter(!sample%in%b$sample)
      a_term_2 <- data.frame("sample"=unique(rgi_term_profile$sample),"drug"=NA,"total"=0,"group"="term")%>%filter(!sample%in%b$sample)
      a_2 <- rbind(a_preterm_2,a_term_2)
      a <- rbind(b,a_2)
      a <- left_join(a,metadata_cb%>%dplyr::select(Study,SequenceID),by=c("sample"="SequenceID"))
      a_fit <- lmerTest::lmer(total ~ group+(1|Study),a)
      a_pval <- anova(a_fit)  
      if (length(unique(b$group))==2) {
         preterm_term_drug_species_pval[nrow(preterm_term_drug_species_pval)+1,] <- c(s,d,b_sum$total[b_sum$group=="preterm"],b_sum$total[b_sum$group=="term"],a_pval)
      }else if (unique(b$group)=="term") {
         preterm_term_drug_species_pval[nrow(preterm_term_drug_species_pval)+1,] <- c(s,d,0,b_sum$total[b_sum$group=="term"],a_pval)
      }else if (unique(b$group)=="preterm") {
         preterm_term_drug_species_pval[nrow(preterm_term_drug_species_pval)+1,] <- c(s,d,b_sum$total[b_sum$group=="preterm"],0,a_pval)
      }
   }
}
preterm_term_drug_species_pval <- preterm_term_drug_species_pval%>%mutate_at(vars(-species,-drug),as.numeric)
preterm_term_drug_species_pval$fdr <- p.adjust(preterm_term_drug_species_pval$pval,method = "fdr")
preterm_term_drug_species_pval$preterm_mean <- preterm_term_drug_species_pval$preterm_total/3300
preterm_term_drug_species_pval$term_mean <- preterm_term_drug_species_pval$term_total/1355
preterm_term_drug_species_pval$diff <- preterm_term_drug_species_pval$preterm_mean-preterm_term_drug_species_pval$term_mean
preterm_term_drug_species_pval_sig <- preterm_term_drug_species_pval%>%filter(pval<0.05)%>%arrange(desc(diff))
preterm_term_drug_species_pval_sig_sum <- preterm_term_drug_species_pval_sig%>%mutate(diff2=sign(diff))%>%
   group_by(drug,diff2)%>%dplyr::summarise(n=sum(diff2))
length(unique(preterm_term_drug_species_pval_sig$species))
length(unique(preterm_term_drug_species_pval_sig$drug))

#to make plot
preterm_drug <- rgi_preterm_profile_drug%>%filter(species%in%shared_list$species)%>%group_by(species,drug)%>%
   dplyr::summarise(total=sum(rpkm))%>%mutate(group=paste0(species,"__preterm"))

term_drug <- rgi_term_profile_drug%>%filter(species%in%shared_list$species)%>%group_by(species,drug)%>%
   dplyr::summarise(total=sum(rpkm))%>%mutate(group=paste0(species,"__term"))

preterm_term_drug <- rbind(preterm_drug%>%ungroup()%>%dplyr::select(-species),term_drug%>%ungroup()%>%dplyr::select(-species))
preterm_term_drug$group2 <- sapply(strsplit(preterm_term_drug$group, split='__', fixed=TRUE), function(x)(x[2]))
preterm_term_drug$group1 <- sapply(strsplit(preterm_term_drug$group, split='__', fixed=TRUE), function(x)(x[1]))

preterm_term_drug_sum <- preterm_term_drug%>%group_by(group2,group1)%>%
   dplyr::summarise(n=sum(total>0))%>%spread(.,group2,n)%>%mutate(diff=preterm-term)

length(unique(preterm_term_drug$drug))#36
preterm_term_drug_spread <- spread(preterm_term_drug%>%dplyr::select(-group1,-group2),drug,total)
names(preterm_term_drug_spread) <- gsub(" antibiotic","",names(preterm_term_drug_spread),fixed = T)
names(preterm_term_drug_spread)[4]
names(preterm_term_drug_spread)[4] <- "AFFA"
names(preterm_term_drug_spread)[10]
names(preterm_term_drug_spread)[10] <- "DAA"
preterm_term_drug_spread <- preterm_term_drug_spread%>%column_to_rownames("group")
preterm_term_drug_spread[is.na(preterm_term_drug_spread)] <- 0

hp_sort <- data.frame("total"=rowSums(preterm_term_drug_spread))%>%rownames_to_column("species_group")
hp_sort$group <- sapply(strsplit(hp_sort$species_group, split='__', fixed=TRUE), function(x)(x[2]))
hp_sort$species <- sapply(strsplit(hp_sort$species_group, split='__', fixed=TRUE), function(x)(x[1]))
hp_sort_preterm <- hp_sort%>%filter(group=="preterm")%>%arrange(desc(total))
hp_sort <- hp_sort%>%arrange(match(species,hp_sort_preterm$species))
preterm_term_drug_spread_2 <- preterm_term_drug_spread%>%rownames_to_column("species_group")%>%
   arrange(match(species_group,hp_sort$species_group))%>%column_to_rownames("species_group")
preterm_term_drug_spread_2[seq(1,56,2),] <- preterm_term_drug_spread_2[seq(1,56,2),]/3300
preterm_term_drug_spread_2[seq(2,56,2),] <- preterm_term_drug_spread_2[seq(2,56,2),]/1355
preterm_term_drug_spread_3 <- preterm_term_drug_spread_2+0.00001
range(preterm_term_drug_spread_3)#0.00001 71427.20458
preterm_term_drug_spread_3_log <- log10(preterm_term_drug_spread_3)
range(preterm_term_drug_spread_3_log)#-5.000000  4.853864

#preterm_term_species_pval is from below code
names(preterm_term_species_pval)[1] <- "mph_species"
preterm_term_species_pval_ncbi <- left_join(preterm_term_species_pval,shared_list_ncbi,by=c("mph_species"="ncbi_species_2"))
preterm_term_species_pval_preterm <- preterm_term_species_pval_ncbi%>%dplyr::select(species,preterm_mean)%>%
   mutate(species_group=paste0(species,"__preterm"))%>%dplyr::select(-species)
names(preterm_term_species_pval_preterm)[1] <- "mean"
preterm_term_species_pval_term <- preterm_term_species_pval_ncbi%>%dplyr::select(species,term_mean)%>%
   mutate(species_group=paste0(species,"__term"))%>%dplyr::select(-species)
names(preterm_term_species_pval_term)[1] <- "mean"
preterm_term_species_pval_term_preterm <- rbind(preterm_term_species_pval_term,preterm_term_species_pval_preterm)
preterm_term_species_pval_term_preterm[nrow(preterm_term_species_pval_term_preterm)+1,] <- c(NA,"Klebsiella variicola__preterm")
preterm_term_species_pval_term_preterm[nrow(preterm_term_species_pval_term_preterm)+1,] <- c(NA,"Klebsiella variicola__term")
preterm_term_species_pval_term_preterm[nrow(preterm_term_species_pval_term_preterm)+1,] <- c(NA,"Streptococcus sp000187445__preterm")
preterm_term_species_pval_term_preterm[nrow(preterm_term_species_pval_term_preterm)+1,] <- c(NA,"Streptococcus sp000187445__term")
preterm_term_species_pval_term_preterm[nrow(preterm_term_species_pval_term_preterm)+1,] <- c(NA,"Veillonella nakazawae__preterm")
preterm_term_species_pval_term_preterm[nrow(preterm_term_species_pval_term_preterm)+1,] <- c(NA,"Veillonella nakazawae__term")
length(unique(preterm_term_species_pval_term_preterm$species_group))#56
preterm_term_species_pval_term_preterm <- preterm_term_species_pval_term_preterm%>%
   arrange(match(species_group,rownames(preterm_term_drug_spread_2)))
preterm_term_species_pval_term_preterm$mean <- as.numeric(preterm_term_species_pval_term_preterm$mean)

term_preterm_species_drug_plus_1 <- preterm_term_drug_species_pval_sig%>%dplyr::select(species,drug,diff)%>%spread(.,drug,diff)
term_preterm_species_drug_plus_2 <- preterm_term_drug%>%filter(!group1%in%term_preterm_species_drug_plus_1$species)%>%distinct(group1)
names(term_preterm_species_drug_plus_2) <- "species"
term_preterm_species_drug_plus <- rbind.fill(term_preterm_species_drug_plus_1,term_preterm_species_drug_plus_2)
term_preterm_species_drug_plus_3 <- preterm_term_drug%>%filter(!drug%in%names(term_preterm_species_drug_plus_1)[-1])%>%distinct(drug)
term_preterm_species_drug_plus_4 <- data.frame(rep(0,28),rep(0,28),rep(0,28),rep(0,28),rep(0,28))#based on the number row of term_preterm_species_drug_plus_3
names(term_preterm_species_drug_plus_4) <- term_preterm_species_drug_plus_3$drug
term_preterm_species_drug_plus_4$species <- term_preterm_species_drug_plus$species
term_preterm_species_drug_plus_preterm <- left_join(term_preterm_species_drug_plus,term_preterm_species_drug_plus_4,by="species")
term_preterm_species_drug_plus_term <- term_preterm_species_drug_plus_preterm
term_preterm_species_drug_plus_term[,-1] <- 0
term_preterm_species_drug_plus_preterm$species <- paste0(term_preterm_species_drug_plus_preterm$species,"__preterm")
term_preterm_species_drug_plus_term$species <- paste0(term_preterm_species_drug_plus_term$species,"__term")
term_preterm_species_drug_plus_term_preterm <- rbind(term_preterm_species_drug_plus_term,term_preterm_species_drug_plus_preterm)
names(term_preterm_species_drug_plus_term_preterm) <- gsub(" antibiotic","",names(term_preterm_species_drug_plus_term_preterm))
names(term_preterm_species_drug_plus_term_preterm)[8]
names(term_preterm_species_drug_plus_term_preterm)[8] <- "DAA"
names(term_preterm_species_drug_plus_term_preterm)[35]
names(term_preterm_species_drug_plus_term_preterm)[35] <- "AFFA"
term_preterm_species_drug_plus_term_preterm[is.na(term_preterm_species_drug_plus_term_preterm)] <- 0
term_preterm_species_drug_plus_term_preterm_order <- term_preterm_species_drug_plus_term_preterm%>%
   dplyr::select(species,names(preterm_term_drug_spread_3_log))%>%
   arrange(match(species,rownames(preterm_term_drug_spread_3_log)))%>%column_to_rownames("species")

row_annotation_right <- 
   rowAnnotation(
      `Types of drug\nclasses` = anno_barplot(rowSums(preterm_term_drug_spread_2>0), border=TRUE, 
                                              gp=gpar(fill='#FEE8C8', col='#FEE8C8'),
                                              axis=TRUE),
      `Mean abundance\nof drug classes\n(log10,rpkm)` = anno_barplot(log(rowSums(preterm_term_drug_spread_2)), border=TRUE, 
                                                                     gp=gpar(fill='#FDBB84', col='#FDBB84'),
                                                                     axis=TRUE),
      `Mean abundance\nof species(%)` = anno_barplot(preterm_term_species_pval_term_preterm$mean,
                                                     border=TRUE, gp=gpar(fill='#FC8D59', col='#FC8D59'),axis=TRUE),
      width=unit(3.8, 'inches'),
      gap = unit(5, "mm"),
      annotation_name_gp=gpar(fontsize = 12))
row_annotation_left <- 
   rowAnnotation(
      `Gestational age` = rep(c("Preterm","Full term"),28),
      col=list(`Gestational age`=c("Preterm"="#4994C4","Full term"="#E09137")),width=unit(0.1, 'inches'),
      annotation_name_gp=gpar(fontsize = 14))
ann <- data.frame("species"=rownames(preterm_term_drug_spread_3_log))
ann$species <- gsub("__term","",ann$species)
ann$species <- gsub("__preterm","",ann$species)
ann$species <- factor(ann$species,levels = unique(ann$species))
ann <- ann%>%group_by(species)%>%dplyr::mutate(species_2=cur_group_id())
col_fun_type <- colorRamp2(c(-5,-2.5,-1,0,1,2,3,4,5),c("#4A7FBF","#83A9D6","#AFC9E8","white",
                                                       "#FFF7EC", "#FDD49E", "#FC8D59", "#EF6548", "#990000"))
#Figure 3
term_preterm_species_drug <- Heatmap(as.matrix(preterm_term_drug_spread_3_log),cluster_rows=F,cluster_columns=T,col = col_fun_type,
                                     show_row_names = T,show_column_names = T,
                                     border = T,column_names_gp = gpar(fontsize =14),
                                     right_annotatio = row_annotation_right,
                                     left_annotation = row_annotation_left,
                                     row_split =ann$species_2, 
                                     row_names_side = "left",
                                     rect_gp = gpar(col = "grey", lwd = 0.7),
                                     # column_names_gp = gpar(fontface="italic"),
                                     border_gp=gpar(col = "black", lwd = 1),column_gap = unit(2, "mm"),
                                     width = unit(26, "cm"), height = unit(22, "cm"),
                                     show_heatmap_legend = TRUE,
                                     heatmap_legend_param = list(title = "Log10(rpkm)",direction = "horizontal", at = c(-5,-2.5, 0, 2.5, 5),
                                                                 legend_hight = unit(1, "cm"),legend_width = unit(3, "cm")),
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                        if(term_preterm_species_drug_plus_term_preterm_order[i,j]>0) {
                                           grid.text("+", x = x, y = y,gp = gpar(fontsize = 10))
                                        } else if(term_preterm_species_drug_plus_term_preterm_order[i,j]<0) {
                                           grid.text("-", x = x, y = y,gp = gpar(fontsize = 10))
                                        } else {
                                           grid.text("", x = x, y = y,gp = gpar(fontsize = 10))
                                        }
                                     }
)

#metaphlan4 and rgi
metaphlan4_preterm <- read.csv("data/metaphlan4_preterm.csv")
metaphlan4_term <- read.csv("data/metaphlan4_term.csv")
metaphlan4_taxa <- left_join(metaphlan4_preterm,metaphlan4_term,by="clade_name")
metaphlan4_taxa[is.na(metaphlan4_taxa)] <- 0
metaphlan4_taxa_2 <- metaphlan4_taxa%>%dplyr::select(clade_name,unique(rgi_preterm_profile$sample),unique(rgi_term_profile$sample))
metaphlan4_species <- metaphlan4_taxa_2%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
metaphlan4_species$species <- sapply(strsplit(metaphlan4_species$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
metaphlan4_species <- metaphlan4_species%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
#compare ncbi and gtdb species: so the species from gtdb match uniquely to species from ncbi
rgi_preterm_profile_term <- rbind(rgi_preterm_profile,rgi_term_profile)
cp_ncbi_gtdb <- rgi_preterm_profile_term%>%group_by(species)%>%dplyr::summarise(n=length(unique(ncbi_species)))
cp_ncbi_gtdb_species <- rgi_preterm_profile_term%>%distinct(species,.keep_all = T)%>%dplyr::select(species,ncbi_species)

preterm_drug_2 <- rgi_preterm_profile_drug%>%filter(species%in%shared_list$species)%>%group_by(sample,species)%>%dplyr::summarise(drug_total=sum(rpkm))
term_drug_2 <- rgi_term_profile_drug%>%filter(species%in%shared_list$species)%>%group_by(sample,species)%>%dplyr::summarise(drug_total=sum(rpkm))
names(metaphlan4_species) <- gsub("_"," ",names(metaphlan4_species))
shared_list_ncbi <- left_join(shared_list,cp_ncbi_gtdb_species,by="species")
shared_list_ncbi$ncbi_species_2 <- shared_list_ncbi$ncbi_species
shared_list_ncbi$ncbi_species_2 <- gsub("[","",shared_list_ncbi$ncbi_species_2,fixed = T)
shared_list_ncbi$ncbi_species_2 <- gsub("]","",shared_list_ncbi$ncbi_species_2,fixed = T)
which(!shared_list_ncbi$species%in%names(metaphlan4_species))#7 11 15 16 18 21 22 25 26 28
which(!shared_list_ncbi$ncbi_species%in%names(metaphlan4_species))#7 21 22 25
which(!shared_list_ncbi$ncbi_species_2%in%names(metaphlan4_species))#21 22 25
#keep the species that are present in gtdb and metaphlan4
metaphlan4_shared_list <- metaphlan4_species%>%dplyr::select(shared_list_ncbi$ncbi_species_2[-c(21,22,25)])

metaphlan4_shared_list_term <- metaphlan4_shared_list%>%filter(rownames(.)%in%term_drug_2$sample)%>%rownames_to_column("sample")
metaphlan4_shared_list_term_g <- gather(metaphlan4_shared_list_term,species,species_abund,-sample)
term_drug_2_ncbi <- left_join(term_drug_2,shared_list_ncbi,by="species")
term_drug_2_species <- left_join(term_drug_2_ncbi,metaphlan4_shared_list_term_g,by=c("sample"="sample","ncbi_species_2"="species"))
term_drug_2_species%>%filter(is.na(species_abund))%>%group_by(ncbi_species_2)%>%dplyr::summarise(n=n())
term_drug_2_species_sum <- term_drug_2_species%>%ungroup()%>%dplyr::select(-sample)%>%group_by(species)%>%
   dplyr::summarise(total_drug=sum(drug_total),total_abund=sum(species_abund))
term_drug_2_species_sum$mean_drug <- term_drug_2_species_sum$total_drug/1355#1355 samples for rgi
term_drug_2_species_sum$mean_abund <- term_drug_2_species_sum$total_abund/1355
term_drug_2_species_sum$group <- "Full term"
cor(term_drug_2_species_sum$mean_drug[-c(17,24,26)],term_drug_2_species_sum$mean_abund[-c(17,24,26)])#0.5439327
trendline_sum(term_drug_2_species_sum$mean_drug[-c(17,24,26)], term_drug_2_species_sum$mean_abund[-c(17,24,26)], model="line2P",Pvalue.corrected = FALSE,summary = TRUE,eDigit = 5)
# F-statistic: 9.6641 on 1 and 23 DF,  p-value: 0.0049444
metaphlan4_shared_list_preterm <- metaphlan4_shared_list%>%filter(rownames(.)%in%preterm_drug_2$sample)%>%rownames_to_column("sample")
metaphlan4_shared_list_preterm_g <- gather(metaphlan4_shared_list_preterm,species,species_abund,-sample)
preterm_drug_2_ncbi <- left_join(preterm_drug_2,shared_list_ncbi,by="species")
preterm_drug_2_species <- left_join(preterm_drug_2_ncbi,metaphlan4_shared_list_preterm_g,by=c("sample"="sample","ncbi_species_2"="species"))
preterm_drug_2_species_sum <- preterm_drug_2_species%>%ungroup()%>%dplyr::select(-sample)%>%group_by(species)%>%
   dplyr::summarise(total_drug=sum(drug_total),total_abund=sum(species_abund))
preterm_drug_2_species_sum$mean_drug <- preterm_drug_2_species_sum$total_drug/3300#3301 samples for rgi
preterm_drug_2_species_sum$mean_abund <- preterm_drug_2_species_sum$total_abund/3300
preterm_drug_2_species_sum$group <- "Preterm"
cor(preterm_drug_2_species_sum$mean_drug[-c(17,24,26)],preterm_drug_2_species_sum$mean_abund[-c(17,24,26)])#0.8874334
trendline_sum(preterm_drug_2_species_sum$mean_drug[-c(17,24,26)], preterm_drug_2_species_sum$mean_abund[-c(17,24,26)], model="line2P",Pvalue.corrected = FALSE,summary = TRUE,eDigit = 5)
preterm_drug_2_species_sum_term <- rbind(term_drug_2_species_sum,preterm_drug_2_species_sum)
#Supplementary Figure 7b
preterm_term_rgi_species <- ggplot(preterm_drug_2_species_sum_term%>%filter(!is.na(total_abund)), aes(x=mean_drug, y=mean_abund,color=group)) + 
   geom_point(alpha=1)+
   geom_smooth(method = lm,size=1,level=0.95)+
   scale_color_manual(values = c("#E09137","#4994C4"))+
   guides(color = guide_legend(reverse=T))+
   labs(x ="Mean abundance of ARG drug classes (rpkm)", y="Mean relative abundance of species (%)")+
   scale_x_continuous(labels=c("0","100,000","200,000","300,000","400,000","500,000"))+
   theme_bw()+
   theme(plot.margin = unit(c(1,5,1,1), "mm"))+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour="black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour="black"),
         axis.title.x = element_text(size=12, colour = "black"))

#procrustes for species and arg families
preterm_taxa_sel_2 <- metaphlan4_taxa%>%dplyr::select(clade_name,unique(rgi_preterm_profile$sample))
preterm_taxa_species_2 <- preterm_taxa_sel_2%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
preterm_taxa_species_2$species <- sapply(strsplit(preterm_taxa_species_2$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
preterm_taxa_species_2 <- preterm_taxa_species_2%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))

rgi_preterm_gene <- rgi_preterm_profile %>%
   dplyr::select(sample,species,cb,rpkm,AMR.Gene.Family)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(sample,species,cb,rpkm,AMR.Gene.Family))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)
rgi_preterm_gene_sum <- rgi_preterm_gene%>%group_by(sample,gene)%>%dplyr::summarise(rpkm=sum(rpkm))%>%spread(.,gene,rpkm)
rgi_preterm_gene_sum[is.na(rgi_preterm_gene_sum)] <- 0
rgi_preterm_gene_sum_sel <- rgi_preterm_gene_sum%>%dplyr::select(sample,which(colSums(.[-1])!=0)+1)
rgi_preterm_gene_sum_sel <- rgi_preterm_gene_sum_sel[which(rowSums(rgi_preterm_gene_sum_sel[-1])!=0),]

preterm_taxa_species_3 <- preterm_taxa_species_2%>%filter(rownames(.)%in%rgi_preterm_gene_sum_sel$sample)
rgi_preterm_gene_sum_1 <- rgi_preterm_gene_sum_sel%>%arrange(match(sample,rownames(preterm_taxa_species_3)))%>%column_to_rownames("sample")
identical(rownames(rgi_preterm_gene_sum_1),rownames(preterm_taxa_species_3))#TRUE

preterm_taxa_species_2_h <- decostand(preterm_taxa_species_3,method = "hellinger")
preterm_taxa_species_2_h_beta <- vegdist(preterm_taxa_species_2_h, method = "bray")
preterm_taxa_species_2_h_beta_pcoa <- cmdscale(preterm_taxa_species_2_h_beta, k = 2, eig = TRUE)
preterm_gene_h <- decostand(rgi_preterm_gene_sum_1,method = "hellinger")
preterm_gene_h_beta <- vegdist(preterm_gene_h, method = "bray")
preterm_gene_h_beta_pcoa <- cmdscale(preterm_gene_h_beta, k = 2, eig = TRUE)
preterm_taxa_gene_proc <- procrustes(X = preterm_taxa_species_2_h_beta_pcoa, Y = preterm_gene_h_beta_pcoa, symmetric = TRUE)
set.seed(1000)
preterm_taxa_gene_proc_p <- protest(X = preterm_taxa_species_2_h_beta_pcoa, Y = preterm_gene_h_beta_pcoa, permutations=9999)
preterm_taxa_gene_proc_p
term_taxa_sel_2 <- metaphlan4_taxa%>%dplyr::select(clade_name,unique(rgi_term_profile$sample))
term_taxa_species_2 <- term_taxa_sel_2%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
term_taxa_species_2$species <- sapply(strsplit(term_taxa_species_2$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
term_taxa_species_2 <- term_taxa_species_2%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
rgi_term_gene <- rgi_term_profile %>%
   dplyr::select(sample,species,cb,rpkm,AMR.Gene.Family)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(sample,species,cb,rpkm,AMR.Gene.Family))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)
rgi_term_gene_sum <- rgi_term_gene%>%group_by(sample,gene)%>%dplyr::summarise(rpkm=sum(rpkm))%>%spread(.,gene,rpkm)
rgi_term_gene_sum[is.na(rgi_term_gene_sum)] <- 0
rgi_term_gene_sum_sel <- rgi_term_gene_sum%>%dplyr::select(sample,which(colSums(.[-1])!=0)+1)
term_taxa_species_3 <- term_taxa_species_2%>%filter(rownames(.)%in%rgi_term_gene_sum_sel$sample)
rgi_term_gene_sum_1 <- rgi_term_gene_sum_sel%>%arrange(match(sample,rownames(term_taxa_species_3)))%>%column_to_rownames("sample")
identical(rownames(rgi_term_gene_sum_1),rownames(term_taxa_species_3))#TRUE

term_taxa_species_2_h <- decostand(term_taxa_species_3,method = "hellinger")
term_taxa_species_2_h_beta <- vegdist(term_taxa_species_2_h, method = "bray")
term_taxa_species_2_h_beta_pcoa <- cmdscale(term_taxa_species_2_h_beta, k = 2, eig = TRUE)
term_gene_h <- decostand(rgi_term_gene_sum_1,method = "hellinger")
term_gene_h_beta <- vegdist(term_gene_h, method = "bray")
term_gene_h_beta_pcoa <- cmdscale(term_gene_h_beta, k = 2, eig = TRUE)
term_taxa_gene_proc <- procrustes(X = term_taxa_species_2_h_beta_pcoa, Y = term_gene_h_beta_pcoa, symmetric = TRUE)
set.seed(1000)
term_taxa_gene_proc_p <- protest(X = term_taxa_species_2_h_beta_pcoa, Y = term_gene_h_beta_pcoa, permutations=9999)
term_taxa_gene_proc_p

#compare bacterial species between term and preterm infants.
term_resistome <- read.csv("data/args_oap_term_resistome.csv")
term_subtype <- term_resistome%>%filter(group=="subtype")%>%dplyr::select(-group)%>%column_to_rownames("names")%>%t()%>%as.data.frame()
nrow(term_subtype)#1356
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_subtype <- preterm_resistome%>%filter(group=="subtype")%>%dplyr::select(-group)%>%column_to_rownames("names")%>%t()%>%as.data.frame()
nrow(preterm_subtype)#3302
term_preterm_species <- metaphlan4_taxa%>%dplyr::select(clade_name,rownames(term_subtype),rownames(preterm_subtype))
term_preterm_species <- term_preterm_species%>%filter(grepl("k__Bacteria",clade_name))%>%filter(grepl("s__",clade_name))%>%filter(!grepl("t__",clade_name))
term_preterm_species$species <- sapply(strsplit(term_preterm_species$clade_name, split='s__', fixed=TRUE), function(x)(x[2]))
term_preterm_species <- term_preterm_species%>%dplyr::select(-clade_name)%>%column_to_rownames("species")%>%
   t()%>%as.data.frame()%>%dplyr::select(which(colSums(.)!=0))
term_preterm_species[is.na(term_preterm_species)] <- 0
names(term_preterm_species) <- gsub("_"," ",names(term_preterm_species))

term_preterm_species_dist <- vegdist(term_preterm_species, method = "bray")
metadata_cb_order <- metadata_cb%>%arrange(match(SequenceID,rownames(term_preterm_species)))
identical(metadata_cb_order$SequenceID,rownames(term_preterm_species))
set.seed(1000)
with(metadata_cb,vegan::adonis2(term_preterm_species_dist ~ Term,data=metadata_cb_order,permutations = 1000,strata = Study))
species_beta_pco <- cmdscale(term_preterm_species_dist, k=2, eig = T)
species_beta_pco_axis.1.title <- paste('PCoA1 [', 
                                       round((species_beta_pco$eig[1]/sum(species_beta_pco$eig))*100,1),
                                       '%]', sep='')
species_beta_pco_axis.2.title <- paste('PCoA2 [', 
                                       round((species_beta_pco$eig[2]/sum(species_beta_pco$eig))*100,1),
                                       '%]', sep='')

species_beta_pco.point <- as.data.frame(species_beta_pco$points)%>%rownames_to_column('SequenceID')
species_beta_pco.tabble <- right_join(metadata_cb, species_beta_pco.point, by="SequenceID")

species_beta_df.plot <- tibble(Axis1 = species_beta_pco.tabble$V1,
                               Axis2 = species_beta_pco.tabble$V2,
                               Sample_ID = species_beta_pco.tabble$SubjectID,
                               Term=species_beta_pco.tabble$Term)
#Supplementary Figure 7a
species_beta_time_p <- ggplot(data=species_beta_df.plot,aes(x=Axis1, y=Axis2, col=Term)) +
   geom_point(size=1.5, alpha=1) + 
   scale_color_manual(values = c("#E09137","#4994C4"),labels=c("Full term","Preterm"),guide = guide_legend(reverse = TRUE))+
   theme_bw()+ 
   xlab(species_beta_pco_axis.1.title) + ylab(species_beta_pco_axis.2.title) +
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"),
         plot.margin = margin(1,5,1,1,unit = "mm"))

term_preterm_species_sum <- term_preterm_species%>%rownames_to_column("SequenceID")%>%
   left_join(.,metadata_cb%>%dplyr::select(SequenceID,Term),by="SequenceID")%>%dplyr::select(-SequenceID)%>%group_by(Term)%>%
   dplyr::summarise_all(sum)%>%column_to_rownames("Term")%>%t()%>%as.data.frame()%>%rownames_to_column("species")
names(term_preterm_species_sum) <- c("species","term_total","preterm_total")
term_preterm_species_sum_n <- term_preterm_species%>%rownames_to_column("SequenceID")%>%
   left_join(.,metadata_cb%>%dplyr::select(SequenceID,Term),by="SequenceID")%>%dplyr::select(-SequenceID)%>%group_by(Term)%>%
   dplyr::summarise_all(~ sum(.>0))%>%column_to_rownames("Term")%>%t()%>%as.data.frame()%>%rownames_to_column("species")
names(term_preterm_species_sum_n) <- c("species","term_n","preterm_n")
term_preterm_species_sum_total_n <- full_join(term_preterm_species_sum,term_preterm_species_sum_n,by="species")
term_preterm_species_sum_total_n_sel <- term_preterm_species_sum_total_n%>%filter(species%in%shared_list_ncbi$ncbi_species_2)
nrow(term_preterm_species_sum_total_n_sel)#22 species
preterm_term_species_pval <- data.frame("species"="x","preterm_mean"=0,"term_mean"=0,"pval"=0,stringsAsFactors = F)[-1,]
count=0
for (s in term_preterm_species_sum_total_n_sel$species) {
   # s="Klebsiella pneumoniae"
   count=count+1
   print(count)
   a <- term_preterm_species%>%rownames_to_column("SequenceID")%>%dplyr::select(SequenceID,s)
   a_meta <- left_join(a,metadata_cb%>%dplyr::select(SequenceID,Term,Study),by="SequenceID")
   names(a_meta)[2] <- "taxa"
   a_meta$Study <- as.factor(a_meta$Study)
   a_meta$Term <- as.factor(a_meta$Term)
   a_meta_sum <- a_meta%>%group_by(Term)%>%dplyr::summarise(mean=mean(taxa))
   a_fit <- lmerTest::lmer(taxa ~ Term+(1|Study),a_meta)
   pval <- anova(a_fit)
   preterm_term_species_pval[nrow(preterm_term_species_pval)+1,] <- c(s,a_meta_sum$mean[2],a_meta_sum$mean[1],pval$`Pr(>F)`)  
}
preterm_term_species_pval <- preterm_term_species_pval%>%mutate_at(vars(-species),as.numeric)
preterm_term_species_pval$fdr <- p.adjust(preterm_term_species_pval$pval,method = "fdr")
preterm_term_species_pval$diff <- preterm_term_species_pval$preterm_mean-preterm_term_species_pval$term_mean
preterm_term_species_pval_sig <- preterm_term_species_pval%>%filter(pval<0.05)%>%arrange(desc(diff))
names(preterm_term_species_pval_sig)[1] <- "mph_species"
preterm_term_species_pval_sig_ncbi <- left_join(preterm_term_species_pval_sig,shared_list_ncbi,by=c("mph_species"="ncbi_species_2"))
preterm_term_species_pval_sig_ncbi%>%dplyr::select(species,preterm_mean,term_mean,pval,fdr,diff)

