library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(coin)
library(stringr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)


#######################################################################################################################################
#ARGs profile
#######################################################################################################################################
std.error <- function(x) sd(x)/sqrt(length(x))
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
#samples summarized by infant age (month)
preterm_metadata$month <- as.numeric(preterm_metadata$month)
#ARG types and subtypes
preterm_resistome <- read.csv("data/args_oap_preterm_resistome.csv")
preterm_type <- preterm_resistome%>%filter(group=="type")%>%dplyr::select(-group)
preterm_type[is.na(preterm_type)] <- 0
preterm_type_richness <- preterm_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"type_richness"=rowSums(.[-1]>0))
preterm_type$names[preterm_type$names=="macrolide-lincosamide-streptogramin"] <- "MLS"
preterm_type$names[preterm_type$names=="other_peptide_antibiotics"] <- "OPA"
preterm_type$names[preterm_type$names=="antibacterial_fatty_acid"] <- "AFA"
preterm_type$names[preterm_type$names=="pleuromutilin_tiamulin"] <- "PT"
preterm_type_sum <- preterm_type%>%
   dplyr::summarise("type"=names,
                    "n"=rowSums(.[-1]>0),
                    "prevalence"=rowSums(.[-1]>0)/(ncol(.)-1),
                    "Mean_copies_per_cell"=rowMeans(.[-1]))%>%
   arrange(desc(prevalence))

type_se <- data.frame("type"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in preterm_type$names) {
   a <- preterm_type%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%setNames("type")
   type_se[nrow(type_se)+1,] <- c(t,std.error(a$type))
}
preterm_type_sum_se <- left_join(preterm_type_sum,type_se,by="type")
preterm_type_sum_se$se <- as.numeric(preterm_type_sum_se$se)
preterm_type_sum_se$type <- factor(preterm_type_sum_se$type,levels = c("multidrug","polymyxin","beta_lactam","MLS","bacitracin","tetracycline",
                                                                       preterm_type_sum$type[7:28]))
preterm_type_sum_se_0.9 <- preterm_type_sum_se%>%filter(prevalence>0.9)
preterm_type_0.9 <- data.frame("total_abund"=colSums(preterm_type[-1]))%>%rownames_to_column("sample")
preterm_type_0.9_1 <- preterm_type%>%filter(names%in%preterm_type_sum_se_0.9$type)
preterm_type_0.9_2 <- data.frame("subset_abund"=colSums(preterm_type_0.9_1[-1]))%>%rownames_to_column("sample")
preterm_type_0.9_3 <- left_join(preterm_type_0.9,preterm_type_0.9_2,by="sample")
preterm_type_0.9_3$prop <- preterm_type_0.9_3$subset_abund/preterm_type_0.9_3$total_abund
preterm_type_sum_abund <- preterm_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"type_total_abund"=rowSums(.[-1]))
preterm_type_sum_se$Mean_copies_per_cell[preterm_type_sum_se$type=="multidrug"]#3.549179
preterm_type_sum_se$Mean_copies_per_cell[preterm_type_sum_se$type=="multidrug"] <- 1.4

#Figure 1b
preterm_type_p <- ggplot() +
   geom_bar(data=preterm_type_sum_se,aes(x=type, y=Mean_copies_per_cell),stat="identity",fill="#4994C4",color="black",width = 0.8)+
   geom_errorbar(data=preterm_type_sum_se,aes(x=type,ymin=Mean_copies_per_cell, ymax=Mean_copies_per_cell+se), width=0.5)+
   geom_line(data=preterm_type_sum_se,aes(x=type, y=prevalence/0.93,group=1),linetype = "dashed",color="#4994C4",linewidth=0.5)+
   geom_point(data=preterm_type_sum_se,aes(x=type, y=prevalence/0.93),size=1.5)+
   geom_hline(yintercept=0.9/0.93,linetype="dashed",size=0.5,color="black")+
   scale_y_continuous(expand = c(0,0),name= "Mean ARG abundance (capc)",limits = c(0,1.5),
                      breaks = c(0,0.5,1.0,1.4),labels = c("0","0.5","1.0","3.5"),
                      sec.axis = sec_axis(~ . *0.93, breaks = c(0,0.25,0.5,0.75,0.9,1)))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         axis.text.y = element_text(size=14, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 90,vjust = 0.4,hjust = 1),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

#ARG type dynamics over infant age
preterm_type_2 <- preterm_type%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
preterm_type_2_meta <- left_join(preterm_metadata%>%dplyr::select(SequenceID,month),preterm_type_2,by="SequenceID")%>%filter(!is.na(month))
preterm_type_2_meta_g <- gather(preterm_type_2_meta,type,abundance,-month,-SequenceID)
preterm_type_2_meta_g <- left_join(preterm_type_2_meta_g,preterm_type_sum_abund,by=c("SequenceID"="sample"))
length(unique(preterm_type_2_meta_g$SequenceID))#3301
preterm_type_2_meta_g$prop <- preterm_type_2_meta_g$abundance/preterm_type_2_meta_g$type_total_abund
preterm_type_2_meta_g_sum <- preterm_type_2_meta_g%>%group_by(month,type)%>%dplyr::summarise(mean_prop=mean(prop),mean_abun=mean(abundance))
preterm_type_2_meta_g_sum_s <- spread(preterm_type_2_meta_g_sum%>%dplyr::select(-mean_abun),month,mean_prop)

#selet types with prevalence >10% and average abundance > 0.1 capc
preterm_type_sum_sel <- preterm_type_sum%>%filter(prevalence>0.1)%>%filter(Mean_copies_per_cell>0.1)
preterm_type_2_meta_g_sum_sel <- preterm_type_2_meta_g_sum%>%filter(type%in%preterm_type_sum_sel$type)
preterm_type_2_meta_g_sum_sel$month <- factor(preterm_type_2_meta_g_sum_sel$month,levels = c("0.25","0.5","0.75","1","2","3","6","12","18"))
preterm_type_2_meta_g_sum_sel_sum <- preterm_type_2_meta_g_sum_sel%>%group_by(type)%>%dplyr::summarise(mean=mean(mean_abun))%>%arrange(desc(mean))
preterm_type_2_meta_g_sum_sel$type <- factor(preterm_type_2_meta_g_sum_sel$type,levels = preterm_type_2_meta_g_sum_sel_sum$type)

#Figure 4d
dynamic_type <- ggplot(preterm_type_2_meta_g_sum_sel, aes(x=month, y=mean_abun, fill=type)) + 
   geom_bar(stat = "identity", width=0.8)+
   scale_fill_manual(values = c("#8DD3C7", "#D0D046", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                                "#B3DE69", "#FCCDE5", "#BBEBF6", "#BC80BD", "#CCEBC5"))+
   theme_bw()+
   labs(x = "Timeponts (month)", y="Mean abundance of ARG types")+
   theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   guides(fill=guide_legend(title="ARG types", ncol=1, override.aes = list(color = "black")))+
   scale_y_continuous(expand = c(0, 0))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"))

#inflence of infant age on the ARG type abundance
#Supplementary Figure 8
type_dynamic_pvalue <- data.frame("type"="x","p_val"=0,stringsAsFactors=F)[-1,]
col_list <- c("#8DD3C7", "#D0D046", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
              "#B3DE69", "#FCCDE5", "#BBEBF6", "#BC80BD", "#CCEBC5")
count=0
for (v in levels(preterm_type_2_meta_g_sum_sel$type)) {
   # v="bacitracin"
   count=count+1
   print(v)
   db_1 <- preterm_type_2_meta_g%>%filter(type==v)
   db_1$month <- factor(db_1$month,levels = c("0.25","0.5","0.75","1","2","3","6","12","18"))
   p_val <- pvalue(kruskal_test(abundance ~ month, data=db_1))
   type_dynamic_pvalue[nrow(type_dynamic_pvalue)+1,] <- c(v,p_val)
   db_1 <- db_1%>%group_by(month)%>%dplyr::mutate(month_2=cur_group_id())
   p <- ggplot(db_1, aes(x=month_2, y=abundance)) +
      geom_smooth(method="loess", level=0.95,colour=col_list[count],linewidth=1)+
      labs(y=v)+
      scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9),
                         labels = c("1"="0.25","2"="0.5","3"="0.75","4"="1","5"="2",
                                    "6"="3","7"="6","8"="12","9"="18"))+
      theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
            panel.border = element_blank())+
      theme(panel.background = element_blank(),
            plot.background = element_blank(),
            legend.position = "none",
            legend.text = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=14, colour = "black"),
            axis.text.x = element_text(size=14, colour = "black"),
            axis.title.y = element_text(size=16, colour = "black"),
            axis.title.x = element_blank())
}

#subtype
preterm_subtype <- preterm_resistome%>%filter(group=="subtype")%>%dplyr::select(-group)
preterm_subtype[is.na(preterm_subtype)] <- 0
length(unique(preterm_subtype$names))#1993 subytpes
preterm_subtype_richness <- preterm_subtype%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"subtype_richness"=rowSums(.[-1]>0))
preterm_subtype_sum <- preterm_subtype%>%
   dplyr::summarise("subtype"=names,
                    "n"=rowSums(.[-1]>0),
                    "prevalence"=rowSums(.[-1]>0)/(ncol(.)-1),
                    "Mean_copies_per_cell"=rowMeans(.[-1]))%>%
   arrange(desc(prevalence))

preterm_subtype_sum_rare_5 <- preterm_subtype_sum%>%filter(prevalence<0.05)
nrow(preterm_subtype_sum_rare_5)/nrow(preterm_subtype_sum)#0.7431009
preterm_subtype_sum_rare_1 <- preterm_subtype_sum%>%filter(prevalence<0.01)
nrow(preterm_subtype_sum_rare_1)/nrow(preterm_subtype_sum)#0.4751631

subtype_se <- data.frame("subtype"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in preterm_subtype$names) {
   a <- preterm_subtype%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%setNames("subtype")
   subtype_se[nrow(subtype_se)+1,] <- c(t,std.error(a$subtype))
}
preterm_subtype_sum_se <- left_join(preterm_subtype_sum,subtype_se,by="subtype")
preterm_subtype_sum_se$se <- as.numeric(preterm_subtype_sum_se$se)

preterm_subtype_sum_se$type <- sapply(strsplit(preterm_subtype_sum_se$subtype, split='__', fixed=TRUE), function(x)(x[1]))
preterm_subtype_sum_se$gene <- sapply(strsplit(preterm_subtype_sum_se$subtype, split='__', fixed=TRUE), function(x)(x[2]))
preterm_subtype_sum_se$type[preterm_subtype_sum_se$type=="macrolide-lincosamide-streptogramin"] <- "MLS"
preterm_subtype_sum_se$type[preterm_subtype_sum_se$type=="other_peptide_antibiotics"] <- "OPA"
preterm_subtype_sum_se$type[preterm_subtype_sum_se$type=="pleuromutilin_tiamulin"] <- "PT"
preterm_subtype_sum_se_2 <- preterm_subtype_sum_se%>%group_by(type)%>%dplyr::summarise(n=n())%>%dplyr::arrange(desc(n))
preterm_subtype_sum_se_2$prop <- preterm_subtype_sum_se_2$n/sum(preterm_subtype_sum_se_2$n)

preterm_subtype_sum_se_0.2 <- preterm_subtype_sum_se%>%filter(Mean_copies_per_cell>0.2)#11 subtypes
preterm_subtype_0.2 <- data.frame("total_abund"=colSums(preterm_subtype[-1]))%>%rownames_to_column("sample")
preterm_subtype_0.2_1 <- preterm_subtype%>%filter(names%in%preterm_subtype_sum_se_0.2$subtype)
preterm_subtype_0.2_2 <- data.frame("subset_abund"=colSums(preterm_subtype_0.2_1[-1]))%>%rownames_to_column("sample")
preterm_subtype_0.2_3 <- left_join(preterm_subtype_0.2,preterm_subtype_0.2_2,by="sample")
preterm_subtype_0.2_3$prop <- preterm_subtype_0.2_3$subset_abund/preterm_subtype_0.2_3$total_abund
#plot the subtype with a prevalence >10%
preterm_subtype_sum_se_sel <- preterm_subtype_sum_se%>%filter(prevalence>0.10)
length(unique(preterm_subtype_sum_se_sel$subtype))
(1993-312)/1993#0.8434521
preterm_subtype_sum_se_sel_sum <- preterm_subtype_sum_se_sel%>%group_by(type)%>%dplyr::summarise(n=n())%>%dplyr::arrange(desc(n))
preterm_subtype_sum_se_nonsel <- preterm_subtype_sum_se%>%filter(prevalence<=0.10)%>%mutate(type="nocolor")
preterm_subtype_sum_se_4 <- rbind(preterm_subtype_sum_se_sel,preterm_subtype_sum_se_nonsel)
preterm_subtype_sum_se_4$type <- factor(preterm_subtype_sum_se_4$type,levels = c(preterm_subtype_sum_se_sel_sum$type,"nocolor"))

#Figure 1c
preterm_subtype_p <- ggplot(preterm_subtype_sum_se_4, aes(x=Mean_copies_per_cell, y=prevalence)) +
   geom_point(aes(fill=type),size=2.5,shape=21,color=alpha("black",0.5))+
   # geom_errorbar(aes(y=prevalence,xmin=Mean_copies_per_cell-se, xmax=Mean_copies_per_cell+se), width=0.02)+
   geom_hline(yintercept=0.1,linetype="dashed",size=0.5,color="black")+
   labs(x ="Mean ARG abundance (capc)", y="ARG prevalence")+
   theme_bw()+
   scale_y_continuous(breaks = c(0,0.10,0.25,0.5,0.75,1))+
   scale_fill_manual(values = c("#43760A", "#BEAED4", "#934B07", "#800404", "#386CB0", "#F3A2CD", 
                                "#DA732D", "#F54646", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                                "#7FC97F", "#E6AB02", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
                                "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#BBBABA"))+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=12, colour = "black"),
         axis.text.y = element_text(size=14, colour = "black"),
         axis.text.x = element_text(size=14, colour="black"),
         axis.title.y = element_text(size=14, colour = "black"),
         axis.title.x = element_text(size=14, colour = "black"))

#rgi results
rgi_profile_coverm <- read.csv("data/rgi_preterm_resistome_coverm_bins_gtdb_stats_plasmid.csv")
#ARO
rgi_profile_coverm_aro <- rgi_profile_coverm%>%group_by(sample,Best_Hit_ARO)%>%dplyr::summarise(rpkm_total=sum(rpkm))
rgi_profile_coverm_aro_richness <- rgi_profile_coverm_aro%>%group_by(sample)%>%dplyr::summarise(n=n())
rgi_profile_coverm_aro_sum <- rgi_profile_coverm_aro%>%
   group_by(Best_Hit_ARO)%>%
   dplyr::summarise("n"=n(),
                    "prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(rpkm_total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))
#compare ARG subtypes and ARGs from rgi
rgi_profile_coverm_aro_sum$Best_Hit_ARO2 <- tolower(rgi_profile_coverm_aro_sum$Best_Hit_ARO)
preterm_subtype_sum$gene <- sapply(strsplit(preterm_subtype_sum$subtype, split='__', fixed=TRUE), function(x)(x[2]))
preterm_subtype_sum$gene_2 <- tolower(preterm_subtype_sum$gene)
length(intersect(rgi_profile_coverm_aro_sum$Best_Hit_ARO2,preterm_subtype_sum$gene_2))#818 genes are overlapped by subtype and aro
818/1639#0.4990848
preterm_subtype_sum_aro_2 <- rgi_profile_coverm_aro_sum%>%filter(Best_Hit_ARO2%in%preterm_subtype_sum$gene_2)%>%
   left_join(.,preterm_subtype_sum%>%dplyr::select(subtype,prevalence,gene,gene_2),by=c("Best_Hit_ARO2"="gene_2"))

cor.test(preterm_subtype_sum_aro_2$prevalence.y,preterm_subtype_sum_aro_2$prevalence.x)#0.6919108, p-value < 2.2e-16

#Supplementary Figure 1d
corr_subtype_rgi_gene <- ggplot(preterm_subtype_sum_aro_2, aes(x=prevalence.y, y=prevalence.x)) +
   geom_point(shape=21,color=alpha("black",0.5),fill=alpha("black",0.5),size=2)+
   geom_smooth(method=lm, formula =  y ~ poly(x,1), level=0,colour="#DA732D",linewidth=2) +
   labs(x ="Pvevalence of ARG subtypes", y="Pvevalence of assembled ARGs")+
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

#organize the function
split_into_multiple <- function(column, pattern, into_prefix){
   cols <- str_split_fixed(column, pattern, n = Inf)
   # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
   cols[which(cols == "")] <- NA
   cols <- tibble::as.tibble(cols)
   # name the 'cols' tibble as 'into_prefix_1', 'into_prefix_2', ..., 'into_prefix_m' 
   # where m = # columns of 'cols'
   m <- dim(cols)[2]
   names(cols) <- paste(into_prefix, 1:m, sep = "_")
   return(cols)
}
#summary for drug
rgi_profile_coverm_2 <- rgi_profile_coverm %>%
   dplyr::select(cb,sample,rpkm,Drug.Class)%>%
   dplyr::bind_cols(split_into_multiple(.$Drug.Class, "; ", "cat"))%>%
   gather(.,cat,drug,-c(cb,sample,rpkm,Drug.Class))%>%
   filter(!is.na(drug))%>%
   dplyr::select(-cat)
rgi_profile_coverm_2_sum <- rgi_profile_coverm_2%>%group_by(sample,drug)%>%dplyr::summarise(rpkm_total=sum(rpkm))
length(unique(rgi_profile_coverm_2_sum$drug))#36 drug class
rgi_drug_richness <- rgi_profile_coverm_2_sum%>%group_by(sample)%>%dplyr::summarise(drug_richness=n())
rgi_profile_coverm_3 <- rgi_profile_coverm_2_sum%>%
   group_by(drug)%>%
   dplyr::summarise("n"=n(),
                    "prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(rpkm_total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))
names(rgi_profile_coverm_3)[1] <- "type"
rgi_profile_coverm_2_sum_1 <- rgi_profile_coverm_2%>%
   group_by(sample,drug)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%filter(rpkm_total>0)%>%spread(.,sample,rpkm_total)
which(colSums(rgi_profile_coverm_2_sum_1[-1])==0)#0
rgi_profile_coverm_2_sum_1[is.na(rgi_profile_coverm_2_sum_1)] <- 0
rgi_drug_se <- data.frame("type"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in rgi_profile_coverm_2_sum_1$drug) {
   a <- rgi_profile_coverm_2_sum_1%>%filter(drug==t)%>%column_to_rownames("drug")%>%t()%>%as.data.frame()%>%setNames("type")
   rgi_drug_se[nrow(rgi_drug_se)+1,] <- c(t,std.error(a$type))
}
rgi_profile_coverm_3_se <- left_join(rgi_profile_coverm_3,rgi_drug_se,by="type")
rgi_profile_coverm_3_se$se <- as.numeric(rgi_profile_coverm_3_se$se)
rgi_profile_coverm_3_0.9 <- rgi_profile_coverm_3%>%filter(prevalence>0.9)
rgi_profile_coverm_3_all <- data.frame("total_abund"=colSums(rgi_profile_coverm_2_sum_1[-1]))%>%rownames_to_column("sample")
rgi_profile_coverm_3_0.9_1 <- rgi_profile_coverm_2_sum_1%>%filter(drug%in%rgi_profile_coverm_3_0.9$type)
rgi_profile_coverm_3_0.9_2 <- data.frame("subset_abund"=colSums(rgi_profile_coverm_3_0.9_1[-1]))%>%rownames_to_column("sample")
rgi_profile_coverm_3_0.9_3 <- left_join(rgi_profile_coverm_3_all,rgi_profile_coverm_3_0.9_2,by="sample")
rgi_profile_coverm_3_0.9_3$prop <- rgi_profile_coverm_3_0.9_3$subset_abund/rgi_profile_coverm_3_0.9_3$total_abund
rgi_profile_coverm_3_se$type <- gsub(" antibiotic","",rgi_profile_coverm_3_se$type,fixed = T)
rgi_profile_coverm_3_se$type <- gsub("disinfecting agents and antiseptics","DAA",rgi_profile_coverm_3_se$type,fixed = T)
rgi_profile_coverm_3_se$type <- gsub("antibacterial free fatty acids","AFFA",rgi_profile_coverm_3_se$type,fixed = T)
rgi_profile_coverm_3_se$type <- factor(rgi_profile_coverm_3_se$type,levels = rgi_profile_coverm_3_se$type)

#Supplementary Figure 1f
rgi_drug_class_p <- ggplot() +
   geom_bar(data=rgi_profile_coverm_3_se, aes(x=type, y=Mean_rpkm),stat="identity",fill="#C48849",color="black",width = 0.8)+
   geom_errorbar(data=rgi_profile_coverm_3_se,aes(x=type,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=rgi_profile_coverm_3_se, aes(x=type, y=prevalence/0.000005,group=1),linetype = "dashed",color="#C48849",size=0.5)+
   geom_point(data=rgi_profile_coverm_3_se, aes(x=type, y=prevalence/0.000005),size=1.5)+
   geom_hline(yintercept=0.9/0.000005,linetype="dashed",size=0.5,color="black")+
   scale_y_continuous(expand = c(0,0),name= "Mean abundance (rpkm)",limits = c(0,280000),
                      breaks = c(0,50000,100000,150000,200000,250000),labels = c("0","50,000","100,000","150,000","200,000","250,000"),
                      sec.axis = sec_axis(~ . *0.000005, breaks = c(0,0.25,0.5,0.75,0.9,1)))+
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
#rgi drug class dynamics
rgi_profile_coverm_2_sum_2 <- rgi_profile_coverm_2_sum_1%>%column_to_rownames("drug")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
preterm_metadata_sel <- preterm_metadata%>%dplyr::select(SequenceID,month)%>%filter(SequenceID%in%rgi_profile_coverm_2_sum_2$SequenceID)
rgi_profile_coverm_2_sum_1_meta <- left_join(preterm_metadata_sel,rgi_profile_coverm_2_sum_2,by="SequenceID")%>%filter(!is.na(month))
rgi_profile_coverm_2_sum_1_meta_g <- gather(rgi_profile_coverm_2_sum_1_meta,drug,abundance,-month,-SequenceID)

rgi_drug_sum_abund <- data.frame("SequenceID"=rgi_profile_coverm_2_sum_2$SequenceID,"drug_total_abund"=rowSums(rgi_profile_coverm_2_sum_2[-1]))
rgi_profile_coverm_2_sum_1_meta_g <- left_join(rgi_profile_coverm_2_sum_1_meta_g,rgi_drug_sum_abund,by="SequenceID")
rgi_profile_coverm_2_sum_1_meta_g$prop <- rgi_profile_coverm_2_sum_1_meta_g$abundance/rgi_profile_coverm_2_sum_1_meta_g$drug_total_abund
rgi_profile_coverm_2_sum_1_meta_g_sum <- rgi_profile_coverm_2_sum_1_meta_g%>%group_by(month,drug)%>%
   dplyr::summarise(mean_prop=mean(prop),mean_abun=mean(abundance))
rgi_profile_coverm_2_sum_1_meta_g_sum_s <- spread(rgi_profile_coverm_2_sum_1_meta_g_sum%>%dplyr::select(-mean_abun),month,mean_prop)

rgi_profile_coverm_3_sel <- rgi_profile_coverm_3%>%filter(prevalence>0.1)%>%filter(Mean_rpkm>50000)
rgi_profile_coverm_3_sel$type <- gsub(" antibiotic","",rgi_profile_coverm_3_sel$type,fixed = T)
rgi_profile_coverm_2_sum_1_meta_g_sum$drug <- gsub(" antibiotic","",rgi_profile_coverm_2_sum_1_meta_g_sum$drug,fixed = T)
rgi_profile_coverm_2_sum_1_meta_g_sum_sel <- rgi_profile_coverm_2_sum_1_meta_g_sum%>%filter(drug%in%rgi_profile_coverm_3_sel$type)

rgi_profile_coverm_2_sum_1_meta_g_sum_sel$month <- factor(rgi_profile_coverm_2_sum_1_meta_g_sum_sel$month,levels = c("0.25","0.5","0.75","1","2","3","6","12","18"))
rgi_profile_coverm_2_sum_1_meta_g_sum_sel_sum <- rgi_profile_coverm_2_sum_1_meta_g_sum_sel%>%group_by(drug)%>%dplyr::summarise(mean=mean(mean_abun))%>%arrange(desc(mean))
rgi_profile_coverm_2_sum_1_meta_g_sum_sel$drug <- factor(rgi_profile_coverm_2_sum_1_meta_g_sum_sel$drug,levels = rgi_profile_coverm_2_sum_1_meta_g_sum_sel_sum$drug)
rgi_profile_coverm_2_sum_1_meta_g_sum_sel$log <- log10(rgi_profile_coverm_2_sum_1_meta_g_sum_sel$mean_abun)

#Supplementary Figure 9a
dynamic_rgi <- ggplot(rgi_profile_coverm_2_sum_1_meta_g_sum_sel, aes(x=month, y=mean_abun, fill=drug)) + 
   geom_bar(stat = "identity", width=0.8)+
   scale_fill_manual(values = c("#8DD3C7", "#D0D046", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                                "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5","#D05046",
                                "#F0A8A3","#DAD5DC","#ECEC5C","#EC695C"))+
   scale_y_continuous(expand = c(0, 0),breaks = c(0,500000,1000000,1500000,2000000),
                      labels = c("0","500,000","1,000,000","1,500,000","2,000,000"))+
   theme_bw()+
   labs(x = "Timeponts (month)", y="Mean abundance of ARG drug classes (rpkm)")+
   theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   guides(fill=guide_legend(title="ARG drug classes", ncol=1, override.aes = list(color = "black")))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "right",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black"),
         axis.title.y = element_text(size=12, colour = "black"),
         axis.title.x = element_text(size=12, colour = "black"))
#summary for gene family
rgi_profile_coverm_gene <- rgi_profile_coverm %>%
   dplyr::select(cb,sample,rpkm,AMR.Gene.Family)%>%
   dplyr::bind_cols(split_into_multiple(.$AMR.Gene.Family, "; ", "cat"))%>%
   gather(.,cat,gene,-c(cb,sample,rpkm,AMR.Gene.Family))%>%
   filter(!is.na(gene))%>%
   dplyr::select(-cat)

rgi_profile_coverm_gene_sum <- rgi_profile_coverm_gene%>%group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))
length(unique(rgi_profile_coverm_gene_sum$gene))#282 gene families
rgi_profile_coverm_gene_sum_richness <- rgi_profile_coverm_gene_sum%>%group_by(sample)%>%dplyr::summarise(drug_richness=n())
rgi_profile_coverm_gene_2 <- rgi_profile_coverm_gene_sum%>%
   group_by(gene)%>%
   dplyr::summarise("n"=n(),
                    "prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(rpkm_total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))
names(rgi_profile_coverm_gene_2)[1] <- "type"

rgi_profile_coverm_gene_sum_1 <- rgi_profile_coverm_gene%>%
   group_by(sample,gene)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%filter(rpkm_total>0)%>%spread(.,sample,rpkm_total)
rgi_profile_coverm_gene_sum_1[is.na(rgi_profile_coverm_gene_sum_1)] <- 0
rgi_gene_se <- data.frame("type"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in rgi_profile_coverm_gene_sum_1$gene) {
   a <- rgi_profile_coverm_gene_sum_1%>%filter(gene==t)%>%column_to_rownames("gene")%>%t()%>%as.data.frame()%>%setNames("type")
   rgi_gene_se[nrow(rgi_gene_se)+1,] <- c(t,std.error(a$type))
}
rgi_profile_coverm_gene_2_se <- left_join(rgi_profile_coverm_gene_2,rgi_gene_se,by="type")
rgi_profile_coverm_gene_2_se$se <- as.numeric(rgi_profile_coverm_gene_2_se$se)

rgi_profile_coverm_gene_2_rare_5 <- rgi_profile_coverm_gene_2%>%filter(prevalence<0.05)
nrow(rgi_profile_coverm_gene_2_rare_5)/nrow(rgi_profile_coverm_gene_2)#0.7765957
rgi_profile_coverm_gene_2_rare_1 <- rgi_profile_coverm_gene_2%>%filter(prevalence<0.01)
nrow(rgi_profile_coverm_gene_2_rare_1)/nrow(rgi_profile_coverm_gene_2)#0.606383

rgi_profile_coverm_gene_2_top <- rgi_profile_coverm_gene_2_se%>%filter(prevalence>0.5)
rgi_profile_coverm_gene_2_top$type <- factor(rgi_profile_coverm_gene_2_top$type,levels = rev(rgi_profile_coverm_gene_2_top$type))
rgi_profile_coverm_gene_2_top$label <- paste0("(n=",formattable::comma(rgi_profile_coverm_gene_2_top$n,digits=0),")")
rgi_profile_coverm_gene_2_top$label[rgi_profile_coverm_gene_2_top$type=="resistance-nodulation-cell division (RND) antibiotic efflux pump"]#"(n=2,860)"
rgi_profile_coverm_gene_2_top$label[rgi_profile_coverm_gene_2_top$type=="resistance-nodulation-cell division (RND) antibiotic efflux pump"] <- NA

#Supplementary Figure 1e
rgi_profile_coverm_gene_2_top_p <- ggplot() +
   geom_bar(data=rgi_profile_coverm_gene_2_top, aes(x=type, y=Mean_rpkm),stat="identity",fill="#C48849",color="black",width = 0.8)+
   geom_errorbar(data=rgi_profile_coverm_gene_2_top,aes(x=type,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=rgi_profile_coverm_gene_2_top, aes(x=type, y=prevalence*230000,group=1),linetype = "dashed",color="#C43A0A",size=1)+
   geom_point(data=rgi_profile_coverm_gene_2_top, aes(x=type, y=prevalence*230000),size=2)+
   geom_text(data=rgi_profile_coverm_gene_2_top,aes(x=type, y=Mean_rpkm+se,label=label),hjust=-0.1,vjust=0.5,size=4)+
   scale_y_continuous(expand = c(0,0),name= "Mean ARG family abundance (rpkm)",limits = c(0,250000),
                      breaks = c(0,50000,100000,150000,200000,250000),labels = c("0","50,000","100,000","150,000","200,000","250,000"),
                      sec.axis = sec_axis(~ . /230000, breaks = c(0,0.3,0.6,0.9)))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black"))+
   theme(legend.position = "none",
         legend.text = element_text(face = "italic"),
         plot.background = element_blank(),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())+
   coord_flip()

#Resistance.Mechanism
rgi_profile_coverm_mecha <- rgi_profile_coverm %>%
   dplyr::select(cb,sample,rpkm,Resistance.Mechanism)%>%
   dplyr::bind_cols(split_into_multiple(.$Resistance.Mechanism, "; ", "cat"))%>%
   gather(.,cat,mecha,-c(cb,sample,rpkm,Resistance.Mechanism))%>%
   filter(!is.na(mecha))%>%
   dplyr::select(-cat)
rgi_profile_coverm_mecha_sum <- rgi_profile_coverm_mecha%>%group_by(mecha)%>%dplyr::summarise(n=n())%>%dplyr::arrange(desc(n))
length(unique(rgi_profile_coverm_mecha_sum$mecha))#6 Resistance.Mechanism
rgi_profile_coverm_mecha_sum$prop <- rgi_profile_coverm_mecha_sum$n/sum(rgi_profile_coverm_mecha_sum$n)
rgi_profile_coverm_mecha_sum
# mecha                                   n   prop
# 1 antibiotic efflux                  143934 0.598 
# 2 antibiotic target alteration        44244 0.184 
# 3 antibiotic inactivation             24096 0.100 
# 4 antibiotic target protection        10935 0.0454
# 5 reduced permeability to antibiotic   9571 0.0398
# 6 antibiotic target replacement        7834 0.0326

#unique ARGs in different studies
preterm_subtype_sum_only1 <- preterm_subtype_sum%>%filter(n==1)
nrow(preterm_subtype_sum_only1)#187
preterm_subtype_only1 <- preterm_subtype%>%filter(names%in%preterm_subtype_sum_only1$subtype)%>%
   dplyr::select(names,which(colSums(.[-1])!=0)+1)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%
   rownames_to_column("SequenceID")%>%dplyr::summarise("SequenceID"=SequenceID,"n"=rowSums(.[-1]>0))%>%
   left_join(.,preterm_metadata%>%dplyr::select(SequenceID,Study,clean_reads),by="SequenceID")
sum(preterm_subtype_only1$n)#187
preterm_subtype_only1$clean_reads_2 <- preterm_subtype_only1$clean_reads/1000000
cor.test(preterm_subtype_only1$clean_reads,preterm_subtype_only1$n)#r = 0.04221598 , p-value = 0.6814

preterm_subtype_only1_sum <- preterm_subtype_only1%>%group_by(Study)%>%dplyr::summarise(total=sum(n))%>%arrange(desc(total))
preterm_subtype_only1_sum$total[1]#111
preterm_subtype_only1_sum$total[1] <- 30

#different from ARG subtypes, ARG familiy is obtained from rgi, not coverm
rgi_profile_coverm_gene_2_only1 <- rgi_profile_coverm_gene_2%>%filter(n==1)
nrow(rgi_profile_coverm_gene_2_only1)#33

preterm_rgi_only1 <- rgi_profile_coverm_gene_sum%>%filter(gene%in%rgi_profile_coverm_gene_2_only1$type)%>%
   dplyr::group_by(sample)%>%dplyr::summarise(n=n())%>%
   left_join(.,preterm_metadata%>%dplyr::select(SequenceID,Study,clean_reads),by=c("sample"="SequenceID"))
sum(preterm_rgi_only1$n)#33
preterm_rgi_only1$clean_reads_2 <- preterm_rgi_only1$clean_reads/1000000

preterm_rgi_only1_sum <- preterm_rgi_only1%>%group_by(Study)%>%dplyr::summarise(total=sum(n))%>%arrange(desc(total))
preterm_rgi_only1_sum_subtype <- rbind(preterm_rgi_only1_sum%>%mutate(group="ARG families"),preterm_subtype_only1_sum%>%mutate(group="ARG subtypes"))
preterm_rgi_only1_sum_subtype$Study <- factor(preterm_rgi_only1_sum_subtype$Study,levels = rev(c(as.character(preterm_subtype_only1_sum$Study),
                                                                                                 "RahmanSF_2018")))
preterm_rgi_only1_sum_subtype$group <- factor(preterm_rgi_only1_sum_subtype$group,levels = c("ARG subtypes","ARG families"))

#Supplementary Figure 1g
subtype_rgi_only_p <- ggplot(data=preterm_rgi_only1_sum_subtype, aes(x=total, y=Study,fill=group)) +
   geom_bar(stat="identity",color="black",width = 0.8)+
   facet_wrap(~group,scales = "free_x")+
   scale_fill_manual(values = c("#4994C4","#C48849"))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=12),
         panel.spacing = unit(5,"mm"),
         plot.background = element_blank(),
         axis.text.y = element_text(size=8, colour = "black"),
         axis.text.x = element_text(size=8, colour="black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())
#Supplementary Figure 1h
subtype_only_p_2 <- ggplot(preterm_subtype_only1, aes(x=clean_reads_2, y=n)) +
   geom_point(size=1,color="#4994C4")+
   geom_smooth(method=lm, formula =  y ~ poly(x,1), level=0,colour="grey",linewidth=1)+
   labs(x ="Number of clean reads (million)", y="Number of ARG subtypes")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         plot.background = element_blank(),
         axis.text.y = element_text(size=8, colour = "black"),
         axis.text.x = element_text(size=8, colour="black"),
         axis.title.y = element_text(size=8, colour="black"),
         axis.title.x = element_blank())
#Supplementary Figure 1h
rgi_only_p_2 <- ggplot(preterm_rgi_only1, aes(x=clean_reads_2, y=n)) +
   geom_point(size=1,color="#C48849")+
   geom_smooth(method=lm, formula =  y ~ poly(x,1), level=0,colour="grey",linewidth=1)+
   labs(x ="Number of clean reads (million)", y="Number of ARG families")+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(panel.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         plot.background = element_blank(),
         axis.text.y = element_text(size=8, colour = "black"),
         axis.text.x = element_text(size=8, colour="black"),
         axis.title = element_text(size=8, colour="black"))
#correlation between arg subytpes and rgi aro
rgi_gene_richness_2 <- rgi_profile_coverm%>%group_by(sample,Best_Hit_ARO)%>%dplyr::summarise(total_rpkm=sum(rpkm))%>%
   ungroup()%>%group_by(sample)%>%
   dplyr::summarise(gene_richness=n())
preterm_subtype_richness_rgi <- left_join(preterm_subtype_richness,rgi_gene_richness_2,by="sample")
preterm_subtype_richness_rgi$gene_richness[is.na(preterm_subtype_richness_rgi$gene_richness)] <- 0
quantile(preterm_subtype_richness_rgi$subtype_richness)
# 0%  25%  50%  75% 100% 
# 9  101  131  169  488 
quantile(preterm_subtype_richness_rgi$gene_richness)
# 0%  25%  50%  75% 100% 
# 0   33   52   71  156 
cor.test(preterm_subtype_richness_rgi$gene_richness,preterm_subtype_richness_rgi$subtype_richness)#0.6080681, p-value < 2.2e-16
#Supplementary Figure 1b
corr_subtype_rgi <- ggplot(preterm_subtype_richness_rgi, aes(x=subtype_richness, y=gene_richness)) +
   geom_point(shape=21,color=alpha("black",0.5),fill=alpha("black",0.5),size=2)+
   geom_smooth(method=lm, formula =  y ~ poly(x,1), level=0,colour="#DA732D",linewidth=2) +
   scale_x_continuous(limits = c(0,510))+
   labs(x ="Number of ARG subtypes", y="Number of assembled ARG")+
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
#Supplementary Figure 1b
density_subtypes <- ggplot(preterm_subtype_richness_rgi, aes(x=subtype_richness)) +
   # geom_histogram(bins = 30,aes(y=..density..))+
   geom_density(fill="#4994C4",alpha=0.5)+
   scale_y_continuous(expand = c(0,0))+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=10, colour = "black"),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())
#Supplementary Figure 1b
density_rgi_gene <- ggplot(preterm_subtype_richness_rgi, aes(x=gene_richness)) +
   geom_density(fill="#C48849",alpha=0.5)+
   theme_bw()+
   scale_y_continuous(expand = c(0,0),position = 'right')+
   scale_x_reverse()+
   theme(panel.grid = element_blank(),axis.line = element_line(colour = "black"),
         panel.border = element_blank())+
   guides(fill=guide_legend(ncol=5))+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         legend.text = element_text(size=10, colour = "black"),
         axis.text.y = element_text(size=10, colour = "black"),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

#######################################################################################################################################
#ARGS host
#######################################################################################################################################
rgi_profile_stats <- read.csv("data/rgi_preterm_resistome_coverm_bins_gtdb_stats_plasmid.csv")
preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")

rgi_profile_stats$domain <- sapply(strsplit(rgi_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[1]))
rgi_profile_stats$phylum <- sapply(strsplit(rgi_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[2]))
rgi_profile_stats$class <- sapply(strsplit(rgi_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[3]))
rgi_profile_stats$order <- sapply(strsplit(rgi_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[4]))
rgi_profile_stats$family <- sapply(strsplit(rgi_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[5]))
rgi_profile_stats$genus <- sapply(strsplit(rgi_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[6]))
rgi_profile_stats$species <- sapply(strsplit(rgi_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[7]))
rgi_profile_stats$domain[is.na(rgi_profile_stats$domain)] <- "unknown"
rgi_profile_stats$phylum[is.na(rgi_profile_stats$phylum)] <- "unknown"
rgi_profile_stats$class[is.na(rgi_profile_stats$class)] <- "unknown"
rgi_profile_stats$order[is.na(rgi_profile_stats$order)] <- "unknown"
rgi_profile_stats$family[is.na(rgi_profile_stats$family)] <- "unknown"
rgi_profile_stats$genus[is.na(rgi_profile_stats$genus)] <- "unknown"
rgi_profile_stats$species[rgi_profile_stats$species=="s__"] <- "unknown"
rgi_profile_stats$species[is.na(rgi_profile_stats$species)] <- "unknown"

rgi_profile_stats_domain <- rgi_profile_stats%>%group_by(domain)%>%dplyr::summarise(n=n())
tax_number <- data.frame("tax"="x","n_tax"=0,"n_args"=0,stringsAsFactors = F)[-1,]
tax_number[1,] <- c("Phylum",length(table(rgi_profile_stats$phylum))-1,sum(table(rgi_profile_stats$phylum))-nrow(rgi_profile_stats%>%filter(phylum=="unknown")))
tax_number[2,] <- c("Class",length(table(rgi_profile_stats$class))-1,sum(table(rgi_profile_stats$class))-nrow(rgi_profile_stats%>%filter(class=="unknown")))
tax_number[3,] <- c("Order",length(table(rgi_profile_stats$order))-1,sum(table(rgi_profile_stats$order))-nrow(rgi_profile_stats%>%filter(order=="unknown")))
tax_number[4,] <- c("Family",length(table(rgi_profile_stats$family))-1,sum(table(rgi_profile_stats$family))-nrow(rgi_profile_stats%>%filter(family=="unknown")))
tax_number[5,] <- c("Genus",length(table(rgi_profile_stats$genus))-1,sum(table(rgi_profile_stats$genus))-nrow(rgi_profile_stats%>%filter(genus=="unknown")))
tax_number[6,] <- c("Species",length(table(rgi_profile_stats$species))-1,sum(table(rgi_profile_stats$species))-nrow(rgi_profile_stats%>%filter(species=="unknown")))
tax_number$n_tax <- as.numeric(tax_number$n_tax)
tax_number$n_args <- as.numeric(tax_number$n_args)
tax_number$tax <- factor(tax_number$tax,levels = c("Phylum", "Class","Order","Family","Genus","Species"))
#Supplementary Figure 2b
tax_number_plot <- ggplot(tax_number, aes(x=tax, y=n_tax)) + 
   geom_bar(stat = "identity",fill="#72A8C4",color="black",width = 0.8)+
   geom_text(aes(label=paste0("(n=",formattable::comma(n_tax,digits=0),"; n=",formattable::comma(n_args,digits=0),")"),hjust=-0.1,size=6))+
   theme_bw() + # remove the backgroud
   labs(x = element_blank(), y=element_blank())+
   theme(axis.title.y = element_text(size = 16))+
   theme(panel.border = element_blank(), axis.line = element_blank()) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the grid
   coord_flip()+
   theme(legend.position="none",
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_text(size=14, colour = "black"),
         axis.text.x = element_blank(),
         axis.ticks = element_blank())
#get the top 10 taxa in each rank
tax=c("phylum","class","order","family","genus")
tax_diver <- data.frame()
for (t in tax) {
   #t="domain"
   a <- rgi_profile_stats%>%dplyr::group_by(get(t))%>%dplyr::summarise(n=n())%>%arrange(desc(n))
   print(sum(a$n))
   b <- a[c(1:11),] 
   b$propor <- b$n/sum(a$n)
   b[nrow(b)+1,] <- list("Others",nrow(rgi_profile_stats)-sum(b$n,na.rm = T),1-sum(b$propor,na.rm = T))
   b$taxonomy <- t
   tax_diver <- rbind(tax_diver,b)
}
names(tax_diver)[1] <- "name" 
tax_diver$name[tax_diver$name=="unknown"] <- "Unknown"
tax_diver$propor <- as.numeric(tax_diver$propor)*100
tax_diver$propor <- round(tax_diver$propor,2)

tax_diver_2 <- tax_diver%>%group_by(taxonomy)%>%dplyr::mutate(new_t=c(1:12))
tax_diver_2$new_t[tax_diver_2$new_t==1] <- 2
tax_diver_2$new_t[tax_diver_2$name=="Unknown"] <- 1
tax_diver_2$new_t <- as.character(tax_diver_2$new_t)
tax_diver_2$name <- gsub(".__","",tax_diver_2$name)
tax_diver_2$taxonomy <- paste(toupper(substr(tax_diver_2$taxonomy, 1, 1)), substr(tax_diver_2$taxonomy, 2, nchar(tax_diver_2$taxonomy)), sep="")
tax_diver_2$taxonomy <- factor(tax_diver_2$taxonomy,levels = c("Domain","Phylum", "Class","Order","Family","Genus"))
tax_diver_2$new_t <- factor(tax_diver_2$new_t, levels = c("1","12","11","10","9","8","7","6","5","4","3","2"))
#Supplementary Figure 2a
rgi_taxa_all <- ggplot(tax_diver_2, aes(x=taxonomy, y=n, fill=new_t)) + 
   geom_bar(stat = "identity", width=0.8)+
   scale_fill_manual(values = c("#CACACA","#8F8F8F","#6AABF5","#BEBADA","#FB8072","#80B1D3","#FDB462",
                                "#B3DE69","#FCCDE5","#CABB57","#BC80BD","#8DD3C7"))+
   theme_bw() +
   labs(x = element_blank(), y="Number of ARG ORFs")+
   theme(axis.title.y = element_text(size = 16))+
   theme(panel.border = element_blank(), axis.line = element_line()) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the grid
   scale_y_continuous(expand = c(0,0),breaks=c(0, 50000, 100000,150000,200000),
                      labels=c("0","50,000", "100,000", "150,000","200,000"))+
   theme(legend.position="none",
         aspect.ratio = 1.3/1,
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 45,vjust = 0.7))

#only species
std.error <- function(x) sd(x)/sqrt(length(x))
rgi_profile_stats_species <- rgi_profile_stats%>%dplyr::select(sample,species,rpkm)
rgi_profile_stats_species$species <- gsub("s__","",rgi_profile_stats_species$species,fixed = T)
rgi_profile_stats_species_sum <- rgi_profile_stats_species%>%group_by(sample,species)%>%dplyr::summarise(total=sum(rpkm))
length(unique(rgi_profile_stats_species_sum$sample))#3300
rgi_profile_stats_species_sum_1 <- rgi_profile_stats_species_sum%>%
   group_by(species)%>%
   dplyr::summarise("prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))

rgi_profile_stats_species_sum_s <- spread(rgi_profile_stats_species_sum,sample,total)
rgi_profile_stats_species_sum_s[is.na(rgi_profile_stats_species_sum_s)] <- 0
rgi_species_se <- data.frame("species"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in rgi_profile_stats_species_sum_s$species) {
   a <- rgi_profile_stats_species_sum_s%>%filter(species==t)%>%column_to_rownames("species")%>%t()%>%as.data.frame()%>%setNames("species")
   rgi_species_se[nrow(rgi_species_se)+1,] <- c(t,std.error(a$species))
}
rgi_profile_stats_species_sum_2 <- left_join(rgi_profile_stats_species_sum_1,rgi_species_se,by="species")
rgi_profile_stats_species_sum_2$se <- as.numeric(rgi_profile_stats_species_sum_2$se)
rgi_profile_stats_species_sum_2 <- rgi_profile_stats_species_sum_2%>%filter(species!="unknown")
rgi_profile_stats_species_sum_3 <- rgi_profile_stats_species_sum_2%>%dplyr::filter(prevalence>0.01)
rgi_profile_stats_species_2 <- rgi_profile_stats_species%>%group_by(species)%>%dplyr::summarise(n_args=n())%>%
   filter(species%in%rgi_profile_stats_species_sum_3$species)%>%arrange(desc(n_args))
rgi_profile_stats_species_sum_3 <- left_join(rgi_profile_stats_species_sum_3,rgi_profile_stats_species_2,by="species")
rgi_profile_stats_species_sum_3$species <- factor(rgi_profile_stats_species_sum_3$species,levels = rev(rgi_profile_stats_species_sum_3$species))
rgi_profile_stats_species_sum_3$label <- paste0("(n=",formattable::comma(rgi_profile_stats_species_sum_3$n_args,digits=0),")")
rgi_profile_stats_species_sum_3$label[rgi_profile_stats_species_sum_3$species=="Escherichia coli"]#"(n=52,569)"
rgi_profile_stats_species_sum_3$label[rgi_profile_stats_species_sum_3$species=="Escherichia coli"] <- NA
rgi_profile_stats_species_sum_3$label[rgi_profile_stats_species_sum_3$species=="Klebsiella pneumoniae"]#"(n=20,907)"
rgi_profile_stats_species_sum_3$label[rgi_profile_stats_species_sum_3$species=="Klebsiella pneumoniae"] <- NA
rgi_profile_stats_species_sum_3$label[rgi_profile_stats_species_sum_3$species=="Enterococcus faecalis"]#"(n=8,744)"
rgi_profile_stats_species_sum_3$label[rgi_profile_stats_species_sum_3$species=="Enterococcus faecalis"] <- NA
rgi_profile_stats_species_sum_3$Mean_rpkm[rgi_profile_stats_species_sum_3$species=="Klebsiella pneumoniae"]#82729.38
rgi_profile_stats_species_sum_3$Mean_rpkm[rgi_profile_stats_species_sum_3$species=="Klebsiella pneumoniae"] <- 55000
rgi_profile_stats_species_sum_3$Mean_rpkm[rgi_profile_stats_species_sum_3$species=="Escherichia coli"]#149144.6
rgi_profile_stats_species_sum_3$Mean_rpkm[rgi_profile_stats_species_sum_3$species=="Escherichia coli"] <- 75000
#Figure 1d
rgi_taxa_species <- ggplot() +
   geom_bar(data=rgi_profile_stats_species_sum_3, aes(x=species, y=Mean_rpkm),stat="identity",fill="#C48849",color="black",width = 0.8)+
   geom_errorbar(data=rgi_profile_stats_species_sum_3,aes(x=species,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=rgi_profile_stats_species_sum_3, aes(x=species, y=prevalence*80000,group=1),linetype = "dashed",color="#C43A0A",size=1)+
   geom_point(data=rgi_profile_stats_species_sum_3, aes(x=species, y=prevalence*80000),size=1.5)+
   geom_text(data=rgi_profile_stats_species_sum_3,aes(x=species, y=Mean_rpkm+se,label=label),hjust=-0.1,vjust=0.5,size=3)+
   scale_y_continuous(expand = c(0,0),name= "Mean abundance (rpkm)",limits = c(0,81000),
                      breaks = c(0,15000,30000,55000,70000),labels = c("0","15,000","30,000","82,729","149,145"),
                      sec.axis = sec_axis(~ . /80000, breaks = c(0,0.15,0.30,0.45)))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(legend.position = "none",
         legend.text = element_text(face = "italic"),
         plot.background = element_blank(),
         panel.background = element_blank(),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         axis.text.y = element_text(size=11, colour = "black",face="italic"),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=12, colour = "black")
   )+
   theme(aspect.ratio=3/1)+
   coord_flip()

#plasmid
rgi_profile_stats$bin[is.na(rgi_profile_stats$bin)] <- "nobin"
rgi_profile_stats_unknown_taxa <- rgi_profile_stats%>%filter(bin=="nobin")
rgi_profile_stats_unknown_taxa_sum <- rgi_profile_stats_unknown_taxa%>%
   group_by(plasmid)%>%dplyr::summarise(n=n())%>%
   mutate(Propotion=n/sum(n))%>%
   mutate(group="Unknown")
rgi_profile_stats_known_taxa <- rgi_profile_stats%>%filter(bin!="nobin")
rgi_profile_stats_known_taxa_sum <- rgi_profile_stats_known_taxa%>%
   group_by(plasmid)%>%dplyr::summarise(n=n())%>%
   mutate(Propotion=n/sum(n))%>%
   mutate(group="Known")
rgi_profile_stats_unknown_known <- rbind(rgi_profile_stats_unknown_taxa_sum,rgi_profile_stats_known_taxa_sum)
chisq.test(rgi_profile_stats$plasmid,rgi_profile_stats$bin,correct = T)
plas_preterm_sum <- rgi_profile_stats%>%group_by(plasmid)%>%dplyr::summarise(n=n())
17546/(205103+17546)#0.07880565
#Supplementary Figure 3a
pie(plas_preterm_sum$n,
    labels = c("Chromosome \n(n=205,103)","Plasmid \n(n=17,546)"),
    col = c("#A8DFFC","#F7A48F"),
    init.angle=360, radius = 1,cex=1)

plas_sum_list <- list()
for (v in c("phylum","class","order","family","genus","species")) {
   # v="genus"
   plas_sum <- data.frame("taxa"="x","pred_n"=0,"prediction"="x","args_n"=0,"args_n_chrom"=0,"args_n_plasmid"=0,stringsAsFactors = F)[-1,]
   count=0
   taxa_list <- rgi_profile_stats%>%filter(get(v)!="unknown")%>%distinct(get(v),.keep_all = F)%>%setNames("taxa")
   for (i in taxa_list$taxa) {
      # i="g__Klebsiella"
      count=count+1
      print(count)
      a <- rgi_profile_stats%>%filter(get(v)==i)
      a_1 <- length(unique(a$plasmid))
      if (a_1==1) {
         plas_sum[nrow(plas_sum)+1,] <- c(i,1,unique(a$plasmid),nrow(a),NA,NA)
      }else{
         plas_sum[nrow(plas_sum)+1,] <- c(i,2,NA,nrow(a),table(a$plasmid)[1],table(a$plasmid)[2])
      }
   }
   plas_sum$pred_n <- as.numeric(plas_sum$pred_n)
   plas_sum$args_n <- as.numeric(plas_sum$args_n)
   plas_sum$args_n_chrom <- as.numeric(plas_sum$args_n_chrom)
   plas_sum$args_n_plasmid <- as.numeric(plas_sum$args_n_plasmid)
   plas_sum_list[[v]] <- plas_sum
}
plas_sum_species_chom <- plas_sum_species%>%filter(prediction=="No")%>%arrange(desc(args_n))
plas_sum_species_chom_sel <- plas_sum_species_chom%>%filter(args_n>1)
length(unique(plas_sum_species_chom_sel$taxa))#375
plas_sum_species_plasmid <- plas_sum_species%>%filter(prediction=="Yes")%>%arrange(desc(args_n))%>%filter(args_n>1)
length(unique(plas_sum_species_plasmid$taxa))#1
plas_sum_species_both <- plas_sum_species%>%filter(pred_n==2)
length(unique(plas_sum_species_both$taxa))#128

plas_tax_number <- data.frame("tax"="x","n_tax"=0,stringsAsFactors = F)[-1,]
plas_tax_number[1,] <- c("Chromosome",375)
plas_tax_number[2,] <- c("Both",128)
plas_tax_number[3,] <- c("Plasmid",1)
plas_tax_number$n_tax <- as.numeric(plas_tax_number$n_tax)
#Supplementary Figure 3b
pie(plas_tax_number$n_tax,
    labels = c("Chromosome \n(n=375)","Both \n(n=128)","Plasmid \n(n=1)"),
    col = c("#A8DFFC","#ACFCA8","#F7A48F"),
    init.angle=30, radius = 1,cex=1)

plas_sum_species_both_sel <- plas_sum_species_both%>%filter(args_n>100)
plas_sum_species_both_sel$taxa <- gsub("s__","",plas_sum_species_both_sel$taxa,fixed = T)
plas_sum_species_both_sel$chrom_prop <- plas_sum_species_both_sel$args_n_chrom/plas_sum_species_both_sel$args_n
plas_sum_species_both_sel$plasmid_prop <- plas_sum_species_both_sel$args_n_plasmid/plas_sum_species_both_sel$args_n
plas_sum_species_both_sel <- plas_sum_species_both_sel%>%arrange(desc(plasmid_prop))
plas_sum_species_both_sel_1 <- plas_sum_species_both_sel%>%dplyr::select(-c(2,3,4,6,8))%>%setNames(c("taxa","n_args","n_prop"))%>%mutate(status="chrom")
plas_sum_species_both_sel_2 <- plas_sum_species_both_sel%>%dplyr::select(-c(2,3,4,5,7))%>%setNames(c("taxa","n_args","n_prop"))%>%mutate(status="plasmid")
plas_sum_species_both_sel_3 <- rbind(plas_sum_species_both_sel_1,plas_sum_species_both_sel_2)
plas_sum_species_both_sel_3$text <- paste0("(n=",formattable::comma(plas_sum_species_both_sel_3$n_args,format = "d"),")")
plas_sum_species_both_sel_3$text <- gsub(" ","",plas_sum_species_both_sel_3$text,fixed = T)
plas_sum_species_both_sel_3$taxa <- factor(plas_sum_species_both_sel_3$taxa,levels = plas_sum_species_both_sel$taxa)
plas_sum_species_both_sel_3 <- plas_sum_species_both_sel_3%>%mutate(position=case_when(status=="chrom" ~ 0.75,
                                                                                       status=="plasmid" ~ 0.01))
#Supplementary Figure 3c
plasmid_proportion <- ggplot(data=plas_sum_species_both_sel_3, aes(x=taxa, y=n_prop, fill=status)) +
   geom_bar(stat="identity")+
   geom_hline(yintercept = 0.5,linetype="dashed", color = "grey")+
   geom_text(aes(label=text,y=position),hjust=-0.05,size=3,angle=90)+
   scale_fill_manual(values = c("#A8DFFC","#F7A48F"),label=c("chrom"="Chromosome","plasmid"="Plasmid"))+
   scale_y_continuous(expand = c(0,0))+
   labs(y = "Proportion")+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),axis.line = element_blank(),panel.grid.major = element_line(linewidth = 1),
         panel.border = element_rect(colour = "black",linewidth = 1.5))+
   theme(legend.text = element_text(size=16),
         aspect.ratio = 1/6,
         legend.position = "top",
         axis.text.y = element_text(size=16, colour = "black"),
         axis.text.x = element_text(size=16, colour = "black",face="italic",angle = 90,hjust = 1,vjust = 0.5),
         axis.title.x = element_blank(),
         plot.margin = unit(c(1,1,1,1), "cm"),
         axis.title.y = element_text(size=16, colour = "black"))

#cluster of ARG ORFs from plasmid
plasme_huaxi <- read.delim("data/plasme_cb_huaxi.csv")
plasme_preterm_repeat <- plasme_huaxi%>%filter(grepl("repeat",contig))
plasme_preterm_repeat$contig2 <- plasme_preterm_repeat$contig
plasme_preterm_repeat$contig2 <- gsub("repeat","",plasme_preterm_repeat$contig2,fixed = T)
rgi_profile_stats_plasmid_yes_orf <- rgi_profile_stats%>%filter(bin!="nobin")%>%filter(plasmid=="Yes")
length(unique(rgi_profile_stats_plasmid_yes_orf$ORF_ID))#3377
rgi_profile_stats_plasmid_yes_orf$run_contig <- gsub("__","_",rgi_profile_stats_plasmid_yes_orf$run_contig,fixed = T)
rgi_profile_stats_plasmid_yes_orf_1 <- rgi_profile_stats_plasmid_yes_orf%>%filter(!run_contig%in%plasme_preterm_repeat$contig2)
rgi_profile_stats_plasmid_yes_orf_2 <- rgi_profile_stats_plasmid_yes_orf%>%filter(run_contig%in%plasme_preterm_repeat$contig2)
rgi_profile_stats_plasmid_yes_orf_2$sample <- paste0(rgi_profile_stats_plasmid_yes_orf_2$sample,"repeat")
rgi_profile_stats_plasmid_yes_orf_cb <- rbind(rgi_profile_stats_plasmid_yes_orf_1,rgi_profile_stats_plasmid_yes_orf_2)
rgi_profile_stats_plasmid_yes_orf_cb$sample_orf <- paste0(rgi_profile_stats_plasmid_yes_orf_cb$sample,"_rgi_card_",rgi_profile_stats_plasmid_yes_orf_cb$ORF_ID)

rgi_cluster_plas <- read.delim("data/rgi_plasmid_orf_nr_cluster.txt")
length(unique(rgi_cluster_plas$clstr))#612 clusters formed
length(unique(rgi_cluster_plas$id))#3377 args
rgi_profile_stats_plas_cluster <- left_join(rgi_profile_stats_plasmid_yes_orf_cb,rgi_cluster_plas%>%dplyr::select(id,clstr),by=c("sample_orf"="id"))

for (v in c("phylum","class","order","family","genus","species")) {
   print(v)
   cluster_sum_1 <- rgi_profile_stats_plas_cluster%>%group_by(clstr)%>%dplyr::summarise(n_taxa_all=length(unique(get(v))),n_args=n())
   cluster_sum_2 <- rgi_profile_stats_plas_cluster%>%filter(get(v)!="unknown")%>%group_by(clstr)%>%dplyr::summarise(n_taxa_sel=length(unique(get(v))),n_args_2=n())
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
#153 clusters with at least 2 species including unknown species; 152 clusters with at least 2 known species
cluster_sum_sel_species <- rgi_profile_stats_plas_cluster%>%filter(clstr%in%cluster_sum_sel$clstr)
length(unique(cluster_sum_sel_species$clstr))#152 clusters
length(unique(cluster_sum_sel_species$species))-1#113 species
cluster_sum_sel_species_cl <- cluster_sum_sel_species%>%group_by(clstr)%>%dplyr::summarise(n_species=length(unique(species)))%>%arrange(desc(n_species))
cluster_sum_sel_species_cl_orf <- cluster_sum_sel_species%>%group_by(clstr,Best_Hit_ARO)%>%dplyr::summarise(n_ARGs=n())
cluster_sum_sel_species_cl_cb <- left_join(cluster_sum_sel_species_cl,cluster_sum_sel_species_cl_orf,by="clstr")

cluster_sum_sel_species_sum <- cluster_sum_sel_species%>%filter(species!="unknown")%>%group_by(clstr,species)%>%dplyr::summarise(n_args=n())
cluster_sum_sel_species_sum_1 <- cluster_sum_sel_species_sum%>%group_by(species)%>%dplyr::summarise(n=n())%>%arrange(desc(n))
cluster_sum_sel_species_sum_2 <- cluster_sum_sel_species_sum%>%group_by(clstr)%>%dplyr::summarise(n=n())%>%filter(n>=10)
cluster_sum_sel_species_sum_sel <- cluster_sum_sel_species_sum%>%
   filter(clstr%in%cluster_sum_sel_species_sum_2$clstr)%>%
   mutate(n_args_2=case_when(n_args>100 ~ "100",
                             TRUE ~ as.character(n_args)))%>%
   dplyr::select(-n_args)
cluster_sum_sel_species_sum_sel$species <- gsub("s__","",cluster_sum_sel_species_sum_sel$species,fixed = T)
cluster_sum_sel_species_sum_sel$n_args_2 <- as.numeric(cluster_sum_sel_species_sum_sel$n_args_2)
quantile(cluster_sum_sel_species_sum_sel$n_args_2)
cluster_sum_sel_species_s <- spread(cluster_sum_sel_species_sum_sel,species,n_args_2)%>%column_to_rownames("clstr")
cluster_sum_sel_species_s[is.na(cluster_sum_sel_species_s)] <- 0
names(cluster_sum_sel_species_s)
col_fun_type <- colorRamp2(c(0,3,10,20,30,50,100),c("white","#A8C6EC","#7FA4D3","#5A87BF","#3768A4","#184B8A","#073876"))
#Supplementary Figure 3d
plasmid_cluster <- Heatmap(as.matrix(cluster_sum_sel_species_s),cluster_rows=T,cluster_columns=T,col = col_fun_type,
                           show_row_names = T,show_column_names = T,
                           border = T,row_names_gp = gpar(fontsize =10),
                           rect_gp = gpar(col = "grey", lwd = 0.7),
                           column_names_gp = gpar(fontface="italic",fontsize=10),
                           border_gp=gpar(col = "black", lwd = 1),column_gap = unit(2, "mm"),
                           width = unit(35, "cm"), height = unit(6, "cm"),
                           show_heatmap_legend = TRUE,
                           heatmap_legend_param = list(title = "Number of ARG ORFs",direction = "horizontal", 
                                                       legend_height = unit(0.3, "cm"),legend_width = unit(4, "cm"),
                                                       at = c(0,5,20,30,50,100), labels = c("0","5","20","30","50","≥100"))
)
cluster_shared <- list()
length(unique(cluster_sum_sel_species_sum$species))#113
for (s in unique(cluster_sum_sel_species_sum$species)) {
   a <- cluster_sum_sel_species_sum%>%filter(species==s)
   b <- cluster_sum_sel_species_sum%>%filter(species!=s)%>%filter(clstr%in%a$clstr)
   b_sum <- b%>%group_by(species)%>%dplyr::summarise(n=n())%>%arrange(desc(n))
   cluster_shared[[s]] <- b_sum
}

#######################################################################################################################################
#ARGs functions from EggNOG-mapper
#######################################################################################################################################
egg_args_sum <- read.csv("data/egg_args_sum.csv")
egg_args_sum%>%group_by(cat)%>%dplyr::summarise(n=n())
egg_args_sum_top30 <- egg_args_sum%>%group_by(cat)%>%dplyr::slice_head(n=30)
#Supplementary Figure 4
#COG
egg_COG <- egg_args_sum_top30%>%filter(cat=="COG_category")
egg_COG$number <- factor(egg_COG$number,levels = rev(egg_COG$number))
egg_COG$text <- paste0("(n=",formattable::comma(egg_COG$n_gene,format = "d"),")")
egg_COG$text <- gsub(" ","",egg_COG$text,fixed = T)
egg_COG_plot <- ggplot(egg_COG, aes(x=n_gene, y=number)) +
   geom_bar(stat="identity", width=0.8,fill="#095A9A")+
   geom_text(aes(label=text),hjust=-0.05,size=3)+
   theme_bw()+
   theme(panel.grid = element_blank())+
   theme(plot.background = element_blank())+
   theme(panel.border = element_rect(colour = "black",linewidth=1),axis.line = element_blank())+
   labs(y="", x= "")+
   scale_x_continuous(expand = c(0,0),limits = c(0,80000),
                      breaks=c(0, 20000,40000,60000),
                      labels=c("0","20,000", "40,000", "60,000"))+
   theme(legend.position="none",
         axis.text.x = element_text(size=15, colour = "black",angle = 90,vjust = 0.5,hjust = 1),
         axis.text.y = element_text(size=15, colour = "black",hjust = 0))

#GOs
egg_GOs <- egg_args_sum_top30%>%filter(cat=="GOs")
egg_GOs$number <- factor(egg_GOs$number,levels = rev(egg_GOs$number))
egg_GOs$text <- paste0("(n=",formattable::comma(egg_GOs$n_gene,format = "d"),")")
egg_GOs_plot <- ggplot(egg_GOs, aes(x=n_gene, y=number)) +
   geom_bar(stat="identity", width=0.7,fill="#099A35")+
   geom_text(aes(label=text),hjust=-0.05,size=3)+
   theme_bw()+
   theme(panel.grid = element_blank())+
   theme(plot.background = element_blank())+
   theme(panel.border = element_rect(colour = "black",linewidth=1),axis.line = element_blank())+
   labs(y="", x= "")+
   scale_x_continuous(expand = c(0,0),limits = c(0,190000),
                      breaks=c(0, 50000, 100000,150000),
                      labels=c("0","50,000", "100,000", "150,000"))+
   theme(legend.position="none",
         axis.text.x = element_text(size=15, colour = "black",angle = 90,vjust = 0.5,hjust = 1),
         axis.text.y = element_text(size=15, colour = "black",hjust = 0))

#ECs
egg_EC <- egg_args_sum_top30%>%filter(cat=="EC")
egg_EC$number <- paste0("EC:",egg_EC$number)
egg_EC$number <- factor(egg_EC$number,levels = rev(egg_EC$number))
egg_EC$text <- paste0("(n=",formattable::comma(egg_EC$n_gene,format = "d"),")")
egg_EC_plot <- ggplot(egg_EC, aes(x=n_gene, y=number)) +
   geom_bar(stat="identity", width=0.7,fill="#9A8F09")+
   geom_text(aes(label=text),hjust=-0.05,size=3)+
   theme_bw()+
   theme(panel.grid = element_blank())+
   theme(plot.background = element_blank())+
   theme(panel.border = element_rect(colour = "black",linewidth=1),axis.line = element_blank())+
   labs(y="", x= "")+
   scale_x_continuous(expand = c(0,0),limits = c(0,13000),
                      breaks=c(0,2500,5000,7500,10000),
                      labels=c("0","2,500","5,000","7,500","10,000"))+
   theme(legend.position="none",
         axis.text.x = element_text(size=15, colour = "black",angle = 90,vjust = 0.5,hjust = 1),
         axis.text.y = element_text(size=15, colour = "black",hjust = 0))

#KEGG_Module
egg_kegg <- egg_args_sum_top30%>%filter(cat=="KEGG_Module")
egg_kegg$number <- factor(egg_kegg$number,levels = rev(egg_kegg$number))
egg_kegg$text <- paste0("(n=",formattable::comma(egg_kegg$n_gene,format = "d"),")")
egg_kegg_plot <- ggplot(egg_kegg, aes(x=n_gene, y=number)) +
   geom_bar(stat="identity", width=0.7,fill="#9A4209")+
   geom_text(aes(label=text),hjust=-0.05,size=3)+
   theme_bw()+
   theme(panel.grid = element_blank())+
   theme(plot.background = element_blank())+
   theme(panel.border = element_rect(colour = "black",linewidth=1),axis.line = element_blank())+
   labs(y="", x= "")+
   scale_x_continuous(expand = c(0,0),limits = c(0,35000),
                      breaks=c(0, 10000,20000,30000),
                      labels=c("0","10,000", "20,000", "30,000"))+
   theme(legend.position="none",
         axis.text.x = element_text(size=15, colour = "black",angle = 90,vjust = 0.5,hjust = 1),
         axis.text.y = element_text(size=15, colour = "black",hjust = 0))


