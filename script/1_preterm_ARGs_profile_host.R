library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(coin)
library(stringr)
library(readxl)

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
                          DOL>360 & DOL<=540 ~ "18",
                          DOL>540 ~ "36"))%>%
   mutate(st_subject=paste0(Study,"_",SubjectID))
preterm_metadata_month <- preterm_metadata%>%filter(!is.na(DOL))%>%group_by(month)%>%dplyr::summarise(n=n())%>%mutate(aa=1)
preterm_metadata_month$ymax <- cumsum(preterm_metadata_month$n)
preterm_metadata_month$aa <- as.factor(preterm_metadata_month$aa)
preterm_metadata_month$month <- factor(preterm_metadata_month$month,levels = rev(preterm_metadata_month$month))

#Supplementary Figure 2a
preterm_month <- ggplot(data=preterm_metadata_month,aes(x=aa, y=n,fill=month)) +
   geom_bar(stat="identity",width = 0.7)+
   theme_bw()+
   scale_fill_manual(values = c("#8DD3C7", "#D0D046", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                                "#B3DE69", "#FCCDE5", "#BBEBF6", "#BC80BD"))+
   scale_x_discrete(expand = c(0.5,0))+
   scale_y_continuous(expand = c(0,0),limits = c(0,5683),breaks = c(0,preterm_metadata_month$ymax))+
   theme(panel.grid = element_blank(),axis.line.y = element_line(color = "black",linewidth = 0.5),
         panel.border = element_blank())+
   theme(panel.background = element_blank(),
         plot.background = element_blank(),
         legend.position = "none",
         strip.background = element_rect(colour="black", linewidth = 1),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

#types and subtypes
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
   arrange(desc(prevalence))
type_se <- data.frame("type"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in preterm_type$names) {
   a <- preterm_type%>%filter(names==t)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%setNames("type")
   type_se[nrow(type_se)+1,] <- c(t,std.error(a$type))
}
preterm_type_sum_se <- left_join(preterm_type_sum,type_se,by="type")
preterm_type_sum_se$se <- as.numeric(preterm_type_sum_se$se)
preterm_type_sum_abund <- preterm_type%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"type_total_abund"=rowSums(.[-1]))
preterm_type_sum_se$Mean_copies_per_cell[preterm_type_sum_se$type=="multidrug"]#3.654255
preterm_type_sum_se$Mean_copies_per_cell[preterm_type_sum_se$type=="multidrug"] <- 2.0
preterm_type_sum_se$type <- factor(preterm_type_sum_se$type,levels = c("multidrug","MLS","beta_lactam","polymyxin","bacitracin","aminoglycoside","tetracycline","trimethoprim",
                                                                       preterm_type_sum$type[9:28]))
#Figure 1b
preterm_type_p <- ggplot() +
   geom_bar(data=preterm_type_sum_se,aes(x=type, y=Mean_copies_per_cell),stat="identity",fill="#4994C4",color="black",width = 0.8)+
   geom_errorbar(data=preterm_type_sum_se,aes(x=type,ymin=Mean_copies_per_cell, ymax=Mean_copies_per_cell+se), width=0.5)+
   geom_line(data=preterm_type_sum_se,aes(x=type, y=prevalence/0.6,group=1),linetype = "dashed",color="#4994C4",size=0.5)+
   geom_point(data=preterm_type_sum_se,aes(x=type, y=prevalence/0.6),size=1.5)+
   geom_hline(yintercept=0.9/0.6,linetype="dashed",size=0.5,color="black")+
   scale_y_continuous(expand = c(0,0),name= "Mean ARG abundance (capc)",limits = c(0,2.2),
                      breaks = c(0,0.5,1.0,1.5,2.0),labels = c("0","0.5","1.0","1.5","3.7"),
                      sec.axis = sec_axis(~ . *0.6, breaks = c(0,0.25,0.5,0.75,0.9,1)))+
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

#type dynamics over days of life
preterm_type_2 <- preterm_type%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
preterm_type_2_meta <- left_join(preterm_metadata%>%dplyr::select(SequenceID,month),preterm_type_2,by="SequenceID")%>%filter(!is.na(month))
preterm_type_2_meta_g <- gather(preterm_type_2_meta,type,abundance,-month,-SequenceID)
preterm_type_2_meta_g <- left_join(preterm_type_2_meta_g,preterm_type_sum_abund,by=c("SequenceID"="sample"))
preterm_type_2_meta_g$prop <- preterm_type_2_meta_g$abundance/preterm_type_2_meta_g$type_total_abund
preterm_type_2_meta_g_sum <- preterm_type_2_meta_g%>%group_by(month,type)%>%dplyr::summarise(mean_prop=mean(prop),mean_abun=mean(abundance))
#selet types with prevalence >10% and average abundance > 0.1 capc
preterm_type_sum_sel <- preterm_type_sum%>%filter(prevalence>0.1)%>%filter(Mean_copies_per_cell>0.1)
preterm_type_2_meta_g_sum_sel <- preterm_type_2_meta_g_sum%>%filter(type%in%preterm_type_sum_sel$type)
preterm_type_2_meta_g_sum_sel$month <- factor(preterm_type_2_meta_g_sum_sel$month,levels = c(0.25,0.5,0.75,1,2,3,6,12,18,36))
preterm_type_2_meta_g_sum_sel_sum <- preterm_type_2_meta_g_sum_sel%>%group_by(type)%>%dplyr::summarise(mean=mean(mean_abun))%>%arrange(desc(mean))
preterm_type_2_meta_g_sum_sel$type <- factor(preterm_type_2_meta_g_sum_sel$type,levels = preterm_type_2_meta_g_sum_sel_sum$type)
#Figure 4d
dynamic_type <- ggplot(preterm_type_2_meta_g_sum_sel, aes(x=month, y=mean_abun, fill=type)) + 
   geom_bar(stat = "identity", width=0.8)+
   scale_fill_manual(values = c("#8DD3C7", "#D0D046", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                                "#B3DE69", "#FCCDE5", "#BBEBF6", "#BC80BD", "#CCEBC5","#FFED6F"))+
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

type_dynamic_pvalue <- data.frame("type"="x","p_val"=0,stringsAsFactors=F)[-1,]
col_list <- c("#8DD3C7", "#D0D046", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
              "#B3DE69", "#FCCDE5", "#BBEBF6", "#BC80BD", "#CCEBC5","#FFED6F")
count=0
#Supplementary Figure 8
for (v in levels(preterm_type_2_meta_g_sum_sel$type)) {
   # v="bacitracin"
   count=count+1
   print(v)
   db_1 <- preterm_type_2_meta_g%>%filter(type==v)
   db_1$month <- factor(db_1$month,levels = c("0.25","0.5","0.75","1","2","3","6","12","18","36"))
   p_val <- pvalue(kruskal_test(abundance ~ month, data=db_1))
   type_dynamic_pvalue[nrow(type_dynamic_pvalue)+1,] <- c(v,p_val)
   db_1 <- db_1%>%group_by(month)%>%dplyr::mutate(month_2=cur_group_id())
   p <- ggplot(db_1, aes(x=month_2, y=abundance)) +
      geom_smooth(method="loess", level=0.95,colour=col_list[count],linewidth=1)+
      labs(y=v)+
      scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10),
                         labels = c("1"="0.25","2"="0.5","3"="0.75","4"="1","5"="2",
                                    "6"="3","7"="6","8"="12","9"="18","10"="36"))+
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
preterm_subtype_richness <- preterm_subtype%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"subtype_richness"=rowSums(.[-1]>0))
preterm_subtype_abund <- preterm_subtype%>%column_to_rownames("names")%>%
   t()%>%as.data.frame()%>%rownames_to_column("sample")%>%dplyr::summarise("sample"=sample,"subtype_total_abund"=rowSums(.[-1]))
preterm_subtype_sum <- preterm_subtype%>%
   dplyr::summarise("subtype"=names,
                    "n"=rowSums(.[-1]>0),
                    "prevalence"=rowSums(.[-1]>0)/(ncol(.)-1),
                    "Mean_copies_per_cell"=rowMeans(.[-1]))%>%
   arrange(desc(prevalence))
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
#plot the subtype with a prevalence >10%
preterm_subtype_sum_se_sel <- preterm_subtype_sum_se%>%filter(prevalence>0.10)
preterm_subtype_sum_se_sel_sum <- preterm_subtype_sum_se_sel%>%group_by(type)%>%dplyr::summarise(n=n())%>%dplyr::arrange(desc(n))
preterm_subtype_sum_se_nonsel <- preterm_subtype_sum_se%>%filter(prevalence<=0.10)%>%mutate(type="ungrouped")
preterm_subtype_sum_se_4 <- rbind(preterm_subtype_sum_se_sel,preterm_subtype_sum_se_nonsel)
preterm_subtype_sum_se_4$type <- factor(preterm_subtype_sum_se_4$type,levels = c(preterm_subtype_sum_se_sel_sum$type,"ungrouped"))
#Figure 1c
preterm_subtype_p <- ggplot(preterm_subtype_sum_se_4, aes(x=Mean_copies_per_cell, y=prevalence)) +
   geom_point(aes(fill=type),size=2.5,shape=21,color=alpha("black",0.5))+
   geom_hline(yintercept=0.1,linetype="dashed",size=0.5,color="black")+
   labs(x ="Mean ARG subtype abundance (capc)", y="Prevalence")+
   theme_bw()+
   scale_y_continuous(breaks = c(0,0.10,0.25,0.5,0.75,1))+
   scale_fill_manual(values = c("#43760A", "#BEAED4", "#934B07", "#800404", "#386CB0", "#F3A2CD", 
                                "#DA732D", "#F54646", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                                "#7FC97F", "#E6AB02", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
                                "#BC80BD", "#CCEBC5", "#FFED6F", "#BBBABA"))+
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

#amr results
amr_profile_coverm <- read.csv("~/Desktop/work/preterm/analysis/revision/amr_preterm_resistome_coverm_bins_gtdb_stat_plasmid.csv")
amr_profile_coverm_zero <- amr_profile_coverm%>%filter(rpkm==0)
#ARGs
amr_profile_coverm_aro <- amr_profile_coverm%>%group_by(sample,Element_symbol)%>%dplyr::summarise(rpkm_total=sum(rpkm))
amr_profile_coverm_aro_richness <- amr_profile_coverm_aro%>%group_by(sample)%>%dplyr::summarise(n=n())
amr_profile_coverm_aro_sum <- amr_profile_coverm_aro%>%
   group_by(Element_symbol)%>%
   dplyr::summarise("n"=n(),
                    "prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(rpkm_total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))
amr_profile_coverm_aro_spread <- amr_profile_coverm_aro%>%spread(.,sample,rpkm_total)
amr_profile_coverm_aro_spread[is.na(amr_profile_coverm_aro_spread)] <- 0
amr_gene_se <- data.frame("Element_symbol"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in amr_profile_coverm_aro_spread$Element_symbol) {
   a <- amr_profile_coverm_aro_spread%>%filter(Element_symbol==t)%>%column_to_rownames("Element_symbol")%>%t()%>%as.data.frame()%>%setNames("Element_symbol")
   amr_gene_se[nrow(amr_gene_se)+1,] <- c(t,std.error(a$Element_symbol))
}
amr_profile_coverm_aro_sum_se <- left_join(amr_profile_coverm_aro_sum,amr_gene_se,by="Element_symbol")
amr_profile_coverm_aro_sum_se$se <- as.numeric(amr_profile_coverm_aro_sum_se$se)
amr_profile_coverm_gene_2_top <- amr_profile_coverm_aro_sum_se%>%filter(prevalence>0.1)
amr_profile_coverm_gene_2_top$Element_symbol <- factor(amr_profile_coverm_gene_2_top$Element_symbol,levels = rev(amr_profile_coverm_gene_2_top$Element_symbol))
amr_profile_coverm_gene_2_top$label <- paste0("(n=",formattable::comma(amr_profile_coverm_gene_2_top$n,digits=0),")")
#Supplementary Figure 2d
amr_profile_coverm_gene_2_top_p <- ggplot() +
   geom_bar(data=amr_profile_coverm_gene_2_top, aes(x=Element_symbol, y=Mean_rpkm),stat="identity",fill="#C48849",color="black",width = 0.8)+
   geom_errorbar(data=amr_profile_coverm_gene_2_top,aes(x=Element_symbol,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=amr_profile_coverm_gene_2_top, aes(x=Element_symbol, y=prevalence*100000,group=1),linetype = "dashed",color="#C43A0A",size=1)+
   geom_point(data=amr_profile_coverm_gene_2_top, aes(x=Element_symbol, y=prevalence*100000),size=1.5)+
   geom_text(data=amr_profile_coverm_gene_2_top,aes(x=Element_symbol, y=57000,label=label),hjust=1,vjust=0.5,size=4)+
   scale_y_continuous(expand = c(0,0),name= "Mean ARG family abundance (rpkm)",limits = c(0,60000),
                      breaks = c(0,20000,40000,60000),labels = c("0","20,000","40,000","60,000"),
                      sec.axis = sec_axis(~ . /100000))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(legend.position = "none",
         legend.text = element_text(face = "italic"),
         plot.background = element_blank(),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         plot.margin = unit(c(0,1,0,0), "cm"),
         axis.text.y = element_text(size=12, colour = "black",face = "italic"),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())+
   coord_flip()

#compare ARG subtypes and ARGs from amr
amr_profile_coverm_aro_sum_modifiy <- read_excel("data/amr_profile_coverm_aro_sum.xlsx", sheet = "Sheet1")
amr_profile_coverm_aro_sum_modifiy$Element_symbol_2 <- tolower(amr_profile_coverm_aro_sum_modifiy$Element_symbol_2)
preterm_subtype_sum$gene <- sapply(strsplit(preterm_subtype_sum$subtype, split='__', fixed=TRUE), function(x)(x[2]))
preterm_subtype_sum$gene_2 <- tolower(preterm_subtype_sum$gene)
preterm_subtype_sum_aro_2 <- amr_profile_coverm_aro_sum_modifiy%>%filter(Element_symbol_2%in%preterm_subtype_sum$gene_2)%>%
   left_join(.,preterm_subtype_sum%>%dplyr::select(subtype,prevalence,gene,gene_2),by=c("Element_symbol_2"="gene_2"))
cor.test(preterm_subtype_sum_aro_2$prevalence.y,preterm_subtype_sum_aro_2$prevalence.x)#0.7680642, p-value < 2.2e-16
#Supplementary Figure 2c
corr_subtype_amr_gene <- ggplot(preterm_subtype_sum_aro_2, aes(x=prevalence.y, y=prevalence.x)) +
   geom_point(shape=21,color=alpha("black",0.5),fill=alpha("black",0.5),size=2)+
   geom_smooth(method=lm, formula =  y ~ poly(x,1), level=0,colour="#DA732D",linewidth=2) +
   labs(x ="Prevalence of ARG subtypes", y="Prevalence of assembled ARGs")+
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
   cols[which(cols == "")] <- NA
   cols <- tibble::as.tibble(cols)
   m <- dim(cols)[2]
   names(cols) <- paste(into_prefix, 1:m, sep = "_")
   return(cols)
}
#summary for drug
amr_profile_coverm_2 <- amr_profile_coverm %>%
   dplyr::select(run_contig,sample,rpkm,Subclass)%>%
   dplyr::bind_cols(split_into_multiple(.$Subclass, "/", "cat"))%>%
   gather(.,cat,drug,-c(run_contig,sample,rpkm,Subclass))%>%
   filter(!is.na(drug))%>%
   dplyr::select(-cat)
amr_profile_coverm_2_sum <- amr_profile_coverm_2%>%group_by(sample,drug)%>%dplyr::summarise(rpkm_total=sum(rpkm))
preterm_amr_drug_richness <- amr_profile_coverm_2_sum%>%group_by(sample)%>%dplyr::summarise(drug_richness=n())
amr_profile_coverm_3 <- amr_profile_coverm_2_sum%>%
   group_by(drug)%>%
   dplyr::summarise("n"=n(),
                    "prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(rpkm_total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))
amr_profile_coverm_2_sum_1 <- amr_profile_coverm_2%>%
   group_by(sample,drug)%>%dplyr::summarise(rpkm_total=sum(rpkm))%>%spread(.,sample,rpkm_total)
which(colSums(amr_profile_coverm_2_sum_1[-1])==0)#0
amr_profile_coverm_2_sum_1[is.na(amr_profile_coverm_2_sum_1)] <- 0
amr_drug_se <- data.frame("drug"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in amr_profile_coverm_2_sum_1$drug) {
   a <- amr_profile_coverm_2_sum_1%>%filter(drug==t)%>%column_to_rownames("drug")%>%t()%>%as.data.frame()%>%setNames("drug")
   amr_drug_se[nrow(amr_drug_se)+1,] <- c(t,std.error(a$drug))
}
amr_profile_coverm_3_se <- left_join(amr_profile_coverm_3,amr_drug_se,by="drug")
amr_profile_coverm_3_se$se <- as.numeric(amr_profile_coverm_3_se$se)
amr_profile_coverm_3_se$drug <- tolower(amr_profile_coverm_3_se$drug)
amr_profile_coverm_3_se$drug[amr_profile_coverm_3_se$drug=="streptogramin b"] <- "streptogramin B"
amr_profile_coverm_3_se$drug[amr_profile_coverm_3_se$drug=="streptogramin a"] <- "streptogramin A"
amr_profile_coverm_3_se$drug[amr_profile_coverm_3_se$drug=="virginiamycin m"] <- "virginiamycin M"
amr_profile_coverm_3_se$drug <- factor(amr_profile_coverm_3_se$drug,levels = amr_profile_coverm_3_se$drug)
#Supplementary Figure 2e
amr_drug_class_p <- ggplot() +
   geom_bar(data=amr_profile_coverm_3_se, aes(x=drug, y=Mean_rpkm),stat="identity",fill="#C48849",color="black",width = 0.8)+
   geom_errorbar(data=amr_profile_coverm_3_se,aes(x=drug,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=amr_profile_coverm_3_se, aes(x=drug, y=prevalence/0.0000055,group=1),linetype = "dashed",color="#C48849",size=0.5)+
   geom_point(data=amr_profile_coverm_3_se, aes(x=drug, y=prevalence/0.0000055),size=1.5)+
   geom_hline(yintercept=0.5/0.0000055,linetype="dashed",size=0.5,color="black")+
   scale_y_continuous(expand = c(0,0),name= "Mean abundance (rpkm)",limits = c(0,150000),
                      breaks = c(0,50000,100000,150000),labels = c("0","50,000","100,000","150,000"),
                      sec.axis = sec_axis(~ . *0.0000055))+
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
         axis.text.x = element_text(size=10, colour="black",angle = 90,vjust = 0.4,hjust = 1),
         axis.title.y = element_blank(),
         axis.title.x = element_blank())

#amr drug class dynamics
amr_profile_coverm_2_sum_2 <- amr_profile_coverm_2_sum_1%>%column_to_rownames("drug")%>%t()%>%as.data.frame()%>%rownames_to_column("SequenceID")
preterm_metadata_sel <- preterm_metadata%>%dplyr::select(SequenceID,month)%>%filter(SequenceID%in%amr_profile_coverm_2_sum_2$SequenceID)
amr_profile_coverm_2_sum_1_meta <- left_join(preterm_metadata_sel,amr_profile_coverm_2_sum_2,by="SequenceID")%>%filter(!is.na(month))
amr_profile_coverm_2_sum_1_meta_g <- gather(amr_profile_coverm_2_sum_1_meta,drug,abundance,-month,-SequenceID)
amr_drug_sum_abund <- data.frame("SequenceID"=amr_profile_coverm_2_sum_2$SequenceID,"drug_total_abund"=rowSums(amr_profile_coverm_2_sum_2[-1]))
amr_profile_coverm_2_sum_1_meta_g <- left_join(amr_profile_coverm_2_sum_1_meta_g,amr_drug_sum_abund,by="SequenceID")
amr_profile_coverm_2_sum_1_meta_g$prop <- amr_profile_coverm_2_sum_1_meta_g$abundance/amr_profile_coverm_2_sum_1_meta_g$drug_total_abund
amr_profile_coverm_2_sum_1_meta_g_sum <- amr_profile_coverm_2_sum_1_meta_g%>%group_by(month,drug)%>%
   dplyr::summarise(mean_prop=mean(prop),mean_abun=mean(abundance))
amr_profile_coverm_2_sum_1_meta_g_sum_s <- spread(amr_profile_coverm_2_sum_1_meta_g_sum%>%dplyr::select(-mean_abun),month,mean_prop)
amr_profile_coverm_3_sel <- amr_profile_coverm_3%>%filter(prevalence>0.5)
amr_profile_coverm_2_sum_1_meta_g_sum_sel <- amr_profile_coverm_2_sum_1_meta_g_sum%>%filter(drug%in%amr_profile_coverm_3_sel$drug)
amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug <- tolower(amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug)
amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug[amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug=="streptogramin b"] <- "streptogramin B"
amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug[amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug=="streptogramin a"] <- "streptogramin A"
amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug[amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug=="virginiamycin m"] <- "virginiamycin M"
amr_profile_coverm_2_sum_1_meta_g_sum_sel$month <- factor(amr_profile_coverm_2_sum_1_meta_g_sum_sel$month,levels = c(0.25,0.5,0.75,1,2,3,6,12,18,36))
amr_profile_coverm_2_sum_1_meta_g_sum_sel_sum <- amr_profile_coverm_2_sum_1_meta_g_sum_sel%>%group_by(drug)%>%dplyr::summarise(mean=mean(mean_abun))%>%arrange(desc(mean))
amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug <- factor(amr_profile_coverm_2_sum_1_meta_g_sum_sel$drug,levels = amr_profile_coverm_2_sum_1_meta_g_sum_sel_sum$drug)
#Supplementary Figure 9a
dynamic_amr <- ggplot(amr_profile_coverm_2_sum_1_meta_g_sum_sel, aes(x=month, y=mean_abun, fill=drug)) + 
   geom_bar(stat = "identity", width=0.8)+
   scale_fill_manual(values = c("#8DD3C7", "#D0D046", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                                "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5","#D05046",
                                "#F0A8A3"))+
   scale_y_continuous(expand = c(0, 0),limits = c(0,1250000),breaks = c(0,250000,500000,750000,1000000,1250000),
                      labels = c("0","250,000","500,000","750,000","1,000,000","1,250,000"))+
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

#unique ARGs in different studies
preterm_subtype_sum_only1 <- preterm_subtype_sum%>%filter(n==1)
preterm_subtype_only1 <- preterm_subtype%>%filter(names%in%preterm_subtype_sum_only1$subtype)%>%
   dplyr::select(names,which(colSums(.[-1])!=0)+1)%>%column_to_rownames("names")%>%t()%>%as.data.frame()%>%
   rownames_to_column("SequenceID")%>%dplyr::summarise("SequenceID"=SequenceID,"n"=rowSums(.[-1]>0))%>%
   left_join(.,preterm_metadata%>%dplyr::select(SequenceID,Study,clean_reads),by="SequenceID")
preterm_subtype_only1$clean_reads_2 <- preterm_subtype_only1$clean_reads/1000000
cor.test(preterm_subtype_only1$clean_reads,preterm_subtype_only1$n)#r = 0.06823991 , p-value = 0.4913
preterm_subtype_only1_sum <- preterm_subtype_only1%>%group_by(Study)%>%dplyr::summarise(total=sum(n))%>%arrange(desc(total))
preterm_subtype_only1_sum$total[1] <- 30
#different from ARG subtypes, ARG gene is obtained from amr, not coverm
amr_profile_coverm_aro_sum_only1 <- amr_profile_coverm_aro_sum%>%filter(n==1)
preterm_amr_only1 <- amr_profile_coverm%>%filter(Element_symbol%in%amr_profile_coverm_aro_sum_only1$Element_symbol)%>%
   dplyr::group_by(sample)%>%dplyr::summarise(n=n())%>%
   left_join(.,preterm_metadata%>%dplyr::select(SequenceID,Study,clean_reads),by=c("sample"="SequenceID"))
preterm_amr_only1$clean_reads_2 <- preterm_amr_only1$clean_reads/1000000
cor.test(preterm_amr_only1$clean_reads,preterm_amr_only1$n)#r = 0.1053792 , p-value = p-value = 0.3818

preterm_amr_only1_sum <- preterm_amr_only1%>%group_by(Study)%>%dplyr::summarise(total=sum(n))%>%arrange(desc(total))
preterm_amr_only1_sum_subtype <- rbind(preterm_amr_only1_sum%>%mutate(group="Assembled ARGs"),preterm_subtype_only1_sum%>%mutate(group="ARG subtypes"))
preterm_amr_only1_sum_subtype$Study <- factor(preterm_amr_only1_sum_subtype$Study,levels = rev(c(as.character(preterm_subtype_only1_sum$Study),
                                                                                                 "RoseG_2017","RahmanSF_2018")))
preterm_amr_only1_sum_subtype$group <- factor(preterm_amr_only1_sum_subtype$group,levels = c("ARG subtypes","Assembled ARGs"))
#Supplementary Figure 2f
subtype_amr_only_p <- ggplot(data=preterm_amr_only1_sum_subtype, aes(x=total, y=Study,fill=group)) +
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

#Supplementary Figure 2g
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
#Supplementary Figure 2g
amr_only_p_2 <- ggplot(preterm_amr_only1, aes(x=clean_reads_2, y=n)) +
   geom_point(size=1,color="#C48849")+
   geom_smooth(method=lm, formula =  y ~ poly(x,1), level=0,colour="grey",linewidth=1)+
   labs(x ="Number of clean reads (million)", y="Number of aARG")+
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

#correlation between arg subytpes and aro
amr_gene_richness_2 <- amr_profile_coverm%>%group_by(sample,Element_symbol)%>%dplyr::summarise(total_rpkm=sum(rpkm))%>%
   ungroup()%>%group_by(sample)%>%
   dplyr::summarise(gene_richness=sum(total_rpkm>0))
preterm_subtype_richness_amr <- left_join(preterm_subtype_richness,amr_gene_richness_2,by="sample")
preterm_subtype_richness_amr$gene_richness[is.na(preterm_subtype_richness_amr$gene_richness)] <- 0
cor.test(preterm_subtype_richness_amr$gene_richness,preterm_subtype_richness_amr$subtype_richness)#0.4285148, p-value < 2.2e-16
#Supplementary Figure 2b
corr_subtype_amr <- ggplot(preterm_subtype_richness_amr, aes(x=subtype_richness, y=gene_richness)) +
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

density_subtypes <- ggplot(preterm_subtype_richness_amr, aes(x=subtype_richness)) +
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

density_amr_gene <- ggplot(preterm_subtype_richness_amr, aes(x=gene_richness)) +
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
amr_profile_stats <- read.csv("data/amr_preterm_resistome_coverm_bins_gtdb_stat_plasmid.csv")
preterm_metadata <- read.csv("data/metadata_cb_infant_preterm.csv")
amr_profile_stats_nobin <- amr_profile_stats%>%filter(is.na(bin))
amr_profile_stats$domain <- sapply(strsplit(amr_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[1]))
amr_profile_stats$phylum <- sapply(strsplit(amr_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[2]))
amr_profile_stats$class <- sapply(strsplit(amr_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[3]))
amr_profile_stats$order <- sapply(strsplit(amr_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[4]))
amr_profile_stats$family <- sapply(strsplit(amr_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[5]))
amr_profile_stats$genus <- sapply(strsplit(amr_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[6]))
amr_profile_stats$species <- sapply(strsplit(amr_profile_stats$GTDB_classification, split=';', fixed=TRUE), function(x)(x[7]))
amr_profile_stats$domain[is.na(amr_profile_stats$domain)] <- "unknown"
amr_profile_stats$phylum[is.na(amr_profile_stats$phylum)] <- "unknown"
amr_profile_stats$class[is.na(amr_profile_stats$class)] <- "unknown"
amr_profile_stats$order[is.na(amr_profile_stats$order)] <- "unknown"
amr_profile_stats$family[is.na(amr_profile_stats$family)] <- "unknown"
amr_profile_stats$genus[is.na(amr_profile_stats$genus)] <- "unknown"
amr_profile_stats$species[amr_profile_stats$species=="s__"] <- "unknown"
amr_profile_stats$species[is.na(amr_profile_stats$species)] <- "unknown"
amr_profile_stats_species <- amr_profile_stats%>%filter(species!="unknown")

amr_profile_stats_domain <- amr_profile_stats%>%group_by(domain)%>%dplyr::summarise(n=n())
tax_number <- data.frame("tax"="x","n_tax"=0,"n_args"=0,stringsAsFactors = F)[-1,]
tax_number[1,] <- c("Phylum",length(table(amr_profile_stats$phylum))-1,sum(table(amr_profile_stats$phylum))-nrow(amr_profile_stats%>%filter(phylum=="unknown")))
tax_number[2,] <- c("Class",length(table(amr_profile_stats$class))-1,sum(table(amr_profile_stats$class))-nrow(amr_profile_stats%>%filter(class=="unknown")))
tax_number[3,] <- c("Order",length(table(amr_profile_stats$order))-1,sum(table(amr_profile_stats$order))-nrow(amr_profile_stats%>%filter(order=="unknown")))
tax_number[4,] <- c("Family",length(table(amr_profile_stats$family))-1,sum(table(amr_profile_stats$family))-nrow(amr_profile_stats%>%filter(family=="unknown")))
tax_number[5,] <- c("Genus",length(table(amr_profile_stats$genus))-1,sum(table(amr_profile_stats$genus))-nrow(amr_profile_stats%>%filter(genus=="unknown")))
tax_number[6,] <- c("Species",length(table(amr_profile_stats$species))-1,sum(table(amr_profile_stats$species))-nrow(amr_profile_stats%>%filter(species=="unknown")))
tax_number$n_tax <- as.numeric(tax_number$n_tax)
tax_number$n_args <- as.numeric(tax_number$n_args)
tax_number$tax <- factor(tax_number$tax,levels = c("Phylum", "Class","Order","Family","Genus","Species"))
#Supplementary Figure 3b
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
   #t="phylum"
   a <- amr_profile_stats%>%dplyr::group_by(get(t))%>%dplyr::summarise(n=n())%>%arrange(desc(n))
   print(sum(a$n))
   b <- a[c(1:11),] 
   b$propor <- b$n/sum(a$n)
   b[nrow(b)+1,] <- list("Others",nrow(amr_profile_stats)-sum(b$n,na.rm = T),1-sum(b$propor,na.rm = T))
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
#Supplementary Figure 3a
amr_taxa_all <- ggplot(tax_diver_2, aes(x=taxonomy, y=n, fill=new_t)) + 
   geom_bar(stat = "identity", width=0.8)+
   scale_fill_manual(values = c("#CACACA","#8F8F8F","#6AABF5","#BEBADA","#FB8072","#80B1D3","#FDB462",
                                "#B3DE69","#FCCDE5","#CABB57","#BC80BD","#8DD3C7"))+
   theme_bw() +
   labs(x = element_blank(), y="Number of ARGs")+
   theme(axis.title.y = element_text(size = 16))+
   theme(panel.border = element_blank(), axis.line = element_line()) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #remove the grid
   scale_y_continuous(expand = c(0,0),breaks=c(0, 20000, 40000,60000,80000),
                      labels=c("0","20,000", "40,000", "60,000","80,000"))+
   theme(legend.position="none",
         aspect.ratio = 1.3/1,
         plot.background = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.text.x = element_text(size=12, colour="black",angle = 45,vjust = 0.7))

#only species
std.error <- function(x) sd(x)/sqrt(length(x))
amr_profile_stats_species <- amr_profile_stats%>%dplyr::select(sample,species,rpkm)
amr_profile_stats_species$species <- gsub("s__","",amr_profile_stats_species$species,fixed = T)
amr_profile_stats_species_sum <- amr_profile_stats_species%>%group_by(sample,species)%>%dplyr::summarise(total=sum(rpkm))
amr_profile_stats_species_sum_1 <- amr_profile_stats_species_sum%>%
   group_by(species)%>%
   dplyr::summarise("prevalence"=n()/length(unique(.$sample)),
                    "Mean_rpkm"=sum(total)/length(unique(.$sample)))%>%
   arrange(desc(prevalence))
amr_profile_stats_species_sum_s <- spread(amr_profile_stats_species_sum,sample,total)
amr_profile_stats_species_sum_s[is.na(amr_profile_stats_species_sum_s)] <- 0
amr_species_se <- data.frame("species"="x","se"=0,stringsAsFactors = F)[-1,]
for (t in amr_profile_stats_species_sum_s$species) {
   a <- amr_profile_stats_species_sum_s%>%filter(species==t)%>%column_to_rownames("species")%>%t()%>%as.data.frame()%>%setNames("species")
   amr_species_se[nrow(amr_species_se)+1,] <- c(t,std.error(a$species))
}
amr_profile_stats_species_sum_2 <- left_join(amr_profile_stats_species_sum_1,amr_species_se,by="species")
amr_profile_stats_species_sum_2$se <- as.numeric(amr_profile_stats_species_sum_2$se)
amr_profile_stats_species_sum_2 <- amr_profile_stats_species_sum_2%>%filter(species!="unknown")
amr_profile_stats_species_sum_3 <- amr_profile_stats_species_sum_2%>%dplyr::filter(prevalence>0.01)
amr_profile_stats_species_2 <- amr_profile_stats_species%>%group_by(species)%>%dplyr::summarise(n_args=n())%>%
   filter(species%in%amr_profile_stats_species_sum_3$species)%>%arrange(desc(n_args))
amr_profile_stats_species_sum_3 <- left_join(amr_profile_stats_species_sum_3,amr_profile_stats_species_2,by="species")
amr_profile_stats_species_sum_3$species <- factor(amr_profile_stats_species_sum_3$species,levels = rev(amr_profile_stats_species_sum_3$species))
amr_profile_stats_species_sum_3$label <- paste0("(n=",formattable::comma(amr_profile_stats_species_sum_3$n_args,digits=0),")")
amr_profile_stats_species_sum_3$label[amr_profile_stats_species_sum_3$species=="Klebsiella pneumoniae"]#"(n=4,346)"
amr_profile_stats_species_sum_3$label[amr_profile_stats_species_sum_3$species=="Klebsiella pneumoniae"] <- NA
amr_profile_stats_species_sum_3$label[amr_profile_stats_species_sum_3$species=="Enterococcus faecalis"]#"(n=4,813)"
amr_profile_stats_species_sum_3$label[amr_profile_stats_species_sum_3$species=="Enterococcus faecalis"] <- NA
#Figure 1d
amr_taxa_species <- ggplot() +
   geom_bar(data=amr_profile_stats_species_sum_3, aes(x=species, y=Mean_rpkm),stat="identity",fill="#C48849",color="black",width = 0.8)+
   geom_errorbar(data=amr_profile_stats_species_sum_3,aes(x=species,ymin=Mean_rpkm, ymax=Mean_rpkm+se), width=0.5)+
   geom_line(data=amr_profile_stats_species_sum_3, aes(x=species, y=prevalence*150000,group=1),linetype = "dashed",color="#C43A0A",size=1)+
   geom_point(data=amr_profile_stats_species_sum_3, aes(x=species, y=prevalence*150000),size=1.5)+
   geom_text(data=amr_profile_stats_species_sum_3,aes(x=species, y=Mean_rpkm+se,label=label),hjust=-0.1,vjust=0.5,size=3)+
   scale_y_continuous(expand = c(0,0),name= "Mean abundance (rpkm)",limits = c(0,80000),
                      breaks = c(0,20000,40000,60000,80000),labels = c("0","20,000","40,000","60,000","80,000"),
                      sec.axis = sec_axis(~ . /150000, breaks = c(0,0.15,0.30,0.45)))+
   theme_bw()+
   theme(panel.grid = element_blank(),axis.line = element_blank(),
         panel.border = element_rect(colour = "black",linewidth = 1))+
   theme(legend.position = "none",
         legend.text = element_text(face = "italic"),
         plot.background = element_blank(),
         panel.background = element_blank(),
         strip.text = element_text(color="black",size=10),
         panel.spacing = unit(5,"mm"),
         plot.margin = unit(c(0,1,0,0), "cm"),
         axis.text.y = element_text(size=11, colour = "black",face="italic"),
         axis.text.x = element_text(size=12, colour = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=12, colour = "black")
   )+
   coord_flip()

#################################################################################################################################################
#plasmid
#################################################################################################################################################
amr_profile_stats$bin[is.na(amr_profile_stats$bin)] <- "nobin"
amr_profile_stats_unknown_taxa <- amr_profile_stats%>%filter(bin=="nobin")
amr_profile_stats_unknown_taxa_sum <- amr_profile_stats_unknown_taxa%>%
   group_by(plasmid)%>%dplyr::summarise(n=n())%>%
   mutate(Propotion=n/sum(n))%>%
   mutate(group="Unknown")
amr_profile_stats_known_taxa <- amr_profile_stats%>%filter(bin!="nobin")
amr_profile_stats_known_taxa_sum <- amr_profile_stats_known_taxa%>%
   group_by(plasmid)%>%dplyr::summarise(n=n())%>%
   mutate(Propotion=n/sum(n))%>%
   mutate(group="Known")
amr_profile_stats_unknown_known <- rbind(amr_profile_stats_unknown_taxa_sum,amr_profile_stats_known_taxa_sum)
plas_preterm_sum <- amr_profile_stats%>%group_by(plasmid)%>%dplyr::summarise(n=n())
amr_profile_stats_plasmid <- amr_profile_stats%>%filter(plasmid=="Yes")
amr_profile_stats_plasmid_nospecies <- amr_profile_stats%>%filter(plasmid=="Yes")%>%filter(domain=="d__Bacteria")%>%filter(species=="unknown")
#Supplementary Figure 4a
pie(plas_preterm_sum$n,
    labels = c("Chromosome \n(n=53,791)","Plasmid \n(n=26,735)"),
    col = c("#A8DFFC","#F7A48F"),
    init.angle=360, radius = 1,cex=1)

plas_sum_list <- list()
for (v in c("phylum","class","order","family","genus","species")) {
   # v="genus"
   plas_sum <- data.frame("taxa"="x","pred_n"=0,"prediction"="x","args_n"=0,"args_n_chrom"=0,"args_n_plasmid"=0,stringsAsFactors = F)[-1,]
   count=0
   taxa_list <- amr_profile_stats%>%filter(get(v)!="unknown")%>%distinct(get(v),.keep_all = F)%>%setNames("taxa")
   for (i in taxa_list$taxa) {
      # i="g__Klebsiella"
      count=count+1
      print(count)
      a <- amr_profile_stats%>%filter(get(v)==i)
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
plas_sum_species <- plas_sum_list$species
plas_sum_species_chom <- plas_sum_species%>%filter(prediction=="No")%>%arrange(desc(args_n))
plas_sum_species_chom_sel <- plas_sum_species_chom%>%filter(args_n>1)
plas_sum_species_plasmid <- plas_sum_species%>%filter(prediction=="Yes")%>%arrange(desc(args_n))%>%filter(args_n>1)
plas_sum_species_both <- plas_sum_species%>%filter(pred_n==2)

#############################################################################################################################
#cluster of arg orfs from plasmid
#############################################################################################################################
amr_profile_stats_plasmid_yes_orf <- amr_profile_stats%>%filter(bin!="nobin")%>%filter(plasmid=="Yes")
amr_profile_stats_plasmid_yes_orf$sample_orf <- paste0(amr_profile_stats_plasmid_yes_orf$sample,"_",amr_profile_stats_plasmid_yes_orf$Gene)
amr_profile_stats_run <- amr_profile_stats%>%distinct(sample,.keep_all = T)%>%mutate(id=paste0(sample,"_gene.fna"))%>%dplyr::select(id)
amr_cluster_plas <- read.delim("data/amr_plasmid_orf_nr_cluster.txt")
amr_cluster_plas$clstr_iden <- gsub("%","",amr_cluster_plas$clstr_iden,fixed = T)
amr_cluster_plas$clstr_cov <- gsub("%","",amr_cluster_plas$clstr_cov,fixed = T)
amr_cluster_plas$clstr_iden <- as.numeric(amr_cluster_plas$clstr_iden)
amr_cluster_plas$clstr_cov <- as.numeric(amr_cluster_plas$clstr_cov)
amr_profile_stats_plas_cluster <- left_join(amr_profile_stats_plasmid_yes_orf,amr_cluster_plas%>%dplyr::select(id,clstr),by=c("sample_orf"="id"))
table(amr_profile_stats_plas_cluster$species)
for (v in c("phylum","class","order","family","genus","species")) {
   print(v)
   cluster_sum_1 <- amr_profile_stats_plas_cluster%>%group_by(clstr)%>%dplyr::summarise(n_taxa_all=length(unique(get(v))),n_args=n())
   cluster_sum_2 <- amr_profile_stats_plas_cluster%>%filter(get(v)!="unknown")%>%group_by(clstr)%>%dplyr::summarise(n_taxa_sel=length(unique(get(v))),n_args_2=n())
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
cluster_sum_sel_species <- amr_profile_stats_plas_cluster%>%filter(clstr%in%cluster_sum_sel$clstr)
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
cluster_sum_sel_species_s <- spread(cluster_sum_sel_species_sum_sel,species,n_args_2)%>%column_to_rownames("clstr")
cluster_sum_sel_species_s[is.na(cluster_sum_sel_species_s)] <- 0
col_fun_type <- colorRamp2(c(0,3,10,20,30,50,100),c("white","#A8C6EC","#7FA4D3","#5A87BF","#3768A4","#184B8A","#073876"))
#Supplementary Figure 4d
plasmid_cluster <- Heatmap(as.matrix(cluster_sum_sel_species_s),cluster_rows=T,cluster_columns=T,col = col_fun_type,
                           show_row_names = T,show_column_names = T,
                           border = T,row_names_gp = gpar(fontsize =10),
                           rect_gp = gpar(col = "grey", lwd = 0.7),
                           column_names_gp = gpar(fontface="italic",fontsize=10),
                           border_gp=gpar(col = "black", lwd = 1),column_gap = unit(2, "mm"),
                           width = unit(35, "cm"), height = unit(8, "cm"),
                           show_heatmap_legend = TRUE,
                           heatmap_legend_param = list(title = "Number of ARG ORFs",direction = "horizontal", 
                                                       legend_height = unit(0.3, "cm"),legend_width = unit(4, "cm"),
                                                       at = c(0,5,20,30,50,100), labels = c("0","5","20","30","50",">100"))
)




