








site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list = c("reshape2","ggplot2","devtools","bindrcpp",
				"ggthemes","agricolae","dplyr")
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

package_list = c("digest")
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		source("https://bioconductor.org/biocLite.R")
		biocLite(p)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}

package_list = c("kassambara/ggpubr")
for(p in package_list){
	q=unlist(strsplit(p,split = "/"))[2]
	if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install_github(p)
		suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}

main_theme = theme(panel.background=element_blank(), panel.grid=element_blank(),
	axis.line.x=element_line(size=.5, colour="black"), axis.line.y=element_line(size=.5, colour="black"),
	axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=),
	legend.position="right", legend.background=element_blank(), legend.key=element_blank(), legend.text= element_text(size=),
	text=element_text(family="sans", size=))



alpha = read.table("result/faprotax/element_tab.txt", header=T, row.names=1, sep="\t", comment.char="") 

if (TRUE){
	alpha = as.data.frame(t(alpha))
}

if (FALSE){
	alpha = as.data.frame(alpha/rowSums(alpha,na=T) * 100) # normalization to 1000
}else{
  alpha=alpha/100
}
design = read.table("/mnt/bai/yongxin/rice/xianGeng/doc/design.txt", header=T, row.names=1, sep="\t")
design$group=design$groupID

if (TRUE){
	sub_design = subset(design, group %in% c("HTEJ","HIND","LTEJ","LIND"))
	sub_design$group  = factor(sub_design$group, levels=c("HTEJ","HIND","LTEJ","LIND"))
}else{
	sub_design = design
}

idx = rownames(sub_design) %in% rownames(alpha)
sub_design=sub_design[idx,]


list = read.table("result/compare_far/IND.list", header=F,  sep="\t")
sub_alpha=alpha[rownames(sub_design), as.vector(list$V1)]
library(dplyr)
sub_alpha$sampleID=rownames(sub_alpha) 
melted = melt(sub_alpha, id.vars = "sampleID")

melted_all = merge(melted, sub_design[,c("subspecies","soiltype")], by.x="sampleID", by.y = "row.names", all.x = T ) 

x = as.data.frame(melted_all %>% group_by(variable) %>%
  summarise(mean = mean(value)))
x = x[order(x$mean, decreasing = T),]
melted_all$variable = factor(melted_all$variable, levels = as.vector(x$variable))

melted_all$soiltype = factor(melted_all$soiltype, levels=c("L", "H"))
p = ggplot(melted_all, aes(x=variable, y=value, color=subspecies)) +
  geom_boxplot(position = "dodge", alpha=1, outlier.size=0.3, size=0.5, width=0.7, fill="transparent") +
  main_theme +
  facet_grid(soiltype ~ .)+
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
p
ggsave(paste("fig1/SF/08.faprotax_IND_boxplot", ".pdf", sep=""), p, width = 8, height = 5)

x = head(x, 9)
melted_all = melted_all[melted_all$variable %in% x$variable,]
melted_all$soiltype = factor(melted_all$soiltype, levels=c("L", "H"))
p = ggplot(melted_all, aes(x=variable, y=value, color=subspecies)) +
  geom_boxplot(position = "dodge", alpha=1, outlier.size=0.3, size=0.5, width=0.7, fill="transparent") +
  main_theme +
  facet_grid(soiltype ~ .)+
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
p
ggsave(paste("fig1/SF/08.faprotax_IND_boxplot_top9", ".pdf", sep=""), p, width = 8, height = 5)


list = read.table("result/compare_far/TEJ.list", header=F,  sep="\t")
sub_alpha=alpha[rownames(sub_design), as.vector(list$V1)]
library(dplyr)
sub_alpha$sampleID=rownames(sub_alpha) 
melted = melt(sub_alpha, id.vars = "sampleID")

melted_all = merge(melted, sub_design[,c("subspecies","soiltype")], by.x="sampleID", by.y = "row.names", all.x = T ) 

x = as.data.frame(melted_all %>% group_by(variable) %>%
                    summarise(mean = mean(value)))
x = x[order(x$mean, decreasing = T),]
melted_all$variable = factor(melted_all$variable, levels = as.vector(x$variable))

melted_all$soiltype = factor(melted_all$soiltype, levels=c("L", "H"))
p = ggplot(melted_all, aes(x=variable, y=value, color=subspecies)) +
  geom_boxplot(position = "dodge", alpha=1, outlier.size=0.3, size=0.5, width=0.7, fill="transparent") +
  main_theme +
  facet_grid(soiltype ~ .)+
  theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
p
ggsave(paste("fig1/SF/08.faprotax_TEJ_boxplot", ".pdf", sep=""), p, width = 8, height = 5)



method = c("nitrate_ammonification","nitrogen_fixation")
for(m in method){
	model = aov(index[[m]] ~ group, data=index)
	Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
	Tukey_HSD_table = as.data.frame(Tukey_HSD$group) 
	write.table(paste(m, "\n\t", sep=""), file=paste("result/faprotax/",m,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
	suppressWarnings(write.table(Tukey_HSD_table, file=paste("result/faprotax/",m,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

	out = LSD.test(model,"group", p.adj="none") # alternative fdr
	stat = out$groups
	index$stat=stat[as.character(index$group),]$groups
	max=max(index[,c(m)])
	min=min(index[,c(m)])
	x = index[,c("group",m)]
	y = x %>% group_by(group) %>% summarise_(Max=paste('max(',m,')',sep=""))
	y=as.data.frame(y)
	rownames(y)=y$group
	index$y=y[as.character(index$group),]$Max + (max-min)*0.05
	
if (FALSE){
	write.table(paste("SampleID\t", sep=""), file=paste("result/faprotax/",m,"_raw.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
	suppressWarnings(write.table(index[,c(m,"group")], file=paste("result/faprotax/",m,"_raw.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
}

	p = ggplot(index, aes(x=group, y=index[[m]], color=group)) +
		geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
		labs(x="Groups", y=paste(m, "index")) + theme_classic() + main_theme +
		geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
		geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
	if (length(unique(sub_design$group))>3){
		p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
	}


	p
	ggsave(paste("result/faprotax/", m, ".pdf", sep=""), p, width = 8, height = 5)
	ggsave(paste("result/faprotax/", m, ".png", sep=""), p, width = 8, height = 5)
	print(paste("Output in result/faprotax/", m, ".txt/pdf finished.", sep = ""))
}

