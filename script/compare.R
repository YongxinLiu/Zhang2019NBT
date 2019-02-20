






rm(list=ls()) 


site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list = c("limma","ggplot2","pheatmap","dplyr","devtools")
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

package_list = c("edgeR")
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


write.table("ID\ttype", file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare//diff.list", sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F)

design = read.table("/mnt/bai/yongxin/rice/xianGeng/doc/design.txt", header=T, row.names=1, sep="\t") # , comment.char=""
design$group = design$groupID

if (TRUE){
	design = subset(design, group %in% c("HTEJ","HIND","LTEJ","LIND","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50LnCp6","A56LnCp6"))
	design$group  = factor(design$group, levels=c("HTEJ","HIND","LTEJ","LIND","V3703HnCp6","ZH11HnCp6","A50LnCp7","A56LnCp7","V3703LnCp6","ZH11LnCp6","A50LnCp6","A56LnCp6"))
}

otutab = read.table(paste("/mnt/bai/yongxin/rice/xianGeng/result/otutab.txt", sep=""), header=T, row.names=1, sep="\t", comment.char="") 

idx = rownames(design) %in% colnames(otutab)
design = design[idx,]
otutab = otutab[,rownames(design)]

if (TRUE){
	norm = t(otutab)/colSums(otutab,na=T)*100
}else{
	norm=t(otutab)/100
}
grp = design[, "groupID2", drop=F]
mat_t2 = merge(grp, norm, by="row.names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=median) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno
filtered = mat_mean_final[apply(mat_mean_final,1,max) > 0.05, ] # select OTU at least one sample > 0.1%
otutab = otutab[rownames(filtered),]

mat_mean_high = mat_mean_final[rownames(filtered),]
write.table(paste("OTUID\t",sep=""), file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", "database.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(round(mat_mean_high*100,3), file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", "database.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T))

print(paste("你正在使用秩和检验！Now, you are using wilcoxon test!", sep=" "))

compare_DA = function(compare){
	group_list = as.vector(as.matrix(compare))
	SampAvsB=paste(group_list[1] ,"-", group_list[2], sep="")
	idx = design$group %in% group_list
	sub_design=design[idx,]
	sub_dat=as.matrix(otutab[,rownames(sub_design)])

	if (TRUE){
		sub_norm = t(t(sub_dat)/colSums(sub_dat,na=T))*100
	}else{
		sub_norm = as.matrix(otutab)/100 # 数据类型一致，计算后矩阵
	}
	idx = sub_design$group %in% group_list[1]
	GroupA = sub_norm[,rownames(sub_design[idx,])]
	idx = sub_design$group %in% group_list[2]
	GroupB = sub_norm[,rownames(sub_design[idx,])]
	
	nrDAO = data.frame(list=rownames(sub_norm), row.names =rownames(sub_norm) )
	for ( i in 1:dim(nrDAO)[1]){
		FC = (mean(GroupA[i,])+0.0001)/(mean(GroupB[i,])+0.0001)
		nrDAO[i,2]=log2(FC)
		nrDAO[i,3]=log2(max(c(GroupA[i,],GroupB[i,]))*10000)
		nrDAO[i,4]= wilcox.test(as.numeric(GroupA[i,]),as.numeric(GroupB[i,]))$p.value
	}	
	nrDAO=nrDAO[,-1]
	colnames(nrDAO)=c("logFC", "logCPM", "PValue")
	nrDAO$FDR = p.adjust(nrDAO$PValue, method="fdr", dim(nrDAO)[1])    #p.adjust就是计算FDR的包，这个可要记得了
	

	nrDAO$logFC=round(nrDAO$logFC,3)
	nrDAO$logCPM=round(nrDAO$logCPM,3)
	nrDAO$level = ifelse(nrDAO$logFC>log2(1.2) & nrDAO$PValue<0.01 & nrDAO$FDR<0.05, "Enriched",ifelse(nrDAO$logFC<log2(1.2)*(-1) & nrDAO$PValue<0.01 & nrDAO$FDR<0.05, "Depleted","NotSig"))
	nrDAO$level=factor(nrDAO$level,levels = c("Enriched","Depleted","NotSig"))

	A_list = subset(sub_design, group %in% group_list[1])
	A_norm = sub_norm[, rownames(A_list)]
	A_mean = as.data.frame(rowMeans(A_norm))
	colnames(A_mean)=c("MeanA")
	B_list = subset(sub_design, group %in% group_list[2])
	B_norm = sub_norm[, rownames(B_list)]
	B_mean = as.data.frame(rowMeans(B_norm))
	colnames(B_mean)=c("MeanB")
	Mean = round(cbind(A_mean, B_mean, A_norm, B_norm),3)
	Mean = Mean[rownames(nrDAO),]   

	if (file.exists("result/taxonomy_8.txt")){
	tax = read.table("result/taxonomy_8.txt", header=T, row.names= 1, sep="\t", comment.char = "") 
	tax = tax[rownames(nrDAO),]
	Mean=cbind(tax, Mean)
	}

	output=cbind(nrDAO,Mean)

	write.table(paste(group_list[1],"_", group_list[2], "\t",sep=""), file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", SampAvsB, "_all.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
	suppressWarnings(write.table(output,file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", SampAvsB, "_all.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T))

	NoE= dim(output[output$level=="Enriched",])[1]
	NoD= dim(output[output$level=="Depleted",])[1]
	NoN= dim(output[output$level=="NotSig",])[1]
	suppressWarnings(write.table(paste( SampAvsB, NoE, NoD, NoN, sep="\t"), file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", "summary.txt",sep=""), append = T, quote = F, sep = '\t', row.names = F, col.names = F))

	output=output[output$level!="NotSig",]
	write.table(paste(group_list[1],"_", group_list[2], "\t",sep=""), file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", SampAvsB, "_sig.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
	suppressWarnings(write.table(output, file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", SampAvsB, "_sig.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T))
	
	if (dim(output)[1]>1){
	write.table(cbind(rownames(output),paste(group_list[1],"_", group_list[2], output$level, sep="")), file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare//diff.list", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
	}
}


write.table("GroupAvsB\tEnriched\tDepleted\tNotSig\n", file=paste("/mnt/bai/yongxin/rice/xianGeng/result/compare/", "summary.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)

if (!file.exists("/mnt/bai/yongxin/rice/xianGeng/doc/compare.txt")) {
	compare_data = as.vector(unique(sub_design$group))
	len_compare_data = length(compare_data)
	for(i in 1:(len_compare_data-1)) {
		for(j in (i+1):len_compare_data) {
			tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
			compare_DA(tmp_compare)
			print(paste(Compared, tmp_compare[,1], "vs",tmp_compare[,2] ,"finished!", sep=" "))
		}
	}
}else {
	compare_data = read.table("/mnt/bai/yongxin/rice/xianGeng/doc/compare.txt", sep="\t", check.names=F, quote='', comment.char="")
	colnames(compare_data) = c("sampA", "sampB")
	for(i in 1:dim(compare_data)[1]){
		compare_DA(compare_data[i,])
		print(paste("Compared", compare_data[i,1], "vs", compare_data[i,2],"finished!", sep=" "))
	}
}

system("sed -i 's/Enriched/_E/;s/Depleted/_D/;' /mnt/bai/yongxin/rice/xianGeng/result/compare//diff.list")

