# The investigation of the immune landscape

# Sankey diagram
options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggalluvial)
library(grid)
library(gtable)
clindf <- read.table("Sankey.txt", head=TRUE, sep="\t", quote="")
dim(clindf)
head(clindf)

clindf <- clindf
for (i in 1:nrow(clindf)) {
    clindf[i, 4:ncol(clindf)] <- clindf[i, 4:ncol(clindf)] /
        sum(clindf[i, 4:ncol(clindf)])
}


newdf1 <- clindf
newlst1 <- list()
for (i in 4:ncol(newdf1)) {
    newdf2 <- data.frame(
        id=rownames(newdf1),
        group=colnames(newdf1)[i],
        dataset=newdf1[, "Type"],
        risk=newdf1[, "Risk"],
        value=newdf1[, i])
    newlst1[[i]] <- newdf2
}
newdf3 <- do.call(rbind, newlst1)

newdf4 <- aggregate(value ~ group, newdf3, mean)
newdf4 <-newdf4[order(newdf4[, "value"]),]
groups <- newdf4[, "group"]

tb <- table(newdf3[, "dataset"])
names <- names(tb[order(-tb)])
datasets <- c()
for (i in 1:length(names)) {
    if (i %% 2 == 1) {
        datasets <- c(names[i], datasets)
    } else {
        datasets <- c(datasets, names[i])
    }
}

risks <- c("high", "low")

newdf3[, "group"] <- factor(newdf3[, "group"], levels=groups)
newdf3[, "dataset"] <- factor(newdf3[, "dataset"], levels=datasets)
newdf3[, "risk"] <- factor(newdf3[, "risk"], levels=risks)
newdf3 <- newdf3[
    order(newdf3[, "risk"], newdf3[, "group"], -newdf3[, "value"]),]
ids <- newdf3[newdf3[, "group"] == groups[length(groups)], "id"]
newdf3[, "id"] <- factor(newdf3[, "id"], levels=ids)


ids2 <- c()
for (i in 1:length(ids)) {
    if (i %% 2 == 1) {
        ids2 <- c(ids[i], ids2)
    } else {
        ids2 <- c(ids2, ids[i])
    }
}
newdf3[, "id"] <- factor(newdf3[, "id"], levels=ids2)


colors1 <- c(
rgb(205,28,0, max=255),       
rgb(202,109,107, max=255),   
rgb(133,109,82, max=255),     
rgb(145,34,30, max=255),     
rgb(209,149,69, max=255),    
rgb(158,109,56, max=255),     
rgb(203,106,57, max=255),     
rgb(58,204,176, max=255),    
rgb(110,110,207, max=255),    
rgb(116,136,211, max=255),   
rgb(77,196,128, max=255),    
rgb(133,180,107, max=255),   
rgb(79,115,47, max=255),      
rgb(130,171,71, max=255),     
rgb(187,204,73, max=255),    
rgb(153,167,29, max=255),     
rgb(188,165,103, max=255),   
rgb(54,71,146, max=255),    
rgb(86,54,143, max=255),    
rgb(0,153,132, max=255),     
rgb(44,101,132, max=255),   
rgb(210,92,161, max=255),   
rgb(212,128,213, max=255),   
rgb(187,87,198, max=255),     
rgb(143,51,144, max=255),    
rgb(202,109,107, max=255),  
rgb(215,97,137, max=255),   
rgb(197,51,115, max=255),   
rgb(227,68,90, max=255),     
rgb(165,0,81, max=255)     
)
colors1


colors2 <- c(
rgb(205,28,0, max=255),     
rgb(202,109,107, max=255),    
rgb(133,109,82, max=255),     
rgb(145,34,30, max=255),     
rgb(209,149,69, max=255),   
rgb(158,109,56, max=255),    
rgb(203,106,57, max=255),     
rgb(58,204,176, max=255),    
rgb(110,110,207, max=255),   
rgb(116,136,211, max=255),   
rgb(77,196,128, max=255),    
rgb(133,180,107, max=255),    
rgb(79,115,47, max=255),      
rgb(130,171,71, max=255),     
rgb(187,204,73, max=255),    
rgb(153,167,29, max=255),    
rgb(188,165,103, max=255),   
rgb(54,71,146, max=255),     
rgb(86,54,143, max=255),     
rgb(0,153,132, max=255),    
rgb(44,101,132, max=255),    
rgb(210,92,161, max=255),    
rgb(212,128,213, max=255),   
rgb(187,87,198, max=255),    
rgb(143,51,144, max=255),   
rgb(202,109,107, max=255),  
rgb(215,97,137, max=255),     
rgb(197,51,115, max=255),    
rgb(227,68,90, max=255),      
rgb(165,0,81, max=255),    
rgb(235,126,107, max=255),     
rgb(84,185,167, max=255),     
rgb(112,132,166, max=255)
)
colors2

colors3 <- c("red", "green") 
colors3 <- rainbow(2, s=0.7, v=0.7)
colors3 <- c(rgb(152,152,220, max=255), rgb(224,167,224, max=255))
colors3

p <- ggplot(newdf3, aes(x=id, y=value, fill=group)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_x_discrete("") +
    scale_y_continuous("", expand=c(0, 0)) +
    scale_fill_manual(values=c(colors1)) +
    facet_grid(~ risk, scales="free", space="free") +
    guides(fill=guide_legend(ncol=1, byrow=TRUE)) +
    theme(
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        legend.title=element_text(),
        legend.position="none",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black"),
        strip.text.x=element_text(color="transparent"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),
        panel.spacing.x=unit(0, "lines"),
        plot.background=element_blank(),
        plot.margin=margin(unit(c(2, 2, 2, 2), "cm")))
p1 <- p

p <- ggplot(newdf3, aes(x=id, y=value, fill=group)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_manual(values=colors1) +
    guides(fill=guide_legend(nrow=3, byrow=TRUE, title.position="left")) +
    theme(
        legend.title=element_text(),
        legend.position="right",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black"),
        plot.margin=margin(unit(c(2, 2, 2, 2), "mm")))
p2 <- p

p <- ggplot(newdf3, aes(x=id, y=value, fill=dataset)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_manual(values=colors2) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE, title.position="left")) +
    theme(
        legend.title=element_text(),
        legend.position="right",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black"),
        plot.margin=margin(unit(c(2, 2, 2, 2), "mm")))
p3 <- p

p <- ggplot(newdf3, aes(x=id, y=value, fill=risk)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_manual(values=colors3) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position="left")) +
    theme(
        legend.title=element_text(),
        legend.position="right",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black"),
        plot.margin=margin(unit(c(2, 2, 2, 2), "mm")))
p4 <- p

clindf2 <- clindf
clindf2[, "Type"] <- factor(clindf2[, "Type"], levels=rev(datasets))
clindf2[, "Risk"] <- factor(clindf2[, "Risk"], levels=rev(risks))
p <- ggplot(
    clindf2, aes(axis1=Risk, axis2=Type)) +
    geom_alluvium(aes(fill=Type), width=1/12, alpha=1) +
    geom_stratum(width=1/12, fill=c(colors3, colors2), color="transparent") +
    scale_fill_manual(values=rev(colors2)) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    theme(
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        legend.position="none",
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank(),
        plot.margin=margin(unit(c(2, 2, -20, 2), "cm"))) + coord_flip()
p5 <- p

g <- ggplotGrob(p1)
strips <- which(grepl('strip-', g$layout$name))
for (i in 1:length(strips)) {
    j <- which(grepl("rect", g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[strips[i]]]$grobs[[1]]$children[[j]]$gp$fill <- colors3[i]
}
g1 <- g

g <- ggplotGrob(p2)
guide <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
g2 <- g$grobs[[guide]]

g2$heights[2] <- unit(0.85, "cm")

g <- ggplotGrob(p3)
guide <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
g3 <- g$grobs[[guide]]

g3$heights[2] <- unit(0.85, "cm")

g <- ggplotGrob(p4)
guide <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
g4 <- g$grobs[[guide]]

g4$heights[2] <- unit(0.85, "cm")

g5 <- ggplotGrob(p5)

t <- gtable(unit(38, "cm"), unit(c(10, 10, 2.9, 2, 1.1), "cm"))

t <- gtable_add_grob(t, g5, 1, 1, 1, 1)
t <- gtable_add_grob(t, g1, 2, 1, 2, 1)
t <- gtable_add_grob(t, g2, 3, 1, 3, 1)
t <- gtable_add_grob(t, g3, 4, 1, 4, 1)
t <- gtable_add_grob(t, g4, 5, 1, 5, 1)
grid.draw(t)

outfile <- "fig_sankey_bar_risk_bottom.pdf"
pdf(outfile, onefile=FALSE, height=28/2.54, width=40/2.54)
grid.draw(t)
dev.off()








# Single-sample gene set enrichment analysis (ssGSEA) 
inputFile="EBPlusPlusAdjustPANCAN_exp_log.txt"                         
gmtFile="immune.gmt"                                         

library(GSVA)
library(limma)
library(GSEABase)
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)







# Comparison of the immune activities estimated by the ssGSEA method between the INB and IB groups
options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)
clindf <- read.table("immune.txt", head=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6, 1:6]

newlst1 <- list()
vars <- colnames(clindf)[4:ncol(clindf)]
for (i in 1:length(vars)) {
    newlst1[[i]] <- data.frame(
        cell=vars[i], group=clindf[, "risk"], value=clindf[, vars[i]]
    )
}
newdf2 <- do.call(rbind, newlst1)
max(newdf2[, "value"])
median(newdf2[, "value"])
mean(newdf2[, "value"])

levels <- levels(factor(newdf2[, "cell"]))
newlst1 <- list()
for (i in 1:length(levels)) {
    newdf3 <- newdf2[newdf2[, "cell"] == levels[i],]
    newdf3[, "group"] <- factor(newdf3[, "group"])
    model <- wilcox.test(value ~ group, newdf3)
    print(model)
    means <- aggregate(value ~ factor(group), newdf3, mean)
    medians <- aggregate(value ~ factor(group), newdf3, median)
    maxs <- aggregate(value ~ factor(group), newdf3, max)
    newlst1[[i]] <- data.frame(
        var=vars[i],
        mean_high=means[1, "value"],
        mean_low=means[2, "value"],
        median_high=medians[1, "value"],
        median_low=medians[2, "value"],
        max_high=maxs[1, "value"],
        max_low=maxs[2, "value"],
        pvalue=model[["p.value"]]
    )
}
newdf4 <- do.call(rbind, newlst1)
write.csv(newdf4, "result_mean_2.csv", row.names=FALSE, na="")

windowsFonts()
library(scales)
colors1 <- rainbow(29, s=1, v=1, alpha = 1)
barplot(1:29,col=rainbow(29, s=1, v=1, alpha = 1))
colors1

colors2 <- rainbow(29, s=1, v=0.6, alpha = 1)
barplot(1:29,col=rainbow(29, s=1, v=0.6, alpha = 1))
colors2

maxvalues <- apply(newdf4[, c("max_high", "max_low")], 1, max)
for (i in 1:length(levels)) {
newdf3 <- newdf2[newdf2[, "cell"] == levels[i],]

colors <- c(colors1[i], colors2[i])
titlex <- levels[i]
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)

labely <- maxvalues[i] *1.15
fig <- ggplot(newdf3, aes(x=group, y=value, group=group)) +

    geom_boxplot(
        aes(colour=group), notch=FALSE,
        outlier.size=0.6, outlier.shape=1, outlier.alpha=0.5,
        size=0.8,

        fatten=0.8
    ) +
    scale_color_manual(values=colors) +

    scale_x_discrete(titlex, labels=c("cold tumor", "hot tumor")) +
    scale_y_continuous("Value") +
    stat_compare_means(
        aes(
            label=ifelse(
                p < 0.001, "P < 0.001", paste0("P = ", ..p.format..)
            )
        ),

        label.x=1.25, label.y=labely, family="sans", size=2.5
    ) +
    stat_compare_means(
        comparisons=comparisons,
        method="wilcox.test",
        symnum.args=symnum.args
    ) +
    theme_bw() +
    theme(
        axis.title.x=element_text(family="sans", size=7),

        axis.title.y=element_blank(),

        axis.text.x=element_text(
            color="black", angle=30, hjust=1, family="sans", size=7),
        axis.text.y=element_text(color="black", family="sans", size=7),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
fig

outfile <- paste0("fig_boxplot_", i, "_", levels[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_boxplot_", i, "_", levels[i], ".pdf")
pdf(outfile, width=2/2.54, height=4/2.54)
print(fig)
dev.off()
}

colors <- c("dodgerblue4", "dodgerblue")
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)
fig <- ggplot(newdf2, aes(x=group, y=value, group=group)) +
    geom_boxplot(
        aes(colour=group), notch=FALSE, outlier.size=1,
        size=1.6, fatten=1
    ) +

    scale_colour_manual(values=colors) +
    scale_x_discrete("") +
    scale_y_continuous("Value", limit=c(0, 20), breaks=seq(0, 18, 3)) +
    stat_compare_means(
        aes(
            label=ifelse(
                p < 0.001, "P < 0.001", paste0("P = ", ..p.format..)
            )
        ), label.x=1.2, label.y=19
    ) +
    stat_compare_means(
        comparisons=comparisons,
        method="wilcox.test",
        symnum.args=symnum.args
    ) +
    facet_wrap(~ cell, ncol=6, strip.position="bottom") +
    theme_bw() +
    theme(
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(color="black"),
        axis.line.y=element_line(),
        legend.title=element_blank(),
        legend.position="top",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA),
        strip.text.x=element_text(color="black", angle=30, hjust=1),
        strip.background=element_rect(colour=NA, fill=NA)
    )
fig

outfile <- "fig_boxplot_all.tiff"
tiff(outfile, width=20, height=20, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- "fig_boxplot_all.pdf"
pdf(outfile, width=20/2.54, height=20/2.54)
print(fig)
dev.off()







# Hierarchical clustering
library(sparcl)
data=read.table("ssgseaOut.txt",sep="\t",header=T,check.names=F,row.names=1)
hc = hclust(dist(t(data)))

y=cutree(hc,2) 

write.table(y,file="cluster.txt",sep="\t",quote=F,col.names=F)

pdf(file="hclust.pdf",width=50,height=20)
ColorDendrogram(hc, y = y, labels = names(y), branchlength = 0.3,xlab=" ",sub=" ",main = " ")
dev.off()





# Heatmap
library(pheatmap)
rt=read.table("ssgseaOut.txt",sep="\t",header=T,row.names=1,check.names=F)    #读取文件
Type=read.table("cluster.txt",sep="\t",check.names=F,header=F)
Type[,2]=paste0("Cluster",Type[,2])
Type=Type[order(Type[,2]),]
rt=rt[,as.vector(Type[,1])]
dim(rt)
rt[1:6,1:6]
max(rt)
min(rt)

cluster=as.data.frame(Type[,2])
row.names(cluster)=Type[,1]
colnames(cluster)="Cluster"

pdf("heatmap.pdf",height=3,width=7)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
bk

pheatmap(rt, 
         annotation=cluster,
         cluster_col=F,
         cluster_rows = TRUE,
         clustering_distance_row="euclidean", 
         clustering_method = "ward.D2",
         fontsize=6,
         fontsize_row=6,
         scale="row",
         show_colnames=F,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         fontsize_col=3)
dev.off()

Cluster1="Immunity_H"
Cluster3="Immunity_L"

a=c()
a[Type[,2]=="Cluster1"]=Cluster1

a[Type[,2]=="Cluster3"]=Cluster3
clusterOut=cbind(Type,a)
write.table(clusterOut,file="cluster.Immunity.txt",sep="\t",quote=F,col.names=F,row.names=F)






# Overlap between immune subtypes
options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)

clusdf <- read.table("clusters.txt", head=TRUE, sep="\t", quote="")
dim(clusdf)
head(clusdf)

CNVCor_METCor_nmf_cluster_table <- table(clusdf[, 2], clusdf[, 3])
CNVCor_METCor_nmf_cluster_table

chisq.test(CNVCor_METCor_nmf_cluster_table)

library(gplots)
outfile <- "fig_ballon.pdf"
pdf(outfile, onefile=FALSE, width=20/2.54, height=20/2.54)
balloonplot(
    CNVCor_METCor_nmf_cluster_table,
    xlab ="", ylab="",label = FALSE, show.margins = FALSE)
chisq.test(CNVCor_METCor_nmf_cluster_table)$p.value->table_pvalue;
if(table_pvalue<1e-5){
    legend("topleft",legend=paste("chisq-p:","<1e-5"))
}else{
    legend("topleft",legend=paste("chisq-p:",round(table_pvalue,5)))
}
dev.off()








# Comparison of the immune activities estimated by the ESTIMATE method between the INB and IB groups
library(limma)
library(estimate)
inputFile="EBPlusPlusAdjustPANCAN_exp_log.txt"                                                 

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
dim(rt)
rt[1:6,1:6]

sum(duplicated(rt[, "gene_id"]))
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dim(data)
data[1:6,1:6]

out=data[rowMeans(data)>0,]
dim(out)
out[1:6,1:6]

out=rbind(ID=colnames(out),out)
dim(out)
out[1:6,1:6]

write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")


scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)







# Bar plot
options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)

clindf <- read.table("estimate.txt", head=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6, 1:6]

newlst1 <- list()
vars <- colnames(clindf)[3:ncol(clindf)]
for (i in 1:length(vars)) {
    newlst1[[i]] <- data.frame(
        type=vars[i], group=clindf[, "risk"], value=clindf[, vars[i]]
    )
}
newdf2 <- do.call(rbind, newlst1)
max(newdf2[, "value"])
median(newdf2[, "value"])
mean(newdf2[, "value"])

levels <- levels(factor(newdf2[, "type"]))
newlst1 <- list()
for (i in 1:length(levels)) {
    newdf3 <- newdf2[newdf2[, "type"] == levels[i],]
    newdf3[, "group"] <- factor(newdf3[, "group"])
    model <- wilcox.test(value ~ group, newdf3)
    print(model)
    means <- aggregate(value ~ factor(group), newdf3, mean)
    medians <- aggregate(value ~ factor(group), newdf3, median)
    maxs <- aggregate(value ~ factor(group), newdf3, max)
    mins <- aggregate(value ~ factor(group), newdf3, min)
    newlst1[[i]] <- data.frame(
        var=vars[i],
        mean_high=means[1, "value"],
        mean_low=means[2, "value"],
        median_high=medians[1, "value"],
        median_low=medians[2, "value"],
        max_high=maxs[1, "value"],
        max_low=maxs[2, "value"],
        min_high=mins[1, "value"],
        min_low=mins[2, "value"],
        pvalue=model[["p.value"]]
    )
}
newdf4 <- do.call(rbind, newlst1)
write.csv(newdf4, "result_mean_2.csv", row.names=FALSE, na="")

windowsFonts()
colors1 <- c(
    rgb(133, 245, 173, max=255),
    rgb(255, 155, 243, max=255),
    rgb(255, 178, 115, max=255),
    rgb(219, 145, 235, max=255),
    "red4", "orange4", "blue4", "green4", "yellow4", "pink4",
    "cyan4", "darkorange4", "seagreen4"
)

colors2 <- c(
    rgb(41, 160, 95, max=255), 
    rgb(196, 28, 134, max=255), 
    rgb(213, 97, 39, max=255), 
    rgb(135, 67, 152, max=255),
    "red1", "orange1", "blue1", "green1", "yellow1", "pink1",
    "cyan1", "darkorange1", "seagreen1"
)

maxvalues <- apply(newdf4[, c("max_high", "max_low")], 1, max)
minvalues <- apply(newdf4[, c("min_high", "min_low")], 1, min)
for (i in 1:length(levels)) {
newdf3 <- newdf2[newdf2[, "type"] == levels[i],]
colors <- c(colors1[i], colors2[i])
titlex <- levels[i]
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)

labely <- maxvalues[i] + (maxvalues[i] - minvalues[i]) * 0.15
fig <- ggplot(newdf3, aes(x=group, y=value, group=group)) +
    geom_boxplot(
        aes(colour=group), notch=FALSE,
        outlier.size=0.6, outlier.shape=1, outlier.alpha=0.5,
        size=1.6,
        fatten=1
    ) +
    scale_color_manual(values=colors) +
    scale_x_discrete(titlex, labels=c("INB", "IB")) +
    scale_y_continuous("Value") +
    stat_compare_means(
        aes(
            label=ifelse(
                p < 0.001, "P < 0.001", paste0("P = ", ..p.format..)
            )
        ),
        label.x=1.25, label.y=labely, family="sans", size=2.8
    ) +
    stat_compare_means(
        comparisons=comparisons,
        method="wilcox.test",
        symnum.args=symnum.args
    ) +
    theme_bw() +
    theme(
        axis.title.x=element_text(family="sans", size=8),
        axis.title.y=element_blank(),
        axis.text.x=element_text(
            color="black", angle=30, hjust=1, family="sans", size=8),
        axis.text.y=element_text(color="black", family="sans", size=8),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
fig

outfile <- paste0("fig_boxplot_", i, "_", levels[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_boxplot_", i, "_", levels[i], ".pdf")
pdf(outfile, width=4/2.54, height=4/2.54)
print(fig)
dev.off()


fig <- ggplot(newdf3, aes(x=group, y=value, group=group)) +
    geom_violin(
        aes(fill=group), color="grey70"
    ) +
    geom_boxplot(
        aes(color=group), notch=FALSE,
        outlier.size=0.3, outlier.shape=1, outlier.alpha=0.5,
        size=0.3, width=0.2,
        fatten=1
    ) +
    geom_boxplot(
        color="grey70", notch=FALSE,
        outlier.colour=NA,
        size=0.3, width=0.2,
        fatten=0.9
    ) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_x_discrete(titlex, labels=c("INB", "IB")) +
    scale_y_continuous("Value") +
    stat_compare_means(
        aes(
            label=ifelse(
                p < 0.001, "P < 0.001", paste0("P = ", ..p.format..)
            )
        ),
        label.x=1.25, label.y=labely, family="sans", size=2.5
    ) +
    stat_compare_means(
        comparisons=comparisons,
        method="wilcox.test",
        symnum.args=symnum.args
    ) +
    theme_bw() +
    theme(
        axis.title.x=element_text(family="sans", size=8),
        axis.title.y=element_blank(),
        axis.text.x=element_text(
            color="black", angle=30, hjust=1, family="sans", size=8),
        axis.text.y=element_text(color="black", family="sans", size=8),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
fig

outfile <- paste0("fig_violin_", i, "_", ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_violin_", i, "_", ".pdf")
pdf(outfile, width=4/2.54, height=4/2.54)
print(fig)
dev.off()
}






# Volcano plot
options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggrepel)
library(ggsci)
clindf <- read.table("immune_activities.txt", head=TRUE, sep="\t", quote="")
dim(clindf)
head(clindf)
colnames(clindf)
sum(clindf[1, 3:ncol(clindf)])
sum(duplicated(substr(clindf[, 1], 1, 12)))


library(limma)
group <- factor(clindf[, "type"])
eset <- t(clindf[, 3:ncol(clindf)])
design <- model.matrix(~ 0 + group)
colnames(design) <- c("high", "low", "normal")
fit <- lmFit(eset, design)

contrast.matrix <- makeContrasts(
    high-normal,
    low-normal,
    levels=design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
vennDiagram(results)


nrow(eset)
x <- topTable(fit2, coef=1, adjust="BH", number=29)
write.csv(x, "result_limma_high.csv")

logFCcut <- 1.5 * 0.01
pvalCut <- 0.05
adjPcut <- 0.05

xmin <- range(x$logFC)[1] * 1.1
xmax <- range(x$logFC)[2] * 1.1
ymin <- 0
ymax <- max(-log10(x$P.Value)) * 1.1

colors <- ifelse(
    x$P.Value < pvalCut & x$logFC > logFCcut,
    "red",
    ifelse(x$P.Value < pvalCut & x$logFC < -logFCcut,"blue", "purple")
)
colors
if (FALSE) {
j <- 1
for (i in 1:length(colors)) {
    if (colors[i] != "purple") {
        colors[i] <- pal_d3("category20b")(sum(colors != "purple"))[j]
        j <- j + 1
    }
}
}
colors

size_values <- c(1, 1.5, 2, 2.5)
size_labels <- c("P>=0.05", "P<0.05", "P<0.01", "P<0.001")
sizes <- ifelse(
   x$P.Value < 0.05,
   ifelse(
       x$P.Value >= 0.01, size_values[2],
       ifelse(
           x$P.Value < 0.01 & x$P.Value >= 0.001, 
           size_values[3], size_values[4]
       )
   ),
   size_values[1]
)

labels <- ifelse(
    x$P.Value < pvalCut & abs(x$logFC) > logFCcut, rownames(x), ""
)
labels

p1 <- ggplot() +
    geom_rect(
        aes(xmin=-Inf, xmax=-logFCcut, ymin=-Inf, ymax=Inf),
        fill="blue", alpha=0.1811
    ) +
    geom_rect(
        aes(xmin=logFCcut, xmax=Inf, ymin=-Inf, ymax=Inf),
        fill="red", alpha=0.1811
    ) +
    geom_rect(
        aes(xmin=-logFCcut, xmax=logFCcut, ymin=-Inf, ymax=Inf),
        fill="purple", alpha=0.1811
    ) +
    geom_point(
        data=x, aes(logFC, -log10(P.Value), size=factor(sizes)),
        color=colors, alpha=0.5
    ) +
    scale_x_continuous(
        bquote(~Log[2]~"(Fold Change)"), limit=c(xmin, -xmin)
    ) +
    scale_y_continuous(
        bquote(~-Log[10]~("P-value")), limit=c(ymin, ymax)
    ) +
    geom_vline(
        xintercept=c(-logFCcut, logFCcut), color=c("blue", "red"), 
        alpha=0.5, linetype="longdash", lwd=1
    ) +
    geom_hline(
        yintercept=-log10(pvalCut), color="purple", alpha=0.5,
        linetype="longdash", lwd=1
    ) +
    theme_bw() +
    theme(
        axis.text=element_text(family="sans", size=12),
        axis.title=element_text(family="sans", size=12),
        legend.position=c(0.8, 0.8),
        legend.background=element_rect(fill=NA),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(family="sans", size=12),
        panel.grid=element_blank()
    )
p2 <- p1 +
    scale_size_manual(values=size_values, labels=size_labels) +
    geom_text_repel(
        data=x, aes(logFC, -log10(P.Value), label=labels),
        family="sans", size=2.5
    ) +
    guides(size=guide_legend(title=NULL))
p2

outfile <- "fig_volcano_high.tiff"
tiff(outfile, width=20, height=20, units="cm", res=350, compression="lzw+p")
p2
dev.off()

outfile <- "fig_volcano_high.pdf"
pdf(outfile, onefile=FALSE, width=15/2.54, height=8/2.54)
print(p2)
dev.off()


nrow(eset)
x <- topTable(fit2, coef=2, adjust="BH", number=29)
write.csv(x, "result_limma_low.csv")

logFCcut <- 1.5 * 0.01
pvalCut <- 0.05
adjPcut <- 0.05

xmin <- range(x$logFC)[1] * 1.1
xmax <- range(x$logFC)[2] * 1.1
ymin <- 0
ymax <- max(-log10(x$P.Value)) * 1.1

colors <- ifelse(
    x$P.Value < pvalCut & x$logFC > logFCcut, "red",
    ifelse(x$P.Value < pvalCut & x$logFC < -logFCcut, "blue", "purple")
)
sum(colors != "purple")

if (FALSE) {
j <- 1
for (i in 1:length(colors)) {
    if (colors[i] != "purple") {
        colors[i] <- pal_d3("category20b")(sum(colors != "purple"))[j]
        j <- j + 1
    }
}
}
colors

size_values <- c(1, 1.5, 2, 2.5)
size_labels <- c("P>=0.05", "P<0.05", "P<0.01", "P<0.001")
sizes <- ifelse(
   x$P.Value < 0.05,
   ifelse(
       x$P.Value >= 0.01, size_values[2],
       ifelse(
           x$P.Value < 0.01 & x$P.Value >= 0.001, 
           size_values[3], size_values[4]
       )
   ),
   size_values[1]
)

labels <- ifelse(
    x$P.Value < pvalCut & abs(x$logFC) > logFCcut, rownames(x), ""
)
labels

p1 <- ggplot() +
    geom_rect(
        aes(xmin=-Inf, xmax=-logFCcut, ymin=-Inf, ymax=Inf),
        fill="blue", alpha=0.1811
    ) +
    geom_rect(
        aes(xmin=logFCcut, xmax=Inf, ymin=-Inf, ymax=Inf),
        fill="red", alpha=0.1811
    ) +
    geom_rect(
        aes(xmin=-logFCcut, xmax=logFCcut, ymin=-Inf, ymax=Inf),
        fill="purple", alpha=0.1811
    ) +
    geom_point(
        data=x, aes(logFC, -log10(P.Value), size=factor(sizes)),
        color=colors, alpha=0.5
    ) +
    scale_x_continuous(
        bquote(~Log[2]~"(Fold Change)"), limit=c(xmin, -xmin)
    ) +
    scale_y_continuous(
        bquote(~-Log[10]~("P-value")), limit=c(ymin, ymax)
    ) +
    geom_vline(
        xintercept=c(-logFCcut, logFCcut), color=c("blue", "red"), 
        alpha=0.5, linetype="longdash", lwd=1
    ) +
    geom_hline(
        yintercept=-log10(pvalCut), color="purple", alpha=0.5,
        linetype="longdash", lwd=1
    ) +
    theme_bw() +
    theme(
        axis.text=element_text(family="sans", size=12),
        axis.title=element_text(family="sans", size=12),
        legend.position=c(0.8, 0.8),
        legend.background=element_rect(fill=NA),
        legend.key=element_rect(fill=NA),
        legend.text=element_text(family="sans", size=12),
        panel.grid=element_blank()
    )
p2 <- p1 +
    scale_size_manual(values=size_values, labels=size_labels) +
    geom_text_repel(
        data=x, aes(logFC, -log10(P.Value), label=labels),
        family="sans", size=2.5
    ) +
    guides(size=guide_legend(title=NULL))
p2

outfile <- "fig_volcano_low.tiff"
tiff(outfile, width=20, height=20, units="cm", res=350, compression="lzw+p")
p2
dev.off()

outfile <- "fig_volcano_low.pdf"
pdf(outfile, onefile=FALSE, width=15/2.54, height=8/2.54)
print(p2)
dev.off()
