# The investigation of extrinsic immune escape mechanisms

# Correlation analysis
options(stringsAsFactors=FALSE)
immudf <- read.table("immune.txt", head=TRUE, sep="\t", quote="")
dim(immudf)
str(immudf)
colnames(immudf)

immudf1 <- immudf[immudf[, "type"] == "high", 3:ncol(immudf)]
cormtx1 <- cor(immudf1)

immudf2 <- immudf[immudf[, "type"] == "low", 3:ncol(immudf)]
cormtx2 <- cor(immudf2)

cor_p <- matrix(0, nrow = ncol(immudf1), ncol = ncol(immudf1))
rownames(cor_p) <- colnames(immudf1)
colnames(cor_p) <- colnames(immudf1)
for (i in 1:ncol(immudf1)){
  for (j in 1:ncol(immudf1)){
      p <- cor.test(immudf1[,i],immudf1[,j])
      cor_p[i,j] <- p$p.value
  }
}
cor_p1 <- cor_p

cor_p <- matrix(0, nrow = ncol(immudf2), ncol = ncol(immudf2))
rownames(cor_p) <- colnames(immudf2)
colnames(cor_p) <- colnames(immudf2)
for (i in 1:ncol(immudf2)){
  for (j in 1:ncol(immudf2)){
      p <- cor.test(immudf2[,i],immudf2[,j])
      cor_p[i,j] <- p$p.value
  }
}
cor_p2 <- cor_p

datR <- cormtx1
for(i in 1:nrow(datR)){
    datR[i,1:i] <- cormtx2[i,1:i]
}
datR[1:3,1:3]

datP <- cor_p1
for (i in 1:nrow(datP)) {
    datP[i,1:i] <- cor_p2[i,1:i]
}
datP[1:3,1:3]

datR[datP > 0.05] <- NA

library(circlize)
colCorRight <-  circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "#ef3b2c"))
colCorLeft <- circlize::colorRamp2(c(-1, 0, 1), c("yellow", "white", "#762a83"))

library(ComplexHeatmap)
p1 <- Heatmap(datR, rect_gp = gpar(type = "none"), 
              show_heatmap_legend = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(col = "grey", fill = NA))
                if(i == j) {
                  grid.circle(x = x, y = y, r = 0.5 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
                  }else if(i > j) {
                    grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                                gp = gpar(fill = colCorLeft(datR[i, j]), col = NA))
                    } else {
                      grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                                  gp = gpar(fill = colCorRight(datR[i, j]), col = NA))
                      }
                },
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = T, show_column_names = T, 
              row_names_side = "right", 
              row_names_rot = 45,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8)
              )

lgdRight <- Legend(col_fun = colCorRight, title = "High", 
                   direction = "horizontal")
lgdLeft <- Legend(col_fun = colCorLeft, title = "Low", 
                  direction = "horizontal")
pd = list(lgdRight, lgdLeft)

pdf("DouleCorPlot.pdf", width = 5, height = 5.5)
draw(p1, annotation_legend_list = pd,
     annotation_legend_side = "top")
dev.off()


newdf1 <- immudf1
vars <- colnames(newdf1)
combdf <- t(combn(vars, 2))
newlst1 <- list()
for (i in 1:nrow(combdf)){
    model <- cor.test(newdf1[, combdf[i, 1]], newdf1[, combdf[i, 2]])
    newlst1[[i]] <- data.frame(
        var1=combdf[i, 1], var2=combdf[i, 2],
        cor=abs(model[["estimate"]]), pvalue=model[["p.value"]]
    )
}
newdf2 <- do.call(rbind, newlst1)
newdf2 <- newdf2[newdf2[, "pvalue"] < 1.05,]
newdf2[, "group"] <- ifelse(newdf2[, "cor"] > 0.4, "r>0.4", "r<=0.4")
tb1 <- table(newdf2[, "group"])
tb1 / nrow(newdf2)
cors1 <- newdf2[, "cor"]
newdf3 <- newdf2
newdf3[, "group"] <- "high"

newdf1 <- immudf2
vars <- colnames(newdf1)
combdf <- t(combn(vars, 2))
newlst1 <- list()
for (i in 1:nrow(combdf)){
    model <- cor.test(newdf1[, combdf[i, 1]], newdf1[, combdf[i, 2]])
    newlst1[[i]] <- data.frame(
        var1=combdf[i, 1], var2=combdf[i, 2],
        cor=abs(model[["estimate"]]), pvalue=model[["p.value"]]
    )
}
newdf2 <- do.call(rbind, newlst1)
newdf2 <- newdf2[newdf2[, "pvalue"] < 1.05,]
newdf2[, "group"] <- ifelse(newdf2[, "cor"] > 0.4, "r>0.4", "r<=0.4")
tb2 <- table(newdf2[, "group"])
tb2 / nrow(newdf2)
cors2 <- newdf2[, "cor"]
newdf4 <- newdf2
newdf4[, "group"] <- "low"


newdf5 <- rbind(newdf3, newdf4)
newdf5 <- newdf5[newdf5[, "cor"] > 0.7,]
write.csv(newdf5, "result_df_cor.csv", row.names=FALSE, na="")

mtx <- matrix(c(tb1, tb2), ncol=2, byrow=TRUE)
mtx
fisher.test(mtx)

wilcox.test(cors1[cors1>0.7], cors2[cors2>0.7])
max(cors1)
min(cors1)
max(cors2)
min(cors2)

library(ggplot2)
library(ggpubr)
df <- data.frame(
    cor=c(cors1[cors1>0.7], cors2[cors2>0.7]),
    group=c(
        rep("high", length(cors1[cors1>0.7])),
        rep("low", length(cors1[cors2>0.7]))
    )
)

windowsFonts()
library(scales)
colors1 <- rainbow(29, s=1, v=1, alpha=1)
barplot(1:29,col=rainbow(29, s=1, v=1, alpha=1))
colors1 <- c(
    rgb(219, 145, 235, max=255),
    rgb(133, 245, 173, max=255),
    rgb(255, 155, 243, max=255),
    rgb(255, 178, 115, max=255),
    "red4", "orange4", "blue4", "green4", "yellow4", "pink4",
    "cyan4", "darkorange4", "seagreen4"
)
colors1

colors2 <- rainbow(29, s=1, v=0.6, alpha=1)
barplot(1:29,col=rainbow(29, s=1, v=0.6, alpha=1))
colors2 <- c(
    rgb(135, 67, 152, max=255),
    rgb(41, 160, 95, max=255),
    rgb(196, 28, 134, max=255),
    rgb(213, 97, 39, max=255),
    "red1", "orange1", "blue1", "green1", "yellow1", "pink1",
    "cyan1", "darkorange1", "seagreen1"
)
colors2

colors <- c(colors1[1], colors2[2])
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)
maxvalues <- max(df[, 1])
minvalues <- min(df[, 1])
labely <- maxvalues[1] *1.05
fig <- ggplot(df, aes(x=group, y=cor, group=group)) +
    geom_boxplot(
        aes(colour=group), notch=FALSE,
        outlier.size=0.6, outlier.shape=1, outlier.alpha=0.5,
        size=0.8,
        fatten=0.9
    ) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_x_discrete("", labels=c("cold tumor", "hot tumor")) +
    scale_y_continuous("Value", limit=c(0.7, 1.05), breaks=seq(0.7, 1, 0.1)) +
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

outfile <- paste0("fig_boxplot", ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_boxplot", ".pdf")
pdf(outfile, width=4/2.54, height=5/2.54)
print(fig)
dev.off()







# Cytolytic activity
options(stringsAsFactors=FALSE)
expdf <- read.table(
    "EBPlusPlusAdjustPANCAN_exp_log.txt", header=TRUE, sep="\t", quote="",
    check.names=FALSE)
dim(expdf)
expdf[1:6, 1:6]

gene_ids <- expdf[, 1]
gene_ids[gene_ids %in% c("GZMA", "PRF1")]

cytdf <- as.data.frame(
    t(expdf[expdf[, 1] %in% c("GZMA", "PRF1"), 2:ncol(expdf)]))
dim(cytdf)
head(cytdf)

colnames(cytdf) <- expdf[expdf[, 1] %in% c("GZMA", "PRF1"), 1]
cytdf[, "CYT"] <- apply(cytdf, 1, function(x) exp(mean(log(x))))
cytdf <- cbind(sample_id=rownames(cytdf), cytdf)
dim(cytdf)
head(cytdf)
write.table(cytdf, "cyt.txt", sep="\t", row.names=FALSE, quote=FALSE)

options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)

clindf <- read.table("group.txt", head=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6, 1:3]
clindf[, "sample_id"] <- clindf[, "Mixture"]

expdf <- read.table("cyt.txt", head=TRUE, sep="\t", quote="")
dim(expdf)
expdf[1:6, 1:3]

clindf <- merge(clindf, expdf, by="sample_id", all=FALSE)
dim(clindf)
clindf[1:6, 1:6]

newlst1 <- list()
vars <- colnames(clindf)[5:ncol(clindf)]
for (i in 1:length(vars)) {
    newlst1[[i]] <- data.frame(
        var=vars[i], group=clindf[, "risk"], value=clindf[, vars[i]]
    )
}
newdf2 <- do.call(rbind, newlst1)
max(newdf2[, "value"])
median(newdf2[, "value"])
mean(newdf2[, "value"])


vars <- colnames(clindf)[5:ncol(clindf)]
outlst <- list()
for (i in 1:length(vars)) {
    newdf3 <- newdf2[newdf2[, "var"] == vars[i],]
    newdf3[, "group"] <- factor(newdf3[, "group"])
    model <- wilcox.test(value ~ group, newdf3)
    print(model)
    means <- aggregate(value ~ factor(group), newdf3, mean)
    medians <- aggregate(value ~ factor(group), newdf3, median)
    maxs <- aggregate(value ~ factor(group), newdf2, max)
    mins <- aggregate(value ~ factor(group), newdf2, min)
    outlst[[i]] <- data.frame(
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
outdf <- do.call(rbind, outlst)
outfile <- "result_mean_2.csv"
write.csv(outdf, outfile, row.names=FALSE, na="")


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

maxvalues <- apply(outdf[, c("max_high", "max_low")], 1, max)
minvalues <- apply(outdf[, c("min_high", "min_low")], 1, min)
for (i in 1:length(vars)) {
newdf3 <- newdf2[newdf2[, "var"] == vars[i],]

windowsFonts()
colors <- c(colors1[i], colors2[i])
titlex <- vars[i]
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)

labely <- maxvalues[i] + (maxvalues[i] - minvalues[i]) * 0.15
fig <- ggplot(newdf3, aes(x=group, y=value, color=group)) +

    geom_boxplot(
        aes(colour=group), notch=FALSE,
        outlier.size=0.6, outlier.shape=1, outlier.alpha=0.5,
        size=1.2,

        fatten=1
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
            color="black", angle=30, hjust=1,
            family="sans", size=8
        ),
        axis.text.y=element_text(color="black", family="sans", size=8),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
fig

outfile <- paste0("fig_boxplot_", i, "_", vars[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_boxplot_", i, "_", vars[i], ".pdf")
pdf(outfile, width=4/2.54, height=5/2.54)
print(fig)
dev.off()
}








# Chemokines
options(stringsAsFactors=FALSE)

expdf <- read.table(
    "EBPlusPlusAdjustPANCAN_exp_log.txt",
    header=TRUE, sep="\t", quote="", check.names=FALSE
)
dim(expdf)
expdf[1:6, 1:6]

dim(expdf)
expdf[1:6, 1:6]

sum(duplicated(expdf[, "gene_id"]))

clindf <- read.table("group.txt", header=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6, 1:3]
colnames(clindf)[1] <- "sample_id"

t_expdf <- t(expdf[, 2:ncol(expdf)])
colnames(t_expdf) <- expdf[, "gene_id"]
t_expdf <- cbind(data.frame(sample_id=rownames(t_expdf)), t_expdf)

clindf2 <- merge(clindf, t_expdf, by="sample_id", all=FALSE)
dim(clindf2)
clindf2[1:6, 1:10]

genedf <- read.table("gene.txt", header=TRUE, sep="\t", quote="")
dim(genedf)
genedf[1:6, 1:2]
genes1 <- genedf[, "gene"]

sum(duplicated(genes1))
as.data.frame(table(genes1))

genes2 <- colnames(clindf2)[3:ncol(clindf2)]
genes1[!genes1 %in% genes2]

genes <- genes1[genes1 %in% genes2]
genes <- genes[!duplicated(genes)]

newlst1 <- list()
for (i in 1:length(genes)) {
    clindf3 <- clindf2[, c("sample_id", "risk", genes[i])]
    colnames(clindf3)[3] <- "value"
    clindf3[, "gene"] <- genes[i]
    newlst1[[i]] <- clindf3
}
clindf4 <- do.call(rbind, newlst1)
clindf4[, "gene"] <- factor(clindf4[, "gene"], levels=genes)
write.table(clindf4, "clin4.txt", sep="\t", row.names=FALSE, na="", quote=FALSE)

newdf1 <- clindf4
newlst1 <- list()
for (i in 1:length(genes)) {
    print(paste0("------- ", genes[i], " --------"))
    newdf2 <- clindf4[clindf4[, "gene"] == genes[i],]
    model <- wilcox.test(value ~ factor(risk), newdf2)
    print(model)
    means <- aggregate(value ~ factor(risk), newdf2, mean)
    medians <- aggregate(value ~ factor(risk), newdf2, median)
    newlst1[[i]] <- data.frame(
        gene=genes[i],
        mean_high=means[1, "value"],
        mean_low=means[2, "value"],
        median_high=medians[1, "value"],
        median_low=medians[2, "value"],
        pvalue=model[["p.value"]])
}
outdf <- do.call(rbind, newlst1)
outfile <- "result_mean_2.csv"
write.csv(outdf, outfile, row.names=FALSE, na="", quote=FALSE)

pvalues <- outdf[, "pvalue"]
pvalues <- ifelse(pvalues < 0.001, "***",
    ifelse(pvalues < 0.01, "**",
    ifelse(pvalues < 0.05, "*", "ns")))

pvaldf <- data.frame(gene=genes, pvalue=pvalues)
write.table(pvaldf, "pval.txt", sep="\t", row.names=FALSE, na="", quote=FALSE)

clindf4 <- read.csv("result_mean_2.csv", header=TRUE, quote="")
dim(clindf4)
clindf4[1:4,1:4]
max(clindf4[,3])
min(clindf4[,3])
median(clindf4[,3])
mean(clindf4[,3])
clindf4 <- merge(genedf, clindf4, by="gene", all=FALSE)
pvaldf <- read.table("pval.txt", header=TRUE, sep="\t", quote="")
clindf4 <- merge(clindf4, pvaldf, by="gene", all=FALSE)
colnames(clindf4)
genes <- clindf4[, "gene"]
types <- clindf4[, "type"]
types <- ifelse(
    types == "Interferons and receptors",
    "Interferons and\n receptors",
    types
)
pvalues <- clindf4[, "pvalue.y"]

library(ComplexHeatmap)
library(circlize)
mtx <- as.matrix(clindf4[, 3:4])
rownames(mtx) <- 1:nrow(clindf4)
colnames(mtx) <- c("INB", "IB")
for (i in 1:nrow(mtx)) {
    mtx[i, ] <- scale(mtx[i,], center=T, scale=F)
}
max(mtx)
min(mtx)
median(mtx)
mean(mtx)


clindf5 <- as.data.frame(t(clindf2[, genes]))
colnames(clindf5) <- clindf2[, 1]
clindf5[, "gene"] <- rownames(clindf5)
clindf5 <- merge(genedf, clindf5, by="gene", all=FALSE)
mtx2 <- as.matrix(clindf5[,3:ncol(clindf5)])
for (i in 1:nrow(mtx2)) {
    mtx2[i, ] <- scale(mtx2[i,], center=T, scale=T)
}
max(mtx2)
min(mtx2)
median(mtx2)
mean(mtx2)

col_fun <- colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3"))
legend_title <- paste0(
    "*   p < 0.05", "\n\n",
    "**  p < 0.01", "\n\n",
    "*** p < 0.001", "\n\n",
    "Expression")

left_annotation <- rowAnnotation(
    gene=anno_text(genes, location=1, just="right",
        width = max_text_width(genes)*1.2)
)

right_annotation <- rowAnnotation(
    pvalue=anno_text(pvalues),
    density=anno_density(
        mtx2, border=FALSE, joyplot_scale=6,
        gp=gpar(fill=rgb(244, 121, 255, max=255)),
        width=unit(2, "cm"), axis=FALSE
    )
)

ht <- Heatmap(
    mtx,
    col=col_fun,
    cluster_columns=FALSE,
    cluster_rows=FALSE,
    show_row_names=FALSE,
    column_names_side="top",
    column_names_rot=0,
    column_names_centered=TRUE,
    left_annotation=left_annotation,
    right_annotation=right_annotation,
    row_split=types,
    row_gap=unit(0.5, "cm"),
    row_title=levels(factor(types)),
    heatmap_legend_param=list(
        title=legend_title,
        title_position="topleft"
    )
)

ht@right_annotation@anno_list$density@name <- ""
draw(ht)

outfile <- "fig_heatmap.tiff"
tiff(outfile, width=25, height=35, unit="cm", res=350, compression="lzw+p")
draw(ht)

decorate_annotation("gene", {
    grid.lines(unit(c(0, 0), "npc"), unit(c(0, 1), "npc"), gp=gpar(lwd=2))
})
decorate_annotation("gene", {
    grid.lines(
        unit(c(0, 0), "npc"), unit(c(-0.19, -0.03), "npc"), gp=gpar(lwd=2)
    )
})
decorate_annotation("gene", {
    grid.lines(
        unit(c(0, 0), "npc"), unit(c(-0.77, -0.22), "npc"), gp=gpar(lwd=2)
    )
})
dev.off()


dim(mtx)
nrow(mtx)
outfile <- "fig_heatmap.pdf"
pdf(outfile, width=12/2.54, height=45/2.54)
draw(ht)
dev.off()









# Fibroblasts
options(stringsAsFactors=FALSE)
expdf <- read.table(
    "EBPlusPlusAdjustPANCAN_exp_log.txt", header=TRUE, sep="\t", quote="",
    check.names=FALSE)
dim(expdf)
expdf[1:6, 1:6]

eset <- expdf[, 2:ncol(expdf)]
rownames(eset) <- expdf[, 1]
eset[1:10, 1:10]

library(MCPcounter)
library(curl)
probesets=read.table("probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character")
genes=read.table("genes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
mcpestimate <- MCPcounter.estimate(
    eset, featuresType="HUGO_symbols", probesets, genes)
heatmap(as.matrix(mcpestimate),col=colorRampPalette(c("blue","white","red"))(100)) 

outdf <- cbind(immu_type=rownames(mcpestimate), mcpestimate)
dim(outdf)
outdf[1:6,1:6]
outdf <- t(outdf)
dim(outdf)
outdf[1:6,1:6]
write.table(outdf, "mcpestimate.txt", sep="\t", row.names=T, col.names=F, quote=FALSE)








# Gene set enrichment analysis
options(stringsAsFactors=FALSE)

expdf <- read.table(
    "EBPlusPlusAdjustPANCAN_exp_log.txt", header=TRUE, sep="\t", quote="", check.names=FALSE)
dim(expdf)
expdf[1:6, 1:6]


clindf <- read.table(
    "group.txt", header=TRUE, sep="\t", quote="", check.names=FALSE)
dim(clindf)
clindf[1:6, 1:18]
rownames(clindf) <- clindf[, 1]

sample_ids1 <- colnames(expdf)[2:ncol(expdf)]
sample_ids2 <- rownames(clindf)
sample_ids <- intersect(sample_ids1, sample_ids2)
sample_ids
groups <- clindf[sample_ids, "risk"]
groups
levels <- groups[!duplicated(groups)]
levels
expmtx <- expdf[, sample_ids]
gene_ids <- expdf[, "gene_id"]
rownames(expmtx) <- gene_ids


outfile <- "risk.cls"
write("", outfile)
conn <- file(outfile, 'r+')
lines <- readLines(conn)
a1 <- paste0(c(length(sample_ids), length(levels), 1), collapse="\t")
a2 <- paste0(c("#", levels), collapse="\t")
a3 <- paste0(groups, collapse="\t")
writeLines(c(a1, a2, a3), con=conn)
close(conn)


outdf <- cbind(NAME=rownames(expmtx), DESCRIPTION="na", expmtx)
outfile <- "risk.gct"
write.table(outdf, outfile, sep="\t", row.names=FALSE, na="", quote=FALSE)
conn <- file(outfile, 'r+')
lines <- readLines(conn)
a1 <- "#1.2"
a2 <- paste0(dim(expmtx), collapse="\t")
writeLines(c(a1, a2, lines), con=conn)
close(conn)



fnames<-Sys.glob("*.xls")
fdataset<-lapply(fnames,read.delim)
names(fdataset) <- fnames
library(plyr)
result <- ldply(fdataset, data.frame)
result$pathway<-unlist(strsplit(result$.id,split = ".xls"))
head(result)



library(ggplot2)
p1<-ggplot(result,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
  geom_point(shape=21) +
  labs(x = "", y = "Enrichment Score", title = "") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),
                     limits =c(min(result$RUNNING.ES-0.02), max(result$RUNNING.ES+0.02))) + 
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) +
  geom_hline(yintercept = 0) +

  guides(fill=guide_legend(title = NULL)) +
  theme(legend.background = element_blank()) +
  theme(legend.key = element_blank())
p1
ggsave(file="point.pdf")


library(RColorBrewer)
p2<-lapply(unique(result$pathway), function(ii) {
    dd <- result[result$pathway == ii,]
    ggplot(dd,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES))+
      geom_bar(stat="identity",colour = "black")+
      labs(x = result$pathway, y = "", title = "") +
      theme_bw() +
      theme(panel.grid =element_blank()) +
      theme(panel.border = element_blank()) + 
      xlim(0,max(result$RANK.IN.GENE.LIST)) +
      theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank())
})

t<-data.frame(a=as.numeric(1:1000),b=as.numeric(1:1000))

p3<-ggplot(t,aes(a,1,fill=b)) +
  geom_tile()+
  theme_bw() +
  labs(x = "high<--------------low", y = "", title = "") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 500)+
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  guides(fill=FALSE)
require(cowplot)
p2[[4]] <- p3
plot_grid(plotlist=p2, ncol=1,rel_heights = c(2,2,2,1))
ggsave(file="bar.pdf")


library(ggplot2)
p4<-ggplot(result,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
  geom_line() + 
  labs(x = "", y = "Enrichment Score", title = "") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),
                     limits =c(min(result$RUNNING.ES-0.02), max(result$RUNNING.ES+0.02))) + 
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) +
  geom_hline(yintercept = 0) + 
  guides(colour=guide_legend(title = NULL)) +
  theme(legend.background = element_blank()) +
  theme(legend.key = element_blank())
p4
ggsave(file="line.pdf")



p5<-ggplot(result,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+
  geom_tile()+
  theme_bw() +
  labs(x = "high<--------------low", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  guides(color=FALSE)
p5


library(gridExtra)
gA=ggplot_gtable(ggplot_build(p4))
gB=ggplot_gtable(ggplot_build(p5))
maxWidth = grid::unit.pmax(gA$widths, gB$widths)
gA$widths <- as.list(maxWidth)
gB$widths <- as.list(maxWidth)

library(grid)
grid.newpage()
grid.arrange(arrangeGrob(gA,gB,nrow=2,heights=c(.8,.6)))

pdf('prettyGSEA.pdf',width=8,height=4)
grid.arrange(arrangeGrob(gA,gB,nrow=2,heights=c(.8,.6)))
dev.off()








# Immunogenomic indicators
options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)

clindf <- read.table("group.txt", head=TRUE, sep="\t", quote="")
clindf[, "group"] <- factor(clindf[, "risk"])
clindf[, "sample_id"] <- substr(clindf[, "Mixture"], 1, 12)
clindf[, "sample_id2"] <- clindf[, "Mixture"]
dim(clindf)
head(clindf)
colnames(clindf)

expdf <- read.csv(
    "The Immune Landscape of Cancer-mmc2.csv", head=TRUE, quote=""
)
dim(expdf)
expdf[1:6, 1:6]

expdf[, "sample_id"] <- expdf[, "TCGA.Participant.Barcode"]
expdf[1:6, 1:6]

sum(duplicated(clindf[, "sample_id"]))
sum(duplicated(expdf[, "sample_id"]))

sample_ids <- intersect(clindf[, "sample_id"], expdf[, "sample_id"])
length(sample_ids)
sample_ids[1:10]
clindf2 <- merge(clindf, expdf, by="sample_id", all=FALSE)
dim(clindf2)
colnames(clindf2)

newdf1 <- clindf2
dim(newdf1)
ncol(newdf1)
for (i in 8:ncol(newdf1)) {
    newdf1[, i] <- log2(newdf1[, i] + 1)
}
vars <- colnames(newdf1)[8:ncol(newdf1)]
outlst <- list()
for (i in 1:length(vars)) {
    print(paste0("------- ", vars[i], " --------"))
    newdf2 <- data.frame(
        value=newdf1[, vars[i]],
        group=newdf1[, "group"])
    newdf2 <- newdf2[complete.cases(newdf2),]
    model <- wilcox.test(value ~ factor(group), newdf2)
    print(model)
    means <- aggregate(value ~ factor(group), newdf2, mean)
    medians <- aggregate(value ~ factor(group), newdf2, median)
    maxs <- aggregate(value ~ factor(group), newdf2, max)
    mins <- aggregate(value ~ factor(group), newdf2, min)
    outlst[[i]] <- data.frame(
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

outdf <- do.call(rbind, outlst)
outfile <- "result_mean_2.csv"
write.csv(outdf, outfile, row.names=FALSE, na="")

windowsFonts()
library(scales)
colors1 <- rainbow(12, s=1, v=1, alpha = 1)
barplot(1:12,col=rainbow(12, s=1, v=1, alpha = 1))
colors1

colors2 <- rainbow(12, s=1, v=0.6, alpha = 1)
barplot(1:12,col=rainbow(12, s=1, v=0.6, alpha = 1))

colors2


maxvalues <- apply(outdf[, c("max_high", "max_low")], 1, max)
minvalues <- apply(outdf[, c("min_high", "min_low")], 1, min)
for (i in 1:length(vars)) {

newdf2 <- data.frame(value=newdf1[, vars[i]], group=newdf1[, "group"])

windowsFonts()
colors <- c(colors1[i], colors2[i])
titlex <- vars[i]
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)

labely <- maxvalues[i] + (maxvalues[i] - minvalues[i]) * 0.15
fig <- ggplot(newdf2, aes(x=group, y=value, color=group)) +

    geom_boxplot(
        aes(colour=group), notch=FALSE,
        outlier.size=0.6, outlier.shape=1, outlier.alpha=0.5,
        size=1.2,

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
            color="black", angle=30, hjust=1,
            family="sans", size=8
        ),
        axis.text.y=element_text(color="black", family="sans", size=8),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
fig

outfile <- paste0("fig_boxplot_", i, "_", vars[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_boxplot_", i, "_", vars[i], ".pdf")
pdf(outfile, width=2/2.54, height=5/2.54)
print(fig)
dev.off()
}
