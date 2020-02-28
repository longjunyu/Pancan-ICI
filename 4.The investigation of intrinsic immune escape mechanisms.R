# The investigation of intrinsic immune escape mechanisms

# Mutation signatures
options(stringsAsFactors=FALSE)
mafdf <- read.table("mc3.v0.2.8.PUBLIC.maf", header=TRUE, sep="\t", quote="")
dim(mafdf)
head(mafdf)
mafdf[1:6,1:6]
head(mafdf[, "Tumor_Sample_Barcode"])
mafdf[, "sample_id"] <- substr(mafdf[, "Tumor_Sample_Barcode"], 1, 16)
head(mafdf[, "sample_id"])
sum(duplicated(mafdf[, "sample_id"]))

silent <- c("Silent", "Intron", "RNA", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR")
keeps <- c(
"Missense_Mutation",
"Frame_Shift_Del",
"Frame_Shift_Ins",
"Nonsense_Mutation",
"Nonstop_Mutation",
"Splice_Site",
"Translation_Start_Site"
)
mafdf <- mafdf[mafdf[, "Variant_Classification"] %in% keeps,]
table(mafdf[, "Variant_Classification"])

sample_ids <- levels(factor(mafdf[, "sample_id"]))
length(sample_ids)
sum(duplicated(sample_ids))

symbols <- levels(factor(mafdf[, "Hugo_Symbol"]))
length(symbols)

mafdf_2 <- mafdf
mafdf_2[, "Hugo_Symbol"] <- factor(mafdf_2[, "Hugo_Symbol"], levels=symbols)
mafdf_2[, "sample_id"] <- factor(mafdf_2[, "sample_id"], levels=sample_ids)
mutdf <- as.data.frame(table(mafdf_2[, "Hugo_Symbol"], mafdf_2[, "sample_id"]))

mutdf[, "mut"] <- ifelse(mutdf[, "Freq"] > 0, 1, 0)
head(mutdf)
genedf <- as.data.frame(
	matrix(
	    mutdf[, "mut"], ncol=length(sample_ids), dimnames=list(symbols, sample_ids)
	)
)
class(genedf)
dim(genedf)
genedf[1:6,1:6]

sample_ids2 <- substr(sample_ids, 1, 15)
dup_ids2 <- sample_ids2[duplicated(sample_ids2)]
sample_ids[sample_ids2 %in% dup_ids2]

keep_ids <- sample_ids[!duplicated(sample_ids2)]
genedf <- genedf[, keep_ids]
dim(genedf)
genedf[1:6,1:6]

clindf <- read.table("group.txt", header=TRUE, sep="\t", quote="")
dim(clindf)
head(clindf)
clindf[, "sample_id"] <- clindf[, "Mixture"]
dim(clindf)
clindf[1:6,1:3]

mafdf2 <- mafdf
mafdf2[, "sample_id"] <- substr(mafdf2[, "Tumor_Sample_Barcode"], 1, 16)
mafdf2 <- merge(mafdf2, clindf, by="sample_id", all=FALSE)
mafdf2[, "mutation_type6"] <- paste0(
	mafdf2[, "Reference_Allele"], ">", mafdf2[, "Tumor_Seq_Allele2"]
)
mafdf2[1:6, "mutation_type6"]
mafdf2[, "mutation_type6"] <- ifelse(
	nchar(mafdf2[, "mutation_type6"]) != 3, "Others", mafdf2[, "mutation_type6"]
)
table(mafdf2[, "mutation_type6"])
mafdf2[, "mutation_type6"] <- ifelse(
	grepl("-", mafdf2[, "mutation_type6"]), "Others", mafdf2[, "mutation_type6"]
)
table(mafdf2[, "mutation_type6"])
mafdf2[, "mutation_type6"] <- ifelse(
	mafdf2[, "mutation_type6"] %in% c("C>A", "G>T"), "C>A",
    ifelse(
    	mafdf2[, "mutation_type6"] %in% c("C>G", "G>C"), "C>G",
		ifelse(
			mafdf2[, "mutation_type6"] %in% c("C>T", "G>A"), "C>T",
		   	ifelse(
		   		mafdf2[, "mutation_type6"] %in% c("T>A", "A>T"), "T>A",
		   	   	ifelse(
		   	   		mafdf2[, "mutation_type6"] %in% c("T>C", "A>G"), "T>C",
		   	   	   	ifelse(mafdf2[, "mutation_type6"] %in% c("T>G", "A>C"), "T>G", NA)
		   	   	)
		   	)
		)
	)
)
table(mafdf2[, "mutation_type6"])
save.image("project_1.RData")




load("maf.RData")
library(ggplot2)
library(grid)
library(gtable)
library(ggsci)


newdf1 <- mafdf2[!is.na(mafdf2[, "mutation_type6"]),]
head(newdf1)
newdf1[, "count"] <- 1
newdf2 <- aggregate(count ~ sample_id + mutation_type6 + risk, newdf1, sum)
newdf3 <- aggregate(count ~ sample_id + risk, newdf1, sum)
colnames(newdf3)[3] <- "total"
newdf4 <- merge(newdf2, newdf3, by=c("sample_id", "risk"), all.x=TRUE)
newdf4[, "percent"] <- newdf4[, "count"] / newdf4[, "total"] * 100
figdf <- newdf4

figdf[, "id"] <- figdf[, "sample_id"]
figdf[, "value"] <- figdf[, "percent"]
figdf[, "group"] <- figdf[, "mutation_type6"]
figdf[, "subgroup"] <- figdf[, "risk"]

sort_group <- c("left", "right", "center")[3]
sort_subgroup <- c("left", "right", "center")[3]

colors_group <- c(
rgb(187,87,198, max=255),
rgb(110,110,207, max=255),
rgb(187,204,73, max=255),
rgb(227,68,90, max=255),
rgb(215,97,137, max=255),
rgb(197,51,115, max=255)
)
colors_group
barplot(1:6,col=colors_group)

colors_subgroup <- pal_jco("default")(2)
colors_subgroup <- c(
rgb(86, 86, 220, max=255),
rgb(224, 77, 224, max=255)
)
barplot(1:2,col=colors_subgroup)

r <- 0.5

title_group <- "Mutation"
title_subgroup <- "Risk"

spacing_subgroup <- unit(0.1 * r, "cm")

width_barplot <- 17 * 2.54 * r
width_group <- 1.5 * 2.54 * r
width_subgroup <- 1.5 * 2.54 * r
width <- width_barplot + width_group + width_subgroup
height <- 10 * 2.54 * r

top_group <- 1 * r
top_subgroup <- 1 * r


if (sort_subgroup == "left") {
    tb <- table(figdf[, "subgroup"])
    subgroups <- names(tb[order(-tb)])
    figdf[, "subgroup"] <- factor(figdf[, "subgroup"], levels=subgroups)
} else if (sort_subgroup == "right") {
    tb <- table(figdf[, "subgroup"])
    subgroups <- names(tb[order(tb)])
    figdf[, "subgroup"] <- factor(figdf[, "subgroup"], levels=subgroups)
} else if (sort_subgroup == "center") {
    tb <- table(figdf[, "subgroup"])
    subgroups <- names(tb[order(-tb)])
    subgroups2 <- c()
    for (i in 1:length(subgroups)) {
        if (i %% 2 == 1) {
            subgroups2 <- c(subgroups[i], subgroups2)
        } else {
            subgroups2 <- c(subgroups2, subgroups[i])
        }
    }
    figdf[, "subgroup"] <- factor(figdf[, "subgroup"], levels=subgroups2)
}

newdf1 <- aggregate(value ~ group, figdf, mean)
groups <- newdf1[order(newdf1[, "value"]), "group"]
figdf[, "group"] <- factor(figdf[, "group"], levels=groups)
if (sort_group == "left") {
    figdf <- figdf[
        order(figdf[, "subgroup"], figdf[, "group"], -figdf[, "value"]),]
    ids <- figdf[figdf[, "group"] == groups[length(groups)], "id"]
    figdf[, "id"] <- factor(figdf[, "id"], levels=ids)
} else if (sort_group == "right") {
    figdf <- figdf[
        order(figdf[, "subgroup"], figdf[, "group"], figdf[, "value"]),]
    ids <- figdf[figdf[, "group"] == groups[length(groups)], "id"]
    figdf[, "id"] <- factor(figdf[, "id"], levels=ids)
} else if (sort_group == "center") {
    figdf <- figdf[
        order(figdf[, "subgroup"], figdf[, "group"], -figdf[, "value"]),]
    ids <- figdf[figdf[, "group"] == groups[length(groups)], "id"]
    ids2 <- c()
    for (i in 1:length(ids)) {
        if (i %% 2 == 1) {
            ids2 <- c(ids[i], ids2)
        } else {
            ids2 <- c(ids2, ids[i])
        }
    }
    figdf[, "id"] <- factor(figdf[, "id"], levels=ids2)
}
head(figdf)

p <- ggplot(figdf, aes(x=id, y=value, fill=group)) +
    geom_bar(position="stack", stat="identity", width=1.1) +
    scale_x_discrete("") +
    scale_y_continuous("", expand=c(0, 0), limit=c(0, 100)) +
    scale_fill_manual(values=colors_group) +
    facet_grid(~ subgroup, scales="free", space="free") +
    guides(fill=guide_legend(ncol=1, byrow=TRUE)) +
    theme(
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_line(),
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
        panel.spacing.x=unit(spacing_subgroup, "cm"),
        plot.background=element_blank()
    )
p1 <- p

p <- ggplot(figdf, aes(x=id, y=value, fill=group)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_manual(values=colors_group) +
    guides(fill=guide_legend(title=title_group, ncol=1, byrow=TRUE)) +
    theme(
        legend.title=element_text(),
        legend.position="right",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black")
    )
p2 <- p

p <- ggplot(figdf, aes(x=id, y=value, fill=subgroup)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_manual(values=colors_subgroup) +
    guides(fill=guide_legend(title=title_subgroup, ncol=1, byrow=TRUE)) +
    theme(
        legend.title=element_text(),
        legend.position="right",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black")
    )
p3 <- p

g <- ggplotGrob(p1)
strips <- which(grepl('strip-', g[["layout"]][["name"]]))
for (i in 1:length(strips)) {
    j <- which(grepl("rect",
        g[["grobs"]][[strips[i]]][["grobs"]][[1]][["childrenOrder"]])
    )
    g[["grobs"]][[strips[i]]][["grobs"]][[1]][[
        "children"]][[j]][["gp"]][["fill"]] <- colors_subgroup[i]
}
g1 <- g

g <- ggplotGrob(p2)
guide <- which(sapply(g[["grobs"]], function(x) x$name) == "guide-box")
g2 <- g[["grobs"]][[guide]]
g2[["heights"]][2] <- unit(top_group, "cm")

g <- ggplotGrob(p3)
guide <- which(sapply(g[["grobs"]], function(x) x$name) == "guide-box")
g3 <- g[["grobs"]][[guide]]
g3[["heights"]][2] <- unit(top_subgroup, "cm")

gt <- gtable(
    unit(c(width_barplot, width_group, width_subgroup), c("cm")),
    unit(height, "cm"))
gt <- gtable_add_grob(gt, g1, 1, 1)
gt <- gtable_add_grob(gt, g2, 1, 2, 1, 2)
gt <- gtable_add_grob(gt, g3, 1, 3, 1, 3)

outfile <- "fig_barplot.pdf"
pdf(outfile, onefile=TRUE, width=width/2.54, height=height/2.54)
grid.draw(gt)
dev.off()




library(ggplot2)
library(ggpubr)

figdf[, "id"] <- figdf[, "sample_id"]
figdf[, "value"] <- figdf[, "percent"]

figdf[, "group"] <- figdf[, "risk"]
figdf[, "var"] <- figdf[, "mutation_type6"]

newdf2 <- figdf
nrow(newdf2)

vars <- levels(factor(newdf2[, "var"]))
newlst1 <- list()
for (i in 1:length(vars)) {
    newdf3 <- newdf2[newdf2[, "var"] == vars[i],]
    newdf3[, "group"] <- factor(newdf3[, "group"])
    model <- wilcox.test(value ~ group, newdf3)
    print(model)
    means <- aggregate(value ~ factor(group), newdf3, mean)
    medians <- aggregate(value ~ factor(group), newdf3, median)
    maxs <- aggregate(value ~ factor(group), newdf3, max)
    mins <- aggregate(value ~ factor(group), newdf3, min)
    newlst1[[i]] <- data.frame(
        var=vars[i],
        mean_1=means[means[, 1] == "high", "value"],
        mean_2=means[means[, 1] == "low", "value"],
        median_1=medians[means[, 1] == "high", "value"],
        median_2=medians[means[, 1] == "low", "value"],
        max_high=maxs[1, "value"],
        max_low=maxs[2, "value"],
        min_high=mins[1, "value"],
        min_low=mins[2, "value"],
        pvalue=model[["p.value"]])
}
newdf4 <- do.call(rbind, newlst1)
write.csv(newdf4, "result_mean_2.csv", row.names=FALSE, na="")


windowsFonts()
library(scales)
colors1 <- rainbow(6, s=1, v=1, alpha=1)
barplot(1:6,col=rainbow(6, s=1, v=1, alpha=1))
colors1


colors2 <- rainbow(6, s=1, v=0.6, alpha=1)
barplot(1:6,col=rainbow(6, s=1, v=0.6, alpha=1))
colors2

maxvalues <- apply(newdf4[, c("max_high", "max_low")], 1, max)
minvalues <- apply(newdf4[, c("min_high", "min_low")], 1, min)

for (i in 1:length(vars)) {
newdf3 <- newdf2[newdf2[, "var"] == vars[i],]

colors <- c(colors1[i], colors2[i])
titlex <- vars[i]
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
    scale_y_continuous("Value", breaks=seq(0, 100, 20)) +
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

outfile <- paste0("fig_boxplot_", i, "_", ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_boxplot_", i, "_", ".pdf")

pdf(outfile, width=5/2.54, height=6/2.54)
print(fig)
dev.off()


fig <- ggplot(newdf3, aes(x=group, y=value, group=group)) +
    geom_violin(
        aes(fill=group)
    ) +

    geom_boxplot(
        colour="black", notch=FALSE,
        outlier.size=0.3, outlier.shape=1, outlier.alpha=0.5,
        size=0.3, width=0.2,

        fatten=0.8
    ) +
    scale_fill_manual(values=colors) +

    scale_x_discrete(titlex, labels=c("cold tumor", "hot tumor")) +
    scale_y_continuous("Value", breaks=seq(0, 100, 20)) +
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

outfile <- paste0("fig_violin_", i, "_", ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_violin_", i, "_", ".pdf")

pdf(outfile, width=5/2.54, height=6/2.54)
print(fig)
dev.off()
}





library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg19)
mut_ref <- mafdf2
mut_ref[1:6, 1:12]
mut_ref[, "Chromosome"] <- paste0("chr", mut_ref[, "Chromosome"])
mut_ref[1:6, 1:12]

sample.mut.ref <- mut_ref
snp_count <- mut.to.sigs.input(
    mut.ref=sample.mut.ref,
    sample.id="sample_id",
    chr="Chromosome",
    pos="Start_Position",
    ref="Reference_Allele",
    alt="Tumor_Seq_Allele2",
    bsg=BSgenome.Hsapiens.UCSC.hg19
)
snp_count[1:6, 1:6]

mut_mtx <- t(snp_count)
dim(mut_mtx)
mut_mtx[1:6, 1:6]

plot_96_profile(mut_mtx[, c(1, 7)])

plot_96_profile(mut_mtx[,c(1, 7)], condensed=TRUE)

counts <- snp_count
colnames(counts) <- substr(colnames(counts), 3, 5)
counts <- t(counts)
types <- levels(factor(row.names(counts)))
lst <- list()
for (i in 1:length(types)) {
    lst[[i]] <- data.frame(
        type=types[i],
        freq=sum(counts[row.names(counts) == types[i],])
    )
}
df <- do.call(rbind, lst)
df$pct <- round(df$freq / sum(df$freq) * 100, 2)
df$type2 <- c(
    "C>A/A>G", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C"
)
df
labels <- paste0(df$type2, "\n(", df$pct, "%)")
colors <- c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
pie(df$pct, labels=labels, col=colors)
outfile <- "fig_1_1.tiff"
tiff(outfile, width=20, height=20, unit="cm", res=350)
pie(df$pct, labels=labels, col=colors)
dev.off()


library(NMF)
mut_mtx2 <- mut_mtx + 0.0001
dim(mut_mtx2)
save.image("project_4_1.RData")


library(MutationalPatterns)
res = nmf(mut_mtx2, rank=4, method="brunet", nrun=30, seed=123456, .opt='p7')

signatures = NMF::basis(res)
contribution = NMF::coef(res)

reconstructed = signatures %*% contribution
nmf_res <- list(signatures = signatures,
			contribution = contribution,
			reconstructed = reconstructed)
str(nmf_res)

colnames(nmf_res$signatures) <- paste0("Signature", LETTERS[1:4])
rownames(nmf_res$contribution) <- paste0("Signature", LETTERS[1:4])


cancer_signatures = read.table("signatures_probabilities.txt", sep="\t", header=TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
colnames(cancer_signatures)
cancer_signatures = as.matrix(cancer_signatures[,4:33])
all(rownames(cancer_signatures) == rownames(nmf_res$signatures))
cos_sim_signatures = cos_sim_matrix(nmf_res$signatures, cancer_signatures)


plot_cosine_heatmap(cos_sim_signatures, cluster_rows=FALSE, plot_values=TRUE)
outfile <- "fig_cosine_heatmap.tiff"
tiff(outfile, width=30, height=20, unit="cm", res=350, compression="lzw+p")
plot_cosine_heatmap(cos_sim_signatures, cluster_rows=FALSE, plot_values=TRUE)
dev.off()

outfile <- "fig_cosine_heatmap.pdf"
pdf(outfile, onefile=FALSE, width=30/2.54, height=20/2.54)
plot_cosine_heatmap(cos_sim_signatures, cluster_rows=FALSE, plot_values=TRUE)
dev.off()

save.image("project_4.RData")
save(nmf_res, file="nmf_res_rank4.RData")




options(stringsAsFactors=FALSE)
library(NMF)
load("mut_mat.RData")
load("nmf_res.RData")
load("mut_mtx.RData")

library(MutationalPatterns)
plot_96_profile(nmf_res$signatures)
outfile <- "fig_barplot_5_1_test.tiff"
tiff(outfile, width=20, height=20, unit="cm", res=350, compression="lzw+p")
plot_96_profile(nmf_res$signatures)
dev.off()



mut_matrix <- nmf_res[["signatures"]]

norm_mut_matrix = apply(mut_matrix, 2, function(x) x / sum(x))

CONTEXTS_96 <- rownames(mut_matrix)
context = CONTEXTS_96
SUBSTITUTIONS <- c(
    "C>A",
    "C>G",
    "C>T",
    "T>A",
    "T>C",
    "T>G"
)
substitution = rep(SUBSTITUTIONS, each=16)

context <- paste0(substr(context, 1, 1), ".", substr(context, 7, 7))

df = data.frame(substitution=substitution, context=context)
rownames(norm_mut_matrix) = NULL
library(reshape2)
df2 = cbind(df, as.data.frame(norm_mut_matrix))
df3 = melt(df2, id.vars=c("substitution", "context"))
head(df3)


df4 <- df3
df4[, "value"] <- round(df4[, "value"] * 100, 1)
outdf <- df4
outfile <- "result_signature_96.csv"
write.csv(outdf, outfile, row.names=FALSE, na="")



library(ggplot2)
library(grid)
library(gtable)
library(scales)

colors <- c(
rgb(204,0,255, max=255),
rgb(255,153,0, max=255),
rgb(51,255,0, max=255),
rgb(197,51,115, max=255),
rgb(160,64,255, max=255),
rgb(255,0,153, max=255)
)
colors



colors_subgroup <- c(
rgb(204,0,255, max=255),
rgb(255,153,0, max=255),
rgb(51,255,0, max=255),
rgb(197,51,115, max=255),
rgb(160,64,255, max=255),
rgb(255,0,153, max=255)
)
colors_subgroup


ymax = 0.35
signatures <- c("SignatureA", "SignatureB", "SignatureC", "SignatureD")
for (i in 1:length(signatures)) {
df4 <- df3[df3[, "variable"] == signatures[i],]
df4[, "substitution2"] <- ""

gplot <- ggplot(
    data=df4,
    aes(
        x=context,
        y=value,
        fill=substitution,
        width=0.6),
    ) +
    geom_rect(
       aes(fill=substitution),
       alpha=0.01, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf
    ) +
    geom_bar(stat="identity", colour="transparent", size=0.2) +
    scale_fill_manual(values=colors) +
    facet_wrap(
        ~ substitution2 + substitution, strip.position="bottom", nrow=1
    ) +
    coord_cartesian(ylim=c(0, ymax)) +
    scale_y_continuous(
        "Relative contribution", breaks=seq(0, ymax, 0.05),
        labels=percent
    ) +
    guides(fill=FALSE) +
    theme_bw() +
    theme(
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_line(),
        legend.title=element_text(),
        legend.position="none",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", color="transparent"),
        legend.background=element_rect(fill="transparent", color="black"),

        strip.background=element_rect(color=NA, fill=NA),
        strip.placement="outside",
        panel.background=element_rect(fill=NA),

        panel.border=element_rect(color=NA),
        panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),

        panel.spacing.x=unit(0.06, "cm"),
        plot.background=element_blank()
    )

grob <- ggplotGrob(gplot)
strips <- which(grepl("strip-", grob[["layout"]][["name"]]))
for (j in 1:length(strips)) {
    k <- which(grepl("rect",
        grob[["grobs"]][[strips[j]]][["grobs"]][[1]][["childrenOrder"]])
    )
    grob[["grobs"]][[strips[j]]][["grobs"]][[1]][[
        "children"]][[k]][["gp"]][["fill"]] <- colors_subgroup[j]
}
grob1 <- grob
grob <- grobTree(
    rectGrob(gp=gpar(fill="black")),
    textGrob(signatures[i], gp=gpar(col="white"))
)
grob2 <- grob
gtab <- gtable(
    unit(c(16, 3, 1), "cm"),
    unit(c(1, 1, 8), "cm")
)

gtab <- gtable_add_grob(gtab, grob1, 1, 1, 3, 3)
gtab <- gtable_add_grob(gtab, grob2, 2, 2, 2, 2)
grid.newpage()
grid.draw(gtab)

outfile <- paste0("fig_barplot_5_1_", i, ".pdf")
pdf(outfile, onefile=FALSE, width=20/2.54, height=10/2.54)
grid.newpage()
grid.draw(gtab)
dev.off()
}


gplot <- ggplot(
    data=df3, aes(
        x=context,
        y=value,
        fill=substitution,
        width=0.6
    )) +
    geom_bar(stat="identity", colour="transparent", size=.2) +
    scale_fill_manual(values=colors) +
    facet_grid(variable ~ substitution) +
    coord_cartesian(ylim=c(0, ymax)) +
    scale_y_continuous("Relative contribution", breaks=seq(0, ymax, 0.1)) +
    guides(fill=FALSE) +
    theme_bw() +
    theme(
        axis.title.y=element_text(size=12,vjust=1),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x=element_text(size=9),
        strip.text.y=element_text(size=9),
        panel.grid.major.x = element_blank()
    )
gplot

outfile <- "fig_barplot_5_1.pdf"
pdf(outfile, onefile=FALSE, width=20/2.54, height=20/2.54)
gplot
dev.off()







options(stringsAsFactors=FALSE)
library(ggplot2)

library(NMF)
load("mut_mat.RData")
load("nmf_res.RData")
load("mut_mtx.RData")
library(MutationalPatterns)
plot_contribution(nmf_res$contribution, nmf_res$signature, mode="absolute")


signatures <- nmf_res[["signatures"]]
contribution <- nmf_res[["contribution"]]

total_signatures = colSums(signatures) 

abs_contribution = contribution * total_signatures

library(reshape2)
m_contribution = melt(abs_contribution)
colnames(m_contribution) = c("Signature", "Sample", "Contribution")
head(m_contribution)

sumdf <- aggregate(Contribution ~ Sample, m_contribution, sum)
ids <- sumdf[order(-sumdf[, "Contribution"]), "Sample"]
m_contribution[, "Sample"] <- factor(m_contribution[, "Sample"], levels=ids)
m_contribution <- m_contribution[order(m_contribution[, "Sample"]),]

m_contribution[, "Contribution"] <- log2(m_contribution[, "Contribution"] + 1)

colors <- rainbow(4)
fig <- ggplot(
    m_contribution, aes(
        x=factor(Sample),
        y=Contribution,
        fill=factor(Signature),
        color=factor(Signature),
        order=Sample
        )
    ) +
    geom_bar(stat="identity", colour = "transparent") +
    scale_color_manual(name="Signature", values=colors) +
    scale_fill_manual(name="Signature", values=colors) +
    scale_x_discrete("") +
    scale_y_continuous(
        "Log2(Absolute contribution + 1) \n Log2(no. mutations + 1)"
    ) +
    theme_bw() +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid=element_blank()
    )
fig


outfile <- "fig_barplot_5_2.pdf"
pdf(outfile, onefile=TRUE, width=20/2.54, height=10/2.54)
fig
dev.off()







options(stringsAsFactors=FALSE)
library(NMF)
library(reshape2)
library(ggplot2)
load("nmf_res.RData")
contribution <- nmf_res$contribution

clindf <- read.table("group.txt", header=TRUE, sep="\t", quote="")
dim(clindf)
head(clindf)
clindf[, "sample_id"] <- clindf[, "Mixture"]
dim(clindf)
clindf[1:6,1:3]

library(MutationalPatterns)
plot_contribution(nmf_res$contribution, nmf_res$signature, mode="relative")


m_contribution = melt(contribution)
colnames(m_contribution) = c("Signature", "Sample", "Contribution")
aggregate(Contribution ~ Signature, m_contribution, sum)

plot = ggplot(m_contribution,
                        aes(x = factor(Sample),
                            y = Contribution,
                            fill = factor(Signature),
                            order = Sample)) +
            geom_bar(position = "fill", stat="identity", colour="transparent")  +

            labs(x = "", y = "Relative contribution") +

            theme_bw() +

            theme(panel.grid.minor.x=element_blank(),
                    panel.grid.major.x=element_blank()) +
            theme(panel.grid.minor.y=element_blank(),
                    panel.grid.major.y=element_blank())
        palette <- rainbow(4)
        plot = plot + scale_fill_manual(name="Signature", values=palette)
        plot = plot + xlim(levels(factor(m_contribution$Sample)))
plot

outfile <- "fig_barplot_5_3_test.pdf"
pdf(outfile, onefile=TRUE, width=30/2.54, height=10/2.54)
plot
dev.off()


library(ggplot2)
library(grid)
library(gtable)
library(ggsci)
newdf1 <- nmf_res[["contribution"]]
newlst1 <- list()
for (i in 1:ncol(newdf1)) {
    newlst1[[i]] <- data.frame(
        sample_id=colnames(newdf1)[i],
        signature=rownames(newdf1),
        contribution=newdf1[, i],
        percent=newdf1[, i] / sum(newdf1[, i]) * 100
    )
}
newdf2 <- do.call(rbind, newlst1)
newdf3 <- merge(newdf2, clindf, by="sample_id", all=FALSE)
figdf <- newdf3
head(figdf)


figdf[, "id"] <- figdf[, "sample_id"]
figdf[, "value"] <- figdf[, "percent"]
figdf[, "group"] <- figdf[, "signature"]
figdf[, "subgroup"] <- figdf[, "risk"]


df4 <- aggregate(value ~ group, figdf, max)
df4[, "value"] <- round(df4[, "value"], 1)
colnames(df4)[2] <- "max"

df5 <- aggregate(value ~ group, figdf, min)
df5[, "value"] <- round(df5[, "value"], 1)
colnames(df5)[2] <- "min"
outdf <- merge(df4, df5, by="group", all=TRUE)
outfile <- "result_signature_sample.csv"
write.csv(outdf, outfile, row.names=FALSE, na="")


sort_group <- c("left", "right", "center")[3]
sort_subgroup <- c("left", "right", "center")[3]


colors_group <- c(
rgb(187,87,198, max=255),
rgb(110,110,207, max=255),
rgb(187,204,73, max=255),
rgb(227,68,90, max=255),
rgb(215,97,137, max=255),
rgb(197,51,115, max=255)
)
barplot(1:6,col=colors_group)
colors_group <- c(
rgb(187,87,198, max=255),
rgb(187,204,73, max=255),
rgb(215,97,137, max=255),
rgb(197,51,115, max=255)
)
colors_group

barplot(1:6,col=colors_group)

colors_subgroup <- pal_jco("default")(2)
colors_subgroup <- c(
rgb(86, 86, 220, max=255),
rgb(224, 77, 224, max=255)
)
barplot(1:2,col=colors_subgroup)

r <- 0.5

title_group <- "Mutation"
title_subgroup <- "Risk"

spacing_subgroup <- unit(0.1 * r, "cm")

width_barplot <- 16.1 * 2.54 * r
width_group <- 2.3 * 2.54 * r
width_subgroup <- 1.6 * 2.54 * r
width <- width_barplot + width_group + width_subgroup
height <- 10 * 2.54 * r

top_group <- 1 * r
top_subgroup <- 1 * r

if (sort_subgroup == "left") {
    tb <- table(figdf[, "subgroup"])
    subgroups <- names(tb[order(-tb)])
    figdf[, "subgroup"] <- factor(figdf[, "subgroup"], levels=subgroups)
} else if (sort_subgroup == "right") {
    tb <- table(figdf[, "subgroup"])
    subgroups <- names(tb[order(tb)])
    figdf[, "subgroup"] <- factor(figdf[, "subgroup"], levels=subgroups)
} else if (sort_subgroup == "center") {
    tb <- table(figdf[, "subgroup"])
    subgroups <- names(tb[order(-tb)])
    subgroups2 <- c()
    for (i in 1:length(subgroups)) {
        if (i %% 2 == 1) {
            subgroups2 <- c(subgroups[i], subgroups2)
        } else {
            subgroups2 <- c(subgroups2, subgroups[i])
        }
    }
    figdf[, "subgroup"] <- factor(figdf[, "subgroup"], levels=subgroups2)
}

newdf1 <- aggregate(value ~ group, figdf, mean)
newdf1

groups <- newdf1[order(newdf1[, "value"]), "group"]
figdf[, "group"] <- factor(figdf[, "group"], levels=groups)
if (sort_group == "left") {
    figdf <- figdf[
        order(figdf[, "subgroup"], figdf[, "group"], -figdf[, "value"]),]
    ids <- figdf[figdf[, "group"] == groups[length(groups)], "id"]
    figdf[, "id"] <- factor(figdf[, "id"], levels=ids)
} else if (sort_group == "right") {
    figdf <- figdf[
        order(figdf[, "subgroup"], figdf[, "group"], figdf[, "value"]),]
    ids <- figdf[figdf[, "group"] == groups[length(groups)], "id"]
    figdf[, "id"] <- factor(figdf[, "id"], levels=ids)
} else if (sort_group == "center") {
    figdf <- figdf[
        order(figdf[, "subgroup"], figdf[, "group"], -figdf[, "value"]),]
    ids <- figdf[figdf[, "group"] == groups[length(groups)], "id"]
    ids2 <- c()
    for (i in 1:length(ids)) {
        if (i %% 2 == 1) {
            ids2 <- c(ids[i], ids2)
        } else {
            ids2 <- c(ids2, ids[i])
        }
    }
    figdf[, "id"] <- factor(figdf[, "id"], levels=ids2)
}
head(figdf)


p <- ggplot(figdf, aes(x=id, y=value, fill=group)) +
    geom_bar(position="stack", stat="identity", width=1.1) +
    scale_x_discrete("") +
    scale_y_continuous("", expand=c(0, 0), limit=c(0, 100)) +
    scale_fill_manual(values=colors_group) +
    facet_grid(~ subgroup, scales="free", space="free") +
    guides(fill=guide_legend(ncol=1, byrow=TRUE)) +
    theme(
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_line(),
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
        panel.spacing.x=unit(spacing_subgroup, "cm"),
        plot.background=element_blank()
    )
p1 <- p

p <- ggplot(figdf, aes(x=id, y=value, fill=group)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_manual(values=colors_group) +
    guides(fill=guide_legend(title=title_group, ncol=1, byrow=TRUE)) +
    theme(
        legend.title=element_text(),
        legend.position="right",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black")
    )
p2 <- p

p <- ggplot(figdf, aes(x=id, y=value, fill=subgroup)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_manual(values=colors_subgroup) +
    guides(fill=guide_legend(title=title_subgroup, ncol=1, byrow=TRUE)) +
    theme(
        legend.title=element_text(),
        legend.position="right",
        legend.justification=c(0, 1),
        legend.key=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="black")
    )
p3 <- p


g <- ggplotGrob(p1)
strips <- which(grepl('strip-', g[["layout"]][["name"]]))
for (i in 1:length(strips)) {
    j <- which(grepl("rect",
        g[["grobs"]][[strips[i]]][["grobs"]][[1]][["childrenOrder"]])
    )
    g[["grobs"]][[strips[i]]][["grobs"]][[1]][[
        "children"]][[j]][["gp"]][["fill"]] <- colors_subgroup[i]
}
g1 <- g

g <- ggplotGrob(p2)
guide <- which(sapply(g[["grobs"]], function(x) x$name) == "guide-box")
g2 <- g[["grobs"]][[guide]]
g2[["heights"]][2] <- unit(top_group, "cm")

g <- ggplotGrob(p3)
guide <- which(sapply(g[["grobs"]], function(x) x$name) == "guide-box")
g3 <- g[["grobs"]][[guide]]
g3[["heights"]][2] <- unit(top_subgroup, "cm")

gt <- gtable(
    unit(c(width_barplot, width_group, width_subgroup), c("cm")),
    unit(height, "cm"))

gt <- gtable_add_grob(gt, g1, 1, 1)
gt <- gtable_add_grob(gt, g2, 1, 2, 1, 2)
gt <- gtable_add_grob(gt, g3, 1, 3, 1, 3)

outfile <- "fig_barplot_5_3.pdf"
pdf(outfile, onefile=TRUE, width=width/2.54, height=height/2.54)
grid.draw(gt)
dev.off()







options(stringsAsFactors=FALSE)

library(NMF)
load("mut_mat.RData")
load("nmf_res.RData")
load("mut_mtx.RData")


clindf <- read.table("group.txt", header=TRUE, sep="\t", quote="")
dim(clindf)
head(clindf)
clindf[, "sample_id"] <- clindf[, "Mixture"]
dim(clindf)
clindf[1:6,1:3]


signatures <- nmf_res[["signatures"]]
contribution <- nmf_res[["contribution"]]

total_signatures = colSums(signatures) 

abs_contribution = contribution * total_signatures

library(reshape2)
m_contribution = melt(abs_contribution)
colnames(m_contribution) = c("Signature", "Sample", "Contribution")
head(m_contribution)
newdf2 <- m_contribution
newdf2[, "sample_id"] <- newdf2[, "Sample"]
newdf2 <- merge(newdf2, clindf, by="sample_id", all=FALSE)
newdf2[, "value"] <- log2(newdf2[, "Contribution"] + 1)
newdf2[, "signature"] <- newdf2[, "Signature"]
newdf2[, "group"] <- newdf2[, "risk"]
nrow(newdf2)


signatures <- c("SignatureA", "SignatureB", "SignatureC", "SignatureD")
newlst1 <- list()
for (i in 1:length(signatures)) {
    newdf3 <- newdf2[newdf2[, "signature"] %in% signatures[i],]
    newdf3[, "group"] <- factor(newdf3[, "group"])
    model <- wilcox.test(value ~ group, newdf3)
    print(model)
    means <- aggregate(value ~ factor(group), newdf3, mean)
    medians <- aggregate(value ~ factor(group), newdf3, median)
    maxs <- aggregate(value ~ factor(group), newdf3, max)
    mins <- aggregate(value ~ factor(group), newdf3, min)
    newlst1[[i]] <- data.frame(
        var=signatures[i],
        mean_1=means[1, "value"],
        mean_2=means[2, "value"],
        median_1=medians[1, "value"],
        median_2=medians[2, "value"],
        max_1=maxs[1, "value"],
        max_2=maxs[2, "value"],
        min_1=mins[1, "value"],
        min_2=mins[2, "value"],
        pvalue=model[["p.value"]])
}
newdf4 <- do.call(rbind, newlst1)
newdf4[, "adj_pvalue"] <- p.adjust(newdf4[, "pvalue"], method="bonferroni")
write.csv(newdf4, "result_mean_2.csv", row.names=FALSE, na="")


library(ggplot2)
library(ggpubr)
library(scales)
windowsFonts()

colors1 <- c(
    rgb(255,0,153, max=255),
    rgb(51,255,0, max=255),
    rgb(255,153,0, max=255),
    rgb(204,0,255, max=255)
)
colors1

colors2 <- c(
    rgb(153,0,92, max=255),
    rgb(31,153,0, max=255),
    rgb(153,92,0, max=255),
    rgb(122,0,153, max=255)
)
colors2

maxvalues <- apply(newdf4[, c("max_1", "max_2")], 1, max)
minvalues <- apply(newdf4[, c("min_1", "min_2")], 1, min)

for (i in 1:length(signatures)) {
newdf3 <- newdf2[newdf2[, "signature"] %in% signatures[i],]

colors <- c(colors1[i], colors2[i])
titlex <- signatures[i]
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

        fatten=0.9
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

outfile <- paste0("fig_boxplot_", i, "_", ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_boxplot_", i, "_", ".pdf")
pdf(outfile, width=4/2.54, height=5/2.54)
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

outfile <- paste0("fig_violin_", i, "_", ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_violin_", i, "_", ".pdf")
pdf(outfile, width=4/2.54, height=5/2.54)
print(fig)
dev.off()
}









# Oncogenic pathways
inputFile="EBPlusPlusAdjustPANCAN_exp_log.txt"                        
gmtFile="pathway.gmt"

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

options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)






clindf <- read.table("pathway.txt", head=TRUE, sep="\t", quote="")
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
colors1 <- rainbow(10, s=1, v=1, alpha = 1)
barplot(1:10,col=rainbow(10, s=1, v=1, alpha = 1))
colors1

colors2 <- rainbow(10, s=1, v=0.6, alpha = 1)
barplot(1:10,col=rainbow(10, s=1, v=0.6, alpha = 1))
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
pdf(outfile, width=3/2.54, height=5/2.54)
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

outfile <- paste0("fig_violin_", i, "_", levels[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig_violin_", i, "_", levels[i], ".pdf")
pdf(outfile, width=3/2.54, height=5/2.54)
print(fig)
dev.off()
}






# MHC molecules, co-stimulators and co-inhibitors
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
