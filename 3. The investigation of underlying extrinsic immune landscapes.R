# 3.1 Sankey diagram

options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggalluvial)
library(grid)
library(gtable)


infile <- "Sankey.txt"
clindf <- read.table(infile, header=TRUE, sep="\t", quote="", check.names=FALSE)
dim(clindf)
clindf[1:6, 1:10]
table(clindf[, "acronym"])
length(table(clindf[, "acronym"]))

clindf <- clindf[complete.cases(clindf[, 7:ncol(clindf)]),]
dim(clindf)
clindf[1:6, 1:10]
table(clindf[, "acronym"])
length(table(clindf[, "acronym"]))
names(table(clindf[, "acronym"]))


newdf1 <- clindf
newlst1 <- list()
for (i in 7:ncol(newdf1)) {
    newdf2 <- data.frame(
        id=rownames(newdf1),
        dataset=newdf1[, "acronym"],
        group=colnames(newdf1)[i],
        risk=newdf1[, "risk"],
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

newdf3[, "group"] <- factor(newdf3[, "group"], levels=groups)
newdf3[, "dataset"] <- factor(newdf3[, "dataset"], levels=datasets)
newdf3 <- newdf3[
    order(newdf3[, "dataset"], newdf3[, "group"], -newdf3[, "value"]),]
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

risks <- c("high", "low")
colors1 <- c(
rgb(145,34,30, max=255),
rgb(202,109,107, max=255),
rgb(143,51,144, max=255),
rgb(212,128,213, max=255),
rgb(86,54,143, max=255),
rgb(187,87,198, max=255),
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
rgb(227,68,90, max=255),
rgb(215,97,137, max=255),
rgb(197,51,115, max=255)
)
barplot(1:length(colors1),col=colors1)


colors2 <- c(
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

colors3 <- c(rgb(152,152,220, max=255), rgb(224,167,224, max=255))
colors3

newdf3[, "dataset2"] <- newdf3[, "dataset"]
library(ggplot2)
library(grid)
library(gtable)
p <- ggplot(newdf3, aes(x=id, y=value, fill=group)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_x_discrete("") +
    scale_y_continuous("", expand=c(0, 0)) +
    scale_fill_manual(values=colors1) +
    facet_grid(~ dataset + dataset2, scales="free", space="free") +
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
        strip.text.x=element_text(color="black", angle=90, hjust=0, size=6),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),
        panel.spacing.x=unit(0, "lines"),
        plot.background=element_blank(),
        plot.margin=margin(unit(c(2, 2, 2, 2), "cm")))
p1 <- p

outfile <- "fig_bar.pdf"
pdf(outfile, onefile=FALSE, height=28/2.54, width=40/2.54)
p1
dev.off()

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
clindf2[, "acronym"] <- factor(clindf2[, "acronym"], levels=rev(datasets))
clindf2[, "risk"] <- factor(clindf2[, "risk"], levels=rev(risks))
p <- ggplot(
    clindf2, aes(axis1=risk, axis2=acronym)) +
    geom_alluvium(aes(fill=acronym), width=1/12, alpha=1) +
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
    g$grobs[[strips[i]]]$grobs[[2]]$children[[j]]$gp$fill <- colors2[i]
    g$grobs[[strips[i]]]$grobs[[2]]$children[[j]]$y <- unit(0, "npc")
    g$grobs[[strips[i]]]$grobs[[2]]$children[[j]]$height <- unit(1.5, "npc")
    g$grobs[[strips[i]]]$grobs[[1]]$children[[j]]$gp$fill <- "white"
    g$grobs[[strips[i]]]$grobs[[1]]$children[[j]]$height <- unit(1, "npc")
    g$grobs[[strips[i]]]$grobs[[1]]$children[[j]]$y <- unit(0, "npc")
    g$grobs[[strips[i]]]$grobs[[1]]$children[[j]]$height <- unit(0, "npc")

    k <- which(grepl("text", g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[strips[i]]]$grobs[[2]]$children[[k]]$children[[1]]$label <- ""
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

t <- gtable(unit(38, "cm"), unit(c(10, 10, 3.3, 2, 1.1), "cm"))
t <- gtable_add_grob(t, g1, 1, 1, 1, 1)
t <- gtable_add_grob(t, g5, 2, 1, 2, 1)
t <- gtable_add_grob(t, g2, 3, 1, 3, 1)
t <- gtable_add_grob(t, g3, 4, 1, 4, 1)
t <- gtable_add_grob(t, g4, 5, 1, 5, 1)
grid.draw(t)

outfile <- "fig_sankey_bar_risk.pdf"
pdf(outfile, onefile=FALSE, height=28/2.54, width=40/2.54)
grid.draw(t)
dev.off()















# 3.2 TIL fraction, leukocyte fraction and lymphocyte fraction analyses

options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)

clindf <- read.table("risk2.txt", head=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6,1:5]
clindf[, "group"] <- factor(clindf[, "risk"])

clindf <- clindf[order(clindf[, "id"]),]
clindf[, "sample_id"] <- substr(clindf[, "id"], 1, 12)
clindf[, "sample_id2"] <- clindf[, "id"]
dim(clindf)
head(clindf)
colnames(clindf)

sum(duplicated(clindf[, "sample_id"]))
clindf <- clindf[!duplicated(clindf[, "sample_id"]),]
sum(duplicated(clindf[, "sample_id"]))


expdf <- read.csv(
    "The Immune Landscape of Cancer-mmc2.csv", head=TRUE, quote=""
)
dim(expdf)
expdf[1:6, 1:5]
expdf[, "sample_id"] <- expdf[, "TCGA.Participant.Barcode"]

sum(duplicated(clindf[, "sample_id"]))
sum(duplicated(expdf[, "sample_id"]))

sample_ids <- intersect(clindf[, "sample_id"], expdf[, "sample_id"])
length(sample_ids)
sample_ids[1:10]
clindf2 <- merge(clindf, expdf, by="sample_id", all=FALSE)
dim(clindf2)
colnames(clindf2)

write.csv(clindf2, "clindf2.csv", row.names=FALSE, na="")

newdf1 <- clindf2
dim(newdf1)
ncol(newdf1)

vars <- colnames(newdf1)[10:ncol(newdf1)]
vars
outlst <- list()
for (i in 1:length(vars)) {
    print(paste0("------- ", vars[i], " --------"))
    newdf2 <- data.frame(
        value=newdf1[, vars[i]],
        group=newdf1[, "group"]
    )
    newdf2 <- newdf2[complete.cases(newdf2),]
    n <- nrow(newdf2)
    model <- wilcox.test(value ~ factor(group), newdf2)
    print(model)
    means <- aggregate(value ~ factor(group), newdf2, mean)
    medians <- aggregate(value ~ factor(group), newdf2, median)
    maxs <- aggregate(value ~ factor(group), newdf2, max)
    mins <- aggregate(value ~ factor(group), newdf2, min)
    outlst[[i]] <- data.frame(
        var=vars[i],
        n=n,
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


dir.create("fig")

windowsFonts()

library(scales)
colors1 <- rainbow(length(vars), s=1, v=1, alpha = 1)
barplot(1:12,col=rainbow(12, s=1, v=1, alpha = 1))
colors1

colors2 <- rainbow(length(vars), s=1, v=0.6, alpha = 1)
barplot(1:12,col=rainbow(12, s=1, v=0.6, alpha = 1))
colors2


maxvalues <- apply(outdf[, c("max_high", "max_low")], 1, max)
minvalues <- apply(outdf[, c("min_high", "min_low")], 1, min)
for (i in 1:length(vars)) {
newdf2 <- data.frame(value=newdf1[, vars[i]], group=newdf1[, "group"])
newdf2 <- newdf2[complete.cases(newdf2),]
newdf2[, "group"] <- factor(
    newdf2[, "group"], levels=c("high", "low"), labels=c("INB", "IB")
)
windowsFonts()
colors <- c(colors1[i], colors2[i])
title_x <- vars[i]
tb <- table(newdf2[, "group"])
title_label_x <- paste0(names(tb), " (n=", tb, ")")
comparisons <- list(c("INB", "IB"))
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
    scale_x_discrete(title_x, labels=title_label_x) +
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
        method="t.test",
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

outfile <- paste0("fig/", "fig_boxplot_", i, "_", vars[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig/", "fig_boxplot_", i, "_", vars[i], ".pdf")
pdf(outfile, width=2/2.54, height=5/2.54)
print(fig)
dev.off()



fig <- ggplot(newdf2, aes(x=group, y=value, group=group)) +
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
    scale_x_discrete(title_x, labels=title_label_x) +
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

outfile <- paste0("fig/", "fig_violin_", i, "_", vars[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig/", "fig_violin_", i, "_", vars[i], ".pdf")
pdf(outfile, width=3/2.54, height=5/2.54)
print(fig)
dev.off()
}













# 3.3 Danaher immune infiltration analysis

library(ggplot2)
library(ggpubr)
library(ggsci)

infile <- "Danaher_immune.txt"
clindf <- read.table(infile, head=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6, 1:6]

newlst1 <- list()
vars <- colnames(clindf)[6:ncol(clindf)]
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
pvalues <- newdf4[, "pvalue"]
pvalue_labels <- ifelse(
    pvalues < 0.001, "P < 0.001",
    paste0("P = ", sprintf("%.3f", pvalues))
)


dir.create("fig")


windowsFonts()
library(scales)
colors1 <- rainbow(29, s=0.3, v=1, alpha=1)
barplot(1:29,col=rainbow(29, s=0.3, v=1, alpha=1))

colors1

colors2 <- rainbow(29, s=1, v=1, alpha=1)
barplot(1:29,col=rainbow(29, s=1, v=1, alpha=1))
colors2

maxvalues <- apply(newdf4[, c("max_high", "max_low")], 1, max)
minvalues <- apply(newdf4[, c("min_high", "min_low")], 1, min)
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

outfile <- paste0("fig/fig_boxplot_", i, "_", levels[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig/fig_boxplot_", i, "_", levels[i], ".pdf")
pdf(outfile, width=2/2.54, height=4/2.54)
print(fig)
dev.off()
}



dir.create("fig2")


colors <- c("dodgerblue4", "dodgerblue")
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)
fig <- ggplot(newdf2, aes(x=group, y=value, group=group)) +
    geom_boxplot(
        aes(colour=group), notch=FALSE, outlier.size=0.6,
        size=1, fatten=1
    ) +

    scale_colour_manual(values=colors, labels=c("cold tumor", "hot tumor")) +
    scale_x_discrete("") +
    scale_y_continuous(
        "Value", limit=c(min(minvalues) * 1.2, max(maxvalues) * 1.2)
    ) +
    stat_compare_means(
        aes(
            label=ifelse(
                p < 0.001, "P < 0.001", paste0("P = ", ..p.format..)
            )
        ), label.x=1.2, label.y=max(maxvalues) * 1.15, family="sans"
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
        axis.text.y=element_text(color="black", family="sans"),
        axis.line.y=element_line(),
        legend.title=element_blank(),
        legend.position="top",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA),
        strip.text.x=element_text(
            color="black", angle=0, hjust=0.5, family="sans"
        ),
        strip.background=element_rect(colour=NA, fill=NA),
        strip.placement="outside"
    )
fig

outfile <- "fig2/fig_boxplot_all.tiff"
tiff(outfile, width=30, height=30, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- "fig2/fig_boxplot_all.pdf"
pdf(outfile, width=30/2.54, height=30/2.54)
print(fig)
dev.off()


colors1 <- rainbow(length(vars), s=0.4, v=1, alpha=1)
colors2 <- rainbow(length(vars), s=1, v=1, alpha=1)
colors <- unlist(
    lapply(1:length(vars), function(x) c(colors1[x], colors2[x]))
)
newdf3 <- newdf2
newdf3[, "id"] <- paste0(newdf3[, "cell"], "_", newdf3[, "group"])
newdf3[, "id"] <- factor(
    newdf3[, "id"],
    levels=paste0(rep(vars, each=2), "_", rep(c("high", "low"), length(vars)))
)
newdf3[, "id"] <- as.numeric(newdf3[, "id"])
newdf3[, "id"] <- ifelse(
    newdf3[, "id"] %% 2 == 1, newdf3[, "id"] + 0.1, newdf3[, "id"]
)
fig <- ggplot() +
    geom_boxplot(
        data=newdf3, aes(x=id, y=value, group=id),
        color=colors, notch=FALSE, outlier.size=0.6, size=1, fatten=1
    ) +
    scale_color_manual(values=colors1, labels=c("cold tumor", "hot tumor")) +
    scale_x_continuous(
        "", breaks=seq(2, length(vars)*2, 2) - 0.45, labels=vars, expand=c(0.01, 0.01)
    ) +
    scale_y_continuous(
        "Value",
        limit=c(min(minvalues) * 1.08, max(maxvalues) * 1.08),
        expand=c(0.01, 0.01)
    ) +
    theme_bw() +
    theme(
        axis.text.x=element_text(family="sans", angle=70, hjust=1),
        axis.text.y=element_text(family="sans"),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="top",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
texts <- ifelse(
    pvalues < 0.0001, "****",
    ifelse(
        pvalues < 0.001, "***",
        ifelse(
            pvalues < 0.01, "**",
            ifelse(pvalues < 0.05, "*", "")
        )
    )
)

labdf <- data.frame(
    x=1:length(vars) * 2 - 0.45, y=max(maxvalues) * 1.06, label=pvalue_labels
)
segdf1 <- data.frame(
    x=1:length(vars) * 2 - 0.45 - 0.46, xend=1:length(vars) * 2 - 0.45 + 0.46,
    y=max(maxvalues) * 1.02, yend=max(maxvalues) * 1.02
)
segdf2 <- data.frame(
    x=1:length(vars) * 2 - 0.45 - 0.45, xend=1:length(vars) * 2 - 0.45 - 0.45,
    y=max(maxvalues) * 1.01, yend=max(maxvalues) * 1.02
)
segdf3 <- data.frame(
    x=1:length(vars) * 2 - 0.45 + 0.45, xend=1:length(vars) * 2 - 0.45 + 0.45,
    y=max(maxvalues) * 1.01, yend=max(maxvalues) * 1.02
)
fig <- fig +

    geom_text(data=labdf, aes(x=x, y=y, label=label), size=2) +

    geom_segment(data=segdf1, aes(x=x, xend=xend, y=y, yend=yend), size=0.2) +
    geom_segment(data=segdf2, aes(x=x, xend=xend, y=y, yend=yend), size=0.2) +
    geom_segment(data=segdf3, aes(x=x, xend=xend, y=y, yend=yend), size=0.2)
fig

outfile <- "fig2/fig_boxplot_all2.pdf"

pdf(outfile, width=40/2.54, height=12/2.54)
print(fig)
dev.off()

outfile <- "fig2/fig_boxplot_all2.tiff"
tiff(outfile, width=30, height=20, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()















# 3.4 Twenty-nine immune signature evaluations
inputFile="EBPlusPlusAdjustPANCAN_exp.txt" 
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

geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())


ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)







# 3.5 Comparison of the 29 immune signatures estimated by the ssGSEA method based on RNA-sequencing data between the high-risk and low-risk groups.

options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggpubr)
library(ggsci)

infile <- "ssGSEA_immune.txt"
clindf <- read.table(infile, head=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6, 1:6]


newlst1 <- list()
vars <- colnames(clindf)[3:ncol(clindf)]
vars

for (i in 1:length(vars)) {
    newlst1[[i]] <- data.frame(
        cell=vars[i], group=clindf[, "type"], value=clindf[, vars[i]]
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
pvalues <- newdf4[, "pvalue"]
pvalue_labels <- ifelse(
    pvalues < 0.001, "P < 0.001",
    paste0("P = ", sprintf("%.3f", pvalues))
)


dir.create("fig")

windowsFonts()
library(scales)
colors1 <- rainbow(29, s=0.3, v=1, alpha=1)
barplot(1:29,col=rainbow(29, s=0.3, v=1, alpha=1))

colors1


colors2 <- rainbow(29, s=1, v=1, alpha=1)
barplot(1:29,col=rainbow(29, s=1, v=1, alpha=1))
colors2

maxvalues <- apply(newdf4[, c("max_high", "max_low")], 1, max)
minvalues <- apply(newdf4[, c("min_high", "min_low")], 1, min)
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

outfile <- paste0("fig/fig_boxplot_", i, "_", levels[i], ".tiff")
tiff(outfile, width=4, height=5, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- paste0("fig/fig_boxplot_", i, "_", levels[i], ".pdf")
pdf(outfile, width=2/2.54, height=4/2.54)
print(fig)
dev.off()
}



dir.create("fig2")

colors <- c("dodgerblue4", "dodgerblue")
comparisons <- list(c("high", "low"))
symnum.args <- list(
    cutpoints=c(0, 0.0001, 0.001, 0.01, 0.05, 1),
    symbols=c("", "", "", "", "")
)
fig <- ggplot(newdf2, aes(x=group, y=value, group=group)) +
    geom_boxplot(
        aes(colour=group), notch=FALSE, outlier.size=0.6,
        size=1, fatten=1
    ) +

    scale_colour_manual(values=colors, labels=c("cold tumor", "hot tumor")) +
    scale_x_discrete("") +
    scale_y_continuous(
        "Value", limit=c(min(minvalues) * 1.2, max(maxvalues) * 1.2)
    ) +
    stat_compare_means(
        aes(
            label=ifelse(
                p < 0.001, "P < 0.001", paste0("P = ", ..p.format..)
            )
        ), label.x=1.2, label.y=max(maxvalues) * 1.15, family="sans"
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
        axis.text.y=element_text(color="black", family="sans"),
        axis.line.y=element_line(),
        legend.title=element_blank(),
        legend.position="top",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA),
        strip.text.x=element_text(
            color="black", angle=0, hjust=0.5, family="sans"
        ),
        strip.background=element_rect(colour=NA, fill=NA),
        strip.placement="outside"
    )
fig

outfile <- "fig2/fig_boxplot_all.tiff"
tiff(outfile, width=30, height=30, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()

outfile <- "fig2/fig_boxplot_all.pdf"
pdf(outfile, width=30/2.54, height=30/2.54)
print(fig)
dev.off()

colors1 <- rainbow(length(vars), s=0.4, v=1, alpha=1)
colors2 <- rainbow(length(vars), s=1, v=1, alpha=1)
colors <- unlist(
    lapply(1:length(vars), function(x) c(colors1[x], colors2[x]))
)
newdf3 <- newdf2
newdf3[, "id"] <- paste0(newdf3[, "cell"], "_", newdf3[, "group"])
newdf3[, "id"] <- factor(
    newdf3[, "id"],
    levels=paste0(rep(vars, each=2), "_", rep(c("high", "low"), length(vars)))
)
newdf3[, "id"] <- as.numeric(newdf3[, "id"])
newdf3[, "id"] <- ifelse(
    newdf3[, "id"] %% 2 == 1, newdf3[, "id"] + 0.1, newdf3[, "id"]
)
fig <- ggplot() +
    geom_boxplot(
        data=newdf3, aes(x=id, y=value, group=id),
        color=colors, notch=FALSE, outlier.size=0.6, size=1, fatten=1
    ) +
    scale_color_manual(values=colors1, labels=c("cold tumor", "hot tumor")) +
    scale_x_continuous(
        "", breaks=seq(2, length(vars)*2, 2) - 0.45, labels=vars, expand=c(0.01, 0.01)
    ) +
    scale_y_continuous(
        "Value",
        limit=c(min(minvalues) * 1.08, max(maxvalues) * 1.08),
        expand=c(0.01, 0.01)
    ) +
    theme_bw() +
    theme(
        axis.text.x=element_text(family="sans", angle=70, hjust=1),
        axis.text.y=element_text(family="sans"),
        axis.line=element_line(),
        legend.title=element_blank(),
        legend.position="top",
        panel.grid=element_blank(),
        panel.border=element_rect(color=NA)
    )
texts <- ifelse(
    pvalues < 0.0001, "****",
    ifelse(
        pvalues < 0.001, "***",
        ifelse(
            pvalues < 0.01, "**",
            ifelse(pvalues < 0.05, "*", "")
        )
    )
)

labdf <- data.frame(
    x=1:length(vars) * 2 - 0.45, y=max(maxvalues) * 1.06, label=pvalue_labels
)
segdf1 <- data.frame(
    x=1:length(vars) * 2 - 0.45 - 0.46, xend=1:length(vars) * 2 - 0.45 + 0.46,
    y=max(maxvalues) * 1.02, yend=max(maxvalues) * 1.02
)
segdf2 <- data.frame(
    x=1:length(vars) * 2 - 0.45 - 0.45, xend=1:length(vars) * 2 - 0.45 - 0.45,
    y=max(maxvalues) * 1.01, yend=max(maxvalues) * 1.02
)
segdf3 <- data.frame(
    x=1:length(vars) * 2 - 0.45 + 0.45, xend=1:length(vars) * 2 - 0.45 + 0.45,
    y=max(maxvalues) * 1.01, yend=max(maxvalues) * 1.02
)
fig <- fig +

    geom_text(data=labdf, aes(x=x, y=y, label=label), size=2) +

    geom_segment(data=segdf1, aes(x=x, xend=xend, y=y, yend=yend), size=0.2) +
    geom_segment(data=segdf2, aes(x=x, xend=xend, y=y, yend=yend), size=0.2) +
    geom_segment(data=segdf3, aes(x=x, xend=xend, y=y, yend=yend), size=0.2)
fig


outfile <- "fig2/fig_boxplot_all2.pdf"

pdf(outfile, width=40/2.54, height=12/2.54)
print(fig)
dev.off()

outfile <- "fig2/fig_boxplot_all2.tiff"
tiff(outfile, width=30, height=20, unit="cm", res=350, compression="lzw+p")
print(fig)
dev.off()










# 3.6 Unsupervised clustering analysis
library(pheatmap)
rt=read.table("ssgseaOut2.txt",sep="\t",header=T,row.names=1,check.names=F)
dim(rt)
rt[1:6,1:6]

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

bk <- c(seq(-1.5,-0.1,by=0.01),seq(0,1.5,by=0.01))

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




# 3.7 The proportion of high immune infiltration and low immune infiltration estimated by 29 immune signatures in the high-risk and low-risk groups.

options(stringsAsFactors=FALSE)
library(ggplot2)

infile <- "stack_barplot.txt"
newdf1 <- read.table(infile, header=TRUE, sep="\t")
newdf1


newlst1 <- list()
for (i in 2:ncol(newdf1)) {
    newdf2 <- newdf1[, c(1, i)]
    colnames(newdf2) <- c("subgroup", "freq")
    newdf2[, "group"] <- colnames(newdf1)[i]
    newlst1[[i]] <- newdf2
}
newdf3 <- do.call(rbind, newlst1)
newdf3


newdf3[, "group"] <- factor(
    newdf3[, "group"], levels=c("high", "low"), labels=c("High", "Low")
)
newdf3[, "subgroup"] <- factor(
    newdf3[, "subgroup"], levels=c("hot", "cold"), labels=c("Hot", "Cold")
)
newdf3 <- newdf3[order(newdf3[, "group"], newdf3[, "subgroup"]),]
newdf3


source("./fig_stack_barplot.R")
fig_stack_barplot(
    newdf3, "freq", "group", "subgroup",
    file_output=TRUE,
    fig_title=NULL, 
    fig_subgroup_color=c("#E3191C", "#309EDE"), 
    fig_xlab_title="",
    fig_ylab_title="Percent(%)", 
    fig_ylab_breaks=seq(0, 1, 0.2), 
    fig_ylab_labels=seq(0, 100, 20), 
    fig_comparison=TRUE, 
    fig_filepath=".", 
    fig_filename="fig_barplot", 
    fig_filetype=c("tiff", "pdf"), 
    fig_width=10, 
    fig_height=10 
)










# 3.8 Volcano plot

options(stringsAsFactors=FALSE)
library(ggplot2)
library(ggrepel)
library(ggsci)

clindf <- read.table("Volcano.txt", head=TRUE, sep="\t", quote="")
dim(clindf)
head(clindf)
colnames(clindf)
sum(clindf[1, 3:ncol(clindf)])
sum(duplicated(substr(clindf[, 1], 1, 12)))


dim(clindf)
table(clindf[, "type"])
sum(table(clindf[, "type"]))
clindf <- clindf[!is.na(clindf[, "type"]),]
dim(clindf)
table(clindf[, "type"])
sum(table(clindf[, "type"]))


sum(!complete.cases(clindf[, "type"]))
eset<-eset[complete.cases(eset),]


library(limma)
table(clindf[, "type"])
sum(table(clindf[, "type"]))

group <- factor(clindf[, "type"])
eset <- t(clindf[, 3:ncol(clindf)])
dim(eset)
eset[1:5,1:6]
design <- model.matrix(~ 0 + group)
colnames(design) <- c("high", "low", "normal")
colnames(design)
dim(eset)
dim(design)
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












# 3.9 Correlation analysis

options(stringsAsFactors=FALSE)
immudf <- read.table("Correlation.txt", head=TRUE, sep="\t", quote="")
dim(immudf)
immudf[1:6,1:6]
table(immudf$type)

immudf <- immudf[!immudf[, "type"] %in% "normal",]
dim(immudf)
immudf[1:6,1:6]
table(immudf$type)

colnames(immudf)

immudf1 <- immudf[immudf[, "type"] == "high", 3:ncol(immudf)]
cormtx1 <- cor(immudf1, use="pairwise.complete.obs")

immudf2 <- immudf[immudf[, "type"] == "low", 3:ncol(immudf)]
cormtx2 <- cor(immudf2, use="pairwise.complete.obs")

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
colors1 <- c(
    rgb(219, 145, 235, max=255),
    rgb(133, 245, 173, max=255),
    rgb(255, 155, 243, max=255),
    rgb(255, 178, 115, max=255),
    "red4", "orange4", "blue4", "green4", "yellow4", "pink4",
    "cyan4", "darkorange4", "seagreen4"
)
colors1

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












#  3.10 Cytolytic activity score assessment

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















# 3.11 Fibroblast assessment

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















# 3.12 Comparison of expression patterns of chemokines between the high-risk and low-risk groups.

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

clindf <- read.table("分组信息.txt", header=TRUE, sep="\t", quote="")
dim(clindf)
clindf[1:6, 1:3]
colnames(clindf)[1] <- "sample_id"


t_expdf <- t(expdf[, 2:ncol(expdf)])
colnames(t_expdf) <- expdf[, "gene_id"]
t_expdf <- cbind(data.frame(sample_id=rownames(t_expdf)), t_expdf)

clindf2 <- merge(clindf, t_expdf, by="sample_id", all=FALSE)
dim(clindf2)
clindf2[1:6, 1:10]

genedf <- read.table("gene-2.txt", header=TRUE, sep="\t", quote="")
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
colnames(mtx) <- c("High", "Low")
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
