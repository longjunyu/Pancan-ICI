# Copy number profiles

options(stringsAsFactors=FALSE)
scores <- read.table(
    "gistic/cnv_results1/scores.gistic", head=TRUE, sep="\t", quote="")
head(scores)
unique(scores$Chromosome)

chrom_extract <- function(BSgenome.hg=NULL) {
    if (is.null(BSgenome.hg)) stop("NULL object !", call.=FALSE)
    obj <- list(
        species=GenomeInfoDb::organism(BSgenome.hg),
        genomebuild=BSgenome::providerVersion(BSgenome.hg))
    df <- data.frame(
        chrom=BSgenome::seqnames(BSgenome.hg),
        chrN=seq_along(BSgenome::seqnames(BSgenome.hg)),
        chr.length=GenomeInfoDb::seqlengths(BSgenome.hg),
        stringsAsFactors=FALSE)
    df <- df[1:24,]
    df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
    df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
    df$middle.chr <- round(diff(c(0, df$chr.length.sum)) / 2)
    df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
    obj$chromosomes <- df
    obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) {
        obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify=FALSE)
    obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) {
       obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify=FALSE)
    names(obj$chr2chrom) <- obj$chromosomes$chrN
    obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length),
        na.rm=TRUE)
    return(obj)
}

BSgenome.hg <- "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)

if ("23" %in% unique(scores$Chromosome)) {
    scores[scores$Chromosome==23, "Chromosome"] <- "X"
}
if ("24" %in% unique(scores$Chromosome)) {
    scores[scores$Chromosome==24, "Chromosome"] <- "Y"
}
chrID <- unname(
    unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))

scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

pdf("fig_cnv_gistic_score1.pdf", onefile=FALSE, width=40/2.54, height=10/2.54)

par(mfrow=c(1, 1), mar=par()$mar + c(0.5, 0.5, 0.5, 0.5))

scores.amp <- scores[scores$Type == "Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type == "Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp, scores.del)

seg.col <- list(
    gain="red", outscale.gain="darkred", loss="blue",
    outscale.red="midnightblue")
ylim <- c(min(scores$G.score) - 0.2, max(scores$G.score) + 0.2)

plot(scores.amp$Start.geno, scores.amp$G.score,
    pch=".", type='h', cex=2, xaxs="i", yaxs="i",
    xlim=c(0,chrom$genome.length), ylim=ylim,
    cex.main=1, ylab="gistic score", xlab=NA,

    cex.lab=1, col=adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8), xaxt="n", lwd=2, las=1)
    
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd=2,

    col=adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.8))

ink <- chrom$chromosomes$chrN %in% chrID
yrange <- abs(diff(ylim))
m.pos <- c(ylim[1] + 0.15, ylim[2] - 0.15)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) + 2
try(text(x=chrom$chromosomes$middle.chr.geno[ink], y=m.pos[m.mod],
    labels=chrom$chromosomes$chrom[ink], cex=1))

abline(h=0.0, col=1, lwd=1, lty=3)
abline(v=c(0, chrom$chromosomes$chr.length.sum), col=1, lty=3, lwd=1)

col1 <- adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8)
col2 <- adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.8)
legend("topleft", c("gain", "loss"), cex=0.6, bty="n", fill=c(col1, col2))
dev.off()

outfile <- "fig_cnv_percentage1.pdf"
pdf(outfile, onefile=FALSE, width=40/2.54, height=10/2.54)
scores.amp <- scores[scores$Type == "Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type == "Del",]
scores.del$frequency <- scores.del$frequency * -100

seg.col <- list(
    gain="red", outscale.gain="darkred", loss="blue",
    outscale.red="midnightblue")
ylim <- c(-90, 90)

plot(scores.amp$Start.geno, scores.amp$frequency,
    pch=".", type='h', cex=2, xaxs="i", yaxs="i", 
    xlim=c(0, chrom$genome.length), ylim=ylim,
    cex.main=1, ylab="gain/loss percentage in cohort", xlab=NA,

    cex.lab=1, col=adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8), xaxt="n", lwd=2, las=1)
    
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd=2,

    col=adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.5))
    
ink <- chrom$chromosomes$chrN %in% chrID
yrange <- abs(diff(ylim))
m.pos <- c(-80, 80)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) + 2
try(text(x=chrom$chromosomes$middle.chr.geno[ink], y=m.pos[m.mod],
    labels=chrom$chromosomes$chrom[ink], cex=1))

abline(h=0.0, col=1, lwd=1, lty=3)
abline(v=c(0, chrom$chromosomes$chr.length.sum), col=1, lty=3, lwd=1)

col1 <- adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8)
col2 <- adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.8)
legend("topleft", c("gain", "loss"), cex=0.6, bty="n", fill=c(col1, col2))
dev.off()









options(stringsAsFactors=FALSE)
scores <- read.table(
    "gistic/cnv_results2/scores.gistic", head=TRUE, sep="\t", quote="")
head(scores)
unique(scores$Chromosome)

chrom_extract <- function(BSgenome.hg=NULL) {
    if (is.null(BSgenome.hg)) stop("NULL object !", call.=FALSE)
    obj <- list(
        species=GenomeInfoDb::organism(BSgenome.hg),
        genomebuild=BSgenome::providerVersion(BSgenome.hg))
    df <- data.frame(
        chrom=BSgenome::seqnames(BSgenome.hg),
        chrN=seq_along(BSgenome::seqnames(BSgenome.hg)),
        chr.length=GenomeInfoDb::seqlengths(BSgenome.hg),
        stringsAsFactors=FALSE)
    df <- df[1:24,]
    df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
    df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
    df$middle.chr <- round(diff(c(0, df$chr.length.sum)) / 2)
    df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
    obj$chromosomes <- df
    obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) {
        obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify=FALSE)
    obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) {
       obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify=FALSE)
    names(obj$chr2chrom) <- obj$chromosomes$chrN
    obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length),
        na.rm=TRUE)
    return(obj)
}

BSgenome.hg <- "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)

if ("23" %in% unique(scores$Chromosome)) {
    scores[scores$Chromosome==23, "Chromosome"] <- "X"
}
if ("24" %in% unique(scores$Chromosome)) {
    scores[scores$Chromosome==24, "Chromosome"] <- "Y"
}
chrID <- unname(
    unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))

scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]


pdf("fig_cnv_gistic_score2.pdf", onefile=FALSE, width=40/2.54, height=10/2.54)

par(mfrow=c(1, 1), mar=par()$mar + c(0.5, 0.5, 0.5, 0.5))

scores.amp <- scores[scores$Type == "Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type == "Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp, scores.del)

seg.col <- list(

    gain=rgb(255, 129, 155, max=255), 

    outscale.gain=rgb(219, 145, 235, max=255), 

    gain=rgb(0, 197, 190, max=255), 

    outscale.loss=rgb(145, 162, 235, max=255)
    )
ylim <- c(min(scores$G.score) - 0.2, max(scores$G.score) + 0.2)

plot(scores.amp$Start.geno, scores.amp$G.score,
    pch=".", type='h', cex=2, xaxs="i", yaxs="i",
    xlim=c(0,chrom$genome.length), ylim=ylim,
    cex.main=1, ylab="gistic score", xlab=NA,

    cex.lab=1, col=adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8), xaxt="n", lwd=2, las=1)

lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd=2,

    col=adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.8))

ink <- chrom$chromosomes$chrN %in% chrID
yrange <- abs(diff(ylim))
m.pos <- c(ylim[1] + 0.15, ylim[2] - 0.15)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) + 2
try(text(x=chrom$chromosomes$middle.chr.geno[ink], y=m.pos[m.mod],
    labels=chrom$chromosomes$chrom[ink], cex=1))

abline(h=0.0, col=1, lwd=1, lty=3)
abline(v=c(0, chrom$chromosomes$chr.length.sum), col=1, lty=3, lwd=1)

col1 <- adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8)
col2 <- adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.8)
legend("topleft", c("gain", "loss"), cex=0.6, bty="n", fill=c(col1, col2))
dev.off()

outfile <- "fig_cnv_percentage2.pdf"
pdf(outfile, onefile=FALSE, width=40/2.54, height=10/2.54)
scores.amp <- scores[scores$Type == "Amp",]
scores.amp$frequency <- scores.amp$frequency * 100
scores.del <- scores[scores$Type == "Del",]
scores.del$frequency <- scores.del$frequency * -100

seg.col <- list(

    gain=rgb(255, 129, 155, max=255), 

    outscale.gain=rgb(219, 145, 235, max=255), 

    gain=rgb(0, 197, 190, max=255), 

    outscale.loss=rgb(145, 162, 235, max=255)
    )
ylim <- c(-90, 90)

plot(scores.amp$Start.geno, scores.amp$frequency,
    pch=".", type='h', cex=2, xaxs="i", yaxs="i", 
    xlim=c(0, chrom$genome.length), ylim=ylim,
    cex.main=1, ylab="gain/loss percentage in cohort", xlab=NA,

    cex.lab=1, col=adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8), xaxt="n", lwd=2, las=1)

lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd=2,

    col=adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.5))

ink <- chrom$chromosomes$chrN %in% chrID
yrange <- abs(diff(ylim))
m.pos <- c(-80, 80)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) + 2
try(text(x=chrom$chromosomes$middle.chr.geno[ink], y=m.pos[m.mod],
    labels=chrom$chromosomes$chrom[ink], cex=1))

abline(h=0.0, col=1, lwd=1, lty=3)
abline(v=c(0, chrom$chromosomes$chr.length.sum), col=1, lty=3, lwd=1)

col1 <- adjustcolor(rgb(219, 145, 235, max=255), alpha.f=.8)
col2 <- adjustcolor(rgb(145, 162, 235, max=255), alpha.f=.8)
legend("topleft", c("gain", "loss"), cex=0.6, bty="n", fill=c(col1, col2))
dev.off()























# Venn diagram

options(stringsAsFactors=FALSE)
infile <- "amp_genes.conf_99-IB.txt"
genedf <- read.table(
    infile, header=FALSE, sep="\t", fill=TRUE, quote="", check.names=FALSE
)
colnames(genedf)
genedf[1:6, 1:6]

genelst <- list()
for (i in 2:ncol(genedf)) {
    genes <- genedf[5:nrow(genedf), i]
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) > 0) {
        genelst[[i]] <- data.frame(
            cytoband=genedf[1, i],
            q_value=genedf[2, i],
            residual_qvalue=genedf[3, i],
            wide_peak_boundaries=genedf[4, i],
            genes_in_wide_peak=genes
        )
    }
}
genedf2 <- do.call(rbind, genelst)
genedf2[1:6,]
genedf_ib <- genedf2
genes_ib <- genedf2[, "genes_in_wide_peak"]

sum(duplicated(genes_ib))
genes_ib <- genes_ib[!duplicated(genes_ib)]

outdf <- genedf2
outfile <- "result_amp_genes.conf_99-IB.txt"
write.table(outdf, outfile, sep="\t", row.names=FALSE, quote=FALSE)
outfile <- "result_amp_genes.conf_99-IB.csv"
write.csv(outdf, outfile, row.names=FALSE, quote=FALSE)



infile <- "amp_genes.conf_99-INB.txt"
genedf <- read.table(
    infile, header=FALSE, sep="\t", fill=TRUE, quote="", check.names=FALSE
)
colnames(genedf)
genedf[1:6, 1:6]

genelst <- list()
for (i in 2:ncol(genedf)) {
    genes <- genedf[5:nrow(genedf), i]
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) > 0) {
        genelst[[i]] <- data.frame(
            cytoband=genedf[1, i],
            q_value=genedf[2, i],
            residual_qvalue=genedf[3, i],
            wide_peak_boundaries=genedf[4, i],
            genes_in_wide_peak=genes
        )
    }
}
genedf2 <- do.call(rbind, genelst)
genedf2[1:6,]
genedf_inb <- genedf2
genes_inb <- genedf2[, "genes_in_wide_peak"]

sum(duplicated(genes_inb))
genes_inb <- genes_inb[!duplicated(genes_inb)]

outdf <- genedf2
outfile <- "result_amp_genes.conf_99-INB.txt"
write.table(outdf, outfile, sep="\t", row.names=FALSE, quote=FALSE)
outfile <- "result_amp_genes.conf_99-INB.csv"
write.csv(outdf, outfile, row.names=FALSE, quote=FALSE)


sum(duplicated(genes_inb))
sum(duplicated(genes_ib))

library(venn)
x <- list("IB"=genes_ib, "INB"=genes_inb)
venn.result <- venn(
    x,
    ilabels=FALSE, zcolor="style", size=25, cexil=1.2, cexsn=1.5
)

outfile <- "fig_venn.pdf"
pdf(outfile, onefile=FALSE, width=10/2.54, height=10/2.54)
venn.result <- venn(
    x,
    ilabels=FALSE, zcolor="style", size=25, cexil=1.2, cexsn=1.5
)
dev.off()

genes_common <- intersect(genes_ib, genes_inb)
sum(duplicated(genes_common))
length(genes_common)
outdf <- data.frame(gene=genes_common)
outfile <- "gene_common.txt"
write.table(outdf, outfile, sep="\t", row.names=FALSE, quote=FALSE)
outfile <- "gene_common.csv"
write.csv(outdf, outfile, row.names=FALSE, quote=FALSE)

genes_only_ib <- genes_ib[!genes_ib %in% genes_common]
length(genes_only_ib)
outdf <- data.frame(gene=genes_only_ib)
outfile <- "genes_only_ib.txt"
write.table(outdf, outfile, sep="\t", row.names=FALSE, quote=FALSE)
outfile <- "genes_only_ib.csv"
write.csv(outdf, outfile, row.names=FALSE, quote=FALSE)

genes_only_inb <- genes_inb[!genes_inb %in% genes_common]
length(genes_only_inb)
outdf <- data.frame(gene=genes_only_inb)
outfile <- "genes_only_inb.txt"
write.table(outdf, outfile, sep="\t", row.names=FALSE, quote=FALSE)
outfile <- "genes_only_inb.csv"
write.csv(outdf, outfile, row.names=FALSE, quote=FALSE)









# Cluster analysis

library(ggplot2)
library(plyr)
library(stringr)
library(ape)
library(GOSemSim)
library(ggtree)
library(scales)
library(cowplot)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

fnames <- Sys.glob("enrichGO*.csv")
fnames

fdataset <- lapply(fnames, function(x){read.csv(x)[,c(2,3,5,10)]})
names(fdataset) <- fnames

ego.all <- ldply(fdataset, data.frame)
ego.all$group <- unlist(strsplit(ego.all$.id, split = ".csv"))
head(ego.all)

dim(ego.all)

ego.ID <- unique(ego.all[,c(2:4)])
head(ego.ID)
dim(ego.ID)

MyMerge <- function(x, y){
  df <- merge(x, y, by= "ID", all.x= TRUE, all.y= TRUE)
  return(df)
}
ego.m <- Reduce(MyMerge, fdataset)
head(ego.m)

ego.m <- ego.m[,c(1,4)]

ego.m <- merge(ego.ID[,1:2], ego.m, by= "ID", all.x= TRUE)
rownames(ego.m) <- ego.m$Description
ego.m$ID <- NULL
ego.m$Description <- NULL

colnames(ego.m) <- paste0("G", seq(1:length(fnames)))

head(ego.m)



hgGO <- godata('org.Hs.eg.db', ont="BP")
save(hgGO, file="hgGO.rdata")
(load("hgGO.rdata"))

ego.sim <- mgoSim(ego.ID$ID, ego.ID$ID, semData=hgGO, measure="Wang", combine=NULL)

ego.sim[1:3, 1:3]

rownames(ego.sim) <- ego.ID$Description
colnames(ego.sim) <- ego.ID$Description
ego.sim[1:3, 1:3]

tree <- nj(as.dist(1-ego.sim))
p <- ggtree(tree) + geom_tiplab() +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  coord_cartesian(xlim=c(-.1,1.3))
p







node <- c(17, 15, 14, 18)
gtree <- groupClade(tree, .node=node)
pbase <- ggtree(gtree, 
                aes(color=group))

fontsize <- 4
offset <- .9
pnode <- pbase + 

  geom_tiplab(size=4, align=TRUE) + 
  geom_cladelabel(node=node[1], align=TRUE, 

                  color = hue_pal()(length(node)+1)[2], 
                  fontsize = fontsize, offset=offset, label="pathway1") +

  geom_cladelabel(node=node[2], align=TRUE, color = hue_pal()(length(node)+1)[3], fontsize = fontsize, offset=offset, label="pathway2") +
  geom_cladelabel(node=node[3], align=TRUE, color = hue_pal()(length(node)+1)[4], fontsize = fontsize, offset=offset, label="pathway3") +
  geom_cladelabel(node=node[4], align=TRUE, color = hue_pal()(length(node)+1)[5], fontsize = fontsize, offset=offset, label="pathway4") +
  geom_cladelabel(node=node[4], align=TRUE, color = hue_pal()(length(node)+1)[6], fontsize = fontsize, offset=offset, label="pathway5") +
  
  coord_cartesian(xlim=c(-.1,1.5)) 

gheatmap(pnode, ego.m, 
         offset=.7,
         width=0.2,
         colnames_angle=90, hjust=0,
         low = rgb(255, 162, 115, max=255), high = rgb(255, 20, 9, max=255)
         ) 

ggsave("GOclustering.pdf", width = 5, height = 8)
