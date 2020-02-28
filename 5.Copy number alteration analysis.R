# Copy number alteration analysis

# Copy number profiles
options(stringsAsFactors=FALSE)
scores <- read.table(
    "GISTIC_2_0_23/cnv_results1/scores.gistic", head=TRUE, sep="\t", quote="")
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





# Gene ontology analysis
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
pvalueFilter=0.05
adjpFilter=1
showNum=30

rt=read.table("id.txt",sep="\t",header=T,check.names=F)
dim(rt)
head(rt)

rt=rt[is.na(rt[,"entrezID"])==F,]
dim(rt)
head(rt)

gene=rt$entrezID
length(gene)
head(gene)

colorSel="p.adjust"
if(adjpFilter>0.05){
	colorSel="pvalue"
}

for(i in c("BP","CC","MF")){

	kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont=i, readable =T)
	GO=as.data.frame(kk)
    GO=GO[(GO$pvalue<pvalueFilter & GO$p.adjust<adjpFilter),]
	write.table(GO,file=paste0(i,".txt"),sep="\t",quote=F,row.names = F)                 #保存富集结果
	if(nrow(GO)<showNum){
		showNum=nrow(GO)
	}
	if(nrow(GO)!=0){

		pdf(file=paste0(i,".barplot.pdf"),width = 11,height = 7)
		bar=barplot(kk, drop = TRUE, showCategory =showNum,color = colorSel)
		print(bar)
		dev.off()

		pdf(file=paste0(i,".bubble.pdf"),width = 11,height = 7)
		bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio", color = colorSel)
		print(bub)
		dev.off()
	}
}






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




node <- c(33, 36, 39, 40, 38)
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
