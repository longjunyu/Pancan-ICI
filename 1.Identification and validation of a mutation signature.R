# Identification and validation of a mutation signature

# Processing of mutation data
options(stringsAsFactors=FALSE)
library(ComplexHeatmap)
packageVersion("ComplexHeatmap")
library(RColorBrewer)
source("./oncoPrint_plus.R")
maf <- read.table("data_mutations_mskcc.txt", header=TRUE, sep="\t", quote="")
head(maf)
dim(maf)
table(maf[, "Variant_Classification"])
maf[, "gene"] <- maf[, "Hugo_Symbol"]
maf[, "mutation_type"] <- maf[, "Variant_Classification"]
maf[, "sample_id"] <- maf[, "Tumor_Sample_Barcode"]
maf[, "mutation_type2"] <- paste0(
    maf[, "Reference_Allele"], ">", maf[, "Tumor_Seq_Allele2"])
maf <- maf[, c("gene", "mutation_type", "sample_id", "mutation_type2")]

clin <- read.table("tmb_mskcc_2018_clinical_data.tsv", header=TRUE, sep="\t", quote="")
head(clin)
dim(clin)
clin[, "sample_id"] <- clin[, "Sample.ID"]

silent <- c("Silent", "Intron", "RNA", "3'UTR", "3'Flank",
            "5'UTR", "5'Flank", "IGR")
maf <- maf[!maf[, "mutation_type"] %in% silent,]
table(maf[, "mutation_type"])

keeps <- c(
"Missense_Mutation",
"Frame_Shift_Del",
"Frame_Shift_Ins",
"Nonsense_Mutation",
"Nonstop_Mutation",
"Splice_Site",
"Translation_Start_Site")
maf <- maf[maf[, "mutation_type"] %in% keeps,]
table(maf[, "mutation_type"])

maf[, "mutation_type2"] <- ifelse(
    nchar(maf[, "mutation_type2"]) != 3, "Others", maf[, "mutation_type2"])
maf[, "mutation_type2"] <- ifelse(
    grepl("-", maf[, "mutation_type2"]), "Others", maf[, "mutation_type2"])
maf[, "mutation_type2"] <-
ifelse(maf[, "mutation_type2"] %in% c("C>A", "G>T"), "C>A",
ifelse(maf[, "mutation_type2"] %in% c("C>G", "G>C"), "C>G",
ifelse(maf[, "mutation_type2"] %in% c("C>T", "G>A"), "C>T",
ifelse(maf[, "mutation_type2"] %in% c("T>A", "A>T"), "T>A",
ifelse(maf[, "mutation_type2"] %in% c("T>C", "A>G"), "T>C",
ifelse(maf[, "mutation_type2"] %in% c("T>G", "A>C"), "T>G", NA))))))
tb <- table(maf[, "mutation_type2"])
tb <- tb[order(-tb)]
tb

cntdf <- as.data.frame(table(maf[, "gene"]), stringsAsFactors=FALSE)
cntdf <- cntdf[order(-cntdf[, "Freq"]),]
genes <- cntdf[1:20, "Var1"]
genes
maf <- maf[maf[, "gene"] %in% genes,]
table(maf[, "mutation_type"])

sample_ids <- intersect(clin[, "sample_id"], maf[, "sample_id"])
sample_ids <- sample_ids[!duplicated(sample_ids)]
head(sample_ids)
clin <- clin[clin[, "sample_id"] %in% sample_ids,]
maf <- maf[maf[, "sample_id"] %in% sample_ids,]

x <- sort(table(factor(clin[,"Cancer.Type"])), dec=TRUE)
y <- numeric(0)
for (i in 1:length(x)) {
	if (i %% 2 == 1) {
		y <- c(x[i], y)
	} else {
		y <- c(y, x[i])
	}
}
y
length(y)
length(clin[,"Cancer.Type"])
clin[, "Cancer.Type"] <- factor(clin[, "Cancer.Type"], levels=names(y))

clin <- clin[order(clin[, "Cancer.Type"]),]
head(clin)
sample_ids <- clin[, "sample_id"]
length(sample_ids)

for (i in 1:length(genes)) {
    maf_2 <- maf[maf["gene"] == genes[i],]
    cntdf_2 <- as.data.frame(table(maf_2[, "gene"], maf_2[, "sample_id"]))
    cntdf_2 <- cntdf_2[, 2:3]
    cntdf_2[, 2] <- ifelse(cntdf_2[, 2] > 1, 1, cntdf_2[, 2])
    colnames(cntdf_2) <- c("sample_id", paste0("gene", i))
    clin <- merge(clin, cntdf_2, by="sample_id", all.x=TRUE)
    clin[, paste0("gene", i)] <- ifelse(
        is.na(clin[, paste0("gene", i)]), 0, clin[, paste0("gene", i)] * -1)
}
ordlst <- list()
ordlst[[1]] <- clin[, "Cancer.Type"]
for (i in 1:length(genes)) {
    ordlst[[i + 1]] <- clin[, paste0("gene", i)]
}
clin <- clin[do.call("order", ordlst),]

sample_ids <- clin[, "sample_id"]

mtx <- matrix(
    rep(0, length(sample_ids) * length(genes)),
    ncol=length(sample_ids),
    dimnames=list(genes, sample_ids))


mutation_types <- maf[!duplicated(maf[, "mutation_type"]), "mutation_type"]
mtxlst <- list()
for (i in 1:length(mutation_types)) {
    mtx2 <- mtx
    maf2 <- maf[maf[, "mutation_type"] == mutation_types[i],]
    cntdf2 <- as.data.frame(
        table(maf2[, "gene"], maf2[, "sample_id"]),
        stringsAsFactors=FALSE)
    for (j in 1:nrow(cntdf2)) {
        if (cntdf2[j, "Freq"] > 0) {
            mtx2[cntdf2[j, "Var1"], cntdf2[j, "Var2"]] <- 1
        }
    }
    mtxlst[[mutation_types[i]]] <- mtx2
}
length(mtxlst)



table(clin[,"Cancer.Type"])
cols_cancer <- c(
	rgb(130, 167, 95, max=255),  
	rgb(136, 255, 177, max=255),
	rgb(255, 203, 119, max=255), 
	rgb(255, 159, 106, max=255),
	rgb(167, 159, 255, max=255),
	rgb(168, 129, 217, max=255),
	rgb(228, 130, 198, max=255),
	rgb(236, 114, 99, max=255),
	rgb(255, 125, 140, max=255),
	rgb(251, 220, 142, max=255), 
	rgb(246, 255, 128, max=255)
)
cols_background <- unlist(lapply(
	as.numeric(clin[, "Cancer.Type"]), function(x) cols_cancer[x]))

cols <- c(
rgb(164, 142, 220, max=255),
rgb(202, 255, 114, max=255),
rgb(255, 110, 199, max=255),
rgb(255, 149, 114, max=255),
rgb(165, 159, 255, max=255),
rgb(255, 100, 83, max=255)
)
names(cols) <- mutation_types

alter_fun <- list(
    background=function(x, y, w, h, j, i) grid.rect(
    	x, y, w, h,
    	gp=gpar(fill=cols_background[j], col=NA, alpha=0.2)),
    Missense_Mutation=function(x, y, w, h) grid.rect(
        x, y, w*0.8, h*0.9,
        gp=gpar(fill=cols["Missense_Mutation"], col=NA)),
    Frame_Shift_Del=function(x, y, w, h) grid.rect(
        x, y, w*0.7, h*0.9,
        gp=gpar(fill=cols["Frame_Shift_Del"], col=NA)),
    Nonsense_Mutation=function(x, y, w, h) grid.rect(
        x, y, w*0.6, h*0.9,
        gp=gpar(fill=cols["Nonsense_Mutation"], col=NA)),
    Splice_Site=function(x, y, w, h) grid.rect(
        x, y, w*0.5, h*0.9,
        gp=gpar(fill=cols["Splice_Site"], col=NA)),
    Frame_Shift_Ins=function(x, y, w, h) grid.rect(
        x, y, w*0.4, h*0.9,
        gp=gpar(fill=cols["Frame_Shift_Ins"], col=NA)),
    Translation_Start_Site=function(x, y, w, h) grid.rect(
        x, y, w*0.3, h*0.9,
        gp=gpar(fill=cols["Translation_Start_Site"], col=NA))
    )
p <- oncoPrint(mtxlst, alter_fun=alter_fun, col=cols)

nrow(clin)
colnames(clin)[1:20]
clin[1:3,1:20]
table(clin[, "Overall.Survival.Status"])
clin[, "Age"] <- ifelse(is.na(clin[, "Age.at.Which.Sequencing.was.Reported..Days."]), NA,
    ifelse(clin[, "Age.at.Which.Sequencing.was.Reported..Days."] <= 65, "<=65", ">65"))
table(clin[, "Age"])
table(clin[, "Drug.Type"])
table(clin[, "Sex"])
max(clin[, "Overall.Survival..Months."])
max(clin[, "TMB.Score"])

cols_os <- unlist(lapply(
    as.numeric(clin[, "Cancer.Type"]), function(x) cols_cancer[x]))
cols_tmb <- unlist(lapply(
    as.numeric(clin[, "Cancer.Type"]), function(x) cols_cancer[x]))

cols_age <- c(
rgb(167, 159, 255, max=255),
rgb(249, 164, 28, max=255),
rgb(111, 144, 201, max=255),
rgb(109, 108, 217, max=255),
rgb(145, 200, 59, max=255)
)
cols_drug <- c(
rgb(181, 255, 249, max=255),
rgb(255, 227, 103, max=255),
rgb(230, 195, 250, max=255)
)
cols_status <- c(
rgb(133, 189, 108, max=255), 
rgb(224, 189, 240, max=255)
)
cols_sex <- c(
rgb(247, 153, 153, max=255),
rgb(164, 180, 220, max=255)
)
cols_mut <- c(
rgb(164, 142, 220, max=255),
rgb(202, 255, 114, max=255),
rgb(255, 110, 199, max=255),
rgb(255, 149, 114, max=255),
rgb(165, 159, 255, max=255),
rgb(255, 100, 83, max=255)
)

top_annotation_col <- list(
	OS_Status=c(
		"LIVING"=cols_status[1],
		"DECEASED"=cols_status[2]),
    Age_Group=c(
        "<=65"=cols_age[1],
        ">65"=cols_age[2]),
    Drug_Type=c(
        "Combo"=cols_drug[1],
        "CTLA4"=cols_drug[2],
        "PD-1/PDL-1"=cols_drug[3]),
    Sex=c(
        "Female"=cols_sex[1],
        "Male"=cols_sex[2]))
top_annotation <- HeatmapAnnotation(
    OS= anno_barplot(
        clin[, "Overall.Survival..Months."], size=unit(1, "mm"),
        border=FALSE, gp=gpar(fill=cols_os, col=cols_os), axis=TRUE, ylim=c(0, 80)),
    TMB_Score=anno_barplot(
        log2(clin[, "TMB.Score"]),
        border=FALSE, gp=gpar(col=cols_tmb, fill=cols_tmb),
        axis=TRUE, ylim=c(0, log2(250))),
    OS_Status=clin[, "Overall.Survival.Status"],
    Age_Group=clin[, "Age"],
    Drug_Type=clin[, "Drug.Type"],
    Sex=clin[, "Sex"],
    annotation_height=unit.c(rep(unit(1.5, "cm"), 2), rep(unit(0.5, "cm"), 4)),
    annotation_legend_param=list(
        labels_gp=gpar(fontsize=9),
        title_gp=gpar(fontsize=10, fontface="bold"), ncol=1),

    gap=unit(c(1, 1, 1, 1, 1), "mm"),
    col=top_annotation_col,
    show_annotation_name=TRUE,
    annotation_name_gp=gpar(fontsize=12))
names(cols_cancer) <- names(table(clin[,"Cancer.Type"]))
bottom_annotation_col <- list(
    Cancer_Type=c(cols_cancer))

maf_2 <- maf
maf_2[, "sample_id"] <- factor(maf_2[, "sample_id"], levels=sample_ids)
pctmtx <- as.matrix(table(maf_2[, "sample_id"], maf_2[, "mutation_type2"]))
pctmtx <- t(apply(pctmtx, 1, function(x) x/sum(x)))
for (i in 1:ncol(pctmtx)) {
    pctmtx[, i] <- ifelse(is.na(pctmtx[, i]), 0, pctmtx[, i])
}
bottom_annotation <- HeatmapAnnotation(
    Cancer_Type=clin[, "Cancer.Type"],
    Mutation_Type=anno_barplot(
    	pctmtx, gp=gpar(col=cols_mut, fill=cols_mut),
    	border=FALSE, axis=TRUE, ylim=c(0, 1), extend=0),

    annotation_height=unit.c(unit(0.2, "cm"),unit(1, "cm")),
    annotation_legend_param=list(
        labels_gp=gpar(fontsize=9),
        title_gp=gpar(fontsize=10, fontface="bold"), ncol=1),

    gap=unit(1, "mm"),
    col=bottom_annotation_col,
    show_annotation_name=TRUE,
    annotation_name_gp=gpar(fontsize=12),
    name="bottom")

lgd_lst <- list(
    Legend(
        labels=c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
        title="Mutation_Type",
        legend_gp=gpar(fill=cols_mut))
)

p2 <- oncoPrint_plus(
    mtxlst, alter_fun=alter_fun, col=cols, alter_fun_is_vectorized=FALSE,
    column_order=sample_ids,
    top_annotation=top_annotation,
    bottom_annotation=bottom_annotation,
    column_split=clin[, "Cancer.Type"],
    column_title=NULL)


pdf("oncoprint.pdf", width=25, height=10)
draw(p2, heatmap_legend_side="bottom", annotation_legend_side="bottom",
     merge_legend=TRUE,
	 heatmap_legend_list=lgd_lst,

	 gap=unit(1.5, "mm"))
dev.off()

tiff("oncoprint.tiff", width=25, height=10, unit="in", res=600)
draw(p2, heatmap_legend_side="bottom", annotation_legend_side="bottom",
	 merge_legend=TRUE,
	 heatmap_legend_list=lgd_lst,

	 gap=unit(1.5, "mm"))
dev.off()





# Identification and validation of a mutation signature
library(survival)
rt=read.table("training.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,"futime"]=rt[,"futime"]/12
rt[1:6,1:6]
dim(rt)
rt <- subset(rt, futime > 0)
dim(rt)

library(survival)
library(glmnet)
set.seed(1000)
x <- as.matrix(rt[,3:ncol(rt)])
y <- Surv(rt$futime, rt$fustat)

model <- glmnet(x, y, family="cox", alpha=1)
plot(model, xvar = "lambda", label = TRUE)
cv.model <- cv.glmnet(x, y, family="cox", alpha=1, nfolds=10)
plot(cv.model)
plot(cv.model$glmnet.fit, xvar="lambda", label=TRUE)
cv.model$lambda.min
cv.model$lambda.1se

library(boot)
lambda_boot <- function(data, indices) {
    d <- data[indices,]
    x <- as.matrix(rt[,3:ncol(d)])
    y <- Surv(d$futime, d$fustat)
    cv.model <- cv.glmnet(x, y, family="cox", alpha=1, nfolds=10)
    return(cv.model$lambda.1se)
    
}

library(boot)
set.seed(1000)
bootcorr <- boot(data=rt, statistic=lambda_boot, R=1000)
bootcorr
summary(bootcorr)
plot(bootcorr)
plot(bootcorr[["t"]])
lambda.1se <- bootcorr[["t"]]
x<-coef(cv.model, s = c(lambda.1se[,1]))
ad <- as.matrix(x)
write.csv(ad, "boot_result-1000.csv")

lambda.1se <- mean(bootcorr[["t"]])
lambda.1se
c(cv.model$lambda.1se,lambda.1se)

Coefficients <- coef(cv.model, s=cv.model$lambda.min)
Coefficients

Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Columns <- rownames(Coefficients)[Active.Index]
cbind(Active.Columns, Active.Coefficients)

predict(cv.model, s=cv.model$lambda.1se, type="coefficients")
predict(cv.model, s=lambda.1se, type="coefficients")

lambda.1se
Coefficients <- coef(cv.model, s=lambda.1se)
Coefficients
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Columns <-rownames(Coefficients)[Active.Index]
cbind(rownames(Coefficients)[Active.Index], Active.Coefficients)
rt[1:6,1:6]


cox <- coxph(Surv(futime, fustat) ~ BRAF+PAK7+PTPRD+ROS1+SETD2+VHL+FAM46C+PTPRT+TET1+ZFHX3, data=rt)
cox
summary(cox)
sum.surv<-summary(cox)
c_index<-sum.surv$concordance
c_index

library(Hmisc)
fp <- predict(cox)
cindex=1-rcorr.cens(fp,Surv(rt$futime,rt$fustat)) [[1]]
cindex

library(survcomp)
cindex <- concordance.index(predict(cox),surv.time = rt$futime, surv.event = rt$fustat,method = "noether")
cindex$c.index; 
cindex$lower; 
cindex$upper; 
cindex$p.value;

riskScore <- predict(cox,type="risk",newdata=rt)
median(riskScore)

risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))
table(risk)
write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="lasso-cox-risk-median(riskScore)1.txt",sep="\t",quote=F,row.names=F)


riskScore <- predict(cox,type="risk",newdata=rt)
risk <- as.vector(ifelse(riskScore > 1.05,"high","low"))
write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="lasso-cox-risk-X-tile1.txt",sep="\t",quote=F,row.names=F)
write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="risk1.txt",sep="\t",quote=F,row.names=F)
cox
summary(cox)
exp(coef(cox))
exp(confint(cox))


getwd() 
rt2=read.table("validation.txt",header=T,sep="\t",check.names=F,row.names=1)
rt2[, "futime"] <- rt2[,"futime"] / 12
riskScore <- predict(cox,type="risk",newdata=rt)
median(riskScore)
riskScore <- predict(cox, type="risk", newdata=rt2)
median(riskScore)
risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt2[,1:2],riskScore,risk)),cbind(rt2[,1:2],riskScore,risk)),file="lasso-cox-risk-median(riskScore)2.txt",sep="\t",quote=F,row.names=F)

riskScore <- predict(cox, type="risk", newdata=rt2)
risk <- as.vector(ifelse(riskScore > 1.05, "high", "low"))
write.table(
    cbind(id=rownames(cbind(rt2[,1:2], riskScore, risk)),
          cbind(rt2[, 1:2], riskScore, risk)),
    file="risk2.txt", sep="\t", quote=F, row.names=F)

cox
summary(cox)
exp(coef(cox))
exp(confint(cox))






# Kaplan-Meier curves
library(survival)
library(survminer)
library(ggfortify)
library(ggplot2)
library(grid)
library(ggsci)
library(scales)
rt <- read.table("risk1.txt", header=T, sep="\t")
dim(rt)
rt[1:5, 1:5]
rt[, "futime"] <- rt[, "futime"] * 365.25 / 30
max(rt[, "futime"])
rt[, "risk"] <- factor(rt[, "risk"], levels=c("low", "high"))

digit <- 3
digit_format <- paste0("%.", digit, "f")

diff <- survdiff(Surv(futime, fustat) ~ risk, data=rt)
diff
print(diff)
pvalue <- 1 - pchisq(diff[["chisq"]], length(diff[["n"]]) - 1)
pvalue <- ifelse(
    pvalue < 0.001, "P < 0.001", 
    paste0("P = ", sprintf(digit_format, round(pvalue, digit))))

cox <- coxph(Surv(futime, fustat) ~ risk, data=rt, ties="breslow")
hr <- summary(cox)[["coefficients"]][, "Pr(>|z|)"]
summary(cox)

print(summary(cox)[["conf.int"]])
hr <- summary(cox)[["conf.int"]][, 1]
hr_lcl <- summary(cox)[["conf.int"]][, 3]
hr_ucl <- summary(cox)[["conf.int"]][, 4]
hr_lcl_ucl <- paste0(
    "Hazard Ratio = ",
    sprintf(digit_format, round(hr, digit)), "\n95% CI: ",
    sprintf(digit_format, round(hr_lcl, digit)), " - ",
    sprintf(digit_format, round(hr_ucl, digit)))

fit <- survfit(Surv(futime, fustat) ~ risk, data=rt)
fit
summary(fit)
summary(fit)[["table"]]

plot(fit, mark.time=TRUE,
    lty=2:3, col=c("red", "blue"),
    xlab="time (year)", ylab="Survival rate",
    main=paste0("Survival curve (", pvalue, ")"))
legend("topright",
    c("low risk", "high risk"), lty=2:3, col=c("red", "blue"))

autoplot(fit, censor.shape='*')

tb <- summary(fit)$table
labels <- paste0(c("Low ", "High"), " ", round(tb[, 7], 2), " mo")
ggsurv <- ggsurvplot(
    fit,
    legend=c(0.8, 0.85),
    legend.title=paste0(c("Risk", rep(" ", 10), "Median"), collapse=""),
    legend.labs=labels,
    palette=pal_jco("default")(2),
    size=c(2, 2),
    pval=paste0(pvalue, "\n", hr_lcl_ucl),
    xlab="Time (months)",
    ylab="Overall Survival (probability)",
    xlim=c(0, 84),
    ylim=c(0, 1),
    break.time.by=24,
    censor=FALSE,
    risk.table=TRUE,
    risk.table.title="No. at risk",
    risk.table.height=0.15,
    tables.theme=theme(
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        axis.title=element_blank(),
        panel.grid.major.y=element_blank()),
    ggtheme=theme(
        axis.text=element_text(size=12),
        axis.line=element_line(color="black"),
        legend.background=element_rect(color="black"),
        legend.key=element_rect(fill=NA),
        legend.key.width=unit(30, "points"),
        panel.background=element_rect(fill=NA),
        panel.grid=element_blank()))
ggsurv$plot <- ggsurv$plot +
    scale_x_continuous(
        expand=c(0, 0), limit=c(0, 84),breaks=seq(24, 84, 24)) +
    scale_y_continuous(
        expand=c(0, 0), limit=c(0, 1.02), breaks=seq(0.2, 1, 0.2)) +
    theme(
        plot.margin=unit(c(10, 20, 10, 30), "points"),
        axis.title.y=element_text(margin=margin(r=-60))) +
    coord_cartesian(clip="off") +
    annotation_custom(textGrob("0", x=-0.03, y=-0.02, hjust=0))
ggsurv$table <- ggsurv$table +
    scale_x_continuous(expand=c(0, 0), limit=c(0, 82)) +
    scale_y_discrete(labels=c("High", "Low")) +
    coord_cartesian(clip="off") + 
    theme(
        plot.title=element_text(hjust=-0.17),
        axis.text.y=element_text(margin=margin(r=40), hjust=0),
        plot.margin=unit(c(10, 20, 10, 30), "points"))
ggsurv
pdf("fig_km1.pdf", onefile=FALSE, width=20/2.54, height=18/2.54)
ggsurv
dev.off()
tiff("fig_km1.tiff", width=20, height=18, unit="cm", res=350)
ggsurv
dev.off()




# receiver operating characteristic (ROC) curves
options(stringsAsFactors=FALSE)
library(timeROC)
library(survival)
tr <- read.table("risk1.txt", header=T, sep="\t")
dim(tr)
tr[1:5,1:5]
tr$time <- tr$futime
tr$status <- tr$fustat

times <- c(3, 6)
roc <- timeROC(
    T=tr$time, delta=tr$status, marker=tr$riskScore,
    cause=1,
    weighting="marginal",
    times=times,
    iid=TRUE)
roc
text1 <- paste0(
    "3 Years: AUC(95%CI),\n",
    round(roc[["AUC"]][1], 3), "(",
    paste0(confint(roc)[["CI_AUC"]][1,], collapse="-"), ")")
text2 <- paste0(
    "5 Years: AUC(95%CI),\n",
    round(roc[["AUC"]][2], 3), "(",
    paste0(confint(roc)[["CI_AUC"]][2,], collapse="-"), ")")

plot(roc, time=3, col="orange", title=FALSE)
plot(roc, time=6, col="blue", add=TRUE, title=FALSE)
text(0.4, 0.6, label=text1, adj=0)
text(0, 0.9, label=text2, adj=0)

library(ggplot2)
df <- data.frame(
    TP=c(roc$TP[, 1], roc$TP[, 2]),
    FP=c(roc$FP[, 1], roc$FP[, 2]),
    group=c(rep("3", nrow(roc$TP)), rep("5", nrow(roc$TP))))
p <- ggplot(data=df, aes(FP, TP, group=group, color=group)) + 
    geom_line() +
    annotate("text", x=0.4, y=0.6, label=text1, color="black", hjust=0) +
    annotate("text", x=0, y=0.9, label=text2, color="black", hjust=0) +
    geom_hline(yintercept=seq(0.2, 1, 0.2), color="gray") +
    scale_color_manual(values=c(
        rgb(42, 83, 95, max=255), rgb(247, 147, 28, max=255))) +
    scale_x_continuous("1 - Specificity", breaks=seq(0, 1, 0.2)) +
    scale_y_continuous("Sensitivity", breaks=seq(0, 1, 0.2), expand=c(0, 0)) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
    theme_bw() +
    theme(
        axis.line=element_line(color="black"),
        legend.position="none",
        panel.grid=element_blank(),
        panel.border=element_blank())
p
pdf("fig_roc.pdf", width=15/2.54, height=15/2.54)
p
dev.off()
tiff("fig_roc.tiff", width=15, height=15, unit="cm", res=350)
p
dev.off()





#  Forest plot
options(stringsAsFactors=FALSE)
library(forestplot)
rawdf <- read.csv("result_hrdf.csv", header=TRUE)
newdf <- rawdf
label <- c(newdf[, "Var1"])
label[!c(1:length(label)) %in% c(1, 2, 3, 5, 8, 10, 12:14, 16, 19, 21)] <- ""
label1 <- label
label <- c(newdf[, "Var2"])
label2 <- label
n1 <- c(newdf[, "n1"])
n2 <- c(newdf[, "n2"])

pct <- c(
    (n1[1:11] + n2[1:11]) / (n1[2] + n2[2]),
    (n1[12:22] + n2[12:22]) / (n1[13] + n2[13]))
pct <- rep(0.12, length(pct))

pct[c(2, 13)] <- 1
hr <- c(newdf[, "hr"])
hr_lcl <- c(newdf[, "hr_lcl"])
hr_ucl <- c(newdf[, "hr_ucl"])
hr_ucl <- ifelse(hr_ucl == "Inf", NA, hr_ucl)
digit <- 2
digit_format <- paste0("%.", digit, "f")
hr_lcl_ucl <- paste0(
    sprintf(digit_format, round(hr, digit)), "(",
    sprintf(digit_format, round(hr_lcl, digit)), "-",
    sprintf(digit_format, round(hr_ucl, digit)), ")")
hr_lcl_ucl <- ifelse(hr_lcl_ucl == "NA(NA-NA)", "", hr_lcl_ucl)
pvalue <- c(newdf[, "pvalue"])
pvalue <- ifelse(pvalue < 0.001, "<0.001", format(round(pvalue, 3), digits=3))
tabletext <- cbind(
    c("Group", "\n", label1),
    c("", "\n", label2),
    c("Risk high", "\n", n1),
    c("Risk low", "\n", n2),
    c("Hazard Ratio(95%CI)", "\n", hr_lcl_ucl),
    c("P Value", "\n", pvalue))
cols <- c(
    rgb(247, 148, 29, max=255),
    rgb(46, 86, 98, max=255),
    rgb(0, 173, 238, max=255))

jpeg("fig_forestplot.jpg", width=30, height=30, unit="cm", res=350)
forestplot(
    labeltext=tabletext, graph.pos=5,
    mean=c(NA, NA, hr),
    lower=c(NA, NA, hr_lcl),
    upper=c(NA, NA, hr_ucl),
    title="",
    xlab=paste0(
         c(rep(" ", 30), "<--Favors Alive--    --Favors Deceased-->"), collapse=""),
    clip=c(0.02, 3),
    xticks=c(0.02, 0.4, 1, 3), xticks.digits=1,
    xlog=T,
    boxsize=c(NA, NA, pct),
    is.summary=c(
        rep(TRUE, 4), rep(FALSE, 9),
        rep(TRUE, 2), rep(FALSE, 9)),
    align=c("l", "l", "r", "r", "r", "r"),
    txt_gp=fpTxtGp(
        label=gpar(cex=1.2),
        ticks=gpar(cex=1.1),
        xlab=gpar(cex=1.2),
        title=gpar(cex=1.2)),
    col=fpColors(
        box=cols[2],
        lines=cols[2],
        zero="gray50",
        summary=cols[1]),
    fn.ci_norm="fpDrawCircleCI",
    zero=1,
    cex=0.9,
    lineheight="auto",
    colgap=unit(6, "mm"),
    lwd.ci=2,
    ci.vertices=FALSE,
    ci.vertices.height=0.4,
    hrzl_lines=list(
        "3"=gpar(lty=1, col="black", cex=3),
        "14"=gpar(lty=1, col="black"),
        "25"=gpar(lty=1, col="black")))
dev.off()

pdf("fig_forestplot.pdf", onefile=FALSE, width=30/2.54, height=30/2.45)
forestplot(
    labeltext=tabletext, graph.pos=5,
    mean=c(NA, NA, hr),
    lower=c(NA, NA, hr_lcl),
    upper=c(NA, NA, hr_ucl),
    title="",
    xlab=paste0(
        c(rep(" ", 30), "<--Favors Alive--    --Favors Deceased-->"), collapse=""),
    clip=c(0.02, 3),
    xticks=c(0.02, 0.3, 1, 3), xticks.digits=1,
    xlog=T,
    boxsize=c(NA, NA, pct),
    is.summary=c(
        rep(TRUE, 4), rep(FALSE, 9),
        rep(TRUE, 2), rep(FALSE, 9)),
    align=c("l", "l", "r", "r", "r", "r"),
    txt_gp=fpTxtGp(
        label=gpar(cex=1.2),
        ticks=gpar(cex=1.1),
        xlab=gpar(cex=1.2),
        title=gpar(cex=1.2)),
    col=fpColors(
        box=cols[2],
        lines=cols[2],
        zero="gray50",
        summary=cols[1]),
    fn.ci_norm="fpDrawCircleCI",
    zero=1,
    cex=0.9,
    lineheight="auto",
    colgap=unit(6, "mm"),
    lwd.ci=2,
    ci.vertices=FALSE,
    ci.vertices.height=0.4,
    hrzl_lines=list(
        "3"=gpar(lty=1, col="black", cex=3),
        "14"=gpar(lty=1, col="black"),
        "25"=gpar(lty=1, col="black")))
dev.off()






# The evaluation of the mutation signature

# The comparison of the C-index
library(survival)
library(survcomp)
library(ggsci)
bc <- read.csv("Nomogram.csv")
dim(bc)
bc[1:6, 1:6]
bc <- bc[complete.cases(bc),]
dim(bc)
bc[1:6, 1:6]
fit1 <- coxph(Surv(futime, fustat) ~ CancerType, data=bc)
ci1 <- concordance.index(
    predict(fit1), surv.time=bc[["futime"]], surv.event=bc[["fustat"]],
    method="noether")
ci1$c.index; 

fit2 <- coxph(Surv(futime, fustat) ~ DrugType, data=bc)
ci2 <- concordance.index(
    predict(fit2), surv.time=bc[["futime"]], surv.event=bc[["fustat"]],
    method="noether")

fit3 <- coxph(Surv(futime, fustat) ~ TMB, data=bc)
ci3 <- concordance.index(
    predict(fit3), surv.time=bc[["futime"]], surv.event=bc[["fustat"]],
    method="noether")
ci3$c.index; 
ci3$lower; 
ci3$upper

fit4 <- coxph(Surv(futime, fustat) ~ Risk, data=bc)
ci4 <- concordance.index(
    predict(fit4), surv.time=bc[["futime"]], surv.event=bc[["fustat"]],
    method="noether")
ci4$c.index; 
ci4$lower; 
ci4$upper

cindex.comp(ci4, ci1)
cindex.comp(ci4, ci2)
cindex.comp(ci4, ci3)
cindex.comp(ci4, ci5)


# Bar plot
cin1 <- ci4
cin2 <- ci1
var1 <- "Risk"
var2 <- "TMB"

comp <- cindex.comp(cin1, cin2)
comp
if (comp$cindex1 <= comp$cindex2) print("C-index排序错误")
df <- data.frame(

    type=c(var1, var2),
    value=c(cin1$c.index, cin2$c.index),
    lower=c(cin1$lower, cin2$lower),
    upper=c(cin1$upper, cin2$upper),
    pvalue=c(cin1$p.value, cin2$p.value))
df[, "pvalue"] <- ifelse(is.na(df[, "pvalue"]), "",
    ifelse(df[, "pvalue"] < 0.001, "P < 0.001",
        paste0("P = ", df[, "pvalue"])))
library(ggplot2)
cols <- pal_jco("default")(4)[c(2, 1)]
fig <- ggplot(data=df, aes(x=type, y=value, fill=type)) +

    geom_col(width=0.5) +

    geom_errorbar(aes(ymin=lower, ymax=upper, width=0.25)) +

    geom_text(aes(y=0.9, label=pvalue)) +

    scale_fill_manual(values=cols) +
    scale_x_discrete("") +

    scale_y_continuous("", limit=c(0, 1)) +
    theme_bw() +
    theme(legend.position="top", legend.title=element_blank(),
          panel.grid=element_blank())
fig
outfile <- paste0("fig_bar_cindex_", var1, "_", var2, ".tiff")
tiff(outfile, width=30, height=30, unit="cm", res=350)
fig
dev.off()
outfile <- paste0("fig_bar_cindex_", var1, "_", var2, ".pdf")
pdf(outfile, onefile=FALSE, width=30/2.54, height=30/2.54)
fig
dev.off()




# The comparisons of the IDI and NRI
options(stringsAsFactors=FALSE)
library(survival)
dat <- read.csv("Nomogram.csv", header=TRUE)
dim(dat)
dat[1:6,1:7]
dat$futime <- as.numeric(dat$futime)
max(dat$futime)
min(dat$futime)

t0 <- 5

indata0 <- as.matrix(dat[, c("futime", "fustat", "TMB")])
indata0[1:6,1:3]

covs0 <- as.matrix(indata0[,c(-1,-2)])
head(covs0)


indata1 <- as.matrix(dat[, c("futime", "fustat", "Risk")])
indata1[1:6,1:3]

covs1 <- as.matrix(indata1[,c(-1,-2)])
head(covs1)

library(survIDINRI)
head(dat[,2:3])
set.seed(123456)

x <- IDI.INF(dat[, c("futime", "fustat")], covs0, covs1, t0, npert=1)

IDI.INF.OUT(x)
x$m1


IDI.INF.GRAPH(x)


shade <- function(x, yl, yu, ...){ 
    nl <- length(x)
    for(i in 1:(nl-1))
        {polygon(c(x[i], x[i], x[i+1], x[i+1]), c(yl[i], yu[i], yu[i+1], yl[i+1]), ...)}
}
IDI.INF.GRAPH2 <- function(
    x, main=NULL, xlab = NULL, ylab = NULL, cex.main=NULL, cex.lab=NULL,...){

    library(ggsci)

    cols <- c(col=rgb(150,85,255, max=255),col=rgb(255,34,177, max=255))
    library(scales)
    show_col(pal_jco("default")(10))

    cc.diff=x$point$cc
    diffs1=1-x$point$FX
    diffs0=1-x$point$GX
   
    Xlim<-c(min(cc.diff[diffs1>0.001 | diffs0>0.001]),max(cc.diff[diffs1<0.999 | diffs0<0.999])) ; Ylim<-c(0, 1)
    plot(cc.diff, diffs1, type="l", lty=1, lwd=3, xlim=Xlim, ylim=Ylim, xlab="", ylab="",...)
    diffs0_a=diffs0; diffs0_a[diffs0>diffs1]=diffs1[diffs0>diffs1]
    diffs1_a=diffs1; diffs1_a[diffs1>diffs0]=diffs0[diffs1>diffs0]

    shade(cc.diff, diffs1, diffs0_a,col=cols[2],border=NA) 
    shade(cc.diff, diffs0, diffs1_a,col=cols[1],border=NA) 

    lines(cc.diff, diffs1, lty=1, lwd=3, col=cols[6]) 
    lines(cc.diff, diffs0, lty=1, lwd=1, col=cols[2]) 
    abline(v=0.0, col="gray", lty=1)
    abline(h=0.5, col="gray", lty=1)

    c0.idx<-which(abs(cc.diff)==min(abs(cc.diff)))

    points(0,diffs1[c0.idx],pch=19, col=cols[3], cex=3) ;
    points(0,diffs0[c0.idx],pch=19, col=cols[3], cex=3) ; 
    m1.idx<-which(abs(diffs1-0.5)==min(abs(diffs1-0.5)))
    m0.idx<-which(abs(diffs0-0.5)==min(abs(diffs0-0.5)))

    points(cc.diff[m1.idx[1]], 0.5, col=cols[4], pch=19, cex=3) ;
    points(cc.diff[m0.idx[1]], 0.5, col=cols[4], pch=19, cex=3) ;
    if(is.null(xlab)){xlab="s"}
    if(is.null(ylab)){ylab=expression(paste("pr(",hat(D)<=s,")"))}
    title(main=main, xlab=xlab, ylab=ylab,cex.lab=cex.lab, cex.main=cex.main)
}
IDI.INF.GRAPH2(x)

tiff("fig_idi_tmb.tiff", width=30, height=30, unit="cm", res=350)
IDI.INF.GRAPH2(x)
dev.off()
pdf("fig_idi_tmb.pdf", onefile=FALSE, width=30/2.54, height=30/2.54)
IDI.INF.GRAPH2(x)
dev.off()



table(dat[, "CancerType"])
covs <- model.matrix(~ CancerType, dat[, c("futime", "fustat", "CancerType")])
ncol(covs)
head(covs)
indata0 <- cbind(as.matrix(dat[, c("futime", "fustat")]), covs)
indata0[1:6,1:4]

covs0<-as.matrix(indata0[,c(-1,-2,-3)])
head(covs0)

indata1 <- as.matrix(dat[, c("futime", "fustat", "Risk")])
indata1[1:6, 1:3]

covs1 <- as.matrix(indata1[,c(-1,-2)])
head(covs1)

library(survIDINRI)
head(dat[,2:3])
set.seed(1234)

x <- IDI.INF(dat[, c("futime", "fustat")], covs0, covs1, t0, npert=1)

IDI.INF.OUT(x)
x$m1

IDI.INF.GRAPH(x)

tiff("fig_idi_cancertype.tiff", width=30, height=30, unit="cm", res=350)
IDI.INF.GRAPH2(x)
dev.off()
pdf("fig_idi_cancertype.pdf", onefile=FALSE, width=30/2.54, height=30/2.54)
IDI.INF.GRAPH2(x)
dev.off()



table(dat[, "DrugType"])
covs <- model.matrix(~ DrugType, dat[, c("futime", "fustat", "DrugType")])
ncol(covs)
head(covs)
indata0 <- cbind(as.matrix(dat[, c("futime", "fustat")]), covs)
indata0[1:6,1:4]

covs0<-as.matrix(indata0[,c(-1,-2,-3)])
head(covs0)

indata1 <- as.matrix(dat[, c("futime", "fustat", "Risk")])
indata1[1:6, 1:3]

covs1 <- as.matrix(indata1[,c(-1,-2)])
head(covs1)

library(survIDINRI)
head(dat[,2:3])
set.seed(123456)

x <- IDI.INF(dat[, c("futime", "fustat")], covs0, covs1, t0, npert=1)

IDI.INF.OUT(x)
x$m1

IDI.INF.GRAPH(x)

tiff("fig_idi_drugtype.tiff", width=30, height=30, unit="cm", res=350)
IDI.INF.GRAPH2(x)
dev.off()
pdf("fig_idi_drugtype.pdf", onefile=FALSE, width=30/2.54, height=30/2.54)
IDI.INF.GRAPH2(x)
dev.off()









# The construction and evaluation of the nomogram

# Nomogram
library(survival)
library(rms)
library(ggsci)


bc <- read.csv("Nomogram1.csv")
dim(bc)
bc[1:6, 1:6]


levels(factor(bc[, "CancerType"]))
bc[, "CancerType"] <- factor(as.numeric(factor(bc[, "CancerType"])))

dd <- datadist(bc)
options(datadist="dd")


cox <- cph(
    Surv(futime, fustat) ~ Risk + CancerType + DrugType + TMB, data=bc,
    x=TRUE, y=TRUE, surv=TRUE)

surv <- Survival(cox)
nom <- nomogram(
    cox, fun=list(
        function(x) surv(1, x), 
        function(x) surv(3, x),
        function(x) surv(5, x)), 
    lp=FALSE,
    funlabel=c("1-year survival", "3-year survival", "5-year survival"),
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01))
plot(nom)


pdf("fig_nom_all.pdf", onefile=FALSE, width=150/2.54, height=30/2.54)
plot(nom)
dev.off()




# Calibration plot
library(survival)
library(rms)
library(ggsci)
bc <- read.csv("Nomogram1.csv")
dim(bc)
bc[1:6, 1:6]
dd<- datadist(bc)
options(datadist="dd")


time <- 1
cox <- cph(
    Surv(futime, fustat) ~ Risk + CancerType + DrugType + TMB, data=bc,
    x=TRUE, y=TRUE, surv=TRUE, time.inc=time)

set.seed(111)
cal <- calibrate(cox, cmethod="KM", method="boot", u=time, m=626, B=1000)

cols <- pal_jco("default")(10)
outfile <- paste0("fig_cal_all", time, ".pdf")
pdf(outfile, onefile=FALSE, width=20/2.54, height=20/2.54)
plot(cal,
    lwd=2, lty=1, errbar.col=cols[2],
    xlim=c(0.445, 0.8), ylim=c(0.38, 0.8),
    xlab="Nomogram-prediced OS (%)", ylab="Observed OS (%)",
    col=cols[1],
    cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal[, c("mean.predicted", "KM")],
      type="b", lwd=1, col=cols[1], pch=16)
mtext("")
dev.off()



time <- 3
cox <- cph(
    Surv(futime, fustat) ~ Risk + CancerType + DrugType + TMB, data=bc,
    x=TRUE, y=TRUE, surv=TRUE, time.inc=time)

set.seed(111)
cal <- calibrate(cox, cmethod="KM", method="boot", u=time, m=626, B=1000)

cols <- pal_jco("default")(10)
outfile <- paste0("fig_cal_all", time, ".pdf")
pdf(outfile, onefile=FALSE, width=20/2.54, height=20/2.54)
plot(cal,
    lwd=2, lty=1, errbar.col=cols[2],
    xlim=c(0.15, 0.54), ylim=c(0.1, 0.6),
    xlab="Nomogram-prediced OS (%)", ylab="Observed OS (%)",
    col=cols[1],
    cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal[, c("mean.predicted", "KM")],
      type="b", lwd=1, col=cols[1], pch=16)
mtext("")
dev.off()



time <- 5
cox <- cph(
    Surv(futime, fustat) ~ Risk + CancerType + DrugType + TMB, data=bc,
    x=TRUE, y=TRUE, surv=TRUE, time.inc=time)

set.seed(111)
cal <- calibrate(cox, cmethod="KM", method="boot", u=time, m=626, B=1000)

cols <- pal_jco("default")(10)
outfile <- paste0("fig_cal_all", time, ".pdf")
pdf(outfile, onefile=FALSE, width=20/2.54, height=20/2.54)
plot(cal,
    lwd=2, lty=1, errbar.col=cols[2],
    xlim=c(0.07, 0.41), ylim=c(0.08, 0.47),
    xlab="Nomogram-prediced OS (%)", ylab="Observed OS (%)",
    col=cols[1],
    cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal[, c("mean.predicted", "KM")],
      type="b", lwd=1, col=cols[1], pch=16)
mtext("")
dev.off()