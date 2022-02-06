# 1.1 Propensity score matching (PSM) weighting algorithm

options(stringsAsFactors=FALSE)
clindf2 <- read.table("Samstein et al. (PMID 30643254)_OS_pancancer_clinic_all_gene.txt", header=TRUE, sep="\t", quote="")
dim(clindf2)
colnames(clindf2)[1:50]
clindf2[1:6, 1:30]

clindf2[, "status"] <- clindf2[, "Os_fustat"]
clindf2[, "time"] <- as.numeric(clindf2[, "Os_futime"])

rownames(clindf2) <- clindf2[, "Sample_id"]

table(clindf2[, "Cancer_type"])
clindf2[, "Cancer_type"] <- gsub(" ", ".", clindf2[, "Cancer_type"])
clindf2[, "Cancer_type"] <- gsub("\\-", ".", clindf2[, "Cancer_type"])
dim(clindf2)
table(clindf2[, "Cancer_type"])

table(clindf2[, "Drug_type"])
clindf2[, "Drug_type"] <- gsub("/", ".", clindf2[, "Drug_type"])
clindf2[, "Drug_type"] <- gsub("\\-", ".", clindf2[, "Drug_type"])
clindf2[, "Drug_type"] <- gsub("\\+", ".", clindf2[, "Drug_type"])
dim(clindf2)
table(clindf2[, "Drug_type"])

colnames(clindf2)[1:30]
clindf3 <- clindf2[complete.cases(clindf2[, c("time", "status")]),]
dim(clindf3)
clindf3[1:6,1:6]
table(clindf3$Cancer_type)

table(clindf3$Cancer_type)
tb <- table(clindf3[, "Cancer_type"])
cancers <- names(tb[tb > 1])
clindf3 <- clindf3[clindf3[, "Cancer_type"] %in% cancers,]
dim(clindf3)
table(clindf3$Cancer_type)

dim(clindf3)
clindf3 <- clindf3[!clindf3[, "Cancer_type"] %in% "Cancer.of.unknown.primary",]
dim(clindf3)

colnames(clindf3)[1:30]
table(clindf3$Age)
table(clindf3$Cancer_type)
table(clindf3$Gender)
table(clindf3$TMB)
table(clindf3$Drug_type)

vars <- c(
    "Age",
    "Drug_type",
    "Cancer_type"
)
dim(clindf3)

clindf33 <- clindf3

sum(!complete.cases((clindf3[, vars])))
clindf3[!complete.cases(clindf3[, vars]),]

clindf3 <- clindf3[complete.cases(clindf3[, vars]),]
dim(clindf3)
clindf3[1:6,1:6]
table(clindf3$Age)
table(clindf3$Gender)
table(clindf3$TMB)
table(clindf3$Drug_type)
table(clindf3$Cancer_type)
sum(table(clindf3$Cancer_type))

colnames(clindf3)[1:30]
colnames(clindf3)[(ncol(clindf3)-6):ncol(clindf3)]
dim(clindf3)
clindf3[1:6,1:6]
ncol(clindf3)
clindf3[1:6,(ncol(clindf3)-6):ncol(clindf3)]

sums <- apply(clindf3[23:(ncol(clindf3)-2)], 2, sum)

genes <- names(sums[sums > 2])
length(genes)
head(genes)
tail(genes)
clindf4 <- clindf3[, c(vars, genes)]
dim(clindf4)

clindf4[1:6,1:6]
table(clindf4$Cancer_type)





colnames(clindf3)[1:50]
clean <- clindf3[, c(colnames(clindf3)[1:22], genes)]
dim(clean)
clean[1:6,1:6]
write.table(cbind(id=rownames(clean), clean), "clean_clinic_gene.txt",sep="\t", row.names=F, quote=FALSE)
write.csv(cbind(id=rownames(clean), clean), "clean_clinic_gene.csv",row.names=F, quote=FALSE)


clean_survival_gene <- clindf3[, c("time", "status", genes)]
dim(clean_survival_gene)
clean_survival_gene[1:6,1:6]
write.table(cbind(id=rownames(clean_survival_gene), clean_survival_gene), "clean_survival_gene.txt",sep="\t", row.names=F, quote=FALSE)
write.csv(cbind(id=rownames(clean_survival_gene), clean_survival_gene), "clean_survival_gene.csv",row.names=F, quote=FALSE)


multiInput_raw <- clindf33[, c("time", "status", genes)]
dim(multiInput_raw)
multiInput_raw[1:6,1:6]
colnames(multiInput_raw)[1:2] <- c("futime", "fustat")
dim(multiInput_raw)
multiInput_raw[1:6,1:6]
write.table(cbind(id=rownames(multiInput_raw), multiInput_raw), "multiInput_raw.txt",sep="\t", row.names=F, quote=FALSE)
write.csv(cbind(id=rownames(multiInput_raw), multiInput_raw), "multiInput_raw.csv",row.names=F, quote=FALSE)



library("dummies")
library("tidyverse")
source("cal.R")

folder <- "PSM_Output"
if (!file.exists(folder)) {
    dir.create(folder)
}

cancerNames <- "all"

dim(clindf3)
clindf3[1:6,1:6]
colnames(clindf3)[1:30]
ImmFeature <- clindf3[, c("time", "status")]
dim(ImmFeature)
ImmFeature[1:6,1:2]


dim(clindf4)
clindf4[1:6,1:6]
clindf4[1:6, c(vars, genes[1])]


analysis <- "All"
sum.ImmFeatureAll <- data.frame()
for (cancer in cancerNames) {

    print(paste0("-------", cancer, "--------"))


dim(clindf4)

fealst <- list()
for (h in 1:length(genes)) {
    print(paste0("-------", genes[h], "--------"))
    data <- clindf4[, c(vars, genes[h])]
    colnames(data)[1] <- "age"

    analysis <- "gene"
    if (analysis == "gene") {

        colnames(data)[which(colnames(data) == genes[h])] <- "Z"
    }

    continuous.feature <- c("Z", "age")
    dummy.feature <- setdiff(colnames(data), continuous.feature)
    if (length(dummy.feature) > 0) {
        data.dum <- dummy.data.frame(data, names=dummy.feature)
        head(data.dum)
        dummy.list <- attr(data.dum, "dummies")
        rm.col <- c()
        for (i in 1:length(dummy.list)) {
            rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
        }
        data.dum <- data.dum[,-rm.col]
        data.dum$X0 <- rep(1, nrow(data.dum))

        exclude.col <- match(c("Z", "X0"), colnames(data.dum))
        colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
        form <- as.formula(
            paste0("Z~", paste0(colnames(data.dum)[-exclude.col], collapse=" + "))
        )
    } else {
        data.dum <- data
        data.dum$X0 <- rep(1, nrow(data.dum))
        exclude.col <- match(c("Z", "X0"), colnames(data.dum))
        colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
        form <- as.formula(
            paste0("Z~", paste0(colnames(data.dum)[-exclude.col], collapse=" + "))
        )
    }

    Feature <- ImmFeature
    Feature <- cbind(barcode=rownames(Feature), Feature)
    Feature.pri <- Feature[2:ncol(Feature)]
    colnames(Feature.pri) <- colnames(Feature)[2:ncol(Feature)]
    rownames(Feature.pri) <- Feature$barcode
    Feature.pri <- rm.zero.col(Feature.pri)

    head(Feature.pri)
    if (!file.exists(folder)) {
        dir.create(folder)
    }
    folder2 <- paste0(folder, "/", cancer, "_", analysis)
    if (!file.exists(folder2)) {
        dir.create(folder2)
    }
    Feature.result <- weight.test(
        data.dum, form, Feature.pri, is.continuous=TRUE, is.survival=TRUE,
        weight=ifelse(analysis == "gene", "MW", "ATT"), mirror.plot=FALSE,
        paste0(folder2, "/", cancer), data.type="Feature",
        outdir=paste0(folder, "/", cancer, "_", analysis), perm=FALSE
    )
    fealst[[h]] <- as.data.frame(Feature.result)
}
length(fealst)
Feature.result <- do.call(rbind, fealst)
length(genes)
nrow(Feature.result)
Feature.result[, "feature"] <- genes
Feature.result[, "fdr"] <- p.adjust(Feature.result[, "pvalue"], "fdr")
head(Feature.result)
Feature.result <- as.list(Feature.result)

    sum.Immune <- summarize.p(Feature.pri, Feature.result, print=TRUE)
    summarize.p(Feature.pri, Feature.result, print=TRUE, cutoff=0.05)
    write.summary(sum.Immune, paste0(folder2, "/", cancer), analysis, "Immune")
    write.result(Feature.result, paste0(folder2, "/", cancer), analysis, "Immune")
    save(Feature.result, file=paste0(folder2, "/", cancer, "_", analysis, "_result.RData"))
    if (length(which(Feature.result$pvalue < 0.05)) > 0) {
        sum.Immune <- data.frame(sum.Immune)
        sum.Immune$class <- rep(cancer, times=nrow(sum.Immune))
        if (nrow(sum.ImmFeatureAll) == 0) {
            sum.ImmFeatureAll <- sum.Immune
        } else {
            sum.ImmFeatureAll <- rbind(sum.ImmFeatureAll, sum.Immune)
        }
    }
}

write.table(sum.ImmFeatureAll, "feature difference.across.Cancer_typesAll.txt",sep="\t", row.names=FALSE, quote=FALSE)

















# 1.2 Lasso-penalized Cox regression analysis

library(survival)
rt=read.table("multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)
dim(rt)
rt[1:6,1:6]
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














# 1.3 Multivariate Cox regression analysis

library(survival)
rt=read.table("multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)
dim(rt)
rt[1:6,1:6]
rt[,"futime"]=rt[,"futime"]/12
rt[1:6,1:6]

dim(rt)
rt[1:6,1:6]
cox <- coxph(Surv(futime, fustat) ~ BRAF+PAK7+PTPRD+PTPRT+ROS1+SETD2+TET1+VHL+FAM46C+RNF43+ZFHX3, data=rt)
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

risk <- as.vector(ifelse(riskScore > 1.07,"high","low"))
table(risk)

write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="lasso-cox-risk-X-tile1.txt",sep="\t",quote=F,row.names=F)
write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="risk1.txt",sep="\t",quote=F,row.names=F)
cox
summary(cox)
exp(coef(cox))
exp(confint(cox))






getwd() 
rt2=read.table("merge_data4.txt",header=T,sep="\t",check.names=F,row.names=1)
dim(rt2)
rt[1:5,1:5]
rt2[1:5,1:5]

rt2[, "futime"] <- rt2[,"futime"] / 12
riskScore <- predict(cox,type="risk",newdata=rt)
median(riskScore)
riskScore <- predict(cox, type="risk", newdata=rt2)
median(riskScore)
risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt2[,1:2],riskScore,risk)),cbind(rt2[,1:2],riskScore,risk)),file="lasso-cox-risk-median(riskScore)2.txt",sep="\t",quote=F,row.names=F)

riskScore <- predict(cox, type="risk", newdata=rt2)
risk <- as.vector(ifelse(riskScore > 1.07, "high", "low"))

head(risk)
table(risk)
write.table(
    cbind(id=rownames(cbind(rt2[,1:2], riskScore, risk)),
          cbind(rt2[, 1:2], riskScore, risk)),
    file="risk2.txt", sep="\t", quote=F, row.names=F)

cox
summary(cox)
exp(coef(cox))
exp(confint(cox))






















# 1.4 Survival analysis

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
max(rt[, "futime"])
rt[, "futime"] <- rt[, "futime"] * 12
rt[1:5, 1:5]
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

    palette=c("#AC88FF", "#FC717F"),

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






rt <- read.table("risk2.txt", header=T, sep="\t")
dim(rt)
rt[1:5, 1:5]
max(rt[, "futime"])

rt[, "futime"] <- rt[, "futime"] * 12
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


cox <- coxph(Surv(futime, fustat) ~ risk, data=rt, ties="breslow") #ָ��ties="breslow"
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


tb <- summary(fit)$table
labels <- paste0(c("Low ", "High"), " ", round(tb[, 7], 2), " mo")
ggsurv <- ggsurvplot(
    fit,
    legend=c(0.8, 0.85),
    legend.title=paste0(c("Risk", rep(" ", 10), "Median"), collapse=""),
    legend.labs=labels,
    palette=c("#AC88FF", "#FC717F"),
    size=c(2, 2),
    pval=paste0(pvalue, "\n", hr_lcl_ucl),
    xlab="Time (months)",
    ylab="Overall Survival (probability)",
    xlim=c(0, 48),
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
        expand=c(0, 0), limit=c(0, 48), breaks=seq(24, 48, 24)) +
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
pdf("fig_km2.pdf", onefile=FALSE, width=20/2.54, height=18/2.54)
ggsurv
dev.off()
tiff("fig_km2.tiff", width=20, height=18, unit="cm", res=350)
ggsurv
dev.off()















# 1.5 Receiver operating characteristic (ROC) curve analysis

options(stringsAsFactors=FALSE)
tr <- read.table("risk1.txt", header=T, sep="\t")
dim(tr)
tr[1:5,1:5]
tr$time <- tr$futime
tr$status <- tr$fustat

times <- c(3, 5)
roc <- timeROC(
    T=tr$time, 
    delta=tr$status, 
    marker=tr$riskScore,
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
plot(roc, time=5, col="blue", add=TRUE, title=FALSE)
text(0.4, 0.6, label=text1, adj=0)
text(0, 0.9, label=text2, adj=0)



smooth.roc.binormal <- function(roc, order) {
se1 <- roc[["FP"]][, order]
sp1 <- roc[["TP"]][, order]
if (length(unique(se1)) <= 3) {

ses <- unique(se1)
ses
sps <- NULL
for (i in 1:length(ses)) {
    sps2 <- sp1[se1 == ses[i]]
    sps <- c(sps, sps2[length(sps2)])
}
sps
selst <- list()
splst <- list()
for (i in 1:(length(ses)-1)) {
    ses2 <- seq(ses[i], ses[i+1], length=10)
    sps2 <- (
        (sps[i] - sps[i+1])/(ses[i] - ses[i+1]) * (ses2 - ses[i+1])
        + sps[i+1]
    )
    selst[[i]] <- ses2
    splst[[i]] <- sps2
}
se1 <- c(0, unlist(selst))
sp1 <- c(0, unlist(splst))
}
n <- 512
df <- data.frame(sp=qnorm(sp1), se=qnorm(se1))
df <- df[apply(df, 1, function(x) all(is.finite(x))),]
if (dim(df)[1] <= 1)
    stop("ROC curve not smoothable (not enough points).")
model <- lm(sp ~ se, df)
if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
se <- qnorm(seq(0, 1, 1/(n-1)))
sp <- predict(model, data.frame(se))
se2 <- pnorm(se)
sp2 <- pnorm(sp)
se_t2 <- se2
sp_t2 <- sp2
return(list(se=se_t2, sp=sp_t2, auc=roc[["AUC"]][order]))
}


outfile <- paste0("fig_roc_smooth_train.pdf")
pdf(outfile, onefile=FALSE, width=15/2.54, height=15/2.54)
par(mar=c(4.5, 4.5, 2, 2))
roc1 <- smooth.roc.binormal(roc, 1)
roc2 <- smooth.roc.binormal(roc, 2)
cols <- c("red", "blue", "darkgreen", "purple", "black")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(3, 3, 3, 3, 3)
legends <- c(
    paste0("3 Years: AUC = ",sprintf("%. 3f", roc1[["auc"]])),
    paste0("5 Years: AUC = ",sprintf("%. 3f", roc2[["auc"]]))
)
plot(
    roc1[["se"]], roc1[["sp"]], type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(roc1[["se"]], 1), c(roc1[["sp"]], 0), col="lightblue")
lines(roc2[["se"]], roc2[["sp"]], col=cols[2], lty=ltys[2], lwd=lwds[2])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA, bg=NA
)
dev.off()


outfile <- paste0("fig_roc_smooth_train.tiff")
tiff(outfile, width=15, height=15, unit="cm", res=600, compression="lzw+p")
par(mar=c(4.5, 4.5, 2, 2))
roc1 <- smooth.roc.binormal(roc, 1)
roc2 <- smooth.roc.binormal(roc, 2)
cols <- c("red", "blue", "darkgreen", "purple", "black")
ltys <- c(1, 1, 1, 1, 1)
lwds <- c(3, 3, 3, 3, 3)
legends <- c(
    paste0("3 Years: AUC = ",sprintf("%. 3f", roc1[["auc"]])),
    paste0("5 Years: AUC = ",sprintf("%. 3f", roc2[["auc"]]))
)
plot(
    roc1[["se"]], roc1[["sp"]], type="l", col=cols[1], lty=ltys[1], lwd=lwds[1],
    xlab="1 - Specificity",
    ylab="Sensitivity"
)
polygon(c(roc1[["se"]], 1), c(roc1[["sp"]], 0), col="lightblue")
lines(roc2[["se"]], roc2[["sp"]], col=cols[2], lty=ltys[2], lwd=lwds[2])
abline(0, 1, col="grey")
legend(
    0.9, 0.4, legend=legends, text.col=cols, adj=1, box.col=NA, bg=NA
)
dev.off()

library(ggplot2)
df <- data.frame(
    TP=c(roc$TP[, 1], roc$TP[, 2]),
    FP=c(roc$FP[, 1], roc$FP[, 2]),
    group=c(rep("4", nrow(roc$TP)), rep("6", nrow(roc$TP))))
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

























# 1.6 Comparison of C-indexes among the mutation-based gene set, frameshift insertion/deletion (indel) mutation burden, tobacco mutation signature, UV signature, APOBEC signature, and DNA damage response pathway mutation. 

options(stringsAsFactors=FALSE)
infile <- "2. frameshift insertiondeletion (indel) mutation burden/data_mutations_mskcc.txt"

mafdf <- read.table(infile, header=TRUE, sep="\t", quote="")
dim(mafdf)
mafdf[1:6,]

mafdf[1:6, "Tumor_Sample_Barcode"]
mafdf[, "sample_id"] <- mafdf[, "Tumor_Sample_Barcode"]
mafdf[1:6, "sample_id"]
mafdf[, "gene"] <- mafdf[, "Hugo_Symbol"]
mafdf[1:6, "gene"]
mafdf[, "mutation_type"] <- mafdf[, "Variant_Classification"]
mafdf[1:6, "mutation_type"]


table(mafdf[, "mutation_type"])

silent <- c(
"3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
"RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA",
"De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
"Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del"
)

nonsyns <- c(
"Missense_Mutation",
"Nonsense_Mutation",
"Nonstop_Mutation",
"Splice_Site",
"Translation_Start_Site",
"Frame_Shift_Del",
"Frame_Shift_Ins",
"In_Frame_Del",
"In_Frame_Ins"
)


mafdf[, "mutation_type2"] <- paste0(
    mafdf[, "Reference_Allele"], ">", mafdf[, "Tumor_Seq_Allele2"]
)
mafdf[, "mutation_type2"] <- ifelse(
    nchar(mafdf[, "mutation_type2"]) != 3, "Others", mafdf[, "mutation_type2"]
)
mafdf[, "mutation_type2"] <- ifelse(
    grepl("-", mafdf[, "mutation_type2"]), "Others", mafdf[, "mutation_type2"]
)
mafdf[, "mutation_type2"] <- (
    ifelse(mafdf[, "mutation_type2"] %in% c("C>A", "G>T"), "C>A",
    ifelse(mafdf[, "mutation_type2"] %in% c("C>G", "G>C"), "C>G",
    ifelse(mafdf[, "mutation_type2"] %in% c("C>T", "G>A"), "C>T",
    ifelse(mafdf[, "mutation_type2"] %in% c("T>A", "A>T"), "T>A",
    ifelse(mafdf[, "mutation_type2"] %in% c("T>C", "A>G"), "T>C",
    ifelse(mafdf[, "mutation_type2"] %in% c("T>G", "A>C"), "T>G", NA))))))
)
tb <- table(mafdf[, "mutation_type2"])
tb <- tb[order(-tb)]
tb

sample_ids <- unique(mafdf[, "sample_id"])
length(sample_ids)

infile <- "2. frameshift insertiondeletion (indel) mutation burden/newdf-ѵ����.csv"
clindf <- read.csv(infile, header=TRUE, check.names=FALSE)
dim(clindf)
clindf[1:6,]
colnames(clindf)
clindf[, "fustat"] <- clindf[, "Os_fustat"]
clindf[, "futime"] <- clindf[, "Os_futime"]

table(clindf[, "Cohort"])
cohort <- "Samstein et al. (PMID 30643254)"
clindf <- clindf[clindf[, "Cohort"] == cohort,]

mafdf2 <- mafdf[mafdf[, "mutation_type"] %in% nonsyns,]
table(mafdf2[, "mutation_type"])
length(table(mafdf2[, "sample_id"]))
freqdf <- as.data.frame(table(mafdf2[, "sample_id"]))
freqdf[1:6,]
colnames(freqdf) <- c("Sample_id", "n_nonsyns")
freqdf[1:6,]
max(freqdf[, "n_nonsyns"])
min(freqdf[, "n_nonsyns"])

clindf <- merge(clindf, freqdf, by="Sample_id", all.x=TRUE)
dim(clindf)
clindf[1:6,]

clindf[, "n_nonsyns"] <- ifelse(
    is.na(clindf[, "n_nonsyns"]), 0, clindf[, "n_nonsyns"]
)
clindf[1:20, c("Sample_id", "TMB", "n_nonsyns")]

table(mafdf[, "Variant_Type"])
mafdf2 <- mafdf[mafdf[, "Variant_Type"] %in% c("DEL", "INS"),]
table(mafdf2[, "Variant_Type"])
length(table(mafdf2[, "sample_id"]))
freqdf <- as.data.frame(table(mafdf2[, "sample_id"]))
freqdf[1:6,]
colnames(freqdf) <- c("Sample_id", "n_indel")
freqdf[1:6,]
max(freqdf[, "n_indel"])
min(freqdf[, "n_indel"])

clindf <- merge(clindf, freqdf, by="Sample_id", all.x=TRUE)
dim(clindf)
clindf[1:6,]

clindf[, "n_indel"] <- ifelse(
    is.na(clindf[, "n_indel"]), 0, clindf[, "n_indel"]
)


table(mafdf[, "Variant_Type"])
mafdf2 <- mafdf[mafdf[, "Variant_Type"] %in% c("SNP"),]
table(mafdf2[, "Variant_Type"])
length(table(mafdf2[, "sample_id"]))
freqdf <- as.data.frame(table(mafdf2[, "sample_id"]))
freqdf[1:6,]
colnames(freqdf) <- c("Sample_id", "n_snv")
freqdf[1:6,]
max(freqdf[, "n_snv"])
min(freqdf[, "n_snv"])

clindf <- merge(clindf, freqdf, by="Sample_id", all.x=TRUE)
clindf[1:6,]
dim(clindf)

clindf[, "n_snv"] <- ifelse(is.na(clindf[, "n_snv"]), 0, clindf[, "n_snv"])



clindf[, "prop_indel"] <- (
    clindf[, "n_indel"] / (clindf[, "n_indel"] + clindf[, "n_snv"])
)

clindf[, "prop_indel"] <- ifelse(
    clindf[, "n_indel"] + clindf[, "n_snv"] == 0, 0, clindf[, "prop_indel"]
)

genelst <- list(
"DNA damage repair"=c(
"MSH2",
"MSH6",
"PMS2",
"POLE",
"BRCA2"
)
)
genes <- unlist(genelst)
genes[!genes %in% mafdf[, "gene"]]

genelst2 <- list(
"DNA damage repair"=c(
"MSH2",
"MSH6",
"PMS2",
"POLE",
"BRCA2"
)
)
for (i in 1:length(genelst2)) {
pathway <- names(genelst2)[[i]]
genes <- genelst2[[i]]
mafdf2 <- mafdf[mafdf[, "gene"] %in% genes,]
dim(mafdf2)
table(mafdf2[, "mutation_type"])
length(table(mafdf2[, "sample_id"]))
freqdf <- as.data.frame(table(mafdf2[, "sample_id"]))
freqdf[1:6,]
colnames(freqdf) <- c("Sample_id", pathway)
freqdf[1:6,]
max(freqdf[, pathway])
min(freqdf[, pathway])


freqdf[, pathway] <- ifelse(freqdf[, pathway] > 1, 1, freqdf[, pathway])


clindf <- merge(clindf, freqdf, by="Sample_id", all.x=TRUE)
dim(clindf)
clindf[1:6,]


clindf[, pathway] <- ifelse(is.na(clindf[, pathway]), 0, clindf[, pathway])
clindf[, pathway] <- factor(clindf[, pathway])
}
clindf[1:6, c("Sample_id", "TMB", pathway)]



infile <- "4.Tobacco mutation signature/data_mutations_mskcc.txt"

mafdf <- read.table(infile, header=TRUE, sep="\t", quote="")
dim(mafdf)
mafdf[1:6,]
mafdf[, "Chromosome"] <- paste0("chr", mafdf[, "Chromosome"])



library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg19)

mafdf[1:6, 1:12]
mut.ref <- mafdf
sigs.input <- mut.to.sigs.input(
    mut.ref=mut.ref,
    sample.id="Tumor_Sample_Barcode",
    chr="Chromosome",
    pos="Start_Position",
    ref="Reference_Allele",
    alt="Tumor_Seq_Allele2",
    bsg=BSgenome.Hsapiens.UCSC.hg19
)
sigs.input[1:6, 1:6]

sample_ids <- rownames(sigs.input)
length(sample_ids)
samplst <- list()
for (i in 1:length(sample_ids)) {
    sigs.sample <- whichSignatures(
	    tumor.ref=sigs.input,
        signatures.ref=signatures.cosmic,
        sample.id=sample_ids[i],
        contexts.needed=TRUE,
        tri.counts.method="default"
    )
    sigs.sample[["weights"]]
    samplst[[i]] <- data.frame(sigs.sample[["weights"]])
}
sampdf <- do.call(rbind, samplst)
sampdf[1:6,]
sampdf[, "Sample_id"] <- rownames(sampdf)

sampdf[, "Smoking_sig"] <- ifelse(sampdf[, "Signature.4"] > 0.06, "1", "0")
sampdf[, "UV_sig"] <- ifelse(sampdf[, "Signature.7"] > 0.06, "1", "0")

sampdf[, "APOBEC_sig"] <- ifelse(sampdf[, "Signature.2"] > 0.06, "1", "0")



clindf <- merge(clindf, sampdf, by="Sample_id", all.x=TRUE)
dim(clindf)
clindf[1:6,]

clindf[, "Smoking_sig"] <- ifelse(is.na(clindf[, "Smoking_sig"]), "0", clindf[, "Smoking_sig"])
clindf[, "UV_sig"] <- ifelse(is.na(clindf[, "UV_sig"]), "0", clindf[, "UV_sig"])
clindf[, "APOBEC_sig"] <- ifelse(is.na(clindf[, "APOBEC_sig"]), "0", clindf[, "APOBEC_sig"])




library(survival)
library(survcomp)
vars <- c(
"risk",
"prop_indel",
"Smoking_sig",
"UV_sig",
"APOBEC_sig",
"DNA damage repair"
)
vars

cilst <- list()
for (i in 1:length(vars)) {
var <- vars[i]
formula <- paste0("Surv(futime, fustat) ~ `", var, "`")
model <- coxph(as.formula(formula), data=clindf)
ci <- concordance.index(
    predict(model),
    surv.time=clindf[["futime"]], surv.event=clindf[["fustat"]],
    method="noether"
)
ci[["c.index"]]
ci[["lower"]]
ci[["upper"]]
cilst[[i]] <- ci
}

cilst2 <- list()
for (i in 1:length(vars)) {
    var <- vars[i]
    ci <- cilst[[i]]
    cidf <- data.frame(
        var=var,
        type=var,
        value=ci$c.index,
        lower=ci$lower,
        upper=ci$upper
    )
    cilst2[[i]] <- cidf
}
cidf <- do.call(rbind, cilst2)
cidf

cmp <- cindex.comp(cilst[[1]], cilst[[2]])
cmp[["p.value"]]


library(ggplot2)
library(ggsci)
newdf1 <- cidf
newdf1 <- newdf1[order(newdf1[, "value"]),]
levels <- newdf1[, "type"]
newdf1[, "type"] <- factor(newdf1[, "type"], levels=levels)
newdf1[, "value"]
ymax <- max(newdf1[, "upper"])

colors <- rainbow(length(levels), s=1, v=0.8, alpha=1)
set.seed(100)
colors <- sample(colors)
plot <- ggplot(data=newdf1, aes(x=type, y=value, fill=type)) +
    geom_hline(yintercept=0.5, linetype=2, color="grey") +
    geom_col(width=0.5) +
    geom_text(aes(y=value*1.1, label=sprintf("%.3f", value))) +

    scale_fill_manual(values=colors) +
    scale_x_discrete("") +
    scale_y_continuous("", limit=c(0, 0.8), breaks=seq(0, 0.8, 0.2)) +
    theme_bw() +
    theme(
        axis.text.x=element_text(angle=0, hjust=0.5),
        legend.position="none",
        legend.title=element_blank(),
        panel.grid=element_blank()
    )
print(plot)

outfile <- "fig_bar_cindex.tiff"
tiff(outfile, width=12, height=10, unit="cm", res=350)
print(plot)
dev.off()

outfile <- "fig_bar_cindex.pdf"
pdf(outfile, onefile=FALSE, width=12/2.54, height=10/2.54)
print(plot)
dev.off()












# 1.7 Comparison of C-indexes among the mutation-based gene set, B2M mutation, JAK1 mutation, JAK2 mutation, KRAS mutation, TP53 mutation, PTEN mutation, STK11 mutation, and BAP1 mutation

options(stringsAsFactors=FALSE)
bc <- read.csv("newdf.csv")

dim(bc)
bc[1:6,]
colnames(bc)

infile <- "Samstein et al. (PMID 30643254)_OS_pancancer_clinic_all_gene.txt"
clindf <- read.table(infile, header=TRUE, sep="\t", check.names=FALSE)
dim(clindf)
clindf[1:6, 1:60]
colnames(clindf)

bc <- merge(bc, clindf[, c(1, 23:ncol(clindf))], by="Sample_id", all.x=TRUE)
dim(bc)
bc[1:6, 1:6]


library(survival)
library(survcomp)
library(ggsci)

fit4 <- coxph(Surv(futime, fustat) ~ risk, data=bc)
ci4 <- concordance.index(
    predict(fit4), surv.time=bc[["futime"]], surv.event=bc[["fustat"]],
    method="noether"
)
ci4$c.index
ci4$lower
ci4$upper


genes <- c(
"KIR3DS1",
"B2M",
"JAK1", 
"JAK2",
"KRAS", 
"TP53",
"PTEN",
"RTK",
"STK11",
"BAP1"
)
for (i in 1:length(genes)) {
    print(genes[i])
    print(colnames(bc)[grepl(genes[i], colnames(bc))])
}

genes <- c(
"B2M",
"JAK1", 
"JAK2",
"KRAS", 
"TP53",
"PTEN",
"STK11",
"BAP1"
)
cilst <- list()
cmplst <- list()
for (i in 1:length(genes)) {
    bc[, genes[i]] <- factor(bc[, genes[i]])
    formula <- paste0("Surv(futime, fustat) ~ ", genes[i])
    fit <- coxph(as.formula(formula), data=bc)
    ci <- concordance.index(
        predict(fit), surv.time=bc[["futime"]], surv.event=bc[["fustat"]],
        method="noether"
    )
    cmp <- cindex.comp(ci4, ci)
    print(genes[i])
    print(ci$c.index)
    print(cmp$p.value)
    cilst[[i]] <- ci
    cmplst[[i]] <- cmp
}


cmplst2 <- list()
for (i in 1:length(genes)) {
    gene <- genes[i]
    ci <- cilst[[i]]
    cmp <- cmplst[[i]]
    cmpdf <- data.frame(
        cindex0=ci4$c.index, cindex_lcl0=ci4$lower, cindex_ucl0=ci4$upper,
        gene=gene,
        cindex=ci$c.index, cindex_lcl=ci$lower, cindex_ucl=ci$upper,
        pvalue=cmp$p.value
    )
    cmplst2[[i]] <- cmpdf
}
cmpdf <- do.call(rbind, cmplst2)
cmpdf


figlst <- list()
for (i in 1:length(genes)) {
    gene <- genes[i]
    ci <- cilst[[i]]
    cmp <- cmplst[[i]]
    figdf <- data.frame(
        gene=gene,
        type=gene,
        value=ci$c.index,
        lower=ci$lower,
        upper=ci$upper,
        pvalue=cmp$p.value
    )
    figlst[[i]] <- figdf
}
figlst[[length(genes)+1]] <- data.frame(
    gene="Best",
    type="Best",
    value=ci4$c.index,
    lower=ci4$lower,
    upper=ci4$upper,
    pvalue=NA
)
figdf <- do.call(rbind, figlst)
figdf <- figdf[order(figdf[, "value"]),]

library(ggplot2)
newdf1 <- figdf
levels <- newdf1[, "type"]
newdf1[, "type"] <- factor(newdf1[, "type"], levels=levels)
newdf1[, "value"]
ymax <- max(newdf1[, "upper"])

colors <- rainbow(length(levels), s=1, v=0.8, alpha=1)
set.seed(100)
colors <- sample(colors)
plot <- ggplot(data=newdf1, aes(x=type, y=value, fill=type)) +
    geom_hline(yintercept=0.5, linetype=2, color="grey") +
    geom_col(width=0.5) +
    geom_text(aes(y=value*1.1, label=sprintf("%.3f", value))) +

    scale_fill_manual(values=colors) +
    scale_x_discrete("") +
    scale_y_continuous("", limit=c(0, 0.8), breaks=seq(0, 0.8, 0.2)) +
    theme_bw() +
    theme(
        legend.position="none",
        legend.title=element_blank(),
        panel.grid=element_blank()
    )
print(plot)

outfile <- "fig_bar_cindex_train.tiff"
tiff(outfile, width=12, height=10, unit="cm", res=350)
print(plot)
dev.off()

outfile <- "fig_bar_cindex_train.pdf"
pdf(outfile, onefile=FALSE, width=12/2.54, height=10/2.54)
print(plot)
dev.off()
