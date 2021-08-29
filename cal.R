library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
hist.mirror2 <- function(data.plot, label) {
    myplot <- ggplot(data.plot, aes(ps, fill=status))+
        geom_histogram(
            data=subset(data.plot, status=="Yes"), aes(ps, y=..count..)
        )+
        geom_histogram(
            data=subset(data.plot,status=="No"), aes(ps, y= -..count..)
        )+
        scale_fill_hue(label)+
        labs(x="Propensity Score")
    print(myplot)
}


weight.test <- function(
    data,form, molecular.pri, is.continuous=TRUE, weight="MW",
    mirror.plot=FALSE, cancer, data.type, outdir=".",
    perm=FALSE, seed=seed, is.survival=FALSE
) {

    rownames(molecular.pri) <- gsub("\\.","\\-",rownames(molecular.pri))
    common <- intersect(rownames(data), rownames(molecular.pri))
    print(paste("Number of samples:", length(common)))

    molecular.common <- as.matrix(
        molecular.pri[match(common, rownames(molecular.pri)),]
    )
    colnames(molecular.common) <- colnames(molecular.pri)
    clinical.common <- data[match(common, rownames(data)),]
    print(summary(factor(clinical.common$Z)))

    write.table(
        t(rbind(rownames(clinical.common), ifelse(clinical.common$Z==1, "Mutation", "Wild"))),
        file=paste0(cancer, "_", data.type,"_samples.txt"),
        sep="\t", row.names=FALSE, quote=FALSE
    )

    print("----------------")
    if("age" %in% colnames(clinical.common)) {
        age.cutoff <- 65
        print(paste(
            "age <", age.cutoff,":",
            length(which(clinical.common$age < age.cutoff))
        ))
        print(paste(
            "age >=", age.cutoff,":",
            length(which(clinical.common$age >= age.cutoff))
        ))
        print("----------------")
    }

    if (perm) {
        n <- nrow(clinical.common)
        set.seed(seed)
        perm <- sample(1:n, n)
        clinical.common$Z <- clinical.common$Z[perm]
    }

    print(paste0("Weighting scheme: ", weight))


    source("GIPW_function_omega.R")

    ans <- GIPW.std.omega(
        dat=clinical.common, form.ps=form, weight=weight, trt.est=FALSE
    )


    if (mirror.plot) {
        ps <- ans$ps
        status <- ifelse(clinical.common$Z == 1, "Yes", "No")
        data.plot <- data.frame(ps, status)
        pdf(paste0(cancer,"_", data.type, ifelse(perm, paste0("_perm_", seed), ""), "_mirror_raw.pdf"))
        
        label <- "Yes"
        if (analysis == "omt") {
            label="OMT"
        }
        if (analysis == "gender") {
            label="Female"
        }
        if (analysis == "race"){
            label="Non-White"
        }
        if (analysis == "tmb_type") {
            label="High"
        }
        if (analysis == "RECIST") {
            label="PD"
        }
        print(hist.mirror2(data.plot, label))
        dev.off()
    }
    wt <- ans$W

    source("check_balance.R")
    index.tr <- which(clinical.common$Z==1)
    index.ctr <- which(clinical.common$Z==0)

    cutoff=0.1
    for (i in 2:(ncol(clinical.common)-1)) {
        print(paste0("check ", colnames(clinical.common)[i]))
        if (nrow(clinical.common[which(clinical.common[,i] == 1),]) >= 3) {

            std.diff <- check.balance(
                index.tr, index.ctr, clinical.common[,i], wt, printout=TRUE
            )
            if (std.diff > cutoff) {
                print(paste0("Fail: standardized difference= ", std.diff))
                print("=======================================")

            } else {
                print("Pass!")
            }
        } else {
            print("too small sample size")
        }
    }

    molecular.pvalues <- c()
    molecular.coefs <- c()
    molecular.0 <- c()
    molecular.1 <- c()
    molecular.0.w <- c()
    molecular.1.w <- c()

    if (is.survival) {
        require(survival)
        for (i in c(1)) {
            tmp <- sapply(split(as.numeric(molecular.common[,i]), clinical.common$Z), mean,na.rm=T)
            molecular.0 <- c(molecular.0, tmp[1])
            molecular.1 <- c(molecular.1, tmp[2])

            molecular.0.w <- c(molecular.0.w, sum(as.numeric(molecular.common[index.ctr,i])*wt[index.ctr],na.rm=T)/sum(wt[index.ctr],na.rm=T))
            molecular.1.w <- c(molecular.1.w, sum(as.numeric(molecular.common[index.tr,i])*wt[index.tr],na.rm=T)/sum(wt[index.tr],na.rm=T))
            #print(i)
            if (i%%1000 == 0) {
                print(i)
            }



            print("-------- cox result ------")
            print(summary(coxph(Surv(molecular.common[,"time"], molecular.common[,"status"])~clinical.common$Z, weights=wt),na.rm=T)$coef)
            pvalue <- summary(coxph(Surv(molecular.common[,"time"], molecular.common[,"status"])~clinical.common$Z, weights=wt),na.rm=T)$coef[5]
            coef <- summary(coxph(Surv(molecular.common[,"time"], molecular.common[,"status"])~clinical.common$Z, weights=wt),na.rm=T)$coef[1]
            if (class(pvalue) == "try-error") {
                pvalue <- NA
                coef <- NA
            }
        }
        molecular.pvalues <- c(molecular.pvalues, pvalue)
        molecular.coefs <- c(molecular.coefs, coef)

        print(paste0("coef:", molecular.coefs))
        print(paste0("pvalue:", molecular.pvalues))
        feature <- "time-stauts"
    } else {
        for (i in 1:ncol(molecular.common)) {
            tmp <- sapply(
                split(as.numeric(molecular.common[,i]), clinical.common$Z),
                mean, na.rm=TRUE
            )
            molecular.0 <- c(molecular.0, tmp[1])
            molecular.1 <- c(molecular.1, tmp[2])

            molecular.0.w <- c(molecular.0.w, sum(as.numeric(molecular.common[index.ctr,i])*wt[index.ctr],na.rm=T)/sum(wt[index.ctr],na.rm=T))
            molecular.1.w <- c(molecular.1.w, sum(as.numeric(molecular.common[index.tr,i])*wt[index.tr],na.rm=T)/sum(wt[index.tr],na.rm=T))

            if (i%%1000 == 0) {
                print(i)
            }

            if (is.continuous) {

                pvalue <- try(summary(lm(molecular.common[,i]~clinical.common$Z, weights=wt))$coef[2,4])
                coef <- try(summary(lm(molecular.common[,i]~clinical.common$Z, weights=wt))$coef[2,1])
                pvalue <- summary(lm(molecular.common[,i]~clinical.common$Z, weights=wt))$coef[2,4]
                if (class(pvalue)=="try-error") {
                    pvalue <- NA
                    coef <- NA
                }
            } else {

                pvalue <- summary(glm(as.numeric(molecular.common[,i])~clinical.common$Z, family=binomial, weights=wt),na.rm=T)$coef[2,4]
                coef <- summary(glm(as.numeric(molecular.common[,i])~clinical.common$Z, family=binomial, weights=wt),na.rm=T)$coef[2,1]
                if (class(pvalue) == "try-error") {
                    pvalue <- NA
                    coef <- NA
                }
            }
            molecular.pvalues <- c(molecular.pvalues, pvalue)
            molecular.coefs <- c(molecular.coefs, coef)
        }
        feature <- colnames(molecular.common)
    }
    molecular.fdr <- p.adjust(molecular.pvalues,"fdr")

    return(list(
        feature=feature, coef=molecular.coefs,
        pvalue=molecular.pvalues, fdr=molecular.fdr,
        mean.0= molecular.0, mean.1= molecular.1,
        mean.0.w=molecular.0.w, mean.1.w=molecular.1.w
    ))
}

summarize.p <- function(
    molecular.data, molecular.result, cutoff=0.05, print=FALSE
) {
    pvalue <- molecular.result$pvalue
    fdr <- molecular.result$fdr
    coef <- molecular.result$coef
    mean.0 <- molecular.result$mean.0
    mean.1 <- molecular.result$mean.1
    mean.0.w <- molecular.result$mean.0.w
    mean.1.w <- molecular.result$mean.1.w
    signif.index <- which(pvalue<cutoff)
    print(paste("Features with p value <", cutoff, "=", length(signif.index)))
    if (print) {
        print(cbind(

            molecular.result[["feature"]][signif.index],
            signif(pvalue[signif.index], 3),
            signif(fdr[signif.index], 3),
            signif(coef[signif.index], 3),
            signif(mean.0[signif.index], 3),
            signif(mean.1[signif.index], 3)
       ))
    }
    return(list(

        feature.sig=molecular.result[["feature"]][signif.index],
        pvalue.sig=pvalue[signif.index], fdr.sig= fdr[signif.index],
        n.sig=length(signif.index), coef.sig=coef[signif.index],
        mean0.sig=mean.0[signif.index], mean1.sig=mean.1[signif.index],
        mean0.sig.w=mean.0.w[signif.index], mean1.sig.w=mean.1.w[signif.index]
    ))
}


rm.zero.col <- function(molecular.data) {

    genes.raw <- colnames(molecular.data)
    genes <- unlist(lapply(genes.raw, function(x) unlist(strsplit(x,"[|]"))[1]))
    ids <- as.numeric(unlist(lapply(genes.raw, function(x) unlist(strsplit(x,"[|]"))[2])))
    if(FALSE)
    {
          rm1 <- grep("?", genes, fixed=TRUE)
          molecular.data <- molecular.data[,-rm1]
          print(paste(length(rm1),"genes were removed."))
    }

    if (FALSE) {
        gene.data <- read.table("protein_coding_gene_annotation.tsv",header=TRUE, sep="\t", quote="")
        name.id <- paste(data$Gene, data$ID, sep="|")
        keep.index <- match(gene.data$ID, ids)
        if (length(which(is.na(keep.index))) > 0) {
            stop("Some protein coding genes are not found. Double check!")
        }
        print(paste("After remove non-coding, total genes:", length(keep.index)))


        molecular.data <- molecular.data[,keep.index]
    }

    zero.col <- which(apply(molecular.data, 2, sd)==0)
    if (length(zero.col) > 0) {
        molecular.data <- molecular.data[, -zero.col]
        print(paste(length(zero.col),"columns were removed because of no variation."))
    }

    return(molecular.data)
}


write.summary <- function(summary, cancer, analysis, type) {
    file <- paste0(cancer,"_",analysis,"_summary_",type,".txt")
    data <- data.frame(
        summary$feature.sig, summary$fdr.sig, summary$coef.sig,
        summary$mean0.sig, summary$mean1.sig, summary$mean0.sig.w,
        summary$mean1.sig.w
    )
    if (analysis == "gender") {
        colnames(data) <- c("feature","fdr","coef","mean_MALE", "mean_FEMALE", "mean_MALE_weighted", "mean_FEMALE_weighted")
    }
    if(analysis == "race") {
        colnames(data) <- c("feature","fdr","coef","mean_nonWHITE", "mean_WHITE","mean_nonWHITE_weighted", "mean_WHITE_weighted")
    }
    if (analysis == "omt") {
        colnames(data) <- c("feature","fdr","coef","mean_nonOMT", "mean_OMT","mean_nonOMT_weighted", "mean_OMT_weighted")
    }
    if (analysis=="tmb_type") {
        colnames(data) <- c("feature","fdr","coef","mean_nonHigh", "mean_High","mean_nonHigh_weighted", "mean_High_weighted")
    }
    write.table(
        data, file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE
    )
}


write.result <- function(result, cancer, analysis, type, fdr.cutoff=0.05) {
    file <- paste0(cancer, "_", analysis, "_result_", type,".txt")
    data <- data.frame(
        result$feature, result$pvalue, result$fdr, result$coef,
        result$mean.0, result$mean.1, result$mean.0.w, result$mean.1.w
    )
    if (analysis == "gender") {
        colnames(data) <- c("feature","pvalue","fdr","coef","mean_MALE", "mean_FEMALE", "mean_MALE_weighted", "mean_FEMALE_weighted")
    }
    if(analysis == "race") {
        colnames(data) <- c("feature","fdr","coef","mean_nonWHITE", "mean_WHITE","mean_nonWHITE_weighted", "mean_WHITE_weighted")
    }
    if (analysis == "omt") {
        colnames(data) <- c("feature","fdr","coef","mean_nonOMT", "mean_OMT","mean_nonOMT_weighted", "mean_OMT_weighted")
    }
    if (analysis == "tmb_type") {
        colnames(data) <- c("feature","fdr","coef","mean_nonHigh", "mean_High","mean_nonHigh_weighted", "mean_High_weighted")
    }
    write.table(
        data, file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE
    )
}


plot.perm <- function(perm, n, cancer, analysis, type, cutoff) {
    p <- length(which(perm>= n))/length(perm)
    print(paste0("Permutation P-value = ", p))
    print(paste0(type,": n.sig = ", n))
    print(paste0("Median perm n.sig = ", median(perm)))
    file <- paste0(cancer, "_", analysis, "_perm_", type, "_", cutoff, ".pdf")
    pdf(file, width=3, height=3)
    par(mgp=c(2,1,0))
    if (type=="mut") {main="Mutation"}
    if (type=="cnv") {main="SCNA"}
    if(type=="methy") {main="Methy"}
    if(type=="mRNAseq") {main="mRNA"}
    if(type=="miRNA") {main="miRNA"}
    if(type=="rppa") {main="Protein"}
    if (n > max(perm)) {
        myhist <- hist(perm, xlim=c(0, max(max(perm),n)), main="", xlab="# Genes" )

    } else {
        myhist <- hist(perm, main="", xlab="# Genes")
    }
    text(n, max(myhist$counts),paste("p-value =", p), pos=ifelse(n> 0.5* max(perm),2,4), col=ifelse(p<=0.05, "red","blue"), cex=1.1, xpd=T)
    abline(v=n, col=ifelse(p<=0.05, "red","blue"), lty=2, lwd=2)
    dev.off()
}


perm.cal <- function(
    cancer, analysis, type, molecular.pri, cutoff=0.05, seedV=1:100
){
    print(paste("Calculate permutation for", cancer,":",type,":"))

    result <- get(paste0(type, ".result"))
    summary <- summarize.fdr(molecular.pri, result, cutoff=cutoff)
    n <- summary$n.sig
    perm <- c()
    print("For permutation:")
    for (seed in seedV) {

        perm.result <- perm.summary <- c()
        if (type == "mut") {
            load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
        } else if (type == "pre") {
            load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_pre_",seed,".RData", sep=""))
        } else {
            load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
        }
        perm.result <- get(paste("perm.",type,".result", sep=""))
        perm.summary <- summarize.fdr(molecular.pri, perm.result, cutoff=cutoff )
        perm <- c(perm, perm.summary$n.sig)
        do.call(rm, list(paste0("perm.",type,".result")))
    }

    plot.perm(perm, n, cancer, analysis, type, cutoff)
}
