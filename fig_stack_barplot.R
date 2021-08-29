get_labely <- function(indf, var, group, adjust=0.5) {
    newdf <- indf
    labely <- NULL
    groups <- levels(newdf[, group])
    for (i in 1:length(groups)) {
        subdf <- newdf[newdf[, group] == groups[i],]
        labely <- c(
            labely, 1 - (cumsum(subdf[, var]) - subdf[, var] * adjust)
        )
    }
    return(labely)
}

fig_stack_barplot <- function(
    indf,
    var,
    group,
    subgroup,
    method="auto",
    n_fisher=10000,
    print_output=TRUE,
    file_output=TRUE,
    fig_title=NULL,
    fig_xlab_title="Group",
    fig_ylab_title="Var",
    fig_subgroup_color=c("#E3191C", "#309EDE"),
    fig_ylab_breaks=seq(0, 1, 0.2),
    fig_ylab_labels=seq(0, 1, 0.2),
    fig_comparison=TRUE,
    fig_filepath=".",
    fig_filename="fig_barplot",
    fig_filetype="tiff",
    fig_width=20,
    fig_height=20
) {
    if (var %in% c("pct", "sum")) {
        print("var, group, subgroup can not be name to pct, sum")
        break
    }

    newdf <- indf[, c(var, group, subgroup)]
    colnames(newdf) <- c("var", "group", "subgroup")

    mtx <- matrix(
        newdf[, "var"],
        ncol=length(levels(newdf[, "group"]))
    )

    n_group <- length(levels(newdf[, "group"]))
    n_subgroup <- length(levels(newdf[, "subgroup"]))
    if (n_group == 2 & n_subgroup == 2 & sum(mtx) < n_fisher) {
        model <- fisher.test(mtx)
        pvalue <- model[["p.value"]]
        pvalue_label <- ifelse(
            pvalue < 0.001, "P < 0.001",
            paste0("P = ", sprintf("%.3f", pvalue))
        )
        auto_title <- paste0("Fisher's exact test ", pvalue_label)
    } else {
        model <- chisq.test(mtx)
        pvalue <- model[["p.value"]]
        pvalue_label <- ifelse(
            pvalue < 0.001, "P < 0.001",
            paste0("P = ", sprintf("%.3f", pvalue))
        )
        auto_title <- paste0("Chi-squared test ", pvalue_label)
    }

    if (print_output) {
        print(mtx)
        print(model)
    }

    if (n_group >= 2) {
        combdf <- as.data.frame(t(combn(n_group, 2)))
        combdf[, "V3"] <- combdf[, 2] - combdf[, 1]
        combdf[, "V4"] <- ifelse(combdf[, "V3"] %% 2 == 1, 1, -1)
        combdf[, "V5"] <- combdf[, "V1"] * combdf[, "V4"]
        print(combdf)
        combdf <- combdf[order(combdf[, 3], combdf[, 5]),]

        pvalues <- NULL
        for (i in 1:nrow(combdf)) {
            mtx2 <- mtx[, as.numeric(combdf[i, 1:2])]
            if (n_subgroup == 2 & sum(mtx2) < n_fisher) {
                model <- fisher.test(mtx2)
            } else {
                model <- chisq.test(mtx2)
            }
            pvalue <- model[["p.value"]]
            pvalues <- c(pvalues, pvalue)

            if (print_output) {
                print(mtx)
                print(model)
            }
        }
        pvalue_labels <- ifelse(
            pvalues < 0.001, "P < 0.001",
            paste0("P = ", sprintf("%.3f", pvalue))
        )

    } else {
        pvalue_labels <- pvalue_label
    }

    newdf2 <- aggregate(var ~ group, newdf, sum)
    colnames(newdf2)[2] <- "sum"
    newdf <- merge(newdf, newdf2, by="group", all.x=TRUE)
    newdf <- newdf[order(newdf[, "group"], newdf[, "subgroup"]),]
    newdf[, "pct"] <- newdf[, "var"] / newdf[, "sum"]
    labelx <- newdf[, "group"]
    labely <- get_labely(newdf, "pct", "group")
    labels <- paste0(sprintf("%.1f", newdf[, "pct"] * 100), "%")

    plot <- ggplot(data=newdf, aes(x=group, y=pct, fill=subgroup)) +
        geom_bar(stat="identity", width=0.7) +
        geom_text(aes(x=labelx, y=labely, label=labels)) +
        scale_fill_manual(values=fig_subgroup_color) +
        xlab(fig_xlab_title) +
        scale_y_continuous(
            fig_ylab_title, limit=c(0, 1),
            breaks=fig_ylab_breaks, labels=fig_ylab_labels
        ) +
        theme_bw() +
        theme(
            legend.position="top",
            legend.title=element_blank(),
            panel.grid=element_blank(),
            axis.line=element_line(color="black"),
            panel.border=element_rect(color=NA, fill=NA),
            plot.title=element_text(hjust=0.5)
        )

    if (fig_comparison) {
        for (i in 1:nrow(combdf)) {
            plot <- plot +
                geom_segment(
                    x=combdf[i, 1], xend=combdf[i, 2],
                    y=1 + 0.1 * i, yend=1 + 0.1 * i
                ) +
                geom_segment(
                    x=combdf[i, 1], xend=combdf[i, 1],
                    y=1 + 0.1 * i - 0.03, yend=1 + 0.1 * i
                ) +
                geom_segment(
                    x=combdf[i, 2], xend=combdf[i, 2],
                    y=1 + 0.1 * i - 0.03, yend=1 + 0.1 * i
                ) +
                geom_text(
                    x=combdf[i, 1] + (combdf[i, 2] - combdf[i, 1]) * 0.5,
                    y=1 + 0.1 * i + 0.05,
                    label=pvalue_labels[i]
                )
        }
        plot <- plot +
            scale_y_continuous(
                fig_ylab_title, limit=c(0, 1 + 0.1 * nrow(combdf) + 0.05),
                breaks=fig_ylab_breaks, labels=fig_ylab_labels
            )
    }
    if (is.null(fig_title)) {
        plot <- plot + ggtitle(auto_title)
    } else {
        plot <- plot + ggtitle(fig_title)
    }
    if (print_output) {
        print(plot)
    }

    if (file_output) {
        if ("jpg" %in% fig_filetype) {
            outfile <- paste0(fig_filepath, "/", fig_filename, ".jpg")
            jpeg(outfile, width=fig_width, height=fig_height, unit="cm",res=350)
            print(plot)
            dev.off()
        }
        if ("tiff" %in% fig_filetype) {
            outfile <- paste0(fig_filepath, "/", fig_filename, ".tiff")
            tiff(
                outfile, width=fig_width, height=fig_height, unit="cm",
                res=350, compression="lzw+p"
            )
            print(plot)
            dev.off()
        }
        if ("pdf" %in% fig_filetype) {
            outfile <- paste0(fig_filepath, "/", fig_filename, ".pdf")
            pdf(
                outfile, onefile=FALSE,
                width=fig_width/2.54, height=fig_height/2.54
            )
            print(plot)
            dev.off()
        }
    }
}
