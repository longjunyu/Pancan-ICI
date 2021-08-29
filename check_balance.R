check.balance <- function(index.tr, index.ctr, X1, wt, printout=FALSE) {
    mean.tr <- (
        sum(wt[index.tr]*X1[index.tr], na.rm=TRUE) /
        sum(wt[index.tr], na.rm=TRUE)
    )
    sd.tr <- (
        sum(wt[index.tr]*(X1[index.tr]-mean(X1[index.tr]))^2, na.rm=TRUE) /
        (sum(wt[index.tr], na.rm=TRUE)-1)
    )
    mean.ctr <- (
        sum(wt[index.ctr]*X1[index.ctr], na.rm=TRUE) /
        sum(wt[index.ctr], na.rm=TRUE)
    )
    sd.ctr <- (
        sum(wt[index.ctr]*(X1[index.ctr]-mean(X1[index.ctr]))^2, na.rm=TRUE) /
        (sum(wt[index.ctr], na.rm=TRUE)-1)
    )

    std.diff <- abs(mean.tr-mean.ctr)/sqrt((sd.tr+sd.ctr)/2)
    if (printout) {
        print(std.diff)
    }

    return(std.diff)
}
