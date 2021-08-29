library(gdata)
library(gtools)
library(gmodels)
library(gplots)
library(formula.tools)
library(Hmisc)
library(Matrix)

GIPW.std.omega <- function(
    dat, form.ps, weight, trt.est=TRUE, delta=0.002, K=4
) {

    n <- nrow(dat)

    out.ps <- ps.model(dat, as.formula(form.ps))
    ps.hat <- out.ps$ps.hat
    dat$ps.hat <- ps.hat
    beta.hat <- as.numeric(coef(out.ps$fm))

    Q <- dat$Z*dat$ps.hat+(1-dat$Z)*(1-dat$ps.hat)
    omega <- sapply(dat$ps, calc.omega, weight=weight, delta=delta, K=K)
    W <- omega/Q
    dat$W <- W

    if (!(trt.est)) {
        ans <- list(weight=weight, W=W, ps=ps.hat)
        return(ans)
    } else {

        mu1.hat <- sum(dat$W*dat$Z*dat$Y)/sum(dat$W*dat$Z)
        mu0.hat <- sum(dat$W*(1-dat$Z)*dat$Y)/sum(dat$W*(1-dat$Z))
        est <- mu1.hat-mu0.hat

        Amat <- Bmat <- 0
        for (i in 1:n) {

            Xi <- as.numeric(dat[i, c("X0", names(coef(out.ps$fm))[-1])])
            Zi <- dat[i, "Z"]
            Yi <- dat[i, "Y"]
            ei <- calc.ps.Xbeta(Xi, beta.hat)
            ei.deriv1 <- calc.ps.deriv1(Xi, beta.hat)
            ei.deriv2 <- calc.ps.deriv2(Xi, beta.hat)
            omega.ei <- omega[i]
            omegaei.deriv <- omega.derive.ei(ei, weight, delta, K)
            Qi <- Q[i]
            Qi.deriv <- 2*Zi-1

            phi <- c(
                Zi*(Yi-mu1.hat)*omega.ei/Qi,
                (1-Zi)*(Yi-mu0.hat)*omega.ei/Qi,
                (Zi-ei)/(ei*(1-ei))*ei.deriv1
            )

            Bmat <- Bmat+outer(phi, phi)

            first.row <- c(
                -Zi*omega.ei/Qi, 0,
                Zi*(Yi-mu1.hat)*ei.deriv1*
                (Qi*omegaei.deriv-omega.ei*Qi.deriv)/Qi^2
            )

            second.row <- c(
                0, -(1-Zi)*omega.ei/Qi,
                (1-Zi)*(Yi-mu0.hat)*ei.deriv1*
                (Qi*omegaei.deriv-omega.ei*Qi.deriv)/Qi^2
            )

            tmp0 <- matrix(0, nrow=length(beta.hat), ncol=2)
            tmp1 <- -ei*(1-ei)*colVec(Xi)%*%Xi
            third.row <- cbind(tmp0, tmp1)

            phi.deriv <- rbind(first.row, second.row, third.row)
            Amat <- Amat+phi.deriv
        }

        Amat <- Amat/n
        Bmat <- Bmat/n
        Amat.inv <- solve(Amat)
        var.mat <- (Amat.inv %*% Bmat %*% t(Amat.inv))/n
        tmp1 <- c(1, -1, rep(0, length(beta.hat)))
        var.est <- rowVec(tmp1) %*% var.mat %*% colVec(tmp1)
        std <- sqrt(as.numeric(var.est))
        CI.lower <- est-1.96*std
        CI.upper <- est+1.96*std

        ans <- list(
            weight=weight, est=est, std=std,
            CI.lower=CI.lower, CI.upper=CI.upper, W=W, ps=ps.hat
        )
        return(ans);
    }
}

rowVec <- function(x) {
    t(x)
}

colVec <- function(x) {
    t(t(x))
}

ps.model <- function(dat, form) {

    fm <- glm(form, data=dat, family=binomial(link="logit"))
    ps.hat <- as.numeric(predict(fm, newdata=dat, type="response"))

    ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999)
    return(list(ps.hat=ps.hat, fm=fm))
}

calc.omega <- function(ps, weight, delta, K) {

    ans <- 0
    if (weight == "IPW") {
        ans <- 1
    } else if (weight == "MW") {
        ans <- calc.omega.MW(ps, delta)
    } else if (weight == "ATT") {
        ans <- ps
    } else if (weight == "ATC") {
        ans <- 1-ps
    } else if (weight == "OVERLAP") {
        ans <- 4*ps*(1-ps)
    } else if (weight == "TRAPEZOIDAL") {
        ans <- calc.omega.trapzd(ps=ps, delta=delta, K=K)
    } else {
        stop("Error in calc.omega: weight method does not exist!")

    }

    ans
}

calc.ps.Xbeta <- function(Xmat, beta) {

    Xmat <- as.matrix(Xmat)
    tmp <- as.numeric(rowVec(beta)%*%Xmat)
    tmp <- exp(tmp)
    names(tmp) <- NULL
    return(tmp/(1+tmp))
}

calc.ps.deriv1 <- function(Xmat, beta) {

    Xmat <- as.matrix(Xmat)
    tmp.ps <- calc.ps.Xbeta(Xmat, beta)

    ans <- tmp.ps*(1-tmp.ps)*t(Xmat)
    names(ans) <- rownames(ans) <- NULL
    return(t(ans))
}

calc.ps.deriv2 <- function(Xi, beta) {

    Xi <- colVec(Xi)
    tmp.ps <- calc.ps.Xbeta(Xi, beta)
    tmp.deriv1 <- calc.ps.deriv1(Xi, beta)
    ans <- Xi%*%rowVec(tmp.deriv1)
    names(ans) <- rownames(ans) <- NULL
    ans <- (1-2*tmp.ps)*ans
    return(ans)
}

calc.omega.MW <- function(ps, delta) {

    ans <- 0
    if (ps <= 0.5-delta) {
        ans <- 2*ps
    } else if (ps >= 0.5+delta) {
        ans <- 2*(1-ps)
    } else {
        ans <- approx.omega.MW(ps, delta)
    }

    ans
}

approx.omega.MW <- function(ps, delta) {

    A <- solve.A.MW(delta)
    ans <- rowVec(A)%*%c(1, ps, ps^2, ps^3)

    ans
}

solve.A.MW <- function(delta) {

    if (delta < 0.00001) {
        stop("*** ERROR in solve.a: delta too small ***")
    }

    tmp1 <- 0.5-delta
    tmp2 <- 0.5+delta

    D <- matrix(c(
        1, tmp1, tmp1^2, tmp1^3,
        0, 1, 2*tmp1, 3*tmp1^2,
        1, tmp2, tmp2^2, tmp2^3,
        0, 1, 2*tmp2, 3*tmp2^2
    ), ncol=4, nrow=4, byrow=TRUE)
    C <- 2*c(tmp1, 1, tmp1, -1)
    A <- solve(D)%*%C

    A
}

calc.omega.trapzd <- function(ps, delta, K) {

    ans <- 0
    if ((0 < ps) & (ps <= 1/K-delta)) {
        ans <- K*ps
    } else if ((1/K+delta <= ps) & (ps <= 1-1/K-delta) ) {
        ans <- 1
    } else if ((1-1/K+delta <= ps) & (ps < 1)) {
        ans <- K*(1-ps)
    } else {
        ans <- approx.omega.trapzd(ps, delta, K)
    }

    ans
}

approx.omega.trapzd <- function(ps, delta, K) {

    A <- 0
    if ((1/K-delta < ps) & (ps < 1/K+delta)) {
        A <- solve.A.trapzd1st(delta, K)
    } else {
        A <- solve.A.trapzd2nd(delta, K)
    }
    ans <- rowVec(A)%*%c(1, ps, ps^2, ps^3)

    ans
}

solve.A.trapzd1st <- function(delta=delta, K=K) {

    if (delta < 0.00001) {
        stop("*** ERROR in solve.a: delta too small ***")
    }

    tmp1 <- 1/K-delta
    tmp2 <- 1/K+delta

    D <- matrix(c(
        1, tmp1, tmp1^2, tmp1^3,
        0, 1, 2*tmp1, 3*tmp1^2,
        1, tmp2, tmp2^2, tmp2^3,
        0, 1, 2*tmp2, 3*tmp2^2
    ), ncol=4, nrow=4, byrow=TRUE)
    C <- 2*c(K*tmp1, K, 1, 0)
    A <- solve(D)%*%C

    A
}

solve.A.trapzd2nd <- function(delta=delta, K=K) {

    if (delta < 0.00001) {
        stop("*** ERROR in solve.a: delta too small ***")
    }

    tmp1 <- 1-1/K-delta
    tmp2 <- 1-1/K+delta

    D <- matrix(c(
        1, tmp1, tmp1^2, tmp1^3,
        0, 1, 2*tmp1, 3*tmp1^2,
        1, tmp2, tmp2^2, tmp2^3,
        0, 1, 2*tmp2, 3*tmp2^2
    ), ncol=4, nrow=4, byrow=TRUE)
    C <- 2*c(1, 0, K*(1/K-delta), -K)
    A <- solve(D)%*%C

    A
}

omega.derive.ei <- function(ps=ei, weight=weight, delta=delta, K=K) {

    ans <- 0
    if (weight == "IPW") {
        ans <- 0
    } else if (weight == "ATT") {
        ans <- 1
    } else if (weight == "ATC") {
        ans <- -1
    } else if (weight == "OVERLAP") {
        ans <- 4*(1-2*ps)
    } else if (weight == "MW") {
        ans <- omega.derive.ei.MW (ps, delta)
    } else if (weight == "TRAPEZOIDAL") {
        ans <- omega.derive.ei.trapzd (ps, delta, K)
    } else {
        stop("User defined first-order derivative of omega function is not provided!")
    }
}

omega.derive.ei.MW <- function(ps, delta) {

    ans <- 0
    if ((0 < ps) & (ps <= 0.5-delta)) {
        ans <- 2*ps
    } else if ((0.5+delta <= ps) & (ps < 1)) {
        ans <- 2*(1-ps)
    } else {
        A <- solve.A.MW(delta)
        ans <- A[2]+2*A[3]*ps+3*A[4]*ps^2
    }

    ans
}

omega.derive.ei.trapzd <- function(ps, delta, K) {

    ans <- 0

    if ((0 < ps) & (ps <= 1/K-delta)) {
        ans <- K
    } else if ((1/K-delta < ps) & (ps < 1/K+delta)) {
        A <- solve.A.trapzd1st(delta=delta, K=K)
        ans <- A[2]+2*A[3]*ps+3*A[4]*ps^2
    } else if ((1/K+delta <= ps) & (ps <= 1-1/K-delta)) {
        ans <- 0
    } else if ((1-1/K+delta <= ps) & (ps<1)) {
        ans <- -K
    } else {
        A <- solve.A.trapzd2nd(delta=delta, K=K)
        ans <- A[2]+2*A[3]*ps+3*A[4]*ps^2
    }

    ans
}
