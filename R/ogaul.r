ogaul<-function (formula, d, data = NULL, na.action, ...) 
{
    d <- as.matrix(d)
    d1 <- d[1L]
    ogaule <- function(formula, d1, data = NULL, na.action) {
        cal <- match.call(expand.dots = FALSE)
        mat <- match(c("formula", "data", "na.action"), names(cal))
        cal <- cal[c(1L, mat)]
        cal[[1L]] <- as.name("model.frame")
        cal <- eval(cal)
        y <- model.response(cal)
        md <- attr(cal, "terms")
        x <- model.matrix(md, cal, contrasts)
        s <- t(x) %*% x
        xin <- solve(s)
        bb <- xin %*% t(x) %*% y
        I <- diag(NCOL(x))
        fd <- solve(s + I) %*% (s + d1 * I)
        bal <- (I - (1 - d1)^2 * solve(s + I) %*% solve(s + I)) %*% 
            bb
balve<-as.vector(bal)
j<-0
sumsq<-0
for (j in 1:NROW(balve))
{
sumsq=(balve[j])^2+sumsq
}
cval<-sumsq
        ev <- (t(y) %*% y - t(bb) %*% t(x) %*% y)/(NROW(x) - 
            NCOL(x))
        ev <- diag(ev)
ahat<-bal%*%t(bal)%*%solve(ev*xin+bal%*%t(bal))
ogaule<-ahat%*%bb
colnames(ogaule) <- c("Estimate")
        dbd <- ev*(ahat%*%xin%*%t(ahat))
        Standard_error <- sqrt(diag(abs(dbd)))
        dbt <- t(ogaule)
        dbd <- ev*(ahat%*%xin%*%t(ahat))
        sdbd_inv <- (sqrt(diag(abs(dbd))))^-1
        sdbd_inv_mat <- diag(sdbd_inv)
        if (NCOL(dbt) == 1L) 
            tbd <- dbt * sdbd_inv
        else tbd <- dbt %*% sdbd_inv_mat
        hggh <- t(tbd)
balm<-as.matrix(bal)
rval<-(1/cval)*balm%*%t(balm)
mse1 <- cval^2*ev*tr(ev*rval*solve(ev*xin+cval*rval)%*%xin%*%solve(ev*xin+cval*rval)%*%rval)+ev^2*t(bal)%*%solve(ev*I+cval*rval%*%s)%*%solve(ev*I+cval*rval%*%s)%*%bal
mse1<-as.vector(mse1)
        mse1 <- round(mse1, digits = 4L)
        names(mse1) <- c("MSE")
        tst <- t(2L * pt(-abs(tbd), df <- (NROW(x) - NCOL(x))))
        colnames(tst) <- c("p_value")
        colnames(hggh) <- c("t_statistic")
        ans1 <- cbind(ogaule, Standard_error, hggh, tst)
        ans <- round(ans1, digits = 4L)
        anw <- list(`*****Ordinary Generalized Almost Unbiased Liu Estimator*****` = ans, 
            `*****Mean square error value*****` = mse1)
        return(anw)
    }
    npt <- ogaule(formula, d1, data, na.action)
    plotogaul <- function(formula, d, data = NULL, na.action, ...) {
        i <- 0
        arr <- 0
        for (i in 1:NROW(d)) {
            if (d[i] < 0L) 
                d[i] <- 0L
            else d[i] <- d[i]
            if (d[i] > 1L) 
                d[i] <- 1L
            else d[i] <- d[i]
            ogaulm <- function(formula, d, data, na.action, ...) {
                cal <- match.call(expand.dots = FALSE)
                mat <- match(c("formula", "data", "na.action"), 
                  names(cal))
                cal <- cal[c(1L, mat)]
                cal[[1L]] <- as.name("model.frame")
                cal <- eval(cal)
                y <- model.response(cal)
                md <- attr(cal, "terms")
                x <- model.matrix(md, cal, contrasts)
                s <- t(x) %*% x
                xin <- solve(s)
                bb <- xin %*% t(x) %*% y
                I <- diag(NCOL(x))
                fd <- solve(s + I) %*% (s + d * I)
                bal <- (I - (1 - d)^2 * solve(s + I) %*% solve(s + 
                  I)) %*% bb
balve<-as.vector(bal)
j<-0
sumsq<-0
for (j in 1:NROW(balve))
{
sumsq=(balve[j])^2+sumsq
}
cval<-sumsq
                ev <- (t(y) %*% y - t(bb) %*% t(x) %*% y)/(NROW(x) - 
                  NCOL(x))
                ev <- diag(ev)
ahat<-bal%*%t(bal)%*%solve(ev*xin+bal%*%t(bal))
                dbd <- ev*(ahat%*%xin%*%t(ahat))
balm<-as.matrix(bal)
rval<-(1/cval)*balm%*%t(balm)
mse1 <- cval^2*ev*tr(ev*rval*solve(ev*xin+cval*rval)%*%xin%*%solve(ev*xin+cval*rval)%*%rval)+ev^2*t(bal)%*%solve(ev*I+cval*rval%*%s)%*%solve(ev*I+cval*rval%*%s)%*%bal
mse1<-as.vector(mse1)
                return(mse1)
            }
            arr[i] <- ogaulm(formula, d[i], data, na.action)
        }
        MSE <- arr
        parameter <- d
        pvl <- cbind(parameter, MSE)
        colnames(pvl) <- c("Parameter", "MSE")
        sval <- pvl
        return(sval)
    }
    paul <- plotogaul(formula, d, data, na.action)
    if (nrow(d) > 1L) 
        val <- paul
    else val <- npt
    val
}
