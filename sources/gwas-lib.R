#-------------------------------------------------------------
# Adjust both pValues and threshold (for Bonferroni)
# Calculate threshold to decide SNPs significance
#-------------------------------------------------------------
adjustPValues <- function (level, pValues, method="FDR") 
{
	m <- length(pValues)
	if (method=="Bonferroni") {
		threshold = -log10(level/m)
	} else if (method=="FDR") {
		tmp <- cbind(pValues,calculateQValue (pValues))
		tmp <- tmp[order(tmp[,2]),]
		if (tmp[1,2] > level) {
			threshold <- -log10(tmp[1,1])*1.2
		} else {
			k <- max(which(tmp[,2] < level))
			threshold <- -log10(mean(tmp[k:(k+1),1]))
		}
	}

	return (list (threshold=threshold, pValues=pValues))
}


calculateQValue <- function(p) {
        smooth.df = 3
        if (min(p) < 0 || max(p) > 1) {
            print("ERROR: p-values not in valid range.")
            return(0)
        }
        lambda = seq(0, 0.9, 0.05)
        m <- length(p)
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }

        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
        pi0 <- min(pi0, 1)
        if (pi0 <= 0) {
            print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
            return(0)
        }
        u <- order(p)
        qvalue.rank <- function(x) {
            idx <- sort.list(x)
            fc <- factor(x)
            nl <- length(levels(fc))
            bin <- as.integer(fc)
            tbl <- tabulate(bin)
            cs <- cumsum(tbl)
            tbl <- rep(cs, tbl)
            tbl[idx] <- tbl
            return(tbl)
        }
        v <- qvalue.rank(p)
        qvalue <- pi0 * m * p/v
        qvalue[u[m]] <- min(qvalue[u[m]], 1)
        for (i in (m - 1):1) {
            qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                1)
        }
        return(qvalue)
}
