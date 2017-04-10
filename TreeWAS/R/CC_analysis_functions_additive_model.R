#' Tree Analysis.
#' Calculate llk grid for the additive genetic model.
#' Scaled version
#'
#' @param prior input
#' @param data input
#' @param do.plot input
#' @param debug input
#' @param scaled input
#'
#' @return a list
#'
#' @export
#'
#' @examples
#'
calculate.llk.grid_scaled_additive <- function(
    prior,
    data,
    do.plot=FALSE,
    debug=FALSE,
    scaled=TRUE
) {

    to.fill <- 0;
    if (scaled) to.fill <- 1;
    op <- array(to.fill,c(length(prior$b.grid)))
    b0.est <- array(0, dim(op));

    if (debug) {
        return(list(op=op, b0.est=b0.est, lmx=0));
    }

    for (i in 1:length(prior$b.grid)) {
        tmp <- uniroot(cc_d.llk_additive, c(-100, 100), b1=prior$b.grid[i], aff=data[4:6], unaf=data[1:3]);
        b0.est[i] <- tmp$root;
        op[i] <- cc_llk_additive(b0=b0.est[i], b1=prior$b.grid[i], aff=data[4:6], unaf=data[1:3]);
    }

    if (do.plot) plot(x=prior$b.grid, y=op, xlab="B.additive", ylab="LLK");

    mx <- max(op);
    if (scaled) op <- exp(op-mx);

    return(list(op=op, b0.est=b0.est, lmx=mx));
}


#' cc_d.llk_additive
#'
#' @param b0 input
#' @param b1 input
#' @param aff input
#' @param unaf input
#' @export
#' @examples
#'
cc_d.llk_additive <- function (b0, b1, aff, unaf)
{

    tmp1 <- (aff[1] + unaf[1])/(1 + exp(-b0))
    tmp2 <- (aff[2] + unaf[2])/(1 + exp(-(b0 + b1)))
    tmp3 <- (aff[3] + unaf[3])/(1 + exp(-(b0 + 2*b1)))
    tmp4 <- -sum(aff)

    return(tmp1+tmp2+tmp3+tmp4)
}


#' cc_llk_additive
#'
#' @param b0 input
#' @param b1 input
#' @param aff input
#' @param unaf input
#' @export
#' @examples
#'
cc_llk_additive <- function(b0, b1, aff, unaf)
{
    return(b0 * sum(aff) + b1 * aff[2] + 2 * b1 * aff[3] - (aff[1] +
        unaf[1]) * log(1 + exp(b0)) - (aff[2] + unaf[2]) * log(1 +
        exp(b0 + b1)) - (aff[3] + unaf[3]) * log(1 + exp(b0 +
        2*b1)))
}


#' calculate.integrated.llk_scaled_additive
#'
#' @param prior input
#' @param llk.surf input
#' @param scaled input
#' @export
#' @examples
#'
calculate.integrated.llk_scaled_additive <- function(prior, llk.surf, scaled = TRUE)
{
    if (scaled) {
        llk.sum <- sum(prior$prior * llk.surf$op)
        return(llk.sum)
    }
    else {
        mx <- max(llk.surf$op)
        op <- llk.surf$op - mx
        op.n <- exp(op)
        llk.sum <- log(sum(op.n * prior$prior))
        return(llk.sum + mx)
    }
}



#' get.posterior.node_additive
#'
#' @param forward TODO
#' @param backward TODO
#' @param prior TODO
#' @param id TODO
#' @param plot TODO
#' @param return.ci TODO
#' @param verbose TODO
#' @param ci.level TODO
#' @param log.plot TODO
#'
#' @export
#'
#' @examples
#'
get.posterior.node_additive <- function(forward, backward, prior, id = 1, plot = FALSE,
    return.ci = TRUE, verbose = FALSE, ci.level = 0.95, log.plot = TRUE)
{
    tmp <- forward[[id]]$op * backward[[id]]$op
    tmp <- tmp/sum(tmp)
    post.null <- tmp[null.id[1]]
    post.active <- 1 - post.null
    tmp[null.id[1]] <- 0
    tmp <- tmp/sum(tmp)
    mx <- arrayInd(which.max(tmp), dim(tmp))
    if (verbose)
        cat("\nNode ", id)
    if (verbose)
        cat("\nMax at b1 = ", prior$b.grid[mx[1]])
    if (verbose)
        cat("\nSummed LLK = ", log(sum(forward[[id]]$op * backward[[id]]$op)) +
            forward[[id]]$lmx + backward[[id]]$lmx)
    if (return.ci) {
        oo <- order(tmp, decreasing = T)
        cs <- cumsum(tmp[arrayInd(oo, dim(tmp))])
        w.ci <- oo[c(1, which(cs <= ci.level))]
        inds <- arrayInd(w.ci, dim(tmp))
        rg.1 <- range(prior$b.grid[inds[, 1]])

        if (verbose)
            cat("\nCI b1(", ci.level, ") = ", paste(rg.1, collapse = " - "),
                sep = "")
    }

    if (verbose)
        cat("\n\n")

    spac <- abs(prior$b.grid[2]-prior$b.grid[1])

    out <- data.frame(
        max_b1 = prior$b.grid[mx[1]],
        summed_llk = log(sum(forward[[id]]$op * backward[[id]]$op)) +
            forward[[id]]$lmx + backward[[id]]$lmx,
        b1_ci_lhs = rg.1[1] - spac/2,
        b1_ci_rhs = rg.1[2] + spac/2,
        POST_ACTIVE = as.numeric(post.active)
    )

    return(out)
}
