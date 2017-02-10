#' Function to generate a 1D prior.
#' For use with GRS and additive genetic effects
#'
#' @param b1.range input
#' @param spac.b1 input
#' @param mn.b1 input
#' @param sd.b1 input
#' @param pi.1 input
#' @export
#' @examples
#' calculate_1d_prior()
#'

calculate_1d_prior <- function (
    b1.range = c(-6, 6),
    spac.b1 = 0.01,
    mn.b1 = 0,
    sd.b1 = 1,
    pi.1 = 0.001
) {
    b1.grid <- seq(b1.range[1], b1.range[2], by = spac.b1)
    prior <- (b1.grid^2/sd.b1^2) * dnorm(b1.grid, mn.b1, sd.b1^2)

    w0 <- which.min(abs(b1.grid))
    prior[w0] <- 0
    s1 <- sum(prior)
    prior <- prior * pi.1 / s1
    prior[w0] <- 1 - pi.1
    
    return(list(b.grid = b1.grid, prior = prior, val.0 = w0))
}


#' Helper function
#'
#' @param b0 input
#' @param b1 input
#' @param x0 input
#' @param x1 input
#' d.llk.d.b0_1d()
#'

d.llk.d.b0_1d <- function(b0,b1,x0,x1) {
    sum(1/(1+exp(b0+b1*x1))) - sum(1/(1+exp(-b0-b1*x0)))
}


#' Helper function
#'
#' @param b0 input
#' @param b1 input
#' @param x input
#' p1.func() 
#'

p1.func <- function(b0,b1,x) {
    1/(1+exp(-(b0+b1*x)))
}


#' Helper function
#'
#' @param b0 input
#' @param b1 input
#' @param x input
#' p0.func() 
#'

p0.func <- function(b0,b1,x) {
    1 - p1.func(b0,b1,x)
}


#' Helper function
#'
#' @param b0 input
#' @param b1 input
#' @param x0 input
#' @param x1 input
#' llk_1d() 
#'

llk_1d <- function(b0,b1,x0,x1) {
    sum(log(p1.func(b0,b1,x1))) + sum(log(p0.func(b0,b1,x0)))    
}


#' Calculate LLK over the grid
#'
#' @param prior input
#' @param x0 input
#' @param x1 input
#' @param scaled input
#' @export
#' @examples
#' calculate_1d_llk_grid_scaled() 
#'

calculate_1d_llk_grid_scaled <- function(
    prior,x0,x1,scaled=TRUE
) {

    to.fill <- 0
    if(scaled)
        to.fill <- 1
    
    op <- array(to.fill,c(length(prior$b.grid)))                
    b0.est <- array(0,dim(op))

    for(i in 1:length(prior$b.grid)) {
        tmp <- uniroot(
            d.llk.d.b0_1d,
            c(-100,100),
            b1=prior$b.grid[i],
            x0=x0,
            x1=x1);
        b0.est[i] <- tmp$root
        op[i] <- llk_1d(b0=b0.est[i],b1=prior$b.grid[i],x0,x1)
    }

    mx <- max(op)
    if(scaled)
        op <- exp(op - mx)

    return(list(op=op, b0.est=b0.est, lmx=mx))
}


#' Wrapper function to use with parallel package
#'
#' @param d input
#' @export
#' @examples
#' calculate.llk.grid.wrapper()
#'

calculate.llk.grid.wrapper <- function(
    d
) {
    
    library(TreeWAS)

    x0 <- d$x0
    x1 <- d$x1
    prior <- d$prior
    
    tmp <- calculate_1d_llk_grid_scaled(
        d$prior,
        d$x0,
        d$x1,
        scaled=TRUE
    )
    
    return(tmp)
}


#' Function to sum over the LLK grid
#'
#' @param prior input
#' @param llk.surf input
#' @param scaled input
#' @export
#' @examples
#' calculate.integrated_1d_llk_scaled()
#'

calculate.integrated_1d_llk_scaled <- function(
    prior,
    llk.surf,
    scaled = TRUE
) {
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


#' Estimate posterior probability and get MAP estimate
#'
#' @param forward input
#' @param backward input
#' @param prior input
#' @param id input
#' @param plot input
#' @param return.ci input
#' @param verbose input
#' @param ci.level input
#' @param log.plot input
#' @export
#' @examples
#' get.posterior.node_1d()
#'

get.posterior.node_1d <- function(
    forward,
    backward,
    prior,
    id = 1,
    plot = FALSE, 
    return.ci = TRUE,
    verbose = FALSE,
    ci.level = 0.95,
    log.plot = TRUE
) {
    tmp <- forward[[id]]$op * backward[[id]]$op
    tmp <- tmp/sum(tmp)
    post.null <- tmp[null.id]
    post.active <- 1 - post.null
    tmp[null.id] <- 0
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
        rg <- range(prior$b.grid[inds[, 1]])

        if (verbose) 
            cat("\nCI b1(", ci.level, ") = ", paste(rg, collapse = " - "), 
                sep = "")      
    }
    if (plot) {
        if (log.plot) 
            tmp <- log(tmp)
        plot(x = prior$b.grid, y = tmp, 
             main = paste("Node", id), xlab = "B1", ylab = "log(post)")
    }
    if (verbose) 
        cat("\n\n")

    out <- data.frame(
        max_b = prior$b.grid[mx[1]],
        summed_llk = log(sum(forward[[id]]$op * backward[[id]]$op)) + 
            forward[[id]]$lmx + backward[[id]]$lmx,
        b_ci_lhs = rg[1], 
        b_ci_rhs = rg[2],
        POST_ACTIVE = post.active)

    return(out)
}
