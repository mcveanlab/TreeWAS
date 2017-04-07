#' Function to generate a 1D non-local prior.
#' For use with GRS and additive genetic effects
#'
#' @param b1.range input
#' @param spac.b1 input
#' @param mn.b1 input
#' @param sd.b1 input
#' @param pi.1 input
#'
#' @return A list with objects b.grid, prior and val.
#'
#' @export
#'
#' @examples
#' calculate_1d_prior(b1.range=c(-4,4),spac.b1=0.01)
#'
calculate_1d_prior <- function (
    b1.range = c(-4,4),
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


#' Function to calculate d(log likelihood)/db0
#'
#' @param b0 input
#' @param b1 input
#' @param x0 input
#' @param x1 input
#'
#' @return A numeric object
#'
#' @examples
#'
d.llk.d.b0_1d <- function(b0=NULL,b1=NULL,x0=NULL,x1=NULL) {
    sum(1/(1+exp(b0+b1*x1))) - sum(1/(1+exp(-b0-b1*x0)))
}


#' Helper function
#'
#' @param b0 input
#' @param b1 input
#' @param x input
#'
#'
p1.func <- function(b0,b1,x) {
    1/(1+exp(-(b0+b1*x)))
}


#' Helper function
#'
#' @param b0 input
#' @param b1 input
#' @param x input
#'
#' @examples
#'
p0.func <- function(b0,b1,x) {
    1 - p1.func(b0,b1,x)
}


#' Likelihood function
#'
#' @param b0 input
#' @param b1 input
#' @param x0 input
#' @param x1 input
#'
#' @examples
#'
llk_1d <- function(b0,b1,x0,x1) {
    sum(log(p1.func(b0,b1,x1))) + sum(log(p0.func(b0,b1,x0)))
}


#' Calculate LLK over the grid
#'
#' @param prior 1D prior
#' @param x0 Values for unaffected individuals
#' @param x1 Values for affected individuals
#' @param scaled Boolean
#'
#' @export
#'
#' @examples
#'
calculate_1d_llk_grid_scaled <- function(
    prior=NULL,x0=NULL,x1=NULL,scaled=TRUE
) {

    if( is.null(prior) | is.null(x0) | is.null(x1) ) {
        msg <- "Missing input data. Check arguments. \n"
        stop(msg)
    }

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
#' @param d List of objects needed by calculate_1d_llk_grid_scaled function
#'
#' @export
#'
#' @examples
#'
#'
calculate.llk.grid.wrapper <- function(
    d
) {

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
#' @param prior Prior
#' @param llk.surf log likelihood surface
#' @param scaled boolean
#'
#' @export
#'
#' @examples
#'
calculate.integrated_1d_llk_scaled <- function(
    prior=prior,
    llk.surf=llk.surf,
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
#'
#' @return A data.frame with the following columns: max_b, summed_llk, b_ci_lhs, b_ci_rhs, POST_ACTIVE
#'
#' @export
#'
#' @examples
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

    spac <- abs(prior$b.grid[2]-prior$b.grid[1])

    out <- data.frame(
        max_b = prior$b.grid[mx[1]],
        summed_llk = log(sum(forward[[id]]$op * backward[[id]]$op)) +
            forward[[id]]$lmx + backward[[id]]$lmx,
        b_ci_lhs = rg[1] - spac/2,
        b_ci_rhs = rg[2] + spac/2,
        POST_ACTIVE = post.active)

    return(out)
}

#' Function to prepare Tree table for analysis
#'
#' @param phens input
#' @param tree input
#' @param NODE.COUNT.LIMIT input
#'
#' @return A data.frame with the tree data sorted for TreeWAS analysis
#'
#' @export
#'
#' @examples
#'

prepare_tree_table_grs <- function(
    phens=phens,
    tree=tree,
    NODE.COUNT.LIMIT=10
) {

    t <- tree
    ## remove 'unclassifiable' node
    t <- t[!t$coding %in% 99999,]

    ## Check if we need to create Top node
    top.nodes <- t[ ! t$parent_id %in% t$node_id,,drop=F]

    if( nrow(top.nodes) > 1 ) {
        ## create TOP node
        top.node <- data.frame(
            coding = -1,
            meaning = "Top node",
            node_id = 99999,
            parent_id = 0,
            selectable = 'N'
        )

        ## assing parent node to category nodes
        t[t$node_id %in% top.nodes$node_id, 'parent_id'] <- '99999'

        ## Add Top node to tree
        t <- rbind(t,top.node)
    } else {

        t <- rbind(t[ ! t$node_id %in% top.nodes$node_id,],
                   top.nodes)
    }

    map <- match(t$coding,colnames(phens))
    colcounts <- colSums(phens)
    t$counts <- colcounts[map]

    ## create dummy nodes for those intermediate nodes that are selectable
    dummy <- which( !is.na(t$counts) &
                       t$node_id %in% t$parent_id &
                           t$counts > 0)

    if( length(dummy) > 0 ) {

        for( i in 1:length(dummy) ) {
            node.id=t[dummy[i],'node_id']
            node.cod=t[dummy[i],'coding']
            node.mea=t[dummy[i],'meaning']
            node.counts=t[dummy[i],'counts']

            new.node <- data.frame(
                coding = node.cod,
                meaning = node.mea,
                node_id = paste(node.id,'_dummy',sep=''),
                parent_id = node.id,
                selectable = "Y",
                counts = node.counts)

            t[dummy[i],'coding'] <- paste(node.cod,"_int",sep='')
            t[dummy[i],'meaning'] <- paste(node.mea,"_int",sep='')
            t[dummy[i],'selectable'] <- "N"
            t[dummy[i],'counts'] <- NA

            t <- rbind(t,new.node)
        }

    }

    ## assign node ID
    t$ID <- 1:nrow(t)
    t$Par <- t[match(t$parent_id,t$node_id),'ID']
    i.ter <- t[which(!(t[,'ID'] %in% t[,'Par'])),'ID'];
    i.par <- setdiff(t[,'ID'], i.ter);

    ## place terminal nodes at the top
    t <- rbind(t[i.ter,],
                  t[i.par,])
    t[nrow(t),'Par'] <- 0

    t$ID <- 1:nrow(t)
    t$Par <- t[match(t$parent_id,t$node_id),'ID']
    t[nrow(t),'Par'] <- 0
    i.ter <- t[which(!(t[,'ID'] %in% t[,'Par'])),'ID'];
    i.par <- setdiff(t[,'ID'], i.ter);

    ## sort tree and make sure it is in ascending order
    t.ter <- t[t$ID %in% i.ter,]
    t2 <- t[nrow(t),,drop=F]
    cont <- TRUE
    while(cont) {

        child.ids <- t[t$Par %in% t2$ID,'ID']
        child.ids <- child.ids[! child.ids %in% t2$ID ]

        t2 <- rbind(t[t$ID %in% child.ids,],
                    t2
                    )

        if(nrow(t2) == nrow(t)) cont <- FALSE
    }

    t <- rbind(t.ter,
               t2[ ! t2$ID %in% t.ter$ID,])

    ## relabel nodes
    tree <- t
    tree$ID <- 1:nrow(tree)
    tree$Par <- tree[match(tree$parent_id,tree$node_id),'ID']
    tree[nrow(tree),'Par'] <- 0
    i.ter <- tree[which(!(tree[,'ID'] %in% tree[,'Par'])),'ID'];
    i.par <- setdiff(tree[,'ID'], i.ter);

    ## Filter nodes with observations below specified limit
    tree$REMOVE <- FALSE

    for( i in 1:length(i.ter) ) {
        if(is.na(tree[i.ter[i],'counts'])) {
            tree[i.ter[i],'REMOVE'] <- TRUE
        } else if(tree[i.ter[i],'counts'] <= NODE.COUNT.LIMIT) {
            tree[i.ter[i],'REMOVE'] <- TRUE
        }
    }

    for( i in 1:length(i.par) ) {
        w.d <- which(tree[,'Par'] == i.par[i])

        if( all(tree[w.d,'REMOVE']) ) {
            tree[i.par[i],'REMOVE'] <- TRUE
        }
    }

    tree <- tree[!tree$REMOVE,]
    tree.id <- 1:nrow(tree)
    tree.par <- tree.id[match(tree$Par,tree$ID)]
    tree$ID <- tree.id
    tree$Par <- tree.par
    tree[nrow(tree),'Par'] <- 0

    i.ter <- tree[which(!(tree[,'ID'] %in% tree[,'Par'])),'ID'];
    i.par <- setdiff(tree[,'ID'], i.ter);

    ## sanity checks
    NUM.NODES <- nrow(tree)
    expected.node.ids <- 1:NUM.NODES

    if(!all(tree[,'ID'] == expected.node.ids)) {
        stop("Node IDs are not compatible. Check that they are numeric and increasing bottom to top\n")
    }

    if(!all(i.ter == 1:length(i.ter))) {
        stop("terminal nodes must be first in the tree table\n")
    }

    tree <- tree[,c("ID","Par","coding","meaning")]

    return(tree)
}


