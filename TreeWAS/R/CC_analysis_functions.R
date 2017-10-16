#' Function to prepare Tree table for analysis
#'
#' @param phens input
#' @param tree input
#' @param samples TODO
#' @param NODE.COUNT.LIMIT input
#' @param NUM.CASES.MIN TODO
#'
#' @return A data.frame with tree information sorted for TreeWAS analysis.
#'
#' @export
#'
#' @examples
#'

prepare_tree_table <- function(
    phens=phens,
    tree=tree,
    samples=samples,
    NODE.COUNT.LIMIT=10,
    NUM.CASES.MIN=10
) {

    ##################################
    ## Sanity checks on data inputs ##
    ##################################

    if(nrow(samples) != nrow(phens)) {
        stop("Sample and Pheno file are not of the same sample size\n")
    }

    if(!all(samples$ID == phens$ID)) {
        stop("Samples must have the same order in both files\n")
    }

    if( NODE.COUNT.LIMIT < NUM.CASES.MIN )
        NODE.COUNT.LIMIT <- NUM.CASES.MIN

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
        t[t$parent_id %in% 0, 'parent_id'] <- '99999'

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

    ## load genotypes to tree
    tree$GT.0.unaf <- NA
    tree$GT.1.unaf <- NA
    tree$GT.2.unaf <- NA
    tree$GT.0.aff <- NA
    tree$GT.1.aff <- NA
    tree$GT.2.aff <- NA
    tree$REMOVE <- FALSE

    tree <- tree[,c("ID","Par","coding","meaning","GT.0.unaf","GT.1.unaf",
                    "GT.2.unaf","GT.0.aff","GT.1.aff","GT.2.aff","REMOVE")]

    for( i in i.ter ) {

        node.coding <- tree[tree$ID %in% i.ter[i],'coding']

        phen.col <- which(colnames(phens) %in% node.coding)

        node.data <- array(c(
            sum(phens[[phen.col]]==0 & samples[,2]==0),
            sum(phens[[phen.col]]==0 & samples[,2]==1),
            sum(phens[[phen.col]]==0 & samples[,2]==2),
            sum(phens[[phen.col]]==1 & samples[,2]==0),
            sum(phens[[phen.col]]==1 & samples[,2]==1),
            sum(phens[[phen.col]]==1 & samples[,2]==2)
        ))

        cond1 <- sum(node.data[4:6]) < NUM.CASES.MIN
        cond2 <- sum(node.data[2:3])==0 & sum(node.data[5:6]) ==0

        if( cond1 | cond2 ) {
            tree[i,'REMOVE'] <- TRUE
        } else {
            tree[i,5:10] <- node.data
        }
    }

    for( i in i.par ) {

        w.d <- which(tree[,'Par'] == i)
        g.sums <- apply(tree[w.d,5:10],2,sum,na.rm=T)

        if(all(tree[w.d,'REMOVE'])) {
            tree[i,'REMOVE'] <- TRUE
        } else {
            tree[i,5:10] <- NA
        }
    }

    tree <- tree[!tree$REMOVE,,drop=F]

    if(nrow(tree) < 1) {
        msg <- paste("tree is empty after filtering.\n\n",sep='')
        stop(msg)
    }

    tree.id <- 1:nrow(tree)
    tree.par <- tree.id[match(tree$Par,tree$ID)]

    tree$ID = tree.id
    tree$Par <- tree.par
    tree[nrow(tree),'Par'] <- 0

    tree <- tree[,c(1:10)]

    return(tree)
}



#' Non-local prior for genetic effects. Scaled and reparametrised
#' Non-local prior of correlated effect sizes
#'
#' @param b1.range input b1 range
#' @param b2.range input b2 range
#' @param spac.b1 input b1 grid density
#' @param spac.b2 input b2 grid density
#' @param mn.b1 input b1 mean
#' @param mn.b2 input b2 mean
#' @param sd.b1 input b1 sd
#' @param sd.b2 input b2 sd
#' @param cor input correlation
#' @param k input k
#' @param pi.1 input pi_1
#' @param f.off input TODO
#' @param do.plot input TODO
#'
#' @return A list
#'
#' @export
#'
#' @examples
#'
calculate_prior_2d <- function(
    b1.range=c(-10,10),
    b2.range=c(-10,10),
    spac.b1=0.1,
    spac.b2=0.1,
    mn.b1=0,
    mn.b2=0,
    sd.b1=2,
    sd.b2=4,
    cor=0.5,
    k=0.5,
    pi.1=1e-3,
    f.off=0.1,
    do.plot=FALSE
) {
    require("mvtnorm")

    b1.grid <- seq(b1.range[1], b1.range[2], by=spac.b1);
    b2.grid <- seq(b2.range[1], b2.range[2], by=spac.b2);
    jt.prior <- array(0,c(length(b1.grid), length(b2.grid)));

    sig <- array(0,c(2,2));
    sig[1,1] <- sd.b1^2;
    sig[2,2] <- sd.b2^2;
    sig[1,2] <- cor*sd.b1*sd.b2;
    sig[2,1] <- sig[1,2];

    for (i in 1:length(b1.grid)) for (j in 1:length(b2.grid)) {
        d1 <- sqrt(b1.grid[i]^2 + b2.grid[j]^2/4)^k;
        e1 <- 1;
        if ((b1.grid[i]*b2.grid[j])<0 | abs(b2.grid[j])<abs(b1.grid[i])) e1 <- f.off;
        jt.prior[i,j] <- mvtnorm::dmvnorm(c(b1.grid[i], b2.grid[j]), c(mn.b1, mn.b2), sig)*d1*e1;
    }

    w0.1 <- which.min(abs(b1.grid));
    w0.2 <- which.min(abs(b2.grid));
    jt.prior[w0.1, w0.2] <- 0;
    s1 <- sum(jt.prior);
    jt.prior <- jt.prior*pi.1/s1;
    jt.prior[w0.1, w0.2] <- 1-pi.1;

    if (do.plot) {
        image(x=b1.grid, y=b2.grid, z=log(jt.prior));
        abline(0,2,lty="dotted");
        abline(h=0);
        abline(v=0);
    }
    return(list(b1.grid=b1.grid, b2.grid=b2.grid, jt.prior=jt.prior, val.0=c(w0.1, w0.2)));
}

#' Tree Analysis.
#' Calculate llk grid.
#' Scaled version
#'
#' @param jt.prior input
#' @param data input
#' @param do.plot input
#' @param debug input
#' @param scaled input
#'
#' @return A list.
#'
#' @export
#'
#' @examples
#'
calculate.llk.grid_scaled <- function(
    jt.prior,
    data,
    do.plot=FALSE,
    debug=FALSE,
    scaled=TRUE
) {

    to.fill <- 0;
    if (scaled) to.fill <- 1;
    op <- array(to.fill,c(length(jt.prior$b1.grid), length(jt.prior$b2.grid)));
    b0.est <- array(0, dim(op));

    if (debug) {
        return(list(op=op, b0.est=b0.est, lmx=0));
    }

    for (i in 1:length(jt.prior$b1.grid)) for (j in 1:length(jt.prior$b2.grid)) {
        tmp <- uniroot(cc_d.llk, c(-100, 100), b1=jt.prior$b1.grid[i], b2=jt.prior$b2.grid[j], aff=data[4:6], unaf=data[1:3]);
        b0.est[i,j] <- tmp$root;
        op[i,j] <- cc_llk(b0=b0.est[i,j], b1=jt.prior$b1.grid[i], b2=jt.prior$b2.grid[j], aff=data[4:6], unaf=data[1:3]);
    }

    if (do.plot) image(x=jt.prior$b1.grid, y=jt.prior$b2.grid, z=op, main="LLK", xlab="B.het", ylab="B.hom");

    mx <- max(op);
    if (scaled) op<-exp(op-mx);

    return(list(op=op, b0.est=b0.est, lmx=mx));
}


#' Tree Analysis.
#' Calculate integrated likelihood
#' Scaled version
#'
#' @param jt.prior input
#' @param llk.surf input
#' @param scaled input
#'
#' @return A numeric value.
#'
#' @export
#'
#' @examples
#'
calculate.integrated.llk_scaled <- function(
    jt.prior,
    llk.surf,
    scaled=TRUE
) {

    if (scaled) {	## Here, max of LLK.surf is 1 - this is NOT logged
        llk.sum <- sum(jt.prior$jt.prior*llk.surf$op);
        return(llk.sum);
    } else {
        mx <- max(llk.surf$op);
        op <- llk.surf$op-mx;
        op.n <- exp(op);
        llk.sum <- log(sum(op.n*jt.prior$jt.prior));
        return(llk.sum+mx);
    }
}


#' Tree Analysis.
#' Function to check B algorithm and
#' get credible sets
#'
#' @param forward TODO
#' @param backward TODO
#' @param jt.prior input
#' @param id TODO
#' @param plot TODO
#' @param return.ci TODO
#' @param verbose TODO
#' @param ci.level TODO
#' @param log.plot TODO
#'
#' @return A data.frame.
#'
#' @export
#'
#' @examples
#'
get.posterior.node <- function(
    forward,
    backward,
    jt.prior,
    id=1,
    plot=FALSE,
    return.ci=TRUE,
    verbose=FALSE,
    ci.level=0.95,
    log.plot=TRUE
) {

    tmp <- forward[[id]]$op*backward[[id]]$op;
    tmp <- tmp/sum(tmp);

    post.null <- tmp[null.id[1],null.id[2]]
    post.active <- 1 - post.null

    ## Posterior | node active
    tmp[null.id[1],null.id[2]] <- 0
    tmp <- tmp/sum(tmp);

    mx <- arrayInd(which.max(tmp), dim(tmp));
    if(verbose) cat("\nNode ", id);
    if(verbose) cat("\nMax at b1 = ", jt.prior$b1.grid[mx[1]], ", b2 = ", jt.prior$b2.grid[mx[2]]);
    if(verbose) cat("\nSummed LLK = ", log(sum(forward[[id]]$op*backward[[id]]$op))+forward[[id]]$lmx+backward[[id]]$lmx);

    ## return the C.I. conditional on node being active
    if (return.ci) {
        oo <- order(tmp, decreasing=T);
        cs <- cumsum(tmp[arrayInd(oo, dim(tmp))]);
        w.ci <- oo[c(1,which(cs<=ci.level))];
        inds <- arrayInd(w.ci, dim(tmp));
        rg.1 <- range(jt.prior$b1.grid[inds[,1]]);
        rg.2 <- range(jt.prior$b2.grid[inds[,2]]);
        if(verbose) cat("\nCI b1(", ci.level, ") = ", paste(rg.1, collapse=" - "), sep="");
        if(verbose) cat("\nCI b2(", ci.level, ") = ", paste(rg.2, collapse=" - "), sep="");
    }

    if (plot) {
        if(log.plot) tmp <- log(tmp);
        image(x=jt.prior$b1.grid, y=jt.prior$b2.grid, z=tmp, main=paste("Node", id), xlab="B1", ylab="B2");
    }

    if(verbose) cat("\n\n");

    spac1 <- abs(jt.prior$b1.grid[2]-jt.prior$b1.grid[1])
    spac2 <- abs(jt.prior$b2.grid[2]-jt.prior$b2.grid[1])

    out <- data.frame(
        max_b1=jt.prior$b1.grid[mx[1]],
        max_b2=jt.prior$b2.grid[mx[2]],
        summed_llk=log(sum(forward[[id]]$op*backward[[id]]$op))+forward[[id]]$lmx+backward[[id]]$lmx,
        b1_ci_lhs=rg.1[1] - spac1/2,
        b1_ci_rhs=rg.1[2] + spac1/2,
        b2_ci_lhs=rg.2[1] - spac2/2,
        b2_ci_rhs=rg.2[2] + spac2/2,
        POST_ACTIVE=post.active
    )

    return(out)

}


#' Function to calculate log likelihood for logistic risk and full genetic
#' model given set of parameters.
#'
#' @param b0 input
#' @param b1 input
#' @param b2 input
#' @param aff input
#' @param unaf input
#'
#' @return numeric value
#'
#' @export
#'
#' @examples
#'
cc_llk <- function(b0=b0, b1=b1, b2=b2, aff=aff, unaf=unaf) {
    return(b0*sum(aff)+b1*aff[2]+b2*aff[3]-(aff[1]+unaf[1])*log(1+exp(b0))-(aff[2]+unaf[2])*log(1+exp(b0+b1))-(aff[3]+unaf[3])*log(1+exp(b0+b2)));
}

#' Function to calculate d(log likelihood)/db0
#' for logistic risk and full genetic model given b0, b1, b2, and data
#'
#' @param b0 input
#' @param b1 input
#' @param b2 input
#' @param aff input
#' @param unaf input
#'
#' @return a nuemeric value
#'
#' @export
#'
#' @examples
#'
cc_d.llk <- function(b0=b0, b1=b1, b2=b2, aff=aff, unaf=unaf) {
    return((aff[1]+unaf[1])/(1+exp(-b0))+(aff[2]+unaf[2])/(1+exp(-(b0+b1)))+(aff[3]+unaf[3])/(1+exp(-(b0+b2)))-sum(aff));
}
