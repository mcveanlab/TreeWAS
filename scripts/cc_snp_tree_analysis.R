#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

#######################
## Collect arguments ##
#######################
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
    args <- c("--help")
}

print.help <- function() {
    cat("
          The R Script
     
          Arguments:
          --sample_file=filename   - string. Sample file
          --tree_file=filename     - string. Tree file
          --pheno_file=filename    - string. Phenotype file
          --outprefix=prefix       - string. Prefix for output files
          --keep=filename          - string. File with IDs to keep
          --theta=num              - numeric
          --p1=num                 - numeric
          --num.cores=num          - numeric
          --b1_max_mag=num         - numeric
          --b2_max_mag=num         - numeric
          --b1_spac=num            - numeric
          --b2_spac=num            - numeric
          --help                   - print this text
          
     
          Example:
          ./cc_snp_tree_analysis.R \
              --sample_file=filename \
              --tree_file=filename \
              --pheno_file=filename \
              --outprefix=text \
              --theta=num \
              --p1=num \
              --num.cores=num

    \n\n")
    
    q(save="no")
    
}

## Help section
if("--help" %in% args) {
    print.help()
}
     
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## Check non-default arguments are parsed
if(is.null(argsL$sample_file)) {
    print.help()
    stop("sample file is missing\n")
} else {
    sample_file <- argsL$sample_file
}

if(is.null(argsL$tree_file)) {
    print.help()
    stop("tree file is missing\n")
} else {
    tree_file <- argsL$tree_file
}

if(is.null(argsL$pheno_file)) {
    print.help()
    stop("pheno file is missing\n")
} else {
    pheno_file <- argsL$pheno_file
}

## Check other arguments and set to default if missing

if(is.null(argsL$outprefix)) {
    outprefix <- 'test'
} else {
    outprefix <- argsL$outprefix
}

if(is.null(argsL$theta)) {
    theta <- 1/3
} else {
    theta <- as.numeric(argsL$theta)
}

if(is.null(argsL$p1)) {
    p1 <- 0.001
} else {
    p1 <- as.numeric(argsL$p1)
}

if(is.null(argsL$num.cores)) {
    num.cores <- 1
} else {
    num.cores <- as.numeric(argsL$num.cores)
}

if(is.null(argsL$keep)) {
    keep.file <- NULL
} else {
    keep.file <- argsL$keep
}

if(is.null(argsL$b1_max_mag)) {
    b1_max_mag <- 4
} else {
    b1_max_mag <- as.numeric(argsL$b1_max_mag)
}

if(is.null(argsL$b2_max_mag)) {
    b2_max_mag <- 4
} else {
    b2_max_mag <- as.numeric(argsL$b2_max_mag)
}

if(is.null(argsL$b1_spac)) {
    b1.spac <- 0.02
} else {
    b1.spac <- as.numeric(argsL$b1_spac)
}

if(is.null(argsL$b2_spac)) {
    b2.spac <- 0.02
} else {
    b2.spac <- as.numeric(argsL$b2_spac)
}

cat("\nArguments:\n")
for( i in 1:length(argsL) ) {
    cat(names(argsL)[i],"\t",argsL[[i]],"\n")
}
cat("\n\n")

##################
## START SCRIPT ##
##################
if( is.na(theta) | ! is.numeric(theta) ) {
    stop("the theta parameters is not numeric\n")
}

use.parallel <- ifelse(num.cores>1,TRUE,FALSE)

if(!require(TreeWAS)) {
    stop("TreeWAS R package not found.\n\n")
}
if( use.parallel ) {
    if(!require(parallel)) {
        stop("parallel R package not found.\n\n")
    }
}
if(!require(data.table)) {
    stop("data.table R package not found.\n\n")
}

#######################
## Start the cluster ##
#######################
if(use.parallel) cl <- makeCluster(num.cores)

#####################
## Read data files ##
#####################
cat("\nReading input files\n")

phens <- fread(pheno_file,header=T,sep=' ')
samples <- read.table(sample_file,header=T)
tree <- read.table(tree_file,header=T,sep="\t",quote='')

ids <- intersect(phens$ID,samples$ID)
if(!is.null(keep.file)) {
   keep.ids <- read.table(keep.file)$V1
   ids <- ids[ids %in% keep.ids]
}
   
phens <- phens[match(ids,phens$ID),]
samples <- samples[match(ids,samples$ID),]

tree <- prepare_tree_table(
    phens=phens,tree=tree,samples=samples,NODE.COUNT.LIMIT=10,NUM.CASES.MIN=10)

########################
## Tree Sanity Checks ##
########################
NUM.NODES <- nrow(tree)
expected.node.ids <- 1:NUM.NODES
if(!all(tree[,'ID'] == expected.node.ids)) {
    stop("Node IDs are not compatible. Check that they are numeric and increasing bottom to top\n")
}

i.ter <- tree[which(!(tree[,'ID'] %in% tree[,'Par'])),'ID'];
i.par <- setdiff(tree[,'ID'], i.ter);
if(!all(i.ter == 1:length(i.ter))) {
    stop("terminal nodes must be first in the tree table\n")
}

## Check if there are empty terminal nodes
aff <- apply(tree[i.ter,8:10],1,sum)
unaf <- apply(tree[i.ter,5:7],1,sum)

if( any(aff==0) | any(unaf==0) ) {
    stop("Remove terminal nodes with 0 affected/unaffected counts\n")
}

######################
## Do tree analysis ##
######################
cat("\nStarting TreeWAS analysis\n")

b1.range = c(-b1_max_mag,b1_max_mag)
b2.range = c(-b2_max_mag,b2_max_mag)

jt.prior <- calculate_prior_2d(
    pi.1     = p1,
    spac.b1  = b1.spac,
    spac.b2  = b2.spac,
    b1.range = b1.range,
    b2.range = b2.range
);

## Transitions
p.stay <- exp(-theta);
p.switch <- 1-p.stay;
p00 <- p.stay+p.switch*(1-p1);
logp00 <- log(p00);
                                                        
debug <- FALSE;
do.backwards <- TRUE;
scaled <- TRUE;

#Find 'null model' entry in beta grid
null.id <- jt.prior$val.0;

#Set up LLK surface for each node and fill in terminal nodes
calculate.llk.grid.wrapper <- function(i,tree,i.ter,jt.prior,debug,scaled) {
    library(TreeWAS)
    data <- as.matrix(tree[i.ter[i],5:10])
    surf <- calculate.llk.grid_scaled(
        jt.prior,data ,
        debug=debug, scaled=scaled);
    return(surf)
}

if(use.parallel) {    
    llk.surf.tree <- parLapply(
        cl,
        1:nrow(tree[i.ter,]),
        calculate.llk.grid.wrapper,
        tree,i.ter,jt.prior,debug,scaled)
} else {
    llk.surf.tree <- list()
    for( i in 1:length(i.ter) ) {
        llk.surf.tree[[i]] <- calculate.llk.grid.wrapper(
            i,tree,i.ter,jt.prior,debug,scaled)
    }
}

## First calculate LLK under model with no active states and prior on this
llk.full.null<-rep(0, nrow(tree));
for (i in i.ter) {
    llk.full.null[i] <- log(llk.surf.tree[[i]]$op[null.id[1], null.id[2]])+llk.surf.tree[[i]]$lmx;
}

for (i in i.par) {
    w.d <- which(tree[,2]==i);
    llk.full.null[i] <- sum(llk.full.null[w.d]) + length(w.d) * logp00;
}

llk.full.null.mrca <- llk.full.null[nrow(tree)] + log(1-p1);
l.p.full.null <- log(1-p1) + (nrow(tree)-1) * logp00;

llk.tree <- array(0, c(nrow(tree), 3));
colnames(llk.tree) <- c("LLK.0", "LLK.1", "Max.LLK.alt");         

## Start with terminal nodes
for (i in i.ter) {
    llk.tree[i,1] <- log(llk.surf.tree[[i]]$op[null.id[1], null.id[2]]) + llk.surf.tree[[i]]$lmx;
    llk.tree[i,2] <- calculate.integrated.llk_scaled(jt.prior, llk.surf.tree[[i]], scaled=scaled);
    llk.tree[i,3] <- llk.surf.tree[[i]]$lmx;
}

## Things needed for backwards version
if (do.backwards) {  
        g.surf.tree <- list();
}

## Now do parents
for (i in i.par) {
    ## Find daughters
    w.d <- which(tree[,2]==i);
    
    tmp <- array(1, dim(llk.surf.tree[[1]]$op));
    s1 <- 0;

    for (j in w.d) {
        tmp.part <- (p.stay * llk.surf.tree[[j]]$op + p.switch * llk.tree[j,2]);
        tmp <- tmp * tmp.part;
        if (do.backwards) {
            gmx <- max(tmp.part);
            g.surf.tree[[j]] <- list(op=tmp.part/gmx, lmx=log(gmx) + llk.surf.tree[[j]]$lmx);
        }
        s1 <- s1 + llk.surf.tree[[j]]$lmx;
    }
    mx <- max(tmp);
    tmp <- tmp/mx;
        
    llk.surf.tree[[i]] <- list(op=tmp, lmx=s1+log(mx));
    llk.tree[i,1] <- log(llk.surf.tree[[i]]$op[null.id[1], null.id[2]]) + llk.surf.tree[[i]]$lmx;
    llk.tree[i,3] <- llk.surf.tree[[i]]$lmx;
    llk.tree[i,2] <- calculate.integrated.llk_scaled(jt.prior, llk.surf.tree[[i]], scaled=scaled);
}

mrca <- nrow(tree);
llk.full.mrca <- log(llk.tree[mrca,2])+llk.tree[mrca,3];


## Get BF
tmp <- c(llk.full.mrca, llk.full.null.mrca);
mx <- max(tmp);
tmp2 <- mx+log(exp(tmp[1]-mx)-exp(tmp[2]-mx))-llk.full.null.mrca;
tmp3 <- l.p.full.null-log(1-exp(l.p.full.null));
log10_treeBF <- (tmp2 + tmp3)/log(10)

cat("\n\nlog10(TreeBF) ",log10_treeBF,"\n\n",sep='')

llk.surf.tree.b <- list();
llk.tree.b <- array(0, c(nrow(tree), 2));
## Col 1 is LLK under node = inactive, Col 2 is integrated LLK under active
colnames(llk.tree.b) <- c("LLK.0", "LLK.1");            

llk.surf.tree.b[[nrow(tree)]] <- list(op=jt.prior$jt.prior, lmx=0);

## Now iterate down tree
for (i in (mrca-1):1) {
    
    ## Find parent of node
    np <- tree[i,2];

    ## Create matrix ratio
    tmp <- log(llk.surf.tree.b[[np]]$op) + log(llk.surf.tree[[np]]$op) - log(g.surf.tree[[i]]$op);

    ## Normalise, exponentiate and sum
    mx <- max(tmp);
    tmp <- exp(tmp-mx);
    val <- p.switch*sum(tmp);
    
    b.tmp <- p.stay*tmp+jt.prior$jt.prior*val;
    mx2 <- max(b.tmp);

    llk.surf.tree.b[[i]] <- list(
        op=b.tmp/mx2, 
        lmx=mx+log(mx2)+llk.surf.tree.b[[np]]$lmx+llk.surf.tree[[np]]$lmx-g.surf.tree[[i]]$lmx);
    
}       

bs.post.dec <- list()
for (i in 1:nrow(tree)) {
    tmp <- get.posterior.node(
        llk.surf.tree,
        llk.surf.tree.b,
        jt.prior, id=i, log.plot=T,plot=F, verbose=F);

    tmp[[4]] <- paste(tmp[[4]],collapse='-')
    tmp[[5]] <- paste(tmp[[5]],collapse='-')
    bs.post.dec[[i]] <- tmp
    
}

bs.post.dec2 <- do.call(rbind,bs.post.dec)

tmp <- bs.post.dec2
tmp <- as.data.frame(tmp)

out.table <- cbind(
    tree[,c(1:10)],
    tmp[,c(1:2,4:8)]
)

out.table <- out.table[order(out.table$POST_ACTIVE,decreasing=T),]

fout <- paste(outprefix,"_TreeWAS_results.txt",sep='')
write.table(out.table,file=fout,
            col.names=T,row.names=F,quote=F,sep="\t")

fout <- paste(outprefix,"_TreeWAS_results_BF.txt",sep='')
write.table(log10_treeBF,file=fout,
            col.names=F,row.names=F,quote=F,sep="\t")


