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
          --pheno_file=filename    - string. Phenotype file
          --tree_file=filename     - string. Tree file
          --outprefix=prefix       - string. Prefix for output files
          --keep=filename          - string. File with IDs to keep
          --theta=num              - numeric
          --p1=num                 - numeric
          --num.cores=num          - numeric
          --b1_max_mag=num         - numeric
          --b1_spac=num            - numeric
          --help                   - print this text
          
     
          Example:
          ./grs_tree_analysis.R \
              --sample_file=filename \
              --tree_file=filename \
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
    stop("input file is missing\n")
} else {
    sample_file <- argsL$sample_file
}

if(is.null(argsL$pheno_file)) {
    print.help()
    stop("input file is missing\n")
} else {
    pheno_file <- argsL$pheno_file
}

if(is.null(argsL$tree_file)) {
    print.help()
    stop("input file is missing\n")
} else {
    tree_file <- argsL$tree_file
}

## Check other arguments and set to default if missing

if(is.null(argsL$outprefix)) {
    outprefix <- "out"
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

if(is.null(argsL$keep)) {
    keep.file <- NULL
} else {
    keep.file <- argsL$keep
}

if(is.null(argsL$b1_max_mag)) {
    b1_max_mag <- 2
} else {
    b1_max_mag <- as.numeric(argsL$b1_max_mag)
}

if(is.null(argsL$b1_spac)) {
    b1_spac <- 0.01
} else {
    b1_spac <- as.numeric(argsL$b1_spac)
}

if(is.null(argsL$num.cores)) {
    num.cores <- 1
} else {
    num.cores <- as.numeric(argsL$num.cores)
}

cat("\nArguments:\n")
for( i in 1:length(argsL) ) {
    cat(names(argsL)[i],"\t",argsL[[i]],"\n")
}
cat("\n\n")

##################
## START SCRIPT ##
##################

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
tree <- read.table(tree_file,header=T,quote='',sep="\t")

## remove samples with missing data
samples <- samples[!is.na(samples[,2]),]

ids <- intersect(phens$ID,samples[,1])
if(!is.null(keep.file)) {
   keep.ids <- read.table(keep.file)$V1
   ids <- ids[ids %in% keep.ids]
}
   
phens <- phens[match(ids,phens$ID),]
samples <- samples[match(ids,samples[,1]),]

tree <- prepare_tree_table_grs(
    phens=phens,
    tree=tree,
    NODE.COUNT.LIMIT=10)

##################################
## Sanity checks on data inputs ##
##################################

if(nrow(samples) != nrow(phens)) {
    stop("Sample and Pheno file are not of the same sample size\n")
}

if(!all(samples[,1] == phens$ID)) {
    stop("Samples must have the same order in both files\n")
}

######################
## Do tree analysis ##
######################
cat("\nStarting TreeWAS analysis\n")

i.ter <- tree[which(!(tree[,'ID'] %in% tree[,'Par'])),'ID'];
i.par <- setdiff(tree[,'ID'], i.ter);

## Prior on effect size (1d)
b1.range <- c(-1*b1_max_mag,b1_max_mag)
prior <- calculate_1d_prior(
    pi.1=p1,
    b1.range=b1.range,
    spac.b1=b1_spac);        

## Transition probabilities
p.stay <- exp(-theta);
p.switch <- 1-p.stay;
p00 <- p.stay+p.switch*(1-p1);
logp00 <- log(p00);

debug <- FALSE;
do.backwards <- TRUE;
scaled <- TRUE;

## Find 'null model' entry in beta grid
null.id <- prior$val.0;

## Calculate the llk surfaces for the terminal nodes in the tree
cat("\nCalculating LLK surfaces\n")
llk.args <- list()
for( i in 1:length(i.ter) ) {
    code <- tree[i.ter[i],'coding']
    col.idx <- which(colnames(phens) %in% code)
    code.p <- phens[,code,with=F]
    code.p$GRS <- samples[,2]
    code.p <- as.data.frame(code.p)
    x0 <- code.p[ code.p[,1] %in% 0,2]
    x1 <- code.p[ code.p[,1] %in% 1,2]

    d <- list(
        code=code,
        x0=x0,
        x1=x1,
        prior=prior
    )

    llk.args[[i]] <- d
}

if(use.parallel) {
    llk.surf.tree <- parLapply(
        cl,
        llk.args,
        calculate.llk.grid.wrapper)
} else {
    llk.surf.tree <- list()
    for( i in 1:length(llk.args) ) {
        llk.surf.tree[[i]] <- calculate.llk.grid.wrapper(llk.args[[i]])
    }
}

## First calculate LLK under model with no active states and prior on this
llk.full.null<-rep(0, nrow(tree));

for (i in i.ter) {
    llk.full.null[i] <- log(llk.surf.tree[[i]]$op[null.id]) + llk.surf.tree[[i]]$lmx
}

for (i in i.par) {
    w.d <- which(tree[,'Par']==i);
    llk.full.null[i]<-sum(llk.full.null[w.d])+length(w.d)*logp00;
}

llk.full.null.mrca <- llk.full.null[nrow(tree)]+log(1-p1);
l.p.full.null <- log(1-p1)+(nrow(tree)-1)*logp00;

## Calculate LLK under model with transitions to-from active states
## Col 1 is LLK under node = inactive, Col 2 is integrated LLK under active, Col 3 is max LLK under non-null beta
llk.tree <- array(0, c(nrow(tree), 3));
colnames(llk.tree) <- c("LLK.0", "LLK.1", "Max.LLK.alt");         

## Start with terminal nodes
for (i in i.ter) {
    llk.tree[i,1] <- log(llk.surf.tree[[i]]$op[null.id[1]])+llk.surf.tree[[i]]$lmx;
    llk.tree[i,2] <- calculate.integrated_1d_llk_scaled(
        prior, llk.surf.tree[[i]], scaled=scaled);
    llk.tree[i,3]<-llk.surf.tree[[i]]$lmx;
}

## Things needed for backwards version
if (do.backwards) {
    g.surf.tree<-list();
}

## Now do parents
for (i in i.par) {
    
    ## Find daughters nodes
    w.d <- which(tree[,'Par']==i);
    
    tmp <- array(1, dim(llk.surf.tree[[1]]$op));
    s1 <- 0;
    
    for (j in w.d) {
        tmp.part <- (p.stay*llk.surf.tree[[j]]$op+p.switch*llk.tree[j,2]);
        tmp <- tmp*tmp.part;
        if (do.backwards) {
            gmx <- max(tmp.part);
            g.surf.tree[[j]] <- list(
                op=tmp.part/gmx, lmx=log(gmx)+llk.surf.tree[[j]]$lmx);
        }
        s1 <- s1+llk.surf.tree[[j]]$lmx;
    }
    mx <- max(tmp);
    tmp <- tmp/mx;
    
    llk.surf.tree[[i]] <- list(op=tmp, lmx=s1+log(mx));
    llk.tree[i,1] <- log(llk.surf.tree[[i]]$op[null.id[1]])+llk.surf.tree[[i]]$lmx;
    llk.tree[i,3] <- llk.surf.tree[[i]]$lmx;
    llk.tree[i,2] <- calculate.integrated_1d_llk_scaled(
        prior, llk.surf.tree[[i]], scaled=scaled);
}

mrca <- nrow(tree);
llk.full.mrca <- log(llk.tree[mrca,2])+llk.tree[mrca,3];

## Get Tree BF
tmp <- c(llk.full.mrca, llk.full.null.mrca);
mx <- max(tmp);
tmp2 <- mx+log(exp(tmp[1]-mx)-exp(tmp[2]-mx))-llk.full.null.mrca;
tmp3 <- l.p.full.null-log(1-exp(l.p.full.null));
log10_treeBF <- (tmp2 + tmp3)/log(10)

cat("\n\nLog10 BF for model with some active nodes = ", log10_treeBF, "\n\n", sep="");

## To do posterior decoding, we need backwards algorithm
## Set up LLK surface for each node and fill in terminal nodes

cat("\nNow doing posterior decoding\n")
llk.surf.tree.b <- list();
llk.tree.b <- array(0, c(nrow(tree), 2));

## Col 1 is LLK under node = inactive, Col 2 is integrated LLK under active
colnames(llk.tree.b) <- c("LLK.0", "LLK.1");            
llk.surf.tree.b[[nrow(tree)]] <- list(op=prior$prior, lmx=0);

## Now iterate down tree
for (i in (mrca-1):1) {
      
    ## Find parent of node
    np <- tree[i,'Par'];
    
    ## Create matrix ratio
    tmp <- log(llk.surf.tree.b[[np]]$op) + log(llk.surf.tree[[np]]$op) - log(g.surf.tree[[i]]$op);
    
    ## Normalise, exponentiate and sum
    mx <- max(tmp);
    tmp <- exp(tmp-mx);
    val <- p.switch*sum(tmp);
    
    b.tmp <- p.stay*tmp + prior$prior*val;
    mx2 <- max(b.tmp);
    
    llk.surf.tree.b[[i]] <- list(
        op=b.tmp/mx2, 
        lmx=mx+log(mx2)+llk.surf.tree.b[[np]]$lmx+llk.surf.tree[[np]]$lmx-g.surf.tree[[i]]$lmx);    
}       

bs.post.dec <- list()
for (i in 1:nrow(tree)) {
    tmp <- get.posterior.node_1d(
        llk.surf.tree,
        llk.surf.tree.b,
        prior, id=i, log.plot=T,plot=F, verbose=F);
    
    tmp[[4]] <- paste(tmp[[4]],collapse='-')
    tmp[[5]] <- paste(tmp[[5]],collapse='-')
    bs.post.dec[[i]] <- tmp
    
}

bs.post.dec2 <- do.call(rbind,bs.post.dec)

tmp <- bs.post.dec2
tmp <- as.data.frame(tmp)
tmp$max_b <- as.numeric(tmp$max_b)
tmp$b_ci_lhs <- as.numeric(tmp$b_ci_lhs)
tmp$b_ci_rhs <- as.numeric(tmp$b_ci_rhs)
tmp$POST_ACTIVE <- as.numeric(tmp$POST_ACTIVE)

out.table <- cbind(
    tree,
    tmp[,c(1,3:5)]
)

out.table <- out.table[order(out.table$POST_ACTIVE,decreasing=T),]

fout <- paste(outprefix,"_TreeWAS_results.txt",sep='')
write.table(out.table,file=fout,
            col.names=T,row.names=F,quote=F,sep="\t")

fout <- paste(outprefix,"_TreeWAS_results_BF.txt",sep='')
write.table(log10_treeBF,file=fout,
            col.names=F,row.names=F,quote=F,sep="\t")


