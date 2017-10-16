context("TreeWAS")

phe_filename <- system.file("extdata", "phenotype_first100.txt",
    package = "TreeWAS")
ge_filename <- system.file("extdata", "sample_file_gen_first100.txt",
    package = "TreeWAS")

phe_table <- read.table(phe_filename, header=TRUE)
phe_table$ID <- NULL
ge <- read.table(ge_filename, header=TRUE)$GE



test_that("Test phenotype input", {
    sparse_phe_M_by_N <- which(phe_table == 1, arr.ind=TRUE)[,c(2,1)]
    count_table = prepare_count_table(173, 100, sparse_phe_M_by_N, ge)

    M <- dim(phe_table)[2]
    cases <- c("0-0", "1-0", "2-0", "0-1", "1-1", "2-1")
    for ( i in 1:M ){
        truth = table(factor(paste(ge, phe_table[,i], sep="-"), levels = cases))
        expect_equal(count_table[i,], as.numeric(truth))
    }
})
