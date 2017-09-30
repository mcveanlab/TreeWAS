#' Preparing count table from phenotype table
#'
#' Let M be the number of disease, and N be the number of individual. We can
#' construct a phenotype indicator matrix to show which person carries which
#' disease
#' @param M number of disease
#' @param N number of individual
#' @param sparse_phe_M_by_N is a two column matrix, to summaries which disease
#' is carried by which patient. The first column indicate the index of the
#' diease, and the second column indicates patient index.
#' @param ge is a indicator vector length N, to identify if the nth patient
#' carries an alternative allele at the target position.
#'
#' @return An M by 6 matrix. Each column show the number of counts of genotype
#' and phenotype pairs of "0-0", "1-0", "2-0", "0-1", "1-1", "2-1"
#'
#' @export
#'
#' @examples
#' N = 1000
#' M = 500
#' N_by_M = N * M
#' phe = matrix(round(round(rexp(N_by_M)/5) > 0), nrow = N)
#' ge = matrix(round(runif(N) * 2), nrow = N)
#' sparse_phe = which(phe == 1, arr.ind=TRUE) # patient idx, phe idx
#' sparse_phe_M_by_N = sparse_phe[,c(2,1)] # phe idx, patient idx
#' count_table = prepare_count_table(M, N, sparse_phe_M_by_N, ge)
#'

prepare_count_table <- function(M, N, sparse_phe_M_by_N, ge){

    ge_idx_of_0 = which(ge==0)
    ge_idx_of_1 = which(ge==1)
    ge_idx_of_2 = which(ge==2)

    ret_tab = matrix(0, nrow = M, ncol = 6)
    for ( i in 1:M ){
        phe_idx_of_0 = 1:N
        phe_idx_of_1 = sparse_phe_M_by_N[,2][sparse_phe_M_by_N[,1] == i]
        if (length(phe_idx_of_1) > 0){
            phe_idx_of_0 = phe_idx_of_0[-phe_idx_of_1]
        }
        ret_tab[i, ] = c(sum(ge_idx_of_0 %in% phe_idx_of_0),
                         sum(ge_idx_of_1 %in% phe_idx_of_0),
                         sum(ge_idx_of_2 %in% phe_idx_of_0),
                         sum(ge_idx_of_0 %in% phe_idx_of_1),
                         sum(ge_idx_of_1 %in% phe_idx_of_1),
                         sum(ge_idx_of_2 %in% phe_idx_of_1))
    }
    return(ret_tab)
}
