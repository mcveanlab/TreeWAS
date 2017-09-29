rm(list=ls())

# M is the number of disease
# N is the number of individual

get_count_table <- function(sparse_phe_M_by_N, ge){

    M_and_N = apply(sparse_phe_M_by_N, 2, max)
    M = M_and_N[1]
    N = max(M_and_N[2], length(ge))

    ge_idx_of_0 = which(ge==0)
    ge_idx_of_1 = which(ge==1)
    ge_idx_of_2 = which(ge==2)

    ret.tab = matrix(0, nrow = M, ncol = 6)
    for ( i in 1:M ){
        phe_idx_of_0 = 1:N
        phe_idx_of_1 = sparse_phe_M_by_N[,2][sparse_phe_M_by_N[,1] == i]
        phe_idx_of_0 = phe_idx_of_0[-phe_idx_of_1]

        ret.tab[i, ] = c(sum(ge_idx_of_0 %in% phe_idx_of_0),
                         sum(ge_idx_of_1 %in% phe_idx_of_0),
                         sum(ge_idx_of_2 %in% phe_idx_of_0),
                         sum(ge_idx_of_0 %in% phe_idx_of_1),
                         sum(ge_idx_of_1 %in% phe_idx_of_1),
                         sum(ge_idx_of_2 %in% phe_idx_of_1))
    }
    return(ret.tab)
}

N = 1000
M = 500
N_by_M = N * M

phe = matrix(round(runif(N_by_M) * 1), nrow = N)
ge = matrix(round(runif(N) * 2), nrow = N)

sparse_phe = which(phe == 1, arr.ind=T) # patient idx, phe idx
sparse_phe_M_by_N = sparse_phe[,c(2,1)] # phe idx, patient idx

a = get_count_table(sparse_phe_M_by_N, ge)
