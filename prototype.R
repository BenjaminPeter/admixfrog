#test for graphical model

#Let: 
# - O be the observed read 0 = ref, 1 = der
# - X be the molecule; 0 = ref, 1 = der
# - G be the genotype; 0 = ref, 1 = der , haploid case
# - C the contamination state 0 = endogenous, 1 = contaminant
# - tau = P(G), to be estimated
# - c = P(C), to be estimated
# - psi = contaminant allele frequency


# 1. generate fake data
fake_data <- function(){
tau1 = 0.7
tau2 = 0.1
c1 = 0.15
c2 = 0.45
e = 0.01
b = 0.01


tau = c(tau1, tau2)
cont = c(c1, c2)

n11 = 1000 # n observations for Z1 L1
n12 = 4500 # n observations for Z1 L1
n21 = 2000 # n observations for Z2 L1
n22 = 2500 # n observations for Z2 L2
n = n11 + n12 + n21 + n22

Z = c(rep(1, n11+n12), rep(2, n21 + n22))
R = c(rep(1, n11), rep(2, n12), rep(1, n21),rep(2, n22))

C <- rep(NA, n)
C[R==1] = rbinom(sum(R==1), 1, c1)
C[R==2] = rbinom(sum(R==2), 1, c2)

G <- rep(NA, n)
G[Z==1] = rbinom(sum(Z==1), 2, tau1)
G[Z==2] = rbinom(sum(Z==2), 2, tau2)

psi <- runif(n, .1, .7)
A <- rbinom(n, 1, psi)
X <- (1-C) * rbinom(n, 1, G/2) + C * A

errors = as.logical(rbinom(n11+n12+n21+n22, 1, e))
O=X;O[errors] = 1-X[errors]
O[O==1] = rbinom(sum(O==1), 1, 1 - b)
return(list(O=O, psi=psi, C=C, G=G, Z=Z, R=R, tau=tau, cont=cont, e=e, b=b))
}

D = fake_data()



#' manual fwd algorithm
#' G[l, j] = Pr(G_l = j | Z_l=k, tau_k, F_k)
#' j = [0, 1, 2]
fwd_p_g <- function(Z, tau){
    pg = matrix((dbinom(0:2, 2, rep(tau[Z], each=3))), ncol=3, byrow=T)
}

#' manual fwd algorithm
#' C[lry, j] = Pr(C_lry = j | c_r)
#' j = [0, 1]
fwd_p_c <- function(R, cont){
    # table [l x 2] [l x 2] of Pr(C_l = 0, 1]
    cbind(1-cont[R], (cont[R])) #col 1 = endo, 2 = cont
}

fwd_p_a <- function(psi){
    cbind(1-psi, psi)
}

#' manual fwd algorithm
#' X[lry, j] = Pr(X_lry = j | c_r, tau_k, F_k, Z_l, psi_l)
#'           = \sum[G_l] \sum_[C_r]  [ Pr(X_lry = j | C_r, g_l, psi_l)  x 
#'                  Pr(C_r | c_r) Pr(G_l | Z_l=k, tau_k, F_k) ] 
#' l: number of sites
#' r: number of read groups
#' y: number of reads in read groups
#' j = [0, 1]
fwd_p_x <- function(pg, pc, pa){
    # table [l x 2] of Pr(X_l = 0, 1)
    x_is_1 = pc[,2] * pa[,2] #is contaminant, and contaminant is alt
    x_is_1 = x_is_1 + pc[,1] * pg[,2] / 2 #is endo, and het, and alt from het
    x_is_1 = x_is_1 + pc[,1] * pg[,3]  #is endo, and homo alt
    return(cbind(1-x_is_1, x_is_1))
}

fwd_p_gc <- function(pg, pc){
    res <- array(NA, c(nrow(pg), 3, 2))
    res[,,1] <- pg * pc[,1]
    res[,,2] <- pg * pc[,2]
    #res[,1,] <- res[,1,] * pc[,1]
    #res[,2,] <- res[,2,] * pc[,2]
    return(res)
}

#' manual fwd algorithm
#' O[lry] = Pr(O_lry | c_r, tau_k, F_k, Z_l, e, b)
#'           = \sum[X_lry] [Pr(O[lry | X_lry, e, b) Pr(X_lry, c_r, tau_k, F_k, Z_l = k)
#' l: number of sites
#' r: number of read groups
#' y: number of reads in read groups

fwd_p_o <- function(o, px, e, b){
    # table [l x 2] of Pr(O_l = 0, 1)
    po = rep(NA, length(o))
    po[o==1] = px[o==1,1] * e     + px[o==1,2] * (1-b)
    po[o==0] = px[o==0,1] * (1-e) + px[o==0,2] * (b)
    return(po)
}

fwd_algorithm = function(
                         obs, #O[lry] Observations
                         Z, # Z[l] SFS category for  SNP[l]
                         R, # R[lr] RG category for a read group
                         tau, # tau[k], conditional derived SFS
                         conts, # conts[r], contamination rates
                         psi,   # psi[l], DAF for SNP l
                         e, b   #error, bias rate)
                         ){
    
    
    fpg = fwd_p_g(Z, tau)
    fpc = fwd_p_c(R, conts)
    fpa = fwd_p_a(psi)
    fpgc = fwd_p_gc(fpg, fpc)
    fpx = fwd_p_x(fpg, fpc, fpa)
    fpo = fwd_p_o(obs, fpx, e, b)
    return(list(G=fpg, C=fpc, GC=fpgc, X=fpx, O=fpo, A=fpa))
}



#' manual bwd algorithm
#' O[lry, j] = Pr(O_lry | X_lry = j, e, b)
bwd_p_o_given_x <- function(O, e, b){
    res = matrix(NA, ncol=2, nrow=length(O))
    res[O==0, 1] = 1 - e #X=0, O=0
    res[O==1, 1] = e #X=0, O=1
    res[O==0, 2] = b #X=1, O=0
    res[O==1, 2] = 1 - b #X=1, O=1
    return(res)
}

#' manual bwd algorithm
#' V[lry, g, c] = Pr(O_lry | G_l=g, c_r = c, psi_l)
#'              = \sum_x Pr(O_lry | X_lry=x) Pr(X_lry=x | G_l=g, c_r=c, psi_l)
bwd_p_o_given_gca <- function(bpx){
    res = matrix(NA, ncol=5, nrow=nrow(bpx)) #5cols: C=1, G=0.2; C=0, A=0,1
    res[,1] = bpx[,1] * 1 + bpx[,2] * 0 #Pr(O | C=0, G=0)
    res[,2] = bpx[,1] * .5 + bpx[,2] * 0.5 #Pr(O | C=0, G=1)
    res[,3] = bpx[,1] * 0 + bpx[,2] * 1 #Pr(O | C=0, G=2)
    res[,4] = bpx[,1] * 1 + bpx[,2] * 0 #Pr(O | C=1, A=0)
    res[,5] = bpx[,1] * 0 + bpx[,2] * 1 #Pr(O | C=1, A=1)
    return(res)
}

#' manual bwd algorithm
#' V[lry, g] = Pr(O_lry | G_l=g, e, b)
#'              = \sum_c Pr(O_lry | X_lry=x, e, b) Pr(X_lry=x | G_l=g, c_r=c)
bwd_p_o_given_g <- function(bpgca, cont, R, psi){
    rowSums(bpgca[,4:5] * cont[R] * fwd_p_a(psi)) + (1-cont[R]) * bpgca[,1:3]
}

#' manual bwd algorithm
#' V[lry, g] = Pr(O_lry | c_r=c, e, b)
#'              = \sum_g Pr(O_lry | X_lry=x, e, b) Pr(X_lry=x | G_l=g, c_r=c)
bwd_p_o_given_c <- function(bpgca, Z, tau, psi){
    c0 = rowSums(bpgca[,1:3] * fwd_p_g(Z, tau))
    c1 = rowSums(bpgca[,4:5] * fwd_p_a(psi))
    cbind(c0, c1)

}

bwd_p_o_given_a <- function(bpgca, cont, R, tau, Z){
    rowSums(bpgca[,1:3] * (1-cont[R]) * fwd_p_g(Z, tau)) +
        cont[R] * bpgca[,4:5]
}

#' manual bwd algorithm
#' V[lry, g, c] = Pr(O_lry | Z_l = k, c_r = c)
#'              = \sum_g Pr(O_lry | G_l=g, c_r =c) Pr(G_l=g | Z_l =k, tau, F)
bwd_p_o_given_z <- function(bpg, tau, Z){
    res  = dbinom(0, 2, tau[Z]) * bpg[,1] +  # c=0=no_cont
      dbinom(1, 2, tau[Z]) * bpg[,2] +  # c=0=no_cont
      dbinom(2, 2, tau[Z]) * bpg[,3]  # c=0=no_cont
    return(res)
}

#' manual bwd algorithm
#' V[lry, g, c] = Pr(O_lry | psi)
#'              = \sum_g Pr(O_lry | G_l=g, c_r =c) Pr(G_l=g | Z_l =k, tau, F)
bwd_p_o_given_z <- function(bpg, tau, Z){
    res  = dbinom(0, 2, tau[Z]) * bpg[,1] +  # c=0=no_cont
      dbinom(1, 2, tau[Z]) * bpg[,2] +  # c=0=no_cont
      dbinom(2, 2, tau[Z]) * bpg[,3]  # c=0=no_cont
    return(res)
}

bwd_p_o_given_psi <- function(bpa, psi){
    bpa[,2] * psi + bpa[,1] * (1-psi)
}

bwd_algorithm = function(
                         obs, #O[lry] Observations
                         Z, # Z[l] SFS category for  SNP[l]
                         R, # R[lr] RG category for a read group
                         tau, # tau[k], conditional derived SFS
                         conts, # conts[r], contamination rates
                         psi,   # psi[l], DAF for SNP l
                         e, b   #error, bias rate)
                         ){
    
    bpx = bwd_p_o_given_x(obs, e, b)
    bpgca = bwd_p_o_given_gca(bpx)
    bpa = bwd_p_o_given_a(bpgca,  conts, R, tau, Z)
    bpg = bwd_p_o_given_g(bpgca, conts, R, psi)
    bpc = bwd_p_o_given_c(bpgca, Z, tau, psi)
    bpz = bwd_p_o_given_z(bpg, tau, Z)
    bppsi = bwd_p_o_given_psi(bpa, psi)
    return(list(GCA=bpgca, C=bpc, G=bpg, X=bpx, Z=bpz, A=bpa, psi=bppsi))
}

FWD = fwd_algorithm(obs=D$O, Z=D$Z, R=D$R, tau=D$tau, conts=D$cont, psi=D$psi, e=0.01, b=0.01)
BWD = bwd_algorithm(obs=D$O, Z=D$Z, R=D$R, tau=D$tau, conts=D$cont, psi=D$psi, e=0.01, b=0.01)

posteriors <- function(FWD, BWD){
    post_x = FWD$X * BWD$X #Pr(X | Z) Pr(O|X) 
    post_x = post_x / rowSums(post_x)

    post_g = FWD$G * BWD$G #Pr(G | Z) Pr(O |G)
    post_g = post_g / rowSums(post_g)

    post_c = FWD$C * BWD$C
    post_c = post_c / rowSums(post_c)

    post_a = FWD$A * BWD$A
    post_a = post_a / rowSums(post_a)


    return(list(G=post_g, X=post_x, C=post_c, A=post_a))
    
}

all_ll <- function(FWD, BWD){
    return(c(
                Z = sum(log(BWD$Z)),
                psi = sum(log(BWD$psi)),
                G = sum(log(rowSums(BWD$G * FWD$G))),
                C = sum(log(rowSums(BWD$C * FWD$C))),
                A = sum(log(rowSums(BWD$A * FWD$A))),
                X = sum(log(rowSums(BWD$X * FWD$X))),
                O = sum(log(FWD$O))))

}

POST = posteriors(FWD, BWD)

em <- function(obs, Z, R, tau0, cont0, psi, e0=0.01, b0=0.01){
    tau <- tau0
    cont <- cont0
    e = e0
    b = b0
    prev_ll <- -Inf
    res <- matrix(ncol=8, nrow=0)
    colnames(res) <- c("iter", "tau0", "tau1", "cont0", "cont1", "e", "b", "ll")
    res = res %>% as_tibble
    row = c(0, tau, cont, e, b, NA)
    names(row) <- names(res)

    for(i in 1:1000){
        FWD = fwd_algorithm(obs=obs, Z=Z, R=R, tau=tau, conts=cont, psi=psi, e=e, b=b)
        BWD = bwd_algorithm(obs=obs, Z=Z, R=R, tau=tau, conts=cont, psi=psi, e=e, b=b)
        POST = posteriors(FWD, BWD)
        tau = sapply(1:2, function(i){
                         mean(POST$G[Z==i,3] + POST$G[Z==i, 2] / 2)
                         })
        cont = sapply(1:2, function(i)mean(POST$C[R==i, 2]))
        #e = sum(POST$X[D$O==1,1]) / sum(POST$X[,1])
        #b = sum(POST$X[D$O==0,2]) / sum(POST$X[,2])
        ll = all_ll(FWD, BWD)
        if(i %% 5 == 0){
        row = c(i, tau, cont, e, b, ll[1])
        names(row) <- names(res)
        res <- bind_rows(res, row)

        print(row)
        print(ll)
        print(c(i, unname(ll[1] - prev_ll)))
        print("---")
        }
        if(ll[1] - prev_ll  < 0.001) return(res)
        prev_ll = ll[1]
    }


    return(res)


}

plot_em <- function(em, tau0, cont0){
    df = em %>% pivot_longer(c(tau0:cont1, ll)) %>% mutate(panel=substr(name, 1, 3)) 
    true_df <- data.frame(panel=c("tau", "tau", "con", "con", "ll"), 
                          name=c("tau0", "tau1", "cont0", "cont1", "ll"), 
                          value=c(tau0, cont0, NA))
    P = df %>% ggplot(aes(x=iter, y=value, color=name)) + 
        facet_grid(panel~., scale='free_y') + geom_line() + geom_point()
    P = P + geom_hline(data=true_df, aes(yintercept=value, color=name), lty=2)
}


