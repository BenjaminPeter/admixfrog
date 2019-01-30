require(VGAM)
LOWER=1e-9

#' throughout this sript, the following conventions are used:
#' 1. ?lomax are lomax distn function
#' 2. ?tlomax is a truncated lomax, meaning that we remove the lower tail of
#' 	the distribution. This reflects the fact that we might have low power to detect short fragments
#' 3. ?tolomax is a truncated and offset lomax, where the offset reflects a minimum age of 
#' 	fragments, here reflecting neandertal extinction
#' 4. *_opt functions are designed to return functions used in `optim`
#' 5. ?texp is a truncated exponential distribution
#' 6. for optimization, we use the log-transform of the parameters, to enforce positivity
#' 	though the _optmv insted are parameterized by mean and sd of he distn
#' 7. d*lomax are density functions, p?lomax are TAIL functions




tlomax_opt <- function(x, trunc=0){
    f <- function(par){
        a=-sum(dtlomax(x, scale=exp(par[1]), shape=exp(par[2]), trunc=trunc))
	return(a)
    }
}
tlomax_opt2 <- function(x, trunc=0){
    f <- function(par){
        a=-sum(dtlomax(x, scale=par[1], shape=par[2], trunc=trunc))
	return(a)
    }
}


#' exponential distn stuff
exp_opt <- function(x){
    f <- function(par){
        -sum(dexp(x, par[1], log=T))
    }
    return(f)
}
texp_opt <- function(x, trunc=3){
    f <- function(par){
        -sum(dtexp(x, exp(par[1]), trunc=trunc))
    }
    return(f)
}
dtexp <- function(x, rate=1, trunc=0, log=T){
    x = x[x>=trunc]
    if(log){
	dexp(x, rate, log=T) - pexp(trunc, rate, lower=F, log=T)
    } else {
	dexp(x, rate, log=F) / pexp(trunc, rate, lower=F, log=F)
    }
}

#' basic lomax stuff
#optimizing functions used in parameter estimates
lomax_opt <- function(x){
    f <- function(par){
        #-sum(dlomax(x, lambda=par[1], kappa=par[2], log=T)) #extradistr::dlomax
        -sum(VGAM::dlomax(x, scale=par[1], shape=par[2], log=T))
    }
    return(f)
}
lomax_opt2 <- function(x){
    f <- function(par){ #pars are omega, s in Kozubowski paper
        par <- pmax(0, par)
        scale = 1/ par[1]
        shape = par[2]
        -sum(VGAM::dlomax(x, scale=scale, shape=shape,log=T))
    }
    return(f)
}


#'truncated lomax stuff
dtlomax <- function(x, scale, shape, trunc, log=T){
    if(log){
	k <- VGAM::dlomax(x, scale=scale, shape=shape, log=T) - 
	VGAM::plomax(trunc, scale=scale, shape=shape, lower=F, log=T)
	k[x<trunc] <- -Inf
	return(k)
    } else {
	k <- VGAM::dlomax(x, scale=scale, shape=shape, log=F) / 
	VGAM::plomax(trunc, scale=scale, shape=shape, lower=F, log=F)
	k[x<trunc] <- NA#0
	return(k)
    }
}
tlomax_optmv <- function(x, trunc){
    f <- function(par){ #pars are omega, s in Kozubowski paper
	mean.par=par[1]
	sd.par=par[2]
	scale = sd.par * sd.par / mean.par
        shape = mean.par / scale
        -sum(dtlomax(x, scale=scale, shape=shape, trunc=trunc))
    }
    return(f)
}

#' truncatd and offset lomax stuff
#truncated and shifted
dtolomax <- function(x, scale, shape, trunc, t0, log=T){
    if(log){
	k <- dolomax(x, scale=scale, shape=shape, t0=t0, log=T) - 
	polomax(trunc, scale=scale, shape=shape, t0=t0, log=T)
	k[x<trunc] <- NA
	return(k)
    } else {
	k <- dolomax(x, scale=scale, shape=shape, t0=t0, log=F) / 
	polomax(trunc, scale=scale, shape=shape, lower=F, t0=t0, log=F)
	k[x<trunc] <- NA#0
	return(k)
    }
}
ptolomax <- function(x, scale, shape, trunc, t0, log=T){
    if(log){
	k <- polomax(x, scale=scale, shape=shape, t0=t0, log=T) - 
	polomax(trunc, scale=scale, shape=shape, t0=t0, log=T)
	k[x<trunc] <- NA#-Inf
	return(k)
    } else {
	k <- polomax(x, scale=scale, shape=shape, t0=t0, log=F) / 
	polomax(trunc, scale=scale, shape=shape, t0=t0, log=F)
	k[x<trunc] <- NA
	return(k)
    }
}
tolomax_opt <- function(x, trunc=0, t0=0){
    f <- function(par){ #pars are scale, shape
	scale = exp(par[1])
	shape = exp(par[2])
#	scale = (par[1])
#	shape = (par[2])
        res <- -sum(dtolomax(x, scale=scale, shape=shape, trunc=trunc, t0=t0, log=T))
#	print(c(res, par))
	return(res)
    }
    return(f)
}


#' offset only lomax
olomax_opt <- function(x, t0=0){
    f <- function(par){ #pars are scale, shape
#	scale = exp(par[1])
#	shape = exp(par[2])
	scale = (par[1])
	shape = (par[2])
        res <- -sum(dolomax(x, scale=scale, shape=shape, t0=t0, log=T))
	#print(c(res, par))
	return(res)
    }
    return(f)
}
#lomax distribution with an offset, to set migration after extinction to zero
#scale theta, shape k according to wiki params
dolomax <- function(x, scale, shape, t0=0, log=F){
    theta=1/scale 
    k=shape
    if(log){
	(-x * t0) - (k+1)*log( (1 + x * theta)) + log(t0+k * theta + x * t0 * theta)
    } else {
	exp(-x * t0) * (1 + x * theta)^-(k+1)*(t0+k * theta + x * t0 * theta)
    }
}
polomax <- function(x, scale, shape, t0=0, log=F){
    theta=1/scale 
    k=shape
    if(log){
	-(x * t0) -k * log(1 + x * theta)
    } else {
	exp(-x * t0) * (1 + x * theta)^-k
    }
}
ptlomax <- function(x, scale, shape, trunc=1, log=T){
    if(log){
	k <- VGAM::plomax(x, scale=scale, shape=shape, lower,F, log=T) - 
	VGAM::plomax(trunc, scale=scale, shape=shape, lower=F, log=T)
	k[x<trunc] <- NA#-Inf
	return(k)
    } else {
	k <- VGAM::plomax(x, scale=scale, shape=shape, log=F, lower=F) / 
	VGAM::plomax(trunc, scale=scale, shape=shape, lower=F, log=F)
	k[x<trunc] <- NA
	return(k)
    }
}


#' other stuff
lomax_ll <- function(x, alpha, lambda){
    n <- length(x)
    n * log(alpha) - n * log(lambda) - (1+alpha) * sum(log(1+x/lambda))
}

comparison <- function(n=10000, lambda=.3, e=1){
    a  <- rexp(n, lambda)
    print(c(length(a), length(a[a>e])))
    f0 <- exp_opt(a)   
    f1 <- lomax_opt(a)    
    f2 <- lomax_opt2(a) 
    f3 <- texp_opt(a[a>e], trunc=1) 
    f4 <- tlomax_opt(a[a>e], trunc=1) 

    o0 <- optim(c(1,1), f0, method="L-BFGS-B", lower=LOWER)
    o1 <- optim(c(1,1), f1, method="L-BFGS-B", lower=LOWER)
    o2 <- optim(c(1,1), f2, method="L-BFGS-B", lower=LOWER)
    o3 <- optim(c(1,1), f3, method="L-BFGS-B", lower=LOWER)
    o4 <- optim(c(1,1), f4, method="L-BFGS-B", lower=LOWER)

    return(list(o0, o1, o2, o3, o4))
}


null_dist <- function(){
    res <- c()
    for(i in 1:1000){
        a  <- rexp(1000, 7)
        f0 <- exp_opt(a)   
        f1 <- lomax_opt(a)    
        f2 <- lomax_opt2(a) 
        o0 <- optim(c(1,1), f0, method="L-BFGS-B", lower=LOWER)
        o2 <- optim(c(1,1), f2, method="L-BFGS-B", lower=LOWER)
        res <- c(res, o0$value - o2$value)
    }

    return(res)
}
alt_dist <- function(){
    res <- c()
    for(i in 1:1000){
        a  <- rlomax(1000, 7)
        f0 <- exp_opt(a)   
        f1 <- lomax_opt(a)    
        f2 <- lomax_opt2(a) 
        o0 <- optim(c(1,1), f0, method="L-BFGS-B", lower=LOWER)
        o2 <- optim(c(1,1), f2, method="L-BFGS-B", lower=LOWER)
        res <- c(res, o0$value - o2$value)
    }

    return(res)
}

mv2ss <- function(o){
	scale = o$par[2] * o$par[2] / o$par[1]
        shape = o$par[1] / scale
	return(c(scale=scale, shape=shape))
}

