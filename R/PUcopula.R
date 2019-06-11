# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' stormflood
#'
#' This is the data set used in Cottin, Pfeifer (2014).
#'     The table contains some original data from an insurance portfolio of storm and flooding losses, observed over a period of 20 years.
#'
#' @name stormflood
#' @docType data
#' @references Cottin, Claudia, and Dietmar Pfeifer. "From Bernstein polynomials to Bernstein copulas." J. Appl. Funct. Anal 9.3-4 (2014): 277-288.
#' @keywords storm and flood losses
NULL

# cf: #http://sbfnk.github.io/mfiidd/mcmc.html
# - target: the target distribution, a function that takes one
#   argument (a number) and returns the (logged) value of a
#   distribution
# - init.theta: the initial value of theta, a number
# - proposal.sd: the standard deviation of (Gaussian) proposal
#   distribution
# - n.iterations: the number of iterations
# The function returns a vector of samples of theta from the target
# distribution
.mcmcMH  <- function(target, init.theta = 0.6, proposal.sd=1, n.iterations=1) {
  # evaluate "target" at "init.theta"
  target.theta.current <- target(init.theta)
  # current value of theta, vector of samples, number of accepted runs
  theta.current <- init.theta
  samples <- theta.current
  accepted <- 0
  # run MCMC for n.iteration interations
  for (i.iteration in seq_len(n.iterations)) {
    # draw new theta from the (Gaussian) proposal distribution
    # Note that this step is vectorized for any arbitratry theta
    # (useful when sampling from a multivariate target distribution)
    theta.proposed <- rnorm(n = length(theta.current),
                            mean = theta.current,
                            sd = proposal.sd)
    # the functions of 'fitmodel' need a named parameter vector:
    names(theta.proposed) <- names(theta.current)
    # evaluate the function target at the proposed theta
    print(target)
    target.theta.proposed <- target(theta.proposed)
    # compute Metropolis-Hastings ratio (acceptance probability). Since
    # the multivariate Gaussian is symmetric no need to consider proposal distribution here
    log.acceptance <- target.theta.proposed - target.theta.current
    print(paste("theta.proposed",theta.proposed))
    print(paste("target.theta.proposed",target.theta.proposed))
    print(paste("target.theta.current",target.theta.current))
    # draw random number number between 0 and 1 using
    r <- runif(1)
    # test acceptance by comparing the random number to the Metropolis-Hastings ratio (acceptance probability)
    # "exp" because we calculated the logarithm of the Metropolis-Hastings ratio before
    if (r < exp(log.acceptance)) { # if accepted:
      # change the current value of theta to the proposed theta
      theta.current <- theta.proposed
      # updated the current value of the target
      target.theta.current <- target.theta.proposed
      # update number of accepted proposals
      accepted <- accepted + 1
    }
    # add current theta to the vector of samples
    # `rbind` in order to deal with multivariate target.
    samples <- rbind(samples, theta.current, deparse.level=0)
    # print current state of chain and acceptance rate
    # use paste() to deal with the case where `theta` is a vector
    message("iteration: ", i.iteration, ", chain:", paste(theta.current, collapse=" "),
            ", acceptance rate:", accepted / i.iteration)
  }
  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
}

.mh_sample <- function(dens, start = 0.5, nreps = 1000, prop_sd = 0.2, ...){ # start 0 oder 0.5?, 1 oder 0.2
  # https://www.nicksolomon.me/post/learn-metropolis-hastings-sampling-with-r/
  theta <- numeric(nreps)
  theta[1] <- start

  for (i in 2:nreps){
    theta_star <- msm::rtnorm(1, mean=theta[i - 1], sd=prop_sd, lower=0, upper=1) #rnorm(1, mean = theta[i - 1], sd = prop_sd)
    # print(paste("theta[i-1]",theta[i-1]))
    # print(paste("theta_star",theta_star))
    # print(paste("dens(theta_star, ...) ",dens(theta_star, ...) ))
    # print(paste("dens(theta[i - 1], ...) ",dens(theta[i - 1], ...) ))
    alpha = dens(theta_star, ...) / dens(theta[i - 1], ...)

    if (runif(1) < alpha) theta[i] <- theta_star
    else theta[i] <- theta[i - 1]
  }

  return(tail(theta,1))
}

#' An S4 class to represent a partition of unity copula.
#'
#' @slot dim A length-one numeric vector; dimension of PU-copula
#' @slot family Character or vector of characters: family determining the desnities fdens
#' @slot par.factor A length-n numeric vector
#' @slot pars.a A length-one numeric vector: defined as integral from 0 to 1 over phi(s,u) for s in r. Will be evauated numerically
#' @slot p A length-one numeric vector
#' @slot phis A list of functions \eqn{\varphi_k(s,u)}{phi_k(s,u)} for \eqn{k=1,dots,d\in N}{k=1,...,d in N}.    Either continuous case: A list of functions which represent Lebesgue densities of distributions over R with a parameter u in (0,1), i.e.
#'       \deqn{\varphi_k(s,u)\geq 0 \text{ and } \int_{-\infty}^\infty \varphi_k(s,u)\, ds = 1 \text{ for } u \in (0,1),}{phi_k(s,u) \geq 0 and \int phi_k(s,u), ds = 1 for u in (0,1),}
#'       Or (discrete case): A list of functions which represent discrete probabilities over Z+ with a parameter u in (0,1).
#' @slot alphs A list of functions. Integral over the elements of phis w.r.t. u, i.e. \deqn{\alpha_k(s):=\int_0^1 \varphi_k(s,u)\, du \in (0,\infty)}{alpha_k(s):=int_0^1 phi_k(s,u) du in (0,\infty)}
#' @slot alphsCDF A list of functions. CDFs corresponding to the densities defined by alphs, i.e. \deqn{A_k(s) := \int_{-\infty}^s \alpha_k(w) \, dw \text{ for } s \in R}{A_k(s) := int_-\infty^s alpha_k(w) dw for s in R}
#' @slot fdenss A list of functions. Densities obtained from normalizing the functions \eqn{\varphi(s,u)}{phi(s,u)} (from slot phis) w.r.t. \eqn{u \in 0,1}{u in (0,1)}, i.e. \deqn{f_k(s,u) := \frac{\varphi_k(s,u)}{\alpha_k(s)}, \, u \in (0,1) \text{ for } s \in R}{f_k(s,u) := phi_k(s,u)/alpha_k(s), u in (0,1) for s in R}
#' @slot PUdens density function of the (continuous) partition of unity copula defined by  \deqn{c(\mathbf{u}) := \int_{-\infty}^\infty \cdots \int_{-\infty}^\infty p(s_1,\cdots,s_d) \prod_{k=1}^d f_k(s_k, u_k) \, ds_1 \cdots ds_d, \mathbf{u} = (u_1,\cdots,u_d) \in (0,1)^d }{c(u) := int ... int p(s_1,...,s_d) \prod_k=1^d f_k(s_k, u_k) ds_1 ... ds_d, u = (u_1,...,u_d) in (0,1)^d }
#'       where \eqn{p(s_1,\cdots,s_d)}{p(s_1,...,s_d)} denotes the density of an arbitrary \eqn{d}-dimensional random vector \eqn{\mathbf{S}=S_1,\cdots,S_d)}{S=S_1,...,S_d)} over \eqn{R^d} with marginal densities \eqn{\alpha(cdot)}{alpha_k(.)} for \eqn{S_k}.
#' @slot fq A length-one numeric vector
#' @slot gq A length-one numeric vector
#' @slot patch A length-one numeric vector
#' @slot data A length-one numeric vector
#' @slot ranks A length-one numeric vector
#' @slot relRanks A length-one numeric vector
#' @slot rpatch A function which generates random samples from the specified patchwork.
#' @slot patchpar Parameter for the chosen patchwork type (used currently for type Gauss)
#' @slot rand random number generator to simulate from PU copula
#'
#' @return An object of class PUcopula
#'
#' @examples
#' # Use dataset stormflood
#' data(stormflood)
#'
#' # example: choose Gamma copula model
#' x <- PUCopula(family="gamma", pars.a=c(40, 43), patch="lFrechet",data=stormflood)
#'
#' # The following plots show the common ranks and several copula drivers;
#' #   these plots do not depend on the above chosen family of the PUcopula.
#' # 2*3 plots in one
#' par(mfrow=c(2,3))
#' # plot the empirical rank vectors
#' plot(x@ranks, sub="emp. rank vectors", xlab="", ylab="")
#' # plot the lower Fréchet driver, i.e. phi=-1
#' plot(x@rpatch(2000,"lFrechet"), sub="lower Fréchet driver, phi=-1", xlab="", ylab="")
#' # plot for phi=-0.8
#' plot(x@rpatch(2000,"Gauss",-0.8), sub="phi = -0.8", xlab="", ylab="")
#' # plot for phi=0 (rook copula)
#' plot(x@rpatch(2000,"rook"), sub="rook copula, phi=0", xlab="", ylab="")
#' # plot for phi=0.9
#' plot(x@rpatch(2000,"Gauss",0.9), sub="phi=0.9", xlab="", ylab="")
#' # plot for phi=1 (upper Fréchet driver)
#' plot(x@rpatch(2000,"uFrechet"), sub="upper Fréchet driver, phi=1", xlab="", ylab="")
#' # single plot window
#' par(mfrow=c(1,1))
#'
#' # plot the densities given by phis (=densities of Gamma distributions)
#' # plot for s=1,3,5,7,...,999 (s must be in (0,inf))
#' palette(rainbow(10))
#' xs=c(seq(0,0.01,length.out=200),seq(0.02,1,length.out=200))
#' # function phi for a fixed s
#' foo <- function(u) x@phi(u,1)
#' ys=foo(xs)
#' # plot for s=1
#' plot(x=xs,y=ys, col=1, ylim=c(0,0.4), xlim=c(0,1),type = "l") #plot(foo, col=1, ylim=c(0,0.94), xlim=c(0,1)) # ylim=c(0,1e-88)
#' # repeat for  s=3,5,7,...,999
#' for (i in 2:500) {
#'   foo <- function(u) x@phi(u,i*2-1)
#'   ys=foo(xs)
#'     lines(x=xs,y=ys, col=i+1) #plot(foo, add=TRUE, col=i+1)
#' }
#' # for a Gamma PUC the alphs are densities of inverse Pareto distributions...
#' plot(x@alphs[[1]], xlim=c(0,500))
#' # the fdenss are basically a normalization of the phis using alphs - for a fixes s they are densities of exponentially transformed Gamma distributions
#' heatplot:
#' us <- seq(0,1,length.out=100)
#' ss <- seq(1,500,length.out=100)
#' zs <- matrix(nrow=length(us),ncol=length(ss))
#' for (s in 1:length(ss)) zs[,s] <- x@fdenss[[1]](us, s)
#' image(x=us, y=ss, z=zs, col=terrain.colors(15))
#'
#' #For the simulation from the Gamma copula, the slots alphsCDF and Qk/corresponnding quantile are used internally
#' # plot the CDF corresponding to density alphs[1]
#' plot(x@alphsCDF[[1]],xlim=c(0,100))
#' x@alphsobjective[[1]](5,0.8)
#' x@alphsquant[[1]](0.8)
#' x@rand(1)
#'
#' # Create a beta copula
#' x <- PUCopula(family="beta", pars.a=c(40, 43), patch="lFrechet",data=stormflood)
#' x@phi(0.5,0.5)
#' #[1] 2.210851e-11
#' x@alph(0.6)
#' #[1] 1.458522e-11
#' #x@fdens(0.5,0.6)
#' #[1] 0.0002381545
#' #x@gdens(0.5,0.6)
#' #[1] 0.0001616076
#' x@fdenss[[1]](0.5,0.6)
#' #[1] 0.0002381545
#' x@fdenss[[2]](0.5,0.6)
#' #[1] 0.0001616076
#' x@rand(5)
#' #does not work
#' # example 2: gamma copula
#' x <- PUCopula(family="gamma", pars.a=c(40, 43), patch="lFrechet",data=stormflood)
#' #plot alpha(s)
#' plot(x@alphs[[1]])
#' #the numerically computed alpha is quite close:
#' s<-10
#' x@pars.a[1]*s^(x@pars.a[1]-1)/((1+s)^(x@pars.a[1]+1))
#' #[1] 0.008034519
#' x@alph(s)
#' #[1] 0.008034519
#' # but of course shows some inaccuracies:
#' s<-2.5
#' x@alph(s)
#' #[1] 6.530307e-06
#' x@pars.a[1]*s^(x@pars.a[1]-1)/((1+s)^(x@pars.a[1]+1))
#' #[1] 6.530261e-06
PUCopula <- setClass("PUCopula",
                     slots = c(
                       dim = "numeric",
                       family = "character",
                       par.factor = "numeric",
                       pars.a = "numeric",
                       p = "matrix",
                         phi = "function",
                         psy = "function",
                         phis = "list",
                          continuous = "logical",
               alph = "function",
               beta = "function",
               alphs = "list",  # this shall replace alph in the future
               alphsCDF = "list",
                  opstart = "list", # alternative, slower, more exact? less robust?
                  cdflim = "list",  #value, up to which cdf works without workaround for tails
                  cdfdom = "list",  #effective domains of the cdfs
                  alphsquantf = "list", # fast version, using spline approximation
                  cdfappx = "list",
                  #cdfast = "list",
               alphsobjective = "list",
               alphsquant = "list",
                       fdens = "function",
                       gdens = "function",
                       fdenss = "list",
               cpatch = "function",
               PUdens = "function",
                       fq = "function",
                       gq = "function",
                       patch = "character",
               dPWork = "function",
                       data = "matrix",
                       ranks = "matrix",
                       relRanks = "matrix",
               ## some methods (density, simulation)
                       rpatch = "function",
                       rand = "function"
                     )
                     #ranks = NULL
                     #validity=function(object){
                     #  if((object@par.a < 0) || (object@par.b < 0)) {
                     #    return("A negative number for one of the coordinates was given.")
                     #  }
                     #  return(TRUE)
                     #},
)

#' @describeIn PUCopula initializes a PUcopula object
setMethod("initialize", "PUCopula", function(.Object, dimension=0, factor=1, family=c("binom","nbinom","poisson","sample","gamma","beta"), pars.a=c(10,10), patch=c("rook","uFrechet","lFrechet","varwc","Gauss"), data, continuous=logical(0)) {
  .Object@dim = dimension # Dimension der Copula
  .Object@family = family # welcher Copula-Typ
  .Object@par.factor = factor #parameter fuer unterteilung patchwork?
  .Object@pars.a = pars.a # Parameter fur Verteilung von phi
  .Object@patch = patch # entweder patchworktyp oder komplett definiertes patchwork
  .Object@data = data # Slot für Daten
  # Ränge berechnen sofern nicht vorhanden; ebenfalls relative Ränge
  if(length(.Object@ranks)==0) {
    .Object@ranks <- as.matrix(apply(.Object@data,2,rank)/.Object@par.factor)
    if(length(.Object@relRanks)==0) {
      par.m <- dim(.Object@ranks)[1]
      .Object@relRanks <- as.matrix(.Object@ranks/(1+par.m))
    }
  }
  # wenn Dimension nicht angegeben wurde (==0), lese Dimension aus .@ranks
  if(dimension==0) {
    .Object@dim <- dim(.Object@ranks)[2]  # replace x@ranks by .Object@
  }

         # binom = {
         #   .Object@phi = function(u,s) {
         #     if (u%%1==0) {
         #       warning(paste("non-integer value for u =",u))
         #       return(0)
         #     }
         #     dbeta(x=u,shape1=s,shape2=pars.a[1]+1-s)/pars.a[1] #pars.a[2]???
         #   }
         # },
         # nbinom = {
         #   .Object@phi = function(u,s) {
         #     if (s%%1==0) {
         #       warning(paste("non-integer value for s =",s))
         #       return(0)
         #     }
         #     dbeta(x=u,shape1=s+1,shape2=pars.a[1]+1)*pars.a[1]/((pars.a[1]+s)*(pars.a[1]+s+1)) #pars.a[2]???
         #   }
         # } ,
         # poisson = {
         #   .Object@phi = function(u,s) {
         #     (1-u)^pars.a[1]*( pars.a[1]^s * ((-log(1-u))^s)/factorial(s)  )  #pars.a[2]???
         #   }
         # },

  fam2phi <- function(family, par.a) {
    force(par.a)
    force(family)
    switch(family,
           binom = {
             phi = function(u,s) {
               if (s%%1!=0) {
                 warning(paste("non-integer value for s =",s))
                 return(rep(0,length(u)) )
               } else
                 dbinom(x=s,size=par.a-1,prob=u)#/par.a
             }
           },
           nbinom = {
             phi = function(u,s) {
               if (s%%1!=0) {
                 warning(paste("non-integer value for s =",s))
                 return(rep(0,length(u))  )
               } else
                 dnbinom(x=s,size=par.a,prob=1-u)
             }
           } ,
           sample = {
             #     rr<-.Object@relRanks     # nicht idp von .Object!!!
             #     m <- .Object@par.factor
             phi = function(u,s) {
               if (s%%1!=0) {
                 warning(paste("non-integer value for s =",s))
                 return(rep(0,length(u)))
               } else {
                 if (s<1) {
                   warning("s must be a positive integer")
                   return(rep(0,length(u)))
                 }
               }
               if (par.a<=1) {
                 warning(paste("par.a too small, must be integer > 1 but is",par.a))
                 return(rep(0,length(u)))
               }

               rr<-.Object@relRanks     # nicht idp von .Object!!!
               #get partition borders for i=s
               m <- par.a #.Object@par.factor
               n <- dim(rr)[1]
               pn1 <- c(s-1,s)*(1/(m-1))
               pn2r <- sum(rr[,1,drop=FALSE]<=pn1[2])/n
               pn2l <- sum(rr[,1,drop=FALSE]<=pn1[1])/n
               #if ((u<=pn2r)&(u>pn2l)) return(1)
               #else return(rep(0,length(u)))
               return((u<=pn2r)*(u>pn2l))
             }
           } ,
           poisson = {
             phi = function(u,s) {
               if (s%%1!=0) {
                 warning(paste("non-integer value for s =",s))
                 return(rep(0,length(u))  )
               } else
                 dpois(x=s, lambda=par.a*(-1)*log(1-u) )
             }
           },
           gamma = {
             phi = function(u,s) {
               if (s<0) { #if ((s<0)|(s>1)) {
                 warning(paste("outside range (0,inf) for s =",s))
                 return(rep(0,length(u))  ) #macht fehler :() ## unnötig? , da dgamma auch 0 zurück gibt?
               }
               # do this in a way which is compatible with vectors
               zero_ind <- c(which(u==0),which(u==1))
               # create 0-vector
               returnVec <- numeric(length(u))
               # compute density where u is not -Inf and not 1   ## aber u > 0!!!!!!!!!!!!!
               nzero_ind <- setdiff(1:length(u),zero_ind) # necessary to define this, bec u[-zero_ind] does not work for empty set of indices
               returnVec[nzero_ind] <- dgamma(x=s, shape=par.a, rate = -log(u[nzero_ind]))
               return(returnVec)
             }
           },
           beta = {
             phi = function(u,s) {
               result <- dbeta(x=s, shape1=(-1)*log(u), shape2 = par.a  )   #shape1=(-1)*log(1-u)# oder: -log(u)
               return(result*(s<=1)*(s>=0))
               # 0<u,s<1, a in N
             }
           }
    )
    return(phi)
  }
  fam2phiV <- function(fam,a,d) {
    fam <- rep_len(fam,d)
    a <- rep_len(a,d)
    result <- list()
    for (i in 1:d) {
      result[[i]] <- fam2phi(fam[i], a[i])
    }
    return(result)
  }
  #Liste fuellen (aber nur wenn nicht schon besetzt?)
  .Object@phis <- fam2phiV(fam=.Object@family,a=.Object@pars.a,d=.Object@dim)
  .Object@psy = fam2phi(rep_len(.Object@family,2)[2], rep_len(.Object@pars.a,2)[2])
  .Object@phi = fam2phi(rep_len(.Object@family,2)[1], rep_len(.Object@pars.a,2)[1])
  fam2cont <- function(family) {
    if ("1" %in% c("binom","nbinom","sample","poisson")) {
      continuous = FALSE
    }
    else {
      if (family %in% c("gamma","beta")) {
        continuous = TRUE
      } else {  ### das sollte eigentlich nur geschehen wenn Wert nicht vorgegeben:
         continuous = .Object@phis[[1]](u=0.5,s=0.9)>0  ## ob das so funktiooniert-....
       }
    }
    return(continuous)
  }
  fam2contV <- function(fam,d) {
    fam <- rep_len(fam,d)
    result <- logical(d)
    for (i in 1:d) {
      result[i] <- fam2cont(fam[i])
    }
    return(result)
  }
  if (length(.Object@continuous)==0) {
    .Object@continuous = fam2contV(fam=.Object@family,d=.Object@dim)
  }


  #.Object@alph = function(s)  return(integrate(f=function(u) {x@phi(u,s)}, lower=0, upper=1, rel.tol = 1e-15)$value)
  .Object@alph = Vectorize(function(s)  return(integrate(f=function(u) {.Object@phi(u,s)}, lower=0, upper=1, rel.tol = 1e-8)$value), vectorize.args=list("s"))
  # und wewnn diskret??
 # Toleranz wegen: > integrate(f=function(u) {x@phi(u,0.021)}, lower=0, upper=1, rel.tol = 1e-15)$value
 # Error in integrate(f = function(u) { : the integral is probably divergent
  .Object@beta = Vectorize(function(s)  return(integrate(f=function(u) {.Object@psy(u,s)}, lower=0, upper=1, rel.tol = 1e-8)$value), vectorize.args=list("s"))
  #alph fuer alle phis:
  phis2alphs <- function(phis) {
    deffun <- function(i) {
      force(i)
      Vectorize(function(s) {
       # return(integrate(f=function(u) {phis[[i]](u,s)}, lower=0, upper=1, rel.tol = 1e-8)$value)   # fehler, weil zu konstant: which(abs(diff(x@alphs[[1]](seq(0,1,length.out=n)[c(-1, -n)])))>1.0e-20)
        ## better use pracma::integral:
        return(pracma::integral(fun=function(u) {phis[[i]](u,s)}, method = "Kronrod", xmin=0, xmax=1, reltol = 1e-5))
        ## wenn zu langsam, dann reltol = 1e-4 statt 1e-8
        ## achtung, fehler immer noch bei: pracma::integral(fun=function(u) {x@phis[[1]](u,0.999999999)}, method = "Kronrod", xmin=0, xmax=1, reltol = 1e-8)
      }  , vectorize.args=list("s"))
    }
   # alphs <- list()
   # for (i in 1:length(phis)) {
  #    alphs[[i]] <- deffun(i)
  #  } #kürzer:
    alphs <- lapply(1:length(phis), deffun)
    return(alphs)
  }
  .Object@alphs = phis2alphs(.Object@phis)

  .Object@cdflim = list()

  alphs2cdfs <- function(alphs) {
    deffun <- function(i) {
      force(i)
      Vectorize(function(s) {
        #return(integrate(f=alphs[[i]], lower=0, upper=s, rel.tol = 1e-4)$value)  #int von 0 oder -Inf
        if (length(.Object@cdflim)>0) if (s>.Object@cdflim[[i]][2] & .Object@cdflim[[i]][2]<1e+100) {
          limup <- .Object@cdflim[[i]][2]
          return(s/(s+((limup/.Object@alphsCDF[[i]](limup))-limup)))
        } #else{
          if (suppressWarnings(alphs[[i]](0)==Inf)) { #war:suppressWarnings(x@alphs[[1]](0)==Inf))
            return(1-pracma::integral(fun=alphs[[i]], method = "Kronrod", xmin=s, xmax=1, reltol = 1e-3) )
          } else {
            return(pracma::integral(fun=alphs[[i]], method = "Kronrod", xmin=0, xmax=s, reltol = 1e-3) )
          }
        #}
      }  , vectorize.args=list("s"))
    }
    cdfs <- lapply(1:length(alphs), deffun)
    return(cdfs)
  }
  .Object@alphsCDF = alphs2cdfs(.Object@alphs) # problematisch für große werte...

  ## ^^ herausfinden ab welchem wert problematisch
  ## 1. startwert für optimierung bestimmen
  nlmstart <- function(alphs) {
    deffun <- function(i) {
      force(i)
      rn <- c(rnorm(100), runif(100,-1,1), rnorm(100,0,100))
      return(rn[which.min(abs(.Object@alphs[[i]](rn)-0.5))])  # als startwert für findquantile
      # Vectorize(function(s) {
      #
      # }  , vectorize.args=list("s"))
    }
    rs <- lapply(1:length(alphs), deffun)
    return(rs)
  }
  .Object@opstart = suppressWarnings(nlmstart(.Object@alphs))

  cdfmima <- function(alphsCDF) {
    seqtails <- function(length.out=100) {
      s <- sort(unique(c(seq(-1,1,length.out=10), #rnorm(100),runif(100,-1,1),rnorm(100,0,100),
                         1*10^(1:10),
                         -1*10^(1:10),
                         1*10^((2:10)*10),
                         -1*10^((2:10)*10))))
      return(s)
    }
    deffun <- function(i) {
      force(i)
      print(paste("please wait for cdfmima to finish",i))
      xs <- seqtails()
      ps <- alphsCDF[[i]](xs)

      xs.fin <- xs
      xs.fin[is.finite(ps)==F] <- NA
      xs.fin[is.nan(ps)==T] <- NA
      ps.fin <- ps
      ps.fin[is.finite(ps)==F] <- NA
      ps.fin[is.nan(ps)==T] <- NA

      return(xs.fin[c(which.min(ps.fin),which.max(ps.fin))])
    }
    rs <- lapply(1:length(alphsCDF), deffun)
    return(rs)
  }
  .Object@cdflim = suppressWarnings(cdfmima(.Object@alphsCDF))

  getcdfdom <- function(alphsCDF) {
    deffun <- function(i) {
      force(i)

      print(paste("please wait for cdfdom to finish",i))
      objective <- function(x, quantile){
        (alphsCDF[[i]](x) - quantile)^2
      }
      find_quantile <- function(quantile){
        result = nlminb(start=.Object@opstart[[i]], objective=objective,
                        quantile = quantile)$par
        return(result)
      }
      itv <- c(find_quantile(quantile = 0.05),
               find_quantile(quantile = 0.95))
      itv2 <- itv + c(-diff(itv),diff(itv))


      return(itv2)
    }
    rs <- lapply(1:length(alphsCDF), deffun)
    return(rs)
  }

  .Object@cdfdom = getcdfdom(.Object@alphsCDF)

  cdfapprox <- function(alphsCDF) {
    deffun <- function(i) {
      force(i)

      print(paste("please wait for cdf being approximated",i))

      pts <- seq(.Object@cdfdom[[i]][1],.Object@cdfdom[[i]][2], length.out = 100)
      pts.int <- alphsCDF[[i]](pts)

      qtfun <- approxfun(c(pts.int,0,1),c(pts,-Inf,Inf))

      return(qtfun)
    }
    rs <- lapply(1:length(alphsCDF), deffun)
    return(rs)
  }

  .Object@cdfappx= cdfapprox(.Object@alphsCDF)

  splineqt <- function(cdfappx) {
    get.mimaq <- function(fun) {
      sequ.tails <- function(length.out=100) {
        s <- sort(unique(c(1*10^(-60:-3), seq(0,1,length.out=length.out) ,(1-1*10^(-60:-3)))))
        return(s)
      }
      qs <- sequ.tails()
      xs <- fun(qs)
      xs.fin <- xs
      xs.fin[is.finite(xs)==F] <- NA
      xs.fin[is.nan(xs)==T] <- NA
      qs.fin <- qs
      qs.fin[is.finite(xs)==F] <- NA
      qs.fin[is.nan(xs)==T] <- NA
      return(qs.fin[c(which.min(qs.fin),which.max(qs.fin))])
    }

    deffun <- function(i) {
      force(i)
      print(paste("please wait for splineqt to finish",i))

      fcdf <- function(x) {
        rs <- numeric(length(x))

        mimax <- get.mimaq(cdfappx[[i]]) # kleiner und großer funktionierender wert
        mimay <- cdfappx[[i]](mimax)
        qm <- cdfappx[[i]](0.5)
        k2 <- (mimay[2]-qm)*(2-2*mimax[2])/(2*mimax[2]-1)
        k1 <- (qm-mimay[1])*2*mimax[1]/(1-2*mimax[1])
        tl2 <- function(x) {
          qm+(((2*x-1)*k2)/(2-2*x))
        }
        tl1 <- function(x) {
          qm-(((1-2*x)*k1)/(2*x))
        }

        rs[x<mimax[2]&x>mimax[1]] <- qtfun(x[x<mimax[2]&x>mimax[1]])
        rs[x>mimax[2]] <- tl2(x[x>mimax[2]])
        rs[x<mimax[1]] <- tl1(x[x<mimax[1]])
        return(rs)
      }
    }
    rs <- lapply(1:length(cdfappx), deffun)
    return(rs)
  }

  .Object@alphsquantf = splineqt(.Object@cdfappx)

  cdfs2objv <- function(cdfs) {
    deffun <- function(i) {
      force(i)
        return( function(x, quantile){ (cdfs[[i]](x) - quantile)^2 })
    }
    objv <- lapply(1:length(cdfs), deffun)
    return(objv)
  }
  .Object@alphsobjective = cdfs2objv(.Object@alphsCDF)

  objv2quant <- function(objv) {
    deffun <- function(i) {
      force(i)
      return( Vectorize( function(quantile){ nlminb(start=0.5, objective=objv[[i]],  quantile = quantile)$par }, vectorize.args = list("quantile")) )
    }
    quant <- lapply(1:length(objv), deffun)
    return(quant)
  }
  .Object@alphsquant = objv2quant(.Object@alphsobjective)

  phis2dens <- function(phis,alphs) {
    deffun <- function(i) {
      force(i)
      function(u,s, log=FALSE) {if (log==FALSE) {phis[[i]](u,s)/alphs[[i]](s)} else {
        log(phis[[i]](u,s)/alphs[[i]](s)) # log(phis[[i]](u,s))-log(alphs[[i]](s)) # wie numerisch stabiler?
        } }
    }
    dens <- lapply(1:length(phis), deffun)
    return(dens)
  }
  .Object@fdens = function(u,s) return(.Object@phi(u,s)/.Object@alph(s))
  .Object@gdens = function(u,s) {.Object@psy(u,s)/.Object@beta(s)}
  .Object@fdenss = phis2dens(.Object@phis,.Object@alphs)
  #einfach mal standardmäßig Gauss(0.8) fuer cpatch festsetzen
  .Object@cpatch = function(s,t) { copula::dCopula(u=c(s,t), copula = copula::normalCopula(0.8, dim=2)) }
  #density of PU copula statt Inf nur 100... :/
  .Object@PUdens = function(u,v) {
    integrate(function(t) {
      sapply(t, function(t) {
        integrate(function(s) {
          sapply(s, function(s) return(.Object@cpatch(.Object@alph(s),.Object@beta(t))*.Object@phi(u,s)*.Object@psy(v,t)) ) },0,100, rel.tol = 1e-7)$value
        #bei 0: -Inf, bei 1: 0
      })
    },0,100, rel.tol = 1e-4)
  }

  # nur bivariat; patch=Gauss
  .Object@dPWork = function(u,v) {
    n <- dim(.Object@data)[1]
    getpatch <- (u<=(.Object@ranks[,1]/n)) * (u>((.Object@ranks[,1]-1)/n)) * (v<=(.Object@ranks[,2]/n)) * (v>((.Object@ranks[,2]-1)/n))
    k<-which(getpatch==1)
    if (length(k)==0) return(0)
    else {
      return(copula::dCopula(c(n*u-x@ranks[k,1]+1,n*v-x@ranks[k,2]+1),copula::normalCopula(0.8, dim=2)) )
    }
  }
  #simulating: patchwork simulation
  .Object@rpatch = function(n=1, patch = .Object@patch, patchpar=NULL){
    # step 1
    rsims.index <- ceiling(runif(n)*dim(.Object@ranks)[1])
    rsims <- as.matrix(.Object@ranks[rsims.index,,drop=FALSE])
    # step 2
    usims <- matrix(runif(.Object@dim*n),nrow=n,ncol=.Object@dim)
    # step 3
    par.m <- dim(.Object@ranks)[1] # oder in obj speichern? #Anz Zeilen/Beobachtungen
    switch(patch,
           none = {Z <- (rsims-1.0)/par.m},
           rook = {Z <- (rsims-usims)/par.m},
           lFrechet = {Z <- cbind(rsims[,1]-usims[,1],rsims[,2]+usims[,1]-1)/par.m}, #nur dim 2
           uFrechet = {Z <- (rsims-usims[,rep(1,.Object@dim)])/par.m},
           Gauss = {
             if (is.null(patchpar)) warning("patchpar must not be NULL when patch is Gauss")
             Z <- (rsims-1+copula::rCopula(n,copula::normalCopula(patchpar, dim=.Object@dim))  )/par.m})
    return(Z)
  }
  #simulating: main function
  .Object@rand =  function(n=1) {
  # step 1 (select random pair of ranks)
  rsims.index <- ceiling(runif(n)*dim(.Object@ranks)[1])
  rsims <- as.matrix(.Object@ranks[rsims.index,,drop=FALSE])
  # step 2 (urvs for each dimension/observation)
  usims <- matrix(runif(.Object@dim*n),nrow=n,ncol=.Object@dim)
  # step 3
  par.m <- dim(.Object@ranks)[1] # oder in obj speichern? #Anz Zeilen/Beobachtungen
  switch(.Object@patch,#match.arg(.Object@patch),
  none = {Z <- (rsims-1.0)/par.m},
  rook = {Z <- (rsims-usims)/par.m},
  lFrechet = {Z <- cbind(rsims[,1]-usims[,1],rsims[,2]+usims[,1]-1)/par.m}, #nur dim 2
  uFrechet = {Z <- (rsims-usims[,rep(1,.Object@dim)])/par.m},
  Gauss = {Z <- (rsims-1+copula::rCopula(n,copula::normalCopula(0.8, dim=obj@dim))  )/par.m})
  # step 4
  # switch(.Object@family,#match.arg(.Object@family),
  # binom = {d<-ceiling(sweep(Z,2,.Object@pars.a,"*"))},
  # nbinom = {d<-floor(sweep(Z/(1-Z),2,.Object@pars.a,"*"))},
  # sample = {
  #           rr<-t(matrixStats::colRanks(Z)-0.5)/dim(Z)[1]  # colRanks(Z)-0.5 wie oben mit stetigkeitskorrektur # var ranks # rel ranks
  #           foo <- function(rrx,n) {table(cut(rrx,breaks=seq(0,1,length.out=n+1), ordered_result=FALSE))}
  #           print(.Object@pars.a) # our "n"; mv?
  #           newpartition <- apply(rr,2,foo, n=.Object@pars.a)/dim(rr)[1]  # ^^verwendet f?r alle spaltengleichen parameter
  #           newpartition <- apply(newpartition,2,cumsum)
  #           p_index <- function(u,part) {      #  function f?r eine Zeile i von newppartition:
  #           test<- matrix(t(u),dim(part)[1],dim(part)[2],TRUE)<part
  #           test[which(!test)]<-NA
  #           return(apply(test*part, 2, which.min))
  #           }
  #           d<-t(apply(Z,1,p_index,part=newpartition))
  #           },
  # poisson = {d<-floor(sweep(-log(1-Z),2,(log(.Object@pars.a+1)-log(.Obect@pars.a)),"/"))} )
            #alternative numerically
            d <- Z # d <- matrix()
            print("Z")
            print(Z)
           # for (i in 1:.Object@dim) d[,i] <- .Object@alphsquant[[i]](Z[,i]) #richtig?
            for (i in 1:.Object@dim) d[,i] <- .Object@alphsquantf[[i]](Z[,i]) #richtiger?

            #step next: result possesses the joint distribution P with desired marginals
            print(paste("d",d))
            rs <- d
            for (i in 1:.Object@dim) for (j in 1:n)  {
              print(paste("i",i,"j",j))
              print(paste("rs[j,i]",rs[j,i]))
              print(paste("d[j,i]",d[j,i]))
           #   print(paste("replacement",.mcmcMH(target = function(theta) { return(.Object@fdenss[[i]](u=theta,s=d[j,i],log=TRUE)) } ) ))
           #   rs[j,i] <- .mcmcMH(target = function(theta) { return(.Object@fdenss[[i]](u=theta,s=d[j,i],log=TRUE)) } )       #ob das funktioniert, index i korrekt benutzt??
             # rs[j,i] <- .mh_sample(dens = function(theta) { return(.Object@fdenss[[i]](u=theta,s=d[j,i],log=F)) } )       #ob das funktioniert, index i korrekt benutzt??
              smpl <- armspp::arms(5000, function(theta) { return(.Object@fdenss[[i]](u=theta,s=d[j,i],log=F)) } , 0, 1)#smpl <- armspp::arms(5000, .Object@fdenss[[i]](u=theta,s=d[j,i],log=T), 0, 1)
              rs[j,i] <- smpl[length(smpl)]
            }
            #^^ ranks ist falsch hier?
            return(rs)
  # step 5
  # switch(.Object@family,#match.arg(.Object@family),
  # binom = {return(qbeta( matrix(runif(n*.Object@dim),nrow=n,ncol=.Object@dim),
  #           d,
  #           matrix(.Object@pars.a+1,nrow=n,ncol=.Object@dim,byrow=TRUE)-d))},
  # nbinom = {return(qbeta( matrix(runif(n*.Object@dim),nrow=n,ncol=.Object@dim),
  #           d+1,
  #           matrix(.Object@pars.a+1,nrow=n,ncol=.Object@dim,byrow=TRUE)))} ,
  # sample = {  unif_i <- function(i,partition) {
  #           npart<-rbind(0,partition)
  #           return(runif(1,min=npart[i],max=npart[i+1]))
  #           }
  #           return( apply(d,1:2,unif_i,partition=newpartition) )} ,
  # poisson = {return(1-exp(-qgamma( matrix(runif(n*.Object@dim),nrow=n,ncol=.Object@dim),
  #           shape=d+1,
  #           scale=1/(1+matrix(.Object@pars.a+1,nrow=n,ncol=.Object@dim,byrow=TRUE))))) }
  # )

}

  return(.Object)
})



