library(maps)
library(dbplyr)
library(leaflet)
library(dse)
library(KFAS)
library(geosphere)
library(foreach)
library(FKF)
library(ggplot2)
library(gridExtra)

# SIMULATION 
# simulate short rate process and zero-coupon rates under Vasicek model
# ----------------------------------------------------------------------------------------
# Model setup

# 1 FACTOR
# kappa     <- c(0.06) 
# theta     <- c(0.05)
# sigma     <- c(0.02)
# lambda    <- c(-0.20)

# 2 FACTOR
# kappa     <- c(0.06, 0.7) 
# theta     <- c(0.05, 0.01)
# sigma     <- c(0.02, 0.05)
# lambda    <- c(-0.20, -0.50)

kappa     <- c(0.06, 0.3, 0.7) 
theta     <- c(0.01, 0.02, 0.04)
sigma     <- c(0.02, 0.05, 0.03)
lambda    <- c(-0.20, -0.50, -0.15)
T         <- 10 # years
delta_t   <- 1/12 # monthly zeros observations
dt        <- 1/12 # monthly r simulations
n         <- T/dt # number of r simulations
r_0       <- theta
measurement_error <- 0.001 # for zero-coupon rates
m         <- length(kappa)   # dimension of state variables
# maturity  <- c(1/12 ,1/4, 1/2, 10) # zeros for 1 factor simulation 
# maturity  <- c(1/12 ,1/4, 1/2, 1, 2, 3, 5, 7, 10) # zeros for 2 factor simulation
maturity  <- c(1/12 ,1/4, 1/2, 1, 2, 3, 4, 5, 7, 10, 12, 15, 20, 30) # zeros for 3 factor simulation
d         <- length(maturity)   # dimension of observations
N         <- T/delta_t # number of observations
# -----------------------------------------------------------------------------------------

# Simulate zero rates driven by short rate process 
# Generate time series of N observations
simulate_zero_rates <- function(null_arg=NULL) 
{
  transition_density_mean <- function(r_t)
  { 
    return((theta)*(1-exp(-kappa*dt)) + exp(-kappa*dt)*r_t);
  }
  transition_density_sd <- sqrt(sigma^2/(2*kappa)*(1-exp(-2*kappa*dt)))
  
  simulated_r <- matrix(rep(0,n*m),nrow=n)
  simulated_r[1,] <- r_0 
  foreach(i=2:n) %do% 
    {
      simulated_r[i,] <- rnorm(m,transition_density_mean(simulated_r[i-1,]),transition_density_sd);
    }
  
  # Zero-coupon rates driven by simulated short rate
  # Affine model: P(t,T) = exp(A(t,T) + B(t,T)r(t))
  kappa <- matrix(kappa,d,m,T)
  theta <- matrix(theta,d,m,T)
  sigma <- matrix(sigma,d,m,T)
  lambda <- matrix(lambda,d,m,T)
  maturity <- matrix(maturity,m,d,T)
  B <- 1/kappa*(1-exp(-kappa*t(maturity)))
  gamma <- kappa^2*(theta-sigma*lambda/kappa) - 0.5*sigma^2 
  A <- t(gamma[1,]/(kappa[1,]^2)) %*% t(B-t(maturity)) + t(-sigma[1,]^2/(4*kappa[1,])) %*% t(B^2)
  A <- t(A)
  
  measurement_noise <- measurement_error * matrix(rnorm(length(maturity[1,])*n),ncol=d)
  simulated_zero_rates <- matrix(rep(-A/maturity[1,],n),ncol=d,byrow=TRUE) + t((B/t(maturity)) %*% t(simulated_r)) + measurement_noise 
  return(list(r=simulated_r, zero_rates=simulated_zero_rates))
}


# PARAMETER ESTIMATION
# ----------------------------------------------------------
# Random initialization of parameters
init_params <- function()
{
  kappa_init <<- runif(m, min=0.0, max=0.5)
  theta_init <<- runif(m, min=0.0, max=0.5)
  sigma_init <<- runif(m, min=0.0, max=1.0)
  lambda_init <<- runif(m, min=-1.0, max=1.0)
  measurement_err_init <<- runif(d, min=0.0, max=0.1)
}
# optimization parameter bounds
upper_bound <- c(rep(c(1.0, 1.0, 1.0, 1.0),each=m), rep(0.1,d))
lower_bound <- c(rep(c(0.0001, 0.0001, 0.0001,-1.0 ),each=m), rep(0.0001,d))
actual_param <- c(kappa=kappa, theta=theta, sigma=sigma, lambda=lambda, err=rep(measurement_error,d))
# ----------------------------------------------------------

vasicek_KF <- function(kappa, theta, sigma, lambda, measurement_errs, observations)
{ 
  # initial state variable (a0: m x 1)
  r_init <- as.vector(theta)
  
  # variance of state variable (P0: m x m)
  P_init <- (sigma^2/(2*kappa))*diag(1,m,m) # unconditional variance of state variable
  
  # intercept of state transition equation (dt: m x 1)
  C <- matrix((theta)*(1-exp(-kappa*delta_t)))
  
  # factor of transition equation (Tt: m x m x 1)
  F_ <- array(exp(-kappa*delta_t)*diag(m),dim=c(m,m,1))
  
  # factor of measurement equation (Zt: d x m x 1)
  B <- array(1/matrix(rep(kappa,d),d,byrow=TRUE)*(1-exp(-matrix(rep(kappa,d),d,byrow=TRUE) * matrix(rep(maturity,m),d))),dim=c(d,m,1))
  
  # intercept of measurement equation (ct: d x 1)
  gamma <- kappa^2*(theta-sigma*lambda/kappa) - 0.5*sigma^2
  A <- t(gamma/(kappa^2)) %*% t(B[,,1]-matrix(rep(maturity,m),d)) + t(-sigma^2/(4*kappa)) %*% t(B[,,1]^2)
  A <- matrix(-A/maturity,nrow=length(maturity))
  
  B <- array(B[,,1]/matrix(rep(maturity,m),d),dim=c(d,m,1))
  
  # variance of innovations of transition (HHt: m x m x 1)
  Q <- array(sigma^2/(2*kappa)*(1-exp(-2*kappa*delta_t))*diag(m),dim=c(m,m,1))
  
  # variance of measurement error (GGt: d x d x 1)
  R <- array(diag(d)*measurement_errs^2,dim=c(d,d,1))
  
  filtered_process <- fkf(a0=r_init, P0=P_init, dt=C, ct=A, Tt=F_, Zt=B, HHt=Q, GGt=R, yt=t(observations))
  return(filtered_process)
}

# Retrieve short rates using Kalman Filter 
retrieve_short_rates <- function(rates, optim_controls, lower_bound=NULL, upper_bound=NULL) 
{
  observations <- rates
  init_params()  
  initial_param <<- c(kappa=kappa_init, theta=theta_init, sigma=sigma_init, lambda=lambda_init, err=measurement_err_init)
  
  vasicek_KF_loglik <- function(x)
  { 
    kappa <- x[1:m]; theta <- x[(m+1):(2*m)]; sigma <- x[(2*m+1):(3*m)]; lambda <- x[(3*m+1):(4*m)]; measurement_errs <- x[(4*m+1):length(x)] 
    return(-vasicek_KF(kappa,theta,sigma,lambda,measurement_errs,observations)$logLik)
  }
  
  # optimization of log likelihood function
  fitted_model <- nlminb(initial_param, vasicek_KF_loglik, control=optim_controls, lower=lower_bound, upper=upper_bound) 
  return(fitted_model)
}


run_simulation <- function(num_simulation)
{
  simulation_results <- vector('list', length=num_simulation)
  simulated_rates <- lapply(simulation_results,FUN=simulate_zero_rates)
  simulated_r <- lapply(simulated_rates,function(x) x$r)
  simulated_zeros <- lapply(simulated_rates, function(x) x$zero_rates)
  return(list(r=simulated_r,zeros=simulated_zeros))
}

run_KF_for_vasicek <- function(rates, num_state_vars, optim_controls, bounds, print_results=FALSE)
{
  num_zeros <- num_state_vars*4 + length(maturity)
  num_params <- num_state_vars*4+d
  fitted_results <- lapply(rates, FUN=retrieve_short_rates, optim_controls,bounds$lower,bounds$upper)
  fitted_parameters <- sapply(fitted_results,function(result) t(cbind(t(result$par),t(result$message))))
  fitted_parameters <- fitted_parameters[,fitted_parameters[length(fitted_parameters[,1]),]=="singular convergence (7)"]
  if(length(rates) && is.na(fitted_parameters[1])) {warning("optimization failed to converge"); return(NULL)}
  fitted_parameters <- apply(matrix(fitted_parameters,num_params+1),2,function(x) x[-(num_params+1)])
  fitted_parameters <- apply(fitted_parameters, 2, sort_param, num_state_vars)
  fitted_parameters <- matrix(as.numeric(fitted_parameters),num_state_vars*4+d)
  return(list(fitted_results,fitted_parameters))
}

sort_param <- function(params,num_state_vars)
{
  idx <- order(params[1:num_state_vars])
  param_seq <- seq(from=1, to=num_state_vars*4, by=num_state_vars)
  foreach(param_type=param_seq) %do%
    {
      params[param_type:(param_type+num_state_vars-1)] <- params[param_type:(param_type+num_state_vars-1)][idx]
    }
  return(params)
}

compute_param_stats <- function(fitted_parameters)
{
  stats <- data.frame(actual=actual_param)
  stats <- cbind(stats,mean=apply(fitted_parameters, 1, mean))
  stats <- cbind(stats, s.d=apply(fitted_parameters, 1, sd))
  return(round(stats, digits=4))
}


# ---------- Functions for printing results ----------------------------------------------------------------------

print_results_and_graph <- function(simulated_r,fitted_model,actual_param,initial_param,fitted_param,observations)
{
  m <- dim(simulated_r)[2]
  print(data.frame(actual=actual_param, start=initial_param, fitted=fitted_param))
  retrieved_r <- with(fitted_model, vasicek_KF(par[1:m], par[(m+1):(2*m)], par[(2*m+1):(3*m)], par[(3*m+1):(4*m)], par[(4*m+1):length(par)],observations))
  
  if(min(dim(simulated_r))>1 & min(dim(simulated_r))<4) plot_rates(simulated_r,KF_rates=retrieved_r$att) 
  plot_rates(rates=apply(simulated_r,1,sum), KF_rates=apply(retrieved_r$att,2,sum))
}

plot_rates <- function(rates, name="Rate", KF_rates, dates=NULL)
{ 
  if(is.null(dates)) dates <- 1:length(rates)
  KF_rates <- t(KF_rates)
  plot_rate <- function(rate, name, KF_rate, dates) {
    data <- data.frame(date=dates, simulated=rate, filtered_rates=KF_rate)
    p <- ggplot(data, aes(date)) +theme_bw() + scale_color_manual("",values=c("steelblue","red")) #61D7A4
    p <- p + geom_line(aes(y=simulated,color='simulated'),alpha=0.6) + geom_line(aes(y=filtered_rates, color='filtered'),size=1.7,alpha=0.3) 
    p <- p + ylab(name) + xlab('Time in Months') + ggtitle(paste('Simulated and Filtered',name))
    return(p)
  }
  
  if(min(dim(rates))>1 & min(dim(rates))<4) {
    graph <- list()
    par(mfrow=c(min(dim(rates)),1))
    for(factor in 1:min(dim(rates))) {
      p <- plot_rate(rates[,factor],paste("Factor",factor),KF_rates[,factor],dates)
      graph[[factor]] <- p
    }
    multiplot(plotlist=graph)
  } else { plot(plot_rate(rates,"Short Rate",t(KF_rates),dates)) }
  
  return()
}

# Multiple plot function
# from the Cookbook for R http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# plot parameter histograms
plot_hist <- function(obs, xlab="Value", plot_density=TRUE)
{
  #par(mfrow=c(4,1))
  obs <- data.frame(obs); names(obs) <- "Value";
  hist <- ggplot(obs, aes(Value)) + 
    geom_histogram(aes(y=..density..), color="steelblue", fill='white') + 
    labs(x=xlab, y="Density")
  
  if(plot_density)
  {
    hist <- hist + stat_function(fun=dnorm, args=list(mean=mean(obs$Value), sd=sd(obs$Value)), colour="orange",size=1)
  }
  return (hist)
}

# ------- End of Functions ----------------------------------------------------------------

# RUN SIMULATION AND KF FILTER 
set.seed(106)
optim_controls <- list(trace=100,eval.max=3000,iter.max=3000,rel.tol=1e-15) 
simulated_rates <- run_simulation(250) # 250 simulations takes long!!!
results <- run_KF_for_vasicek(rates=simulated_rates$zeros, num_state_vars=m, optim_controls=optim_controls, 
                              bounds=list(upper=upper_bound,lower=lower_bound), print_results=FALSE)
fitted_param <- results[[2]]
k <- 1 # 1 factor k=1; 2 factors k=3; 3 factors k=40
model <- results[[1]][[k]]
print_results_and_graph(simulated_rates$r[[k]], model, actual_param,initial_param,fitted_param[,k],simulated_rates$zeros[[k]])
compute_param_stats(fitted_param)


# Plot Parameter Histograms 
hists <- list()
hist_num <- 1; param <- 1;
while(param <= (4*m+d)) {
  hists[[hist_num]] <- plot_hist(fitted_param[param,],names(initial_param)[param])
  if(hist_num==6 | param==(4*m+d)) {
    multiplot(hists[[1]],hists[[2]],hists[[3]],hists[[4]],
              hists[[5]],hists[[6]],cols=2)
    hists <- list('','','','','','')
    hist_num <- hist_num-6
  }
  hist_num <- hist_num + 1
  param <- param + 1
}

# END


mapF <- leaflet() %>% addTiles() %>%
  addMarkers(lat=44.5884, lng=33.5224, popup = "1") %>%
  addMarkers(lat=44.5888, lng=33.5227, popup = "2") %>%
  addMarkers(lat=44.5890, lng=33.5237, popup = "3") %>%
  addMarkers(lat=44.5893, lng=33.5245, popup = "4") %>%
  addMarkers(lat=44.5899, lng=33.5240, popup = "5") %>%
  addMarkers(lat=44.5912, lng=33.5245, popup = "6") %>%
  addMarkers(lat=44.5919, lng=33.5252, popup = "7") %>%
  addMarkers(lat=44.5927, lng=33.5265, popup = "8") %>%
  addMarkers(lat=44.5933, lng=33.5272, popup = "9") %>%
  addMarkers(lat=44.5945, lng=33.5277, popup = "10") %>%
  addMarkers(lat=44.5944, lng=33.5289, popup = "11") 
mapF
#map %>% addProviderTiles(providers$Stamen.Toner)


m1b.dse <- dse::SS(F = matrix(1, 1, 1), Q = matrix(40, 1, 1),
                      H = matrix(1, 1, 1), R = matrix(130, 1, 1),
                      constants = list(Q = matrix(TRUE, 1, 1), P0 = matrix(TRUE, 1, 1)),
                      z0 = matrix(0, 1, 1), P0 = matrix(10^5, 1, 1))


data("Nile", package = "datasets")
m1b.dse.est <- estMaxLik(m1b.dse, TSdata(output = Nile))
data("Nile", package = "datasets")
m1.sspir <- sspir::SS(Fmat = function(tt, x, phi) {
  return(matrix(1))
   }, Gmat = function(tt, x, phi) {
     return(matrix(1))
     }, Vmat = function(tt, x, phi) {
       return(matrix(exp(phi[1])))
       }, Wmat = function(tt, x, phi) {
         return(matrix(exp(phi[2])))
         }, y = as.matrix(Nile, ncol = 1))
str(m1.sspir)




# 
# inter <- gcIntermediate(addMarkers_1, addMarkers_2,addMarkers_3, n=50, addStartEnd=TRUE, breakAtDateLine=F)
# lines(inter, col="slateblue", lwd=2)
# map('world',
#     col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
#     mar=rep(0,4),border=0, ylim=c(-80,80) 
# )
# 
# # Dot for cities
# points(x=data$long, y=data$lat, col="slateblue", cex=3, pch=20)
# 
# # Compute the connection between Buenos Aires and Paris
# inter <- gcIntermediate(addMarkers_1, addMarkers_2, n=50, addStartEnd=TRUE, breakAtDateLine=F)
# 
# # Show this connection
# lines(inter, col="slateblue", lwd=2)
# 
# # Between Paris and Melbourne
# inter <- gcIntermediate(addMarkers_2,addMarkers_3,  addMarkers_1, n=50, addStartEnd=TRUE, breakAtDateLine=F)             
# lines(inter, col="slateblue", lwd=2)

logistG <- function(r, p, k, t){
  k * p * exp(r*t) / (k + p * (exp(r*t) - 1))
}

k <- 100
p0 <- 0.1*k
r <- 0.2
deltaT <- 0.1

# Let's create some sample data:
set.seed(12345)

obsVariance <- 25
nObs = 250
nu <- rnorm(nObs, mean=0, sd=sqrt(obsVariance)) 
pop <- c(p0, logistG(r, p0, k, (1:(nObs-1))*deltaT)) + nu

Estimate <- data.frame(Rate=rep(NA, nObs),
                       Population=rep(NA,nObs))

library(numDeriv)
a <- function(x, k, deltaT){
  c(r=x[1],
    logistG(r=x[1], p=x[2], k, deltaT))
}
G <- t(c(0, 1))

# Evolution error
Q <- diag(c(0, 0))
# Observation error
R <-  obsVariance
# Prior
x <- c(r, p0)
Sigma <-  diag(c(144, 25))

for(i in 1:nObs){
  # Observation
  xobs <- c(0, pop[i])
  y <- G %*% xobs
  # Filter  
  SigTermInv <- solve(G %*% Sigma %*% t(G) + R)
  xf <- x + Sigma %*% t(G) %*%  SigTermInv %*% (y - G %*% x)
  Sigma <- Sigma - Sigma %*% t(G) %*% SigTermInv %*% G %*% Sigma 
  
  A <- jacobian(a, x=x, k=k, deltaT=deltaT)   
  K <- A %*% Sigma %*% t(G) %*% solve(G %*% Sigma %*% t(G) + R)
  Estimate[i,] <- x
  
  # Predict
  x <- a(x=xf, k=k, deltaT=deltaT) + K %*% (y - G %*% xf)
  Sigma <- A %*% Sigma %*% t(A) - K %*% G %*% Sigma %*% t(A) + Q
}

# Plot output
op <- par(mfrow=c(2,1))
time <- c(1:nObs)*deltaT

m <- leaflet() %>% 
  addTiles() %>% 
  setView( lng = 2.34, lat = 48.85, zoom = 4 ) %>% 
  addProviderTiles("NASAGIBS.ViirsEarthAtNight2012", group="????? 1") %>%
  addTiles(options = providerTileOptions(noWrap = TRUE), group="background 2") %>%
  
  addCircleMarkers(data=df, lng=df$Long , lat=df$Lat, radius=8 , color="black",
                   fillColor="red", stroke = TRUE, fillOpacity = 0.8, group="Red") %>%
  
  addLayersControl(overlayGroups ="Red" , baseGroups = c("????? 1","????? 2"), 
                   options = layersControlOptions(collapsed = FALSE))

mapT <- leaflet(df) %>% addTiles() %>%
  addMarkers(lat=44.5884, lng=33.5224, popup = "1") %>%
  addMarkers(lat=44.5886, lng=33.5227, popup = "2") %>%
  addMarkers(lat=44.5893, lng=33.5237, popup = "3") %>%
  addMarkers(lat=44.5897, lng=33.5240, popup = "4") %>%
  addMarkers(lat=44.5900, lng=33.5242, popup = "5") %>%
  addMarkers(lat=44.5912, lng=33.5247, popup = "6") %>%
  addMarkers(lat=44.5919, lng=33.5252, popup = "7") %>%
  addMarkers(lat=44.5927, lng=33.5268, popup = "8") %>%
  addMarkers(lat=44.5933, lng=33.5272, popup = "9") %>%
  addMarkers(lat=44.5942, lng=33.5277, popup = "10") %>%
  addMarkers(lat=44.5944, lng=33.5289, popup = "11") 

my_Clayton <- function(Clayton_copula,data,method = 'ml'){
  tryCatch(
    {
      y = fitCopula(Clayton_copula,data,method)
      return(y)
    },error=function(error_message) {
      return(NA)
    }
  )
}
marginal_straight<-function(vecX,distrType,paramVector){
  switch(distrType,
         normal=pnorm(vecX,mean=paramVector[1],sd=paramVector[2]),
         gamma=pgamma(vecX,shape=paramVector[1],rate=paramVector[2]),
         exponential=pgamma(vecX,shape=paramVector[1],rate=paramVector[2]),
         lognormal=plnorm(vecX,mean=paramVector[1],sd=paramVector[2]),
         logistic=plogis(vecX,location=paramVector[1],scale=paramVector[2]))
}

marginal_inverse<-function(vecX,distrType,paramVector){
  switch(distrType,
         normal=qnorm(vecX,mean=paramVector[1],sd=paramVector[2]),
         gamma=qgamma(vecX,shape=paramVector[1],rate=paramVector[2]),
         exponential=qgamma(vecX,shape=paramVector[1],rate=paramVector[2]),
         lognormal=qlnorm(vecX,mean=paramVector[1],sd=paramVector[2]),
         logistic=qlogis(vecX,location=paramVector[1],scale=paramVector[2]))
}


myvec <- c(Normal_fit@loglik, T_fit@loglik,loglic, Frank_fit@loglik)
(copulanum <- which (myvec == max(myvec)))
if (copulanum==1)
  best_parameters <- param1
if (copulanum==2)
  best_parameters <- param2
if (copulanum==3)
  best_parameters <- param3
if (copulanum==4)
  best_parameters <- param4

mapT

