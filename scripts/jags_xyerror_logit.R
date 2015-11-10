###Script to fit a logistic regression with a quasi-gaussian functional form and
###error-in-variables model in JAGS on simulated data.
###Ian Breckheimer
###November 7th 2015

####Simulate Data####
n <- 1000

##First covariate without measurement error
meanx1 <- 1
sdx1 <- 1
x1 <- rnorm(n,meanx1,sdx1)

##Second covariate with measurement error
x2mean <- 0
sdx2 <- 1
pop_taux2 <- 1 / (sdx2 * sdx2)
sdobs <- 0.1
taux2 <- 1 / (sdobs * sdobs)
truex2 <- rnorm(n, x2mean, sdx2)
errorx2 <- rnorm(n, 0, sdobs)
obsx2 <- truex2 + errorx2

##Grouping variable.
ngroups <- n/20
group_sd <- 0.5
groups <- 1:ngroups
nphotos <- rpois(n=ngroups,lambda = n/(ngroups-5))
photos_groups <- as.numeric(factor(rep(groups,times=nphotos)[1:n]))
means_groups <- rnorm(ngroups,0,group_sd)
grpmeans_y <- rep(means_groups,times=nphotos)[1:n]

height_int <- 0
height_slope <- 0.25
width_int <- 2.5
width_slope <- 0.8
opt_int <- 0.5
opt_slope <- 0.7
height <- height_int + height_slope * truex2 + grpmeans_y
width <- exp(width_int + width_slope * truex2) * -1
opt <- opt_int + opt_slope * truex2
  
antilogit <- function(x){
  exp(x) / (1 + exp(x))
}
alpha <- antilogit(width * (x1 - opt)^2 + height)

# process error on y
obsy <- rbinom(n=length(alpha),size=1,prob=alpha)

# plots simulated "true" and observed data.
obsy_t <- obsy[obsy==1]
x1_t <- x1[obsy==1]
obsx2_t <- obsx2[obsy==1]
obsy_n <- obsy[obsy==0]
x1_n <- x1[obsy==0]
obsx2_n <- obsx2[obsy==0]

par(mfrow=c(1,1))
plot(obsx2_n,x1_n,col="grey80")
points(obsx2_t,x1_t,col="red")


# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, 
               group = photos_groups, n = length(x1),
               ngroups = ngroups)

#### Ordinary logistic regression. ####

# write model
cat("
    model{
    ## Priors
    height_int ~ dunif(-20,20)
    height_slope ~ dunif(-5,5)
    opt_int ~ dunif(-20,20)
    opt_slope ~ dunif(-5,5)
    width_int ~ dunif(-20,20)
    width_slope ~ dnorm(0,0.1)T(-5,5)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)

    ## Likelihood
    for (j in 1:ngroups){
      height.g[j] ~ dnorm(0,group_tau)
    }
    for (i in 1:n){
      height[i] <- height_int + height_slope * x2[i] + height.g[group[i]] 
      opt[i] <- opt_int + opt_slope * x2[i]
      width[i] <- exp(width_int + width_slope * x2[i]) * -1
      logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
      y[i] ~ dbern(alpha[i])
    }
    }
    ",
    fill=TRUE, file="yerror.txt")

require(rjags)
load.module("glm") 

# initiate model
mod1 <- jags.model("yerror.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)

# simulate posterior
out1 <- coda.samples(mod1, n.iter=10000, thin=10,
                    variable.names=c("height_int","height_slope", 
                                     "opt_int","opt_slope", 
                                     "width_int","width_slope",
                                     "height.g[1]"))

# store parameter estimates
require(ggmcmc)
ggd <- ggs(out1)
mod1_height_int <- ggd$value[which(ggd$Parameter == "height_int")]
mod1_height_slope <- ggd$value[which(ggd$Parameter == "height_slope")]
mod1_width_int <- ggd$value[which(ggd$Parameter == "width_int")]
mod1_width_slope <- ggd$value[which(ggd$Parameter == "width_slope")]
mod1_opt_int <- ggd$value[which(ggd$Parameter == "opt_int")]
mod1_opt_slope <- ggd$value[which(ggd$Parameter == "opt_slope")]

d <- data.frame(height_int=mod1_height_int,
                height_slope=mod1_height_slope,
                width_int=mod1_width_int,
                width_slope=mod1_width_slope,
                opt_int=mod1_opt_int,
                opt_slope=mod1_opt_slope)

#### Model with observation error in x.

# specify model
cat("
    model {
    ## Priors
    height_int ~ dunif(-20,20)
    height_slope ~ dunif(-5,5)
    opt_int ~ dunif(-20,20)
    opt_slope ~ dunif(-5,5)
    width_int ~ dunif(-50,20)
    width_slope ~ dnorm(0,0.1)T(-10,10)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)
    err_sd ~ dunif((sd_obs - 0.1),(sd_obs + 0.1))
    tau_obs <- 1 / (err_sd * err_sd)
    
    ## Random effects priors
    for (j in 1:ngroups){
      height.g[j] ~ dnorm(0,group_tau)
    }

    ## Likelihood
    for (i in 1:n){
      x2true[i] ~ dnorm(0,pop_taux)
      x2[i] ~ dnorm(x2true[i],tau_obs)
      height[i] <- height_int + height_slope * x2true[i] + height.g[group[i]] 
      opt[i] <- opt_int + opt_slope * x2true[i]
      width[i] <- exp(width_int + width_slope * x2true[i]) * -1
      logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
      y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="xerror.txt")

# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,group = photos_groups,
               ngroups = ngroups, n = length(x1))

# initiate model
mod2 <- jags.model("xerror.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)

# simulate posterior
out2 <- coda.samples(mod2, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_int","opt_slope", 
                                      "width_int","width_slope",
                                      "height.g[1]","group_sd"))

# store parameter estimates
ggd2 <- ggs(out2)
mod2_height_int <- ggd2$value[which(ggd$Parameter == "height_int")]
mod2_height_slope <- ggd2$value[which(ggd$Parameter == "height_slope")]
mod2_width_int <- ggd2$value[which(ggd$Parameter == "width_int")]
mod2_width_slope <- ggd2$value[which(ggd$Parameter == "width_slope")]
mod2_opt_int <- ggd2$value[which(ggd$Parameter == "opt_int")]
mod2_opt_slope <- ggd2$value[which(ggd$Parameter == "opt_slope")]

d2 <- data.frame(height_int=mod2_height_int,
                height_slope=mod2_height_slope,
                width_int=mod2_width_int,
                width_slope=mod2_width_slope,
                opt_int=mod2_opt_int,
                opt_slope=mod2_opt_slope)

####Plots both models####

ggplot() +
  geom_point(aes(x=obsx, obsy),size=2,color="grey60") +
  geom_abline(aes(intercept=a, slope=b), data=d, color="red", alpha=0.01) +
  geom_abline(aes(intercept=a2, slope=b2), data=d2, color="blue", alpha=0.01) +
  geom_abline(aes(intercept=a3, slope=b3), data=d3, color="green", alpha=0.01) +
  geom_abline(aes(intercept=alpha, slope=beta),
              color="black", size=1.5, linetype="dashed") +
  theme_bw() +
  geom_point(aes(x=truex,y=truey),size=2,color="black")+
  xlab("X") + ylab("Y") +
  guides(color=guide_legend(title="Legend"))+
  ggtitle("Model results with and without modeling error in X and Y")