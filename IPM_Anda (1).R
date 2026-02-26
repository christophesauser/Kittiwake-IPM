
library(jagsUI)

# DATA

setwd("C:/Users/christophe.sauser/Documents/Paper3/R")

CH <- read.csv("CMR_Anda.csv", header = FALSE)
BS <- read.csv("BS_Anda.csv")
Counts <- read.csv("Counts_Anda.csv")

mat <- as.matrix(CH)
colnames(mat) <- NULL
mat <- apply(mat, 2, as.numeric)
CH <- mat
nind <- nrow(CH)
n.occasions <- ncol(CH)

# TEMPORAL ALIGNMENT

yr.cmr <- 2005:(2005 + n.occasions - 1)
yr.bs <- BS$Year
yr.counts <- Counts$Year

idx.prod <- match(yr.bs, yr.cmr)
idx.counts <- match(yr.counts, yr.cmr)

cat("CMR:", range(yr.cmr), "| BS:", range(yr.bs), "| Counts:", range(yr.counts), "\n")
cat("idx.prod:", idx.prod, "\n")
cat("idx.counts:", idx.counts, "\n")

# CMR PREPARATION

get.first <- function(x) min(which(x != 0))
f <- apply(CH, 1, get.first)

known.state.cjs <- function(ch) {
  state <- ch
  for (i in 1:nrow(ch)) {
    n1 <- min(which(ch[i,] == 1))
    n2 <- max(which(ch[i,] == 1))
    state[i, n1:n2] <- 1
    state[i, n1] <- NA
  }
  state[state == 0] <- NA
  return(state)
}

cjs.init.z <- function(ch, f) {
  for (i in 1:nrow(ch)) {
    if (sum(ch[i,]) == 1) next
    n2 <- max(which(ch[i,] == 1))
    ch[i, f[i]:n2] <- NA
  }
  for (i in 1:nrow(ch)) {
    ch[i, 1:f[i]] <- NA
  }
  return(ch)
}

m <- CH[, 1:(n.occasions - 1)]
m[m == 0] <- 2

# INITIAL POPULATION STRUCTURE

AP <- 0.842
JP <- 0.6
FEC <- 0.548
b <- 0.25

M <- matrix(c(
  0,       0,         0,         FEC*JP,
  JP,      0,         0,         0,
  0,       AP*(1-b),  AP*b,      0,
  0,       AP*b,      AP*(1-b),  AP
), 4, 4, byrow = TRUE)

ev <- eigen(M)$vectors[, 1]
ss1 <- abs(ev / sum(ev))

C1 <- Counts$count[1]

# JAGS MODEL

cat(file = "IPM_anda.txt", "
model {

# --- POPULATION STATE-SPACE ---

for (t in 2:T) {
  log.Nimm[t] ~ dnorm(log.mu.imm, tau.imm)
  Nimm[t] <- exp(log.Nimm[t])
}
log.mu.imm ~ dunif(0, 10)
sigma.imm ~ dunif(0.01, 2)
tau.imm <- pow(sigma.imm, -2)

delta.sj <- 0.6

nb ~ dunif(lower.nb, upper.nb)
ntot <- nb / ss1[4]
N1[1] <- round(ntot * ss1[1])
N2[1] <- round(ntot * ss1[2])
N3[1] <- round(ntot * ss1[3])
N4[1] <- round(nb)

for (t in 2:T) {
  prod[t-1] <- (fe[t-1] / 2) * delta.sj
  N1[t] ~ dbin(prod[t-1], N4[t-1])
  N2[t] ~ dbin(delta.sj, N1[t-1])

  n.prebreed[t-1] <- round(0.75 * N2[t-1] + 0.25 * N3[t-1])
  N3[t] ~ dbin(phi.a[t-1], n.prebreed[t-1])

  n.breed[t-1] <- round(0.25 * N2[t-1] + 0.75 * N3[t-1] + N4[t-1])
  N4.res[t] ~ dbin(phi.a[t-1], n.breed[t-1])
  N4[t] <- N4.res[t] + round(Nimm[t])
}

sigma.obs ~ dunif(0.001, 0.5)
tau.obs <- pow(sigma.obs, -2)
for (t in 1:T.counts) {
  Ntot[t] <- N4[idx.counts[t]]
  log.Ntot[t] <- log(max(Ntot[t], 1))
  C[t] ~ dlnorm(log.Ntot[t], tau.obs)
}

# --- PRODUCTIVITY ---

mean.fe ~ dunif(0, 3)
mu.fe <- log(mean.fe)
sigma.fe ~ dunif(0, 3)
tau.fe <- pow(sigma.fe, -2)
for (t in 1:T) {
  log.fe[t] ~ dnorm(mu.fe, tau.fe)
  fe[t] <- exp(log.fe[t])
}
for (t in 1:T.prod) {
  J[t] ~ dpois(B[t] * fe[idx.prod[t]])
}

# --- CJS (no transience, trap-dependence) ---

mean.phi ~ dunif(0, 1)
mu.phi <- log(mean.phi / (1 - mean.phi))
for (t in 1:(T-1)) {
  epsilon.phi[t] ~ dnorm(0, tau.phi)
  logit(phi.a[t]) <- mu.phi + epsilon.phi[t]
}
tau.phi <- pow(sigma.phi, -2)
sigma.phi ~ dunif(0.001, 1)

for (t in 1:(T-1)) {
  epsilon.p[t] ~ dnorm(0, tau.p)
}
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0.001, 1)

for (j in 1:2) {
  alpha.p[j] ~ dunif(0, 1)
  beta.p[j] <- log(alpha.p[j] / (1 - alpha.p[j]))
}

for (i in 1:nind) {
  for (t in f[i]:(T-1)) {
    logit(phi[i,t]) <- mu.phi + epsilon.phi[t]
    logit(p[i,t]) <- beta.p[m[i,t]] + epsilon.p[t]
  }
}

for (i in 1:nind) {
  z[i, f[i]] <- 1
  for (t in (f[i]+1):T) {
    z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
    y[i,t] ~ dbern(p[i,t-1] * z[i,t])
  }
}

# --- DERIVED ---

for (t in 1:(T-1)) {
  lambda[t] <- N4[t+1] / N4[t]
}
for (t in 2:T) {
  omega[t] <- round(Nimm[t]) / N4[t]
}
for (j in 1:2) {
  rec.p[j] <- 1 / (1 + exp(-beta.p[j]))
}

# --- GOF ---

for (t in 1:T.counts) {
  C.new[t] ~ dlnorm(log.Ntot[t], tau.obs)
  sq.C[t] <- pow(log(C[t]) - log.Ntot[t], 2)
  sq.CN[t] <- pow(log(C.new[t]) - log.Ntot[t], 2)
}
fit.C <- sum(sq.C[])
fitN.C <- sum(sq.CN[])

for (t in 1:T.prod) {
  J.exp[t] <- B[t] * fe[idx.prod[t]]
  J.new[t] ~ dpois(J.exp[t])
  chi2.J[t] <- pow(J[t] - J.exp[t], 2) / (J.exp[t] + 0.5)
  chi2.JN[t] <- pow(J.new[t] - J.exp[t], 2) / (J.exp[t] + 0.5)
}
fit.J <- sum(chi2.J[])
fitN.J <- sum(chi2.JN[])

}
")

# BUNDLE DATA

jags.data <- list(
  y = CH, f = f, nind = nind, T = n.occasions,
  z = known.state.cjs(CH), m = m, ss1 = ss1,
  B = BS$AON, J = BS$chicks,
  C = as.numeric(Counts$count),
  T.prod = nrow(BS), T.counts = nrow(Counts),
  idx.prod = idx.prod, idx.counts = idx.counts,
  lower.nb = C1 * 0.95, upper.nb = C1 * 1.05
)

inits <- function() { list(z = cjs.init.z(CH, f)) }

parameters <- c(
  "phi.a", "mean.phi", "sigma.phi",
  "rec.p",
  "mean.fe", "fe", "sigma.fe",
  "N4", "lambda", "omega", "Nimm", "sigma.obs", "sigma.imm", "log.mu.imm",
  "fit.C", "fitN.C", "fit.J", "fitN.J"
)

# TEST RUN

cat("--- Test run ---\n")
test <- jags(jags.data, inits, parameters, "IPM_anda.txt",
             n.chains = 2, n.thin = 1, n.iter = 1000,
             n.burnin = 500, n.adapt = 500, parallel = TRUE)
cat("Test OK. Launching full run...\n")

# MCMC

ni <- 100000; nb <- 20000; nc <- 3; nt <- 5; na <- 5000

out <- jags(jags.data, inits, parameters, "IPM_anda.txt",
            n.chains = nc, n.thin = nt, n.iter = ni,
            n.burnin = nb, n.adapt = na, parallel = TRUE)

print(out, digits = 3)

# DIAGNOSTICS

max.rhat <- max(unlist(out$Rhat), na.rm = TRUE)
min.neff <- min(unlist(out$n.eff), na.rm = TRUE)
cat("\nMax Rhat:", max.rhat, "| Min n.eff:", min.neff, "\n")
if (max.rhat > 1.1) cat("CONVERGENCE ISSUES - consider longer run\n")

pval.C <- mean(out$sims.list$fit.C > out$sims.list$fitN.C)
pval.J <- mean(out$sims.list$fit.J > out$sims.list$fitN.J)
cat("Bayesian p-value (Counts):", pval.C, "\n")
cat("Bayesian p-value (Productivity):", pval.J, "\n")

# SAVE

options(max.print = 1e5)
sink("IPM_Anda_out.txt")
print(out, digits = 3)
sink()
save(out, file = "IPM_Anda.Rdata")
