
library(jagsUI)

# DATA

setwd("C:/Users/chris/OneDrive/Documents/Papier3/Data")

CH <- read.csv("CMR_SorGjeslingan.csv", header = FALSE)
BS <- read.csv("BS_SorGjeslingan.csv")
Counts <- read.csv("Counts_SorGjeslingan.csv")

mat <- as.matrix(CH)
colnames(mat) <- NULL
mat <- apply(mat, 2, as.numeric)
CH <- mat
nind <- nrow(CH)
n.occasions <- ncol(CH)

# TEMPORAL ALIGNMENT

yr.cmr <- 2011:(2011 + n.occasions - 1)
yr.bs <- BS$Year
yr.counts <- Counts$Year

idx.prod <- match(yr.bs, yr.cmr)
idx.counts <- match(yr.counts, yr.cmr)

cat("CMR:", range(yr.cmr), "(", n.occasions, "occ)\n")
cat("BS:", range(yr.bs), "(", nrow(BS), "years, gaps 2018-2019)\n")
cat("Counts:", range(yr.counts), "(", nrow(Counts), "years)\n")
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

AP <- 0.85
S0 <- 0.656 * AP
S1 <- 0.937 * AP
FEC <- 0.311
b <- 0.25

M <- matrix(c(
  0,       0,         0,         FEC*S0,
  S1,      0,         0,         0,
  0,       AP*(1-b),  AP*b,      0,
  0,       AP*b,      AP*(1-b),  AP
), 4, 4, byrow = TRUE)

ev <- eigen(M)$vectors[, 1]
ss1 <- abs(ev / sum(ev))

C1 <- Counts$count[1]

# JAGS MODEL

cat(file = "IPM_sorgjeslingan.txt", "
model {

# --- PRIORS ---

mean.phi ~ dunif(0, 1)
mu.phi <- log(mean.phi / (1 - mean.phi))
sigma.phi ~ dunif(0.001, 1)
tau.phi <- pow(sigma.phi, -2)

for (j in 1:2) {
  alpha.p[j] ~ dunif(0, 1)
  beta.p[j] <- log(alpha.p[j] / (1 - alpha.p[j]))
}
sigma.p ~ dunif(0.001, 1)
tau.p <- pow(sigma.p, -2)

mean.fe ~ dunif(0, 3)
mu.fe <- log(mean.fe)
sigma.fe ~ dunif(0, 3)
tau.fe <- pow(sigma.fe, -2)

log.mu.imm ~ dunif(0, 10)
sigma.imm ~ dunif(0.01, 2)
tau.imm <- pow(sigma.imm, -2)

ratio.s0 <- 0.656
ratio.s1 <- 0.937

nb ~ dunif(lower.nb, upper.nb)
ntot <- nb / ss1[4]
N1[1] <- round(ntot * ss1[1])
N2[1] <- round(ntot * ss1[2])
N3[1] <- round(ntot * ss1[3])
N4[1] <- round(nb)

# --- CJS (trap-dependence, no transience) ---

for (t in 1:(T-1)) {
  epsilon.phi[t] ~ dnorm(0, tau.phi)
  logit(phi.a[t]) <- mu.phi + epsilon.phi[t]
  epsilon.p[t] ~ dnorm(0, tau.p)
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

# --- PRODUCTIVITY ---

for (t in 1:T) {
  log.fe[t] ~ dnorm(mu.fe, tau.fe)
  fe[t] <- exp(log.fe[t])
}
for (t in 1:T.prod) {
  J[t] ~ dpois(B[t] * fe[idx.prod[t]])
}

# --- POPULATION ---

for (t in 2:T) {
  log.Nimm[t] ~ dnorm(log.mu.imm, tau.imm)
  Nimm[t] <- exp(log.Nimm[t])

  prod[t-1] <- (fe[t-1] / 2) * delta.s0[t-1]
  delta.s0[t-1] <- ratio.s0 * phi.a[t-1]
  delta.s1[t-1] <- ratio.s1 * phi.a[t-1]
  N1[t] ~ dbin(prod[t-1], N4[t-1])
  N2[t] ~ dbin(delta.s1[t-1], N1[t-1])

  n.prebreed[t-1] <- round(0.75 * N2[t-1] + 0.25 * N3[t-1])
  N3[t] ~ dbin(phi.a[t-1], n.prebreed[t-1])

  n.breed[t-1] <- round(0.25 * N2[t-1] + 0.75 * N3[t-1] + N4[t-1])
  N4.res[t] ~ dbin(phi.a[t-1], n.breed[t-1])
  N4[t] <- N4.res[t] + round(Nimm[t])
}

# --- OBSERVATION (Poisson) ---

for (t in 1:T.counts) {
  C[t] ~ dpois(N4[idx.counts[t]])
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
  C.new[t] ~ dpois(N4[idx.counts[t]])
  sq.C[t] <- pow(C[t] - N4[idx.counts[t]], 2)
  sq.CN[t] <- pow(C.new[t] - N4[idx.counts[t]], 2)
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
  y = CH, f = f, nind = nind, T = n.occasions, m = m, ss1 = ss1,
  z = known.state.cjs(CH),
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
  "N4", "lambda", "omega", "Nimm", "sigma.imm", "log.mu.imm",
  "delta.s0", "delta.s1",
  "fit.C", "fitN.C", "fit.J", "fitN.J"
)

# TEST RUN

cat("--- Test run ---\n")
test <- jags(jags.data, inits, parameters, "IPM_sorgjeslingan.txt",
             n.chains = 2, n.thin = 1, n.iter = 1000,
             n.burnin = 500, n.adapt = 500, parallel = TRUE)
cat("Test OK. Launching full run...\n")

# MCMC

ni <- 300000; nb <- 50000; nc <- 6; nt <- 10; na <- 10000

out <- jags(jags.data, inits, parameters, "IPM_sorgjeslingan.txt",
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
sink("../output/IPM_SorGjeslingan_out.txt")
print(out, digits = 3)
sink()
save(out, file = "../output/IPM_SorGjeslingan.Rdata")
