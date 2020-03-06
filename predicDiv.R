

# Wrapper for prfreq to get divergence for the "Complex" model (see Fig. 1A of Ying et al. 2020)

# Put here the folder to the executable prfreq versions (i.e. prfreq and prfreq_KIRK)
# Download prfreq from here: https://bustamantelab.stanford.edu/lab-developed-software

pathToPrfreq = "~/myBins/"
Sys.setenv(PATH = paste(Sys.getenv()["PATH"], paste0(":", pathToPrfreq), sep=""))


# Return SFS and divergence for neutral single population model

prfreq_returnNeutralSFS = function(sampleSize = 100, nSizeChange = 1, Ne_grid = 2000, fixed=3, poisson = 0, theta = 20000, tdiv = 6.7992, tau = 0.0346, omega = 0.3890, tau_b = 0.01, omega_b = 0.01) {
  
  if (! nSizeChange %in% c(0, 1, 2)) stop("Only one, two or three epochs!")
  if (! poisson %in% c(0,1)) stop("either multinomial or poisson")

  wd = getwd()
  status = TRUE
  rdir = paste("./randdir", floor(runif(1)*1000000), sep="_")
  status = dir.create(rdir)
  
  while (status == FALSE) {
    rdir = paste("./randdir", floor(runif(1)*1000000), sep="_")
    status = dir.create(rdir)
  }
  setwd(rdir)
  
  par = c()
  par = c(par, paste("fixed =", fixed))
  if (poisson == 0) par = c(par, "poisson = 0") else par = c(par, "poisson = 1")
  par = c(par, paste("samplesize = ", sampleSize))
  par = c(par, paste("Ne =", Ne_grid))   # more or less a grid parameter...
  par = c(par, paste("THETA =", theta))
  par = c(par, paste("tdiv =", tdiv))
  if (nSizeChange==0) par = c(par, "num_epochs = 1")
  if (nSizeChange==1) par = c(par, "num_epochs = 2")
  if (nSizeChange==2) par = c(par, "num_epochs = 3")
  if (nSizeChange>0) par = c(par, paste("TAU =", tau, "-1", "-1"))
  if (nSizeChange>0) par = c(par, paste("OMEGA =", omega, "-1", "-1"))
  if (nSizeChange>1) par = c(par, paste("TAU_B =", tau_b, "-1", "-1"))
  if (nSizeChange>1) par = c(par, paste("OMEGA_B =", omega_b, "-1", "-1"))
  par = c(par, "distrib = 0")
  par = c(par, "read_in_snps = 0")
  par = c(par, "run_type = 2")
  par = c(par, "input_grids = 0")
  par = c(par, "output_grids = 0")
  par = c(par, "num_iters = 0")
  
  writeLines(par, con="setting_tmp.txt")
  
  console = system("prfreq setting_tmp.txt temp_out", intern=TRUE, ignore.stderr = TRUE)
  res = readLines(con = "temp_out")
  
  suppressWarnings ({
    div = as.numeric(strsplit(tail(res, 1), split="=")[[1]][2])
    sfs = sapply(strsplit(res[3:(length(res)-1)], split = "="), function(x) as.numeric(x[2]))
    sfs_count = as.numeric(gsub("[^0-9]*", "", sapply(strsplit(res[3:(length(res)-1)], split = "="), function(x) x[1])))
    
    setwd(wd)
    system(paste("rm -r", rdir))
    
    return(list(divergence = div, sfs = data.frame(sfs_count = sfs_count, sfs = sfs), console = console))
  })
}


# Return SFS and divergence for single population model with selection

prfreq_returnSelectionSFS = function(sfs_neutral = sfs_neutral, sampleSize = 100, nSizeChange = 1, Ne_grid = 2000, theta = 20000, tdiv = 6.7992, tau = 0.0346, omega = 0.3890, tau_b = 0.01, omega_b = 0.01, distr = 6, 
                                     P = c(0.25, 1000), I = list(c(-13634658, -1363466), c(-1363466, -136346), c(-136346, -8000), c(-8000, -800), c(-800, -8), c(-8, -0.8), c(-0.8, -0.08), c(-0.08, -0.00001))) {
  
  fixed = 1
  poisson = 0
  if (! nSizeChange %in% c(0, 1, 2)) stop("Only one, two or three epochs!")
  if (! poisson %in% c(0,1)) stop("either multinomial or poisson")
  
  wd = getwd()
  status = TRUE
  rdir = paste("./randdir", floor(runif(1)*1000000), sep="_")
  status = dir.create(rdir)
  
  while (status == FALSE) {
    rdir = paste("./randdir", floor(runif(1)*1000000), sep="_")
    status = dir.create(rdir)
  }
  setwd(rdir)
  
  write.table(sfs_neutral, file="temp_sfs", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  par = c()
  par = c(par, paste("fixed =", fixed))
  if (poisson == 0) par = c(par, "poisson = 0") else par = c(par, "poisson = 1")
  par = c(par, paste("samplesize = ", sampleSize))
  par = c(par, paste("Ne =", Ne_grid))   # more or less a grid parameter...
  par = c(par, paste("THETA =", theta))
  par = c(par, paste("tdiv =", tdiv))
  if (nSizeChange==0) par = c(par, "num_epochs = 1")
  if (nSizeChange==1) par = c(par, "num_epochs = 2")
  if (nSizeChange==2) par = c(par, "num_epochs = 3")
  if (nSizeChange>0) par = c(par, paste("TAU =", tau, "-1", "-1"))
  if (nSizeChange>0) par = c(par, paste("OMEGA =", omega, "-1", "-1"))
  if (nSizeChange>1) par = c(par, paste("TAU_B =", tau_b, "-1", "-1"))
  if (nSizeChange>1) par = c(par, paste("OMEGA_B =", omega_b, "-1", "-1"))
  par = c(par, paste("distrib =", distr))
  for (i in 1:length(P)) par = c(par, paste("P", P[i], "-1", "-1"))
  for (i in 1:length(I)) par = c(par, paste("I 1", paste(I[[i]], collapse=" ")))
  
  par = c(par, "read_in_snps = 0")
  par = c(par, "run_type = 3")
  par = c(par, "input_grids = 0")
  par = c(par, "output_grids = 0")
  par = c(par, "num_iters = 0")
  
  writeLines(par, con="setting_tmp.txt")
  
  console = system("prfreq_KIRK setting_tmp.txt temp_out temp_sfs", intern=TRUE, ignore.stderr = TRUE)
  res = readLines(con = "temp_out")
  
  suppressWarnings ({
    div = as.numeric(strsplit(tail(res, 1), split="=")[[1]][2])
    sfs = sapply(strsplit(res[3:(length(res)-1)], split = "="), function(x) as.numeric(x[2]))
    sfs_count = as.numeric(gsub("[^0-9]*", "", sapply(strsplit(res[3:(length(res)-1)], split = "="), function(x) x[1])))
    
    setwd(wd)
    system(paste("rm -r", rdir))
    
    return(list(divergence = div, sfs = data.frame(sfs_count = sfs_count, sfs = sfs), console = console, output = res))
  })
  
}


sfsFromPrfreq_gammaDFE = function(sample_size, grid_size=2000, theta, alpha, beta) {
  fixed = 1
  poisson = 1
  distr = 6
  IntegrationRanges = list(c(-13634658, -1363466), c(-1363466, -136346), c(-136346, -8000), c(-8000, -800), c(-800, -8), c(-8, -0.8), c(-0.8, -0.08), c(-0.08, -0.00001))
  sfs_neutral = prfreq_returnNeutralSFS(sampleSize = sample_size, nSizeChange = 0, Ne_grid = grid_size, theta = theta, tdiv = 6)$sfs$sfs
  sfs = prfreq_returnSelectionSFS(sfs_neutral = sfs_neutral, sampleSize = sample_size, nSizeChange = 0, Ne_grid = grid_size, theta = theta, tdiv = 6, distr = distr, P = c(alpha, beta), I = IntegrationRanges)$sfs$sfs
  return(sfs)
}


sfsFromPrfreq_singleS = function(sample_size, grid_size=2000, theta, S) {
  fixed = 1
  poisson = 1
  distr = 1
  IntegrationRanges = list(c(-13634658, -1363466), c(-1363466, -136346), c(-136346, -8000), c(-8000, -800), c(-800, -8), c(-8, -0.8), c(-0.8, -0.08), c(-0.08, -0.00001))
  sfs_neutral = prfreq_returnNeutralSFS(sampleSize = sample_size, nSizeChange = 0, Ne_grid = grid_size, theta = theta, tdiv = 6)$sfs$sfs
  sfs = prfreq_returnSelectionSFS(sfs_neutral = sfs_neutral, sampleSize = sample_size, nSizeChange = 0, Ne_grid = grid_size, theta = theta, tdiv = 6, distr = distr, P = c(S), I = IntegrationRanges)$sfs$sfs
  return(sfs)
}

# Computing fixed differences between two species according to Sawyer & Hartl (1992) eqn 13:

Dval = function(tdiv, theta, S) {
  retval_neutral = tdiv*theta
  retval = ifelse (S < -1.0e-7, retval_neutral*2.0*S*exp(2.0*S)/(exp(2.0*S)-1.0), ifelse(S > 1.0e-7, retval_neutral*2.0*S/(1.0-exp(-2.0*S)), retval_neutral))
  return(retval)
}

# Finding divergence by Monte Carlo integration assuming gamma DFE:

Divergence_gammaDFE = function(tdiv, theta, alpha, beta, nrep=1e6) {
  
  gammaDFE_2Nes = -rgamma(nrep, shape = alpha, scale = beta)
  Dval_gammaDFE = Dval(tdiv, theta, gammaDFE_2Nes)
  
  return(sum(Dval_gammaDFE)/nrep)
}

# Function for divergence under the demographic model in Ying et al. 2020, Fig. 1A

# Here, the size change at tPopSizeChange is different from the population size in the ancestral population.
# I.e., we assume that there is an ancestral population size of the ingroup that is different from the population size of the ancestral pop of both ingroup and outgroup

# Note however that this is not modeling the recent population size change (e.g. this relates to very recent increases in humans etc.)
# Repeated estimates might change slightly because of the random nature of the Monte Carlo integration

predicDiv_gammaDFE = function(theta = 20000, tdiv = 6.7992, tau = 0.0346, omega = 0.3890, omega2=NA, omega_ancestral = 1, omega_outgroup = NA, alpha=NA, beta=NA, tPopSizeChange = tdiv) {
  
  nrep = 1e6
  if (is.na(alpha) | is.na(beta)) {
    alpha = 0; beta = 0; nrep = 1
  }
  
  if (alpha == 0 | beta == 0) {
    alpha = 0; beta = 0; nrep = 1
  }
  
  if (is.na(omega_outgroup)) omega_outgroup = omega_ancestral
  
  if (is.na(omega2)) omega2 = omega_ancestral
  
  if (tPopSizeChange < 0 | tPopSizeChange > tdiv) stop("The population size change in the ingroup has to be positive and younger than the divergence time!")
  
  print("Computing ingroup divergence")
  ingroupDivergence1 = Divergence_gammaDFE(tPopSizeChange/omega/2, theta*omega, alpha, beta*omega, nrep=nrep)
  ingroupDivergence2 = Divergence_gammaDFE((tdiv-tPopSizeChange)/omega2/2, theta*omega2, alpha, beta*omega2, nrep=nrep)
  ingroupDivergence = ingroupDivergence1 + ingroupDivergence2
  print("Computing outgroup divergence")
  outgroupDivergence = Divergence_gammaDFE(tdiv/omega_outgroup/2, theta*omega_outgroup, alpha, beta*omega_outgroup, nrep=nrep)
  print("Computing ancestral divergence")
  ancestralDivergence = sum(sfsFromPrfreq_gammaDFE(theta = theta*omega_ancestral, alpha = alpha, beta = beta*omega_ancestral, sample_size = 2))
  
  return(data.table(ingroupDivergence = ingroupDivergence, outgroupDivergence = outgroupDivergence, ancestralDivergence = ancestralDivergence, totalDivergence = ingroupDivergence+outgroupDivergence+ancestralDivergence))
}
  


predicDiv_singleS = function(theta = 20000, tdiv = 6.7992, tau = 0.0346, omega = 0.3890, omega2 = NA, omega_ancestral = 1, omega_outgroup = NA, S=NA, tPopSizeChange = tdiv) {
  
  if (is.na(S)) {
    S = 0
  }
  
  if (is.na(omega_outgroup)) omega_outgroup = omega_ancestral
  
  if (tPopSizeChange < 0 | tPopSizeChange > tdiv) stop("The population size change in the ingroup has to be positive and younger than the divergence time!")
  
  if (is.na(omega2)) omega2 = omega_ancestral
  
  print("Computing ingroup divergence")
  ingroupDivergence1 = Dval(tdiv = tPopSizeChange/omega/2, theta = theta*omega, S = S*omega)
  ingroupDivergence2 = Dval(tdiv = (tdiv-tPopSizeChange)/omega2/2, theta = theta*omega2, S = S*omega2)
  ingroupDivergence = ingroupDivergence1 + ingroupDivergence2
  print("Computing outgroup divergence")
  outgroupDivergence = Dval(tdiv = tdiv/omega_outgroup/2, theta = theta*omega_outgroup, S = S*omega_outgroup)
  print("Computing ancestral divergence")
  ancestralDivergence = sum(sfsFromPrfreq_singleS(theta = theta*omega_ancestral, S = S*omega_ancestral, sample_size = 2))
  
  return(data.table(ingroupDivergence = ingroupDivergence, outgroupDivergence = outgroupDivergence, ancestralDivergence = ancestralDivergence, totalDivergence = ingroupDivergence+outgroupDivergence+ancestralDivergence))
}


# Example:

predicDiv_gammaDFE(theta = 20000, tdiv = 6.7992, tau = 0.0346, omega = 0.3890, omega2=2.3, omega_ancestral = 0.8, omega_outgroup = 1.5, alpha=0.19, beta=1200, tPopSizeChange = 4.5)

# [1] "Computing ingroup divergence"
# [1] "Computing outgroup divergence"
# [1] "Computing ancestral divergence"
# ingroupDivergence outgroupDivergence ancestralDivergence totalDivergence
# 1:          18633.46           16000.19            5293.396        39927.05

