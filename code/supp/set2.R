args = (commandArgs(TRUE))
eval(parse(text = args[[1]]))
cat("Iteration Number = ", seedN)

source("/home/bc758/research/decals/code/jasa/code/DECALS.R")
########################################################
################### Simulation 2 #######################
########################################################
path="/home/bc758/research/decals/code/jasa/"
load(file=paste0(path,"data/sim2/sig.RData"))
load(file=paste0(path,"data/sim2/frac0.RData"))
load(file=paste0(path,"data/sim2/sigma_ols0.RData"))

sim_data2 = data_gen2(frac0,sig,sigma_ols0,seedN=seedN)

bulk=sim_data2$bulk
sig=sim_data2$sig
CTS_proportion=constraint_ols(sig=sig,bulk=bulk)
lambda_ols <- tune_select(y_all=bulk,W=t(sig),propor_est=t(CTS_proportion),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")
sigma_ols0<-sigma_upt0(propor_all=t(CTS_proportion),y_all=bulk,W=t(sig),lambda_all=lambda_ols[[1]])

CTS_proportion_var=propor_cov0(W=t(sig),sigma_all=sigma_ols0,propor_all=t(CTS_proportion))
dim(CTS_proportion_var) # CTS_proportions_var * samples
decals_res=list(CTS_proportion=CTS_proportion,CTS_proportion_var=CTS_proportion_var)


save(sim_data2,file=paste0("/home/bc758/research/decals/result/sim/sim2/sim_data2_",seedN,".RData"))
save(decals_res,file=paste0("/home/bc758/research/decals/result/sim/sim2/decals_res_",seedN,".RData"))

