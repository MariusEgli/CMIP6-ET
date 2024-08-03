
library(ncdf4)
library(data.table)
library(magrittr)
source("/net/h2o/climphys1/sippels/_projects/_small_projects/aer_vs_ghg_LE/code/_functions_CMIP5_extr.R")
source("/net/h2o/climphys1/sippels/_projects/_small_projects/aer_vs_ghg_LE/code/_functions_CMIP6.R")
source("/net/h2o/climphys1/maegli/CMIP6/Pr_vs_ET_functions.R")


#====================pr585===========================

CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025"
target_dir_5d00_monthly = "/net/h2o/climphys1/maegli/CMIP6/data/NetCDFs/"
CMIP6.files_pr       = get.CMIP6.file.list(vari= "pr", 
                                           temp.res = "mon", 
                                           scen = c("ssp585"), 
                                           CMIP6.dir = CMIP6.orig_dir)


CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025"
CMIP6.dir = "/net/ch4/data/cmip6-Next_Generation/"
CMIP6.files_hfls       = get.CMIP6.file.list(vari= "hfls", 
                                             temp.res = "mon", 
                                             scen = c("ssp585"), 
                                             CMIP6.dir = CMIP6.orig_dir)



CMIP6_LEs <- c("ACCESS-ESM1-5","CanESM5","EC-Earth3","MIROC-ES2L","MIROC6","MPI-ESM1-2-LR")

#modall_int <- intersect(CMIP6.files_pr$modall, CMIP6.files_hfls$modall)
CMIP6.files_pr <- CMIP6.files_pr[which(CMIP6.files_pr$mod %in% CMIP6_LEs),]
CMIP6.files_hfls <- CMIP6.files_hfls[which(CMIP6.files_hfls$mod %in% CMIP6_LEs),]
CMIP6.files_psl <-  CMIP6.files_psl[which(CMIP6.files_psl$mod %in% CMIP6_LEs),]


modall_int <- intersect(CMIP6.files_pr$modall, CMIP6.files_hfls$modall)
CMIP6.files_pr <- CMIP6.files_pr[which(CMIP6.files_pr$modall %in% modall_int),]
CMIP6.files_hfls <- CMIP6.files_hfls[which(CMIP6.files_hfls$modall %in% modall_int),]


meta_DT_pr <- stringr::str_split_fixed(CMIP6.files_pr$file.name, "_", n = 6) %>% data.table() %>% 
  setnames( c("vari", "res", "mod", "scen", "ens.mem", "rest")) %>% 
  .[, c("file.name", "modall") :=  list(CMIP6.files_pr$file.name,paste0(mod, "_",scen , "_", ens.mem))]

cmip6_pr_585_mon_2d50 = read.CMIP6_novar(file.name = CMIP6.files_pr$file.name, 
                                                  scen = CMIP6.files_pr$scen, var="pr", res = "mon", 
                                                  CMIP5.dir = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025/", time.count = rep(-1, 1579))

CMIP6_pr_mon_585_2d50_XAX = get.XAX_mon(X = cmip6_pr_585_mon_2d50, M = meta_DT_pr,  rm.HIST.from.rcp26 = F, ncol = 144*72)
CMIP6_pr_mon_585_2d50_XAX$X = CMIP6_pr_mon_585_2d50_XAX$X* (3600*24*30)
time_pr <- lapply(CMIP6.files_pr$file.name, get_proper_time, wd = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025/") %>% rbindlist()
CMIP6_pr_mon_585_2d50_XAX$M$year <- time_pr$year
CMIP6_pr_mon_585_2d50_XAX$M$mon <- time_pr$mon


save(list = "CMIP6_pr_mon_585_2d50_XAX", 
     file = "/net/h2o/climphys1/maegli/CMIP6/data/XAX/CMIP6_pr_mon_585_2d50_XAX.RData")



#=============hfsl585==================


meta_DT_hfls <- stringr::str_split_fixed(CMIP6.files_hfls$file.name, "_", n = 6) %>% data.table() %>% 
  setnames( c("vari", "res", "mod", "scen", "ens.mem", "rest")) %>% 
  .[, c("file.name", "modall") :=  list(CMIP6.files_hfls$file.name,paste0(mod, "_",scen , "_", ens.mem))]


cmip6_hfls_585_mon_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hfls$file.name, 
                                         scen = CMIP6.files_hfls$scen, var="hfls", res = "mon", 
                                         CMIP5.dir = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025/", time.count = rep(-1, 1579))


CMIP6_hfls_mon_585_2d50_XAX = get.XAX_mon(X = cmip6_hfls_585_mon_2d50, M = meta_DT_hfls,  rm.HIST.from.rcp26 = F, ncol = 144*72)
time_hfls <- lapply(CMIP6.files_hfls$file.name, get_proper_time, wd = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025/") %>% rbindlist()
CMIP6_hfls_mon_585_2d50_XAX$X <- CMIP6_hfls_mon_585_2d50_XAX$X * (3600*24*30 /(2.5*10^6))
CMIP6_hfls_mon_585_2d50_XAX$M$year <- time_hfls$year
CMIP6_hfls_mon_585_2d50_XAX$M$mon <- time_hfls$mon

save(list = "CMIP6_hfls_mon_585_2d50_XAX", 
     file = "/net/h2o/climphys1/maegli/CMIP6/data/XAX/CMIP6_hfls_mon_585_2d50_XAX.RData")
#=================PiC===============================


CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025"
target_dir_5d00_monthly = "/net/h2o/climphys1/maegli/CMIP6/data/NetCDFs/"
CMIP6.files_pr_PiC       = get.CMIP6.file.list(vari= "pr", 
                                           temp.res = "mon", 
                                           scen = "piControl", 
                                           CMIP6.dir = CMIP6.orig_dir)

CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025"
CMIP6.dir = "/net/ch4/data/cmip6-Next_Generation/"
CMIP6.files_hfls_PiC       = get.CMIP6.file.list(vari= "hfls", 
                                             temp.res = "mon", 
                                             scen = "piControl", 
                                             CMIP6.dir = CMIP6.orig_dir)

#modall_int_PiC <- intersect(CMIP6.files_pr_PiC$modall, CMIP6.files_hfls_PiC$modall)
CMIP6.files_pr_PiC <- CMIP6.files_pr_PiC[CMIP6.files_pr_PiC$mod %in% CMIP6_LEs & CMIP6.files_pr_PiC$ens.mem %in% c("r1i1p1f1","r1i1p1f2"),]
CMIP6.files_hfls_PiC <- CMIP6.files_hfls_PiC[which(CMIP6.files_hfls_PiC$mod %in% CMIP6_LEs) & CMIP6.files_hfls_PiC$ens.mem %in% c("r1i1p1f1","r1i1p1f2"),]


meta_DT_pr_PiC <- stringr::str_split_fixed(CMIP6.files_pr_PiC$file.name, "_", n = 6) %>% data.table() %>% 
  setnames( c("vari", "res", "mod", "scen", "ens.mem", "rest")) %>% 
  .[, c("file.name", "modall") :=  list(CMIP6.files_pr_PiC$file.name,paste0(mod, "_",scen , "_", ens.mem))]



cmip6_pr_PiC_mon_2d50 = read.CMIP6_novar(file.name = CMIP6.files_pr_PiC$file.name, 
                                         scen = CMIP6.files_pr_PiC$scen, var="pr", res = "mon", 
                                         CMIP5.dir = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025/", time.count = rep(-1, 24000),
                                         rm.years = 0)

CMIP6_pr_mon_PiC_2d50_XAX = get.XAX_mon(X = cmip6_pr_PiC_mon_2d50, M = meta_DT_pr_PiC,  rm.HIST.from.rcp26 = F, ncol = 144*72)
CMIP6_pr_mon_PiC_2d50_XAX$X = CMIP6_pr_mon_PiC_2d50_XAX$X *24*3600*30
time_pr <- lapply(CMIP6.files_pr_PiC$file.name, get_proper_time, wd = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025/") %>% rbindlist()
CMIP6_pr_mon_PiC_2d50_XAX$M$year <- time_pr$year
CMIP6_pr_mon_PiC_2d50_XAX$M$mon <- time_pr$mon

save(list = "CMIP6_pr_mon_PiC_2d50_XAX", 
     file = "/net/h2o/climphys1/maegli/CMIP6/data/XAX/CMIP6_pr_mon_PiC_2d50_XAX.RData")





meta_DT_hfls_PiC <- stringr::str_split_fixed(CMIP6.files_hfls_PiC$file.name, "_", n = 6) %>% data.table() %>% 
  setnames( c("vari", "res", "mod", "scen", "ens.mem", "rest")) %>% 
  .[, c("file.name", "modall") :=  list(CMIP6.files_hfls_PiC$file.name,paste0(mod, "_",scen , "_", ens.mem))]

cmip6_hfls_PiC_mon_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hfls_PiC$file.name, 
                                         scen = CMIP6.files_hfls_PiC$scen, var="hfls", res = "mon", 
                                         CMIP5.dir = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025/", time.count = rep(-1, 24000),
                                         rm.years = 0)

CMIP6_hfls_mon_PiC_2d50_XAX = get.XAX_mon(X = cmip6_hfls_PiC_mon_2d50, M = meta_DT_hfls_PiC,  rm.HIST.from.rcp26 = F, ncol = 144*72)
CMIP6_hfls_mon_PiC_2d50_XAX$X = CMIP6_hfls_mon_PiC_2d50_XAX$X * 3600*24*30 /(2.5*10^6)
time_pr <- lapply(CMIP6.files_hfls_PiC$file.name, get_proper_time, wd = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025/") %>% rbindlist()
CMIP6_hfls_mon_PiC_2d50_XAX$M$year <- time_pr$year
CMIP6_hfls_mon_PiC_2d50_XAX$M$mon <- time_pr$mon

save(list = "CMIP6_hfls_mon_PiC_2d50_XAX", 
     file = "/net/h2o/climphys1/maegli/CMIP6/data/XAX/CMIP6_hfls_mon_PiC_2d50_XAX.RData")


CMIP6_hfls_mon_PiC_2d50_XAX$X[1,] %>% mean()
CMIP6_hfls_mon_585_2d50_XAX$X[1,] %>% mean()


rbind(data.table(var = "ET", scen = "PiC", value = CMIP6_hfls_mon_PiC_2d50_XAX$X[1,]),
data.table(var = "ET", scen = "SSP585", value =CMIP6_hfls_mon_585_2d50_XAX$X[1,] ),
data.table(var = "Pr" , scen = "PiC", value = CMIP6_pr_mon_PiC_2d50_XAX$X[1,] ),
data.table(var = "Pr" , scen = "SSP585", value =  CMIP6_pr_mon_585_2d50_XAX$X[1,] )) %>% 
  ggplot()+
  geom_density(aes(value, group = interaction(scen, var), fill = interaction(var,scen)), alpha = 0.5)



#===============hist===========================

CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025"
target_dir_5d00_monthly = "/net/h2o/climphys1/maegli/CMIP6/data/NetCDFs/"
CMIP6.files_hist_pr       = get.CMIP6.file.list(vari= "pr", 
                                           temp.res = "mon", 
                                           scen = "historical", 
                                           CMIP6.dir = CMIP6.orig_dir)


CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025"
CMIP6.dir = "/net/ch4/data/cmip6-Next_Generation/"
CMIP6.files_hist_hfls       = get.CMIP6.file.list(vari= "hfls", 
                                             temp.res = "mon", 
                                             scen = "historical", 
                                             CMIP6.dir = CMIP6.orig_dir)




CMIP6_LEs <- c("ACCESS-ESM1-5","CanESM5","EC-Earth3","MIROC-ES2L","MIROC6","MPI-ESM1-2-LR")

#modall_int <- intersect(CMIP6.files_pr$modall, CMIP6.files_hfls$modall)
CMIP6.files_hist_pr <- CMIP6.files_hist_pr[which(CMIP6.files_hist_pr$mod %in% CMIP6_LEs),]
CMIP6.files_hist_hfls <- CMIP6.files_hist_hfls[which(CMIP6.files_hist_hfls$mod %in% CMIP6_LEs),]

modall_int <- intersect(CMIP6.files_hist_pr$modall, CMIP6.files_hist_hfls$modall)
CMIP6.files_hist_pr <- CMIP6.files_hist_pr[which(CMIP6.files_hist_pr$modall %in% modall_int),]
CMIP6.files_hist_hfls <- CMIP6.files_hist_hfls[which(CMIP6.files_hist_hfls$modall %in% modall_int),]


meta_DT_pr_hist <- stringr::str_split_fixed(CMIP6.files_hist_pr$file.name, "_", n = 6) %>% data.table() %>% 
  setnames( c("vari", "res", "mod", "scen", "ens.mem", "rest")) %>% 
  .[, c("file.name", "modall") :=  list(CMIP6.files_hist_pr$file.name,paste0(mod, "_",scen , "_", ens.mem))]

cmip6_pr_hist_mon_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hist_pr$file.name, 
                                         scen = CMIP6.files_hist_pr$scen, var="pr", res = "mon", 
                                         CMIP5.dir = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025/", time.count = rep(-1, 1579))


CMIP6_pr_mon_hist_2d50_XAX = get.XAX_mon(X = cmip6_pr_hist_mon_2d50, M = meta_DT_pr_hist,  rm.HIST.from.rcp26 = F, ncol = 144*72)
CMIP6_pr_mon_hist_2d50_XAX$X = CMIP6_pr_mon_hist_2d50_XAX$X * 3600*24*30
time_pr <- lapply(CMIP6.files_hist_pr$file.name, get_proper_time, wd = "/net/ch4/data/cmip6-Next_Generation/pr/mon/g025/") %>% rbindlist()
CMIP6_pr_mon_hist_2d50_XAX$M$year <- time_pr$year
CMIP6_pr_mon_hist_2d50_XAX$M$mon <- time_pr$mon

save(list = "CMIP6_pr_mon_hist_2d50_XAX", 
     file = "/net/h2o/climphys1/maegli/CMIP6/data/XAX/CMIP6_pr_mon_hist_2d50_XAX.RData")



meta_DT_hfls_hist <- stringr::str_split_fixed(CMIP6.files_hist_hfls$file.name, "_", n = 6) %>% data.table() %>% 
  setnames( c("vari", "res", "mod", "scen", "ens.mem", "rest")) %>% 
  .[, c("file.name", "modall") :=  list(CMIP6.files_hist_hfls$file.name,paste0(mod, "_",scen , "_", ens.mem))]

cmip6_hfls_hist_mon_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hist_hfls$file.name, 
                                          scen = CMIP6.files_hist_hfls$scen, var="hfls", res = "mon", 
                                          CMIP5.dir = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025/", time.count = rep(-1, 1579))

CMIP6_hfls_mon_hist_2d50_XAX = get.XAX_mon(X = cmip6_hfls_hist_mon_2d50, M = meta_DT_hfls_hist,  rm.HIST.from.rcp26 = F, ncol = 144*72)
CMIP6_hfls_mon_hist_2d50_XAX$X = CMIP6_hfls_mon_hist_2d50_XAX$X * 3600*24*30 /(2.5*10^6)
time_pr <- lapply(CMIP6.files_hist_hfls$file.name, get_proper_time, wd = "/net/ch4/data/cmip6-Next_Generation/hfls/mon/g025/") %>% rbindlist()
CMIP6_hfls_mon_hist_2d50_XAX$M$year <- time_pr$year
CMIP6_hfls_mon_hist_2d50_XAX$M$mon <- time_pr$mon

save(list = "CMIP6_hfls_mon_hist_2d50_XAX", 
     file = "/net/h2o/climphys1/maegli/CMIP6/data/XAX/CMIP6_hfls_mon_hist_2d50_XAX.RData")






