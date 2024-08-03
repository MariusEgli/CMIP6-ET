library(ncdf4)
library(stringr)
library(data.table)
source("/net/xenon/climphys_backedup/maegli/CMIP6/data/process_CMIP6_fun.R")

mon_dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/"


CMIP6.files_pr       = data.table(get.CMIP6.file.list(vari= "pr", 
                                           temp.res = "mon", 
                                           scen = "hist585", 
                                           CMIP6.dir = paste0(mon_dir, "pr/hist585/")))[str_detect(file.name, "mmday")]


CMIP6.files_hfls = data.table(get.CMIP6.file.list(vari= "hfls", 
                                                      temp.res = "mon", 
                                                      scen = "hist585", 
                                                      CMIP6.dir = paste0(mon_dir, "hfls/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_tas = data.table(get.CMIP6.file.list(vari= "tas", 
                                                  temp.res = "mon", 
                                                  scen = "hist585", 
                                                  CMIP6.dir = paste0(mon_dir, "tas/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_huss = data.table(get.CMIP6.file.list(vari= "huss", 
                                                 temp.res = "mon", 
                                                 scen = "hist585", 
                                                 CMIP6.dir = paste0(mon_dir, "huss/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_psl = data.table(get.CMIP6.file.list(vari= "psl", 
                                                  temp.res = "mon", 
                                                  scen = "hist585", 
                                                  CMIP6.dir = paste0(mon_dir, "psl/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_rnet = data.table(get.CMIP6.file.list(vari= "rsds", 
                                                 temp.res = "mon", 
                                                 scen = "hist585", 
                                                 CMIP6.dir = paste0(mon_dir, "rnet/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_netlw_hist585   = data.table(get.CMIP6.file.list(vari= "netlw", 
                                                             temp.res = "mon", 
                                                             scen = c("hist585"), 
                                                             CMIP6.dir = paste0(mon_dir, "netlw/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_netsw_hist585   = data.table(get.CMIP6.file.list(vari= "netsw", 
                                                             temp.res = "mon", 
                                                             scen = c("hist585"), 
                                                             CMIP6.dir = paste0(mon_dir, "netsw/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_hurs_hist585   = data.table(get.CMIP6.file.list(vari= "hurs", 
                                                             temp.res = "mon", 
                                                             scen = c("hist585"), 
                                                             CMIP6.dir = paste0(mon_dir, "hurs/hist585/")))[!str_detect(file.name, "anom")]

CMIP6.files_dtot_hist585 = data.table(get.CMIP6.file.list(vari= "dtot", 
                                                          temp.res = "mon", 
                                                          scen = "hist585", 
                                                          CMIP6.dir = paste0(mon_dir, "dtot/hist585/")))[!str_detect(file.name, "anom")]


CMIP6.files_rsds_hist585 = data.table(get.CMIP6.file.list(vari= "rsds", 
                                                          temp.res = "mon", 
                                                          scen = "hist585", 
                                                          CMIP6.dir = paste0(mon_dir, "rsds/hist585/")))[!str_detect(file.name, "anom")]
CMIP6.files_rlds_hist585 = data.table(get.CMIP6.file.list(vari= "rlds", 
                                                          temp.res = "mon", 
                                                          scen = "hist585", 
                                                          CMIP6.dir = paste0(mon_dir, "rlds/hist585/")))[!str_detect(file.name, "anom")]




modall_intersect <- intersect(CMIP6.files_pr$modall, CMIP6.files_hfls$modall) %>% 
  intersect(CMIP6.files_tas$modall) %>% 
  intersect(CMIP6.files_huss$modall) %>% 
  intersect(CMIP6.files_psl$modall) %>% 
  intersect(CMIP6.files_rnet$modall)

modall_intersect <- modall_intersect[!modall_intersect %in% c("MIROC6_hist585_r29i1p1f1", "MIROC6_hist585_r30i1p1f1")]


mod_table <- CMIP6.files_hfls[modall %in% modall_intersect, table(mod)] 
models <- c(names(mod_table[mod_table >=5]), "FGOALS-g3", "NorESM2-MM", "HadGEM3-GC31-LL")


CMIP6.files_hfls <- CMIP6.files_hfls[modall %in% modall_intersect][mod %in% models]
CMIP6.files_huss <- CMIP6.files_huss[modall %in% modall_intersect][mod %in% models]
CMIP6.files_tas <- CMIP6.files_tas[modall %in% modall_intersect][mod %in% models]
CMIP6.files_pr <- CMIP6.files_pr[modall %in% modall_intersect][mod %in% models]
CMIP6.files_psl <- CMIP6.files_psl[modall %in% modall_intersect][mod %in% models]
CMIP6.files_rnet <- CMIP6.files_rnet[modall %in% modall_intersect][mod %in% models]
CMIP6.files_netlw <- CMIP6.files_netlw_hist585[modall %in% modall_intersect][mod %in% models]
CMIP6.files_netsw <- CMIP6.files_netsw_hist585[modall %in% modall_intersect][mod %in% models]
CMIP6.files_hurs <- CMIP6.files_hurs_hist585[modall %in% modall_intersect][mod %in% models]
CMIP6.files_dtot <- CMIP6.files_dtot_hist585[modall %in% modall_intersect][mod %in% models]
CMIP6.files_rsds <- CMIP6.files_rsds_hist585[modall %in% modall_intersect][mod %in% models]
CMIP6.files_rlds <- CMIP6.files_rlds_hist585[modall %in% modall_intersect][mod %in% models]



CMIP6.files_dtot_hist585_anom <- CMIP6.files_dtot[, lapply(.SD, function(x)  str_replace(x, "_g025", "anom_g025"))]
CMIP6.files_dtot_hist585_JJA_anom <- CMIP6.files_dtot_hist585_anom[, lapply(.SD, function(x)  str_replace(x, "mon", "JJA"))]

CMIP6.files_rsds_hist585_anom <- CMIP6.files_rsds[, lapply(.SD, function(x)  str_replace(x, "2d50", "anom_2d50"))]
CMIP6.files_rsds_hist585_JJA <- CMIP6.files_rsds[, lapply(.SD, function(x)  str_replace(x, "mon ", "JJA"))]
CMIP6.files_rsds_hist585_JJA_anom <- CMIP6.files_rsds_hist585_anom[, lapply(.SD, function(x)  str_replace(x, "mon", "JJA"))]

CMIP6.files_rlds_hist585_anom <- CMIP6.files_rlds[, lapply(.SD, function(x)  str_replace(x, "2d50", "anom_2d50"))]
CMIP6.files_rlds_hist585_JJA <- CMIP6.files_rlds[, lapply(.SD, function(x)  str_replace(x, "mon ", "JJA"))]
CMIP6.files_rlds_hist585_JJA_anom <- CMIP6.files_rlds_hist585_anom[, lapply(.SD, function(x)  str_replace(x, "mon", "JJA"))]


#hfls
in_names  <- sapply(CMIP6.files_hfls[,file.name], FUN = function(x) paste0(mon_dir,"hfls/hist585/",x), USE.NAMES = F)
out_names  <- sapply(CMIP6.files_hfls[,file.name], FUN = function(x) paste0(mon_dir,"hfls/seas/",x), USE.NAMES = F) %>% str_replace("mon", "JJA")
anom_names <- str_replace(out_names, "2d50", "anom_2d50")


for(i in 1:length(in_names)){
  system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", out_names[i]))
  system(paste0("cdo -ymonsub ", out_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", out_names[i], " ", anom_names[i]))
  }

#huss
in_names  <- sapply(CMIP6.files_huss[,file.name], FUN = function(x) paste0(mon_dir,"huss/hist585/",x), USE.NAMES = F)
out_names  <- sapply(CMIP6.files_huss[,file.name], FUN = function(x) paste0(mon_dir,"huss/seas/",x), USE.NAMES = F) %>% str_replace("mon", "JJA")
anom_names <- str_replace(out_names, "2d50", "anom_2d50")

for(i in 1:length(in_names)){
  #system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", out_names[i]))
  system(paste0("cdo -ymonsub ", out_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", out_names[i], " ", anom_names[i]))
  
}

#tas
in_names  <- sapply(CMIP6.files_tas[,file.name], FUN = function(x) paste0(mon_dir,"tas/hist585/",x), USE.NAMES = F)
out_names  <- sapply(CMIP6.files_tas[,file.name], FUN = function(x) paste0(mon_dir,"tas/seas/",x), USE.NAMES = F) %>% str_replace("mon", "JJA")
anom_names <- str_replace(out_names, "2d50", "anom_2d50")


for(i in 1:length(in_names)){
  system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", str_replace(out_names[i], "mon", "JJA")))
  system(paste0("cdo -ymonsub ", out_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", out_names[i], " ", anom_names[i]))
}

#psl
in_names  <- sapply(CMIP6.files_psl[,file.name], FUN = function(x) paste0(mon_dir,"psl/hist585/",x), USE.NAMES = F)
out_names  <- sapply(CMIP6.files_psl[,file.name], FUN = function(x) paste0(mon_dir,"psl/seas/",x), USE.NAMES = F) %>% str_replace("mon", "JJA")
anom_names <- str_replace(out_names, "2d50", "anom_2d50")


for(i in 1:length(in_names)){
  #system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", str_replace(out_names[i], "mon", "JJA")))
  system(paste0("cdo -ymonsub ", out_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", out_names[i], " ", anom_names[i]))
  
}

#rnet
in_names  <- sapply(CMIP6.files_rnet[,file.name], FUN = function(x) paste0(mon_dir,"rnet/hist585/",x), USE.NAMES = F)
out_names  <- sapply(CMIP6.files_rnet[,file.name], FUN = function(x) paste0(mon_dir,"rnet/seas/",x), USE.NAMES = F) %>% str_replace("mon", "JJA")
anom_names <- str_replace(out_names, "g025", "anom_2d50")

for(i in 1:length(in_names)){
  system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", str_replace(out_names[i], "mon", "JJA")))
  system(paste0("cdo -ymonsub ", out_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", out_names[i], " ", anom_names[i]))
  
}

#pr
in_names  <- sapply(CMIP6.files_pr[,file.name], FUN = function(x) paste0(mon_dir,"pr/hist585/",x), USE.NAMES = F)
out_names  <- sapply(CMIP6.files_pr[,file.name], FUN = function(x) paste0(mon_dir,"pr/seas/",x), USE.NAMES = F) %>% str_replace("mon", "JJA")
anom_names <- str_replace(out_names, "2d50", "anom_2d50")

for(i in 1:length(in_names)){
  #system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", str_replace(out_names[i], "mon", "JJA")))
  system(paste0("cdo -ymonsub ", out_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", out_names[i], " ", anom_names[i]))
  
}

#netlw
in_names  <- sapply(CMIP6.files_netlw[,file.name], FUN = function(x) paste0(mon_dir,"netlw/hist585/",x), USE.NAMES = F)
anom_names <- str_replace(in_names, "g025", "anom_2d50")
seasanom_names  <- sapply(CMIP6.files_netlw[,file.name], FUN = function(x) paste0(mon_dir,"netlw/seas/",x), USE.NAMES = F) %>% 
  str_replace("mon", "JJA") %>% 
  str_replace("g025", "anom_2d50")
seas_names  <- sapply(CMIP6.files_netlw[,file.name], FUN = function(x) paste0(mon_dir,"netlw/seas/",x), USE.NAMES = F) %>% 
  str_replace("mon", "JJA") 
for(i in 1:length(in_names)){
  system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", seas_names[i]))
  system(paste0("cdo -ymonsub ", in_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", in_names[i], " ", anom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", anom_names[i], " ", seasanom_names[i]))
  
}

#netsw
in_names  <- sapply(CMIP6.files_netsw[,file.name], FUN = function(x) paste0(mon_dir,"netsw/hist585/",x), USE.NAMES = F)
anom_names <- str_replace(in_names, "g025", "anom_2d50")
seasanom_names  <- sapply(CMIP6.files_netsw[,file.name], FUN = function(x) paste0(mon_dir,"netsw/seas/",x), USE.NAMES = F) %>% 
  str_replace("mon", "JJA") %>% 
  str_replace("g025", "anom_2d50")
seas_names  <- sapply(CMIP6.files_netsw[,file.name], FUN = function(x) paste0(mon_dir,"netsw/seas/",x), USE.NAMES = F) %>% 
  str_replace("mon", "JJA") 
for(i in 1:length(in_names)){
  system(paste0("cdo -ymonsub ", in_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", in_names[i], " ", anom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", anom_names[i], " ", seasanom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", seas_names[i]))
  
}

#hurs
in_names  <- sapply(CMIP6.files_hurs[,file.name], FUN = function(x) paste0(mon_dir,"hurs/hist585/",x), USE.NAMES = F)
anom_names <- str_replace(in_names, "2d50", "anom_2d50")
seasanom_names  <- sapply(CMIP6.files_hurs[,file.name], FUN = function(x) paste0(mon_dir,"hurs/seas/",x), USE.NAMES = F) %>% 
  str_replace("mon", "JJA") %>% 
  str_replace("2d50", "anom_2d50")
seas_names  <- sapply(CMIP6.files_hurs[,file.name], FUN = function(x) paste0(mon_dir,"hurs/seas/",x), USE.NAMES = F) %>% 
  str_replace("mon", "JJA") 

for(i in 1:length(in_names)){
  system(paste0("cdo -ymonsub ", in_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", in_names[i], " ", anom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", anom_names[i], " ", seasanom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", in_names[i], " ", seas_names[i]))
}



#dtot
for(i in 1:nrow(CMIP6.files_dtot_hist585_anom)){
  infile <- paste0(mon_dir, "dtot/hist585/",  CMIP6.files_dtot_hist585_anom$file.name[i])
  outfile <- paste0(mon_dir, "dtot/seas/",  CMIP6.files_dtot_hist585_JJA_anom$file.name[i])
  #system(paste0("cdo -ymonsub ", in_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", in_names[i], " ", anom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", infile, " ", outfile))
}

#rsds
for(i in 1:nrow(CMIP6.files_rsds_hist585_anom)){
  infile <- paste0(mon_dir, "rsds/hist585/",  CMIP6.files_rsds_hist585$file.name[i])
  outfile <- paste0(mon_dir, "rsds/seas/",  CMIP6.files_rsds_hist585_JJA$file.name[i])
  #system(paste0("cdo -ymonsub ", in_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", in_names[i], " ", anom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", infile, " ", outfile))
}

#rlds
for(i in 1:nrow(CMIP6.files_rlds_hist585_anom)){
  infile <- paste0(mon_dir, "rlds/hist585/",  CMIP6.files_rlds_hist585$file.name[i])
  outfile <- paste0(mon_dir, "rlds/seas/",  CMIP6.files_rlds_hist585_JJA$file.name[i])
  #system(paste0("cdo -ymonsub ", in_names[i], " -ymonmean -seldate,1950-01-01,1999-12-31 ", in_names[i], " ", anom_names[i]))
  system(paste0("cdo -seasmean -select,season=JJA ", infile, " ", outfile))
}


CMIP6.files_pr[, file.name := str_replace(file.name, "mon", "MAMJJA")]
CMIP6.files_hfls[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_huss[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_tas[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_psl[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_rnet[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_netlw[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_netsw[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_hurs[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_rlds[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_rsds[, file.name := str_replace(file.name, "mon", "JJA")]

#pr
cmip6_pr_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_pr$file.name, 
                                                  scen = CMIP6.files_pr$scen, var="pr", res = "ann", 
                                                  CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/pr/seas", time.count = rep(-1, 1579))

CMIP6_pr_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_pr_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_pr_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_pr, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/pr/seas/", mons = c(6,7,8))

save(list = "CMIP6_pr_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_pr_JJA_hist585_2d50_XAX.RData")

rm(cmip6_pr_hist585_JJA_2d50)


#hfls
cmip6_hfls_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hfls$file.name, 
                                             scen = CMIP6.files_hfls$scen, var="hfls", res = "ann", 
                                             CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/seas", time.count = rep(-1, 1579))

CMIP6_hfls_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_hfls_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_hfls_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_hfls, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/seas/", mons = c(6,7,8))

save(list = "CMIP6_hfls_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_hfls_JJA_hist585_2d50_XAX.RData")
rm(cmip6_hfls_hist585_JJA_2d50)

#tas
cmip6_tas_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_tas$file.name, 
                                               scen = CMIP6.files_tas$scen, var="tas", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/tas/seas", time.count = rep(-1, 1579))

CMIP6_tas_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_tas_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_tas_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_tas, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/tas/seas/", mons = c(6,7,8))

save(list = "CMIP6_tas_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_tas_JJA_hist585_2d50_XAX.RData")

rm(cmip6_tas_hist585_JJA_2d50)

#huss
cmip6_huss_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_huss$file.name, 
                                              scen = CMIP6.files_huss$scen, var="huss", res = "ann", 
                                              CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/huss/seas", time.count = rep(-1, 1579))

CMIP6_huss_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_huss_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_huss_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_huss, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/huss/seas/", mons = c(6,7,8))

save(list = "CMIP6_huss_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_huss_JJA_hist585_2d50_XAX.RData")

rm(cmip6_huss_hist585_JJA_2d50)

#psl
cmip6_psl_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_psl$file.name, 
                                               scen = CMIP6.files_psl$scen, var="psl", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/psl/seas", time.count = rep(-1, 1579))

CMIP6_psl_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_psl_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_psl_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_psl, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/psl/seas/", mons = c(6,7,8))

save(list = "CMIP6_psl_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_psl_JJA_hist585_2d50_XAX.RData")

rm(cmip6_psl_hist585_JJA_2d50)
gc()


#rnet
cmip6_rnet_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_rnet$file.name, 
                                              scen = CMIP6.files_rnet$scen, var="rsds", res = "ann", 
                                              CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/seas", time.count = rep(-1, 1579))

CMIP6_rnet_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_rnet_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_rnet_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_rnet, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/seas/", mons = c(6,7,8))

save(list = "CMIP6_rnet_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_rnet_JJA_hist585_2d50_XAX.RData")

rm(cmip6_rnet_hist585_JJA_2d50)
gc()

#netsw
cmip6_netsw_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_netsw$file.name, 
                                               scen = CMIP6.files_netsw$scen, var="rsds", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netsw/seas", time.count = rep(-1, 1579))

CMIP6_netsw_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_netsw_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_netsw_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_netsw, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netsw/seas/", mons = c(6,7,8))

save(list = "CMIP6_netsw_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_netsw_JJA_hist585_2d50_XAX.RData")

rm(cmip6_netsw_hist585_JJA_2d50)
gc()

#netlw
cmip6_netlw_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_netlw$file.name, 
                                               scen = CMIP6.files_netlw$scen, var="rlds", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netlw/seas", time.count = rep(-1, 1579))

CMIP6_netlw_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_netlw_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_netlw_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_netlw, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netlw/seas/", mons = c(6,7,8))

save(list = "CMIP6_netlw_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_netlw_JJA_hist585_2d50_XAX.RData")

rm(cmip6_netlw_hist585_JJA_2d50)
gc()


#hurs
cmip6_hurs_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hurs$file.name, 
                                                scen = CMIP6.files_hurs$scen, var="hurs", res = "ann", 
                                                CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hurs/seas", time.count = rep(-1, 1579))

CMIP6_hurs_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_hurs_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_hurs_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_hurs, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hurs/seas/", mons = c(6,7,8))

save(list = "CMIP6_hurs_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_hurs_JJA_hist585_2d50_XAX.RData")

rm(cmip6_hurs_hist585_JJA_2d50)
gc()
#===========anom=================

CMIP6.files_pr[, file.name := str_replace(file.name, "2d50", "anom_2d50")]
CMIP6.files_hfls[, file.name := str_replace(file.name, "2d50", "anom_2d50")]
CMIP6.files_huss[, file.name := str_replace(file.name, "2d50", "anom_2d50")]
CMIP6.files_tas[, file.name := str_replace(file.name, "2d50", "anom_2d50")]
CMIP6.files_psl[, file.name := str_replace(file.name, "2d50", "anom_2d50")]
CMIP6.files_rnet[, file.name := str_replace(file.name, "g025", "anom_2d50")]
CMIP6.files_netlw[, file.name := str_replace(file.name, "g025", "anom_2d50")]
CMIP6.files_netsw[, file.name := str_replace(file.name, "g025", "anom_2d50")]
CMIP6.files_hurs[, file.name := str_replace(file.name, "2d50", "anom_2d50")]
CMIP6.files_rlds[, file.name := str_replace(file.name, "2d50", "anom_2d50")]
CMIP6.files_rsds[, file.name := str_replace(file.name, "2d50", "anom_2d50")]


CMIP6.files_pr_MAMJJA_anom <-  CMIP6.files_pr[, file.name := str_replace(file.name, "2d50", "anom_2d50")][, file.name := str_replace(CMIP6.files_pr$file.name, "_mmday", "")]
#pr
cmip6_pr_hist585_MAMJJA_anom_2d50 = read.CMIP6_novar(file.name = str_replace(CMIP6.files_pr$file.name, "_mmday", ""), 
                                             scen = CMIP6.files_pr$scen, var="pr", res = "ann", 
                                             CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/pr/seas", time.count = rep(-1, 1579))

CMIP6_pr_MAMJJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_pr_hist585_MAMJJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_pr_MAMJJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_pr, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/pr/seas/", mons = c(6,7,8))

save(list = "CMIP6_pr_MAMJJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_pr_MAMJJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_pr_hist585_JJA_anom_2d50)


#hfls
cmip6_hfls_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hfls$file.name, 
                                               scen = CMIP6.files_hfls$scen, var="hfls", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/seas", time.count = rep(-1, 1579))

CMIP6_hfls_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_hfls_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_hfls, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/seas/", mons = c(6,7,8))

save(list = "CMIP6_hfls_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hfls_JJA_hist585_anom_2d50_XAX.RData")
rm(cmip6_hfls_hist585_JJA_anom_2d50)

#tas
cmip6_tas_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_tas$file.name, 
                                              scen = CMIP6.files_tas$scen, var="tas", res = "ann", 
                                              CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/tas/seas", time.count = rep(-1, 1579))

CMIP6_tas_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_tas_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_tas_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_tas, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/tas/seas/", mons = c(6,7,8))

save(list = "CMIP6_tas_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_tas_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_tas_hist585_JJA_anom_2d50)

#huss
cmip6_huss_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_huss$file.name, 
                                               scen = CMIP6.files_huss$scen, var="huss", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/huss/seas", time.count = rep(-1, 1579))

CMIP6_huss_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_huss_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_huss_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_huss, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/huss/seas/", mons = c(6,7,8))

save(list = "CMIP6_huss_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_huss_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_huss_hist585_JJA_anom_2d50)

#psl
cmip6_psl_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_psl$file.name, 
                                              scen = CMIP6.files_psl$scen, var="psl", res = "ann", 
                                              CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/psl/seas", time.count = rep(-1, 1579))

CMIP6_psl_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_psl_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_psl_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_psl, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/psl/seas/", mons = c(6,7,8))

save(list = "CMIP6_psl_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_psl_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_psl_hist585_JJA_anom_2d50)
gc()


#rnet
cmip6_rnet_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_rnet$file.name, 
                                               scen = CMIP6.files_rnet$scen, var="rsds", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/seas", time.count = rep(-1, 1579))

CMIP6_rnet_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_rnet_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_rnet, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/seas/", mons = c(6,7,8))

save(list = "CMIP6_rnet_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_rnet_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_rnet_hist585_JJA_anom_2d50)
gc()


#netlw
cmip6_netlw_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_netlw$file.name, 
                                                    scen = CMIP6.files_netlw$scen, var="rlds", res = "ann", 
                                                    CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netlw/seas", time.count = rep(-1, 1579))

CMIP6_netlw_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_netlw_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_netlw_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_netlw, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netlw/seas/", mons = c(6,7,8))

save(list = "CMIP6_netlw_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netlw_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_netlw_hist585_JJA_anom_2d50)
gc()

#netsw
cmip6_netsw_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_netsw$file.name, 
                                                     scen = CMIP6.files_netsw$scen, var="rsds", res = "ann", 
                                                     CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netsw/seas", time.count = rep(-1, 1579))

CMIP6_netsw_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_netsw_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_netsw_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_netsw, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/netsw/seas/", mons = c(6,7,8))

save(list = "CMIP6_netsw_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netsw_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_netsw_hist585_JJA_anom_2d50)
gc()


#hurs
cmip6_hurs_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hurs$file.name, 
                                                     scen = CMIP6.files_hurs$scen, var="hurs", res = "ann", 
                                                     CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hurs/seas", time.count = rep(-1, 1579))

CMIP6_hurs_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_hurs_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_hurs_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_hurs, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hurs/seas/", mons = c(6,7,8))

save(list = "CMIP6_hurs_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hurs_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_hurs_hist585_JJA_anom_2d50)
gc()

#dtot
cmip6_dtot_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_dtot_hist585_JJA_anom$file.name, 
                                                    scen = CMIP6.files_dtot_hist585_JJA_anom$scen, var="rlds", res = "ann", 
                                                    CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/dtot/seas", time.count = rep(-1, 1579))

CMIP6_dtot_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_dtot_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_dtot_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_dtot_hist585_JJA_anom, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/dtot/seas/", mons = c(6,7,8))

save(list = "CMIP6_dtot_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_dtot_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_dtot_hist585_JJA_anom_2d50)
gc()


#rsds
cmip6_rsds_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_rsds_hist585_JJA_anom$file.name, 
                                                    scen = CMIP6.files_rsds_hist585_JJA_anom$scen, var="rsds", res = "ann", 
                                                    CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rsds/seas", time.count = rep(-1, 1579))

CMIP6_rsds_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_rsds_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_rsds_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_rsds_hist585_JJA_anom, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rsds/seas/", mons = c(6,7,8))

save(list = "CMIP6_rsds_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_rsds_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_rsds_hist585_JJA_anom_2d50)
gc()

#rlds
cmip6_rlds_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_rlds_hist585_JJA_anom$file.name, 
                                                    scen = CMIP6.files_rlds_hist585_JJA_anom$scen, var="rlds", res = "ann", 
                                                    CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rlds/seas", time.count = rep(-1, 1579))

CMIP6_rlds_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_rlds_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_rlds_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_rlds_hist585_JJA_anom, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rlds/seas/", mons = c(6,7,8))

save(list = "CMIP6_rlds_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_rlds_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_rlds_hist585_JJA_anom_2d50)
gc()
#==================eq1============


lambda  <- 2.5008*10^6 #J kg^-1
cp <- 1005 #J kg^−1 K^−1)
Rv <- 461.5 #(J kg−1 K−1) 

cp*Rv / lambda^2

CMIP6.files_huss[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_tas[, file.name := str_replace(file.name, "mon", "JJA")]
CMIP6.files_q_T <- CMIP6.files_tas[, lapply(.SD, function(x) str_replace(x,"tas", "q_T"))]

for(i in 1:nrow(CMIP6.files_tas)){
system(paste0("cdo -O merge /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/huss/seas/", CMIP6.files_huss$file.name[i], " ",
       " /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/tas/seas/", CMIP6.files_tas$file.name[i], " ",
       "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", CMIP6.files_q_T$file.name[i]))
}

CMIP6.files_B <- CMIP6.files_tas[, lapply(.SD, function(x) str_replace(x,"tas", "B"))]
CMIP6.files_B_anom <- CMIP6.files_B[,.(file.name = str_replace(file.name, "2d50", "anom_2d50"))]

for(i in 1:nrow(CMIP6.files_tas)){
  system(paste0("cdo exprf,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/mccoll_eq1 ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", CMIP6.files_q_T$file.name[i], " ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", CMIP6.files_B$file.name[i]))
  }


  for(i in 1:nrow(CMIP6.files_tas)){
  
  system(paste0("cdo -selvar,B -ymonsub /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", CMIP6.files_B$file.name[i],
                " -ymonmean -seldate,1950-01-01,1999-12-31 /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", CMIP6.files_B$file.name[i],
                " /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/",CMIP6.files_B_anom$file.name[i]))
}


eq1_files <- CMIP6.files_tas[, lapply(.SD, function(x) str_replace(x,"tas", "eq1"))]



for(i in 1:nrow(CMIP6.files_tas)){
  system(paste0("cdo -O merge /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/seas/", CMIP6.files_rnet$file.name[i], " ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", CMIP6.files_B$file.name[i], " ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", eq1_files$file.name[i]))

}


CMIP6.files_collET <- CMIP6.files_tas[, lapply(.SD, function(x) str_replace(x,"tas", "collET"))][, lapply(.SD, function(x) str_replace(x,"mon", "JJA"))]
CMIP6.files_collET_anom <- CMIP6.files_collET[, lapply(.SD, function(x) str_replace(x,"mon", "JJA"))]
CMIP6.files_collET_anom[, file.name := str_replace(file.name, "2d50", "anom_2d50")]

for(i in 1:nrow(CMIP6.files_tas)){
  system(paste0("cdo -selvar,ET -exprf,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/bowen_to_et ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", eq1_files$file.name[i], " ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", CMIP6.files_collET$file.name[i]))
}

for(i in 1:nrow(CMIP6.files_tas)){
  system(paste0("cdo -ymonsub ", "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", CMIP6.files_collET$file.name[i],
                " -ymonmean -seldate,1950-01-01,1999-12-31 /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", CMIP6.files_collET$file.name[i]," ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", CMIP6.files_collET_anom$file.name[i]))
  

}


#McColl ET XAX
cmip6_MCET_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_collET$file.name, 
                                                   scen = CMIP6.files_collET$scen, var="ET", res = "ann", 
                                                   CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll", time.count = rep(-1, 1579))

CMIP6_MCET_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_MCET_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_MCET_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_collET, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", mons = c(6,7,8))

save(list = "CMIP6_MCET_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_MCET_JJA_hist585_2d50_XAX.RData")

rm(cmip6_MCET_hist585_JJA_2d50)
gc()

#McColl ET ANOM XAX
cmip6_MCET_hist585_anom_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_collET_anom$file.name, 
                                               scen = CMIP6.files_collET_anom$scen, var="ET", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll", time.count = rep(-1, 1579))

CMIP6_MCET_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_MCET_hist585_anom_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_MCET_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_collET_anom, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", mons = c(6,7,8))

save(list = "CMIP6_MCET_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_MCET_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_MCET_hist585_anom_JJA_2d50)
gc()


#McColl B XAX anom 
cmip6_B_hist585_JJA_2d50_anom = read.CMIP6_novar(file.name = CMIP6.files_B_anom$file.name, 
                                               scen = "hist585", var="B", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas", time.count = rep(-1, 1579))

CMIP6_B_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_B_hist585_JJA_2d50_anom ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_B_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_B_anom, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", mons = c(6,7,8))

save(list = "CMIP6_B_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_B_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_B_hist585_JJA_2d50_anom)
gc()


cmip6_B_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_B$file.name, 
                                            scen = "hist585", var="B", res = "ann", 
                                            CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas", time.count = rep(-1, 1579))

CMIP6_B_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_B_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_B_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_B, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", mons = c(6,7,8))

save(list = "CMIP6_B_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_B_JJA_hist585_2d50_XAX.RData")

rm(cmip6_B_hist585_JJA_2d50)
gc()



#===========R net as LH + SH============
setwd("/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs")

CMIP6.files_hfls
CMIP6.files_hfss <- CMIP6.files_hfls[, lapply(.SD, function(x) str_replace(x,"hfls", "hfss"))]
CMIP6.files_LHSH <- CMIP6.files_hfls[, lapply(.SD, function(x) str_replace(x,"hfls", "LHSH"))]
CMIP6.files_LHSH_JJA <- CMIP6.files_LHSH[, lapply(.SD, function(x) str_replace(x,"mon", "JJA"))]


target_dir <- "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/LHSH/"

CMIP6.files_hfss <- CMIP6.files_hfls[, lapply(.SD, function(x) str_replace(x,"hfls", "hfss"))][, lapply(.SD, function(x) str_replace(x,"JJA", "mon"))]
CMIP6.files_hfss_anom <- CMIP6.files_hfss[, lapply(.SD, function(x) str_replace(x,"2d50", "anom_2d50"))]
CMIP6.files_hfss_anom_JJA <-  CMIP6.files_hfss_anom[, lapply(.SD, function(x) str_replace(x,"mon", "JJA"))]
CMIP6.files_hfss_JJA <-  CMIP6.files_hfss[, lapply(.SD, function(x) str_replace(x,"mon", "JJA"))]


for(i in 1:nrow(CMIP6.files_hfss)){
  system(paste0("cdo seasmean -select,season=JJA /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/hist585/", CMIP6.files_hfss_anom$file.name[i],
                              " /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/seas/", CMIP6.files_hfss_anom_JJA$file.name[i]))
  
}

for(i in 1:nrow(CMIP6.files_hfss)){
  system(paste0("cdo seasmean -select,season=JJA /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/hist585/", CMIP6.files_hfss$file.name[i],
                " /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/seas/", CMIP6.files_hfss_JJA$file.name[i]))
  
}


for(i in 1:nrow(CMIP6.files_hfls)){
  system(paste0("cdo add ./hfls/orig/", CMIP6.files_hfls$file.name[i], " ./hfss/hist585/", CMIP6.files_hfss$file.name[i], " ./rnet/LHSH/", CMIP6.files_LHSH$file.name[i]))
  system(paste0("cdo seasmean -select,season=JJA /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/LHSH/", CMIP6.files_LHSH$file.name[i],
                " /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/LHSH/", CMIP6.files_LHSH_JJA$file.name[i]))

  
}

eq1_files_LHSH <- CMIP6.files_LHSH_JJA[, lapply(.SD, function(x) str_replace(x,"LHSH", "eq1_LHSH"))]

for(i in 1:nrow(CMIP6.files_tas)){
  system(paste0("cdo -O merge /net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/rnet/LHSH/", CMIP6.files_LHSH_JJA$file.name[i], " ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/q_T/seas/", CMIP6.files_B$file.name[i], " ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", eq1_files_LHSH$file.name[i]))
}


CMIP6.files_collSHLH <- CMIP6.files_LHSH_JJA[, lapply(.SD, function(x) str_replace(x,"LHSH", "collSHLH"))]



for(i in 1:nrow(CMIP6.files_tas)){
  system(paste0("cdo -selvar,ET -exprf,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/bowen_to_et_LHSH ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", eq1_files_LHSH$file.name[i], " ",
                "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", CMIP6.files_collSHLH$file.name[i]))
}




cmip6_MCLHSH_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_collSHLH$file.name, 
                                               scen = CMIP6.files_collSHLH$scen, var="ET", res = "ann", 
                                               CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll", time.count = rep(-1, 1579))

CMIP6_MCLHSH_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_MCLHSH_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_MCLHSH_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_collSHLH, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/mccoll/", mons = c(6,7,8))

save(list = "CMIP6_MCLHSH_JJA_hist585_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_MCLHSH_JJA_hist585_2d50_XAX.RData")

rm(cmip6_MCLHSH_hist585_JJA_2d50)
gc()

##================HFSS=================

cmip6_hfss_hist585_JJA_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hfss_JJA$file.name, 
                                                    scen = CMIP6.files_hfss_JJA$scen, var="hfss", res = "ann", 
                                                    CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/seas", time.count = rep(-1, 1579))

CMIP6_hfss_JJA_hist585_2d50_XAX = XAX_mon(X = cmip6_hfss_hist585_JJA_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_hfss_JJA_hist585_2d50_XAX$M <- get_MAX(CMIP6.files_hfss_JJA, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/seas/", mons = c(6,7,8))

save(list = "CMIP6_hfss_JJA_hist585_2d50_XAX",
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_hfss_JJA_hist585_2d50_XAX.RData")

rm(cmip6_hfss_hist585_JJA_2d50)
gc()





cmip6_hfss_hist585_JJA_anom_2d50 = read.CMIP6_novar(file.name = CMIP6.files_hfss_anom_JJA$file.name, 
                                                 scen = CMIP6.files_hfss_anom_JJA$scen, var="hfss", res = "ann", 
                                                 CMIP5.dir = "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/seas", time.count = rep(-1, 1579))

CMIP6_hfss_JJA_hist585_anom_2d50_XAX = XAX_mon(X = cmip6_hfss_hist585_JJA_anom_2d50 ,rm.HIST.from.rcp26 = F,  ncol = 144*72)
CMIP6_hfss_JJA_hist585_anom_2d50_XAX$M <- get_MAX(CMIP6.files_hfss_anom_JJA, "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfss/seas/", mons = c(6,7,8))

save(list = "CMIP6_hfss_JJA_hist585_anom_2d50_XAX", 
     file = "/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hfss_JJA_hist585_anom_2d50_XAX.RData")

rm(cmip6_hfss_hist585_JJA_anom_2d50)
gc()




#===========SW down====================




