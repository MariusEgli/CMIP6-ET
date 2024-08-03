library(ncdf4)

for (year in 1950:2021){
  #T
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
               "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.t2m.025deg.1m.",year,".nc ",
               "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/T/orig/era5_deterministic_recent.t2m.2.5deg.1m",year,".nc"
               ))
  
  #q
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
                 "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.sp.025deg.1m.",year,".nc ",
                 "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/Q/orig/era5_deterministic_recent.sp.2.5deg.1m",year,".nc"
  ))
  
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
                "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.d2m.025deg.1m.",year,".nc ",
                "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/Q/orig/era5_deterministic_recent.d2m.2.5deg.1m",year,".nc"
  ))
  
  #SLP
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
                 "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.msl.025deg.1m.",year,".nc ",
                 "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/SLP/orig/era5_deterministic_recent.msl.2.5deg.1m",year,".nc"
  ))
}

for (year in 1950:2021){
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
                "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.e.025deg.1m.",year,".nc ",
                "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/ET/orig/era5_deterministic_recent.e.2.5deg.1m",year,".nc"
  ))
  }

for (year in 1950:2021){
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
                "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.tp.025deg.1m.",year,".nc ",
                "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/P/orig/era5_deterministic_recent.tp.2.5deg.1m",year,".nc"
  ))
}

for (year in 1950:2021){
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
                "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.msnswrf.025deg.1m.",year,".nc ",
                "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/sw_orig/era5_deterministic_recent.msnswrf.2.5deg.1m",year,".nc"
  ))
}

for (year in 1950:2021){
  system(paste0("cdo remapcon,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ",
                "/net/atmos/data/ERA5_deterministic/recent/0.25deg_lat-lon_1m/original/era5_deterministic_recent.msnlwrf.025deg.1m.",year,".nc ",
                "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/lw_orig/era5_deterministic_recent.msnlwrf.2.5deg.1m",year,".nc"
  ))
}





system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/T/orig/*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021.nc")
system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/Q/orig/era5_deterministic_recent.sp.*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/Q/era5_deterministic_recent.sp.2.5deg.1m.1950-2021.nc")
system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/Q/orig/era5_deterministic_recent.d2m.*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/Q/era5_deterministic_recent.d2m.2.5deg.1m.1950-2021.nc")
system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/SLP/orig/era5_deterministic_recent.msl*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/SLP/era5_deterministic_recent.msl.2.5deg.1m.1950-2021.nc")
system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/ET/orig/era5_deterministic_recent.e*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc")
system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/P/orig/era5_deterministic_recent.tp*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/P/era5_deterministic_recent.tp.2.5deg.1m.1950-2021.nc")

system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/sw_orig/era5_deterministic_recent.msnswrf*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/era5_deterministic_recent.msnswrf.2.5deg.1m.1950-2021.nc")
system("cdo -b F64 mergetime /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/lw_orig/era5_deterministic_recent.msnlwrf*.nc /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/era5_deterministic_recent.msnlwrf.2.5deg.1m.1950-2021.nc")

system(paste0("cdo add /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/era5_deterministic_recent.msnlwrf.2.5deg.1m.1950-2021.nc ",
              "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/era5_deterministic_recent.msnswrf.2.5deg.1m.1950-2021.nc ",
              "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/era5_deterministic_recent.rnet.2.5deg.1m.1950-2021.nc"))


system(paste0("cdo sub /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/P/era5_deterministic_recent.tp.2.5deg.1m.1950-2021.nc ",
       "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc ",
       "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/P/era5_deterministic_recent.PmE.2.5deg.1m.1950-2021.nc"))

system(paste0("cdo mulc,1000 /net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/P/era5_deterministic_recent.PmE.2.5deg.1m.1950-2021.nc ",
              "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/P/era5_deterministic_recent.PmE_mm.2.5deg.1m.1950-2021.nc"))


setwd("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/Q")

#e = 6.112 * exp((17.67 * Td) / (Td + 243.5))
system("cdo exprf,dewpoint_to_e era5_deterministic_recent.d2m.2.5deg.1m.1950-2021.nc era5_deterministic_recent.eandstuff.2.5deg.1m.1950-2021.nc")
system("cdo select,name=e era5_deterministic_recent.eandstuff.2.5deg.1m.1950-2021.nc era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc")

#q = 0.622 * e / (P - e)
system("cdo mulc,0.622 era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc 0.622e.nc")
system("cdo sub era5_deterministic_recent.sp.2.5deg.1m.1950-2021.nc era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc SPminusvaporpressure.nc")
system("cdo div 0.622e.nc SPminusvaporpressure.nc era5_deterministic_recent.q.2.5deg.1m.1950-2021.nc")

#RH = 100 × {exp[17.625 × Dp/(243.04 + Dp)]/exp[17.625 × T/(243.04 + T)]}
setwd("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5")


system("cdo merge ./Q/era5_deterministic_recent.d2m.2.5deg.1m.1950-2021.nc ./T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021.nc ./RH/orig/era5_deterministic_recent.d_t2m.2.5deg.1m.1950-2021.nc")
system("cdo exprf,./RH/orig/d2mtoRH ./RH/orig/era5_deterministic_recent.d_t2m.2.5deg.1m.1950-2021.nc ./RH/orig/era5_deterministic_recent.rh.2.5deg.1m.1950-2021.nc")

#anomaly calculation
system("cdo -ymonsub ./T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021.nc ./T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021_anom.nc")
system("cdo -ymonsub ./Q/era5_deterministic_recent.q.2.5deg.1m.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./Q/era5_deterministic_recent.q.2.5deg.1m.1950-2021.nc ./Q/era5_deterministic_recent.q.2.5deg.1m.1950-2021_anom.nc")
system("cdo -ymonsub ./SLP/era5_deterministic_recent.msl.2.5deg.1m.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./SLP/era5_deterministic_recent.msl.2.5deg.1m.1950-2021.nc ./SLP/era5_deterministic_recent.msl.2.5deg.1m.1950-2021_anom.nc")
system("cdo -ymonsub ./ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc ./ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021_anom.nc")
system("cdo -ymonsub ./R/era5_deterministic_recent.rnet.2.5deg.1m.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./R/era5_deterministic_recent.rnet.2.5deg.1m.1950-2021.nc ./R/era5_deterministic_recent.rnet.2.5deg.1m.1950-2021_anom.nc")
system("cdo -ymonsub ./P/era5_deterministic_recent.PmE_mm.2.5deg.1m.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./P/era5_deterministic_recent.PmE_mm.2.5deg.1m.1950-2021.nc ./P/era5_deterministic_recent.PmE_mm.2.5deg.1m.1950-2021_anom.nc")
system("cdo -ymonsub ./RH/orig/era5_deterministic_recent.rh.2.5deg.1m.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./RH/orig/era5_deterministic_recent.rh.2.5deg.1m.1950-2021.nc ./RH/era5_deterministic_recent.rh.2.5deg.1m.1950-2021_anom.nc")


#seas mean
system("cdo -seasmean -select,season=JJA ./T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021.nc ./T/era5_deterministic_recent.t2m.2.5deg.JJA.1950-2021.nc")
system("cdo -seasmean -select,season=JJA ./Q/era5_deterministic_recent.q.2.5deg.1m.1950-2021.nc ./Q/era5_deterministic_recent.q.2.5deg.JJA.1950-2021.nc")
system("cdo -seasmean -select,season=JJA ./SLP/era5_deterministic_recent.msl.2.5deg.1m.1950-2021.nc ./SLP/era5_deterministic_recent.msl.2.5deg.JJA.1950-2021.nc")
system("cdo -seasmean -select,season=JJA ./R/era5_deterministic_recent.rnet.2.5deg.1m.1950-2021.nc ./R/era5_deterministic_recent.rnet.2.5deg.JJA.1950-2021.nc")
system("cdo -seasmean -select,season=JJA ./P/era5_deterministic_recent.PmE.2.5deg.1m.1950-2021.nc ./P/era5_deterministic_recent.PmE.2.5deg.JJA.1950-2021.nc")
system("cdo -seasmean -select,season=JJA ./ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc ./ET/era5_deterministic_recent.e.2.5deg.JJA.1950-2021.nc")
system("cdo -seasmean -select,season=JJA ./RH/orig/era5_deterministic_recent.rh.2.5deg.1m.1950-2021.nc ./RH/era5_deterministic_recent.rh.2.5deg.JJA.1950-2021.nc")


#seas anom
system("cdo -ymonsub ./T/era5_deterministic_recent.t2m.2.5deg.JJA.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./T/era5_deterministic_recent.t2m.2.5deg.JJA.1950-2021.nc ./T/era5_deterministic_recent.t2m.2.5deg.JJA.1950-2021_anom.nc")
system("cdo -ymonsub ./Q/era5_deterministic_recent.q.2.5deg.JJA.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./Q/era5_deterministic_recent.q.2.5deg.JJA.1950-2021.nc ./Q/era5_deterministic_recent.q.2.5deg.JJA.1950-2021_anom.nc")
system("cdo -ymonsub ./SLP/era5_deterministic_recent.msl.2.5deg.JJA.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./SLP/era5_deterministic_recent.msl.2.5deg.JJA.1950-2021.nc ./SLP/era5_deterministic_recent.msl.2.5deg.JJA.1950-2021_anom.nc")
system("cdo -ymonsub ./ET/era5_deterministic_recent.e.2.5deg.JJA.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./ET/era5_deterministic_recent.e.2.5deg.JJA.1950-2021.nc ./ET/era5_deterministic_recent.e.2.5deg.JJA.1950-2021_anom.nc")
system("cdo -ymonsub ./R/era5_deterministic_recent.rnet.2.5deg.JJA.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./R/era5_deterministic_recent.rnet.2.5deg.JJA.1950-2021.nc ./R/era5_deterministic_recent.rnet.2.5deg.JJA.1950-2021_anom.nc")
system("cdo -ymonsub ./P/era5_deterministic_recent.PmE.2.5deg.JJA.1950-2021.nc -ymonmean -seldate,1950-01-01,1999-12-31 ./P/era5_deterministic_recent.PmE.2.5deg.JJA.1950-2021.nc ./P/era5_deterministic_recent.PmE.2.5deg.JJA.1950-2021_anom.nc")
system("cdo -seasmean -select,season=JJA ./RH/era5_deterministic_recent.rh.2.5deg.1m.1950-2021_anom.nc ./RH/era5_deterministic_recent.rh.2.5deg.JJA.1950-2021_anom.nc")


get_obs_XAX <- function(file){
  file_nc <- nc_open(file)
  obs_data  <- ncvar_get(nc = file_nc)
  obs_time <- data.table(stringr::str_split_fixed(ncdf4.helpers::nc.get.time.series(file_nc), "-", n= 3))
  nc_close(file_nc)
  
  XAX <- list()
  XAX$X <- t(apply(obs_data, FUN = c, MARGIN = 3))
  XAX$M$year <- as.integer(obs_time$V1)
  XAX$M$mon <- as.integer(obs_time$V2)
  return(XAX)
}

#XAX JJA
ERA5_T_JJA_2d50_XAX <- get_obs_XAX("./T/era5_deterministic_recent.t2m.2.5deg.JJA.1950-2021.nc")
ERA5_q_JJA_2d50_XAX <- get_obs_XAX("./Q/era5_deterministic_recent.q.2.5deg.JJA.1950-2021.nc")
ERA5_q_JJA_2d50_XAX$X <- ERA5_q_JJA_2d50_XAX$X *100
ERA5_SLP_JJA_2d50_XAX <- get_obs_XAX("./SLP/era5_deterministic_recent.msl.2.5deg.JJA.1950-2021.nc")
ERA5_Rnet_JJA_2d50_XAX <- get_obs_XAX("./R/era5_deterministic_recent.rnet.2.5deg.JJA.1950-2021.nc")
ERA5_ET_JJA_2d50_XAX <- get_obs_XAX("./ET/era5_deterministic_recent.e.2.5deg.JJA.1950-2021.nc")
ERA5_ET_JJA_2d50_XAX$X <- ERA5_ET_JJA_2d50_XAX$X*-1000



save(list = "ERA5_T_JJA_2d50_XAX", file = "./T/ERA5_T_JJA_2d50_XAX.RData")
save(list = "ERA5_q_JJA_2d50_XAX", file = "./Q/ERA5_q_JJA_2d50_XAX.RData")
save(list = "ERA5_SLP_JJA_2d50_XAX", file = "./SLP/ERA5_SLP_JJA_2d50_XAX.RData")
save(list = "ERA5_Rnet_JJA_2d50_XAX", file = "./R/ERA5_Rnet_JJA_2d50_XAX.RData")
save(list = "ERA5_ET_JJA_2d50_XAX", file = "./ET/ERA5_ET_JJA_2d50_XAX.RData")


#XAX JJA anom
ERA5_T_JJA_2d50_anom_XAX <- get_obs_XAX("./T/era5_deterministic_recent.t2m.2.5deg.JJA.1950-2021_anom.nc")
ERA5_q_JJA_2d50_anom_XAX <- get_obs_XAX("./Q/era5_deterministic_recent.q.2.5deg.JJA.1950-2021_anom.nc")
ERA5_q_JJA_2d50_anom_XAX$X <- ERA5_q_JJA_2d50_anom_XAX$X *100
ERA5_SLP_JJA_2d50_anom_XAX <- get_obs_XAX("./SLP/era5_deterministic_recent.msl.2.5deg.JJA.1950-2021_anom.nc")
ERA5_Rnet_JJA_2d50_anom_XAX <- get_obs_XAX("./R/era5_deterministic_recent.rnet.2.5deg.JJA.1950-2021_anom.nc")
ERA5_ET_JJA_2d50_anom_XAX <- get_obs_XAX("./ET/era5_deterministic_recent.e.2.5deg.JJA.1950-2021_anom.nc")
ERA5_ET_JJA_2d50_anom_XAX$X <- ERA5_ET_JJA_2d50_anom_XAX$X*-1000
ERA5_RH_JJA_2d50_anom_XAX <- get_obs_XAX("./RH/era5_deterministic_recent.rh.2.5deg.JJA.1950-2021_anom.nc")


save(list = "ERA5_T_JJA_2d50_anom_XAX", file = "./T/ERA5_T_JJA_2d50_anom_XAX.RData")
save(list = "ERA5_q_JJA_2d50_anom_XAX", file = "./Q/ERA5_q_JJA_2d50_anom_XAX.RData")
save(list = "ERA5_SLP_JJA_2d50_anom_XAX", file = "./SLP/ERA5_SLP_JJA_2d50_anom_XAX.RData")
save(list = "ERA5_Rnet_JJA_2d50_anom_XAX", file = "./R/ERA5_Rnet_JJA_2d50_anom_XAX.RData")
save(list = "ERA5_ET_JJA_2d50_anom_XAX", file = "./ET/ERA5_ET_JJA_2d50_anom_XAX.RData")
save(list = "ERA5_RH_JJA_2d50_anom_XAX", file = "./RH/ERA5_RH_JJA_2d50_anom_XAX.RData")



#XAX T
ERA5_T_mon_2d50_XAX <- get_obs_XAX("./T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021.nc")
save(list = "ERA5_T_mon_2d50_XAX", file = "./T/ERA5_T_mon_2d50_XAX.RData")

#XAX T anom
ERA5_T_mon_2d50_anom_XAX <- get_obs_XAX("./T/era5_deterministic_recent.t2m.2.5deg.1m.1950-2021_anom.nc")

#XAX q 
ERA5_q_mon_2d50_XAX <- get_obs_XAX("./Q/era5_deterministic_recent.q.2.5deg.1m.1950-2021.nc")
ERA5_q_mon_2d50_XAX$X <- ERA5_q_mon_2d50_XAX$X *100

save(list = "ERA5_q_mon_2d50_XAX", file = "./Q/ERA5_q_mon_2d50_XAX.RData")

#XAX q anom
ERA5_q_mon_2d50_anom_XAX <- get_obs_XAX("./Q/era5_deterministic_recent.q.2.5deg.1m.1950-2021_anom.nc")
ERA5_q_mon_2d50_anom_XAX$X <- ERA5_q_mon_2d50_anom_XAX*100

#XAX SLP
ERA5_SLP_mon_2d50_XAX <- get_obs_XAX("./SLP/era5_deterministic_recent.msl.2.5deg.1m.1950-2021.nc")
save(list = "ERA5_SLP_mon_2d50_XAX", file = "./SLP/ERA5_SLP_mon_2d50_XAX.RData")


#XAX SLP anom
ERA5_SLP_mon_2d50_anom_XAX <- get_obs_XAX("./SLP/era5_deterministic_recent.msl.2.5deg.1m.1950-2021_anom.nc")


save(list = "ERA5_SLP_mon_2d50_anom_XAX", file = "./SLP/ERA5_SLP_mon_2d50_anom_XAX.RData")
save(list = "ERA5_q_mon_2d50_anom_XAX", file = "./Q/ERA5_q_mon_2d50_anom_XAX.RData")
save(list = "ERA5_T_mon_2d50_anom_XAX", file = "./T/ERA5_T_mon_2d50_anom_XAX.RData")


#XAX ET
ERA5_ET_mon_2d50_XAX <- get_obs_XAX("./ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021.nc")
ERA5_ET_mon_2d50_XAX$X <- ERA5_ET_mon_2d50_XAX$X*-1000
save(list = "ERA5_ET_mon_2d50_XAX", file = "./ET/ERA5_ET_mon_2d50_XAX.RData")

#XAX ET anom
ERA5_ET_mon_2d50_anom_XAX <- get_obs_XAX("./ET/era5_deterministic_recent.e.2.5deg.1m.1950-2021_anom.nc")
ERA5_ET_mon_2d50_anom_XAX$X <- ERA5_ET_mon_2d50_anom_XAX*-1000
save(list = "ERA5_ET_mon_2d50_anom_XAX", file = "./ET/ERA5_ET_mon_2d50_anom_XAX.RData")

#XAX P
ERA5_P_mon_2d50_XAX <- get_obs_XAX("./P/era5_deterministic_recent.tp.2.5deg.1m.1950-2021.nc")
ERA5_P_mon_2d50_XAX$X <- ERA5_P_mon_2d50_XAX$X *1000
save(list = "ERA5_P_mon_2d50_XAX", file = "./P/ERA5_P_mon_2d50_XAX.RData")

#XAX Rnet
ERA5_Rnet_mon_2d50_XAX <- get_obs_XAX("./R/era5_deterministic_recent.rnet.2.5deg.1m.1950-2021.nc")
save(list = "ERA5_Rnet_mon_2d50_XAX", file = "./R/ERA5_Rnet_mon_2d50_XAX.RData")

#XAX Rnet anom
ERA5_Rnet_mon_2d50_anom_XAX <- get_obs_XAX("./R/era5_deterministic_recent.rnet.2.5deg.1m.1950-2021_anom.nc")
save(list = "ERA5_Rnet_mon_2d50_anom_XAX", file = "./R/ERA5_Rnet_mon_2d50_anom_XAX.RData")

#XAX PmE
ERA5_PmE_mon_2d50_anom_XAX <- get_obs_XAX("./P/era5_deterministic_recent.PmE_mm.2.5deg.1m.1950-2021_anom.nc")
save(list = "ERA5_PmE_mon_2d50_anom_XAX", file = "./P/ERA5_PmE_mon_2d50_anom_XAX.RData")






