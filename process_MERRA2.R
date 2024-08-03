library(ncdf4)
library(data.table)


setwd("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/orig/")
#t2m
system(paste0("/usr/local/cdo-2.1.1/bin/cdo -O mergetime -apply,-selvar,T2M [ ./bin_t2m/MERRA2_*.tavgM_2d_slv_Nx.*.nc4 ] MERRA2_t2m.tavgM_1980-2023.nc"))

remap_file <- "/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc"
  
system(paste0("cdo remapcon2,/net/xenon/climphys_backedup/maegli/CMIP6/data/NetCDFs/hfls/hist585/hfls_mon_ACCESS-CM2_hist585_r1i1p1f1_2d50.nc ", 
              "./MERRA2_t2m.tavgM_1980-2023.nc ",
              "./MERRA2_t2m.tavgM_1980-2023_2.5.nc"))

system("cdo -ymonsub MERRA2_t2m.tavgM_1980-2023_2.5.nc -ymonmean -seldate,1980-01-01,1999-12-31 MERRA2_t2m.tavgM_1980-2023_2.5.nc MERRA2_t2m.tavgM_1980-2023_2.5_anom.nc")
system("cdo -seasmean -select,season=JJA MERRA2_t2m.tavgM_1980-2023_2.5_anom.nc MERRA2_t2m.tavg_JJA_1980-2023_2.5_anom.nc")

#RH
system(paste0("/usr/local/cdo-2.1.1/bin/cdo -O mergetime -apply,'-remapcon,", remap_file,
              " -selname,T2M,T2MDEW' [ ./bin_t2m/MERRA2_*.tavgM_2d_slv_Nx.*.nc4 ] MERRA2_allvar.tavgM_1980-2023_2.5.nc"))

system("cdo exprf,expr_d2mtoRH MERRA2_allvar.tavgM_1980-2023_2.5.nc MERRA2_RH.tavgM_1980-2023_2.5.nc")
system("cdo -ymonsub MERRA2_RH.tavgM_1980-2023_2.5.nc -ymonmean -seldate,1980-01-01,1999-12-31 MERRA2_RH.tavgM_1980-2023_2.5.nc MERRA2_RH.tavgM_1980-2023_2.5_anom.nc")
system("cdo -seasmean -select,season=JJA MERRA2_RH.tavgM_1980-2023_2.5_anom.nc MERRA2_RH.tavg_JJA_1980-2023_2.5_anom.nc")

#SLP
system(paste0("/usr/local/cdo-2.1.1/bin/cdo -O mergetime -apply,'-remapbil,", remap_file,
              " -selname,SLP' [ ./bin_t2m/MERRA2_*.tavgM_2d_slv_Nx.*.nc4 ] MERRA2_SLP.tavgM_1980-2023_2.5.nc"))
system("cdo -ymonsub MERRA2_SLP.tavgM_1980-2023_2.5.nc -ymonmean -seldate,1980-01-01,1999-12-31 MERRA2_SLP.tavgM_1980-2023_2.5.nc MERRA2_SLP.tavgM_1980-2023_2.5_anom.nc")
system("cdo -seasmean -select,season=JJA MERRA2_SLP.tavgM_1980-2023_2.5_anom.nc MERRA2_SLP.tavg_JAA_1980-2023_2.5_anom.nc")

#Rnet
system(paste0("/usr/local/cdo-2.1.1/bin/cdo -O mergetime -apply,'-remapbil,", remap_file,
              " -selname,LWGNT,SWGNT' [ ./bin_rad/MERRA2_*.tavgM_2d_rad_Nx.*.nc4 ] MERRA2_rad.tavgM_1980-2023_2.5.nc"))

system("cdo expr,netR=LWGNT+SWGNT MERRA2_rad.tavgM_1980-2023_2.5.nc MERRA2_netrad.tavgM_1980-2023_2.5.nc")
system("cdo -ymonsub MERRA2_netrad.tavgM_1980-2023_2.5.nc -ymonmean -seldate,1980-01-01,1999-12-31 MERRA2_netrad.tavgM_1980-2023_2.5.nc MERRA2_netrad.tavgM_1980-2023_2.5_anom.nc")
system("cdo -seasmean -select,season=JJA MERRA2_netrad.tavgM_1980-2023_2.5_anom.nc MERRA2_netrad.tavg_JJA_1980-2023_2.5_anom.nc")


#ET
system(paste0("/usr/local/cdo-2.1.1/bin/cdo -O mergetime -apply,'-remapbil,", remap_file,
              " -selname,EFLUX' [ ./bin_flx/MERRA2_*.tavgM_2d_flx_Nx.*.nc4 ] MERRA2_ET.tavgM_1980-2023_2.5.nc"))
system("cdo -ymonsub MERRA2_ET.tavgM_1980-2023_2.5.nc -ymonmean -seldate,1980-01-01,1999-12-31 MERRA2_ET.tavgM_1980-2023_2.5.nc MERRA2_ET.tavgM_1980-2023_2.5_anom.nc")
system("cdo -seasmean -select,season=JJA MERRA2_ET.tavgM_1980-2023_2.5_anom.nc MERRA2_ET.tavg_JJA_1980-2023_2.5_anom.nc")





MERRA_T_JJA_2d50_anom_XAX <- get_obs_XAX("MERRA2_t2m.tavg_JJA_1980-2023_2.5_anom.nc")
MERRA_RH_JJA_2d50_anom_XAX <- get_obs_XAX("MERRA2_RH.tavg_JJA_1980-2023_2.5_anom.nc")
MERRA_SLP_JJA_2d50_anom_XAX <-  get_obs_XAX("MERRA2_SLP.tavg_JAA_1980-2023_2.5_anom.nc")
MERRA_Rnet_JJA_2d50_anom_XAX <-  get_obs_XAX("MERRA2_netrad.tavg_JJA_1980-2023_2.5_anom.nc")
MERRA_ET_JJA_2d50_anom_XAX <-  get_obs_XAX("MERRA2_ET.tavg_JJA_1980-2023_2.5_anom.nc")
MERRA_ET_JJA_2d50_anom_XAX$X <- MERRA_ET_JJA_2d50_anom_XAX$X * 3600 * 24 /(2.5*10^6)
MERRA_ET_JJA_2d50_anom_XAX$M$obs <- "MERRA2"

save(list = "MERRA_T_JJA_2d50_anom_XAX", file = "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_T_JJA_2d50_anom_XAX.RData")
save(list = "MERRA_RH_JJA_2d50_anom_XAX", file = "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_RH_JJA_2d50_anom_XAX.RData")
save(list = "MERRA_SLP_JJA_2d50_anom_XAX", file = "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_SLP_JJA_2d50_anom_XAX.RData")
save(list = "MERRA_Rnet_JJA_2d50_anom_XAX", file = "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_Rnet_JJA_2d50_anom_XAX.RData")
save(list = "MERRA_ET_JJA_2d50_anom_XAX", file = "/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_ET_JJA_2d50_anom_XAX.RData")

