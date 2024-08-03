library(data.table)
library(sf)
library(ggplot2)
library(ggpubr)

load("/net/xenon/climphys_backedup/maegli/CMIP6/data/dist_mat.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_tas_JJA_hist585_anom_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hurs_JJA_hist585_anom_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hfls_JJA_hist585_anom_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_rnet_JJA_hist585_anom_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/land_mask.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_pr_JJA_hist585_anom_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_psl_JJA_hist585_anom_2d50_XAX.RData")

load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_MCET_JJA_hist585_anom_2d50_XAX.RData")


load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netlw_JJA_hist585_anom_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netsw_JJA_hist585_anom_2d50_XAX.RData")


source("/net/xenon/climphys_backedup/maegli/CMIP6/Pr_vs_ET_functions.R")
source("/net/xenon/climphys_backedup/maegli/ET_Adj/Scripts/00_dynadj_functions.R")
renice_session()

latlonDT <- get_lonlatDT(SREX = T)
latlonDT[, ix := .I]
CMIP6_models <- unique(CMIP6_MCET_JJA_hist585_anom_2d50_XAX$M$mod)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
SREX_sel <- c("WCE", "MED", "EAS", "SAS", "WNA", "NSA")

FR_mat_tas <- get_FR_mat(CMIP6_tas_JJA_hist585_anom_2d50_XAX)
FR_mat_hurs <- get_FR_mat(CMIP6_hurs_JJA_hist585_anom_2d50_XAX)
FR_mat_rnet <- get_FR_mat(CMIP6_rnet_JJA_hist585_anom_2d50_XAX)
FR_mat_hfls <- get_FR_mat(CMIP6_hfls_JJA_hist585_anom_2d50_XAX)
FR_mat_netlw <- get_FR_mat(CMIP6_netlw_JJA_hist585_anom_2d50_XAX)
FR_mat_netsw <- get_FR_mat(CMIP6_netsw_JJA_hist585_anom_2d50_XAX)
FR_mat_pr <- get_FR_mat(CMIP6_pr_JJA_hist585_anom_2d50_XAX)
FR_mat_mccoll <- get_FR_mat(CMIP6_MCET_JJA_hist585_anom_2d50_XAX)


#===============================
#==============OLS==============
#===============================

#==========LMO CV================
train_start <- 1950
train_end <- 2100

test_start <- 1950
test_end <- 2100


lm_list <- list()
beta_list <- list()
lmo_cv_list <- list()
lmo_cv <- list()

CMIP6_models_all <- CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$mod


cv_folds_GCMs <- data.table(mod = unique(CMIP6_models_all),
                            #fold = c("ACCESS", "ACCESS", "CanESM", "CESM", "CNRM", "CNRM", "EC-Earth", "EC-Earth", "GISS", "GISS", "IPSL", "MIROC", "MIROC", "MPI", "MRI", "UKESM"),
                            fold = c("UKESM", "UKESM", "CanESM", "CESM", "CNRM", "CNRM", "EC-Earth", "EC-Earth", "notuse", "GISS", "GISS","notuse", "CNRM", "MIROC", "MIROC", "MPI", "MRI","notuse","UKESM"),
                            weights = as.numeric(19/table(CMIP6_models_all)) /sum(19/table(CMIP6_models_all)))




M_full <- CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M[cv_folds_GCMs, on = "mod"]


folds <- unique(cv_folds_GCMs[fold != "notuse", fold])

for(k in seq_along(folds)){
  test_mods <-  cv_folds_GCMs[mod %in% CMIP6_models][fold == folds[k]]$mod
  train_mods <- cv_folds_GCMs[mod %in% CMIP6_models][fold != folds[k]]$mod
  
  train_ix <- which(CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$year %in% train_start:train_end & CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$mod %in% train_mods)
  train_mccoll <- which(CMIP6_MCET_JJA_hist585_anom_2d50_XAX$M$year %in% train_start:train_end & CMIP6_MCET_JJA_hist585_anom_2d50_XAX$M$mod %in% train_mods)
  
  
  test_ix <- which(CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$year %in% test_start:test_end & CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$mod %in% test_mods)
  test_mccoll <- which(CMIP6_MCET_JJA_hist585_anom_2d50_XAX$M$year %in% test_start:test_end & CMIP6_MCET_JJA_hist585_anom_2d50_XAX$M$mod %in% test_mods)
  
  
  for(i in seq_along(land_mask)){
  
    x_data_train <- data.table(ET = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         Rnet = CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         RH = CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         SAT = CMIP6_tas_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         Pr = CMIP6_pr_JJA_hist585_anom_2d50_XAX$X[train_mccoll, land_mask[i]],
                         SLP = CMIP6_psl_JJA_hist585_anom_2d50_XAX$X[train_mccoll, land_mask[i]],
                         weights  = M_full$weights[train_ix])
    
    lm_base <- lm(ET ~ Rnet + RH +SAT, data = x_data_train, weights = weights)
    #lm_pr <- lm(ET ~ Rnet + RH +SAT + Pr, data = x_data_train, weights = weights)
    lm_noSAT <- lm(ET ~ Rnet + RH, data = x_data_train, weights = weights)
    lm_noRnet<- lm(ET ~ SAT + RH, data = x_data_train, weights = weights)
    lm_noRH <- lm(ET ~ Rnet + SAT, data = x_data_train, weights = weights)
    #lm_SLP <- lm(ET ~ Rnet + RH +SAT + SLP, data = x_data_train, weights = weights)
    
    x_test <- data.table(Rnet = CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X[test_ix, land_mask[i]],
                             RH = CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X[test_ix, land_mask[i]],
                             SAT = CMIP6_tas_JJA_hist585_anom_2d50_XAX$X[test_ix, land_mask[i]],
                             Pr = CMIP6_pr_JJA_hist585_anom_2d50_XAX$X[test_mccoll, land_mask[i]],
                         SLP = CMIP6_psl_JJA_hist585_anom_2d50_XAX$X[test_mccoll, land_mask[i]]
                              )
    
    
    pred_base <- predict(lm_base, x_test)
    #pred_pr <- predict(lm_pr, x_test)
    #pred_SLP <-  predict(lm_SLP, x_test)
    noSAT <- predict(lm_noSAT, x_test)
    noRnet <- predict(lm_noRnet, x_test)
    noRH <- predict(lm_noRH, x_test)
    
    lmo_cv_list[[i]] <- data.table(base = pred_base,
                                   #pr = pred_pr,
                                   #SLP = pred_SLP,
                                   noSAT = noSAT,
                                   noRnet = noRnet,
                                   noRH = noRH,
                                   McColl = CMIP6_MCET_JJA_hist585_anom_2d50_XAX$X[test_mccoll, land_mask[i]],
               year = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$year[test_ix],
               mod = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$mod[test_ix],
               mem = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$ens.mem[test_ix],
               fold = M_full$fold[test_ix],
               true = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[test_ix, land_mask[i]],
               latlonDT[land_mask[i]])
    
  }
  lmo_cv[[k]] <- rbindlist(lmo_cv_list)
}

lmo_cv_all <- rbindlist(lmo_cv) %>% melt(id.vars = c("year", "fold", "mod", "true", "lat", "lon", "lonlat", "SREX", "ix", "mem"))
lmo_cv_ensmean <- rbindlist(lmo_cv) %>% melt(id.vars = c("year", "fold", "mod", "lat", "lon", "lonlat", "SREX", "ix", "mem")) %>% 
  .[, .(value = mean(value)), .(fold, SREX, year, variable)]


  
bias <- function(y, y_pred)  abs(mean(y) - mean(y_pred))

lmo_cv_mse <- lmo_cv_all[, .(rmse = sqrt(mse(value, true)), bias = bias(true, value)), .(lat, lon, lonlat, variable, fold)]

lmo_cv_mse[, .(rmse = mean(rmse)), .(lat, lon, variable)] %>% dcast(lat + lon ~ variable, value.var = "rmse") %>% 
  .[, rmse := noRH/base] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = rmse))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  #facet_wrap(~fold)+
  scale_fill_gradientn(name = "Ratio", rescaler = ~ scales::rescale_mid(.x, mid = 1),
                       colours=c("dodgerblue3", "white", "firebrick2"))+
  theme(panel.grid.major = element_blank()) +
  ggtitle(paste0("Train ", train_start, "-", train_end, " | Test ",test_start,"-",test_end," | RMSE(SAT + Rnet) / RMSE(SAT + RH + Rnet)"))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  ylim(-55,80) 

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision//",
       filename = paste0("lmo_noRH_",train_start, "-",train_end,"_",test_start,"-", test_end, ".png"),
       height = 7, width = 12,
       bg = "white")


lmo_cv_mse[, .(bias = mean(bias)), .(lat, lon, variable)] %>% dcast(lat + lon ~ variable, value.var = "bias") %>% 
  .[, bias := McColl/base] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = bias))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  #facet_wrap(~fold)+
  scale_fill_gradientn(name = "RMSE", rescaler = ~ scales::rescale_mid(.x, mid = 1),
                       colours=c("blue", "white", "red"))+
  theme(panel.grid.major = element_blank()) +
  ggtitle(paste0("Train ", train_start, "-", train_end, " / Test ",test_start,"-",test_end," | Difference RMSE(SAT + RH + Rnet) - RMSE(SAT + RH + Rnet + Pr)"))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  ylim(-55,80) 


lmo_cv_ensmean[SREX %in% SREX_sel][, .(value = mean(value), sd = sd(value)), .(SREX, year, variable)] %>% 
  ggplot()+
  geom_line(aes(year, value, color = variable))+
  geom_ribbon(aes(year, ymin = value -sd, ymax = value + sd, fill = variable), alpha = 0.3)+
  facet_wrap(~SREX)+
  ylab("ET anomaly [mm/day]") + xlab("Year") +
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("goldenrod", "mediumorchid4", "forestgreen", "grey20"), labels = c("SAT + RH + Rnet","SAT + RH + Rnet +Pr", "SFE", "CMIP6 ET"), name = "")+
  scale_fill_manual(values = c("goldenrod", "mediumorchid4", "forestgreen", "grey20"), labels = c("SAT + RH + Rnet","SAT + RH + Rnet +Pr", "SFE", "CMIP6 ET"), name = "")+
  ggtitle(paste0("Train ", train_start, "-", train_end, " / Test ",test_start,"-",test_end,""))
  
ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision//",
       filename = paste0("SREX_ts_lmo_",train_start, "-",train_end,"_",test_start,"-", test_end, ".png"),
       height = 7, width = 12,
       bg = "white")


ggplot(lmo_cv_mse[, .(V1 = mean(V1)), .(lat, lon, variable)])+
  geom_tile(aes(lon, lat, fill = V1))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~variable)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9") +
  theme(panel.grid.major = element_blank()) 



#==================single model====================
sds_tas <- matrixStats::colSds(CMIP6_tas_JJA_hist585_anom_2d50_XAX$X)
sds_hurs <- matrixStats::colSds(CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X)
sds_rnet <- matrixStats::colSds(CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X)


FRintoIAV_lm <- list()
beta_list_RHRnetSAT <- list()
lm_list_RHRnetSAT <- list()
lm_list_RnetSAT <- list()
lm_list_RHRnet <- list()
lm_list_RHSAT <- list()

train_ix_all <- which(CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$year %in% 1980:2100 & CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$mod %in% CMIP6_models)
test_ix_all <-  which(FR_mat_hfls$M$mod %in% CMIP6_models)


for(i in seq_along(land_mask)){
  
  x_data_RHRnetSAT <- data.table(ET = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[train_ix_all, land_mask[i]],
                              RH = CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X[train_ix_all, land_mask[i]],
                              Rnet = CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X[train_ix_all, land_mask[i]],
                              SAT = CMIP6_tas_JJA_hist585_anom_2d50_XAX$X[train_ix_all, land_mask[i]],
                              weights = M_full$weights[train_ix_all]) 
  
  
  lm_list_RHRnetSAT[[i]] <- lm(ET~ RH + Rnet + SAT, data = x_data_RHRnetSAT, weights = weights)
  lm_list_RnetSAT[[i]] <- lm(ET~  Rnet + SAT, data = x_data_RHRnetSAT, weights = weights)
  lm_list_RHRnet[[i]] <- lm(ET~ RH + Rnet , data = x_data_RHRnetSAT, weights = weights)
  lm_list_RHSAT[[i]] <- lm(ET~ RH  + SAT, data = x_data_RHRnetSAT, weights = weights)
  
  
  beta_list_RHRnetSAT[[i]] <- data.table(beta = lm_list_RHRnetSAT[[i]]$coefficients,
                                     beta_scaled = lm_list_RHRnetSAT[[i]]$coefficients * c(1, 
                                                                                        sds_hurs[land_mask[i]], 
                                                                                        sds_rnet[land_mask[i]], 
                                                                                        sds_tas[land_mask[i]]),
                               var = names(lm_list_RHRnetSAT[[i]]$coefficients),
                               latlonDT[land_mask[i]])
  

  
  x_new_RHRnet <- data.table(RH = FR_mat_hurs$X[test_ix_all, land_mask[i]],
                             Rnet = FR_mat_rnet$X[test_ix_all, land_mask[i]],
                             SAT = FR_mat_tas$X[test_ix_all, land_mask[i]])
  
 lm_frame <- data.table(pred_RHRnetSAT = predict(lm_list_RHRnetSAT[[i]], x_new_RHRnet),
                                  year = FR_mat_hfls$M$year[test_ix_all],
                                  mod = FR_mat_hfls$M$mod[test_ix_all],
                                  true = FR_mat_hfls$X[test_ix_all, land_mask[i]],
                                  latlonDT[land_mask[i]] 
                                  ) %>% melt(id.vars = c("year", "mod", "lat", "lon", "lonlat", "SREX", "ix"))
  
  mccoll_frame <- data.table(value =FR_mat_mccoll$X[, land_mask[i]], 
                             variable = "McColl",
                             year = FR_mat_mccoll$M$year,
                             mod = FR_mat_mccoll$M$mod,
                             latlonDT[land_mask[i]])
  FRintoIAV_lm[[i]] <- rbind(lm_frame, mccoll_frame)
}


FRintoIAV_lm_frame <- rbindlist(FRintoIAV_lm)
FRintoIAV_lm_frame[cv_folds_GCMs]

FRintoIAV_lm_frame %>% 
  ggplot()+
  geom_hex(aes(pred, true))+
  geom_abline(xintercept = 0 , slope = 1)+
  facet_wrap(~mod)+
  scale_fill_viridis_b(trans = "log")
  
ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231003_lm_on_obs/",
       filename = "betas_SAT_vs_SATonly.png",
       height = 4, width = 12)


FRintoIAV_lm_frame[SREX %in% SREX_sel] %>%
  .[, .(value = mean(value)), .(mod, year, SREX, variable)] %>% 
  .[, .(value = mean(value), sd = sd(value)), .(year, SREX, variable)] %>% 
  ggplot()+
  geom_hline(yintercept = 0, linewidth = 0.2)+
  annotate("rect", xmin = 1980, xmax = 2020, ymin = - Inf, ymax = Inf,
           alpha = .3,fill = "grey60")+
  geom_line(aes(year, value, color = variable))+
  geom_ribbon(aes(year, ymin = value -sd, ymax = value + sd, fill = variable), alpha = 0.3)+
  facet_wrap(~SREX)+
  scale_color_manual(values = c("mediumorchid4", "grey20", "goldenrod"), labels = c("LM Prediction","CMIP6 ET", "McColl"), name = "")+
  scale_fill_manual(values = c("mediumorchid4", "grey20", "goldenrod"), labels = c("LM Prediction","CMIP6 ET", "McColl"), name = "")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  xlab("Year") +ylab("ET anomaly [mm/day]")

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240502_draft_update2/",
       filename = "train1980-2014.png",
       height = 7, width = 12,
       bg = "white")


FRintoIAV_lm_frame[SREX == "WCE"] %>%
  melt(id.vars = c("year", "mod", "lat", "lon", "lonlat", "SREX", "ix")) %>% 
  .[variable == "true"] %>% 
  .[, .(value = mean(value)), .(mod, year, SREX, variable)] %>% 
  ggplot()+
  geom_line(aes(year, value))+
  facet_wrap(~mod)+
  #scale_color_manual(values = c("grey20", "goldenrod"), labels = c("CMIP6 ET", "LM Prediction"), name = "")+
  theme_minimal()+
  ggtitle("WCE") +
  xlab("Year") +ylab("ET anomaly [mm/day]")


ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231119_groupmeeting/",
       filename = "WCE_FR_hist585.png",
       height = 7, width = 12,
       bg = "white")


p1 <- rbindlist(beta_list)[var != "(Intercept)"] %>% 
ggplot()+
  geom_tile(aes(lon, lat, fill = beta))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  facet_wrap(~var)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9", name = "") +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  ylim(-55,80) +
  ggtitle("Coefficients on original scale")

p2 <- rbindlist(beta_list)[var != "(Intercept)"] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = beta_scaled))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  facet_wrap(~var)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9", name = "") +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  ylim(-55,80) +
  ggtitle("Coefficients scaled with variance of variable")

ggpubr::ggarrange(p1, p2, ncol = 1)


ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240227_draft_update//",
       filename = "beta_lwsw_scaled_plot.png",
       height = 8, width = 12,
       bg = "white")

p1 <- rbindlist(beta_list_RHRnetSAT)[var != "(Intercept)"] %>% 
ggplot()+
  geom_tile(aes(lon, lat, fill = beta))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  facet_wrap(~var)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9", name = "") +
  theme(panel.grid.major = element_blank())+
  ggtitle("Coefficients on original scale")

p2 <- rbindlist(beta_list_RHRnetSAT)[var != "(Intercept)"] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = beta_scaled))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  facet_wrap(~var)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9", name = "") +
  theme(panel.grid.major = element_blank())+
  ggtitle("Coefficients scaled with variance of variable")

ggpubr::ggarrange(p1, p2, ncol = 1)
ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240227_draft_update//",
       filename = "beta_rnet_scaled_plot.png",
       height = 5, width = 11,
       bg = "white")



#==========contribtuions RH + SAT + Rnet======
hist585_var_RHSATRnet_list <- list()

for(i in seq_along(land_mask)){
  hist585_newdata_rh <- data.table(RH = FR_mat_hurs$X[, land_mask[i]],
                                   SAT = 0,
                                   Rnet = 0)
  
  hist585_newdata_rnet <- data.table(RH = 0,
                                    Rnet = FR_mat_rnet$X[, land_mask[i]],
                                    SAT = 0)
  
  hist585_newdata_tas <- data.table(RH = 0,
                                     Rnet = 0,
                                     SAT = FR_mat_tas$X[, land_mask[i]])
  
  hist585_var_RHSATRnet_list[[i]] <- data.table(RH = predict(lm_list_RHRnetSAT[[i]], hist585_newdata_rh),
                                           Rnet = predict(lm_list_RHRnetSAT[[i]], hist585_newdata_rnet),
                                           SAT = predict(lm_list_RHRnetSAT[[i]], hist585_newdata_tas),
                                           ET = FR_mat_hfls$X[, land_mask[i]],
                                           year = FR_mat_tas$M$year,
                                           mod = FR_mat_tas$M$mod,
                                           latlonDT[land_mask[i]])
  
}

hist585_RHRnetSAT_sens <- rbindlist(hist585_var_RHSATRnet_list) %>% melt(id.vars= c("year", "mod", "lat", "lon", "lonlat", "SREX", "ix"))
hist585_RHRnetSAT_sens[, value_smooth := predict(loess(value ~ year, span = 0.3, degree = 2), year), .(mod, lat, lon, lonlat, SREX, variable)]



sens_baseline <- hist585_RHRnetSAT_sens[year %in% 1850:1900, .(baseline = mean(value)), .(mod, lat, lon, lonlat, variable)]

hist585_RHRnetSAT_sens <- hist585_RHRnetSAT_sens[sens_baseline,on =  c("mod", "lat", "lon", "lonlat", "variable")]

hist585_RHRnetSAT_sens[,value_smooth := value_smooth - baseline]

#var.labs <- c(bquote(ET[RH]), bquote(ET[R[net]]), bquote(ET[SAT]), "ET")
#names(var.labs) <- c("RH", "Rnet", "SAT", "ET")
#levels(hist585_RHRnetSAT_sens$variable) <- c(bquote(ET[RH]), bquote(ET[R[net]]), bquote(ET[SAT]), "ET")



hist585_RHRnetSAT_sens[SREX %in% SREX_sel
                       ][, .(value_smooth = mean(value_smooth)), .(year, variable, SREX, mod)] %>% 
  ggplot()+
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.5)+
  geom_line(aes(year, value_smooth, color = variable, group = mod))+ 
  facet_grid(SREX~variable)+
  scale_color_manual(values = c("steelblue", "goldenrod","forestgreen",  "grey30"), name  = "", labels = c(bquote(ET[RH]), bquote(ET[R[net]]), bquote(ET[SAT]), "ET")) +
  theme_classic()+
  theme(panel.grid = element_blank(), strip.background = element_blank())+ 
  ylab("ET anomaly [mm/day]") + xlab("Year")


ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240502_draft_update2/",
       filename = "hist585_sens.png",
       height = 10, width = 10, bg = "white")


#==============================



hist585_var_sens_smooth_FR[ix %in% latlonDT[SREX == "WCE", ix], .(sd = sd(value), mean = mean(value)), .(year, variable, mod)
                           ] %>% 
  ggplot()+ 
  geom_line(aes(year, mean, color = variable)) +
  facet_wrap(~mod)

hist585_var_sens_smooth_FR[ix %in% latlonDT[SREX == "WCE", ix], .(sd = sd(value), mean = mean(value)), .(year, variable, mod)
                           ][, .(ymin = min(mean), ymax = max(mean), mean = mean(mean)), .(year, variable)] %>% 
  ggplot()+
  geom_ribbon(aes(year, ymin = ymin, ymax = ymax, fill = variable), alpha = 0.3)+
  geom_line(aes(year, mean, color= variable))+
  scale_fill_manual(values = c("steelblue","goldenrod", "tomato2"), name  = "")+
  scale_color_manual(values = c("steelblue","goldenrod", "tomato2"), name  = "") +
  ylab("ET [mm/day]")

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231214_SWLW_and_constr/",
       filename = "FR_SWLW.png",
       height = 4, width = 6, bg = "white")

hist585_var_sens[variable != "ET"][SREX == "NEU", .(sd = sd(value_IAV), mean = mean(value_IAV)), .(year, variable, mod)] %>% 
  ggplot()+ 
  geom_line(aes(year, mean, color = variable)) +
  facet_wrap(~mod)+
  ylab("ET [mm/day]")+
  scale_color_manual(values = c("steelblue","goldenrod", "tomato2"), name  = "")
  

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231214_SWLW_and_constr/",
       filename = "IAV_SWLW.png",
       height = 8, width = 12, bg = "white")

#==============================
var_change_lonlat <- var_sens_summary %>% dcast(lonlat+lon+lat~eval, value.var = "variable") %>% 
  .[FR == IAV] %>% .[FR == "RH"] %>%  .$lonlat


var_change_ix <- latlonDT[lonlat %in% var_change_lonlat]$ix
hist585_var_sens[variable != "RHRnet"][ix %in%  var_change_ix
                   ][, .(value = mean(value)), .(year, lonlat, variable, mod)
                     ][, .(value = mean(value), sd = sd(value)), .(year, variable, mod)] %>% 
  ggplot()+
  geom_ribbon(aes(year, ymax = value+sd, ymin=value-sd, fill = variable), alpha = 0.3, color = NA)+
  geom_line(aes(year, value, color = variable))+
  scale_color_manual(values = c("steelblue", "goldenrod2","tomato2", "grey30"))+ 
  scale_fill_manual(values = c("steelblue", "goldenrod2","tomato2", "grey30"))+
  facet_wrap(~mod)

#MSE of trend
FR_trend_variable_hist585 <- hist585_var_sens[year %in% 2015:2100, .(value = mean(value), ET = mean(ET)), .(year, mod, variable, lonlat, lat, lon)
                                    ][, .(trend_pred = lm(value~year)$coefficients[2],
      trend_ET = lm(ET~year)$coefficients[2]), .(mod, variable, lonlat, lat, lon)
      ][, mse(trend_pred, trend_ET), .(variable, mod, lonlat, lat, lon)
        ][, .SD[which.min(V1)], .(mod, lonlat, lat, lon)]


#correlation of FR
FR_cor_sens <- hist585_var_sens[year %in% 2015:2100, .(value = mean(value), ET = mean(ET)), .(year, mod, variable, lonlat, lat, lon)
                 ][, .(value = loess(value ~ year)$fitted, ET = loess(ET ~ year)$fitted), .(year, mod, variable, lonlat, lon, lat)
                 ][,cor(value, ET), .(mod, variable, lonlat, lat, lon)
                   ][, .SD[which.max(V1)], .(mod, lonlat, lat, lon)]



IAV_corr_variable_hist585 <- hist585_var_sens[year %in% 2015:2100, cor(value_detrend, ET_detrend), .(mod, variable, lonlat, lat, lon)][, .SD[which.max(V1)], .(mod, lonlat, lat, lon)]

FR_trend_variable_hist585[, eval := "FR"]
IAV_corr_variable_hist585[, eval := "IAV"]
FR_cor_sens[, eval := "FR"]
FRvsIAV_trend_variable_hist585 <- rbind(FR_cor_sens,IAV_corr_variable_hist585)

FR_trend_variable_hist585 %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = v<ariable))+
  facet_wrap(~mod)+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  scale_color_manual(values = c("steelblue", "goldenrod2", "tomato2"), name  = "")+
  theme(panel.grid.major = element_blank())

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231024_mccoll_vs_lm/",
       filename = "IAV_sensitivity_hist585.png",
       height = 8, width = 12, bg = "white")



FRvsIAV_trend_variable_hist585[, .N, .(lonlat, lon, lat, variable, eval)][, .SD[which.max(N)], .(lonlat, lon, lat, eval)] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = variable, alpha = N))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  scale_fill_manual(values = c("tomato2", "steelblue", "goldenrod2"), name  = "")+
  theme(panel.grid.major = element_blank())+
  facet_wrap(~eval)


ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231024_mccoll_vs_lm/",
       filename = "IAVvsFR_sensitivity_hist585_summary.png",
       height = 8, width = 12, bg = "white")



SREX_avg_var_sens <- rbindlist(hist585_var_sens_list)[, .(ET = mean(ET), RH = mean(RH), Rnet = mean(Rnet)), .(year, mod, mem, lat, lon, lonlat, SREX)]
var_sens_IAV <- SREX_avg_var_sens[year %in% 1980:2020, .(slope_IAV = lm(Rnet~RH)$coefficients[2],
                                         intercept_IAV = lm(Rnet~RH)$coefficients[1]), .(mod, lat, lon, lonlat, SREX)]

var_sens_FR <- SREX_avg_var_sens[year %in% 1850:2100,.(ET = mean(ET), RH = mean(RH), Rnet = mean(Rnet)), .(year, mod, lat, lon, lonlat, SREX)
                                  ][, .(RH =  predict(loess(RH ~year), year),
                                        Rnet =  predict(loess(Rnet ~year), year),
                                        year), .(mod, lat, lon, lonlat, SREX)]


ggplot(var_sens_FR) +
  geom_path(aes(RH, Rnet, color = year, group = mod))+
  geom_abline(data = var_sens_FR[year %in% 2020:2100, .(trend = lm(Rnet ~ RH)$coefficient[2],
                                                        intercept = lm(Rnet ~RH)$coefficient[1]), .(mod)], 
                aes(slope =trend, intercept = intercept, group = mod))+
  scale_color_viridis_c()

SREX_sel <- "SES"
ggplot(SREX_avg_var_sens[SREX == SREX_sel, .(ET = mean(ET), RH = mean(RH), Rnet = mean(Rnet)), .(mod, year, mem)])+
  geom_point(aes(RH, Rnet, color = year), size = 0.2) +
  facet_wrap(~mod)+
  geom_abline(data = var_sens_IAV[SREX == SREX_sel, .(slope_IAV = mean(slope_IAV), intercept_IAV= mean(intercept_IAV)), .(mod)], aes(intercept = intercept_IAV, slope = slope_IAV, linetype = "IAV (1850-1950)"))+
  geom_path(data = var_sens_FR[SREX == SREX_sel, .(RH = mean(RH), Rnet = mean(Rnet)), .(mod, year)], aes(RH, Rnet, linetype = "FR (1850-2100)"))+ 
  scale_color_viridis_c()+
  ggtitle(SREX_sel) +
  coord_equal()+
  theme_minimal()

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231012_decomp/",
       filename = paste0("RHvsRnet_contribution_scatter_", SREX_sel , ".png"),
       height = 8, width = 12, bg = "white")

var_sens_FR_trend <- var_sens_FR[year %in% 1980:2020, .(trend = lm(Rnet ~ RH)$coefficient[2]), .(mod, lat, lon, lonlat, SREX)][abs(trend)>1, var := "Rnet"][abs(trend)<1,var := "RH"]

var_sens_IAV[abs(slope_IAV)>1, var := "Rnet"][abs(slope_IAV)<1,var := "RH"]
var_sens_IAV_summary <- var_sens_IAV[, .N, .(var,lon, lat, lonlat)][, .SD[which.max(N)], .(lat, lon, lonlat)]
var_sens_FR_trend_summary <- var_sens_FR_trend[, .N, .(var,lon, lat, lonlat)][, .SD[which.max(N)], .(lat, lon, lonlat)]

var_sens_IAV_summary[,eval := "IAV"]
var_sens_FR_trend_summary[, eval := "FR"]

rbind(var_sens_IAV_summary, var_sens_FR_trend_summary) %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = var, alpha = N))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_grid(~eval)+
  theme(panel.grid.major = element_blank())

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20231012_decomp/",
       filename = "RHvsRnet_contribution_FR_vs_IAV.png",
       height = 6, width = 12, bg = "white")


#================Reanalysis==================
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/ET/ERA5_ET_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/RH/ERA5_RH_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/T/ERA5_T_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5/R/ERA5_Rnet_JJA_2d50_anom_XAX.RData")



load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5-Land/ERA5Land_tas_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5-Land/ERA5Land_RH_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5-Land/ERA5Land_rnet_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5-Land/ERA5Land_ET_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_Rnet_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_RH_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_T_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/MERRA2/MERRA_ET_JJA_2d50_anom_XAX.RData")

obs_eval_lm_list <- list()

for(i in seq_along(land_mask)){
  era5_newdata <- data.table(SAT = ERA5_T_JJA_2d50_anom_XAX$X[, land_mask[i]],
                             RH = ERA5_RH_JJA_2d50_anom_XAX$X[, land_mask[i]],
                             Rnet =ERA5_Rnet_JJA_2d50_anom_XAX$X[, land_mask[i]])

era5land_newdata <- data.table(SAT = ERA5Land_tas_JJA_2d50_anom_XAX$X[, land_mask[i]],
                           RH = ERA5Land_RH_JJA_2d50_anom_XAX$X[, land_mask[i]],
                           Rnet = ERA5Land_rnet_JJA_2d50_anom_XAX$X[, land_mask[i]])

merra2_newdata <- data.table(SAT = MERRA_T_JJA_2d50_anom_XAX$X[, land_mask[i]],
                           RH = MERRA_RH_JJA_2d50_anom_XAX$X[, land_mask[i]],
                           Rnet = MERRA_Rnet_JJA_2d50_anom_XAX$X[, land_mask[i]])


obs_eval_lm_list[[i]] <- rbind(data.table("SAT + RH + Rnet" = predict(lm_list_RHRnetSAT[[i]], merra2_newdata),
                                          "RH + SAT" = predict(lm_list_RHSAT[[i]], merra2_newdata),
                                          "Rnet + SAT" = predict(lm_list_RnetSAT[[i]], merra2_newdata),
                                          "RH + Rnet" = predict(lm_list_RHRnet[[i]], merra2_newdata),
                                          year = MERRA_T_JJA_2d50_anom_XAX$M$year,
                                          ET = MERRA_ET_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                          obs = "MERRA2",
                                          latlonDT[land_mask[i]]),
                               
data.table("SAT + RH + Rnet" = predict(lm_list_RHRnetSAT[[i]], era5land_newdata),
           "RH + SAT" = predict(lm_list_RHSAT[[i]], era5land_newdata),
           "Rnet + SAT" = predict(lm_list_RnetSAT[[i]], era5land_newdata),
           "RH + Rnet" = predict(lm_list_RHRnet[[i]], era5land_newdata),
           year = ERA5Land_ET_JJA_2d50_anom_XAX$M$year,
           ET = ERA5Land_ET_JJA_2d50_anom_XAX$X[, land_mask[i]],
           obs = "ERA5 Land",
           latlonDT[land_mask[i]]),

data.table("SAT + RH + Rnet" = predict(lm_list_RHRnetSAT[[i]], era5_newdata),
           "RH + SAT" = predict(lm_list_RHSAT[[i]], era5_newdata),
           "Rnet + SAT" = predict(lm_list_RnetSAT[[i]], era5_newdata),
           "RH + Rnet" = predict(lm_list_RHRnet[[i]], era5_newdata),
           year = ERA5Land_ET_JJA_2d50_anom_XAX$M$year,
           ET = ERA5Land_ET_JJA_2d50_anom_XAX$X[, land_mask[i]],
           obs = "ERA5",
           latlonDT[land_mask[i]]))

           
}

obs_eval_lm <- rbindlist(obs_eval_lm_list) %>% melt(id.vars = c("year","ET","obs","lat","lon","lonlat","SREX","ix"))
obs_eval_lm_cor <- obs_eval_lm[, .(cor = cor(value, ET)), .(obs, lat, lon, variable)]
  

obs_eval_lm_cor %>% 
  ggplot()+
  geom_density(aes(cor, color = variable), linewidth = 0.6) +
  scale_color_brewer(palette = "Set2", name ="Predictors")+
  scale_fill_brewer(palette = "Set2", name ="Predictors")+
  xlab("Correlation distribution")+ ylab("Density")+
  theme_minimal()+
  facet_wrap(~obs)#+
  theme(panel.grid.minor = element_blank()) 

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision/",
       filename = "cor_distribution_era5vsera5land.png",
       height = 6, width = 14, bg = "white")


pal <-  RColorBrewer::brewer.pal(6, "Spectral") %>% rev()

scale_fill_fermenter_custom <- function(pal, na.value = "grey50", guide = "coloursteps", aesthetics = "fill", ...) {
  binned_scale("fill", "fermenter", ggplot2:::binned_pal(scales::manual_pal(unname(pal))), na.value = na.value, guide = guide, ...)  
}


cormap_obs <- obs_eval_lm_cor[cor < 0, cor := -0.1] %>% 
ggplot()+
  geom_tile(aes(lon, lat, fill = cor))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~obs)+
  scale_fill_fermenter_custom(pal, n.break = 7,  name = "Correlation")+
  theme_minimal()+
  ylim(-55,80)+
  theme(panel.grid.major = element_blank(), axis.text = element_blank()) 





obs_eval_lm_cor %>% dcast(value.var = "cor", variable+lat+lon~obs) %>% setnames("ERA5 Land", "ERA5Land") %>% .[, diff := ERA5Land - ERA5] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = diff))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  scale_fill_fermenter_custom(pal, n.break = 7, name = "Difference")+
  theme_minimal()+
  ggtitle("Correlation difference (ERA5 Land - ERA5)")+
  ylim(-55,80)+
  theme(panel.grid.major = element_blank(), axis.text = element_blank()) 

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision/",
       filename = "cor_diff_era5vsera5land.png",
       height = 5, width = 8, bg = "white")

cormap_gcms <- lmo_cv_mse[variable != "pr", .(V1 = mean(V1)), .(lat, lon, fold, variable)
           ][V1 < 0, V1:= -0.1
             ][, variable := "Linear Model"] %>% 
ggplot()+
  geom_tile(aes(lon, lat, fill = V1))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~variable)+
  scale_fill_fermenter_custom(pal, n.break = 7, name = "Correlation")+
  theme_minimal()+
  ylim(-55,80)+
  theme(panel.grid.major = element_blank(), axis.text = element_blank()) 

ggpubr::ggarrange(cormap_gcms, cormap_obs, nrow = 2, common.legend = T, legend = "right")

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240502_draft_update2//",
       filename = "cor_lm_obs.png",
       height = 6, width = 12, bg = "white")

obs_eval_lm[, .(cor = cor(value, ET)), .(obs, variable)] 

obs_eval_lm[variable == "pred_lm"][SREX == "NEU"][obs == "ERA5"] %>% 
  ggplot()+
  geom_line(aes(year, value, color = obs))+
  geom_line(aes(year, ET))+
  facet_wrap(~lonlat)

smooth_obs <- function(XAX){
  obs_smooth <- data.table(XAX$X, 
             year = XAX$M$year)[, lapply(.SD, function(x) loess(x~year)$fitted)]
  
  obs_smooth[, year := NULL]
  obs_smooth <- as.matrix(obs_smooth)
  return(obs_smooth)
}

obs_smooth_lm_list <- list()
ERA5_T_smooth <- smooth_obs(ERA5_T_JJA_2d50_anom_XAX)
ERA5_RH_smooth <- smooth_obs(ERA5_RH_JJA_2d50_anom_XAX)
ERA5_Rnet_smooth <- smooth_obs(ERA5_Rnet_JJA_2d50_anom_XAX)
ERA5_ET_smooth <- smooth_obs(ERA5_ET_JJA_2d50_anom_XAX)

MERRA2_T_smooth <- smooth_obs(MERRA_T_JJA_2d50_anom_XAX)
MERRA2_RH_smooth <- smooth_obs(MERRA_RH_JJA_2d50_anom_XAX)
MERRA2_Rnet_smooth <- smooth_obs(MERRA_Rnet_JJA_2d50_anom_XAX)
MERRA2_ET_smooth <- smooth_obs(MERRA_ET_JJA_2d50_anom_XAX)

obs_smooth_decomp_list <- list()
for(i in seq_along(land_mask)){

era5_smooth_newdata_T <- data.table(SAT = ERA5_T_JJA_2d50_anom_XAX$X[, land_mask[i]],
                           RH = 0,
                           Rnet = 0)

era5_smooth_newdata_RH <- data.table(SAT = 0,
                                  RH = ERA5_RH_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                  Rnet = 0)


era5_smooth_newdata_Rnet <- data.table(SAT = 0,
                                  RH = 0,
                                  Rnet = ERA5_Rnet_JJA_2d50_anom_XAX$X[, land_mask[i]])

era5_smooth_newdata_all <- data.table(SAT = ERA5_T_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                       RH = ERA5_RH_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                       Rnet = ERA5_Rnet_JJA_2d50_anom_XAX$X[, land_mask[i]])


obs_smooth_decomp_list[[i]] <- data.table(pred_smooth_T = predict(lm_list[[i]], era5_smooth_newdata_T),
                                          pred_smooth_RH = predict(lm_list[[i]], era5_smooth_newdata_RH),
                                          pred_smooth_Rnet = predict(lm_list[[i]], era5_smooth_newdata_Rnet),
                                          #pred_all = predict(lm_list[[i]], era5_smooth_newdata_all),
                                          #ET_smooth = ERA5_ET_smooth[,land_mask[i]],
                                          year = ERA5_RH_JJA_2d50_anom_XAX$M$year,
                                          ET = ERA5_ET_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                          obs = "ERA5",
                                          latlonDT[land_mask[i]])

}

obs_smooth_decomp <- rbindlist(obs_smooth_decomp_list) %>% melt(id.vars = c("year", "obs", "lon", "lat", "lonlat", "SREX", "ix"))

obs_smooth_decomp[SREX == "EAS"][, .(value = mean(value)), .(year, SREX, variable)] %>% 
  ggplot()+
  geom_line(aes(year, value, color = variable))

obs_smooth_decomp %>%  dcast(value.var = "value", variable ~ year + obs + lon + lat + lonlat + SREX + ix)[, cor(value, ET), .(variable, lonlat, lat, lon)] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = V1))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~variable)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9") +
  theme(panel.grid.major = element_blank())


obs_smooth_decomp_trends <- obs_smooth_decomp[year %in% 1980:2020, .(trend = lm(value~year)$coefficients[2], lat = lat, lon= lon), .(lonlat, variable, obs)]

obs_smooth_decomp_trends %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = trend))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~variable)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9") +
  theme(panel.grid.major = element_blank())

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20230927_lm_vs_ridge/",
       filename = "MERRA2_smooth_lm.png",
       height = 7, width = 12)


#MERRA
obs_smooth_decomp_list_merra <- list()
for(i in seq_along(land_mask)){
  
  merra_smooth_newdata_T <- data.table(SAT = MERRA_T_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                      RH = 0,
                                      Rnet = 0)
  
  merra_smooth_newdata_RH <- data.table(SAT = 0,
                                       RH = MERRA_RH_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                       Rnet = 0)
  
  merra_smooth_newdata_Rnet <- data.table(SAT = 0,
                                         RH = 0,
                                         Rnet = MERRA_RH_JJA_2d50_anom_XAX$X[, land_mask[i]])
  
  obs_smooth_decomp_list_merra[[i]] <- data.table(pred_T = predict(lm_list[[i]], merra_smooth_newdata_T),
                                            pred_RH = predict(lm_list[[i]], merra_smooth_newdata_RH),
                                            pred_Rnet = predict(lm_list[[i]], merra_smooth_newdata_Rnet),
                                            year = MERRA_RH_JJA_2d50_anom_XAX$M$year,
                                            ET = MERRA_ET_JJA_2d50_anom_XAX$X[, land_mask[i]],
                                            obs = "MERRA2",
                                            latlonDT[land_mask[i]])
  
}

merra2_smooth_decomp <- rbindlist(obs_smooth_decomp_list_merra) %>% melt(id.vars = c("year", "obs", "lon", "lat", "lonlat", "SREX", "ix"))
merra2_smooth_decomp_trends <- merra2_smooth_decomp[year %in% 1980:2020, .(trend = lm(value~year)$coefficients[2], lat = lat, lon= lon), .(lonlat, variable, obs)]

merra2_smooth_decomp_trends %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = trend))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~variable)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9") +
  theme(panel.grid.major = element_blank())


merra2_smooth_decomp[SREX == "NEU"][, .(value = mean(value)), .(year, SREX, variable)] %>% 
  ggplot()+
  geom_line(aes(year, value, color = variable))

#================actual OBS===================
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/CERES/CERES_Rnet_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/CRUTS/CruTS_RH_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/GLEAM/GLEAM_3.6b_JJA_2d50_anom_XAX.RData")

era_ix <- ERA5_Rnet_JJA_2d50_anom_XAX$M$year %in% 2003:2021
merra_ix <- MERRA_Rnet_JJA_2d50_anom_XAX$M$year %in% 2003:2021
cru_ix <- CruTS_RH_JJA_2d50_anom_XAX$M$year %in% 2003:2021
ceres_ix <- CERES_Rnet_JJA_2d50_anom_XAX$M$year %in% 2003:2021
GLEAM_ix <- GLEAM_3.6b_JJA_2d50_anom_XAX$M$year %in% 2003:2021


CruTS_pred_list <- list() 
for (i in seq_along(land_mask)){

cru_merra_newdata <- data.table(
                                RH = ERA5_RH_JJA_2d50_anom_XAX$X[era_ix, land_mask[i]],
                                Rnet = CERES_Rnet_JJA_2d50_anom_XAX$X[ceres_ix, land_mask[i]])

cru_era_newdata <- data.table(
                                    RH = MERRA_RH_JJA_2d50_anom_XAX$X[merra_ix, land_mask[i]],
                                    Rnet = CERES_Rnet_JJA_2d50_anom_XAX$X[ceres_ix, land_mask[i]])

cru_ceres_newdata <- data.table(
                                RH = CruTS_RH_JJA_2d50_anom_XAX$X[cru_ix, land_mask[i]],
                                Rnet = CERES_Rnet_JJA_2d50_anom_XAX$X[ceres_ix, land_mask[i]])


CruTS_pred_list[[i]] <- data.table(cru_merra= predict(lm_list_RHRnet[[i]], cru_merra_newdata),
                   cru_era = predict(lm_list_RHRnet[[i]], cru_era_newdata),
                   cru_ceres = predict(lm_list_RHRnet[[i]], cru_ceres_newdata),
                   ET_ERA5 = ERA5_ET_JJA_2d50_anom_XAX$X[era_ix,land_mask[i]],
                   ET_MERRA2 = MERRA_ET_JJA_2d50_anom_XAX$X[merra_ix,land_mask[i]],
                   ET_GLEAM = GLEAM_3.6b_JJA_2d50_anom_XAX$X[GLEAM_ix, land_mask[i]],
                   year = ERA5_RH_JJA_2d50_anom_XAX$M$year[era_ix],
                   latlonDT[land_mask[i]])
}


rbindlist(CruTS_pred_list) %>% melt(id.vars = c("year", "lat", "lon", "lonlat", "SREX", "ix")) %>% 
  .[, .(value = mean(value)), .(SREX, year, variable)] %>% 
  .[SREX == "SAS"] %>% 
  ggplot()+
  geom_line(aes(year, value, color = variable))



rbindlist(CruTS_pred_list) %>% melt(id.vars = c("year", "lat", "lon", "lonlat", "SREX", "ix", "ET_ERA5", "ET_MERRA2", "ET_GLEAM")) %>% 
  .[, .(cor_ERA5 = cor(value, ET_ERA5), cor_MERRA2 = cor(value, ET_MERRA2), cor_GLEAM= cor(value, ET_GLEAM)), .(ix,lat, lon, variable)] %>% 
  melt(id.vars = c("lat", "lon", "ix", "variable"), variable.name = "eval_data") %>% 
  ggplot() +
  geom_tile(aes(lon, lat, fill = value))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_grid(variable~eval_data)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9") +
  theme(panel.grid.major = element_blank())





#==========fitting LM to every GCM===============
beta_singe_GCM_list_temp <- list()
beta_singe_GCM_list <- list()

for(mod in seq_along(CMIP6_models)){
  
  train_ix <- which(CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$year %in% 1980:2020 & CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$mod == CMIP6_models[mod])
  
  #train_ix <- which(CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$year %in% 1980:2020 &
  #                    CMIP6_rnet_JJA_hist585_anom_2d50_XAX$M$mod == CMIP6_models[mod] & 
  #                    CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$ens.mem == CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$ens.mem[train_ix][1])
  


  for(i in seq_along(land_mask)){
    
    x_data <- data.table(ET = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         RH = CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         Rnet = CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]])

    lm <-       lm(ET ~., data = x_data)
    
    pred = predict(lm, data.table(ET = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                        RH = CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                        Rnet = CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]]))

    beta_singe_GCM_list_temp[[i]] <-  data.table(beta = lm$coefficients[2:3], 
                                                 beta_scaled = lm$coefficients[2:3] * c(sds_hurs[land_mask[i]], 
                                                                                           sds_rnet[land_mask[i]]),
                                                 var = c("RH", "Rnet"),
                                                 latlonDT[land_mask[i]],
                                                 mod = CMIP6_models[mod])
    
  }
  beta_singe_GCM_list[[mod]] <- rbindlist(beta_singe_GCM_list_temp)
}

beta_singe_GCM <- beta_singe_GCM_list %>% rbindlist()
beta_singe_GCM_ratio <- beta_singe_GCM %>% dcast(lat + lon + lonlat +SREX +ix + mod~var, value.var = "beta") %>% 
  .[, ratio := RH/Rnet]

beta_singe_GCM[, .(beta_var = var(beta)), .(lat, lon, lonlat, var)] %>% 
  ggplot() +
  geom_tile(aes(lon, lat, fill = beta_var))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~var)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9") +
  theme(panel.grid.major = element_blank())


beta_singe_GCM[var == "RH"] %>% 
  ggplot() +
  geom_tile(aes(lon, lat, fill = beta))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~mod)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9") +
  theme(panel.grid.major = element_blank())



beta_singe_GCM_ratio %>% 
  ggplot() +
  geom_tile(aes(lon, lat, fill = ratio))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  facet_wrap(~mod)+
  scale_fill_gradient2(high = "#c90202", mid = "white", low ="#0237c9", limits = c(-4,4)) +
  theme(panel.grid.major = element_blank())
#============= hierarcical partitioning ======


library(hier.part)
hierpart_list <- list()
hierpart_table_mods <- list()
CMIP6_models_all <- CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$mod %>% unique()

for (mod in seq_along(CMIP6_models_all)){
  for (i in 1:length(land_mask)){

    train_ix <-which(CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$mod == CMIP6_models_all[mod] &
                       CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$year %in% 1950:2014)
    x_data <- data.table(ET = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         #SWnet = CMIP6_netsw_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         #LWnet = CMIP6_netlw_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         Rnet = CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         RH = CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         SAT = CMIP6_tas_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                         year = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$year[train_ix],
                         mem = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$ens.mem[train_ix])
  
    
    x_data_FR <- x_data[,!"mem"][, lapply(.SD, mean), year][, lapply(.SD, function(x) predict(loess(x ~ year,degree = 2, span = 0.3), year))] 
    temp_FR <- x_data_FR %>% melt(id.vars = "year", value.name = "FR")
    temp_FR[, year := round(year, 4)]
    
    temp_IAV <- melt(x_data, id.vars = c("year", "mem")
                     )[temp_FR, on = c("year", "variable")
                       ][, value_IAV := value-FR
                         ][, .(year, variable, value_IAV, mem)]
    

    
    x_data_IAV <- dcast(temp_IAV, year+mem ~variable, value.var= "value_IAV")
    

    fit_hierpart_IAV <- hier.part(y = x_data_IAV$ET, xcan = x_data_IAV[,.(Rnet, RH, SAT)], barplot = F)
    fit_hierpart_FR <- hier.part(y = x_data_FR$ET, xcan = x_data_FR[,.(Rnet, RH, SAT)], barplot = F)
    
  
  
    hierpart_list[[i]] <- data.table(IAV = fit_hierpart_IAV$I.perc,
                                     FR = fit_hierpart_FR$I.perc,
                                     variable = row.names(fit_hierpart_FR$I.perc),
                                     latlonDT[land_mask[i]],
                                     type = "hier_part",
                                     mod = CMIP6_models_all[mod]) %>% 
      setnames(c("IAV.ind.exp.var", "FR.ind.exp.var"), c("IAV", "FR"))
  }
  hierpart_table_mods[[mod]] <- rbindlist(hierpart_list) 

}

HP_importance_FR <- rbindlist(hierpart_table_mods)[,.SD[which.max(FR)], .(lonlat, lat, lon, mod)
                               ][, .N, .(lonlat, lon, lat, variable)
                                 ][, .SD[which.max(N)], .(lonlat, lon, lat)
                                   ][,eval := "Long term change"
                                     ][, seas := "JJA"]


HP_importance_IAV <- rbindlist(hierpart_table_mods)[,.SD[which.max(IAV)], .(lonlat, lat, lon, mod)
                                                   ][, .N, .(lonlat, lon, lat, variable)
                                                     ][, .SD[which.max(N)], .(lonlat, lon, lat)
                                                       ][,eval := "Short term variability"
                                                         ][, seas := "JJA"]

weights <- get_area_weights()[land_mask]
weights <- weights * 1/sum(weights)
importance_IAV_temp <- rbindlist(hierpart_table_mods)[,.SD[which.max(IAV)], .(lonlat, lat, lon, mod)
                               ][, .N, .(lonlat, lon, lat, variable)
                                 ][, .SD[which.max(N)], .(lonlat, lon, lat)
                                   ]

importance_IAV_bar <- importance_IAV_temp[, weights := weights][, .(area = sum(weights)), variable][,eval := "Short term variability"]
importance_FR_bar <- rbindlist(hierpart_table_mods)[,.SD[which.max(FR)], .(lonlat, lat, lon, mod)
                                                    ][, .N, .(lonlat, lon, lat, variable)
                                                      ][, .SD[which.max(N)], .(lonlat, lon, lat)
                                                        ][, weights := weights
                                                          ][, .(area = sum(weights)), variable
                                                            ][,eval := "Long term change"]







importance_map <- rbind(HP_importance_FR, HP_importance_IAV) %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = variable, alpha = N/19*100))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  scale_fill_manual(values = c("steelblue","goldenrod" ,"forestgreen"), name  = "Most important variable")+
  scale_alpha_continuous(labels = c("40% GCM agreement", "60% GCM agreement", "80% GCM agreement", "100% GCM agreement"), name = "Opacity", )+
  theme_minimal()+
  facet_wrap(~eval, ncol = 1)+
  ylim(-55,80)+
  theme(panel.grid.major = element_blank(), axis.text = element_blank()) 

importance_bars <- rbind(importance_IAV_bar, importance_FR_bar) %>% 
  ggplot(aes(variable))+
  geom_bar(aes(y = area, fill = variable), stat = "identity",)+
  theme_minimal()+
  theme(panel.grid = element_blank(), strip.text = element_blank(), legend.position = "none")+
  facet_wrap(~eval, ncol = 1)+
  xlab("")+ylab("Fraction of area")+
  scale_fill_manual(values = c("steelblue","goldenrod" ,"forestgreen"))
  
ggarrange(importance_bars, importance_map,ncol = 2, widths = c(1,7))

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240502_draft_update2/",
       filename = "FR_vs_IAV_sens_hierpart.png",
       height = 6, width = 12,
       bg = "white")
#=========Shapely===========

shap_val_list <- list()
for (i in seq_along(land_mask)){
  
  x_data_test <- data.table(ET = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                            Rnet = CMIP6_rnet_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                            RH = CMIP6_hurs_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]],
                            SAT = CMIP6_tas_JJA_hist585_anom_2d50_XAX$X[train_ix, land_mask[i]])
  
  
  
  val <- shapleyvalue(x_data_test[,ET], x_data_test[,-1, with = F])
  shap_val_list[[i]] <- data.table(val[2,],
                                   latlonDT[land_mask[i]])
}

shapely_val <- rbindlist(shap_val_list) %>% melt(id.vars = c("lat", "lon", "lonlat", "SREX", "ix"))
shapely_val %>% 
  .[,.SD[which.max(value)], .(lonlat, lat, lon)] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = variable))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylim(-52,80)+
  ylab("") + xlab("") +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())


shapely_val[,eval_type := "Shapely"]
hierpart_table[,eval_type := "HierPart"]

rbind(shapely_val, hierpart_table) %>% 
  .[,.SD[which.max(value)], .(lonlat, lat, lon, type)] %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = variable))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylim(-52,80)+
  facet_wrap(~eval_type)+
  ylab("") + xlab("") +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())


#==========VIF=============

vif_rnet <- lapply(lm_list_RHRnetSAT, car::vif)
vif_swlw <- lapply(lm_list, car::vif)

rbindlist(lapply(vif_swlw, function(x) data.table(values = x, variable = names(x))))

vif_rnet_list <- list()
vif_lwsw_list <- list()
for (i in seq_along(land_mask)){
  vif_rnet_list[[i]] <- data.table(RH = vif_rnet[[i]][1],
           Rnet =  vif_rnet[[i]][2],
           SAT =  vif_rnet[[i]][3],
           latlonDT[land_mask][i])
  
  vif_lwsw_list[[i]] <- data.table(SWnet = vif_swlw[[i]][1],
                                   LWnet =  vif_swlw[[i]][2],
                                   SAT =  vif_swlw[[i]][3],
                                   RH =  vif_swlw[[i]][4],
                                   latlonDT[land_mask][i])
  
}
vif_rnet <- rbindlist(vif_rnet_list) %>% melt(id.vars = c("lat", "lon", "lonlat", "SREX", "ix"))

vif_rnet_summary <- vif_rnet[value > 5, .(N= round(.N/1799*100, 1)), variable][, lat := -40][, lon := -120]

ggplot(vif_rnet)+
  geom_tile(aes(lon, lat, fill = value))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylim(-55,80)+
  geom_text(data = vif_rnet_summary, aes(lon,lat, label = paste0(N, "%")))+
  facet_wrap(~variable, ncol = 1)+
  scale_fill_viridis_c(name = "VIF")+
  theme_minimal()+
  ylab("") + xlab("") +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240502_draft_update2/",
       filename = "SATRnetRH_VIF.png",
       height = 10, width = 8,
       bg = "white")



rbindlist(vif_lwsw_list) %>% melt(id.vars = c("lat", "lon", "lonlat", "SREX", "ix")) %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = value))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylim(-55,80)+
  facet_wrap(~variable)+
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())


