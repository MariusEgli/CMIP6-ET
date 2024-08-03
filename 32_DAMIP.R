load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netlw_JJA_histGHG_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netlw_JJA_histAER_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netsw_JJA_histGHG_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_netsw_JJA_histAER_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hfls_JJA_histAER_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hfls_JJA_histGHG_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hurs_JJA_histAER_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_hurs_JJA_histGHG_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_tas_JJA_histAER_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_tas_JJA_histGHG_2d50_XAX.RData")


SREX_sel <- c("WCE", "MED", "EAS", "SAS", "WNA", "NSA")


histGHG_var_sens_list <- list()

for(i in seq_along(land_mask)){
  histghg_newdata_rh <- data.table(RH = CMIP6_hurs_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                   SAT = 0,
                                   LWnet = 0,
                                   SWnet = 0)
  
  histghg_newdata_swd <- data.table(SWnet = CMIP6_netsw_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                    SAT = 0,
                                    LWnet = 0,
                                    RH = 0)
  
  histghg_newdata_lwd <- data.table(SWnet = 0,
                                    SAT = 0,
                                    LWnet = CMIP6_netlw_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                    RH = 0)
  
  histghg_newdata_tas <- data.table(SWnet = 0,
                                    SAT = CMIP6_tas_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                    LWnet = 0,
                                    RH = 0)
  
  histGHG_var_sens_list[[i]] <- data.table(RH = predict(lm_list[[i]], histghg_newdata_rh),
                                           SWnet = predict(lm_list[[i]], histghg_newdata_swd),
                                           LWnet = predict(lm_list[[i]], histghg_newdata_lwd),
                                           SAT = predict(lm_list[[i]], histghg_newdata_tas),
                                           ET = CMIP6_hfls_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                           year = CMIP6_hfls_JJA_histGHG_2d50_XAX$M$year,
                                           mem = CMIP6_hfls_JJA_histGHG_2d50_XAX$M$ens.mem,
                                           mod = CMIP6_hfls_JJA_histGHG_2d50_XAX$M$mod,
                                           latlonDT[land_mask[i]])
  
}

histGHG_var_sens <- rbindlist(histGHG_var_sens_list) %>% melt(id.vars= c("year", "mod", "lat", "lon", "lonlat", "SREX", "mem", "ix"))

histGHG_var_sens[, ens.mean := mean(value), .(year, mod, lat, lon, lonlat, SREX, variable)]
histGHG_var_sens[, value_smooth := loess(ens.mean ~ year, degree = 2, span = 0.5 )$fitted, .(mod, lat, lon, lonlat, SREX, variable, mem)]




histGHG_var_sens[SREX  %in% SREX_sel
                 ][,.(value_smooth = mean(value_smooth)), .(mod, year, variable, SREX)] %>% 
  ggplot()+
  geom_line(aes(year, value_smooth, color = variable, group = mod))+ 
  facet_grid(variable~SREX)+
  scale_color_manual(values = c("steelblue", "goldenrod","tomato2", "forestgreen",  "grey30"), name  = "", guide = guide_legend(reverse = TRUE)) +
  theme_minimal()+
  ylab("ET anomaly [mm/day]") + xlab("Year")




histAER_var_sens_list <- list()

for(i in seq_along(land_mask)){
  histAER_newdata_rh <- data.table(RH = CMIP6_hurs_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                   SAT = 0,
                                   LWnet = 0,
                                   SWnet = 0)
  
  histAER_newdata_swd <- data.table(SWnet = CMIP6_netsw_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                    SAT = 0,
                                    LWnet = 0,
                                    RH = 0)
  
  histAER_newdata_lwd <- data.table(SWnet = 0,
                                    SAT = 0,
                                    LWnet = CMIP6_netlw_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                    RH = 0)
  
  histAER_newdata_tas <- data.table(SWnet = 0,
                                    SAT = CMIP6_tas_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                    LWnet = 0,
                                    RH = 0)
  
  histAER_var_sens_list[[i]] <- data.table(RH = predict(lm_list[[i]], histAER_newdata_rh),
                                           SWnet = predict(lm_list[[i]], histAER_newdata_swd),
                                           LWnet = predict(lm_list[[i]], histAER_newdata_lwd),
                                           SAT = predict(lm_list[[i]], histAER_newdata_tas),
                                           ET = CMIP6_hfls_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                           year = CMIP6_hfls_JJA_histAER_2d50_XAX$M$year,
                                           mem = CMIP6_hfls_JJA_histAER_2d50_XAX$M$ens.mem,
                                           mod = CMIP6_hfls_JJA_histAER_2d50_XAX$M$mod,
                                           latlonDT[land_mask[i]])
  
}

histAER_var_sens <- rbindlist(histAER_var_sens_list) %>% melt(id.vars= c("year", "mod", "lat", "lon", "lonlat", "SREX", "mem", "ix"))

histAER_var_sens[, ens.mean := mean(value), .(year, mod, lat, lon, lonlat, SREX, variable)]
histAER_var_sens[, value_smooth := loess(ens.mean ~ year, degree = 2, span = 0.5 )$fitted, .(mod, lat, lon, lonlat, SREX, variable, mem)]

histAER_var_sens[SREX   %in% SREX_sel
                 ][,.(value_smooth = mean(value_smooth)), .(year, variable, SREX, mod)] %>% 
  ggplot()+
  geom_line(aes(year, value_smooth, color = variable, group = mod))+ 
  facet_grid(variable~SREX) +
  scale_color_manual(values = c("steelblue", "goldenrod","tomato2", "forestgreen",  "grey30"), name  = "") +
  theme_minimal() +
  ylab("ET anomaly [mm/day]") + xlab("Year")+
  xlim(1950, 2020)


histAER_var_sens_trends <- histAER_var_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, SREX, mod)
                 ][,.(histAER = lm(value~year)$coefficients[2]), .(mod, SREX, variable)]

histGHG_var_sens_trend <- histGHG_var_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, SREX, mod ,mem)
                                           ][ ,.(value = mean(value), .N), .(year, variable, SREX, mod )
                                              ][,.(histGHG = lm(value~year)$coefficients[2]), .(mod, SREX, variable)]



histGHG_var_sens[mod == "GISS-E2-1-G"]
histGHG_var_sens[year %in% 1950:2014, .(value = mean(value)), .(mod ,mem)
                 ][ ,.(value = mean(value), .N), .(mod )]


hist585_var_sens_trend <- hist585_var_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, SREX, mod)
                 ][,.(histAll = lm(value~year)$coefficients[2]), .(mod, SREX, variable)]

DAMIP_trends <- hist585_var_sens_trend[histGHG_var_sens_trend, on = .(mod, SREX, variable)
                                       ][histAER_var_sens_trends, on = .(mod, SREX, variable)
                                         ][SREX %in% SREX_sel]

histGHG_cor <- DAMIP_trends[!is.na(histAll), .(round(cor(histAll, histGHG, method = "spearman"), 2)), .(SREX)][,x := -0.005][,y:= 0.005]
histAER_cor <- DAMIP_trends[!is.na(histAll), .(round(cor(histAll, histAER, method = "spearman"), 2)), .(SREX)][,x := -0.005][,y:= 0.005]

load("/net/xenon/climphys_backedup/maegli/Misc/SREX/IPCC-WGI-reference-regions-v4_R.rda")


SREX_name_table <- data.table(SREX_long = IPCC_WGI_reference_regions_v4$Name,
           SREX=IPCC_WGI_reference_regions_v4$Acronym)[SREX %in% SREX_sel]

SREX_name_table[, SREX_long := c("Western North-America", "")]

p1 <- ggplot(DAMIP_trends)+
  geom_hline(yintercept = 0, linewidth = 0.2)+
  geom_vline(xintercept = 0, linewidth = 0.2)+
  geom_point(aes(histAll, histGHG, color = variable))+
  geom_abline(intercept = 0, slope = 1, color = "grey20", linetype = 3)+
  facet_wrap(~SREX)+
  scale_color_manual(values = c("steelblue", "goldenrod","tomato2", "forestgreen",  "grey30"), name  = "Contribution") +
  geom_text(data = histGHG_cor, aes(x,y, label = paste0("r = ", V1)))+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  coord_equal()+
  ylim(-0.015,0.006)+
  labs(
    x = bquote(HistAll~ET~trend~ "["~mm ~day^-1~year^-1~"]"),
    y = bquote(HistGHG~ET~trend~ "["~mm ~day^-1~year^-1~"]"))



p2 <- ggplot(DAMIP_trends)+
  geom_hline(yintercept = 0, linewidth = 0.2)+
  geom_vline(xintercept = 0, linewidth = 0.2)+
  geom_point(aes(histAll, histAER, color = variable))+
  geom_abline(intercept = 0, slope = 1, color = "grey20", linetype = 3)+
  facet_wrap(~SREX)+
  scale_color_manual(values = c("steelblue", "goldenrod","tomato2", "forestgreen",  "grey30"), name  = "") +
  geom_text(data = histAER_cor, aes(x,y, label = paste0("r = ", V1)))+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  ylim(-0.015,0.006)+
  coord_equal()+
  labs(
    x = bquote(HistAll~ET~trend~ "["~mm ~day^-1~year^-1~"]"),
    y = bquote(HistAER~ET~trend~ "["~mm ~day^-1~year^-1~"]"))


library(ggpubr)
ggarrange(p1, p2, nrow = 2)
ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240227_draft_update/",
       filename = "damip_trends.png",
       height = 12, width = 8,
       bg = "white")


DAMPI_trends[, unique(paste(mod, mem))]

rbind(histAER_var_sens_trends, histGHG_var_sens_trend, hist585_var_sens_trend)[SREX %in% SREX_sel] %>% 
  ggplot()+
  geom_jitter(aes(scen, trend, fill = variable), width = 0.1)+
  facet_grid(variable~SREX) +
  scale_fill_manual(values = c("steelblue", "goldenrod","tomato2",  "grey30"), name  = "") 
  
  


table(CMIP6_netlw_JJA_histGHG_2d50_XAX$M$mod) /171
CMIP6_netlw_JJA_histGHG_2d50_XAX$M$year %in% 1950:2014 %>% 

CMIP6_netlw_JJA_histGHG_2d50_XAX$M$year[CMIP6_netlw_JJA_histGHG_2d50_XAX$M$mod == "CESM2"]




#=======Rnet====================
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_rnet_JJA_histGHG_2d50_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/anom/CMIP6_rnet_JJA_histAER_2d50_XAX.RData")


histGHG_Rnet_sens_list <- list()

for(i in seq_along(land_mask)){
  histghg_newdata_rh <- data.table(RH = CMIP6_hurs_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                   SAT = 0,
                                   Rnet = 0)
  
  histghg_newdata_rnet <- data.table(Rnet = CMIP6_rnet_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                    SAT = 0,
                                    RH = 0)
  
  histghg_newdata_tas <- data.table(Rnet = 0,
                                    SAT = CMIP6_tas_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                    RH = 0)
  
  histGHG_Rnet_sens_list[[i]] <- data.table(RH = predict(lm_list_RHRnetSAT[[i]], histghg_newdata_rh),
                                           Rnet = predict(lm_list_RHRnetSAT[[i]], histghg_newdata_rnet),
                                           SAT = predict(lm_list_RHRnetSAT[[i]], histghg_newdata_tas),
                                           ET = CMIP6_hfls_JJA_histGHG_2d50_XAX$X[, land_mask[i]],
                                           year = CMIP6_hfls_JJA_histGHG_2d50_XAX$M$year,
                                           mem = CMIP6_hfls_JJA_histGHG_2d50_XAX$M$ens.mem,
                                           mod = CMIP6_hfls_JJA_histGHG_2d50_XAX$M$mod,
                                           latlonDT[land_mask[i]])
  
}

histGHG_Rnet_sens <- rbindlist(histGHG_Rnet_sens_list) %>% melt(id.vars= c("year", "mod", "lat", "lon", "lonlat", "SREX", "mem", "ix"))
histGHG_Rnet_sens[, ens.mean := mean(value), .(year, mod, lat, lon, lonlat, SREX, variable)]


histAER_Rnet_sens_list <- list()

for(i in seq_along(land_mask)){
  histAER_newdata_rh <- data.table(RH = CMIP6_hurs_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                   SAT = 0,
                                   Rnet = 0)
  
  histAER_newdata_rnet <- data.table(Rnet = CMIP6_rnet_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                     SAT = 0,
                                     RH = 0)
  
  histAER_newdata_tas <- data.table(Rnet = 0,
                                    SAT = CMIP6_tas_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                    RH = 0)
  
  histAER_Rnet_sens_list[[i]] <- data.table(RH = predict(lm_list_RHRnetSAT[[i]], histAER_newdata_rh),
                                            Rnet = predict(lm_list_RHRnetSAT[[i]], histAER_newdata_rnet),
                                            SAT = predict(lm_list_RHRnetSAT[[i]], histAER_newdata_tas),
                                            ET = CMIP6_hfls_JJA_histAER_2d50_XAX$X[, land_mask[i]],
                                            year = CMIP6_hfls_JJA_histAER_2d50_XAX$M$year,
                                            mem = CMIP6_hfls_JJA_histAER_2d50_XAX$M$ens.mem,
                                            mod = CMIP6_hfls_JJA_histAER_2d50_XAX$M$mod,
                                            latlonDT[land_mask[i]])
  
}

histAER_Rnet_sens <- rbindlist(histAER_Rnet_sens_list) %>% melt(id.vars= c("year", "mod", "lat", "lon", "lonlat", "SREX", "mem", "ix"))
histAER_Rnet_sens[, ens.mean := mean(value), .(year, mod, lat, lon, lonlat, SREX, variable)]


histAER_Rnet_sens_trends <- histAER_Rnet_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, SREX, mod)
                                            ][,.(histAER = lm(value~year)$coefficients[2]), .(mod, SREX, variable)]

histGHG_Rnet_sens_trend <- histGHG_Rnet_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, SREX, mod)
                                           ][ ,.(value = mean(value), .N), .(year, variable, SREX, mod )
                                              ][,.(histGHG = lm(value~year)$coefficients[2]), .(mod, SREX, variable)]



hist585_Rnet_sens_trend <- hist585_RHRnetSAT_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, SREX, mod)
                                           ][,.(histAll = lm(value~year)$coefficients[2]), .(mod, SREX, variable)]


DAMIP_trends_rnet <- hist585_Rnet_sens_trend[histGHG_Rnet_sens_trend, on = .(mod, SREX, variable)
][histAER_Rnet_sens_trends, on = .(mod, SREX, variable)
][SREX %in% SREX_sel]

histGHG_cor_rnet <- DAMIP_trends_rnet[mod != "FGOALS-g3"
                                      ][mod != "GISS-E2-1-G"
                                        #][variable == "ET"
                                          ][!is.na(histAll), .(round(cor(histAll, histGHG, method = "spearman"), 2)), .(SREX)
                                            ][,x := -0.005][,y:= 0.005]
histAER_cor_rnet <- DAMIP_trends_rnet[mod != "FGOALS-g3"
                                      ][mod != "GISS-E2-1-G"
                                        #][variable == "ET"
                                          ][!is.na(histAll), .(round(cor(histAll, histAER, method = "spearman"), 2)), .(SREX)
                                            ][,x := -0.005][,y:= 0.005]
histBoth_cor_rnet <- DAMIP_trends_rnet[mod != "FGOALS-g3"][mod != "GISS-E2-1-G"][!is.na(histAll), .(round(cor(histAll, histAER + histGHG, method = "spearman"), 2)), .(SREX)][,x := -0.005][,y:= 0.005]


DAMIP_trends_rnet[variable == "Rnet" & histGHG < -0.005]

p1_rnet <- ggplot(DAMIP_trends_rnet[mod != "FGOALS-g3"][mod != "GISS-E2-1-G"])+
  geom_hline(yintercept = 0, linewidth = 0.2)+
  geom_vline(xintercept = 0, linewidth = 0.2)+
  geom_point(aes(histAll, histGHG, color = variable))+
  geom_abline(intercept = 0, slope = 1, color = "grey20", linetype = 3)+
  facet_wrap(~SREX)+
  scale_color_manual(values = c("dodgerblue", "goldenrod","forestgreen", "grey30"), name  = "Contribution",
                     labels = c(bquote(ET[RH]), bquote(ET[R[net]]), bquote(ET[SAT]), "ET")) +
  geom_text(data = histGHG_cor_rnet, aes(x,y, label = paste0("r = ", V1)))+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  coord_equal()+
  ylim(-0.008,0.008)+
  labs(
    x = bquote(HistAll~ET~trend~ "["~mm ~day^-1~year^-1~"]"),
    y = bquote(HistGHG~ET~trend~ "["~mm ~day^-1~year^-1~"]"))


p2_rnet <- ggplot(DAMIP_trends_rnet[mod != "FGOALS-g3"][mod != "GISS-E2-1-G"])+
  geom_hline(yintercept = 0, linewidth = 0.2)+
  geom_vline(xintercept = 0, linewidth = 0.2)+
  geom_point(aes(histAll, histAER, color = variable))+
  geom_abline(intercept = 0, slope = 1, color = "grey20", linetype = 3)+
  facet_wrap(~SREX)+
  scale_color_manual(values = c("dodgerblue", "goldenrod","forestgreen", "grey30"), name  = "Contribution", 
                     labels = c(bquote(ET[RH]), bquote(ET[R[net]]), bquote(ET[SAT]), "ET")) +
  geom_text(data = histAER_cor_rnet, aes(x,y, label = paste0("r = ", V1)))+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  ylim(-0.008,0.008)+
  coord_equal()+
  labs(
    x = bquote(HistAll~ET~trend~ "["~mm ~day^-1~year^-1~"]"),
    y = bquote(HistAER~ET~trend~ "["~mm ~day^-1~year^-1~"]"))


ggarrange(p1_rnet, p2_rnet, nrow = 2)

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision/",
       filename = "damip_trends_rnet.png",
       height = 12, width = 8,
       bg = "white")

ggplot(DAMIP_trends_rnet[mod != "FGOALS-g3"][mod != "GISS-E2-1-G"])+
  geom_hline(yintercept = 0, linewidth = 0.2)+
  geom_vline(xintercept = 0, linewidth = 0.2)+
  geom_point(aes(histAll, histAER + histGHG, color = variable))+
  geom_abline(intercept = 0, slope = 1, color = "grey20", linetype = 3)+
  facet_wrap(~SREX)+
  scale_color_manual(values = c("steelblue", "goldenrod","forestgreen", "grey30"), name  = "Contribution") +
  geom_text(data = histBoth_cor_rnet, aes(x,y, label = paste0("r = ", V1)))+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  ylim(-0.008,0.008)+
  coord_equal()+
  labs(
    x = bquote(HistAll~ET~trend~ "["~mm ~day^-1~year^-1~"]"),
    y = bquote(HistAER+HistGHG~ET~trend~ "["~mm ~day^-1~year^-1~"]"))


histAER_Rnet_sens_trends_lonlat <- histAER_Rnet_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, mod, lon, lat)
                                                     ][,.(histAER = lm(value~year)$coefficients[2]), .(mod, lon, lat, variable)]

histGHG_Rnet_sens_trend_lonlat <- histGHG_Rnet_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, mod, lon, lat)
                                             ][,.(histGHG = lm(value~year)$coefficients[2]), .(mod, lat, lon, variable)]

hist585_Rnet_sens_trend_lonlat <- hist585_RHRnetSAT_sens[year %in% 1950:2014, .(value = mean(value)), .(year, variable, lon, lat, mod)
                                                         ][,.(histAll = lm(value~year)$coefficients[2]), .(mod, lon, lat, variable)]


DAMIP_trends_rnet_lonlat <- histAER_Rnet_sens_trends_lonlat[histGHG_Rnet_sens_trend_lonlat, on = .(mod, lon, lat, variable)
                                ][hist585_Rnet_sens_trend_lonlat, on = .(mod, lon, lat, variable)][!is.na(histGHG)]


scen_lonlat_mod <- DAMIP_trends_rnet_lonlat[mod != "FGOALS-g3"][mod != "GISS-E2-1-G"] %>% melt(id.vars = c("mod", "lon", "lat", "variable", "histAll"), variable.name = "scen") %>% 
  .[, .(value = mse(histAll, value)), .(mod, lon, lat, scen)] %>% 
  .[, .SD[which.min(value)], .(mod, lat, lon)]


ggplot(scen_lonlat_mod)+
  geom_tile(aes(lon, lat, fill = scen))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  facet_wrap(~mod)+
  scale_fill_manual(values= c("dodgerblue2", "tomato2"), name = "Scenario")+
  theme(panel.grid.major = element_blank()) +
  ylim(-55,80)


scen_lonlat_mod_summary <- scen_lonlat_mod[, .N, .(lon, lat, scen)][, .SD[which.max(N)], .(lon, lat)]


dominant_variabel_map <- ggplot(scen_lonlat_mod_summary)+
  geom_tile(aes(lon, lat, fill = scen, alpha = N))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  theme_minimal()+
  scale_fill_manual(values= c("dodgerblue2", "tomato2"), name = "Scenario")+
  theme(panel.grid.major = element_blank(),
        axis.text =element_blank(), 
        axis.ticks=element_blank())+
  ylim(-55,80)+ 
  scale_alpha_continuous(labels = c("60% GCM agreement", "70% GCM agreement", "85% GCM agreement", "100% GCM agreement"), name = "Opacity")
  

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision/",
       filename = "damip_map.png",
       height = 6, width = 10,
       bg = "white")


weights <- get_area_weights()[land_mask]


scen_lonlat_mod[, weights := rep(weights, length(unique(mod)))]


hist585_var_sens_sum <- hist585_RHRnetSAT_sens[variable == "ET"
                                         ][mod != "FGOALS-g3"
                                           ][mod != "GISS-E2-1-G"
                                             ][mod %in% unique(DAMIP_trends_rnet_lonlat$mod)
                                               ][, .(ens.mean = mean(value)), .(mod, year, lat, lon, variable)]

hist585_var_sens_sum[, value_smooth := predict(loess(ens.mean ~ year,degree = 2, span = 0.3), year), .(mod, lat, lon, variable)]

reg_ET_variance <- hist585_var_sens_sum[year %in% 1950:2014, .(var = var(value_smooth)), .(mod, lat, lon)
                     ][scen_lonlat_mod, on = c("lat", "lon", "mod")]



reg_ET_variance[, total_var := sum(weights * var), mod]
reg_ET_variance[, .(var = sum(weights * var), total_var), .(scen, mod)] %>% unique() %>% 
  ggplot()+
  geom_bar(aes (y = var/total_var, x = mod, fill = scen), stat='identity')+
  #geom_text(aes(x = scen, y = var/total_var, label = mod, color = scen), alpha = 0.6) +
  #geom_point(aes(x = scen, y = var/total_var, fill = scen, color = scen), alpha = 0.6) +
  scale_fill_manual(values= c("dodgerblue2", "tomato2"), name = "Scenario")+
  scale_color_manual(values= c("dodgerblue2", "tomato2"), name = "Scenario")+
  ylab("Fraction of total ET variance over land") + xlab("")+
  theme_minimal()+
  theme(panel.grid = element_blank())

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision/",
       filename = "total_ET_fraction_DAMIP.png",
       height = 4, width = 8,
       bg = "white")


  reg_ET_variance_sum <- reg_ET_variance[, .(var = sum(weights * var)), .(scen, mod)]
 
dominat_variable_bars <- dcast(reg_ET_variance_sum, mod ~scen,  value.var = "var") %>%  .[, tot_var := histAER+histGHG] %>% melt(id.vars = c("mod", "tot_var")) %>% 
    ggplot()+
    geom_bar(aes(x = mod, y = value/tot_var, fill= variable), stat = "identity", alpha = 0.9) +
    theme_minimal() +
    theme(panel.grid = element_blank(), legend.position = "none")+
    scale_fill_manual(values= c("dodgerblue2", "tomato2"), name = "Scenario") +
    ylab("Fraction of total long-term ET varicance")+ xlab("") + 
  coord_flip()
    
ggarrange(dominat_variable_bars, dominant_variabel_map, nrow = 1, widths = c(1,3))

ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision/",
       filename = "damip_map_bars.png",
       height = 3.5, width = 12,
       bg = "white")



#=======abs ET===========

load("/net/xenon/climphys_backedup/maegli/CMIP6/data/XAX/seas/orig/CMIP6_hfss_JJA_hist585_2d50_XAX.RData")

hfls_FR <- get_FR_mat(CMIP6_hfss_JJA_hist585_2d50_XAX)

abs_ET_FR <- data.table(hfls_FR$X[, land_mask], 
                        year = hfls_FR$M$year,
                        mod = hfls_FR$M$mod) %>% 
  melt(id.vars = c("mod", "year"))

abs_ET_FR[, ix := as.integer(substr(variable, 2,5))][,variable := NULL]
abs_ET_FR_table <- latlonDT[abs_ET_FR, on = "ix"]

abs_ET_FR_table_scen <- abs_ET_FR_table[mod %in% unique(DAMIP_trends_rnet_lonlat$mod)
                ][mod != "FGOALS-g3"
                  ][mod != "GISS-E2-1-G"
                    ][year %in% 1950:2014, .(value = mean(value)), .(mod, lat, lon, lonlat)][scen_lonlat_mod[, 1:4], on = c("mod", "lat", "lon")]

abs_ET_FR_table_scen[!is.na(value), .(mean(value)),.(scen, mod)][, V1 := V1/sum(V1), (mod)] %>% 
  ggplot+
  geom_bar(aes(x= mod, y = V1, fill = scen), stat = "identity")+    
  scale_fill_manual(values= c("dodgerblue2", "tomato2"), name = "Scenario")+
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("") +ylab("Fraction of total 1950-2014 land ET")
  
  ggsave(path = "/net/xenon/climphys_backedup/maegli/ET_Adj/Figures/20240711_revision/",
         filename = "damip_total_ET_bars.png",
         height = 5, width = 8,
         bg = "white")

