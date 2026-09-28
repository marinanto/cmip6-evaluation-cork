# =============================================================================
# County Cork CMIP6 temperature analysis
# Reproducible manuscript workflow
# =============================================================================
# Reference period: 1985-2014
# Validation period: 2015-2025
# GCMs: CNRM, HadGEM, EC-Earth
# Scenarios: SSP2-4.5, SSP5-8.5
# Signed errors: reference - GCM
#   Gridded evaluation: ERA5-Land - GCM
#   Station validation: station observation - GCM
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
required <- c("loadeR", "transformeR", "downscaleR", "climate4R.climdex",
              "climate4R.UDG", "ncdf4", "readxl", "dplyr", "tidyr",
              "lubridate", "ggplot2", "RColorBrewer", "openair", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(loadeR); library(transformeR); library(downscaleR)
  library(climate4R.climdex); library(climate4R.UDG); library(ncdf4)
  library(readxl); library(dplyr); library(tidyr); library(lubridate)
  library(ggplot2); library(RColorBrewer); library(openair); library(scales)
})
options(stringsAsFactors = FALSE)

# ---- 1. Paths and settings --------------------------------------------------
repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_root, "data")
fig_dir <- file.path(repo_root, "outputs", "figures")
tab_dir <- file.path(repo_root, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

ref_years <- 1985:2014
future_years <- 2015:2025
cork_lon <- -8.486
cork_lat <- 51.847
colours <- rev(brewer.pal(4, "Spectral"))

files <- list(
  era5 = file.path(data_dir, "new_eraland", "era5_land_cork_1985_2025_final.nc"),
  cnrm_ref = file.path(data_dir, "new_eraland", "cnrm_ref_final.nc"),
  hadgem_ref = file.path(data_dir, "new_eraland", "hadgem_ref_final.nc"),
  ec_earth_ref = file.path(data_dir, "new_eraland", "ecearth_ref_final.nc"),
  cnrm45 = file.path(data_dir, "new_eraland", "cnrm45_final.nc"),
  cnrm85 = file.path(data_dir, "new_eraland", "cnrm85_final.nc")
)

stations <- tibble(
  station = c("MOORE PARK", "SHERKIN ISLAND", "ROCHES POINT", "CORK AIRPORT"),
  lon = c(-8.264, -9.428, -8.244, cork_lon),
  lat = c(52.164, 51.476, 51.793, cork_lat),
  file = file.path(data_dir, "stations", c("mly575.xlsx", "mly775.xlsx", "mly1075.xlsx", "mly3904.xlsx"))
)
required_files <- c(unlist(files), stations$file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing input files:\n", paste(missing_files, collapse = "\n"))

# ---- 2. Data loading --------------------------------------------------------
load_grid <- function(file, var, years, lonLim = c(-13, -3), latLim = c(49, 58)) {
  x <- loadGridData(file, var = var, lonLim = lonLim, latLim = latLim, years = years)
  x$Dates$start <- as.POSIXct(format(as.Date(x$Dates$start), "%Y-%m-01"), tz = "UTC")
  x$Dates$end <- as.POSIXct(format(as.Date(x$Dates$end), "%Y-%m-01"), tz = "UTC")
  x
}
era5_land_9km_mon <- load_grid(files$era5, "t2m", 1985:2025)
cnrm_9km_ref_mon <- load_grid(files$cnrm_ref, "tas", ref_years)
hadgem_9km_ref_mon <- load_grid(files$hadgem_ref, "tas", ref_years)
ec_earth_9km_ref_mon <- load_grid(files$ec_earth_ref, "tas", ref_years)
cnrm45_9km <- load_grid(files$cnrm45, "tas", 2015:2100)
cnrm85_9km <- load_grid(files$cnrm85, "tas", 2015:2100)

# ---- 3. Common functions ----------------------------------------------------
spatial_monthly <- function(x) {
  d <- as.Date(x$Dates$start); a <- x$Data
  tibble(date = d, year = year(d), month = month(d),
         temp = vapply(seq_len(dim(a)[1]), function(i) mean(as.numeric(a[i,,]), na.rm = TRUE), numeric(1)))
}
annual_mean <- function(x) spatial_monthly(x) %>% group_by(year) %>% summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")
compute_metrics <- function(obs, sim) {
  ok <- complete.cases(obs, sim); obs <- obs[ok]; sim <- sim[ok]; err <- obs - sim
  q05 <- quantile(obs, .05, na.rm = TRUE); q95 <- quantile(obs, .95, na.rm = TRUE)
  tibble(RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)), MBE = mean(err),
         Correlation = cor(obs, sim), R_squared = cor(obs, sim)^2,
         Cold_bias = mean(err[obs <= q05], na.rm = TRUE),
         Warm_bias = mean(err[obs >= q95], na.rm = TRUE),
         Extreme_bias = mean(abs(c(mean(err[obs <= q05], na.rm = TRUE), mean(err[obs >= q95], na.rm = TRUE))), na.rm = TRUE))
}
get_season_year <- function(date) {
  m <- month(date); y <- year(date)
  tibble(season = case_when(m %in% c(12,1,2) ~ "DJF", m %in% 3:5 ~ "MAM", m %in% 6:8 ~ "JJA", TRUE ~ "SON"),
         season_year = if_else(m == 12, y + 1L, y))
}
nearest_idx <- function(x, lon, lat) {
  g <- expand.grid(lon = x$xyCoords$x, lat = x$xyCoords$y)
  which.min((g$lon-lon)^2 + (g$lat-lat)^2)
}
extract_point <- function(x, lon, lat) {
  i <- nearest_idx(x, lon, lat); nx <- length(x$xyCoords$x)
  li <- ((i-1) %% nx) + 1; la <- ((i-1) %/% nx) + 1; d <- as.Date(x$Dates$start)
  tibble(date=d, year=year(d), month=month(d), temp=as.numeric(x$Data[,la,li]))
}
annualise <- function(df) df %>% mutate(year=year(date)) %>% group_by(year, dataset) %>% summarise(temp=mean(temp, na.rm=TRUE), .groups="drop")
load_station <- function(file) {
  read_excel(file) %>% rename(temp=meant) %>% filter(year >= min(future_years), year <= max(future_years)) %>%
    mutate(date=as.Date(sprintf("%04d-%02d-01", year, month))) %>% select(date,year,month,temp)
}
raw_point <- function(x, lon, lat) extract_point(x,lon,lat) %>% transmute(date,year,month,temp,dataset="CNRM raw") %>% filter(year %in% future_years)
bc_point <- function(pred, base, lon, lat) {
  i <- nearest_idx(base,lon,lat); nx <- length(base$xyCoords$x); li <- ((i-1) %% nx)+1; la <- ((i-1) %/% nx)+1; d <- as.Date(base$Dates$start)
  tibble(date=d,year=year(d),month=month(d),temp=as.numeric(pred$Data[,la,li]),dataset="CNRM bias-corrected") %>% filter(year %in% future_years)
}

# ---- 4. Baseline annual/monthly evaluation ---------------------------------
era5_a <- annual_mean(era5_land_9km_mon) %>% filter(year %in% ref_years)
cnrm_a <- annual_mean(cnrm_9km_ref_mon) %>% filter(year %in% ref_years)
hadgem_a <- annual_mean(hadgem_9km_ref_mon) %>% filter(year %in% ref_years)
ec_a <- annual_mean(ec_earth_9km_ref_mon) %>% filter(year %in% ref_years)

annual_metrics <- bind_rows(
  compute_metrics(era5_a$temp, cnrm_a$temp) %>% mutate(Model="CNRM"),
  compute_metrics(era5_a$temp, hadgem_a$temp) %>% mutate(Model="HadGEM"),
  compute_metrics(era5_a$temp, ec_a$temp) %>% mutate(Model="EC-Earth")
)

m_era5 <- spatial_monthly(era5_land_9km_mon) %>% filter(year %in% ref_years)
m_cnrm <- spatial_monthly(cnrm_9km_ref_mon) %>% filter(year %in% ref_years)
m_hadgem <- spatial_monthly(hadgem_9km_ref_mon) %>% filter(year %in% ref_years)
m_ec <- spatial_monthly(ec_earth_9km_ref_mon) %>% filter(year %in% ref_years)
monthly_metrics <- bind_rows(
  compute_metrics(m_era5$temp, m_cnrm$temp) %>% mutate(Model="CNRM"),
  compute_metrics(m_era5$temp, m_hadgem$temp) %>% mutate(Model="HadGEM"),
  compute_metrics(m_era5$temp, m_ec$temp) %>% mutate(Model="EC-Earth")
)
write.csv(annual_metrics,file.path(tab_dir,"Table_2_annual_metrics.csv"),row.names=FALSE)
write.csv(monthly_metrics,file.path(tab_dir,"Table_3_monthly_metrics.csv"),row.names=FALSE)

# ---- 5. Figures 2 and 3 -----------------------------------------------------
base_plot <- function(df, annual=FALSE) {
  ds <- if (annual) bind_rows(era5_a %>% mutate(dataset="ERA5-Land"),cnrm_a %>% mutate(dataset="CNRM"),hadgem_a %>% mutate(dataset="HadGEM"),ec_a %>% mutate(dataset="EC-Earth")) else
       bind_rows(m_era5 %>% mutate(dataset="ERA5-Land"),m_cnrm %>% mutate(dataset="CNRM"),m_hadgem %>% mutate(dataset="HadGEM"),m_ec %>% mutate(dataset="EC-Earth"))
  p <- ggplot(ds,aes(if(annual) year else date,temp,color=dataset)) + geom_line(linewidth=if(annual) .6 else .4)
  if(annual) p <- p + geom_point(size=1) + scale_x_continuous(breaks=seq(1985,2015,5))
  p + scale_color_manual(values=setNames(colours[1:4],c("ERA5-Land","CNRM","HadGEM","EC-Earth"))) +
    labs(x=if(annual) "Year" else "Date", y=if(annual) "Mean Annual Temperature (°C)" else "Temperature (°C)", color="Dataset") + theme_minimal()
}
ggsave(file.path(fig_dir,"Figure_2_annual_timeseries.png"),base_plot(NULL,TRUE),width=9,height=6,dpi=300)
ggsave(file.path(fig_dir,"Figure_3_monthly_timeseries.png"),base_plot(NULL,FALSE),width=10,height=6,dpi=300)

# ---- 6. Seasonal metrics and Figure 5 --------------------------------------
seasonal_fun <- function(model,name) {
  a <- spatial_monthly(model) %>% bind_cols(get_season_year(spatial_monthly(model)$date)) %>% filter(season_year %in% ref_years)
  o <- spatial_monthly(era5_land_9km_mon) %>% bind_cols(get_season_year(spatial_monthly(era5_land_9km_mon)$date)) %>% filter(season_year %in% ref_years)
  a <- a %>% group_by(season,season_year) %>% summarise(temp_model=mean(temp,na.rm=TRUE),n=sum(!is.na(temp)),.groups="drop") %>% filter(n==3)
  o <- o %>% group_by(season,season_year) %>% summarise(temp_obs=mean(temp,na.rm=TRUE),n=sum(!is.na(temp)),.groups="drop") %>% filter(n==3)
  inner_join(o,a,by=c("season","season_year")) %>% group_by(season) %>% summarise(RMSE=sqrt(mean((temp_obs-temp_model)^2)),MBE=mean(temp_obs-temp_model),Correlation=cor(temp_model,temp_obs),n_years=n(),.groups="drop") %>% mutate(Model=name)
}
seasonal_metrics <- bind_rows(seasonal_fun(cnrm_9km_ref_mon,"CNRM"),seasonal_fun(hadgem_9km_ref_mon,"HadGEM"),seasonal_fun(ec_earth_9km_ref_mon,"EC-Earth"))
write.csv(seasonal_metrics,file.path(tab_dir,"Table_4_seasonal_metrics.csv"),row.names=FALSE)
Figure_5 <- ggplot(seasonal_metrics,aes(factor(Model,levels=c("CNRM","HadGEM","EC-Earth")),factor(season,levels=c("DJF","MAM","JJA","SON")),fill=MBE))+geom_tile(color="white")+geom_text(aes(label=sprintf("%.2f",MBE)),size=3.5)+scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,name="MBE (°C)")+labs(x="",y="Season",title="Seasonal Temperature Bias by GCM (Cork County, 1985-2014)")+theme_minimal()
ggsave(file.path(fig_dir,"Figure_5_seasonal_bias_heatmap.png"),Figure_5,width=8,height=5.5,dpi=300)

# ---- 7. Split-period stability and GCM ranking ------------------------------
period_fun <- function(model,name,a,b) {
  x <- spatial_monthly(model) %>% filter(year>=a,year<=b); o <- spatial_monthly(era5_land_9km_mon) %>% filter(year>=a,year<=b)
  j <- inner_join(o,x,by="date",suffix=c("_era5","_model"))
  tibble(Model=name,Period=paste0(a,"-",b),RMSE=sqrt(mean((j$temp_era5-j$temp_model)^2)),MBE=mean(j$temp_era5-j$temp_model),Correlation=cor(j$temp_model,j$temp_era5))
}
stability <- bind_rows(period_fun(cnrm_9km_ref_mon,"CNRM",1985,1999),period_fun(cnrm_9km_ref_mon,"CNRM",2000,2014),period_fun(hadgem_9km_ref_mon,"HadGEM",1985,1999),period_fun(hadgem_9km_ref_mon,"HadGEM",2000,2014),period_fun(ec_earth_9km_ref_mon,"EC-Earth",1985,1999),period_fun(ec_earth_9km_ref_mon,"EC-Earth",2000,2014))
write.csv(stability,file.path(tab_dir,"Table_5_split_period_metrics.csv"),row.names=FALSE)
write.csv(stability %>% select(Model,Period,MBE) %>% pivot_wider(names_from=Period,values_from=MBE) %>% mutate(Delta_MBE=`2000-2014`-`1985-1999`),file.path(tab_dir,"split_period_MBE_change.csv"),row.names=FALSE)
ranking <- annual_metrics %>% mutate(RMSE_s=1-rescale(RMSE),MAE_s=1-rescale(MAE),MBE_s=1-rescale(abs(MBE)),COR_s=rescale(Correlation),EXT_s=1-rescale(Extreme_bias),Final_score=100*rowMeans(cbind(RMSE_s,MAE_s,MBE_s,COR_s,EXT_s),na.rm=TRUE)) %>% arrange(desc(Final_score))
write.csv(ranking,file.path(tab_dir,"GCM_composite_ranking.csv"),row.names=FALSE)
selected_gcm <- "CNRM"

# ---- 8. Taylor diagram ------------------------------------------------------
taylor <- function(x) spatial_monthly(x) %>% filter(year %in% ref_years) %>% select(date,temp)
tdf <- taylor(era5_land_9km_mon)%>%rename(ERA5=temp)%>%left_join(taylor(cnrm_9km_ref_mon)%>%rename(CNRM=temp),by="date")%>%left_join(taylor(hadgem_9km_ref_mon)%>%rename(HadGEM=temp),by="date")%>%left_join(taylor(ec_earth_9km_ref_mon)%>%rename(`EC-Earth`=temp),by="date")%>%drop_na()
tlong <- pivot_longer(tdf,cols=c(CNRM,HadGEM,`EC-Earth`),names_to="Model",values_to="Model_temp")
png(file.path(fig_dir,"Figure_4_taylor_diagram.png"),width=2400,height=2000,res=300); TaylorDiagram(tlong,obs="ERA5",mod="Model_temp",group="Model",normalise=FALSE,main="Taylor Diagram: GCMs vs ERA5 (1985-2014)"); dev.off()

# ---- 9. GLM bias correction -------------------------------------------------
# Training is 1985-2014. spatial.predictors = 0.95, Gaussian GLM.
downscale_data <- prepareData(x=cnrm_9km_ref_mon,y=era5_land_9km_mon,spatial.predictors=list(v.exp=0.95))
regression <- downscaleTrain(downscale_data,method="GLM",family=gaussian)
cnrm45_9km$xyCoords <- attr(downscale_data,"xyCoords")
cnrm85_9km$xyCoords <- attr(downscale_data,"xyCoords")
p45 <- prepareNewData(newdata=cnrm45_9km,data.structure=downscale_data)
p85 <- prepareNewData(newdata=cnrm85_9km,data.structure=downscale_data)
regression_predictions45 <- downscalePredict(newdata=p45,model=regression,simulate=FALSE)
regression_predictions85 <- downscalePredict(newdata=p85,model=regression,simulate=FALSE)

# ---- 10. Future domain validation and Table 6 --------------------------------
future_spatial <- function(x) spatial_monthly(x) %>% filter(year %in% future_years)
E <- future_spatial(era5_land_9km_mon); R45 <- future_spatial(cnrm45_9km); B45 <- future_spatial(regression_predictions45); R85 <- future_spatial(cnrm85_9km); B85 <- future_spatial(regression_predictions85)
plus <- function(obs,sim) { tibble(RMSE=sqrt(mean((obs-sim)^2)),MBE=mean(obs-sim),MAE=mean(abs(obs-sim)),Correlation=cor(obs,sim)) }
tail <- function(obs,sim) { q05=quantile(obs,.05); q95=quantile(obs,.95); tibble(Cold=mean((obs-sim)[obs<=q05]),Warm=mean((obs-sim)[obs>=q95])) }
row_future <- function(obs,sim,scenario,label) bind_cols(plus(obs$temp,sim$temp),tail(obs$temp,sim$temp)) %>% mutate(Scenario=scenario,Model=label)
table_6 <- bind_rows(row_future(E,R45,"SSP2-4.5","CNRM raw"),row_future(E,B45,"SSP2-4.5","CNRM bias-corrected"),row_future(E,R85,"SSP5-8.5","CNRM raw"),row_future(E,B85,"SSP5-8.5","CNRM bias-corrected")) %>% select(Scenario,Model,RMSE,MBE,MAE,Correlation,Cold,Warm)
write.csv(table_6,file.path(tab_dir,"Table_6_bias_correction_metrics.csv"),row.names=FALSE)

# ---- 11. Figures 7 and 8 ----------------------------------------------------
err <- bind_rows(R45%>%select(date,CNRM=temp)%>%left_join(E%>%select(date,ERA5=temp),by="date")%>%mutate(Scenario="SSP2-4.5",Type="Raw"),B45%>%select(date,CNRM=temp)%>%left_join(E%>%select(date,ERA5=temp),by="date")%>%mutate(Scenario="SSP2-4.5",Type="Bias corrected"),R85%>%select(date,CNRM=temp)%>%left_join(E%>%select(date,ERA5=temp),by="date")%>%mutate(Scenario="SSP5-8.5",Type="Raw"),B85%>%select(date,CNRM=temp)%>%left_join(E%>%select(date,ERA5=temp),by="date")%>%mutate(Scenario="SSP5-8.5",Type="Bias corrected"))%>%mutate(Error=ERA5-CNRM,Dataset=paste(Scenario,Type))
Figure_7 <- ggplot(err,aes(Dataset,Error,fill=Type))+geom_boxplot(width=.65)+geom_hline(yintercept=0,linetype="dashed")+scale_x_discrete(limits=c("SSP2-4.5 Raw","SSP2-4.5 Bias corrected","SSP5-8.5 Raw","SSP5-8.5 Bias corrected"))+labs(x="Scenario and dataset",y="Temperature error (°C)",fill="Dataset type",title="Temperature error distribution: CNRM before and after bias correction")+theme_bw()+theme(legend.position="bottom")
ggsave(file.path(fig_dir,"Figure_7_error_distribution.png"),Figure_7,width=9,height=6,dpi=300)

bias <- B45%>%select(date,year,month,temp_bc=temp)%>%left_join(E%>%select(date,temp_era5=temp),by="date")%>%mutate(bias=temp_era5-temp_bc,month_label=factor(month.abb[month],levels=month.abb))
Figure_8 <- ggplot(bias,aes(factor(year),month_label,fill=bias))+geom_tile(color="white")+scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,limits=c(-4,4),name="Bias (°C)")+labs(title="Monthly Bias: CNRM bias-corrected vs ERA5 Mean (2015-2025)",x="Year",y="Month")+theme_minimal()
ggsave(file.path(fig_dir,"Figure_8_monthly_bias_heatmap.png"),Figure_8,width=9,height=6,dpi=300)

# ---- 12. Station validation, Tables 7 and Figures 9-10 ----------------------
safe_station <- function(obs,sim,label) { j<-inner_join(obs,sim,by="date",suffix=c("_obs","_mod")); if(nrow(j)<12)return(tibble(Model=label,RMSE=NA_real_,MBE=NA_real_,Correlation=NA_real_)); compute_metrics(j$temp_obs,j$temp_mod)%>%select(RMSE,MBE,Correlation)%>%mutate(Model=label) }
station_monthly <- function(st,scenario,base,pred) {
  o<-load_station(st$file)%>%select(date,temp); r<-raw_point(base,st$lon,st$lat)%>%select(date,temp); b<-bc_point(pred,base,st$lon,st$lat)%>%select(date,temp)
  bind_rows(safe_station(o,r,"CNRM raw"),safe_station(o,b,"CNRM bias-corrected"))%>%mutate(Station=st$station,Scenario=scenario)
}
station_annual <- function(st,scenario,base,pred) {
  o<-load_station(st$file)%>%mutate(dataset="Station")%>%annualise()%>%select(year,temp); r<-raw_point(base,st$lon,st$lat)%>%annualise()%>%select(year,temp); b<-bc_point(pred,base,st$lon,st$lat)%>%annualise()%>%select(year,temp)
  jr<-inner_join(o,r,by="year",suffix=c("_obs","_mod")); jb<-inner_join(o,b,by="year",suffix=c("_obs","_mod"))
  bind_rows(compute_metrics(jr$temp_obs,jr$temp_mod),compute_metrics(jb$temp_obs,jb$temp_mod))%>%select(RMSE,MBE,Correlation)%>%mutate(Model=c("CNRM raw","CNRM bias-corrected"),Station=st$station,Scenario=scenario)
}
monthly_station_metrics <- bind_rows(lapply(seq_len(nrow(stations)),function(i)bind_rows(station_monthly(stations[i,],"SSP2-4.5",cnrm45_9km,regression_predictions45),station_monthly(stations[i,],"SSP5-8.5",cnrm85_9km,regression_predictions85))) )
annual_station_metrics <- bind_rows(lapply(seq_len(nrow(stations)),function(i)bind_rows(station_annual(stations[i,],"SSP2-4.5",cnrm45_9km,regression_predictions45),station_annual(stations[i,],"SSP5-8.5",cnrm85_9km,regression_predictions85))) )
write.csv(monthly_station_metrics,file.path(tab_dir,"station_monthly_metrics.csv"),row.names=FALSE)
write.csv(annual_station_metrics,file.path(tab_dir,"station_annual_metrics.csv"),row.names=FALSE)
table_7 <- annual_station_metrics%>%select(Scenario,Station,Model,MBE)%>%pivot_wider(names_from=c(Scenario,Model),values_from=MBE)
write.csv(table_7,file.path(tab_dir,"Table_7_station_MBE.csv"),row.names=FALSE)

fig9data <- annual_station_metrics%>%mutate(Type=ifelse(Model=="CNRM raw","Raw","Bias-Corrected"),Station=case_when(Station=="CORK AIRPORT"~"Cork Airport",Station=="MOORE PARK"~"Moore Park",Station=="ROCHES POINT"~"Roches Point",TRUE~"Sherkin Island"))
Figure_9 <- ggplot(fig9data,aes(Station,RMSE,fill=Type))+geom_col(position="dodge")+facet_wrap(~Scenario,ncol=1)+scale_fill_manual(values=c("Raw"="orange","Bias-Corrected"="darkgreen"))+labs(x="Station",y="RMSE (°C)",fill="Model Type",title="Future Performance: Raw vs Bias-Corrected CNRM")+theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="bottom")
ggsave(file.path(fig_dir,"Figure_9_station_RMSE.png"),Figure_9,width=8,height=7,dpi=300)

station_ann_plot <- bind_rows(lapply(seq_len(nrow(stations)),function(i){st<-stations[i,]; bind_rows(load_station(st$file)%>%mutate(dataset="Station",Station=st$station)%>%annualise(),extract_point(era5_land_9km_mon,st$lon,st$lat)%>%filter(year%in%future_years)%>%mutate(dataset="ERA5",Station=st$station)%>%annualise(),raw_point(cnrm45_9km,st$lon,st$lat)%>%annualise()%>%mutate(Scenario="SSP2-4.5",Station=st$station),bc_point(regression_predictions45,cnrm45_9km,st$lon,st$lat)%>%annualise()%>%mutate(Scenario="SSP2-4.5",Station=st$station),raw_point(cnrm85_9km,st$lon,st$lat)%>%annualise()%>%mutate(Scenario="SSP5-8.5",Station=st$station),bc_point(regression_predictions85,cnrm85_9km,st$lon,st$lat)%>%annualise()%>%mutate(Scenario="SSP5-8.5",Station=st$station))}))
station_ann_plot <- station_ann_plot%>%mutate(Scenario=ifelse(is.na(Scenario),"Observed/reference",Scenario))
Figure_10 <- ggplot(station_ann_plot,aes(year,temp,color=dataset))+geom_line(linewidth=.8)+facet_grid(Station~Scenario,scales="free_y")+labs(x="Year",y="Temperature (°C)",color="Dataset",title="Annual Mean Temperature by Station and Scenario")+theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust=.5),legend.position="bottom")
ggsave(file.path(fig_dir,"Figure_10_station_annual_comparison.png"),Figure_10,width=11,height=9,dpi=300)

# ---- 13. Completion ---------------------------------------------------------
message("Analysis completed. Figures: ", normalizePath(fig_dir), " | Tables: ", normalizePath(tab_dir))
