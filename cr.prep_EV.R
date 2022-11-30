
##Load Spatial packages
library(here)
library(janitor)
library(sf)
library(lubridate)
library(openCR) 
library(hrbrthemes)
library(units)
library(tidyverse)
library(tidylog)


####LOOK AT 2018 DETECTIONS HISTORIES

##################################################################
####LOAD CAPTURE DATA
##################################################################

df <- read_csv(here::here("data","caprecap","SRGBP_SPI_Export_6Nov2022.csv"))%>%
  rbind(read_csv(here::here("data","caprecap","FH2007_SPI_Export_13May2022.csv")))%>%
  clean_names()

##Add in "Session Column" for SECR
df <- df%>%
  mutate(year=year(visit_end_date),
         occassion=str_sub(session,-1,-1)%>%as.numeric)%>%
  mutate(occassion_text=paste0("v",occassion))

##clean up
df<-df%>%
  drop_na(longitude)%>%##Cull any missing spatial data
  filter(!sample_is_incidental%in%"Y", ##remove incidental records
         occassion>0) ##session 0 is rub tree first check, has hair from previous year
  

##traps
df <- df%>%
  mutate(trap.type=case_when(trap_type_code%in%c("HR-BA-BS")~"bs",
                             trap_type_code%in%c("HR-BA-RS","HR-NB-RS","HR-BA-TR","HR-DC")~"rt")) ## pool trail and fence with rt's, similar lambda0 estimates in early runs

  # mutate(trap.type=case_when(type_of_device_code%in%c("HR-BA-BS")~"bs",
  #                            type_of_device_code%in%c("HR-BA-RS","HR-NB-RS")~"rt",
  #                            type_of_device_code%in%c("HR-BA-TR","HR-DC")~"other"))


df <- df%>%filter(site!="s0810") ##in a town, wrong location


sessions <- df%>%distinct(year)%>%arrange(year)%>%pull(year)

##summary stats
##how many individual bears?
df%>%
  filter(species%in%"M-URAR")%>%
  distinct(animal_id)%>%
  nrow()

##how many detections
df%>%
  filter(species%in%"M-URAR")%>%
  distinct(animal_id,site,year, occassion)%>%
  nrow()


a <- df%>%
  group_by(trap.type,year)%>%
  summarize(n_distinct(occassion))
##clip to EV only

##load EV study area
sa <-st_read(here::here("data","studyarea","EV_grizz_sa.shp"))%>%
  st_transform(3005)%>%
  mutate(area="EV")%>%
  dplyr::select(area)

##average hr size of animals is ~220 sq.km, so buffer by -8.5km radius
sa.small <- st_buffer(sa, (-5*1000))
st_area(sa.small)/1E6
#library(mapview)
#mapview(sa)+mapview(sa.small)

df <-df%>%
  st_as_sf(coords=c("longitude","latitude"), crs=4326)%>%
  st_transform(3005)

df <- df%>%st_intersection(sa.small)
#mapview(sa.small)+mapview(df)

##back to tibble
df <- df%>%
  cbind(st_coordinates(.))%>%
  tibble()


##summary stats

##how many individual bears?
df%>%
  filter(species%in%"M-URAR")%>%
  distinct(animal_id)%>%
  nrow()

##how many detections
df%>%
  filter(species%in%"M-URAR")%>%
  distinct(animal_id,site,year, occassion)%>%
  nrow()

##detections across occasions
df%>%
  filter(species%in%"M-URAR")%>%
  distinct(animal_id,site,year, occassion,trap.type)%>%
  group_by(trap.type,occassion)%>%
  summarize(n=n())

########################################
### put into SECR format
########################################

  
##Retain columns needed for Capture Data,
  CapData<-df%>%
  filter(species%in%"M-URAR")%>%
    select(`#session`=year,
           id=animal_id,
           occassion,
           detector=site,
           sex)
  
  ##Remove all non-capture data
  CapData<-CapData%>%drop_na(id)
  
  ##Remove any duplicated catches at the same detector
  CapData<-CapData%>%
    distinct(.keep_all=TRUE)
  
  ##Order by bear and year
  ##Re-Order dataframe
  CapData<-CapData%>%
    arrange(`#session`,
            id,
            occassion)

  
  ##Export
  for(i in 1:length(sessions)){
    write.table(CapData%>%dplyr::filter(`#session`==sessions[i]),file=here::here("data/caprecap/clean/cap",paste0("cap_",sessions[i],".txt")), sep=",", row.names=F, quote=FALSE, col.names=TRUE)
  }
  

  #write.table(CapData,file=here::here("/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/SR_Demography/data/caprecap/clean/", paste(df$Project[1], ".txt",sep="")), sep=",", row.names=F, quote=FALSE, col.names=TRUE)
  
  
  ########################################
  ########################################
  ## PREP trap data
  ########################################
  ########################################
  #st_transform(3005)
  
  ###### Prepare Detector Layout  ##########
  
  ##Retain the columns that will be of use for Detector Layout
  dfdetectcull<-df%>%
    select(detector=site,
           X,
           Y,
           trap=trap.type,
           trap_nights=nights_deployed,
           occassion,
           occassion_text,
           session=year)
    
  
  
  ##Prepare Usage Matrix which shows when each site was run
usage <-   dfdetectcull%>% 
    distinct(detector,session,occassion, .keep_all=TRUE)%>%
    group_by(detector,session,occassion)%>%
    summarise(n=sum(trap_nights))%>%
    ungroup%>%
    arrange(detector,session,occassion)%>%
    pivot_wider(names_from=occassion, values_from=n)%>%
    mutate_if(is.numeric , replace_na, replace = 0)%>%
    select(detector,session,`1`,`2`,`3`,`4`,`5`)%>%
    unite("ch",-c("detector","session"), sep=" ")

  ## Add Usage to data
  dfdetectcull<-dfdetectcull%>%
    left_join(usage,by=c("detector","session"))
  
  
  ##convert to BC albers
  # dfdetectcull<-  dfdetectcull%>%
  #   st_as_sf(coords=c("longitude","latitude"), crs=4326)%>%
  #   st_transform(3005)
  
  # library(mapview)
  # mapview(dfdetectcull)
  
  
  ##turn back into tibble
  # dfdetectcull<-  dfdetectcull%>%
  #   cbind(st_coordinates(.))%>%
  #   tibble()
 
   ##Rename Columns
  DetectLayout<-dfdetectcull%>%
    select(`#detector`=detector,
           X,Y,
           usage=ch,
           trap,
           session)

  
  
  ##remove any duplicates
  DetectLayout<-DetectLayout%>%
    distinct()
  
  
  ###throw error if
  if (sum(!CapData$detector %in% DetectLayout$`#detector`)>0){ 
    stop("failed to match some capture locations to detector sites")
  }
  
  ##round coordinates to nearest meter
  DetectLayout$X<- round(DetectLayout$X,0)
  DetectLayout$Y <- round(DetectLayout$Y,0)
  
  DetectLayout$trap <- paste0("/",  DetectLayout$trap)
  
  ##Export 
  for(i in 1:length(sessions)){
  write.table(DetectLayout%>%dplyr::filter(session==sessions[i])%>%dplyr::select(-session),file=here::here("data/caprecap/clean/trap",paste0("trap_",sessions[i],".txt")), sep=",", row.names=F, quote=FALSE, col.names=TRUE)
  }



  
  
  

###PREP DATA
  
  ##load models if running again
  # m1 <- readRDS("mods/m1.rds")
  # m2 <- readRDS("mods/m2.rds")
  # m3 <- readRDS("mods/m3.rds")
  # m4 <- readRDS("mods/m4.rds")

##list to populate
  ch.list <- list()
  
  for(i in 1:length(sessions)){
    ##Read capture history files
    a <- read.capthist(here::here("data/caprecap/clean/cap",paste0("cap_",sessions[i],".txt")), 
                       here::here("data/caprecap/clean/trap",paste0("trap_",sessions[i],".txt")),
                       detector ="proximity",
                       trapcovnames=c("trap"),
                       covnames=c("sex"),
                       sep = ",",
                       binary.usage=FALSE)
    ch.list[[i]] <- a
    names(ch.list)[i]<- paste0("sess",sessions[i])
  }


##get factors consistent across sessions
traps(ch.list) <- shareFactorLevels(traps(ch.list), columns="trap")

##bind into one CH
grizzCH <- MS.capthist(ch.list)

##fix names
names(grizzCH) <- paste0("sess",sessions)


##have a look
#covariates(traps(grizzCH))
summary(grizzCH)

##make mask
GBmask.o <- make.mask (rbind(traps(grizzCH)), buffer = 30000, type = 'trapbuffer', spacing = 4000)
plot(GBmask.o)


# 
# ##Find a top detection model
# trap.sex <- secr.fit(capthist = grizzCH, 
#                      mask=GBmask.o, 
#                      list(lambda0~trap + g, sigma~g, D~1), 
#                      detectfn="HHN",
#                      trace=TRUE, 
#                      ncores=8,
#                      groups = "sex",
#                      verify=FALSE,
#                      start=readRDS("mods/starts/trap.sex.rds"))
# 
# saveRDS(trap.sex, "mods/trap.sex.rds")
# 
# trap.bk <- secr.fit(capthist = grizzCH, 
#                     mask=GBmask.o, 
#                     list(lambda0~trap*bk + g, sigma~g, D~1), 
#                     detectfn="HHN",
#                     trace=TRUE, 
#                     ncores=8,
#                     groups = "sex",
#                     verify=FALSE,
#                     start=readRDS("mods/starts/trap.bk.rds"))
# 
# saveRDS(trap.bk, "mods/trap.bk.rds")
# 
# trap.pbk <- secr.fit(capthist = grizzCH, 
#                      mask=GBmask.o, 
#                      list(lambda0~trap + bk + g, sigma~g, D~1), 
#                      detectfn="HHN",
#                      trace=TRUE, 
#                      ncores=8,
#                      groups = "sex",
#                      verify=FALSE,
#                      start=readRDS("mods/starts/trap.pbk.rds"))
# 
# saveRDS(trap.pbk, "mods/trap.pbk.rds")
# 
# trap.bk.t <- secr.fit(capthist = grizzCH, 
#                     mask=GBmask.o, 
#                     list(lambda0~trap*bk + g + t, sigma~g, D~1), 
#                     detectfn="HHN",
#                     trace=TRUE, 
#                     ncores=8,
#                     groups = "sex",
#                     verify=FALSE,
#                     start=readRDS("mods/starts/trap.bk.t.rds"))
# 
# saveRDS(trap.bk.t, "mods/trap.bk.t.rds")
# 
# AIC(trap.sex, trap.bk,trap.pbk,trap.bk.t)
# collate(trap.sex, trap.bk,trap.pbk,trap.bk.t,
#         realnames='D', perm=c(2,3,4,1))
# 
# ##run trap.bk.t and trap.pbk, basically get same answer, use simple trap.pbk
# 
# 
# #### Detection fn
# trap.pbk.hn <- secr.fit(capthist = grizzCH, 
#                      mask=GBmask.o, 
#                      list(g0~trap + bk + g, sigma~g, D~1), 
#                      detectfn="HN",
#                      trace=TRUE, 
#                      ncores=8,
#                      groups = "sex",
#                      verify=FALSE,
#                      start=readRDS("mods/starts/trap.pbk.rds"))
# 
# saveRDS(trap.pbk, "mods/trap.pbk.hn.rds")
# 
# trap.pbk.hex <- secr.fit(capthist = grizzCH, 
#                         mask=GBmask.o, 
#                         list(lambda0~trap + bk + g, sigma~g, D~1), 
#                         detectfn="HEX",
#                         trace=TRUE, 
#                         ncores=8,
#                         groups = "sex",
#                         verify=FALSE,
#                         start=readRDS("mods/starts/trap.pbk.rds"))
# 
# saveRDS(trap.pbk, "mods/trap.pbk.hex.rds")
# 
# 
# AIC(trap.pbk,trap.pbk.hn,trap.pbk.hex)
# 
# collate(trap.pbk,trap.pbk.hn,trap.pbk.hex,
#         realnames='D', perm=c(2,3,4,1)) ##no effect on D
# 
# collate(trap.pbk,trap.pbk.hn,trap.pbk.hex,
#         realnames='sigma', perm=c(2,3,4,1)) ##large effect on sigma, but HHN most plausable. Stick with that. Murray said caution needed with hex anyways.
# 
# 
# ##sigma and lambda0 changes through time?
# trap.pbk.b <- secr.fit(capthist = grizzCH%>%subset(sessions=11:16), 
#                      mask=GBmask.o, 
#                      list(lambda0~trap + bk + g, sigma~g + Session, D~1), 
#                      detectfn="HHN",
#                      trace=TRUE, 
#                      ncores=8,
#                      groups = "sex",
#                      verify=FALSE)
# 
# 
# trap.pbk.c <- secr.fit(capthist = grizzCH%>%subset(sessions=11:16), 
#                        mask=GBmask.o, 
#                        list(lambda0~trap + bk + g + Session, sigma~g + Session, D~1), 
#                        detectfn="HHN",
#                        trace=TRUE, 
#                        ncores=8,
#                        groups = "sex",
#                        verify=FALSE)
# 
# 
# AIC(trap.pbk,trap.pbk.b,trap.pbk.c)
# 
# collate(trap.pbk,trap.pbk.b,trap.pbk.c,
#         realnames='D', perm=c(2,3,4,1)) ##effect on D



##"lambda0~trap + bk + g, sigma~g + Session" is best detection model






##OPENCR modelling
##run model

# m1.open.starts<- openCR::openCR.fit(capthist = grizzCH, 
#                         mask=GBmask.o, 
#                         type = 'JSSAsecrbCL', 
#                         list(lambda0~trap + bk + sex, sigma~sex, lambda~1), 
#                         detectfn="HHN",
#                         ncores=8)
# 
# saveRDS(m1.open.starts, "mods/starts/m1.open.starts.rds")


m1o<- openCR::openCR.fit(capthist = grizzCH, 
                         mask=GBmask.o, 
                         type = 'JSSAsecrlCL', 
                         list(lambda0~trap + bk + sex, sigma~sex, lambda~1), 
                         detectfn="HHN",
                         ncores=8,
                         start=c(-4.7,-1.7, 1, -0.4,  1.8, 0.02, 8.2,  1.1))
saveRDS(m1o, "mods/m1o.rds")


m1o.EVstudy <-  openCR::openCR.fit(capthist = grizzCH, 
                                       mask=GBmask.o, 
                                       type = 'JSSAsecrlCL', 
                                       list(lambda0~trap + bk + sex, sigma~sex, lambda~scov), 
                                       detectfn="HHN",
                                       sessioncov=c(rep("Early",times=10),rep("Late", times=6)),
                                       ncores=8,
                                       start=c(-4.7,-1.7, 1, -0.4,  1.8,  0.02, -0.04,  8.2,  1.1),
                                       trace=TRUE)
saveRDS(m1o.EVstudy, "mods/m1o.EVstudy.rds")

AIC(m1o,m1o.EVstudy) ##no support for splitting trend, stick with simpler m1o model



##extract estimates
ests <- predict(m1o)$lambda%>%
  tibble%>%
  slice(1)%>%
  select(estimate:ucl)%>%
  mutate(ucl90=estimate-(SE.estimate*1.645),
         lcl90=estimate+(SE.estimate*1.645))


write_csv(ests, here::here("data/caprecap/lambda.opencr.csv"))
##SECR modelling


##Linear trend through time
# m2<- secr.fit(capthist = grizzCH%>%subset(sessions=11:16), 
#               mask=GBmask.o, 
#               list(lambda0~trap + bk + g, sigma~g, D~Session), 
#               detectfn="HHN",
#               trace=TRUE, 
#               ncores=8,
#               groups = "sex",
#               verify=FALSE,
#               start=readRDS("mods/starts/m2.rds"))
# 
# 
# saveRDS(m2, "mods/m2.rds")
# 
# secr.pred.temp <- predict(m2)
# secr.pred <-tibble()
# for(i in 1:length(secr.pred.temp)){
#   a <- tibble(
#     year=str_split(names(secr.pred.temp[i]),",", simplify = TRUE)[,1]%>%str_sub(-4,-1)%>%as.numeric,
#     sex=str_split(names(secr.pred.temp[i]),",", simplify = TRUE)[,2]%>%str_sub(-1,-1),
#     D=secr.pred.temp[i][[1]]$estimate[1]*1E5,
#     D.lcl=secr.pred.temp[i][[1]]$lcl[1]*1E5,
#     D.ucl=secr.pred.temp[i][[1]]$ucl[1]*1E5,
#     N=secr.pred.temp[i][[1]]$estimate[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
#     N.lcl=secr.pred.temp[i][[1]]$lcl[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
#     N.ucl=secr.pred.temp[i][[1]]$ucl[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric)
#   
#   secr.pred <- bind_rows(secr.pred,a)
# }
# secr.pred <- secr.pred%>%
#   distinct(year,sex, .keep_all=TRUE)%>%
#   group_by(year)%>%
#   dplyr::summarise_if(is.numeric,sum)
# 
# 
# ggplot(secr.pred, aes(x=year, y=N, ymin=N.lcl, ymax=N.ucl))+
#   geom_ribbon(alpha=0.2)+
#   geom_line()+
#   theme_ipsum()+
#   labs(x="Year", y="N", color="Sex",
#        title="B) Abundance", subtitle="1462 genetic capture-recaptures of 291 individuals (2006-2021)")+
#   expand_limits(y=0)+
#   theme_ipsum()+
#   theme(axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         strip.text.x = element_text(size=15),
#         strip.text.y = element_text(size=15),
#         legend.text = element_text(size=13),
#         legend.title=element_text(size=15),
#         legend.position = "bottom")
# 
# 
# write_csv(secr.pred, here::here("data/caprecap/dens.pred.csv"))


##Individual D each year 
m3<- secr.fit(capthist = grizzCH%>%subset(sessions=11:16), 
              mask=GBmask.o, 
              list(lambda0~trap + bk + g, sigma~g+Session, D~session), 
              detectfn="HHN",
              trace=TRUE, 
              ncores=8,
              groups = "sex",
              verify=FALSE)

saveRDS(m3, "mods/m3.rds")
m3 <- readRDS(here::here("mods/m3.rds"))

secr.pred.temp <- predict(m3)
secr.pred <-tibble()
for(i in 1:length(secr.pred.temp)){
  a <- tibble(
    year=str_split(names(secr.pred.temp[i]),",", simplify = TRUE)[,1]%>%str_sub(-4,-1)%>%as.numeric,
    sex=str_split(names(secr.pred.temp[i]),",", simplify = TRUE)[,2]%>%str_sub(-1,-1),
    D=secr.pred.temp[i][[1]]$estimate[1]*1E5,
    D.se=secr.pred.temp[i][[1]]$SE.estimate[1]*1E5,
    D.lcl=secr.pred.temp[i][[1]]$lcl[1]*1E5,
    D.ucl=secr.pred.temp[i][[1]]$ucl[1]*1E5,
    N=secr.pred.temp[i][[1]]$estimate[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.se=secr.pred.temp[i][[1]]$SE.estimate[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.lcl=secr.pred.temp[i][[1]]$lcl[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.ucl=secr.pred.temp[i][[1]]$ucl[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric)
  
  
  secr.pred <- bind_rows(secr.pred,a)
}

secr.pred <- secr.pred%>%
  distinct(year,sex, .keep_all=TRUE)%>%
  group_by(year)%>%
  dplyr::summarise_if(is.numeric,sum)


ggplot(secr.pred, aes(x=year, y=N, ymin=N.lcl, ymax=N.ucl))+
  geom_ribbon(alpha=0.2)+
  geom_line()+
  theme_ipsum()+
  labs(x="Year", y="N", color="Sex",
       title="B) Abundance", subtitle="1462 genetic capture-recaptures of 291 individuals (2006-2021)")+
  expand_limits(y=0)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "bottom")

write_csv(secr.pred, here::here("data/caprecap/dens.pred.annual.csv"))



##pooled D for across years during study
m4<- secr.fit(capthist = grizzCH%>%subset(sessions=11:16), 
              mask=GBmask.o, 
              list(lambda0~trap + bk + g, sigma~g+Session, D~1), 
              detectfn="HHN",
              trace=TRUE, 
              ncores=8,
              groups = "sex",
              verify=FALSE,
              start=readRDS("mods/starts/m4.rds"))

saveRDS(m4, "mods/m4.rds")


secr.pred.temp <- predict(m4)
secr.pred <-tibble()
for(i in 1:length(secr.pred.temp)){
  a <- tibble(
    year=str_split(names(secr.pred.temp[i]),",", simplify = TRUE)[,1]%>%str_sub(-4,-1)%>%as.numeric,
    sex=str_split(names(secr.pred.temp[i]),",", simplify = TRUE)[,2]%>%str_sub(-1,-1),
    D=secr.pred.temp[i][[1]]$estimate[1]*1E5,
    D.se=secr.pred.temp[i][[1]]$SE.estimate[1]*1E5,
    D.lcl=secr.pred.temp[i][[1]]$lcl[1]*1E5,
    D.ucl=secr.pred.temp[i][[1]]$ucl[1]*1E5,
    D.lcl90=((secr.pred.temp[i][[1]]$estimate[1])-(secr.pred.temp[i][[1]]$SE.estimate[1]*1.645))*1E5,
    D.ucl90=((secr.pred.temp[i][[1]]$estimate[1])+(secr.pred.temp[i][[1]]$SE.estimate[1]*1.645))*1E5,
    
    N=secr.pred.temp[i][[1]]$estimate[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.se=secr.pred.temp[i][[1]]$SE.estimate[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.lcl=secr.pred.temp[i][[1]]$lcl[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.ucl=secr.pred.temp[i][[1]]$ucl[1]*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.lcl90=(((secr.pred.temp[i][[1]]$estimate[1])-(secr.pred.temp[i][[1]]$SE.estimate[1]*1.645)))*st_area(sa.small)%>%set_units("ha")%>%as.numeric,
    N.ucl90=((secr.pred.temp[i][[1]]$estimate[1])+(secr.pred.temp[i][[1]]$SE.estimate[1]*1.645))*st_area(sa.small)%>%set_units("ha")%>%as.numeric)
  
  secr.pred <- bind_rows(secr.pred,a)
}

secr.pred <- secr.pred%>%
  distinct(year,sex, .keep_all=TRUE)%>%
  group_by(year)%>%
  dplyr::summarise_if(is.numeric,sum)

write_csv(secr.pred, here::here("data/caprecap/dens.pred.combined.csv"))




###figure out low last year
dens <-tibble()
for(i in 1:16){
a<- secr.fit(capthist = grizzCH%>%subset(sessions=i),
                 mask=GBmask.o,
                 list(lambda0~trap, sigma~1, D~1),
                 detectfn="HHN",
                 trace=TRUE,
                 ncores=8,
                 verify=FALSE,
             start=c(-8,-6,0.08,8.3))

secr.pred.temp <- predict(a)

dens <- dens%>%
  rbind(
    tibble(
  year=i,
  param=c("D","lambda0","sigma"),
  secr.pred.temp))
}


dens%>%
  mutate(year=rep(2006:2021, each=3))%>%
  filter(year>2006)%>%
  ggplot(aes(x=year,y=estimate,ymin=lcl,ymax=ucl))+
  geom_point()+
  geom_line()+
  geom_linerange()+
  facet_wrap(vars(param), scales="free_y")+
  expand_limits(y=0)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "bottom")

