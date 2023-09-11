Elk Valley Grizzly Demography 2016-2022
================
Clayton T. Lamb
05 September, 2023

## Load Data

``` r
library(here)
library(ggnewscale)
library(tidyterra)
library(survival)
library(hrbrthemes)
library(lubridate)
library(survminer)
library(forcats)
library(ggdist)
library(patchwork)
library(ggmap)
library(RStoolbox)
library(sf)
library(raster)
library(terra)
library(basemaps)
library(amt)
library(RColorBrewer)
library(demogR)
library(mapview)
library(lwgeom)
library(units)
library(ggallin)
library(readxl)
library(scales)
library(knitr)
library(tidyverse)
library(tidylog)

##set params
study.end <- ymd("2022-11-06") ##median den entry date from Lamb et al. 2022 WCB
lamba.opencr <- read_csv(here::here("data/caprecap/lambda.opencr.csv"))$estimate

## Load Data and trim to <study end
telem <- st_read(here::here("data", "telem","shp","EVcollar_Relocs_raw.shp"))%>%
  filter(!Name%in%c("Sheena","Honey"),
         DateTime<study.end)%>%
  st_transform(3005)

cap.raw <- read_csv(here::here("data", "Grizzly Capture Form","Grizzly Capture Form.csv"))%>%
  filter(!`Capture ID`%in%c("Sheena","Honey"),
         !`Collar Number`%in%21015,# Sid's Alberta capture
         `Date and Time`<study.end)

mort.raw <- read_csv(here::here("data", "Grizzly Mortality Form.csv"))%>%
  filter(!`Capture ID`%in%c("Sheena","Honey"),
         `Suspected Date of Mortality`<study.end)  

cubs.raw <- read_csv(here::here("data", "Cubs","Cubs.csv"))%>%
  left_join(
    read_csv(here::here("data", "Cubs","LinkTo-Cubs Observations.csv")),
    by=c("form_record_id"="parent_record_id")
    )%>%
  filter(!`Capture ID`%in%c("Sheena","Honey"),
         Date<study.end)


conflict.raw <- read_csv(here::here("data","conflict","encounter-356181-176470.csv"))%>%
  filter(encounter_date<study.end)

#####Load some other spatial data#####
##Pacific northwest for inset
pnw<- st_read(here::here("data/administrative/North_America.shp"))%>%
  filter(FID_usa%in%c(1,2,6,8,11)|FID_canada%in%c(8,9))%>%
  mutate(prov=c("WA","MT","WY","ID","OR","AB","BC","BC"),
         country=c(rep("USA",times=5),rep("CAN",times=3)))%>%
  st_make_valid()%>%
  group_by(prov, country)%>%
  summarise(id=mean(FID_canada))%>%
  ungroup%>%
  st_transform(3005)


##grizz range
range <- st_read(here::here("data/range/Current_LambEdits_poly.shp"))%>%
  st_transform(3005)%>%
  st_make_valid()%>%
  st_intersection(st_make_valid(pnw))
```

## Delineate study area

``` r
##calculate UD for each individual 

telem.amt <-  telem%>%
  # left_join(cap.raw%>%select(Name=`Capture ID`,Sex)%>%distinct())%>%
  # filter(Sex=="F")%>%
  filter(DateTime>ymd("2016-08-14"))%>%
  cbind(st_coordinates(.))%>%
  tibble()%>%
  mutate(DateTime=ymd_hms(DateTime),
         Name2=paste0(Name, year(DateTime)))%>%
  group_by(Name)%>%
  #group_by(Name2)%>%
  # mutate(dur=(max(DateTime)-min(DateTime))%>%as.numeric("days"))%>%
  # filter(dur>120)%>%
  ungroup%>%
  make_track(X, Y,DateTime, id = Name, 
             crs=sp::CRS("+init=epsg:3005"))


hr <-telem.amt%>%
  nest(data =c(x_, y_, t_))%>%  
  mutate(kde = map(data, ~ hr_isopleths(hr_kde(.,level = c(0.95)))))%>%
  select(id,kde)%>%
  unnest(kde)

hr <- hr%>%
  pull(geometry)%>%
  st_as_sf%>%
  mutate(id=hr$id,
         area=st_area(.)%>%as.numeric/1E6)

hr.center <- hr%>%st_centroid()

##export hr centers for designating DNA area
st_write(hr.center, here::here("data","studyarea","hr_centers.shp"), delete_dsn = TRUE)
```

    ## Deleting source `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/studyarea/hr_centers.shp' using driver `ESRI Shapefile'
    ## Writing layer `hr_centers' to data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/studyarea/hr_centers.shp' using driver `ESRI Shapefile'
    ## Writing 64 features with 2 fields and geometry type Point.

``` r
##do again for all animals pooled
sa <-  telem%>%
  filter(DateTime>ymd("2016-08-14"))%>%
  cbind(st_coordinates(.))%>%
  tibble()%>%
  mutate(DateTime=ymd_hms(DateTime),
         ID="A")%>% ##make a single individual, don't do individual HR's
  amt::make_track(X, Y,DateTime, id = ID, 
                  crs=3005)%>%
  amt::hr_kde(level = c(0.99))%>%
  hr_isopleths()

st_write(sa, here::here("data","studyarea","EV_grizz_sa.shp"), delete_dsn = TRUE)
```

    ## Deleting source `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/studyarea/EV_grizz_sa.shp' using driver `ESRI Shapefile'
    ## Writing layer `EV_grizz_sa' to data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/studyarea/EV_grizz_sa.shp' using driver `ESRI Shapefile'
    ## Writing 1 features with 3 fields and geometry type Multi Polygon.

``` r
##map
ggplot()+
  geom_sf(data=sa, fill=NA)+
  geom_sf(data=sa%>%st_buffer(-5000), fill="grey")+
  geom_sf(data=hr.center)+
  theme_classic()+
  labs(title="Study Area Extents", subtitle="Larger area is primary study area defined by 99th percentile of telemetry\nShaded is interior 5 km buffer for DNA capture-recapture encompassing homerange centers")
```

![](README_files/figure-gfm/Delineate%20study%20area-1.png)<!-- -->

``` r
ggsave(here::here("plots","study_extents.png"), width=8, height=10, bg="white")
```

## Capture data

``` r
## Clean up capture data
cap <- cap.raw%>%
  filter(Species%in%"Grizzly Bear")%>%
  select(BearID=`Bear #`,
         CapID=`Capture ID`,
         ColBrand=`Collar Brand`,
         ColNum=`Collar Number`,
         Sex,
         Age_est=`Age Estimate`,
         Age=`Tooth Age`,
         Weight=`Weight (kg)`,
         CapReason=`Capture Reason`,
         Date=`Date and Time`,
         lat=`Location (latitude)`,
         long=`Location (longitude)`,
         CapType=`Capture Type`,
         Neck=`Neck Size`,
         svl.nat=`Snout to Vent Length (natural contours)`,
         svl.strt=`Snout to Vent Length (straight line)`,
         zoo.l=`Zoo Length`,
         Chest=`Chest Girth`,
         Zyg.Width=`Zygomatic Width`,
         Zyg.Leng=`Zygomatic Length`,
         Fat=`Percent Body Fat (grizzly only)`,
         GPSend=`GPS Transmit End Date`,
         GPSendReas=`GPS Transmit End Reason`,
         VHFend=`VHF Transmit End Date`,
         VHFendReas=`VHF Transmit End Reason`,
         Detectafter=`Other Detection Date`,
         DectectafterType=`Other Detection Reason`,
         `Fate unknown?`
         )

## use age guess if no known age,
##pull out year and seasons
cap <- cap%>%
  mutate(age.combine=case_when(is.na(Age)~Age_est,
                               TRUE~Age),
         year=year(Date),
         month=month(Date),
         Season=case_when(month(Date,label=TRUE) %in% c("Apr", "May", "Jun")~ "Spring",
                          month(Date,label=TRUE) %in% c("Jul", "Aug")~ "Summer",
                          month(Date,label=TRUE) %in% c("Sep", "Oct", "Nov")~ "Fall"))%>%
  mutate(AgeClass=case_when(age.combine%in%0:1~"Dependent (0-1)",
                            age.combine%in%2:6~"Subadult (2-6)",
                            age.combine>6~"Adult (>6)"))

##Capture Data Plots

##Weight
set.seed(6789);ggplot(cap, aes(x=age.combine, y=Weight, color=Sex%>%fct_relevel("M","F")))+
  stat_smooth(geom='line', alpha=0.8, se=FALSE, linetype="dashed",size=1, method="loess")+
  geom_jitter(alpha=0.9,width=0.5)+
  scale_color_manual(values=RColorBrewer::brewer.pal(5, "Accent")[c(1,5)])+
  labs(x="Age", y="Weight (Kg)", color="Sex",
       title="c) Weight by Age")+
  expand_limits(y=0)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = "bottom")->weight.plot
weight.plot
```

![](README_files/figure-gfm/Capture%20data-1.png)<!-- -->

``` r
##Fat
fat.plot <- ggplot(cap%>%
         arrange(desc(CapReason))%>%
         mutate(date2000=as.Date(paste0("2000-",format(Date%>%ymd_hm, "%j")), "%Y-%j"))
         )+
  stat_smooth( aes(x=date2000, y=Fat), color=brewer.pal(3, "Accent")[c(3)], geom='line', alpha=1, se=FALSE, linetype="dashed",size=1, method="loess", span=1)+
  geom_point( aes(x=date2000, y=Fat, color=CapReason, shape=Sex), alpha=0.95,width=0.15, size=2)+
  scale_color_manual(values=RColorBrewer::brewer.pal(3, "Accent")[c(2:3)])+
  labs(x="Month", y="Fat (%)", color="Capture reason",shape="Sex",
       title="d) % Fat by Month")+
  expand_limits(y=0)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = "bottom",
        legend.box="vertical", legend.margin=margin())+
  scale_x_date(labels = date_format("%b"), breaks = date_breaks("1 month"))+
    guides(color = guide_legend(nrow = 2))
fat.plot
```

![](README_files/figure-gfm/Capture%20data-2.png)<!-- --> \## Summarize
capture data

``` r
cap%>%
  summarise(captured_inds=n_distinct(BearID),
            total_caps=n(),
            management=sum(CapReason=="Management")
            )

cap%>%
  drop_na(ColNum)%>%
  summarise(captured_inds=n_distinct(BearID),
            total_caps=n(),
            management=sum(CapReason=="Management")
  )


cap%>%
  count(CapType)

cap%>%
  drop_na(ColNum)%>%
  summarise(collared_inds=n_distinct(BearID),
            total_collars=n_distinct(ColNum)
  )

cap%>%
  #drop_na(ColNum)%>%
  group_by(Sex,AgeClass)%>%
  summarise(weight=mean(Weight, na.rm=TRUE),
            age=mean(age.combine,na.rm=TRUE),
            age.max=max(age.combine,na.rm=TRUE),
            fat=mean(Fat,na.rm=TRUE),
            fat.max=max(Fat,na.rm=TRUE),
            n=n())

cap%>%
  group_by(Sex)%>%
  summarise(weight=mean(Weight, na.rm=TRUE),
            age=mean(age.combine,na.rm=TRUE),
            age.max=max(age.combine,na.rm=TRUE),
            fat=mean(Fat,na.rm=TRUE),
            fat.max=max(Fat,na.rm=TRUE))

cap%>%
  summarise(weight=mean(Weight, na.rm=TRUE),
            age=mean(age.combine,na.rm=TRUE),
            age.max=max(age.combine,na.rm=TRUE),
            fat=mean(Fat,na.rm=TRUE),
            fat.max=max(Fat,na.rm=TRUE))

cap%>%
  drop_na(ColNum)%>%
  group_by(AgeClass,Sex)%>%
count()

cap <- cap%>%
  mutate(dif=Age_est-Age)

cap%>%
  drop_na(ColNum)%>%
  count(`Fate unknown?`)

cap%>%
  group_by(`Fate unknown?`)%>%
  drop_na(ColNum)%>%
  count(VHFendReas)

cap.summary <- cap%>%
  group_by(Sex, AgeClass)%>%
  summarise(Individuals=n_distinct(BearID),
            Captures=n(),
            Management=sum(CapReason=="Management"),
            Weight.mean=mean(Weight, na.rm=TRUE),
            Weight.min=min(Weight, na.rm=TRUE),
            Weight.max=max(Weight, na.rm=TRUE),
            Age.mean=mean(age.combine,na.rm=TRUE),
            Age.min=min(age.combine,na.rm=TRUE),
            Age.max=max(age.combine,na.rm=TRUE),
            fat.mean=mean(Fat,na.rm=TRUE),
            fat.min=min(Fat,na.rm=TRUE),
            fat.max=max(Fat,na.rm=TRUE),
            neck.mean=mean(Neck,na.rm=TRUE),
            neck.min=min(Neck,na.rm=TRUE),
            neck.max=max(Neck,na.rm=TRUE),
            length.mean=mean(svl.strt,na.rm=TRUE),
            length.min=min(svl.strt,na.rm=TRUE),
            length.max=max(svl.strt,na.rm=TRUE),
            chest.mean=mean(Chest,na.rm=TRUE),
            chest.min=min(Chest,na.rm=TRUE),
            chest.max=max(Chest,na.rm=TRUE)
            )%>%
  mutate_if(is.numeric, round,0)%>%
  ungroup%>%
  mutate(AgeClass=fct_relevel(AgeClass,"Dependent (0-1)","Subadult (2-6)","Adult (>6)"),
         Weight=paste(Weight.mean," (",Weight.min,"-",Weight.max,")", sep=""),
         Age=paste(Age.mean," (",Age.min,"-",Age.max,")", sep=""),
         Fat=case_when(is.na(fat.mean)~"-",
                       TRUE~paste(fat.mean," (",fat.min,"-",fat.max,")", sep="")),
         Neck=paste(neck.mean," (",neck.min,"-",neck.max,")", sep=""),
         Length=paste(length.mean," (",length.min,"-",length.max,")", sep=""),
         Chest=paste(chest.mean," (",chest.min,"-",chest.max,")", sep=""))%>%
  arrange(Sex,Age.mean)%>%
  select(Sex:Management, Age, Weight, Fat, Neck, Chest, Length)

write_csv(cap.summary, here::here("tables/cap.summary.csv"))

kable(cap.summary)
```

| Sex | AgeClass        | Individuals | Captures | Management | Age       | Weight        | Fat        | Neck          | Chest         | Length        |
|:----|:----------------|------------:|---------:|-----------:|:----------|:--------------|:-----------|:--------------|:--------------|:--------------|
| F   | Dependent (0-1) |           4 |        4 |          2 | 0 (0-1)   | 34 (34-34)    | \-         | NaN (Inf–Inf) | NaN (Inf–Inf) | NaN (Inf–Inf) |
| F   | Subadult (2-6)  |          18 |       23 |          4 | 4 (2-6)   | 95 (46-130)   | 24 (5-34)  | 54 (43-61)    | 93 (82-108)   | 142 (118-155) |
| F   | Adult (\>6)     |          20 |       30 |          1 | 11 (7-20) | 125 (84-163)  | 25 (16-39) | 59 (48-68)    | 100 (86-117)  | 148 (130-162) |
| M   | Dependent (0-1) |           4 |        5 |          0 | 1 (0-1)   | 59 (30-83)    | 35 (35-35) | 42 (36-48)    | 72 (60-84)    | 112 (96-128)  |
| M   | Subadult (2-6)  |          17 |       21 |          3 | 4 (2-6)   | 112 (76-156)  | 23 (6-30)  | 56 (48-67)    | 96 (86-111)   | 147 (128-170) |
| M   | Adult (\>6)     |          17 |       27 |          2 | 12 (7-27) | 210 (139-269) | 26 (16-39) | 75 (61-86)    | 121 (104-136) | 172 (154-183) |

# Survival

## Prep survival data

``` r
##Prep survival data
surv.raw <- cap%>%
  drop_na(ColBrand)%>% #keep only collared bears
  mutate(end=pmax(VHFend,GPSend),
         start=date(Date),
         id_period = paste0(CapID,date(Date)))%>%
  select(id=CapID,
         id_period,
         sex=Sex,
         age=age.combine,
         start,
         end,
         reason=VHFendReas,
         fate_unknown=`Fate unknown?`)%>%
  mutate(end=case_when(is.na(end)~study.end,
                       TRUE~end),
         reason=case_when(is.na(reason)~"Alive",
                          TRUE~reason),
         event=case_when(reason%in% "Mortality"~1,
                         TRUE~0))


##make annual
source(here::here("fn", "seasonalsurvival_fn.r"))
surv.yr <- surv.raw%>%
  dplyr::select(id,id_period,sex,start,end,event)%>%
   dplyr::ungroup()%>%
  stretch_survival_data('1 day')%>%
  mutate(year=year(start))%>%
  group_by(id, id_period,sex,year)%>%
  summarise(enter.date=min(start),
            exit.date=max(end)-1,
            event=max(dead))%>%
  mutate(enter=month(enter.date)-1,
         exit=month(exit.date),
         time=(exit.date-enter.date)%>%as.numeric)

##add age
surv.yr <-surv.yr%>%
  left_join(cap%>%
              mutate(year=year(Date))%>%
              select(id=CapID,
                     age=age.combine,
                     year)%>%
              mutate(birthyear=year-age)%>%
              select(id,birthyear)%>%
              distinct(),
            by="id")%>%
  mutate(age=year-birthyear,
         ageclass=case_when(age%in%0:1~"0-1",
                                   age%in%2:6~"2-6",
                                   age>6~">6"))%>%
  dplyr::ungroup()


##summarise
surv.yr%>%
  ungroup%>%
  distinct(id, year,.keep_all=TRUE)%>%
  summarise(n=n(),
            inds=n_distinct(id))
```

## Estimate survival

``` r
##rough by year for R2
surv.yr%>%
  filter(age>2)%>%
  group_by(age,sex)%>%
  summarise(m=sum(event)/sum(time),
            s=(1-m)^365,
            n=sum(time)/365,
            mort=sum(event))%>%
    ggplot(aes(x=age, y=s, color=sex, label=mort))+
    geom_point()+
    geom_path()+
  geom_text(vjust=-1)+
  labs(x="Age", y="Annual survival")+
  theme_ipsum()+
  ylim(0.4,1.05)+
  scale_x_continuous(n.breaks=14)+
    theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/surv%20est-1.png)<!-- -->

``` r
fit <- survfit(Surv(time, event) ~ sex + ageclass, data = surv.yr, conf.type = "log-log")
ggsurvplot(fit, 
           pval = TRUE,
           #risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           risk.table.height = 0.5,
           ggtheme = theme_ipsum(), # Change ggplot2 theme
           ylim = c(0.5,1),
           xlim=c(0,360))
```

![](README_files/figure-gfm/surv%20est-2.png)<!-- -->

``` r
out<-summary(
  survfit(
    Surv(enter, exit, event)~sex + ageclass, 
    conf.type = "log-log",
    data=surv.yr),
  times=12,
  extend = TRUE)

surv.yr.est <-data.frame(strata=out$strata,
                                est=out$surv,
                                se=out$std.err,
                                lower=out$lower,
                                upper=out$upper,
                                n=out$n,
                                sd=out$std.err)%>% ##consistent with Eacker App
  mutate(sex=str_split(strata," ",simplify = TRUE)[,1]%>%str_sub(5,5),
         ageclass=str_split(strata," ",simplify = TRUE)[,2]%>%str_remove("ageclass=")%>%fct_relevel("0-1","2-6",">6"))



##PLOT survival

ggplot(surv.yr.est%>%filter(ageclass%in%c(">6","2-6")) , aes(x=sex, y=est,ymin=lower, ymax=upper, color=ageclass))+
  geom_point(position = position_dodge(width = 0.3))+
  geom_linerange(position = position_dodge(width = 0.3))+
  geom_text(position = position_dodge(width = 0.3),hjust=2, aes(label=n))+
  labs(x="Sex", y="Annual survival", color="Age class")+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/surv%20est-3.png)<!-- -->

## bootstrap to get error distribtion

``` r
surv.boot <- data.frame()
mort.add <- tibble(dplyr::slice(surv.yr%>%ungroup,1))%>%
  dplyr::mutate(event=1,
                enter=0,
                exit=12,
                sex=c("M"),
                ageclass=c(">6"))

set.seed(2022)
for(i in 1:5000){
dat.i <- surv.yr%>%
  bind_rows(mort.add)%>% ### add a mort record to 100% survival
  dplyr::group_by(sex,ageclass)%>%
  dplyr::sample_frac(1,replace=TRUE)

out.i<-summary(
    survfit(
      Surv(enter, exit, event)~sex + ageclass, 
      conf.type = "log-log",
      data=dat.i),
    times=12,
    extend = TRUE)
  
  surv.est.i <-data.frame(strata=out.i$strata,
                           est=out.i$surv)%>% 
    dplyr::mutate(sex=str_split(strata," ",simplify = TRUE)[,1]%>%str_sub(5,5),
           ageclass=str_split(strata," ",simplify = TRUE)[,2]%>%str_remove("ageclass=")%>%fct_relevel("0-1","2-6",">6"),
           iter=i)

  surv.boot <- rbind(surv.boot,surv.est.i)  
}


survival.plot <- ggplot()+
  # stat_dots(data=surv.boot%>%filter(ageclass!="0-2"), 
  #           aes(fill =ageclass, y =est, x=sex),position = position_dodge(width = 0.3),slab_alpha=0.5,justification = 3,quantiles = 50)+
  stat_slab(data=surv.boot%>%filter(ageclass!="0-1"), 
            aes(fill =ageclass, y =est, x=sex),position = position_dodge(width = 0.5),slab_alpha=0.7,justification = 0)+
  geom_pointinterval(data=surv.yr.est%>%filter(ageclass!="0-1"),
                     aes(fill =ageclass, y =est,ymin=est-se,ymax=est+se, x=sex),position = position_dodge(width = 0.5),justification = 0.3)+
  geom_text(data=surv.yr.est%>%filter(ageclass!="0-1"),
            position = position_dodge(width = 0.5),hjust=1.5, aes(label=paste0("n=",n), y =est, x=sex, group=ageclass), size=5)+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Accent")[c(3:2)])+
  labs(x="Sex", y="Annual survival", fill="Age class", title="a) Survival"
       #, subtitle = "by age class, with bootstrapped error distribution"
       )+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = c(0.4,0.2))
survival.plot 
```

![](README_files/figure-gfm/surv.boot-1.png)<!-- -->

# Reproduction

## Prep reproduction data

``` r
cubs <- cubs.raw%>%
  select(id=`Capture ID`,
         Date,
         n.cubs=`Cubs#`,
         cub.ageclass=`Cub Age Class`,
         litter=LitterID,
         n.cubs.end=`Cubs observed at end`)%>%
  mutate(year=year(Date),
         litter=as.character(litter))%>%
  group_by(id,year,litter,cub.ageclass)%>%
  summarise(n.cubs=max(n.cubs, na.rm=TRUE),
            n.cubs.end=mean(n.cubs.end, na.rm=TRUE))%>%
  mutate(cub.age=case_when(cub.ageclass%in%"COY"~0,
                           cub.ageclass%in%"1YO"~1,
                           cub.ageclass%in%"2YO"~2,
                           cub.ageclass%in%"3YO"~3))%>%
  left_join(cap%>%
              mutate(year=year(Date))%>%
              select(id=CapID,
                     age=age.combine,
                     year)%>%
              mutate(birthyear=year-age)%>%
              select(id,birthyear)%>%
              distinct(),
            by="id")%>%
  mutate(mom.age=year-birthyear,
         mom.ageclass=case_when(mom.age%in%0:4~"0-4",
                            mom.age%in%5:6~"5-6",
                            mom.age>6~">6"),
         mom.agecub0=case_when(n.cubs>0~mom.age-cub.age,
                               n.cubs==0~mom.age))

##summarise
cubs%>%
  ungroup%>%
  distinct(id, year,.keep_all=TRUE)%>%
  summarise(n=n(),
            inds=n_distinct(id))

cubs%>%
  ungroup%>%
  distinct(id, year,.keep_all=TRUE)%>%
  drop_na(litter)%>%
  summarise(n=n(),
            inds=n_distinct(id),
            litters=n_distinct(paste0(id,litter)))

cubs%>%
  ungroup%>%
  distinct(id, year,.keep_all=TRUE)%>%
  count(cub.ageclass)

cubs%>%
  ungroup%>%
  distinct(id, year,.keep_all=TRUE)%>%
  group_by(cub.ageclass)%>%
  summarise(n=mean(n.cubs))


cubs%>%
  ungroup%>%
  group_by(id, litter)%>%
  summarise(n=max(n.cubs,na.rm=TRUE))%>%
  ungroup%>%
  summarise(n=sum(n,na.rm=TRUE))


cubs%>%
  ungroup%>%
  drop_na(n.cubs.end)%>%
  filter(cub.ageclass%in%c("COY","1YO"))%>%
  group_by(id, litter)%>%
  summarise(n=max(n.cubs,na.rm=TRUE))%>%
  ungroup%>%
  summarise(n=sum(n,na.rm=TRUE))

cubs%>%
  ungroup%>%
  drop_na(n.cubs.end)%>%
  group_by(cub.ageclass)%>%
  mutate(dead=n.cubs.end-n.cubs)%>%
  summarise(sum(dead), n=sum(n.cubs))

##plot

cubs%>%
  ungroup%>%
  filter(litter==1|is.na(litter),
         mom.agecub0%in%c(2:10))%>%
  mutate(cub=as.numeric(n.cubs>0)/n())%>%
  group_by(mom.agecub0)%>%
  summarise(prop=sum(cub))%>%
  ggplot(aes(x=mom.agecub0%>%as.integer,y=prop))+
  geom_col()+
  labs(x="Age", y="Proportion", title="Age of First Reproduction for Elk Valley Grizzly Bears", subtitle = "2016-2021")+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "bottom")
```

![](README_files/figure-gfm/prep%20repro-1.png)<!-- -->

``` r
cubs%>%
  ungroup%>%
  mutate(cub.ageclass=case_when(is.na(cub.ageclass)~"No cubs",
                                TRUE~cub.ageclass))%>%
  group_by(mom.age)%>%
  mutate(cub=1/n())%>%
  group_by(mom.age,cub.ageclass)%>%
  summarise(prop=sum(cub),
            n=n())%>%
  ungroup%>%
  ggplot(aes(x=mom.age,y=prop, fill=cub.ageclass,label = n))+
  geom_col()+
  geom_text(position = position_stack(vjust = 0.5),size=6)+
  scale_fill_manual(values=RColorBrewer::brewer.pal(6, "Accent")[c(2:6)])+
  labs(x="Age", y="Proportion", fill="Cub age", title="Reproduction by age", subtitle = "2016-2022")+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15),
        legend.position = "bottom")
```

![](README_files/figure-gfm/prep%20repro-2.png)<!-- -->

``` r
ggsave(here::here("plots","repro_age.png"), width=10, height=8, bg="white")


##expand by individual cub
cub.expand <- tibble()
cubs.surv <- cubs%>%drop_na(n.cubs.end)
for(i in 1:nrow(cubs.surv)){
  a <- cubs.surv[i,]
  b <- bind_rows(replicate(n=a[1,"n.cubs"][[1]], a, simplify = FALSE))
  b$fate <- c(rep(1,a[1,"n.cubs.end"]), rep(0,(a[1,"n.cubs"]-a[1,"n.cubs.end"])))
  cub.expand <- bind_rows(cub.expand,b)
}

cub.summary <- cub.expand%>%
  group_by(cub.age)%>%
  summarise(survival=mean(fate),
            n=n(),
            survived=sum(fate))%>%
  mutate(se= sqrt( survival*( 1 - survival) / n ),
         ageclass=as.factor(cub.age))

kable(cub.summary)
```

| cub.age |  survival |   n | survived |        se | ageclass |
|--------:|----------:|----:|---------:|----------:|:---------|
|       0 | 0.7307692 |  26 |       19 | 0.0869893 | 0        |
|       1 | 0.7333333 |  15 |       11 | 0.1141798 | 1        |
|       2 | 0.6666667 |   3 |        2 | 0.2721655 | 2        |

## bootstrap cub surv

``` r
cubsurv.boot<- data.frame()
set.seed(2022)
for(i in 1:5000){
  
cubsurv.boot <-cubsurv.boot%>%
  rbind(
    cub.expand%>%
    dplyr::ungroup()%>%
    dplyr::filter(cub.age<2)%>%
    dplyr::sample_frac(1, replace=TRUE)%>%
    dplyr::summarise(est=mean(fate))%>%
    dplyr::mutate(iter=i)
    )
}

#repro.boot$ageclass <- as.factor(repro.boot$cub.age)
cubsurv.boot$ageclass <- "0-1" ##0-1 similar, pool survival across years then. 


##combine  survival
surv.boot.all <-
cubsurv.boot%>%
  mutate(sex="M")%>%
  rbind(
    cubsurv.boot%>%
          mutate(sex="F")
    )%>%
  rbind(
  surv.boot%>%
    filter(ageclass!="0-1")%>%
    select(est,iter,ageclass,sex)
    )%>%
  mutate(param="survival")


surv.boot.all.summary <-surv.boot.all%>%
  group_by(ageclass,sex,param)%>%
  summarise(median.boot=median(est),
            sd.boot=sd(est),
            lower=quantile(est,0.05),
            upper=quantile(est, 0.95))

write_csv(surv.boot.all.summary, here::here("tables/surv.csv"))

kable(surv.boot.all.summary)
```

| ageclass | sex | param    | median.boot |   sd.boot |     lower |     upper |
|:---------|:----|:---------|------------:|----------:|----------:|----------:|
| 0-1      | F   | survival |   0.7317073 | 0.0680960 | 0.6097561 | 0.8292683 |
| 0-1      | M   | survival |   0.7317073 | 0.0680960 | 0.6097561 | 0.8292683 |
| 2-6      | F   | survival |   0.7152902 | 0.1021722 | 0.5354051 | 0.8770043 |
| 2-6      | M   | survival |   0.5965909 | 0.1358203 | 0.3717150 | 0.8181818 |
| \>6      | F   | survival |   0.9600000 | 0.0285896 | 0.9074074 | 1.0000000 |
| \>6      | M   | survival |   0.9444444 | 0.0579261 | 0.8333333 | 1.0000000 |

## Estimate reproduction rate

``` r
cubs%>%
  filter(mom.ageclass!="0-4")%>%
  mutate(repro=case_when(cub.ageclass=="COY"~n.cubs,
                       TRUE~0))%>%
  group_by(mom.ageclass)%>%
  summarise(repro=mean(repro),
            n=n())%>%
  mutate(se= sqrt( repro*( 1 - repro) / n ))

repro.summary <- cubs%>%
  dplyr::filter(mom.ageclass!="0-4")%>%
  dplyr::mutate(repro=case_when(cub.ageclass=="COY"~n.cubs,
                                TRUE~0))%>%
  dplyr::group_by(mom.ageclass)%>%
  dplyr::summarise(est=mean(repro),
                   n=n())%>%
  mutate(se= sqrt( est*( 1 - est) / n ),)

##bootstrap repro
##0-1 similar, pool survival across years then. 
repro.boot<- data.frame()
set.seed(2022)
for(i in 1:5000){
  
  repro.boot <-repro.boot%>%
    rbind(
      cubs%>%
        dplyr::filter(mom.ageclass!="0-4")%>%
        dplyr::mutate(repro=case_when(cub.ageclass=="COY"~n.cubs,
                               TRUE~0))%>%
        dplyr::group_by(mom.ageclass)%>%
        dplyr::sample_frac(1, replace=TRUE)%>%
        dplyr::summarise(est=mean(repro))%>%
        dplyr::mutate(iter=i)
    )
}

repro.boot<-repro.boot%>%
  mutate(param="reproduction")%>%
  select(est,iter,ageclass=mom.ageclass,param)

repro.boot.summary <-repro.boot%>%
  group_by(ageclass,param)%>%
  summarise(median.boot=median(est/2),
            sd.boot=sd(est/2),
            lower=quantile(est/2,0.05),
            upper=quantile(est/2, 0.95))

write_csv(repro.boot.summary, here::here("tables/repro.csv"))

repro.plot <- ggplot()+
  # stat_dots(data=surv.boot%>%filter(ageclass!="0-2"), 
  #           aes(fill =ageclass, y =est, x=sex),position = position_dodge(width = 0.3),slab_alpha=0.5,justification = 3,quantiles = 50)+
  stat_slab(data=repro.boot, 
            aes(fill =ageclass, y =est/2, x=ageclass%>%fct_relevel("5-6",">6")),slab_alpha=0.7, adjust=2)+
  geom_pointinterval(data=repro.summary,
                     aes(fill =mom.ageclass, y =est/2,ymin=(est-se)/2,ymax=(est+se)/2, x=mom.ageclass))+
  geom_text(data=repro.summary,
            hjust=1.5, aes(label=paste0("n=",n), y =est/2, x=mom.ageclass,  group=mom.ageclass), size=5)+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Accent")[c(2:3)])+
  labs(x="Age class", y="Reproductive rate", fill="", title="b) Reproduction"
       #, subtitle = "Female cubs per year by mother age class, with bootstrapped error distribution"
       )+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = "none")

repro.plot
```

![](README_files/figure-gfm/cub.repro.boot-1.png)<!-- -->

``` r
kable(repro.boot.summary)
```

| ageclass | param        | median.boot |   sd.boot |     lower |     upper |
|:---------|:-------------|------------:|----------:|----------:|----------:|
| 5-6      | reproduction |   0.1538462 | 0.0982936 | 0.0000000 | 0.3076923 |
| \>6      | reproduction |   0.2348485 | 0.0528544 | 0.1515152 | 0.3257576 |

## Age of primiparity estimation from Garshelis et al. 1998

``` r
primipar <- cubs%>%
  ungroup%>%
  filter(mom.agecub0%in%5:9,
         #!id%in%c("Cali","Emma"),
         litter==1 | is.na(litter))

primipar%>%distinct(id)%>%nrow()
```

    ## [1] 16

``` r
##impute records for F observed post COY, but likely first litter
input.par <- primipar%>%
  drop_na(cub.ageclass)%>%
  group_by(id)%>%
  mutate(cub.ageclass=as.factor(cub.ageclass))%>%
  count(cub.ageclass,.drop = FALSE)%>%
  filter(cub.ageclass%in%"COY"& n==0)

primipar.add <-
  primipar%>%
  filter(id%in%input.par$id)%>%
  group_by(id)%>%
  slice(1)%>%
  mutate(cub.ageclass="COY",mom.age=mom.age-cub.age)

primipar <- primipar%>%
  rbind(primipar.add)

##split into those that did and didn't have cubs during period
primipar.cubs <- primipar%>%
  filter(cub.ageclass=="COY")

primipar.end <- primipar%>%
  filter(!id%in%primipar.cubs$id)%>%
  group_by(id)%>%
  filter(mom.age==max(mom.age))%>%
  rbind(primipar.cubs)

##prep datat
primipar.together <- tibble()
bears<- unique(primipar.end$id)
for(i in 1:length(bears)){
  end <- primipar.end%>%
    filter(id%in%bears[i])
  
  cub <- primipar.cubs%>%
    filter(id%in%bears[i])%>%
    mutate(cub=1)%>%
    select(id,mom.age,cub)
  
b <- tibble(id=end$id,
            mom.age=5:end$mom.age)%>%
  dplyr::left_join(cub)%>%
  mutate(cub=replace_na(cub,0))
  
primipar.together <- rbind(primipar.together,b)
}

# ##make into matrix
# primipar.matrix <- primipar.together%>%
#   ungroup()%>%
#   tidyr::spread(mom.age,cub)

##prep Garshelis 1998 stats
primip.summary <- primipar.together%>%
  filter(mom.age<=9)%>%
  group_by(mom.age)%>%
  summarise(null=n(),
          first=sum(cub))%>%
  mutate(prop.cubs=first/null,
         perc.avail=c(100,NA,NA,NA,NA))

primip.summary[2,"perc.avail"]<-primip.summary[1,"perc.avail"]*(1-primip.summary[1,"prop.cubs"])
primip.summary[3,"perc.avail"]<-primip.summary[2,"perc.avail"]*(1-primip.summary[2,"prop.cubs"])
primip.summary[4,"perc.avail"]<-primip.summary[3,"perc.avail"]*(1-primip.summary[3,"prop.cubs"])
primip.summary[5,"perc.avail"]<-primip.summary[4,"perc.avail"]*(1-primip.summary[4,"prop.cubs"])
  
primip.summary$perc.prod <- primip.summary$perc.avail*primip.summary$prop.cubs 

primip.summary <- primip.summary%>%
  mutate(wt=(perc.prod*mom.age)/sum(perc.prod),
         mom.age=as.character(mom.age))%>%
  add_row(mom.age="Total",wt=sum(.$wt))


kable(primip.summary)
```

| mom.age | null | first | prop.cubs | perc.avail | perc.prod |        wt |
|:--------|-----:|------:|----------:|-----------:|----------:|----------:|
| 5       |   16 |     2 | 0.1250000 |  100.00000 |  12.50000 | 0.8940397 |
| 6       |   14 |     1 | 0.0714286 |   87.50000 |   6.25000 | 0.5364238 |
| 7       |    9 |     2 | 0.2222222 |   81.25000 |  18.05556 | 1.8079470 |
| 8       |    7 |     2 | 0.2857143 |   63.19444 |  18.05556 | 2.0662252 |
| 9       |    3 |     1 | 0.3333333 |   45.13889 |  15.04630 | 1.9370861 |
| Total   |   NA |    NA |        NA |         NA |        NA | 7.2417219 |

# Population growth

``` r
##FEMALE ONLY
Px <- c(rep(surv.boot.all.summary[3,"median.boot"][[1]],times=2),rep(surv.boot.all.summary[5,"median.boot"][[1]], times=5),rep(surv.boot.all.summary[1,"median.boot"][[1]], times=19),0)
Fx <- c(rep(0,times=5),rep(repro.boot.summary[2,"median.boot"][[1]], times=2),rep(repro.boot.summary[1,"median.boot"][[1]], times=21))
L <- odiag(Px, -1)
L[1, ] <- Fx

eigen(L)$values[1]
```

    ## [1] 0.9164652+0i

``` r
##boot
lambda.boot <- data.frame()
imi.boot <- data.frame()
set.seed(2022)
for(i in 1:5000){
  a <-surv.boot.all%>%
    dplyr::filter(iter==i & sex=="F")
  
  b <- repro.boot%>%
    dplyr::filter(iter==i & ageclass=="5-6")
  
  c <- repro.boot%>%
    dplyr::filter(iter==i & ageclass==">6")
  
  Px <- c(rep(a%>%dplyr::filter(ageclass=="0-1")%>%pull(est),times=2),rep(a%>%dplyr::filter(ageclass=="2-6")%>%pull(est), times=5),rep(a%>%dplyr::filter(ageclass==">6")%>%pull(est), times=19),0)
  Fx <- c(rep(0,times=5),rep(b$est*0.5, times=2),rep(c$est*0.5, times=21)) ##divided by two to represent females only
  L <- odiag(Px, -1)
  L[1, ] <- Fx
  
  lambda.boot <- rbind(lambda.boot, tibble(l=eigen(L)$values[1], iter=i))
  imi.boot  <-rbind(imi.boot, tibble(imi=lamba.opencr-eigen(L)$values[1], iter=i)) ##lamba.opencr =1.007 and comes from openCR lambda model of SRGBP data 2016-2020
}

lambda.summary <- lambda.boot%>%
  mutate(l=as.numeric(l))%>%
  summarise(median.boot=median(l),
            sd.boot=sd(l),
            upper=quantile(l, 0.95),
            lower=quantile(l,0.05),
            below1=mean(l<1))

kable(lambda.summary)
```

| median.boot |   sd.boot |    upper |     lower | below1 |
|------------:|----------:|---------:|----------:|-------:|
|   0.9389419 | 0.0442223 | 1.008867 | 0.8623848 |  0.923 |

``` r
imi.boot%>%
  mutate(l=as.numeric(imi))%>%
  summarise(median.boot=median(l),
            sd.boot=sd(l),
            upper=quantile(l, 0.95),
            lower=quantile(l,0.05),
            below1=mean(l>0))%>%
  kable()
```

| median.boot |   sd.boot |   upper |      lower | below1 |
|------------:|----------:|--------:|-----------:|-------:|
|   0.0689829 | 0.0442223 | 0.14554 | -0.0009419 | 0.9484 |

``` r
lambda.plot <-ggplot()+
  stat_slab(data=lambda.boot,aes(x=l%>%as.numeric%>%round(2)), justification=-0.05)+
    geom_linerange(data=lambda.boot%>%summarise(upper=quantile(l%>%as.numeric,0.975),
                                                lower=quantile(l%>%as.numeric,0.025)),
                   aes(xmin=lower,xmax=upper,y=0),size=0.7)+
  geom_linerange(data=lambda.boot%>%summarise(upper=quantile(l%>%as.numeric,0.9),
                                              lower=quantile(l%>%as.numeric,0.1)),
                     aes(xmin=lower,xmax=upper,y=0), size=1.5)+
    geom_point(data=lambda.boot%>%summarise(est=median(l)%>%as.numeric),
                   aes(x =est,y=0))+
  geom_vline(xintercept = 1, linetype="dashed")+
  labs(x="Population growth rate", y="Sample density", title="b) Intrinsic Population Growth Rate", subtitle = "without immigration")+
  annotate("text", x = 1.05, y = 0.95, label = "Increasing\npopulation",size=5) +
  annotate("text", x = 0.83, y = 0.95, label = "Decreasing\npopulation", size=5) +
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
        plot.title=element_text(size=20),
        plot.subtitle = element_text(size=15),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = "bottom")
lambda.plot
```

![](README_files/figure-gfm/pop%20growth-1.png)<!-- -->

## CI Data

``` r
##There were two kills from MU 423 without a location in CI data. 
##One was EVGM103 which said “Alpine Trails” and the other was an illegal kill in Brule Creek
##CL added these manually to data
ci <- read_excel(here::here("data","CI","CI_MASTER_25Apr2022.XLSX"))

##add in 0's, code kill types, make spatial
ci <- ci%>%
  filter(GRID_EAST>0 & ZONE_NO>0 & ZONE_NO!="00")%>%
  mutate(x=paste0(GRID_EAST,"000")%>%as.numeric,
         y=paste0(GRID_NORTH,"000")%>%as.numeric,
         date=ymd(KILL_DATE),
         source=case_when(KILL_CODE%in%1~"Hunter",
                          KILL_CODE%in%2~"Human-bear conflict",
                          KILL_CODE%in%3~"Picked-up",
                          KILL_CODE%in%4~"Poaching",
                          KILL_CODE%in%5~"Road",
                          KILL_CODE%in%6~"Rail",
                          KILL_CODE%in%7~"Unknown",
                          KILL_CODE%in%9~"Trapped"))%>%
  filter(!KILL_CODE%in%c(3,9))

##remove bears who were unreported but added into CI database once found via collar
ci<-ci%>%
  filter(!CID_NO%in%(mort.raw%>%filter(`Mortality Reported?`%in%"No")%>%drop_na(`CI Number`)%>%pull(`CI Number`)%>%as.character()))


####Make Bear Data Spatial
Zones<-unique(ci$ZONE_NO)
all <- list()
##Loop to get each zone into correct UTM's and then project to lat long
for(i in 1:length(Zones)){
  a<-ci %>%
    filter(ZONE_NO %in% Zones[i])%>%
    st_as_sf(coords=c("x", "y"), crs=paste("+proj=utm +zone=",Zones[i]," +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",sep=""))%>%
    st_transform(crs=st_crs(sa))
  
  
  all[[i]] <- a
}

##Join all CI Zones into a single layer
ci.spat<-sf::st_as_sf(data.table::rbindlist(all))
#mapview(ci.spat)

##assign in study area
ci.spat <- ci.spat%>%
  st_join(sa%>%mutate(region="Elk Valley")%>%select(region))%>%
  mutate(region=case_when(is.na(region)~"Rest of BC", TRUE~region))

##BC grizz range area
range%>%st_make_valid()%>%
  st_intersection(pnw%>%filter(prov=="BC")%>%st_make_valid())%>%st_area()%>%set_units("km^2") ##764,330 km2
```

    ## Units: [km^2]
    ## [1]      0.0 764330.2      0.0      0.0      0.0

``` r
##EV study area
sa%>%st_area()%>%set_units("km^2") ##5,073 km2
```

    ## 5073.642 [km^2]

``` r
##null expected proportion=5125/(764330)=0.67%

##compare prevelance
ci.compare <- ci.spat%>%
  tibble%>%
  filter(year(date)>2000)%>%
  filter(!CID_NO%in%c("134689","134688","150302","150301","127132","137208"))%>% ##remove collared bears that would have been unreported if not collared
  mutate(source=case_when(source%in%c("Poaching")~"Human-bear conflict",
                             TRUE~source))%>%
  group_by(region)%>%
  count(source)%>%
  ungroup%>%
  pivot_wider(names_from=region, values_from=n)%>%
  mutate(`Provincial total`=`Rest of BC`+`Elk Valley`,
         ev.d=(`Elk Valley`/5.135),
         rbc.d=(`Rest of BC`/(764.33-5.135)),
         bc.d=(`Provincial total`/(764.33)))%>%
  mutate_if(is.numeric, ~replace_na(., 0))%>%
  mutate(`Elk Valley share (%)`=((`Elk Valley`/`Provincial total`)*100)%>%round(0),
         `Excess (x higher than expected)`=(ev.d/rbc.d)%>%round(0))%>%
  mutate(`Elk Valley`=paste0(ev.d%>%round(2), " (",`Elk Valley`,")"),
         `Rest of BC`=paste0(rbc.d%>%round(2), " (",`Rest of BC`,")"))%>%
  select(source:`Rest of BC`,`Elk Valley share (%)`:`Excess (x higher than expected)`)

write_csv(ci.compare, here::here("tables","ci.compare.csv"))

kable(ci.compare)
```

| source              | Elk Valley | Rest of BC  | Elk Valley share (%) | Excess (x higher than expected) |
|:--------------------|:-----------|:------------|---------------------:|--------------------------------:|
| Human-bear conflict | 13.44 (69) | 1.25 (947)  |                    7 |                              11 |
| Hunter              | 14.22 (73) | 5.72 (4340) |                    2 |                               2 |
| Rail                | 3.7 (19)   | 0.03 (26)   |                   42 |                             108 |
| Road                | 3.51 (18)  | 0.05 (37)   |                   33 |                              72 |

## Mortality sources and reporting rates

``` r
###sources of mortality
mort<- mort.raw%>%
  mutate(cause_grouped=case_when(`Suspected Cause of Mortality?`%in%c("Collision (road OR rail), Predation","Rail","Road")~"Road/Rail",
                                 `Suspected Cause of Mortality?`%in%c("Defence of life or property","Poaching")~"Conflict", ##Including Poaching here in Conflicts**
                                 `Suspected Cause of Mortality?`%in%c("unknown-human suspected")~"Unk-human suspected",
                                 TRUE~`Suspected Cause of Mortality?`),
         cause_grouped_cos=case_when(cause_grouped=="Conflict"&`Killed by COS`=="Yes"~"Conflict-COS",
                                     TRUE~cause_grouped))%>%
  
           select(id=`Capture ID`,
                  monitored=`GPS/VHF Functioning?`,
                  fat=`Rump fat thickness (cm)`,
                  lat=`Location (latitude)`,
                  long=`Location (longitude)`,
                  mortdate=`Suspected Date of Mortality`,
                  cause=`Suspected Cause of Mortality?`,
                  cause_grouped,
                  cause_grouped_cos,
                  cos.kill=`Killed by COS`,
                  reported=`Mortality Reported?`
                  )

##all recorded morts of tagged bears, some not actively monitored when killed
mort%>%
  count(cause_grouped)

mort%>%
  count(cause_grouped_cos)

mort%>%
  group_by(cause_grouped)%>%
  drop_na(fat)%>%
  summarise(n=n(),
            mean=mean(fat),
            min=min(fat),
            max=max(fat))
         

##morts that were tagged and summarize reporting       
mort.table <- mort%>%
  filter(monitored=="Yes")%>%         
  count(cause_grouped)%>%
  left_join(
    mort%>%
      filter(monitored=="Yes")%>%
      count(cause_grouped, reported)%>%
      pivot_wider(names_from="reported", values_from="n")
  )%>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0))%>%
  mutate(unreporting.rate=No/(Yes+No))

(sum(mort.table$No)-1)/(sum(mort.table$No+mort.table$Yes)-1)
```

    ## [1] 0.5384615

``` r
##reporting rate with error
report.boot <- c()

for(i in 1:5000){
report.boot[i] <-mort%>%
  dplyr::filter(monitored=="Yes" & cause_grouped!="Natural")%>%
  dplyr::sample_frac(1, replace=TRUE)%>%
  dplyr::summarise(unreported=sum(reported=="No")/n())%>%
  dplyr::pull(unreported)
}
  
mean(report.boot)
```

    ## [1] 0.5373385

``` r
mort%>%
  dplyr::filter(monitored=="Yes" & cause_grouped!="Natural")%>%
  dplyr::summarise(unreported=sum(reported=="No")/n())%>%
  dplyr::pull(unreported)
```

    ## [1] 0.5384615

``` r
quantile(report.boot,0.95)
```

    ##       95% 
    ## 0.7692308

``` r
quantile(report.boot,0.05)
```

    ##        5% 
    ## 0.3076923

``` r
##another way of estimating unreported by cause using ear tags
co.eartag.rate <- mort%>%
  filter(cause_grouped_cos!="Natural" & id!="EVGF111")%>% ##remove "EVGF111" she was killed with Melissa, not independent event
  count(cause_grouped_cos, name="n.all")%>%
  left_join(
    mort%>%
      filter(monitored=="Yes"& cause_grouped_cos!="Natural")%>%   
      group_by(cause_grouped_cos)%>%
      summarise(n.monitored=n(),
                n.reported=sum(reported=="Yes"),
                n.unreported=sum(reported=="No"))
    
  )%>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))%>%
  filter(cause_grouped_cos=="Conflict-COS")%>%
  mutate(cos.rate=n.monitored/(n.all-n.unreported))%>%
  pull(cos.rate)

unreported.v2 <- mort%>%
  filter(cause_grouped_cos!="Natural" & id!="EVGF111")%>% ##remove "EVGF111" she was killed with Melissa, not independent event
  count(cause_grouped_cos, name="n.all")%>%
  left_join(
    mort%>%
      filter(monitored=="Yes"& cause_grouped_cos!="Natural")%>%   
      group_by(cause_grouped_cos)%>%
      summarise(n.monitored=n(),
                n.reported=sum(reported=="Yes"),
                n.unreported=sum(reported=="No"))
                
  )%>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))%>%
  mutate(expected=n.monitored/co.eartag.rate) %>%
  mutate(unreported.rate.v1=(n.unreported/n.monitored),
         unreported.rate.v2=1-((n.all-n.unreported)/expected))

1-((unreported.v2$n.all[5]-unreported.v2$n.unreported[5])/(unreported.v2$n.monitored[5]/co.eartag.rate ))
```

    ## [1] 0.7606838

``` r
##boot
report.boot2 <- tibble()

set.seed(2023)
for(i in 1:5000){
  dat.i  <- mort%>%
    dplyr::group_by(cause_grouped_cos)%>%
    dplyr::filter(cause_grouped_cos!="Natural" & id!="EVGF111")%>% ##remove "EVGF111" she was killed with Melissa, not independent event
    dplyr::sample_frac(1, replace=TRUE)%>%
    dplyr::ungroup()
  
  ratio.i <-dat.i%>%
    dplyr::filter(cause_grouped_cos=="Conflict-COS")%>%
    dplyr::summarise(prop=sum(monitored=="Yes")/n())%>%
    dplyr:: pull(prop)
  
report.boot2 <- dat.i %>%
  dplyr::summarise(n.all=n(),
            n.monitored=sum(monitored=="Yes"),
            n.unreported=sum(reported=="No"))%>%
  dplyr::mutate(expected.all=n.monitored/ratio.i,
                expected=n.monitored/co.eartag.rate,
         unreported.bootall=1-((n.all-n.unreported)/expected.all),
         unreported.boot=1-((n.all-n.unreported)/expected))%>%
  rbind(report.boot2)

}


median(report.boot2$unreported.bootall)
```

    ## [1] 0.7619048

``` r
quantile(report.boot2$unreported.bootall,0.95)
```

    ## 95% 
    ##   1

``` r
quantile(report.boot2$unreported.bootall,0.05)
```

    ##        5% 
    ## 0.5424837

``` r
# median(report.boot2$unreported.boot)
# quantile(report.boot2$unreported.boot,0.9)
# quantile(report.boot2$unreported.boot,0.1)


##3rd way, using CI data and McLellan 2018 method
ev.ci <- ci.spat%>%
  tibble%>%
  filter(region=="Elk Valley", year(date)>2015)%>%
  left_join(read_excel(here::here("data","CI","ci.EV.conf.LS.xlsx"))%>%select(CID_NO,`Killed By`))%>%
  mutate(source.co=case_when(`Killed By`%in%"COS"~"Conflict-COS",TRUE~source),
         type="CI")%>%
  select(source.co,type)%>%
  rbind(mort%>%mutate(type="Collar")%>%filter(!cause_grouped_cos %in%"Natural", monitored=="Yes")%>%select(source.co=cause_grouped_cos,type))%>%
  mutate(source=case_when(source.co=="Conflict-COS"~"COS",TRUE~"human-caused"))

ev.ci.summary <- ev.ci%>%group_by(type)%>%count(source)%>%ungroup
cos.ci <- ev.ci.summary%>%filter(type=="CI"&source=="COS")%>%pull(n)
cos.col <- ev.ci.summary%>%filter(type=="Collar"&source=="COS")%>%pull(n)
hum.ci <- ev.ci.summary%>%filter(type=="CI"&source=="human-caused")%>%pull(n)
hum.col <- ev.ci.summary%>%filter(type=="Collar"&source=="human-caused")%>%pull(n)

#Mclellan 2018 #'s to check
# cos.ci <- 71
# cos.col <- 10
# hum.ci <- 10
# hum.col <- 12

ci.unreported.n <- ((cos.ci/(cos.col/hum.col))-hum.ci)

ci.unreported <- ci.unreported.n/(ci.unreported.n+hum.ci)

##boot approach 3
report.boot3 <- c()

for(i in 1:5000){
  ev.ci.summary.i <- ev.ci%>%dplyr::group_by(type)%>%dplyr::sample_frac(1, replace=TRUE)%>%dplyr::count(source)
  cos.ci <- ev.ci.summary.i%>%dplyr::filter(type=="CI"&source=="COS")%>%pull(n)
  cos.col <- ev.ci.summary.i%>%dplyr::filter(type=="Collar"&source=="COS")%>%pull(n)
  hum.ci <- ev.ci.summary.i%>%dplyr::filter(type=="CI"&source=="human-caused")%>%pull(n)
  hum.col <- ev.ci.summary.i%>%dplyr::filter(type=="Collar"&source=="human-caused")%>%pull(n)

  
  if(nrow(ev.ci.summary.i)==4){
    ci.unreported.n <- ((cos.ci/(cos.col/hum.col))-hum.ci)
    
    report.boot3[i] <- ci.unreported.n/(ci.unreported.n+hum.ci)
  }
  
  if(nrow(ev.ci.summary.i%>%dplyr::filter(type=="Collar"&source=="COS"))==0){ ##no collars killed by COS, add one in there so not infinity
    ci.unreported.n <- ((cos.ci/(1/hum.col))-hum.ci)
    report.boot3[i] <- ci.unreported.n/(ci.unreported.n+hum.ci)
  }
  
}

##anything <0 to 0
report.boot3[report.boot3<0] <- 0

quantile(report.boot3,0.5, na.rm=TRUE)
```

    ##       50% 
    ## 0.6727273

``` r
quantile(report.boot3,0.95, na.rm=TRUE)
```

    ##       95% 
    ## 0.8888889

``` r
quantile(report.boot3,0.05, na.rm=TRUE)
```

    ## 5% 
    ##  0

``` r
##approach 3 for each type
types <- c(unreported.v2$cause_grouped_cos[1:3],"Hunter")
unreported.rate.v3 <- tibble()
ev.ci.summary.grp <- ev.ci%>%
  mutate(source.co=case_when(source.co%in%c("Animal Control","Human-bear conflict","Poaching")~"Conflict",
                             source.co%in%c("Rail","Road")~"Road/Rail",
                             TRUE~source.co))%>%
  group_by(type)%>%count(source.co)%>%ungroup

for(i in 1:length(types)){
cos.ci.grp <- ev.ci.summary.grp%>%filter(type=="CI"&source.co=="Conflict-COS")%>%pull(n)
cos.col.grp <- ev.ci.summary.grp%>%filter(type=="Collar"&source.co=="Conflict-COS")%>%pull(n)
hum.ci.grp <- ev.ci.summary.grp%>%filter(type=="CI"&source.co==types[i])%>%pull(n)
hum.col.grp <- ev.ci.summary.grp%>%filter(type=="Collar"&source.co==types[i])%>%pull(n)
if(sum(hum.col.grp)==0){
  hum.col.grp <-0
}

ci.unreported.n.grp <- ((cos.ci.grp/(cos.col.grp/hum.col.grp))-hum.ci.grp)

ci.unreported.grp <- ci.unreported.n.grp/(ci.unreported.n.grp+hum.ci.grp)

if(ci.unreported.grp<0){ci.unreported.grp <- 0}

unreported.rate.v3 <- rbind(unreported.rate.v3, tibble(Cause=types[i], 
                                                       `CI reported`=hum.ci.grp,
                                                       `unreported (CI method)`=round(ci.unreported.grp,2)))
}


##unrecorded mort plot
uncrecorded.plot.dat <- tibble(Method=rep(c("Collar fates","Return ratio","Return ratio, parial error","CI"),each=5000)%>%as.character(),
       value=c(report.boot,
               report.boot2$unreported.bootall,
               report.boot2$unreported.boot,
               report.boot3),
       i=rep(1:5000,times=4))%>%
  filter(!Method%in%"Return ratio, parial error")

##ensemble
ensemble <- uncrecorded.plot.dat%>%
  group_by(i)%>%
  summarise(mean=mean(value))%>%
  pull(mean)
quantile(ensemble,0.5, na.rm=TRUE)
```

    ##       50% 
    ## 0.6482399

``` r
quantile(ensemble,0.05, na.rm=TRUE)
```

    ##        5% 
    ## 0.4212767

``` r
quantile(ensemble,0.95, na.rm=TRUE)
```

    ##      95% 
    ## 0.793512

``` r
uncrecorded.plot.dat <- uncrecorded.plot.dat%>%
  rbind(tibble(Method="Ensemble",
               value=ensemble,
        i=1:5000))
  
uncrecorded.plot <-ggplot(data=uncrecorded.plot.dat%>%mutate(Method=fct_relevel(Method,"Collar fates", "Return ratio", "CI", "Ensemble")), aes(fill = Method, color = Method, x = value*100)) +
  stat_slab(alpha = .2, adjust=2) +
  stat_pointinterval(position = position_dodge(width = .4, preserve = "single"),.width = c(.66, .95))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Dark2"))+
  scale_color_manual(values=RColorBrewer::brewer.pal(4, "Dark2"))+
  scale_y_continuous(breaks = NULL)+
  labs(x="Unrecorded mortalities (%)", y="Sample density", title="d) Estimated Unrecorded Mortalities"
       #, subtitle = "5000 bootstrapped samples"
       )+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
        strip.text.x = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = c(0.3,0.9),
        legend.direction = "horizontal")+
  xlim(0,100)+
    guides(fill = guide_legend(nrow = 4))


##put into cleaner table summarizing both
unreported.table <- unreported.v2%>%
  mutate(monitored=paste0(n.monitored, " (",n.unreported,")"),
         `tagged returned (reported)`=n.all-n.unreported,
         `tagged expected`=expected)%>%
  select(Cause=cause_grouped_cos,
         monitored,
         `tagged returned (reported)`, 
         `tagged expected`,
         `unreported (collar method)`=unreported.rate.v1,
         `unreported (eartag method)`=unreported.rate.v2)%>%
    add_row(Cause="Hunter",monitored="0 (0)")%>%
  mutate(Cause=factor( as.character(Cause), levels= c("Conflict","Conflict-COS","Road/Rail","Unk-human suspected","Hunter","Total") ))%>%
  arrange(Cause)%>%
  left_join(unreported.rate.v3%>%
              bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))%>%
              mutate(`unreported (CI method)`=case_when(Cause=="Total"~ci.unreported,
                                                        TRUE~`unreported (CI method)`)))%>%
  mutate_if(is.numeric, round,2)%>%
  relocate(`CI reported`, .after =  `tagged expected`)%>%
  mutate_if(is.numeric, ~replace_na(.x, 0))%>%
    select(Cause,
         monitored,
         `CI reported`,
         `tagged returned (reported)`, 
         `tagged expected`,
         `unreported (collar method)`,
         `unreported (CI method)`,
         `unreported (eartag method)`)

write_csv(unreported.table, here::here("tables/unreported.csv"))

kable(unreported.table)
```

| Cause               | monitored | CI reported | tagged returned (reported) | tagged expected | unreported (collar method) | unreported (CI method) | unreported (eartag method) |
|:--------------------|:----------|------------:|---------------------------:|----------------:|---------------------------:|-----------------------:|---------------------------:|
| Conflict            | 4 (2)     |          14 |                          2 |            18.0 |                       0.50 |                   0.50 |                       0.89 |
| Conflict-COS        | 2 (0)     |          14 |                          9 |             9.0 |                       0.00 |                   0.00 |                       0.00 |
| Road/Rail           | 6 (4)     |          11 |                          3 |            27.0 |                       0.67 |                   0.74 |                       0.89 |
| Unk-human suspected | 1 (1)     |           0 |                          0 |             4.5 |                       1.00 |                   0.00 |                       1.00 |
| Hunter              | 0 (0)     |           3 |                          0 |             0.0 |                       0.00 |                   0.00 |                       0.00 |
| Total               | 13 (7)    |          42 |                         14 |            58.5 |                       0.54 |                   0.64 |                       0.76 |

## Plot Others’ Survival

``` r
##boot yellowstone survival
##mort and month counts taken from Schwartz et al. 2006
# gye <- tibble(sex=c("F","F","M","M"),
#               ageclass=c("Subadult","Adult","Subadult","Adult"),
#               mort=c(3,4,2,17),
#               months=c(388,1998,491,1304))%>%
#   mutate(surv=(1-(mort/months))^12,
#          se= sqrt(surv*( 1 - surv) / months),
#          lcl=surv-(1.96*se),
#          ucl=surv+(1.96*se))

surv.others<-read_csv(here::here("data","demog_studies","grizz_survival_others.csv"))%>%
  rbind(surv.boot.all.summary%>%
          ungroup%>%
          filter(ageclass%in%c(">6","2-6"))%>%
          mutate(Ageclass=case_when(ageclass==">6"~"Adult",
                                    ageclass=="2-6"~"Subadult"),
                 Region="Elk Valley, BC",
                 citation="this study")%>%
          select(Region, Sex=sex, Ageclass, S=median.boot, S.lcl=lower, S.ucl=upper, S.se=sd.boot,citation))


##order by lowest subadult survival
order <- surv.others%>%
  filter(Ageclass=="Subadult")%>%
  group_by(Region)%>%
  summarise(S=mean(S))%>%
  arrange(-S)%>%
  pull(Region)

order <- c("Barren-ground, NT",order)  #add barren ground to end, didnt have SA survival

surv.others<-surv.others%>%
  mutate(Region.plot=fct_relevel(Region,order))

compare.plot <- ggplot(surv.others, aes(x=S,y=Region.plot,xmin=S.lcl,xmax=S.ucl,color=Sex))+
  geom_linerange(position = position_dodge(width = 0.4))+
  geom_point(position = position_dodge(width = 0.4))+
  facet_wrap(vars(Ageclass))+
  labs(x="Survival", y="", 
       title="c) Survival comparison"
       #, subtitle="Elk Valley has lowest documented subadult survival in North America"
       )+
  expand_limits(y=0)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.x = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = c(0.15,0.2),
        legend.direction = "vertical")+
  scale_color_manual(values=RColorBrewer::brewer.pal(7, "Accent")[c(5,1)])

##plot together
fig3 <-(survival.plot+repro.plot)/(compare.plot+uncrecorded.plot)
ggsave(plot=fig3, here::here("plots/demog_plate.png"), width=13, height=12, bg="white")
fig3 
```

![](README_files/figure-gfm/others%20surv-1.png)<!-- -->

## Maps

``` r
###Map captures, morts, and telemetry data

##prep some spatial data
##Cities
cities <- st_read(here::here("data/administrative/places.shp"))%>%
  filter(NAME%in%c("SPARWOOD", "FERNIE",  "Jaffray", "CRANBROOK"))%>%
  mutate(Name=str_to_sentence(NAME))%>%
  select(Name,geometry)%>%
  st_transform(4326)%>%
  rbind(tibble(Name="Elkford", Y=50.02, X=-114.921)%>%st_as_sf(coords=c("X","Y"),crs=4326))%>%
  st_transform(3005)
```

    ## Reading layer `places' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/administrative/places.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 787 features and 3 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: -134.9996 ymin: 48.37586 xmax: -114.6224 ymax: 59.92801
    ## Geodetic CRS:  NAD27

``` r
##disturbances
hwy <- st_read(here::here("data/disturbance/hwy.shp"))
```

    ## Reading layer `hwy' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/disturbance/hwy.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 23904 features and 2 fields
    ## Geometry type: LINESTRING
    ## Dimension:     XY
    ## Bounding box:  xmin: 459196.8 ymin: 5099605 xmax: 1021389 ymax: 5768250
    ## Projected CRS: NAD83 / UTM zone 11N

``` r
mines <- st_read(here::here("data/disturbance/mine_dist_ann.shp"))%>%
  filter(year==max(year))
```

    ## Reading layer `mine_dist_ann' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/disturbance/mine_dist_ann.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 105 features and 2 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 647821.7 ymin: 5482934 xmax: 671106 ymax: 5568438
    ## Projected CRS: NAD83 / UTM zone 11N

``` r
##building density
# bd <- rast("/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/ElkValley_Grizzly_Demography_22/data/buildingdensity/bldgdens_30m.tiff")%>%
#   crop(sa%>%st_buffer(30000)%>%st_transform(26911)%>%vect)
# bd[bd==0]<-NA
# bd <- project(bd, sa%>%st_transform(cust.crs))
# writeRaster(bd,"/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/ElkValley_Grizzly_Demography_22/data/buildingdensity/bldgdens_30m_subset.tiff", overwrite=TRUE)

bd <- rast(here::here("data/buildingdensity/bldgdens_30m_subset.tiff"))
  
##inset
##custom crs to keep things straight
cust.crs <- "+proj=aea +lat_0=50 +lon_0=-114.9 +lat_1=49 +lat_2=50.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"
inset <- ggplot() +
  geom_sf(data = pnw%>%st_transform(cust.crs),fill="grey80", color=NA) +
  geom_sf(data = range%>%st_transform(cust.crs), fill="grey60", color=NA) +
  geom_sf(data = pnw%>%st_transform(cust.crs),size=0.5, fill=NA, color="grey30") +
  geom_sf(data=sa%>%st_buffer(30000)%>%st_transform(cust.crs)%>%st_bbox()%>%st_as_sfc, inherit.aes = FALSE, fill=NA, color="white", linewidth=1)+
  geom_sf_text(data=pnw%>%st_transform(cust.crs), aes(label=prov),inherit.aes = FALSE, size=3)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "transparent"),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position = c(0.65,0.075),
        plot.background = element_rect(fill = "transparent",colour = NA))





##get basemap
register_google("AIzaSyCOwGx2D77XOqRgGhKmcb5F4Kt_S61tCLI")
#set_defaults(map_service = "osm", map_type = "terrain_bg")

bmap.big <- basemap_raster(ext=sa%>%st_buffer(30000)%>%st_transform(3857),
                        map_res=1)%>%projectRaster(crs=cust.crs)
```

    ## Loading basemap 'terrain' from map service 'osm_stamen'...

``` r
bmap.small <- basemap_raster(ext=sa%>%st_buffer(10000)%>%st_transform(3857),
                    map_res=1)%>%projectRaster(crs=cust.crs)
```

    ## Loading basemap 'terrain' from map service 'osm_stamen'...

``` r
ggRGB(bmap.big , r=1, g=2, b=3)+
  theme_bw()+
  geom_spatraster(data = bd, na.rm=TRUE)+
  scale_fill_hypso_c()+
  new_scale_fill()+
  geom_sf(data = hwy%>%st_transform(cust.crs), aes(color="Highway"),size=0.7) +
  geom_sf(data = pnw%>%st_transform(cust.crs),size=1, fill=NA, linetype="dashed") +
 # geom_sf(data=sa%>%st_cast("POLYGON")%>%slice(1)%>%st_transform(cust.crs), aes(color="Study area"), size=1,inherit.aes = FALSE, fill=NA,show.legend = F)+
  geom_sf(data=sa%>%st_cast("POLYGON")%>%st_transform(cust.crs), aes(color="Study area"), linewidth=1,inherit.aes = FALSE, fill=NA,show.legend = F)+
    geom_sf(data = mines%>%st_transform(cust.crs),aes(fill="Coal mine"),alpha=0.8, color=NA) +
    scale_color_manual(values = c("Highway"="black", "Study area"="grey90"))+
  scale_fill_manual(values = c("Coal mine"="black"))+
  geom_sf(data=cities%>%st_transform(cust.crs), inherit.aes = FALSE, size=2,pch=21, fill="white", color="black")+
  geom_sf_label(data=cities%>%filter(!Name%in%c("Sparwood","Elkford"))%>%st_transform(cust.crs), aes(label=Name),inherit.aes = FALSE, size=4,hjust = -0.1, vjust = 1)+
  geom_sf_label(data=cities%>%filter(Name%in%c("Sparwood","Elkford"))%>%st_transform(cust.crs), aes(label=Name),inherit.aes = FALSE, size=4, vjust = 1,hjust=1.1)+
  annotation_custom(ggplotGrob(inset), xmin =92.2E4, xmax = 98.5E4, ymin = -0.5E4, ymax = 6.5E4)+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=15, color="white"),
        legend.title=element_blank(),
        legend.spacing.y = unit(0.05, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.position = c(0.17,0.19))+
  ggsn::scalebar(x.min=93E4, x.max=105E4,y.min=-12.3E4,y.max=5E4, dist = 20,  height=0.02, dist_unit="km",transform=FALSE, location="bottomleft", st.color = "black",st.bottom = FALSE)+
  scale_y_continuous(expand = c(0,NA), limits = c(-12.5E4,6E4))+
  scale_x_continuous(expand = c(0,0),limits = c(92.5E4,105.5E4))+
  annotate("text",x=95E4,y=-9.6E4,label="Building\ndensity\n/3 sq.km", size=4.7, color="white",hjust = 0,lineheight = 1)+
  guides(fill = guide_legend(order = 2),
         color = guide_legend(order = 1, reverse=TRUE))
```

![](README_files/figure-gfm/maps-1.png)<!-- -->

``` r
ggsave(here::here("plots/sa_map.boot.png"), width=8, height=10, bg="white")



ggRGB(bmap.small, r=1, g=2, b=3)+
  theme_bw()+
  geom_sf(data=sa%>%st_transform(cust.crs), inherit.aes = FALSE, fill=NA, color="grey90", linewidth=1)+
  geom_sf(data=telem%>%filter(DateTime>ymd("2016-08-14"))%>%st_transform(cust.crs),inherit.aes = FALSE, alpha=0.1, size=0.1)+
  geom_sf(data=cities%>%st_transform(cust.crs), inherit.aes = FALSE, size=3,pch=21, fill="white", color="black")+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=7),
        legend.title=element_text(size=0),
        legend.position = "bottom")+
  labs(title="Telemetry",
       subtitle = "2016-2021, 63 individuals, 93,623 relocations",
       fill="")+
  scale_x_continuous(expand = c(0,0), limits = c(943588.1,1029124))+
  scale_y_continuous(expand = c(0,0), limits = c(-104054.3,50787.66))
```

![](README_files/figure-gfm/maps-2.png)<!-- -->

``` r
cap.plot <- ggRGB(bmap.small, r=1, g=2, b=3)+
  theme_bw()+
  geom_sf(data=sa%>%st_transform(cust.crs), inherit.aes = FALSE, fill=NA, color="grey90", linewidth=0.8)+
  geom_sf(data=cap%>%st_as_sf(coords=c("long","lat"), crs=4326)%>%st_transform(cust.crs),aes(fill=CapReason), pch=21,inherit.aes = FALSE, size=5, alpha=0.8)+
  geom_sf(data=cities%>%st_transform(cust.crs), inherit.aes = FALSE, size=3,pch=21, fill="white", color="black")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Accent")[c(2:3)])+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
          strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = "bottom")+
  labs(title="a) Captures",
       fill="Capture reason")+
  scale_x_continuous(expand = c(0,0), limits = c(943588.1,1029124))+
  scale_y_continuous(expand = c(0,0), limits = c(-104054.3,50787.66))+
    guides(fill = guide_legend(nrow = 2))

mort.plot <- ggRGB(bmap.small, r=1, g=2, b=3)+
  theme_bw()+
  geom_sf(data=sa%>%st_transform(cust.crs), inherit.aes = FALSE, fill=NA, color="grey90", linewidth=0.8)+
  geom_sf(data=telem%>%filter(DateTime>ymd("2016-08-14"))%>%st_transform(cust.crs),inherit.aes = FALSE, alpha=0.4, size=0.2)+
  geom_sf(data=mort%>%st_as_sf(coords=c("long","lat"), crs=4326)%>%st_transform(cust.crs)%>%st_jitter(1000),aes(color=cause_grouped,shape=monitored),inherit.aes = FALSE, size=5,alpha=0.9)+
  geom_sf(data=cities%>%st_transform(cust.crs), inherit.aes = FALSE, size=3,pch=21, fill="white", color="black")+
  scale_color_manual(values=RColorBrewer::brewer.pal(5, "Accent")[c(3,2,5,4,6)])+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position="bottom", 
        legend.box="vertical",
        legend.margin=margin())+
  labs(title="b) Telemetry and Mortality",
       color="Cause",
       shape="Monitored at time?")+
  scale_x_continuous(expand = c(0,0), limits = c(943588.1,1029124))+
  scale_y_continuous(expand = c(0,0), limits = c(-104054.3,50787.66))+
  guides(title.position = "top",
         color = guide_legend(nrow = 2))

##plot together
fig2 <- cap.plot+mort.plot+(weight.plot/fat.plot)+plot_layout(widths = c(1.6,1.6,1.3))
ggsave(plot=fig2, here::here("plots/cap_plate.png"), width=15, height=10.5, bg="white")



#####Conflict data#####
conflict <- conflict.raw%>%
  filter(encounter_lat>47,encounter_lat<60, encounter_lng>(-140))%>%
  st_as_sf(coords=c("encounter_lng","encounter_lat"), crs=4326)%>%
  st_transform(26911)

grid <- conflict %>% 
  st_buffer(40000)%>%
  st_make_grid(cellsize = 100000, what = "polygons")%>%
  st_sf()%>%
  st_intersection(pnw%>%filter(prov=="BC")%>%st_transform(26911)%>%st_make_valid())%>%
  mutate(count=lengths(st_intersects(.,conflict)),
         area=st_area(.)%>%set_units(km^2),
         conf.dens=as.numeric((count/6)/(area/10000)))

#mapview(grid["conf.dens"])

mean(grid$conf.dens)
```

    ## [1] 5.785677

``` r
conflict.ev <- conflict%>%
  st_transform(st_crs(sa))%>%
  st_intersection(sa)

conf.seas <- tibble(month=1:12, lab=month.abb[month])%>%
  left_join(conflict.ev%>%
              tibble%>%
  mutate(month=month(encounter_date))%>%
  group_by(month)%>%
  count)%>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0))%>%
  mutate(lab=fct_reorder(lab,month))%>%
  ggplot(aes(x=lab,y=n/(2021-2016)))+ ##scale by N years observed
  geom_col()+
  labs(x="Month", y="Human-bear conflict reports", title="b) Seasonal conflict trend")+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=16, angle=45,hjust=1),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = "bottom")


conf.ev <- ggRGB(bmap.small, r=1, g=2, b=3)+
  theme_bw()+
  geom_sf(data=sa%>%st_transform(cust.crs), inherit.aes = FALSE, fill=NA, color="grey90", linewidth=0.8)+
  geom_sf(data=conflict%>%st_transform(cust.crs)%>%st_jitter(1000),inherit.aes = FALSE, size=2,alpha=0.6)+
  geom_sf(data=cities%>%st_transform(cust.crs), inherit.aes = FALSE, size=3,pch=21, fill="white", color="black", alpha=0.6)+
  scale_color_manual(values=RColorBrewer::brewer.pal(5, "Accent")[c(3,2,4,5)])+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position="bottom", 
        legend.box="vertical", 
        legend.margin=margin())+
  labs(title="a) Human-bear conflicts")+
  scale_x_continuous(expand = c(0,0), limits = c(943588.1,1029124))+
  scale_y_continuous(expand = c(0,0), limits = c(-104054.3,35787.66))

conf.bc <- ggplot()+
  theme_bw()+
  #geom_sf(data=sa%>%st_transform(cust.crs), inherit.aes = FALSE, fill=NA, color="grey90",size=0.8)+
  geom_sf(data=grid%>%st_transform(cust.crs), aes(fill=conf.dens))+
  scale_fill_viridis_c()+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
       strip.text.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position="bottom", 
        legend.box="vertical", 
        legend.margin=margin())+
  labs(title="c) Provincial conflict hotspots",
       fill="Conflicts per 10,000 sq.km/year")

fig3 <- conf.ev+(conf.seas/conf.bc)
ggsave(plot=fig3, here::here("plots/conf_plate.png"), width=13, height=11, bg="white")






#####Known Dispersers#####

#mort data
mort.disp <-mort%>%
  filter(id%in% c("Matt", "Donnie", "Sid"))%>%
  st_as_sf(coords=c("long","lat"),
           crs = 4326)%>%
  select(id)%>%
  st_transform(cust.crs)%>%
  cbind(st_coordinates(.))


###known far dispersers  
disp <-telem%>%
  filter(Name%in% c("Matt", "Donnie", "Sid"))%>%
  st_transform(crs = 4326)%>%
  mutate(type="dispersed",
         Date=ymd_hms(DateTime))%>%
  select(id=Name, Date, type)

##Donnie AB data
disp <- disp%>%
  rbind(
    st_read(here::here("data","telem","dispersbears","Bear163", "locs_2016.shp"))%>%
      mutate(id="Donnie",
             Date=ymd_hms(DateObs),
             type="natal")%>%
      st_transform(crs = 4326)%>%
      select(id, Date, type))%>%
  rbind(
    st_read(here::here("data","telem","dispersbears","Bear163", "locs_2017.shp"))%>%
      mutate(id="Donnie",
             Date=ymd_hms(DateObs),
             type="natal")%>%
      st_transform(crs = 4326)%>%
      select(id, Date, type))
```

    ## Reading layer `locs_2016' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/telem/dispersbears/Bear163/locs_2016.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 252 features and 48 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: -115.1735 ymin: 50.59842 xmax: -115.0693 ymax: 50.92693
    ## Geodetic CRS:  WGS 84
    ## Reading layer `locs_2017' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/telem/dispersbears/Bear163/locs_2017.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 62 features and 62 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: -115.1842 ymin: 50.62482 xmax: -115.109 ymax: 50.71597
    ## Geodetic CRS:  WGS 84

``` r
#Add Sid AB data
disp <- disp%>%
  rbind(
    st_read(here::here("data","telem","dispersbears","Bear162", "Bear162.shp"))%>%
      mutate(id="Sid",
             Date=ymd(LMT_Date),
             type="natal")%>%
      st_transform(crs = 4326)%>%
      cbind(st_coordinates(.))%>%
  mutate(type=case_when(Y<50.5~"dispersed",
                        TRUE~type))%>%
    select(id, Date, type)
  )
```

    ## Reading layer `Bear162' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/telem/dispersbears/Bear162/Bear162.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 896 features and 32 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: -1.797693e+308 ymin: -1.797693e+308 xmax: -114.8253 ymax: 50.78919
    ## Geodetic CRS:  WGS 84

``` r
#Add Matt USA data
disp <- disp%>%
  rbind(
    st_read(here::here("data","telem","dispersbears","GB810_MattsMom_GPSVHFLocs", "GB810_MattsMom_GPSVHFLocs.shp"))%>%
      mutate(id="Matt",
             Date=ymd(Date),
             type="natal")%>%
      st_transform(crs = 4326)%>%
      select(id, Date, type))%>%
  rbind(
    st_read(here::here("data","telem","dispersbears","GB18986_Matt_DNAHits", "GB18986_Matt_DNAHits.shp"))%>%
      mutate(id="Matt",
             Date=ymd(Date),
             type="natal")%>%
      st_transform(crs = 4326)%>%
      select(id, Date, type))
```

    ## Reading layer `GB810_MattsMom_GPSVHFLocs' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/telem/dispersbears/GB810_MattsMom_GPSVHFLocs/GB810_MattsMom_GPSVHFLocs.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 1408 features and 18 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: -116.1643 ymin: 48.78457 xmax: -115.9272 ymax: 48.95825
    ## Geodetic CRS:  WGS 84
    ## Reading layer `GB18986_Matt_DNAHits' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/PDF/PDF Analyses/Published/ElkValley_Grizzly_Demography_22/data/telem/dispersbears/GB18986_Matt_DNAHits/GB18986_Matt_DNAHits.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 16 features and 16 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: -116.1262 ymin: 48.81362 xmax: -116.02 ymax: 48.935
    ## Geodetic CRS:  WGS 84

``` r
disp <- disp%>%
  st_transform(cust.crs)%>%
  cbind(st_coordinates(.))%>%
  drop_na(X)%>%
  mutate(id=case_when(id%in%"Donnie"~"EVGM61",
                      id%in%"Matt"~"EVGM79",
                      id%in%"Sid"~"EVGM58"))


bear.line <- disp%>%group_by(id)%>%arrange(Date)%>%summarize(do_union=FALSE) %>% st_cast("LINESTRING")
bmap.imi <- basemap_raster(ext=disp%>%st_buffer(30000)%>%st_transform(3857),
                           map_res=1)%>%projectRaster(crs=cust.crs)
```

    ## Loading basemap 'terrain' from map service 'osm_stamen'...

``` r
imi.map <- ggRGB(bmap.imi, r=1, g=2, b=3)+
  theme_bw()+
  geom_sf(data=sa%>%st_transform(cust.crs), inherit.aes = FALSE, fill=NA, color="grey90", linewidth=0.8)+
  geom_sf(data=disp%>%st_transform(cust.crs),inherit.aes = FALSE, alpha=0.1, size=0.1,aes(color=id))+
  geom_sf(data=bear.line%>%st_transform(cust.crs),inherit.aes = FALSE, alpha=0.8, size=1,aes(color=id))+
  geom_sf(data=cities%>%st_transform(cust.crs), inherit.aes = FALSE, size=3,pch=21, fill="white", color="black")+
  theme_ipsum()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title=element_text(size=20),
        plot.subtitle = element_text(size=15),
        legend.title=element_blank(),
        legend.spacing.y = unit(0.05, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.position = "none")+
  labs(title="a) Known immigrants into study area",
       subtitle="3 male grizzly bears, 77-95 km displacements")+
  scale_y_continuous(expand = c(0,NA))+
  ggrepel::geom_label_repel(
    data = disp%>%filter(type=="natal")%>%as_tibble()%>%select(id,X,Y)%>%group_by(id)%>%summarise_all(mean)%>%mutate(lab="Natal Range")%>%
      rbind(mort.disp%>%as_tibble()%>%select(id,X,Y)%>%mutate(lab="Mortality")),
    aes(label = lab, x=X, y=Y, color=lab),
    segment.size  = 1,
    nudge_x=25000,
    size=4)+
  scale_color_manual(values=c(RColorBrewer::brewer.pal(7, "Accent")[c(3,4,6)],"firebrick3","salmon4"))+
  ggsn::scalebar(x.min=97E4, x.max=102E4,y.min=-16E4,y.max=13.7E4, dist = 20,  height=0.02, dist_unit="km",transform=FALSE, location="bottomright", st.color = "black",st.bottom = FALSE)
imi.map
```

![](README_files/figure-gfm/maps-3.png)<!-- -->

## Plot SCR Density and simulate trends

``` r
start.year <- 2016
secr.pred <- read_csv(here::here("data","caprecap","dens.pred.annual.csv"))%>%
  mutate(se=D.se)%>%
  select(year,N,se,N.lcl,N.ucl)%>%
  mutate(method="DNA",
         type="observed")%>%
  filter(year>=start.year)



col.projection <-tibble()
for(i in 1:5000){
col.projection <- tibble(year=start.year:2040)%>%
  dplyr::mutate(Ncol=(secr.pred[1,2][[1]]*(lambda.boot[i,1][[1]]^(year-start.year)))%>%as.numeric)%>%
  dplyr::select(year,Ncol)%>%
  dplyr::mutate(method="Collar",
                type=case_when(year<2023~"observed",TRUE~"projected"))%>%
  rbind(col.projection)
}

#add in collar
secr.pred <-secr.pred%>%
  rbind(col.projection%>%
          group_by(year,method,type)%>%
          summarise(N=median(Ncol),
                    se=sd(Ncol),
                    N.lcl=quantile(Ncol,0.05),
                    N.ucl=quantile(Ncol,0.95))%>%
          select(year, N, se,N.lcl, N.ucl, method, type))

dna.projection <-tibble()
for(i in 1:5000){
  dna.lambda <- rnorm(1,mean=1,sd=0.01)
  dna.projection <- tibble(year=2022:2040)%>%
    dplyr::mutate(N=(secr.pred%>%dplyr::filter(year==2021 & method=="DNA")%>%pull(N)*(dna.lambda^(year-start.year)))%>%as.numeric)%>%
    dplyr::select(year,N)%>%
    dplyr::mutate(method="DNA",
                  type="projected")%>%
    rbind(dna.projection)
}

#add in DNA projections
secr.pred <-secr.pred%>%
  rbind(dna.projection%>%
          group_by(year,method,type)%>%
          summarise(N.lcl=quantile(N,0.05),
                    se=sd(N),
                    N.ucl=quantile(N,0.95),
                    N=median(N))%>%
          select(year, N, se,N.lcl, N.ucl, method, type))

##immigrants/yr
imi.peryear <- read_csv(here::here("data","caprecap","dens.pred.combined.csv"))%>%
  slice(1)%>%
  mutate(imi=N*median(imi.boot$imi))%>%
  pull(imi)%>%
  as.numeric()%>%
  round(0)

abund.proj <- ggplot(secr.pred)+
  geom_cloud(aes(x=year, y=N, ymin=N-se, ymax=N+se, fill=method),steps=20,alpha=0.9)+
  geom_line(aes(x=year, y=N, ymin=N-se, ymax=N+se, fill=method, linetype=type))+
  theme_ipsum()+
  labs(x="Year", y="N", fill="Method", linetype="Type",
       title="c) Abundance", subtitle="Genetic capture recapture shows abundance is stable,\ncollars show population would decline without immigration")+
  expand_limits(y=0)+
  theme_ipsum()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17),
        strip.text.y = element_text(size=20),
            plot.title=element_text(size=20),
        legend.text = element_text(size=18),
        legend.title=element_text(size=19),
        legend.position = c(0.4,0.1),
        legend.direction = "horizontal")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Accent")[c(2:3)])+
  coord_cartesian(ylim=c(0,110))+
  annotate("segment", x = 2038, xend = 2038, y =secr.pred%>%filter(year==2038, method=="DNA")%>%pull(N), yend = 65,
           colour = "black")+
  annotate("segment", x = 2038, xend = 2038, y =45, yend = secr.pred%>%filter(year==2038, method=="Collar")%>%pull(N),
           colour = "black")+
  annotate("segment", x = 2037.5, xend = 2038.5, y =secr.pred%>%filter(year==2038, method=="DNA")%>%pull(N), yend = secr.pred%>%filter(year==2038, method=="DNA")%>%pull(N),
           colour = "black")+
  annotate("segment", x = 2037.5, xend = 2038.5, y =secr.pred%>%filter(year==2038, method=="Collar")%>%pull(N), yend = secr.pred%>%filter(year==2038, method=="Collar")%>%pull(N),
           colour = "black")+
annotate("text", x = 2037,  y =60,
         label = paste0("Persistence subsidy from\nimmigration & connectivity\n~",imi.peryear, " grizzly bears\nper year"),
         colour = "black",
         hjust = 1)+
  guides(linetype = "none")

fig6 <- imi.map+(lambda.plot/abund.proj)+plot_layout(widths=c(1.5,1))
ggsave(plot=fig6, here::here("plots/imi_plate.png"), width=12, height=11, bg="white")
fig6
```

![](README_files/figure-gfm/DNA%20sim-1.png)<!-- -->
