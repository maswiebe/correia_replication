library(tidyverse)
library(skimr)
library(fixest)
library(specr)
library(tabulator)
library(robomit)
library(vtable)

setwd("/home/michael/Dropbox/replications/air_pollution/correia-pollution/data/")

df <- read_csv('data.csv') %>% 
  filter(is.na(pm25)==0)

#-------------------------------------------------------------------------------

dict = c(
  loginc = 'Log income', logpop = 'Log population', hisp = 'Hispanic share', blk = 'Black share', edu = 'Education', pm25 = 'PM2.5', mort_cr = 'Mortality rate: cardio-respiratory deaths', mort_non_cr = 'Mortality rate: non-cardio-respiratory deaths', mort_all = 'Mortality rate: all causes', mort_canc = 'Mortality rate: cancers', mort_inj = 'Mortality rate: injuries', mort_acc = "Mortality rate: accidents", fips = 'County', stfips = 'State', year = 'Year', region = 'Region', d_urban = 'Urban',
  note = dsb('*Note*: *Mortality rate* is per 100,000 population. *PM2.5* is a county-year average, measured in $\\mu \\text{g} / \\text{m}^{3}$. *Log income* is inflation-adjusted per-capita county income. *Education* is the share with a high school degree or above. *Region* is the 8 BEA regions. Standard errors (in parentheses) are clustered at the state level. Significance levels: $\\text{* } p<0.1, \\text{** } p<0.05, \\text{*** } p<0.01$.'))
setFixest_dict(dict)

my_style = style.tex('aer', tpt=T, notes.tpt.intro = "\\footnotesize", signif.code=c("***"=0.01, "**"=0.05, "*"=0.10)) # signif.code seems to be disabled for 'aer'
setFixest_etable(style.tex = my_style, page.width = 'fit', fitstat = ~ n + r2, view.cache=T)

# note: need to install packages tinytex and pdftools to view png
# https://lrberge.github.io/fixest/articles/etable_new_features.html

#-------------------------------------------------------------------------------

# region-year FEs
df <- df %>% 
  mutate(
    regions = case_when(
      region==1 ~ 'New England',
      region==2 ~ 'Mideast',
      region==3 ~ 'Great Lakes',
      region==4 ~ 'Plains',
      region==5 ~ 'Southeast',
      region==6 ~ 'Southwest',
      region==7 ~ 'Rocky Mountain',
      region==8 ~ 'Far West',
    )
  )

# interaction with urban; heterogeneous TE by urban/rural
# note; urban is from 2010 census, not time-varying.

# correia does: subsample with urban rate >90%; n=169 counties

df <- df %>% 
  mutate(
    d_urban = ifelse(urban > 0.9,1,0)
  )

# separate datasets
df_urban <- df %>%
  filter(d_urban==1)
df_rural <- df %>%
  filter(d_urban==0)


### long difference: keep every 5th year
years5 <- c(1999,2004,2009,2014)
df_long <- df %>% 
  filter(year %in% years5)

df_long_urban <- df_urban %>% 
  filter(year %in% years5)
df_long_rural <- df_rural %>% 
  filter(year %in% years5)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### summary statistics
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sumtable(df, vars=c('mort_all', 'mort_cr', 'mort_canc', 'mort_acc', 'pm25', 'loginc', 'logpop', 'hisp', 'blk', 'edu', 'urban'))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### full sample results
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# all-cause
est_all <- feols(mort_all ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df)
est_all_ry <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df)
est_all_sy <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)

etable(est_all, est_all_ry, est_all_sy,
       notes = 'note',  
       view=TRUE)

# cardio-respiratory
est_cr <- feols(mort_cr ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df)
est_cr_ry <- feols(mort_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df)
est_cr_sy <- feols(mort_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)

etable(est_cr, est_cr_ry, est_cr_sy,
       notes = 'note',
       view=TRUE)

# non-cardio-respiratory
est_ncr <- feols(mort_non_cr ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df)
est_ncr_ry <- feols(mort_non_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df)
est_ncr_sy <- feols(mort_non_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)

etable(est_ncr, est_ncr_ry, est_ncr_sy,
       notes = 'note',
       view=TRUE)

# cancer
est_canc <- feols(mort_canc ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df)
est_canc_ry <- feols(mort_canc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df)
est_canc_sy <- feols(mort_canc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)

etable(est_canc, est_canc_ry, est_canc_sy,
       notes = 'note',
       view=TRUE)

# accidents
est_acc <- feols(mort_acc ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df)
est_acc_ry <- feols(mort_acc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df)
est_acc_sy <- feols(mort_acc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)

etable(est_acc, est_acc_ry, est_acc_sy,
       notes = 'note',
       view=TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### all-cause mortality, urban vs rural: urban rate
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#----------
# urban
est_all_urb <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + year, cluster=~stfips, df, split=df$d_urban)

etable(est_all_urb,
       notes = 'note',
       view=TRUE)

feols(mort_all ~ pm25 + factor(d_urban):pm25 + loginc + logpop + hisp + blk + edu | fips + year, cluster=~stfips, df)
feols(mort_all ~ pm25 + factor(d_urban):pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)
feols(mort_all ~ i(d_urban, pm25) + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)
# when restricting coefficients on covariates (using full sample), get similar effects: -2 in rural, 5 in urban

### add covariates
est_all_cov_urb <- feols(mort_all ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_urban)
est_all_ry_cov_urb <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_urban)
est_all_sy_cov_urb <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_urban)

etable(est_all_cov_urb, est_all_ry_cov_urb, est_all_sy_cov_urb,
       notes = 'note',
       view=TRUE)

df_urban %>% skim(pm25)
df_urban %>% skim(mort_all)
# interpret: increasing pm25 by 10mg (100% of mean) increases mortality rate by ~5% (=41.42/852)

#----------
# rural

### add covariates
est_all_cov_rur <- feols(mort_all ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_rural)
est_all_ry_cov_rur <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_rural)
est_all_sy_cov_rur <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_rural)

etable(est_all_cov_rur, est_all_ry_cov_rur, est_all_sy_cov_rur,
       notes = 'note',
       view=TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### cancer mortality, urban vs rural: urban rate
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#----------
# urban
est_canc_urb <- feols(mort_canc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + year, cluster=~stfips, df, split=df$d_urban)

etable(est_canc_urb,
       notes = 'note',
       view=TRUE)

### add covariates
est_canc_cov_urb <- feols(mort_canc ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_urban)
est_canc_ry_cov_urb <- feols(mort_canc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_urban)
est_canc_sy_cov_urb <- feols(mort_canc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_urban)

etable(est_canc_cov_urb, est_canc_ry_cov_urb, est_canc_sy_cov_urb,
       notes = 'note',
       view=TRUE)

#----------
# rural

### add covariates
est_canc_cov_rur <- feols(mort_canc ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_rural)
est_canc_ry_cov_rur <- feols(mort_canc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_rural)
est_canc_sy_cov_rur <- feols(mort_canc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_rural)

etable(est_canc_cov_rur, est_canc_ry_cov_rur, est_canc_sy_cov_rur,
       notes = 'note',
       view=TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### cardio-respiratory mortality, urban vs rural: urban rate
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#----------
# urban
est_cr_urb <- feols(mort_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + year, cluster=~stfips, df, split=df$d_urban)

etable(est_cr_urb,
       notes = 'note',
       view=TRUE)

### add covariates
est_cr_cov_urb <- feols(mort_cr ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_urban)
est_cr_ry_cov_urb <- feols(mort_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_urban)
est_cr_sy_cov_urb <- feols(mort_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_urban)

etable(est_cr_cov_urb, est_cr_ry_cov_urb, est_cr_sy_cov_urb,
       notes = 'note',
       view=TRUE)

#----------
# rural

### add covariates
est_cr_cov_rur <- feols(mort_cr ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_rural)
est_cr_ry_cov_rur <- feols(mort_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_rural)
est_cr_sy_cov_rur <- feols(mort_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_rural)

etable(est_cr_cov_rur, est_cr_ry_cov_rur, est_cr_sy_cov_rur,
       notes = 'note',
       view=TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### non-cardio-respiratory mortality, urban vs rural: urban rate
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#----------
# urban
est_ncr_urb <- feols(mort_non_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + year, cluster=~stfips, df, split=df$d_urban)

etable(est_ncr_urb,
       notes = 'note',
       view=TRUE)

### add covariates
est_ncr_cov_urb <- feols(mort_non_cr ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_urban)
est_ncr_ry_cov_urb <- feols(mort_non_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_urban)
est_ncr_sy_cov_urb <- feols(mort_non_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_urban)

etable(est_ncr_cov_urb, est_ncr_ry_cov_urb, est_ncr_sy_cov_urb,
       notes = 'note',
       view=TRUE)

#----------
# rural

### add covariates
est_ncr_cov_rur <- feols(mort_non_cr ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_rural)
est_ncr_ry_cov_rur <- feols(mort_non_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_rural)
est_ncr_sy_cov_rur <- feols(mort_non_cr ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_rural)

etable(est_ncr_cov_rur, est_ncr_ry_cov_rur, est_ncr_sy_cov_rur,
       notes = 'note',
       view=TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### accident mortality, urban vs rural: urban rate
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#----------
# urban
est_acc_urb <- feols(mort_acc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + year, cluster=~stfips, df, split=df$d_urban)

etable(est_acc_urb,
       notes = 'note',
       view=TRUE)

feols(mort_acc ~ pm25 + factor(d_urban):pm25 + loginc + logpop + hisp + blk + edu | fips + year, cluster=~stfips, df)
feols(mort_acc ~ pm25 + factor(d_urban):pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df)

### add covariates
est_acc_cov_urb <- feols(mort_acc ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_urban)
est_acc_ry_cov_urb <- feols(mort_acc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_urban)
est_acc_sy_cov_urb <- feols(mort_acc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_urban)

etable(est_acc_cov_urb, est_acc_ry_cov_urb, est_acc_sy_cov_urb,
       notes = 'note',
       view=TRUE)

#----------
# rural

### add covariates
est_acc_cov_rur <- feols(mort_acc ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_rural)
est_acc_ry_cov_rur <- feols(mort_acc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_rural)
est_acc_sy_cov_rur <- feols(mort_acc ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_rural)

etable(est_acc_cov_rur, est_acc_ry_cov_rur, est_acc_sy_cov_rur,
       notes = 'note',
       view=TRUE)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### 'long difference': use every fifth year
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# urban
est_long_all_cov_urb <- feols(mort_all ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_long_urban)
est_long_all_ry_cov_urb <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_long_urban)
est_long_all_sy_cov_urb <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_long_urban)

etable(est_long_all_cov_urb, est_long_all_ry_cov_urb, est_long_all_sy_cov_urb,
       notes = 'note', 
       view=TRUE)

# rural
est_long_all_cov_rur <- feols(mort_all ~ pm25 + csw0(loginc, logpop, hisp, blk, edu) | fips + year, cluster=~stfips, df_long_rural)
est_long_all_ry_cov_rur <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + region^year, cluster=~stfips, df_long_rural)
est_long_all_sy_cov_rur <- feols(mort_all ~ pm25 + loginc + logpop + hisp + blk + edu | fips + stfips^year, cluster=~stfips, df_long_rural)

etable(est_long_all_cov_rur, est_long_all_ry_cov_rur, est_long_all_sy_cov_rur,
       notes = 'note', 
       view=TRUE)
