library(tidyverse)
library(data.table)
library(skimr)
library(zoo)

setwd("/home/michael/Dropbox/replications/air_pollution/correia-pollution/data/")

#-------------------------------------------------------------------------------
### Mortality data
# IHME: https://ghdx.healthdata.org/record/ihme-data/united-states-mortality-rates-county-1980-2014
#-------------------------------------------------------------------------------

df_mort <- list.files(path = "mortality/ihme/ihme_county/", full.names=T) %>%
  map_dfr(fread) %>%  # fread from data.table, should be faster
  tibble()
  
# note: this data is age-standardized
# drop national and states: FIPS<1000
# df_mort %>% 
#   filter(FIPS<1000) %>% 
#   distinct(location_name) %>% 
#   print(n=52)

df_mort <- df_mort %>% 
  rename(fips = FIPS, year = year_id, county=location_name, mortality = mx)

# split by sex
df_mort_m <- df_mort %>% 
  filter(fips>=1000 & sex=="Male" & cause_name %in% c('All causes', 'Cardiovascular diseases', 'Chronic respiratory diseases', 'Neoplasms', 'Injuries', 'Unintentional injuries'))
df_mort_f <- df_mort %>% 
  filter(fips>=1000 & sex=="Female" & cause_name %in% c('All causes', 'Cardiovascular diseases', 'Chronic respiratory diseases', 'Neoplasms', 'Injuries', 'Unintentional injuries'))

# combined sexes
df_mort <- df_mort %>% 
  filter(fips>=1000 & sex=="Both" & cause_name %in% c('All causes', 'Cardiovascular diseases', 'Chronic respiratory diseases', 'Neoplasms', 'Injuries', 'Unintentional injuries'))

df_mort_m <- df_mort_m %>% 
  pivot_wider(id_cols=c(fips, year), names_from=cause_id, values_from=mortality) %>% 
  rename(mort_all=`294`, mort_card=`491`, mort_resp=`508`, mort_canc=`410`, mort_inj=`687`, mort_acc=`696`) %>% 
  mutate(
    mort_cr = mort_card + mort_resp,
    mort_non_cr = mort_all - mort_cr) 

df_mort_f <- df_mort_f %>% 
  pivot_wider(id_cols=c(fips, year), names_from=cause_id, values_from=mortality) %>% 
  rename(mort_all=`294`, mort_card=`491`, mort_resp=`508`, mort_canc=`410`, mort_inj=`687`, mort_acc=`696`) %>% 
  mutate(
    mort_cr = mort_card + mort_resp,
    mort_non_cr = mort_all - mort_cr) 

df_mort <- df_mort %>% 
  pivot_wider(id_cols=c(fips, year), names_from=cause_id, values_from=mortality) %>% 
  rename(mort_all=`294`, mort_card=`491`, mort_resp=`508`, mort_canc=`410`, mort_inj=`687`, mort_acc=`696`) %>% 
  mutate(
    mort_cr = mort_card + mort_resp,
    mort_non_cr = mort_all - mort_cr) 

df_mort <- df_mort %>% select(fips, year, mort_all, mort_cr, mort_non_cr, mort_canc, mort_inj, mort_acc)
df_mort_m <- df_mort_m %>% select(fips, year, mort_all_m = mort_all, mort_cr_m = mort_cr,  mort_non_cr_m = mort_non_cr, mort_canc_m = mort_canc, mort_inj_m = mort_inj, mort_acc_m = mort_acc)
df_mort_f <- df_mort_f %>% select(fips, year, mort_all_f = mort_all, mort_cr_f = mort_cr, mort_non_cr_f = mort_non_cr,  mort_canc_f = mort_canc, mort_inj_f = mort_inj, mort_acc_f = mort_acc)

df_mort <- df_mort %>% 
  left_join(df_mort_m, by=c('year', 'fips')) %>% 
  left_join(df_mort_f, by=c('year', 'fips'))
  
#-------------------------------------------------------------------------------
### Population data
# https://www.nber.org/research/data/survey-epidemiology-and-end-results-seer-us-state-and-county-population-data-age-race-sex-hispanic
# https://data.nber.org/data/seer_u.s._county_population_data.html
  # county-level, adjusted, 1990-, 19 age groups, 4 races
    # 1969 data doesn't have hispanic
# total, share by ethnicity

#-------------------------------------------------------------------------------

df_pop_full <- read_csv("population/usrace19agesadj.csv")

# df_pop_full %>% distinct(hispanic) # 0/1
# df_pop_full %>% distinct(race) # 1=white, 2=black, 3=american indian, 4=asian/pacific islander
# note: hispanic and race are orthogonal (hispanic origin)

# rename fips=county, convert to numeric (drop leading 0)
df_pop_full <- df_pop_full %>% 
  rename(fips=county) %>% 
  mutate(fips = as.numeric(fips)) %>% 
  filter(fips<60000) # have state=KR in 2005; hurricane katrina displacement?

# aggregate across age, sex
df_pop_full <- df_pop_full %>% 
  group_by(year, st, stfips, fips, race, hispanic) %>% 
  summarize(pop = sum(pop)) %>% 
  arrange(year, fips)

# total county-year population
df_pop <- df_pop_full %>% 
  group_by(year, fips, st, stfips) %>% 
  summarize(pop_tot = sum(pop)) %>% 
  mutate(logpop = log(pop_tot))

# total hispanic population
  # rows are dropped entirely if count=0
df_pop_hisp <- df_pop_full %>% 
  filter(hispanic==1) %>% 
  group_by(year, fips) %>% 
  summarize(pop_hisp = sum(pop))

# total black population
df_pop_blk <- df_pop_full %>% 
  filter(race==2) %>% 
  group_by(year, fips) %>% 
  summarize(pop_blk = sum(pop))

# merge
df_pop <- df_pop %>% 
  left_join(df_pop_hisp, by=c('year', 'fips')) %>% 
  left_join(df_pop_blk, by=c('year', 'fips'))

# share hispanic, share black
  # have NAs
df_pop <- df_pop %>% 
  mutate(
    hisp = pop_hisp/pop_tot,
    hisp = replace_na(hisp,0),
    blk = pop_blk/pop_tot,
    blk = replace_na(blk,0))

#-------------------------------------------------------------------------------
### Pollution data
# https://www.caces.us/data

# EPA: annual
  # https://aqs.epa.gov/aqsweb/airdata/download_files.html
    # Annual summary data, Concentration by monitor, 1999-2014
  # https://aqs.epa.gov/aqsweb/airdata/FileFormats.html
#-------------------------------------------------------------------------------
df_pol <- read_csv('pollution/caces.csv') %>% 
  # convert fips to numeric (drop leading 0)
  rename(pm25=pred_wght) %>% 
  mutate(fips = as.numeric(fips)) %>% 
  select(fips,year,pm25) %>% 
  arrange(year, fips)

#-------------------------------------------------------------------------------
# EPA: state-county-monitor-instrument-standard-pollutant level data
  # have multiple measurement standards per pollutant
  # have multiple instruments per monitor

df_epa_raw  <- list.files(path = "pollution/epa/annual/", full.names=T) %>%
  map_dfr(fread, colClasses=list(character=c('State Code', 'County Code'))) %>% 
  tibble()
# some file has 'state code' as int

# create fips from concatenating state code, county code; strip leading zero
df_epa_raw <- df_epa_raw %>% 
  mutate(fips = as.numeric(paste0(`State Code`, `County Code`))) %>% 
  filter(is.na(fips)==0 & fips<60000) # keep only contiguous US

# df_epa_raw %>% select(`State Code`, `County Code`, fips, `State Name`, `Address`, `Local Site Name`) %>% filter(is.na(fips)) %>% View()
  # somehow have monitors in Canada; State Code == CC, but numeric County Code

# rename columns: 
df_epa_raw <- df_epa_raw %>% 
  rename(monitor = `Site Num`, 
         parameter=`Parameter Code`,
         standard=`Pollutant Standard`, 
         pm25=`Arithmetic Mean`, 
         year=`Year`)

# filter rows with PM25: parameter==88101
  # pick a measurement standard: `Pollutant Standard` == "PM25 Annual 2012"
    # or could average across standards
df_epa_raw <- df_epa_raw %>% 
  filter(parameter==88101 & standard =="PM25 Annual 2012")

# aggregate to county-year level:
df_epa <- df_epa_raw %>% 
  group_by(fips, year) %>% 
  summarize(pm25 = mean(pm25))

df_epa <- df_epa %>% 
  rename(pm25_epa = pm25)

#-------------------------------------------------------------------------------
### real per-capita county income
# income: https://apps.bea.gov/regional/downloadzip.cfm
  # Personal Income (State and Local)
  # CAINC1: Annual Personal Income by County
# inflation: https://data.worldbank.org/indicator/FP.CPI.TOTL.ZG
#-------------------------------------------------------------------------------

df_inc <- read_csv('controls/income/income.csv')
df_inc <- df_inc %>% 
  rename(fips=GeoFIPS, region=Region) %>% 
  mutate(fips = as.numeric(fips)) %>% 
  filter(LineCode==3) # per capita income

# df_inc %>% filter(fips>90000 ) %>%  select(region, GeoName)
# regions to use for heterogeneity
# 1 New England   
# 2 Mideast       
# 3 Great Lakes   
# 4 Plains        
# 5 Southeast     
# 6 Southwest     
# 7 Rocky Mountain
# 8 Far West  

df_inc <- df_inc %>% 
  filter(fips%%1000!=0 & fips<60000) %>% 
  select(fips, region, c(`1999`:`2014`))

df_inc <- df_inc %>% 
  pivot_longer(c(`1999`:`2014`), names_to='year', values_to='income')

#---------------
df_inf <- read_csv('controls/income/inflation.csv', skip=4) %>% 
  filter(`Country Code`=='USA') %>% 
  select(c(`2000`:`2014`)) # use 1999 as baseline

df_inf <- df_inf %>% 
  pivot_longer(c(`2000`:`2014`), names_to='year', values_to='inflation')

df_inf <-  df_inf %>% 
  mutate(compound=1, # prev year = 100
         inflation = 1 + inflation/100,
         compound = cumprod(inflation))

df_inf <- df_inf %>% add_row(year='1999', inflation=1, compound=1, .before=1)

#---------------

df_inc <- df_inc %>% 
  left_join(df_inf, by='year')

df_inc <- df_inc %>% 
  mutate(income = as.numeric(income)) %>% 
  rename(income_nom = income) %>% 
  mutate(income = income_nom/compound,
         year = as.numeric(year))

df_inc <- df_inc %>% 
  filter(is.na(income)==0) %>%  # have NAs in original data; eg Alaska
  mutate(loginc = log(income))

df_inc <- df_inc %>% select(fips, region, year, loginc)

#-------------------------------------------------------------------------------
### High school education
# proportion with high school or above
# https://www.ers.usda.gov/data-products/county-level-data-sets/download-data/
# excel file, save as csv
#-------------------------------------------------------------------------------
df_edu <- read_csv('controls/education/education.csv', skip=4) %>% 
  filter(!(State %in% c('US', 'AK', 'PR'))) %>% 
  select(fips=`Federal Information Processing Standards (FIPS) Code`, 
         #state = State, 
         edu1990= `Percent of adults with less than a high school diploma, 1990`, edu2000= `Percent of adults with less than a high school diploma, 2000`, edu2007 = `Percent of adults with less than a high school diploma, 2007-11`, edu2016 = `Percent of adults with less than a high school diploma, 2016-20`) %>% 
  mutate(fips = as.numeric(fips)) %>% 
  filter(fips%%1000!=0 & fips<60000)

# ACS data is 5-year average
df_edu <- df_edu %>% 
  mutate(
    edu2008 = edu2007,
    edu2009 = edu2007,
    edu2010 = edu2007,
    edu2011 = edu2007,
    )

# pivot longer
df_edu <- df_edu %>% 
  pivot_longer(cols=starts_with('edu'), names_to='year', names_prefix='edu', values_to = 'edu') %>% 
  arrange(fips, year) %>% 
  mutate(year = as.numeric(year))

# add rows
df_edu <- df_edu %>% 
  # group_by(fips) %>% 
  complete(year = full_seq(year,1), fips) %>% 
  arrange(fips, year)

# interpolate
df_edu <- df_edu %>% 
  group_by(fips) %>% 
  mutate(edu = na.approx(edu, rule=20)) # some counties have no data in 1990

# transform from 'less than high school' to 'high school or above: (100-x)/100, 
df_edu <- df_edu %>% 
  mutate(edu = (100-edu)/100) %>% 
  filter(between(year, 1999, 2014))

#-------------------------------------------------------------------------------
# urban rate
  # percent of population living in urban areas 
  # US Census 2010
# https://www.census.gov/programs-surveys/geography/guidance/geo-areas/urban-rural/2010-urban-rural.html

df_urb <- read_csv('controls/urban/urban.csv') %>% 
  filter(!(STATENAME %in% c('Alaska', 'Puerto Rico'))) %>% 
  mutate(fips = as.numeric(paste0(STATE, COUNTY))) %>% 
  select(fips, area=AREA_COU, urban=POPPCT_URBAN) %>% 
  mutate(area = area/1e6 * 0.386102, # convert from square meters to square miles
         urban = urban/100
         ) 

#-------------------------------------------------------------------------------
### final dataset

df <- df_mort %>% 
  left_join(df_pop, by=c('fips', 'year')) %>% 
  left_join(df_pol, by=c('fips', 'year')) %>% 
  left_join(df_epa, by=c('fips', 'year')) %>% 
  left_join(df_inc, by=c('fips', 'year')) %>% 
  left_join(df_edu, by=c('fips', 'year')) %>% 
  left_join(df_urb, by=c('fips')) %>% 
  filter(between(year, 1999, 2014)) %>% 
  mutate(popden = pop_tot / area)

df %>% 
  filter(is.na(pm25)==0) %>% 
  write_csv('data.csv')
#-------------------------------------------------------------------------------
