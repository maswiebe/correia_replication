# correia_replication
Replicating [Correia et al. (2013)](https://journals.lww.com/epidem/Fulltext/2013/01000/Effect_of_Air_Pollution_Control_on_Life_Expectancy.4.aspx) on the correlation between PM2.5 and mortality, using US county-level data and `R`.

Clean the raw data using `correia_clean.r`, and run regressions with `correia_results.r`.

The cleaned dataset is `data.csv`.

To re-run the results using restricted NCHS mortality data, use `correia_results_nchs.r`.

### Data sources
- [Mortality](https://ghdx.healthdata.org/record/ihme-data/united-states-mortality-rates-county-1980-2014): download and extract the five CSV files for 'Mortality Rates by County 1980-2014'.
- [Population](https://www.nber.org/research/data/survey-epidemiology-and-end-results-seer-us-state-and-county-population-data-age-race-sex-hispanic): download the CSV file for All States Combined (Adjusted), County, 19 Age Groups, 1990-, 4 Expanded Races by Origin (usrace19agesadj.csv).
- [Modeled PM2.5](https://www.caces.us/data): enter your email and organization, select LUR, National, County, PM2.5, 1999-2014. Data is from [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0228535).
- [Raw PM2.5 (EPA)](https://aqs.epa.gov/aqsweb/airdata/download_files.html#Annual): Annual summary data, Concentration by monitor, download and extract files for 1999-2014.
- [Income](https://apps.bea.gov/regional/downloadzip.cfm): under 'Personal Income (State and Local)', download 'CAINC1: Annual Personal Income by County'. Adjust for inflation using US data in the CSV [here](https://data.worldbank.org/indicator/FP.CPI.TOTL.ZG).
- [Education](https://www.ers.usda.gov/data-products/county-level-data-sets/download-data/): download the excel file for 'Educational attainment', save as CSV.
- [Urban rate](https://www.census.gov/programs-surveys/geography/guidance/geo-areas/urban-rural/2010-urban-rural.html): download the excel file for 'Percent Urban and Rural In 2010 by State and County', save as CSV.
