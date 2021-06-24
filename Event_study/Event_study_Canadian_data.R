# ==================================================
# Event_study_Canadian_data.R
# 
# Download the Canadian data for the event study.
# ==================================================

# == Preliminaries ========================
rm(list = ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Load useful packages
library(cansim)
library(fredr)
library(dplyr)
library(zoo)
library(quantmod)
library(zoo)
library(mFilter)
library(reshape2)
require(broom)
require(knitr)

# === Set API key ===
# NB: Must install your own Fred API key to run the code
fredr_set_key("079901da86c80566465535381740061c")


# === Download US 2-year yield data ===
# Get data from FRED
seriescode <- c("DGS2")
seriesname <- c("US_twoyear")

t.start <- '2009-01-01'
t.end   <- '2015-12-31'

for (i_ in 1 : length(seriescode))
{
  temp <-  fredr(series_id = seriescode[i_], 
                 observation_start = as.Date(t.start),
                 observation_end = as.Date(t.end),
                 frequency = "d") %>% filter(!is.na(value)) %>%
    xts(x = .$value, order.by = .$date)
  names(temp) <- seriesname[i_]
  temp$dd <- (as.numeric((temp[,c(1)]))) - (as.numeric(stats::lag(temp[,c(1)])))
  temp$seriesname[i_] <- as.numeric(temp$seriesname[i_])
  assign(seriesname[i_],temp)
  
}

rm(temp)

# === Download exchange rate data ===
# Get data from FRED
seriescode <- c("DEXCAUS")
seriesname <- c("Exchange_Rate")

# Note exchange rate is already in format up = depreciation of CAD
t.start <- '2009-01-01'
t.end   <- '2015-12-31'

# Download exchange rate and format
for (i_ in 1 : length(seriescode))
{
  temp <-  fredr(series_id = seriescode[i_], 
                 observation_start = as.Date(t.start),
                 observation_end = as.Date(t.end),
                 frequency = "d") %>% filter(!is.na(value)) %>%
    xts(x = .$value, order.by = .$date)
  names(temp) <- seriesname[i_]
  temp$dd <- 100*(log(dplyr::lead(as.numeric((temp[,c(1)])))) - log((as.numeric((temp[,c(1)]))))) 
  temp$seriesname[i_] <- as.numeric(temp$seriesname[i_])
  assign(seriesname[i_],temp)
  
}

rm(temp)

# === Download Canadian interest rate data ===
# Vectors of the Canadian interest rate series that I use
# 3, month, 6 month, 1 year, 2, year, 5 year, 10 year bond yields
cirvec <- c("v39065", "v39067", "v39051", "v39053", "v39055")
cirnam <- c("threemonth", "oneyear", "twoyear", "fiveyear", "tenyear") 

# Loop over the interest rate series and calculate daily change
for (i_ in 1 : length(cirvec))
{
  temp <- get_cansim_vector(cirvec[i_], start_time ="2009-01-01", end_time = "2015-12-31")
  temp <- temp[,c(2,3)]
  names(temp) <- c(cirnam[i_],"Date")
  temp$Date <- as.Date(temp$Date)
  temp <- xts(temp, order.by = temp$Date)
  temp <- temp[,c(1)]
  temp$dd <- as.numeric(temp[,c(1)]) - as.numeric(stats::lag(temp[,c(1)]))
  temp$cirnam[i_] <- as.numeric(temp$cirnam[i_])
  assign(cirnam[i_],temp)
  
}

# Tidy up the workspace
rm(temp)

# ==============================================================================
# Save data series
save(Exchange_Rate, US_twoyear, threemonth, oneyear, twoyear, fiveyear, tenyear, file = "Event_study_data.RData")
