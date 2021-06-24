# Event_study.R
#
# Runs an event study of effect of US FG shocks on Canadian
# interest and exchange rates
# ==============================================================

rm(list = ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Install useful packages
library(openxlsx)
library(reshape2)
library(zoo)
library(mFilter)
library(quantmod)
library(tidyverse)
library(dplyr)
require(broom)
require(knitr)
library(pixiedust)
library(lmtest)
library(sandwich)

# === Load FG shocks ======================================================
fgshock_file <- "FG_Shocks.xlsx"

fg.shocks <- read.xlsx(xlsxFile = fgshock_file,
                     sheet = "Sheet1", startRow = 1, colNames = TRUE )

# Multiply by -1 to make plus an expansionary FG shock
fg.shocks$FG.shock = fg.shocks$FG.shock*-1
fg.shocks$Date <- as.Date(fg.shocks$Date, origin = "1899-12-30", tz = "UTC")

# FG shocks time series
fg.shocks.ts = fg.shocks[,2] 
fg.shocks.ts = xts(fg.shocks.ts, order.by=fg.shocks$Date)
# =========================================================================

# === Load interest and exchange rate data ================================

load("Event_study_data.RData")

# =========================================================================

# === Estimate event study regression for Canadian interest rates ===

cirnam <- c("US_twoyear", "threemonth", "oneyear", "twoyear", 
            "fiveyear", "tenyear", "Exchange_Rate") 


for (i_ in 1 : length(cirnam))
{
  temp <- merge(fg.shocks.ts, eval(parse(text=paste(cirnam[i_],"$dd",sep=""))),join='left')
  
  # Run regression and calculate robust standard errors
  tempmod <- lm(dd ~ fg.shocks.ts, data = temp)
  robust_t <- coeftest(tempmod,vcov=vcovHC(tempmod,type = "HC0"))
  
  assign( paste(cirnam[i_],".estimates",sep=""), robust_t )
  
  show(paste("Estimates for: ",cirnam[i_]))
  show(robust_t)
  
  # Tidy up the workspace
  rm(temp, temp_plot, temp_melt, temp_df, tempmod, robust_t)
  
}



