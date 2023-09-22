# This file computes several measures of a wholesale price index. It is only used as an
# input for the counterfactuals that compute indifference markups, as well as some statistics
# in Online Appendix B.1.

rm(list=ls())
cat("\014")
# Load project path definitions.
source("project_paths.R")
library(tidyverse)
library(ggplot2)
library(lubridate)
# On some Mac machines this code snippet may have to be run first.
# If you encounter an error, please run this code manually.
# library(remotes)
# remotes::install_version("Rttf2pt1", version = "1.3.8")
# library(extrafont)
# extrafont::font_import()
# loadfonts(device="pdf")
# # On Windows systems.
# # loadfonts(device="win")
# fonts()
# extrafont::font_import()
extrafont::loadfonts()
        
wholesale_spot <- read_delim(str_c(PATH_IN_DATA, "/", "wholesale_spot.csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                             grouping_mark = "."), trim_ws = TRUE) %>%
                    filter(YearMonth>=201201 & YearMonth<201607) %>%
                    rename(spot_price = BPXBBPXPRC) %>%
                    mutate(YearMonth = parse_date(as.character(YearMonth),format='%Y%m')) %>%
                    mutate(part_day = cut(quarter_num,breaks=c(1,29,81,100),labels=c('Early','Peak','Late')))
# Define split of day.
# Early: quarter_num 1-30 (midnight -7:00am)
# During day: quarter_num 30-80 (7:00am-8:00pm)
# Late evening: quarter_num > 80 (8:00pm to midnight)
summary(wholesale_spot)

# Group data by month to be comparable to what we actually need in the model counterfactual for markups.
wholesale_month <- group_by(wholesale_spot, YearMonth, part_day) %>%
    filter(is.na(part_day)==FALSE) %>%
    summarize(mean_spot = mean(spot_price,na.rm=TRUE))

# Create graph similar to ACER on year level..
wholesale_avg_by_day <- group_by(wholesale_spot,year,part_day) %>%
    filter(is.na(part_day)==FALSE) %>%
    summarize(mean_spot = mean(spot_price,na.rm=TRUE))
wholesale_avg_by_partday <- group_by(wholesale_avg_by_day,part_day) %>%
    summarize(mean_spot = mean(mean_spot,na.rm=TRUE))
wholesale_avg_by_partday

ggplot(data=wholesale_avg_by_day,mapping=aes(x=year,y=mean_spot,color=part_day)) +
    geom_line() +
    ggtitle('Average wholesale price over years by part of day')
ggsave(str_c(PATH_OUT_FIGURES, "/", "wspot_partofday.pdf"))

# Plot average retail price (over days) as function of time of day.
wholesale_intraday <- group_by(wholesale_spot,quarter_num) %>%
    summarize(spot_price = mean(spot_price,na.rm=TRUE))
ggplot(data=wholesale_intraday,mapping=aes(x=quarter_num,y=spot_price)) +
    geom_line() +
    ggtitle("Evolution of wholesale price during average day")
ggsave(str_c(PATH_OUT_FIGURES, "/", "wspot_intraday.pdf"))
# Plot different types of wholesale prices
ggplot(data=wholesale_month,mapping=aes(x=YearMonth,y=mean_spot,color=part_day)) +
    geom_line() +
    ggtitle('Evolution of monthly wholesale prices by part of day')
ggsave(str_c(PATH_OUT_FIGURES, "/", "wspot_monthly_part_of_day.pdf"))

# Prepare wholesale price time series for use in graphs and counterfactuals.
# For sample period January 2012 to June 2016.
wholesale_spot_export <- spread(wholesale_month,key=part_day,value=mean_spot)
write_csv(wholesale_spot_export,str_c(PATH_OUT_DATA,"/", "wholesale_corr_export.csv"))
