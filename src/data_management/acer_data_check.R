# This file generates the statistics for Online Appendix A.1
# with statistics about the evolution of wholesale prices, retail prices, and markups
# in several European countries.
rm(list=ls())
cat("\014")
# Load project path definitions.
source("project_paths.R")
library(hrbrthemes)
library(tidyverse)
library("RColorBrewer")
library(ggplot2)
library(lubridate)
library(readxl)
library(xtable)
# This only has to be run once on each system, afterwards, restart RStudio
# Please see the separate file that should be run manually on a Mac.
#extrafont::font_import()
extrafont::loadfonts()

options(xtable.floating = FALSE)
options(xtable.timestamp = "")
acer_mu <- read_excel(str_c(PATH_IN_DATA,"/","acer_markups.xlsx"), skip = 11, col_names= c("year","country","pr","pw","mar")) %>% 
    filter(year>=2012) %>%
    #filter(year<=2016) %>%
    mutate(muw = 100 * mar / pw ) %>%
    mutate(mur = 100 * mar / pr) %>%
    mutate(country=ifelse(country=="The Netherlands","Netherlands",country)) %>%
    mutate(country=ifelse(country=="Great Britain","United Kingdom",country))
# Subset of countries with deregulated markets.
countries_comp <- c("Belgium","Germany", "United Kingdom","Italy","Luxembourg","Netherlands","Sweden","Norway")
# Subset of countries with regulated markets.
countries_reg <- c("Denmark","Spain","France", "Poland","Portugal")
# Select relevant subsets of countries.
acer_mu_comp <- filter(acer_mu,country %in% countries_comp)
acer_mu_reg <- filter(acer_mu,country %in% countries_reg)
################################################################################
# Look only at Belgium (this is not identical to Flanders, but should be reasonably comaprable)
mu_bel <- filter(acer_mu,country=="Belgium") %>%
    arrange(year)
mu_bel_mat <- matrix(c(mu_bel$pr,mu_bel$pw,mu_bel$muw),nrow=7)
rownames(mu_bel_mat) <- c("2012","2013","2014","2015","2016","2017","2018")
colnames(mu_bel_mat) <- c("Retail Price","Wholesale Price","Markup")
# Print table with markups based on ACER.
mu_bel_table <- xtable(mu_bel_mat)
align(mu_bel_table) <- "r|lll"
print(mu_bel_table, file=str_c(PATH_OUT_TABLES,"/","mu_bel_table.tex"))

################################################################################
# Set fixed color scheme to that schemes are comparable across different graphs.
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc_muw <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 10))

# Create graph with markups relative to wholesale price.
ggplot(data = acer_mu_comp, aes(country, year, fill= muw)) + 
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Markups over time (deregulated prices)") +
    labs(x = "", y="") +
    guides(fill = guide_colourbar(title = "Markup in %"))
ggsave(str_c(PATH_OUT_FIGURES, "/", "markups_time_deregulated.png"))

ggplot(data = acer_mu_reg, aes(country, year, fill= muw)) + 
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    theme_ipsum() +
    labs(x = "", y="") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Markups over time (regulated prices)") +
    guides(fill = guide_colourbar(title = "Markup (in %)"))
ggsave(str_c(PATH_OUT_FIGURES, "/", "markups_time_regulated.png"))

################################################################################
# Create graph with retail prices.
ggplot(data = acer_mu_comp, aes(country, year, fill= pr)) + 
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Retail prices over time (deregulated prices)") +
    labs(x = "", y="") +
    guides(fill = guide_colourbar(title = "Price (EUR per mwh)"))
ggsave(str_c(PATH_OUT_FIGURES, "/", "retailprices_time_deregulated.png"))

ggplot(data = acer_mu_reg, aes(country, year, fill= pr)) + 
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Retail prices over time (regulated prices)") +
    labs(x = "", y="") +
    guides(fill = guide_colourbar(title = "Price (EUR per mwh)"))
ggsave(str_c(PATH_OUT_FIGURES, "/", "retailprices_time_regulated.png"))

################################################################################
# Create graph with wholesale prices.
ggplot(data = acer_mu_comp, aes(country, year, fill= pw)) + 
    geom_tile() +
    #scale_fill_distiller(palette = "RdPu") +
    scale_fill_gradient(low="white", high="blue") +
    theme_ipsum() +
    labs(x = "", y="") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Wholesale prices over time (deregulated prices)") +
    guides(fill = guide_colourbar(title = "Price (EUR per mwh)"))
ggsave(str_c(PATH_OUT_FIGURES, "/", "wholesaleprices_time_deregulated.png"))

ggplot(data = acer_mu_reg, aes(country, year, fill= pw)) + 
    geom_tile() +
    #scale_fill_distiller(palette = "RdPu") +
    scale_fill_gradient(low="white", high="blue") +
    theme_ipsum() +
    labs(x = "", y="") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Wholesale prices over time (regulated prices)") +
    guides(fill = guide_colourbar(title = "Price (EUR per mwh)"))
ggsave(str_c(PATH_OUT_FIGURES, "/", "wholesaleprices_time_regulated.png"))


################################################################################
# Illustration of role of other charges and taxes.
pc_all <- read_excel(str_c(PATH_IN_DATA,"/","acer_price_breakdown.xlsx"), skip = 11, 
                         col_names= c("year","country_code","country","energy","network","res","tac","vat","potb")) %>%
    select(-country_code) %>%
    arrange(year) %>%
    mutate(energy = energy *100, network = network *100, res = res *100, tac = tac * 100, vat = vat * 100)
# Split data by countries with competition and regulation
pc_comp <- pc_all %>%
    filter(country %in% countries_comp)
pc_reg <- pc_all %>%
    filter(country %in% countries_reg)
pc_bel <- filter(pc_all,country=="Belgium") %>%
    select(-country)
# Create overview of average costs for consumers in different countries.
comp_levels_small <- c(Energy = "energ_cost", Charges = "charges")

pc_total_comp <- group_by(pc_comp,country) %>%
    summarize(cost = mean(potb),energy=mean(energy)) %>%
    mutate(country = factor(country)) %>%
    mutate(country = fct_reorder(country,cost)) %>%
    arrange(desc(cost)) %>%
    mutate(energ_cost = cost * energy/100) %>%
    mutate(charges = cost - energ_cost) %>%
    gather(charges,energ_cost,key=component,value=eur_cost) %>%
    mutate(component = factor(component))
# Recode component factor variable.
pc_total_comp$component <-    fct_recode(pc_total_comp$component, !!!comp_levels_small)
# Repeat procedure for regulated markets.
pc_total_reg <- group_by(pc_reg,country) %>%
    summarize(cost = mean(potb),energy=mean(energy)) %>%
    mutate(country = factor(country)) %>%
    mutate(country = fct_reorder(country,cost)) %>%
    arrange(desc(cost)) %>%
    mutate(energ_cost = cost * energy/100) %>%
    mutate(charges = cost - energ_cost) %>%
    gather(charges,energ_cost,key=component,value=eur_cost) %>%
    mutate(component = factor(component))
# Recode component factor variable.
pc_total_reg$component <-    fct_recode(pc_total_reg$component, !!!comp_levels_small)

# Overview of total retail electricity costs for consumers.
ggplot(data = pc_total_comp, mapping = aes(x = country,y=eur_cost,fill=component)) + 
    geom_bar(stat="identity", width=0.65) +
    ggtitle('Average total yearly electricity costs - Dereg. prices') +
    theme_ipsum() +
    scale_fill_brewer() +
    geom_text(aes(label = round(eur_cost)), size = 2.5, hjust = 1.5, vjust = 0, position ="stack") +
    labs(x = "", y="Total yearly cost (in EUR)",fill="Components") +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
    coord_flip()
ggsave(str_c(PATH_OUT_FIGURES, "/", "price_composition_deregulated.png"))

# Overview of total retail electricity costs for consumers.
ggplot(data = pc_total_reg, mapping = aes(x = country,y=eur_cost,fill=component)) + 
    geom_bar(stat="identity", width=0.65) +
    ggtitle('Average total yearly electricity costs - Reg. prices') +
    theme_ipsum() +
    scale_fill_brewer() +
    geom_text(aes(label = round(eur_cost)), size = 2.5, hjust = 1.5, vjust = 0, position ="stack") +
    labs(x = "", y="Total yearly cost (in EUR)",fill="Components") +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
    coord_flip()
ggsave(str_c(PATH_OUT_FIGURES, "/", "price_composition_regulated.png"))

# Create table for Belgium.
pc_bel_mat <- matrix(c(pc_bel$energy,pc_bel$network,pc_bel$res,pc_bel$tac,pc_bel$vat,pc_bel$potb),nrow=4)
rownames(pc_bel_mat) <- c("2015","2016","2017","2018")
comp_names <- c("Energy","Network","RES","TandC","VAT","POTB")
colnames(pc_bel_mat) <- comp_names
comp_levels <- c(Energy = "energy", Network = "network", RES = "res", TaC = "tac",VAT = "vat")
# Print table with markups based on ACER.
pc_bel_table <- xtable(pc_bel_mat,digits=c(0,0,0,0,0,0,0))
align(pc_bel_table) <- "r|lllll|l"
pc_bel_table
print(pc_bel_table, file=str_c(PATH_OUT_TABLES,"/","pc_bel_table.tex"))


# Create plot for regulated markets.
pc_reg_plot <- pc_reg %>%
    gather(energy,network,res,tac,vat, key = "component", value = "share") %>%
    mutate(component = factor(component)) %>%
    group_by(country,component) %>%
    mutate(factor(component)) %>%
    summarize(share=mean(share),potb=mean(potb))
# Recode component factor variable.
pc_reg_plot$component <-    fct_recode(pc_reg_plot$component, !!!comp_levels)
ggplot(pc_reg_plot, aes(x = country, y = share , fill = fct_rev(component))) + 
    geom_bar(stat = "identity") +
    theme_ipsum() +
    ggtitle("Total price composition (regulated prices)") +
    scale_fill_brewer() +
    geom_text(aes(label = round(share)), size = 2, hjust = 0.5, vjust = 3, position ="stack") +
    labs(x = "", y="Share of component (in %)",fill="Components") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
ggsave(str_c(PATH_OUT_FIGURES, "/", "price_breakdown_deregulated.png"))

# Create plot for deregulated markets.
pc_comp_plot <- pc_comp %>%
    gather(energy,network,res,tac,vat, key = "component", value = "share") %>%
    mutate(component = factor(component)) %>%
    group_by(country,component) %>%
    mutate(factor(component)) %>%
    summarize(share=mean(share),potb=mean(potb))
# Recode component factor variable.
pc_comp_plot$component <-    fct_recode(pc_comp_plot$component, !!!comp_levels)
ggplot(pc_comp_plot, aes(x = country, y = share , fill = fct_rev(component))) + 
    geom_bar(stat = "identity") +
    theme_ipsum() +
    ggtitle("Total price composition (deregulated prices)") +
    scale_fill_brewer() +
    geom_text(aes(label = round(share)), size = 2, hjust = 0.5, vjust = 3, position ="stack") +
    labs(x = "", y="Share of component (in %)",fill="Components") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
ggsave(str_c(PATH_OUT_FIGURES, "/", "price_breakdown_regulated.png"))