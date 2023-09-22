# If you encounter problems when running some of the R code (which could happen on some Mac systems),
# please run this file before trying to nicely format graphs in R.
# It fixes some issues with saving particular types of graphs as pdf and png.
# Please note that this program might downgrade your Rttf2pt1 package.
library(remotes)
remotes::install_version("Rttf2pt1", version = "1.3.8")
library(extrafont)
loadfonts(device = "pdf")
font_import()
fonts()
