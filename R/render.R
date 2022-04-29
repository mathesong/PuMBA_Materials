library(here)
library(rmarkdown)

setwd(here::here("R"))

render("PuMBA_Demonstration.Rmd", output_format = "all")

file.rename("PuMBA_Demonstration.md", "README.md")

setwd(here::here())
