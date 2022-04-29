library(here)
library(rmarkdown)

render("PuMBA_Demonstration.Rmd", output_format = "all")

file.rename("PuMBA_Demonstration.md", "README.md")

setwd(here::here())
