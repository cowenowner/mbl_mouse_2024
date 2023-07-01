# CowenLib.R
# Core R functions used by the cowen group.
# To use this, just use source('.../CowenLib.R')
#
# Cowen 2021

# Include libraries that are required for these functions.

# Fisher Z transform.
FishZ <- function(r){
  result <- 0.5*(log(1+r)-log(1-r))
  return(result)
}

# PLotting functions ------
PlotItSVG <- function(p, fname, width = 4.5, height = 4.5){
  svg(fname, width = width , height = height) # This can be exported into powerpoint.
  print(p)
  dev.off()
  p
}


library(officer)
library(rvg)


PlotItPPT <- function(plt, fname){
  ## NOTE: the height and width parameters in ph_with/ph_locaion do not work at all so don't bother.
  read_pptx() %>% add_slide(layout = "Two Content", master = "Office Theme") %>%
    ph_with(dml(ggobj = plt), location = ph_location_type()) %>% print(target = fname)
  print(plt)
  return(plt)
}
