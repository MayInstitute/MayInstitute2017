
require(testthat)

context("ProtExp")

test_that("sanity check", {

	expect_equal(1 + 1, 2)

})

test_that("loading data works", {

	data(twin)
  
  expect_is(twin_dia, "data.frame")
  
  expect_is(twin_srm, "data.frame")

	expect_equal(twin_dia$protein[1], "A1AG_BOVINE")

})

test_that("ProtExp methods works", {
  
  library(tidyverse)
  
  data(twin)
  
  twin_dia2 <- twin_dia %>%
    rename(heavy = intensity_h, light = intensity_l) %>% 
    gather(label, intensity, heavy, light)
  
  data1 <- ProtExp(protein=twin_dia2$protein,
                  feature=twin_dia2$feature,
                  run=twin_dia2$run,
                  intensity=twin_dia2$intensity,
                  label=twin_dia2$label)

  expect_is(data1, "ProtExp")
    
  expect_is(data1, "data.frame")

  data1norm <- normalize(data1, by="heavy")
  
  data1summ <- summary(data1norm)
  
  expect_equal(data1summ["R001","A1AG_BOVINE","light"],
               18.29013, tolerance=0.01)
    
})

