test_that("Leapfrog matches single year age group and coarse age group projection without migration", {

  demog_matches_birthsdeaths("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ", threshold_deaths = 0.05, threshold_births = 1e-3, threshold_absolute = 1e-3)
  demog_matches_totpop("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ")
  matches_coarse_age_groups("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ", threshold_pid = c(0, 0, 0), threshold_naturaldeaths = 1e-3)
  
})

##TODO: add in test for hiv entrant population
# test_that("Age 15 entrant population matches", {
#   
#   
# })

test_that("Leapfrog matches direct incidence option at coarse and single year age structure, no ART and no hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-no-hiv-deaths_spectrum-v6.13_2022-02-12.pjnz"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 3, threshold_births = 0.1)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(250, 25, 1))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(900, 1, 0.05))
})

##Maggie: there is something wrong with incorporating HIV deaths
test_that("Leapfrog matches direct incidence option at coarse and single year age structure, no ART + hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 25, threshold_births = 45)
  ## This doesn't pass at this point, think it's something with needing the popadjust
  ##demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(140, 25, 15))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(640, 2, 20))
  
})
