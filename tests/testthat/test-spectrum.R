test_that("Leapfrog matches single year age group and coarse age group projection without migration", {

  demog_matches_birthsdeaths("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ", threshold_deaths = 0.05, threshold_births = 1e-3)
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
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(230, 25, 1))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(705, 0.01, 0.05))
})

##Maggie: there is something wrong with incorporating HIV deaths
test_that("Leapfrog matches direct incidence option at coarse and single year age structure, no ART + hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(145, 25, 15))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(640, 4, 20), threshold_naturaldeaths = 3)

})

test_that("Leapfrog matches direct incidence option in EPP-ASM, no ART + hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-hiv-fert_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 0.2, threshold_births = 0.13, threshold_absolute = 1e-3)
  ## demog_matches_totpop(pjnz1)
  trans_matches(pjnz1, threshold_absolute_pid = c(355, 50, 30))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(700, 50, 35))
  
})

test_that("Input child transmission matches", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(230, 25, 1))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(705, 0.01, 0.05))
})


test_that("Distribution of child infections matches", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-child-input-distribute-infections_spectrum-v6.13_2022-02-12.PJNZ"
  demog_matches_birthsdeaths(pjnz1, threshold_deaths = 1e-3, threshold_births = 0.01)
  demog_matches_totpop(pjnz1)
  transmission_matches(pjnz1, threshold_absolute_pid = c(230, 25, 1))
  matches_coarse_age_groups(pjnz1, threshold_pid = c(705, 0.01, 0.05))
})















