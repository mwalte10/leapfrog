test_that("DemProj only matches EPP-ASM", {

  ## Check that population age 15:79 matches between
  ## Note: the open 80+ population does not match because EPP-ASM did
  ##   not handle survivorship of the open age group correctly. This
  ##   is corrected in leapfrog.

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
  
  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)

  ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
  ## in EPP-ASM preparation
  demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
  demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
  
  lmod <- leapfrogR(demp, hivp)

  expect_warning(fp <- eppasm::prepare_directincid(pjnz1),
                 "no non-missing arguments to min; returning Inf")
  fp$tARTstart <- 61L

  ## Replace ASFR because demp$asfr is normalised, but fp$asfr is not
  fp$asfr <- demp$asfr
    
  mod <- eppasm::simmod(fp)

  expect_equal(lmod$totpop1[16:80,,], mod[1:65,,1,])
})

test_that("Leapfrog matches DemProj projection without migration", {

  demo_matches_birthsdeaths("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ")
  demo_matches_totpop("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ")
  
})

test_that("Leapfrog matches direct incidence option in EPP-ASM, no ART and no hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-no-hiv-deaths_spectrum-v6.13_2022-02-12.pjnz"
  demo_matches_birthsdeaths(pjnz1)
  demo_matches_totpop(pjnz1)
  trans_matches(pjnz1)
})

test_that("Leapfrog matches direct incidence option in EPP-ASM, no ART + hiv mort", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-no-art-plus-hiv-deaths_spectrum-v6.13_2022-02-12.PJNZ"
  demo_matches_birthsdeaths(pjnz1)
  demo_matches_totpop(pjnz1)
  trans_matches(pjnz1)
})


test_that("Leapfrog matches direct incidence option in EPP-ASM, ART", {
  ## Check that prevalence, deaths and incidence  matches between
  ## the two models
  pjnz1 <- "../testdata/spectrum/v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ"
  demo_matches_birthsdeaths(pjnz1)
  demo_matches_totpop(pjnz1)
  trans_matches(pjnz1)
})
