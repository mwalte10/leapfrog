test_that("DemProj only matches EPP-ASM", {

  ## Check that population age 15:79 matches between
  ## Note: the open 80+ population does not match because EPP-ASM did
  ##   not handle survivorship of the open age group correctly. This
  ##   is corrected in leapfrog.

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
  
  expect_warning(fp <- eppasm::prepare_directincid(pjnz1),
                 "no non-missing arguments to min; returning Inf")
  fp$tARTstart <- 61L
  mod <- eppasm::simmod(fp)

  demp <- prepare_leapfrog_demp(pjnz1)
  hivp <- prepare_leapfrog_projp(pjnz1)
  lmod <- leapfrogR(demp, hivp)

  expect_equal(lmod$totpop1[16:80,,], mod[1:65,,1,])
})

test_that("Leapfrog matches DemProj projection without migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff <- lmod1$totpop[,,2] - demp1$basepop[,,2]

  specres <- eppasm::read_hivproj_output(pjnz1)

  ## For age 0, there is a difference with Spectrum; unsure why
  diff[1,]

  expect_true(all(abs(diff[-1,]) < 0.001))

  ## deaths by sex/age
  expect_true(all(abs(lmod1$natdeaths[,,-1] - specres$natdeaths[,,-1]) < 0.003))

  ## births by age
  ## !!! QUITE LENIENT THRESHOLD
  ## TO DO: review births calculation wiht Avenir
  expect_true(all(abs(lmod1$births[-1] - specres$births[-1]) < 0.2))
  
})

test_that("Leapfrog matches DemProj projection with migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff <- lmod1$totpop1[,,2] - demp1$basepop[,,2]

  ## With migration, difference goes through age 5; suspect some Beers
  ## step in Spectrum that is missing?
  diff[1:10,]

  ## Age 80+ group also slightly larger diff now
  diff[81,]

  expect_true(all(abs(diff[-(1:6),]) < 0.005))
})
