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

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only-no-mig_spectrum-v6.13_2022-02-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff <- lmod1$totpop[,,2] - demp1$basepop[,,2]

  specres <- eppasm::read_hivproj_output(pjnz1)

  expect_true(all(abs(diff < 0.001)))

  ## deaths by sex/age
  expect_true(all(abs(lmod1$natdeaths[,,-1] - specres$natdeaths[,,-1]) < 0.003))

  ## births by age
  expect_true(all(abs(lmod1$births[-1] - specres$births[-1]) < 0.002))
  
})

test_that("Leapfrog matches DemProj projection with migration", {

  pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
  demp1 <- prepare_leapfrog_demp(pjnz1)
  hivp1 <- prepare_leapfrog_projp(pjnz1)
  lmod1 <- leapfrogR(demp1, hivp1)

  diff <- lmod1$totpop1[,,2:6] - demp1$basepop[,,2:6]

  expect_true(all(abs(diff) < 0.01))
})

# test_that("Leapfrog matches direct incidence option in EPP-ASM, no ART", {
#   ## Check that prevalence, deaths and incidence  matches between
#   ## the two models
#   pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_aim-adult-no-art_spectrum-v6.13_2022-02-12.pjnz")
# 
#   demp <- prepare_leapfrog_demp(pjnz1)
#   hivp <- prepare_leapfrog_projp(pjnz1)
# 
#   ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
#   ## in EPP-ASM preparation
#   demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
#   demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
# 
#   lmod <- leapfrogR(demp, hivp)
# 
#   expect_warning(fp <- eppasm::prepare_directincid(pjnz1),
#                  "no non-missing arguments to min; returning Inf")
#   fp$tARTstart <- 61L
# 
#   ## Replace ASFR because demp$asfr is normalised, but fp$asfr is not
#   fp$asfr <- demp$asfr
# 
#   mod <- eppasm::simmod(fp)
#   
#   ##EPPASM and leapfrog are done at different age specifications
#   age_mapping <- data.frame(eppasm_age = c(rep(1,2), rep(2,2), unlist(lapply(c(3:8), rep, times = 5)), rep(9, 32)), leapfrog_age = 1:66)
# 
#   ## PREVALENCE
#   mod_prev <- attr(mod, 'hivpop')
#   dimnames(mod_prev) <- list('CD4' = c(1:7), 'eppasm_age' = c(1:9), 'sex' = c(1:2), 'year' = c(1:61))
#   mod_prev <- mod_prev %>%
#     cubelyr::as.tbl_cube(met_name = "EPPASM") %>%
#     as_tibble %>%
#     group_by(year, eppasm_age, CD4) %>% 
#     summarise(across(`EPPASM`, sum))  
#   
#   lmod_prev <- lmod$hivstrat_adult
#   dimnames(lmod_prev) <- list('CD4' = c(1:7), 'leapfrog_age' = c(1:66), 'sex' = c(1:2), 'year' = c(1:61))
#   lmod_prev <- lmod_prev %>%
#     cubelyr::as.tbl_cube(met_name = "leapfrog") %>%
#     as_tibble
#   lmod_prev <- left_join(lmod_prev, age_mapping, by = 'leapfrog_age')
#   lmod_prev <- lmod_prev %>%
#     group_by(year, eppasm_age, CD4) %>% 
#     summarise(across(`leapfrog`, sum))  
#   
#   prev_comparison <- full_join(lmod_prev, mod_prev)  %>% 
#     group_by(year, eppasm_age, CD4) %>% 
#     mutate(pct_diff = abs(EPPASM - leapfrog)  / leapfrog, leapfrog = leapfrog, EPPASM = EPPASM) %>%
#     mutate(pct_diff = ifelse(is.nan(pct_diff), 0, pct_diff), leapfrog = leapfrog, EPPASM = EPPASM)
# 
#   ##percent difference as less than 1%, this is failing rn bc
#   ##the transmission model hasn't been fully implemented
#   ##write.csv(prev_comparison, "./tests/testdata/plotting_diagnostics/prev_comparison.csv", row.names = F)
#   expect_true(all(abs(prev_comparison$pct_diff) < 0.01))
# 
#   ##INCIDENCE
#   mod_inc <- attr(mod, 'infections')
#   mod_inc <- mod_inc[-66,,]
#   dimnames(mod_inc) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   mod_inc <- mod_inc %>%
#     cubelyr::as.tbl_cube(met_name = "EPPASM") %>%
#     as_tibble 
#   
#   lmod_inc <- lmod$infections
#   lmod_inc <- lmod_inc[16:80,,]
#   dimnames(lmod_inc) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   lmod_inc <- lmod_inc %>%
#     cubelyr::as.tbl_cube(met_name = "leapfrog") %>%
#     as_tibble 
#   inc_comparison <- full_join(lmod_inc, mod_inc)  %>% 
#     group_by(age, sex, year) %>% 
#     mutate(pct_diff = abs(EPPASM - leapfrog)  / leapfrog, leapfrog = leapfrog, EPPASM = EPPASM) %>%
#     mutate(pct_diff = ifelse(is.nan(pct_diff), 0, pct_diff), leapfrog = leapfrog, EPPASM = EPPASM)
#   
#   ##assumed less than 1% difference
#   expect_true(all(abs(inc_comparison$pct_diff) < 0.01))
#   
#   ##HIV DEATHS
#   mod_deaths <- attr(mod, 'hivdeaths')
#   mod_deaths <- mod_deaths[-66,,]
#   dimnames(mod_deaths) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   mod_deaths <- mod_deaths %>%
#     cubelyr::as.tbl_cube(met_name = "EPPASM") %>%
#     as_tibble 
#   
#   lmod_deaths <- lmod$hivdeaths
#   lmod_deaths <- lmod_deaths[16:80,,]
#   dimnames(lmod_deaths) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   lmod_deaths <- lmod_deaths %>%
#     cubelyr::as.tbl_cube(met_name = "leapfrog") %>%
#     as_tibble 
#   deaths_comparison <- full_join(lmod_deaths, mod_deaths)  %>% 
#     group_by(age, sex, year) %>% 
#     mutate(pct_diff = abs(EPPASM - leapfrog)  / leapfrog, leapfrog = leapfrog, EPPASM = EPPASM) %>%
#     mutate(pct_diff = ifelse(is.nan(pct_diff), 0, pct_diff), leapfrog = leapfrog, EPPASM = EPPASM)
# 
#   ##percent difference as less than 1%, this is failing rn bc
#   ##the transmission model hasn't been fully implemented
#   ##write.csv(deaths_comparison, "./tests/testdata/plotting_diagnostics/deaths_comparison.csv", row.names = F)
#   expect_true(all(abs(deaths_comparison$pct_diff) < 0.01))
#   
# 
# })
# 
# test_that("Leapfrog matches direct incidence option in EPP-ASM, ART", {
#   ## Check that prevalence, deaths and incidence  matches between
#   ## the two models
#   pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ")
#   
#   demp <- prepare_leapfrog_demp(pjnz1)
#   hivp <- prepare_leapfrog_projp(pjnz1)
#   
#   ## Replace netmigr with unadjusted age 0-4 netmigr, which are not
#   ## in EPP-ASM preparation
#   demp$netmigr <- read_netmigr(pjnz1, adjust_u5mig = FALSE)
#   demp$netmigr_adj <- adjust_spectrum_netmigr(demp$netmigr)
#   
#   lmod <- leapfrogR(demp, hivp)
#   
#   expect_warning(fp <- eppasm::prepare_directincid(pjnz1),
#                  "no non-missing arguments to min; returning Inf")
#   fp$tARTstart <- 61L
#   
#   ## Replace ASFR because demp$asfr is normalised, but fp$asfr is not
#   fp$asfr <- demp$asfr
#   
#   mod <- eppasm::simmod(fp)
#   
#   ##EPPASM and leapfrog are done at different age specifications
#   age_mapping <- data.frame(eppasm_age = c(rep(1,2), rep(2,2), unlist(lapply(c(3:8), rep, times = 5)), rep(9, 32)), leapfrog_age = 1:66)
#   
#   ## PREVALENCE
#   mod_prev <- attr(mod, 'hivpop')
#   dimnames(mod_prev) <- list('CD4' = c(1:7), 'eppasm_age' = c(1:9), 'sex' = c(1:2), 'year' = c(1:61))
#   mod_prev <- mod_prev %>%
#     cubelyr::as.tbl_cube(met_name = "EPPASM") %>%
#     as_tibble %>%
#     group_by(year, eppasm_age, CD4) %>% 
#     summarise(across(`EPPASM`, sum))  
#   
#   lmod_prev <- lmod$hivstrat_adult
#   dimnames(lmod_prev) <- list('CD4' = c(1:7), 'leapfrog_age' = c(1:66), 'sex' = c(1:2), 'year' = c(1:61))
#   lmod_prev <- lmod_prev %>%
#     cubelyr::as.tbl_cube(met_name = "leapfrog") %>%
#     as_tibble
#   lmod_prev <- left_join(lmod_prev, age_mapping, by = 'leapfrog_age')
#   lmod_prev <- lmod_prev %>%
#     group_by(year, eppasm_age, CD4) %>% 
#     summarise(across(`leapfrog`, sum))  
#   
#   prev_comparison <- full_join(lmod_prev, mod_prev)  %>% 
#     group_by(year, eppasm_age, CD4) %>% 
#     mutate(pct_diff = abs(EPPASM - leapfrog)  / leapfrog, leapfrog = leapfrog, EPPASM = EPPASM) %>%
#     mutate(pct_diff = ifelse(is.nan(pct_diff), 0, pct_diff), leapfrog = leapfrog, EPPASM = EPPASM)
#   
#   ##percent difference as less than 1%, this is failing rn bc
#   ##the transmission model hasn't been fully implemented
#   ##write.csv(prev_comparison, "./tests/testdata/plotting_diagnostics/prev_comparison.csv", row.names = F)
#   expect_true(all(abs(prev_comparison$pct_diff) < 0.01))
#   
#   ##INCIDENCE
#   mod_inc <- attr(mod, 'infections')
#   mod_inc <- mod_inc[-66,,]
#   dimnames(mod_inc) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   mod_inc <- mod_inc %>%
#     cubelyr::as.tbl_cube(met_name = "EPPASM") %>%
#     as_tibble 
#   
#   lmod_inc <- lmod$infections
#   lmod_inc <- lmod_inc[16:80,,]
#   dimnames(lmod_inc) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   lmod_inc <- lmod_inc %>%
#     cubelyr::as.tbl_cube(met_name = "leapfrog") %>%
#     as_tibble 
#   inc_comparison <- full_join(lmod_inc, mod_inc)  %>% 
#     group_by(age, sex, year) %>% 
#     mutate(pct_diff = abs(EPPASM - leapfrog)  / leapfrog, leapfrog = leapfrog, EPPASM = EPPASM) %>%
#     mutate(pct_diff = ifelse(is.nan(pct_diff), 0, pct_diff), leapfrog = leapfrog, EPPASM = EPPASM)
#   
#   ##assumed less than 1% difference
#   expect_true(all(abs(inc_comparison$pct_diff) < 0.01))
#   
#   ##HIV DEATHS
#   mod_deaths <- attr(mod, 'hivdeaths')
#   mod_deaths <- mod_deaths[-66,,]
#   dimnames(mod_deaths) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   mod_deaths <- mod_deaths %>%
#     cubelyr::as.tbl_cube(met_name = "EPPASM") %>%
#     as_tibble 
#   
#   lmod_deaths <- lmod$hivdeaths
#   lmod_deaths <- lmod_deaths[16:80,,]
#   dimnames(lmod_deaths) <- list(age = c(1:65), sex = c(1,2), year = c(1:61))
#   lmod_deaths <- lmod_deaths %>%
#     cubelyr::as.tbl_cube(met_name = "leapfrog") %>%
#     as_tibble 
#   deaths_comparison <- full_join(lmod_deaths, mod_deaths)  %>% 
#     group_by(age, sex, year) %>% 
#     mutate(pct_diff = abs(EPPASM - leapfrog)  / leapfrog, leapfrog = leapfrog, EPPASM = EPPASM) %>%
#     mutate(pct_diff = ifelse(is.nan(pct_diff), 0, pct_diff), leapfrog = leapfrog, EPPASM = EPPASM)
#   
#   ##percent difference as less than 1%, this is failing rn bc
#   ##the transmission model hasn't been fully implemented
#   ##write.csv(deaths_comparison, "./tests/testdata/plotting_diagnostics/deaths_comparison.csv", row.names = F)
#   expect_true(all(abs(deaths_comparison$pct_diff) < 0.01))
#   
#   
# })
# 
# 
# test_that("Dimensions for ages are all equal",{
#   pjnz1 <- test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ")
#   demp1 <- prepare_leapfrog_demp(pjnz1)
#   hivp1 <- prepare_leapfrog_projp(pjnz1)
# 
#   ##full pop
#   basepop_ages <- dim(demp1$basepop)[1]
#   mx_ages <- dim(demp1$mx)[1]
#   netmigr_ages <- dim(demp1$netmigr)[1]
#   Sx_ages <- dim(demp1$Sx)[1] -1
#   netmigr_adg_ages <- dim(demp1$netmigr_adj)[1]
# 
#   age_dim <- c(basepop_ages, mx_ages, netmigr_ages, Sx_ages, netmigr_adg_ages)
#   expect_equal(length(unique(age_dim)), 1)
# 
#   ##adult pop
#   asfr_ages <- dim(demp1$asfr)[1]
#   asfd_ages <- dim(demp1$asfd)[1]
# 
#   expect_equal(asfr_ages, asfd_ages)
# })

## check a direct incidence with 0 input incidence
## check that incidence pattern looks normal, following incrr_age