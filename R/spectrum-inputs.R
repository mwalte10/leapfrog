
read_sx <- function(pjnz, use_ep5=FALSE){

  if(use_ep5) {
    dpfile <- grep(".ep5$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  } else {
    dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  }

  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)

  exists_dptag <- function(tag, tagcol=1){tag %in% dp[,tagcol]}
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  version <- paste("Spectrum", dp[which(dp[,1] == "<ValidVers MV>")+2, 4])

  ## projection parameters
  yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
  
  ## mx
  Sx <- dpsub("<SurvRate MV2>", 3+c(0:81, 82+0:81), timedat.idx)
  Sx <- array(as.numeric(unlist(Sx)), c(82, 2, length(proj.years)))
  dimnames(Sx) <- list(age=c(0:80, "80+"), sex=c("male", "female"), year=proj.years)

  Sx
}

#' Prepare demographic inputs from Spectrum PJNZ
#'
#' @param pjnz path to PJNZ file
#'
#' @return list of demographic input parameters
#' 
#' @export
prepare_leapfrog_demp <- function(pjnz) {

  demp <- eppasm::read_specdp_demog_param(pjnz)
  demp$Sx <- read_sx(pjnz)

  ## Spectrum adjusts net-migration to occur half in
  ## current age group and half in next age group
  
  demp$netmigr_adj <- demp$netmigr
  demp$netmigr_adj[-1,,] <- (demp$netmigr[-1,,] + demp$netmigr[-81,,])/2
  demp$netmigr_adj[1,,] <- demp$netmigr[1,,]/2
  demp$netmigr_adj[81,,] <- demp$netmigr_adj[81,,] + demp$netmigr[81,,]/2

  demp$births_sex_prop <- rbind(male = demp$srb, female = 100) / (demp$srb + 100)

  demp
}

#' Prepare adult HIV projection parameters from Spectrum PJNZ
#'
#' @param pjnz path to PJNZ file
#' @param hiv_steps_per_year number of Euler integration steps per year; default 10
#' @param hTS number of HIV treatment stages; default 3 (0-5 months,
#'   6-11 months, 12+ months)
#'
#' @return list of HIV projection parameters
#'
#' @export
prepare_leapfrog_projp <- function(pjnz, hiv_steps_per_year = 10L, hTS = 3) {

  projp <- eppasm::read_hivproj_param(pjnz)

  ## Hard coded to expand age groups 15-24, 25-34, 35-44, 45+ to
  ## single-year ages 15:80.
  ## Requires extension for coarse HIV age group stratification
  idx_expand_full <- rep(1:4, times = c(10, 10, 10, 36))
  idx_expand_coarse <- rep(1:4, times = c(3, 2, 2, 2))

  v <- list()
  v$incidinput <- eppasm::read_incid_input(pjnz)
  v$incidpopage <- attr(v$incidinput, "incidpopage")    
  v$incrr_sex <- projp$incrr_sex
  
  ## Use Beer's coefficients to distribution IRRs by age/sex
  Amat <- eppasm:::create_beers(17)
  v$incrr_age <- apply(projp$incrr_age, 2:3, function(x)  Amat %*% x)[16:81, , ] ## !! Hard coded
  v$incrr_age[v$incrr_age < 0] <- 0
  
  v$cd4_initdist_full <- projp$cd4_initdist[ , idx_expand_full, ]
  v$cd4_prog_full <- (1-exp(-projp$cd4_prog[ , idx_expand_full, ] / hiv_steps_per_year)) * hiv_steps_per_year
  v$cd4_mort_full <- projp$cd4_mort[ ,idx_expand_full, ]
  v$art_mort_full <- projp$art_mort[c(1, 2, rep(3, hTS - 2)), , idx_expand_full, ]

  v$cd4_initdist_coarse <- projp$cd4_initdist[ , idx_expand_coarse, ]
  v$cd4_prog_coarse <- (1-exp(-projp$cd4_prog[ , idx_expand_coarse, ] / hiv_steps_per_year)) * hiv_steps_per_year
  v$cd4_mort_coarse <- projp$cd4_mort[ ,idx_expand_coarse, ]
  v$art_mort_coarse <- projp$art_mort[c(1, 2, rep(3, hTS - 2)), , idx_expand_coarse, ]
  
  v$artmx_timerr <- projp$artmx_timerr[c(1, 2, rep(3, hTS - 2)), ]

  ## ## ART eligibility and numbers on treatment

  v$art15plus_num <- projp$art15plus_num
  v$art15plus_isperc <- projp$art15plus_numperc == 1

  ## convert percentage to proportion
  v$art15plus_num[v$art15plus_isperc] <- v$art15plus_num[v$art15plus_isperc] / 100

  
  ## eligibility starts in projection year idx
  ## ## !! NOTE: from EPP-ASM; not yet implemented  
  ## v$specpop_percelig <- rowSums(with(projp$artelig_specpop[-1,], mapply(function(elig, percent, year) rep(c(0, percent*as.numeric(elig)), c(year - proj_start, proj_end - year + 1)), elig, percent, year)))
  v$artcd4elig_idx <- findInterval(-projp$art15plus_eligthres, -c(999, 500, 350, 250, 200, 100, 50))

  ## Update eligibility threshold from CD4 <200 to <250 to account for additional
  ## proportion eligible with WHO Stage 3/4.
  v$artcd4elig_idx <- replace(v$artcd4elig_idx, v$artcd4elig_idx==5L, 4L)

  v$pw_artelig <- with(projp$artelig_specpop["PW",], rep(c(0, elig), c(year - projp$yr_start, projp$yr_end - year + 1)))  # are pregnant women eligible (0/1)

  ## ## !! NOTE: from EPP-ASM; not yet implemented
  ## ## percentage of those with CD4 <350 who are based on WHO Stage III/IV infection
  ## v$who34percelig <- who34percelig

  v$art_dropout <- projp$art_dropout/100

  proj_years <- as.integer(projp$yr_end - projp$yr_start + 1L)
  v$t_ART_start <- min(c(unlist(apply(v$art15plus_num > 0, 1, which)), proj_years))

  ## New ART patient allocation options
  v$art_alloc_method <- projp$art_alloc_method
  v$art_alloc_mxweight <- projp$art_prop_alloc[1]

  ## Scale mortality among untreated population by ART coverage
  v$scale_cd4_mort <- projp$scale_cd4_mort

  ## State space dimensions
  v$hAG_SPAN_full <- rep(1L, 66L)
  v$hAG_SPAN_coarse <- c(2L, 3L, 5L, 5L, 5L, 5L, 5L, 5L, 31L)

  v
}
