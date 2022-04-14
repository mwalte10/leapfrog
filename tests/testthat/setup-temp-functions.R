read_sx <- function(pjnz, use_ep5=FALSE){

  if(use_ep5) {
    dpfile <- grep(".ep5$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  } else {
    dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  }

  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)

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
