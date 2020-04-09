getData <- function(connectionDetails,
                    cohortId,
                    outcomeIds,
                    cdmDatabaseSchema,
                    cdmDatabaseName,
                    cohortDatabaseSchema,
                    cohortTable,
                    oracleTempSchema,
                    standardCovariates,
                    endDay,
                    firstExposureOnly,
                    sampleSize,
                    cdmVersion){
  
  pathToCustom <- system.file("settings", 'CustomCovariates.csv', package = "CovidSimpleSurvival")
  cohortVarsToCreate <- utils::read.csv(pathToCustom)
  covSets <- list()
  length(covSets) <- nrow(cohortVarsToCreate)+1
  covSets[[1]] <- standardCovariates
  
  for(i in 1:nrow(cohortVarsToCreate)){
    covSets[[1+i]] <- createCohortCovariateSettings(covariateName = as.character(cohortVarsToCreate$cohortName[i]),
                                                    covariateId = cohortVarsToCreate$cohortId[i]*1000+456,
                                                    cohortDatabaseSchema = cohortDatabaseSchema,
                                                    cohortTable = cohortTable,
                                                    cohortId = cohortVarsToCreate$atlasId[i],
                                                    startDay=cohortVarsToCreate$startDay[i], 
                                                    endDay=endDay,
                                                    count= as.character(cohortVarsToCreate$count[i]), 
                                                    ageInteraction = as.character(cohortVarsToCreate$ageInteraction[i]))
  }
  
  result <- PatientLevelPrediction::getPlpData(connectionDetails = connectionDetails,
                                               cdmDatabaseSchema = cdmDatabaseSchema,
                                               oracleTempSchema = oracleTempSchema, 
                                               cohortId = as.double(as.character(cohortId)), 
                                               outcomeIds = as.double(as.character(outcomeIds)), 
                                               cohortDatabaseSchema = cohortDatabaseSchema, 
                                               outcomeDatabaseSchema = cohortDatabaseSchema, 
                                               cohortTable = cohortTable, 
                                               outcomeTable = cohortTable, 
                                               cdmVersion = cdmVersion, 
                                               firstExposureOnly = firstExposureOnly, 
                                               sampleSize =  sampleSize, 
                                               covariateSettings = covSets)
  
  return(result)
  
}


getModel <- function(model = 'SimpleModel.csv'){
  pathToCustom <- system.file("settings", model , package = "CovidSimpleSurvival")
  coefficients <- utils::read.csv(pathToCustom)
  return(coefficients)
}

getPredict <- function(model){
  predictExisting <- function(plpData, population){
    coefficients <- model
    
    prediction <- merge(plpData$covariates, ff::as.ffdf(coefficients), by = "covariateId")
    prediction$value <- prediction$covariateValue * prediction$points
    prediction <- PatientLevelPrediction:::bySumFf(prediction$value, prediction$rowId)
    colnames(prediction) <- c("rowId", "value")
    prediction <- merge(population, prediction, by ="rowId", all.x = TRUE)
    prediction$value[is.na(prediction$value)] <- 0
    
    # add any final mapping here (e.g., add intercept and mapping)
    prediction$value <- prediction$value + model$points[model$covariateId==0]
    prediction$value <- prediction$value/10
    prediction$value <- 1/(1+exp(-1*prediction$value))
    
    scaleVal <- max(prediction$value)
    if(scaleVal>1){
      prediction$value <- prediction$value/scaleVal
    }
    
    attr(prediction, "metaData") <- list(predictionType = 'binary', scale = scaleVal)
    
    return(prediction)
  }
  return(predictExisting)
}


getSettings <- function(predictShortTermSurvival,
                        predictLongTermSurvival, 
                        usePackageCohorts){
  
  settingR <- c()
  
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "CovidSimpleSurvival")
  cohortsToCreate <- utils::read.csv(pathToCsv, stringsAsFactors = F)
  
  toList<- expand.grid(cohortsToCreate$cohortId[cohortsToCreate$type=='target'], cohortsToCreate$cohortId[cohortsToCreate$type=='outcome'])
  colnames(toList) <- c("cohortId", "outcomeId")
  
  if(predictShortTermSurvival){
    settingR <- rbind(data.frame(settingR), 
                      data.frame(cohortId = ifelse(usePackageCohorts, toList$cohortId, 1), 
                                 outcomeId = ifelse(usePackageCohorts, toList$outcomeId, 2),
                                 analysisId = ifelse(usePackageCohorts, rep(6,nrow(toList)), 6), 
                                 model = ifelse(usePackageCohorts, rep('ShortTermSurvival.csv',nrow(toList)), 'ShortTermSurvival.csv'), 
                                 riskWindowStart = ifelse(usePackageCohorts, rep(0,nrow(toList)), 0),
                                 startAnchor = rep('cohort start',nrow(toList)),
                                 riskWindowEnd = ifelse(usePackageCohorts, rep(60,nrow(toList)), 60),
                                 endAnchor = rep('cohort start',nrow(toList))
                      )
    )
  }
  if(predictLongTermSurvival){
    settingR <- rbind(settingR, c(cohortId = ifelse(usePackageCohorts, 1, 1), 
                                  outcomeId = ifelse(usePackageCohorts, 2, 2),
                                  analysisId = 7, 
                                  model = 'LongTermSurvival.csv',
                                  riskWindowStart = 0,
                                  startAnchor = 'cohort start',
                                  riskWindowEnd = 365,
                                  endAnchor = 'cohort start'))
    settingR <- rbind(data.frame(settingR), 
                      data.frame(cohortId = ifelse(usePackageCohorts, toList$cohortId, 1), 
                                 outcomeId = ifelse(usePackageCohorts, toList$outcomeId, 2),
                                 analysisId = ifelse(usePackageCohorts, rep(7,nrow(toList)), 7), 
                                 model = ifelse(usePackageCohorts, rep('LongTermSurvival.csv',nrow(toList)), 'LongTermSurvival.csv'), 
                                 riskWindowStart = ifelse(usePackageCohorts, rep(0,nrow(toList)), 0),
                                 startAnchor = ifelse(usePackageCohorts,rep('cohort start',nrow(toList)), 'cohort start'),
                                 riskWindowEnd = ifelse(usePackageCohorts, rep(365,nrow(toList)), 365),
                                 endAnchor = rep('cohort start',nrow(toList))
                      )
    )
  }
  
  settingR <- as.data.frame(settingR)
  return(settingR)
  
}


