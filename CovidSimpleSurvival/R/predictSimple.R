predictSimple <- function(connectionDetails,
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
                                           cdmVersion,
                                           riskWindowStart,
                                           startAnchor,
                                           riskWindowEnd,
                                           endAnchor,
                                           removeSubjectsWithPriorOutcome,
                                           priorOutcomeLookback,
                                           requireTimeAtRisk,
                                           minTimeAtRisk,
                                           includeAllOutcomes,
                                           model = 'ShortTermSurvival.csv',
                                           analysisId = 1){
  ParallelLogger::logInfo("Extracting data for predicting Simple Model")
  plpData <- getData(connectionDetails = connectionDetails,
                     cohortId = cohortId,
                     outcomeIds = outcomeIds,
                     cdmDatabaseSchema = cdmDatabaseSchema,
                     cdmDatabaseName = cdmDatabaseName,
                     cohortDatabaseSchema = cohortDatabaseSchema,
                     cohortTable = cohortTable,
                     oracleTempSchema = oracleTempSchema,
                     standardCovariates = standardCovariates,
                     endDay = endDay,
                     firstExposureOnly = firstExposureOnly,
                     sampleSize = sampleSize,
                     cdmVersion = cdmVersion)
  
  ParallelLogger::logInfo("Creating population")
  ParallelLogger::logInfo(paste0("TAR:", startAnchor,'+',riskWindowStart,'-', endAnchor,'+',riskWindowEnd))
  population <- PatientLevelPrediction::createStudyPopulation(plpData = plpData, 
                                                              outcomeId = outcomeIds,
                                                              riskWindowStart = riskWindowStart,
                                                              startAnchor = startAnchor,
                                                              riskWindowEnd = riskWindowEnd,
                                                              endAnchor = endAnchor,
                                                              firstExposureOnly = firstExposureOnly,
                                                              removeSubjectsWithPriorOutcome = removeSubjectsWithPriorOutcome,
                                                              priorOutcomeLookback = priorOutcomeLookback,
                                                              requireTimeAtRisk = requireTimeAtRisk,
                                                              minTimeAtRisk = minTimeAtRisk,
                                                              includeAllOutcomes = includeAllOutcomes)
  
  
  # apply the model:
  plpModel <- list(model = getModel(model = model),
                   analysisId = paste0('Analysis_',analysisId),
                   hyperParamSearch = NULL,
                   index = NULL,
                   trainCVAuc = NULL,
                   modelSettings = list(model = 'score', modelParameters = NULL),
                   metaData = NULL,
                   populationSettings = attr(population, "metaData"),
                   trainingTime = NULL,
                   varImp = NULL,
                   dense = T,
                   cohortId = cohortId,
                   outcomeId = outcomeIds,
                   endDay = endDay,
                   covariateMap = NULL,
                   predict = getPredict(getModel(model = model))
  )
  class(plpModel) <- 'plpModel'
  ParallelLogger::logInfo("Applying and evaluating model")
  result <- PatientLevelPrediction::applyModel(population = population,
                                               plpData = plpData,
                                               plpModel = plpModel)
  
  result$inputSetting$database <- cdmDatabaseName
  result$executionSummary  <- list()
  result$model <- plpModel
  result$analysisRef <- list()
  result$covariateSummary <- PatientLevelPrediction:::covariateSummary(plpData = plpData, population = population)
  
  result$covariateSummary <- merge(result$covariateSummary, result$model$model[,c('covariateId', 'points')], by ='covariateId')
  result$covariateSummary$covariateValue = result$covariateSummary$points
  
 return(result)
  }

  