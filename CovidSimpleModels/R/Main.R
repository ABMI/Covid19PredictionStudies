# Copyright 2020 Observational Health Data Sciences and Informatics
#
# This file is part of CovidVulnerabilityIndex
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Execute the Study
#'
#' @details
#' This function executes the CovidVulnerabilityIndex Study.
#' 
#' @param connectionDetails    An object of type \code{connectionDetails} as created using the
#'                             \code{\link[DatabaseConnector]{createConnectionDetails}} function in the
#'                             DatabaseConnector package.
#' @param usePackageCohorts    Whether to use the T and O cohorts in the package
#' @param newTargetCohortId    The cohort id of a new target cohort (must be generated already)
#' @param newTargetCohortName  The name of the new target cohort
#' @param newOutcomeCohortId   The cohort id of a new outcome (must be generated already)  
#' @param newOutcomeCohortName  The name of the new outcome cohort
#' @param newCohortDatabaseSchema Schema name where the new cohorts are stored. Note that for SQL Server, 
#'                             this should include both the database and schema name, for example 'cdm_data.dbo'.
#' @param newCohortTable       The name of the table that contains the target population and outcome
#'                             cohorts.                      
#' @param cdmDatabaseSchema    Schema name where your patient-level data in OMOP CDM format resides.
#'                             Note that for SQL Server, this should include both the database and
#'                             schema name, for example 'cdm_data.dbo'.
#' @param cdmDatabaseName      Shareable name of the database 
#' @param cohortDatabaseSchema Schema name where the cohorts and covariates are created. You will need to have
#'                             write priviliges in this schema. Note that for SQL Server, this should
#'                             include both the database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTable          The name of the table that will be created in the work database schema.
#'                             This table will hold the target, outcome and variable cohorts used in this study.
#' @param oracleTempSchema     Should be used in Oracle to specify a schema where the user has write
#'                             priviliges for storing temporary tables.
#' @param sampleSize           How many patients to sample from the target population                             
#' @param riskWindowStart      The start of the risk window (in days) relative to the startAnchor.                           
#' @param startAnchor          The anchor point for the start of the risk window. Can be "cohort start" or "cohort end".
#' @param riskWindowEnd        The end of the risk window (in days) relative to the endAnchor parameter
#' @param endAnchor            The anchor point for the end of the risk window. Can be "cohort start" or "cohort end".
#' @param firstExposureOnly    Should only the first exposure per subject be included? Note that this is typically done in the createStudyPopulation function,
#' @param removeSubjectsWithPriorOutcome Remove subjects that have the outcome prior to the risk window start?
#' @param priorOutcomeLookback How many days should we look back when identifying prior outcomes?
#' @param requireTimeAtRisk    Should subject without time at risk be removed?
#' @param minTimeAtRisk        The minimum number of days at risk required to be included
#' @param includeAllOutcomes   (binary) indicating whether to include people with outcomes who are not observed for the whole at risk period
#' @param standardCovariates   Use this to add standard covariates such as age/gender
#' @param endDay               The end day relative to index for the custom covariates (default is -1)
#' @param outputFolder         Name of local folder to place results; make sure to use forward slashes
#'                             (/). Do not use a folder on a network drive since this greatly impacts
#'                             performance.
#' @param createCohorts        Create the cohortTable table with the cohorts
#' @param predictSevereAtOutpatientVisit   Run the model that predicts 30 day risk of hospitalization during outpatient visit with covid or flu
#' @param predictCriticalAtOutpatientVisit  Run the model that predicts 30 day risk of hospitalization requring critical care during outpatient visit with covid or flu
#' @param predictDeathAtOutpatientVisit   Run the model that predicts 30 day risk of death during outpatient visit with covid or flu
#' @param predictCriticalAtInpatientVisit  Run the model that predicts requring critical care during initial inpatient visit
#' @param predictDeathAtInpatientVisit  Run the model that predicts 30 day risk of death at initial inpatient visit with pneumonia
#' @param packageResults       Should results be packaged for later sharing?     
#' @param minCellCount         The minimum number of subjects contributing to a count before it can be included 
#'                             in packaged results.
#' @param verbosity            Sets the level of the verbosity. If the log level is at or higher in priority than the logger threshold, a message will print. The levels are:
#'                                         \itemize{
#'                                         \item{DEBUG}{Highest verbosity showing all debug statements}
#'                                         \item{TRACE}{Showing information about start and end of steps}
#'                                         \item{INFO}{Show informative information (Default)}
#'                                         \item{WARN}{Show warning messages}
#'                                         \item{ERROR}{Show error messages}
#'                                         \item{FATAL}{Be silent except for fatal errors}
#'                                         }                              
#' @param cdmVersion           The version of the common data model                             
#'
#' @examples
#' \dontrun{
#' connectionDetails <- createConnectionDetails(dbms = "postgresql",
#'                                              user = "joe",
#'                                              password = "secret",
#'                                              server = "myserver")
#'
#' execute(connectionDetails,
#'         cdmDatabaseSchema = "cdm_data",
#'         cdmDatabaseName = 'shareable name of the database'
#'         cohortDatabaseSchema = "study_results",
#'         cohortTable = "cohort",
#'         outcomeId = 1,
#'         oracleTempSchema = NULL,
#'         riskWindowStart = 1,
#'         startAnchor = 'cohort start',
#'         riskWindowEnd = 365,
#'         endAnchor = 'cohort start',
#'         standardCovariates = FeatureExtraction::createCovariateSettings(useDemographicsAgeGroup = T, useDemographicsGender = T),
#'         outputFolder = "c:/temp/study_results", 
#'         createCohorts = T,
#'         predictSevereAtOutpatientVisit = T,
#'         predictCriticalAtOutpatientVisit = T,
#'         predictDeathAtOutpatientVisit = T,
#'         predictCriticalAtInpatientVisit = T,
#'         predictDeathAtInpatientVisit = T,
#'         packageResults = F,
#'         minCellCount = 10,
#'         verbosity = "INFO",
#'         cdmVersion = 5)
#' }
#'
#' @export
execute <- function(connectionDetails,
                    usePackageCohorts = T,
                    newTargetCohortId = NULL,
                    newTargetCohortName = '',
                    newOutcomeCohortId = NULL,
                    newOutcomeCohortName = '',
                    newCohortDatabaseSchema = NULL,
                    newCohortTable = NULL,
                    cdmDatabaseSchema,
                    cdmDatabaseName = 'friendly database name',
                    cohortDatabaseSchema = cdmDatabaseSchema,
                    cohortTable = "cohort",
                    oracleTempSchema = cohortDatabaseSchema,
                    sampleSize = NULL,
                    riskWindowStart = 1,
                    startAnchor = 'cohort start',
                    riskWindowEnd = 365,
                    endAnchor = 'cohort start',
                    firstExposureOnly = F,
                    removeSubjectsWithPriorOutcome = F,
                    priorOutcomeLookback = 99999,
                    requireTimeAtRisk = F,
                    minTimeAtRisk = 1,
                    includeAllOutcomes = T,
                    standardCovariates,
                    endDay = -1,
                    outputFolder,
                    createCohorts = F,
                    predictSevereAtOutpatientVisit = F,
                    predictCriticalAtOutpatientVisit = F,
                    predictDeathAtOutpatientVisit = F,
                    predictCriticalAtInpatientVisit = F,
                    predictDeathAtInpatientVisit = F,
                    packageResults = F,
					          minCellCount = 10,
                    verbosity = "INFO",
                    cdmVersion = 5) {
  
  if(usePackageCohorts){
    if(!is.null(newTargetCohortId) | !is.null(newOutcomeCohortId)){
      warning('Input new ids but also said to use package cohorts - new ids will be ignored...')
    }
  }
  
  if (!file.exists(file.path(outputFolder,cdmDatabaseName)))
    dir.create(file.path(outputFolder,cdmDatabaseName), recursive = TRUE)
  
  ParallelLogger::addDefaultFileLogger(file.path(outputFolder,cdmDatabaseName, "log.txt"))
  
  if (createCohorts) {
    ParallelLogger::logInfo("Creating cohorts")
    createCohorts(connectionDetails = connectionDetails,
                  cdmDatabaseSchema = cdmDatabaseSchema,
                  cohortDatabaseSchema = cohortDatabaseSchema,
                  cohortTable = cohortTable,
                  usePackageCohorts = usePackageCohorts,
                  newTargetCohortId = newTargetCohortId,
                  newOutcomeCohortId = newOutcomeCohortId,
                  newCohortDatabaseSchema = newCohortDatabaseSchema,
                  newCohortTable = newCohortTable,
                  oracleTempSchema = oracleTempSchema,
                  outputFolder = file.path(outputFolder, cdmDatabaseName)) 
  }
  
  studyAnalyses <- getSettings(predictSevereAtOutpatientVisit = predictSevereAtOutpatientVisit,
                               predictCriticalAtOutpatientVisit = predictCriticalAtOutpatientVisit,
                               predictDeathAtOutpatientVisit = predictDeathAtOutpatientVisit,
                               predictCriticalAtInpatientVisit = predictCriticalAtInpatientVisit,
                               predictDeathAtInpatientVisit = predictDeathAtInpatientVisit,
                               usePackageCohorts = usePackageCohorts)
  if( nrow(studyAnalyses)!=0){
    for(k in 1:nrow(studyAnalyses)){
      sa <- studyAnalyses[k,]
      result <- predictSimple(connectionDetails = connectionDetails,
                    cohortId = sa$cohortId,
                    outcomeIds = sa$outcomeId,
                    model = sa$model,
                    analysisId = sa$analysisId,
                    cdmDatabaseSchema = cdmDatabaseSchema,
                    cdmDatabaseName = cdmDatabaseName,
                    cohortDatabaseSchema = cohortDatabaseSchema,
                    cohortTable = cohortTable,
                    oracleTempSchema = oracleTempSchema,
                    standardCovariates = standardCovariates,
                    endDay = endDay,
                    sampleSize = sampleSize,
                    cdmVersion = cdmVersion,
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
      
      if(!dir.exists(file.path(outputFolder,cdmDatabaseName, paste0('Analysis_',sa$analysisId)))){
        dir.create(file.path(outputFolder,cdmDatabaseName,paste0('Analysis_',sa$analysisId)), recursive = T)
      }
      ParallelLogger::logInfo("Saving results")
      saveRDS(result, file.path(outputFolder,cdmDatabaseName,paste0('Analysis_',sa$analysisId),'validationResults.rds'))
      ParallelLogger::logInfo(paste0("Results saved to:",file.path(outputFolder,cdmDatabaseName,paste0('Analysis_',sa$analysisId),'validationResults.rds')))
      
    }
  }
  
  if (packageResults) {
    ParallelLogger::logInfo("Packaging results")
    packageResults(outputFolder = outputFolder,
                   cdmDatabaseName = cdmDatabaseName,
                   minCellCount = minCellCount)
  }
   
  return(TRUE)
}




