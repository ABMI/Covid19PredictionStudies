library(CovidSimpleSurvival)
# USER INPUTS
#=======================
# Specify where the temporary files (used by the ff package) will be created:
options(fftempdir = "location with space to save big data")

# The folder where the study intermediate and result files will be written:
outputFolder <- "./CovidSimpleSurvivalResults"

# Details for connecting to the server:
dbms <- "you dbms"
user <- 'your username'
pw <- 'your password'
server <- 'your server'
port <- 'your port'

connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = dbms,
                                                                server = server,
                                                                user = user,
                                                                password = pw,
                                                                port = port)

# Add the database containing the OMOP CDM data
cdmDatabaseSchema <- 'cdm database schema'
# Add a database with read/write access as this is where the cohorts will be generated
cohortDatabaseSchema <- 'work database schema'
oracleTempSchema <- NULL

# table name where the cohorts will be generated
cohortTable <- 'CovidSimpleSurvivalCohort'

#============== Pick Study Parts To Run: ===========
createCohorts = TRUE
predictShortTermSurvival = TRUE
predictLongTermSurvival = TRUE
packageResults = TRUE

minCellCount <- 5
sampleSize <- NULL


#============== Pick T and O cohorts ===========

# [option 1] use default cohorts (same as those used to develop the models)
usePackageCohorts <- TRUE

# [option 2] use your own specified cohorts (these must be created already)
# uncomment the line below to set usePackageCohorts as False
## usePackageCohorts <- False
newTargetCohortId = 'the cohort definition id for the new target' 
newOutcomeCohortId = 'the cohort definition id for the new outcome' 
newCohortDatabaseSchema = 'the database schema containing the cohorts' 
newCohortTable = 'the table name of the table containing the cohorts' 

#=======================
# TAR settings - recommended to not edit
#=======================
riskWindowStart <- NULL # uses default of 0
startAnchor <- NULL # uses default of 'cohort start'
riskWindowEnd <- NULL # uses default of 60 for short term and 365 for long term
endAnchor <- NULL # uses default of 'cohort start'
firstExposureOnly <- F
removeSubjectsWithPriorOutcome <- F
priorOutcomeLookback <- 99999
requireTimeAtRisk <- F
minTimeAtRisk <- 1
includeAllOutcomes <- T

standardCovariates <- FeatureExtraction::createCovariateSettings(useDemographicsAge= T, 
                                                                 useDemographicsGender = T, 
                                                                 excludedCovariateConceptIds = 8532 )

execute(connectionDetails = connectionDetails,
        usePackageCohorts = usePackageCohorts,
        newTargetCohortId = newTargetCohortId,
        newOutcomeCohortId = newOutcomeCohortId,
        newCohortDatabaseSchema = newCohortDatabaseSchema,
        newCohortTable = newCohortTable,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cdmDatabaseName = cdmDatabaseName,
        cohortDatabaseSchema = cohortDatabaseSchema,
        cohortTable = cohortTable,
        sampleSize = sampleSize,
        riskWindowStart = riskWindowStart,
        startAnchor = startAnchor,
        riskWindowEnd = riskWindowEnd,
        endAnchor = endAnchor,
        firstExposureOnly = firstExposureOnly,
        removeSubjectsWithPriorOutcome = removeSubjectsWithPriorOutcome,
        priorOutcomeLookback = priorOutcomeLookback,
        requireTimeAtRisk = requireTimeAtRisk,
        minTimeAtRisk = minTimeAtRisk,
        includeAllOutcomes = includeAllOutcomes,
        standardCovariates = standardCovariates,
        outputFolder = outputFolder,
        createCohorts = createCohorts,
        predictShortTermSurvival = predictShortTermSurvival,
        predictLongTermSurvival = predictLongTermSurvival,
        packageResults = packageResults,
        minCellCount = minCellCount,
        verbosity = "INFO",
        cdmVersion = 5)