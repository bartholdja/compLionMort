load("/Users/Viktualia/Dropbox/Projects/003_Lions/RScripts/dataCleaning/workspaceLions.rdata")

lion <- lion[!duplicated(lion$ID), ]
lion <- lion[lion$location== "serengeti", ]
lion <- lion[!is.na(lion$birth), ]
lion <- lion[!is.na(lion$lastSeen), ]
lion <- lion[(lion$sex == "f" | lion$sex == "m" | lion$sex == "u"), ]
lion <- lion[!is.na(lion$ID), ]

# downsize the data
idKeep <-rbinom(nrow(lion), 1, 0.2)
lion <- lion[idKeep == 1, ]

# Extract variables:
study <- julian(as.Date(c("1967-01-01", "2010-05-15")), 
                origin = as.POSIXct("1970-01-01"))
n <- nrow(lion)
birth <- julian(as.Date(lion[, "birth"]), 
                origin = as.POSIXct("1970-01-01"))
last <- julian(as.Date(lion[, "lastSeen"]), 
               origin = as.POSIXct("1970-01-01"))
death <- rep(NA, length(last))
death[lion$lsinfo == "d" | lion$lsinfo == "d*" | lion$lsinfo == "**d"] <- 1
idNoDeath <- which(is.na(death))
first <- rep(NA, n)
sex <- as.character(lion[, 'sex'])
idLeftTr <- which(lion$birth <=  "1967-01-01")
ageTrunc <- apply(cbind(study[1] - birth, 0), 1, max) / 365.25
ageToLast <- (last - birth) / 365.25
# in Pusey & Packer 1987 all males dispersed by the age of 4.2, minimum age 1.8 (1 ind out of 12)
idMigr <- which(sex == 'm' & is.na(death) &
                  ageToLast >= 1.75 & ageToLast <= 4.25)  
idNonMigr <- which(!is.na(birth))
idNonMigr <- idNonMigr[!(idNonMigr %in% idMigr)] # n = sum(!is.na(death)) + length(idMigr) + length(idNonMigr), everyone accounted for
idNoSex <- which(sex == "u")
probFem <- 0.45


