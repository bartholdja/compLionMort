rm(list = ls())

setwd("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions")

serenM <- read.csv("data/serengeti/serenMal05Nov.csv",
                   header = T, na.strings = "")
serenM[ , which(names(serenM) == "lsDate")] <- factor(serenM[ , which(names(serenM) == "lsDate")], levels = c(levels(serenM[ , which(names(serenM) == "lsDate")]), "2013-08-01"))

serenM$alive <- rep(0, nrow(serenM))
for (i in 1:nrow(serenM)) {
    if (serenM[i, which(names(serenM) == "lsDate")] == "ok") {
      serenM$alive[i] <- 1
      serenM[i, which(names(serenM) == "lsDate")] <- "2013-08-01" 
  }
}
rm(i)
serenM[ ,8] <- factor(serenM[ ,8], levels = levels(serenM[ ,8])[!levels(serenM[ ,8]) == "ok"])
serenM$birthDate <- as.Date(serenM$birthDate)
serenM$firstSeenRes <- as.Date(serenM$firstSeenRes)
serenM$firstSeenNomNevRes <- as.Date(serenM$firstSeenNomNevRes)
serenM$firstSeenNomLatRes <- as.Date(serenM$firstSeenNomLatRes)
serenM$firstSeenIm <- as.Date(serenM$firstSeenIm)
serenM$lsDate <- as.Date(serenM$lsDate)
#load("data/serengeti/serenMalCens01Aug2013.Rdata")
serenM$ageYrs <- as.numeric((serenM$lsDate - serenM$birthDate)/365.25)
serenM$firstSeenDate <- c("1900-01-01", rep(NA, (nrow(serenM)-1)))
serenM$firstSeenDate <- as.Date(serenM$firstSeenDate)

# make cumulative firstSeenDate (from the 4 different FSdate categories)
for (i in 1:nrow(serenM)){
  if(sum(!is.na(c(serenM$firstSeenRes[i]
                  ,serenM$firstSeenNomNevRes[i],
                  serenM$firstSeenNomLatRes[i], serenM$firstSeenIm[i]))) > 0){
  serenM$firstSeenDate[i] <-  c(serenM$firstSeenRes[i], serenM$firstSeenNomNevRes[i], 
                             serenM$firstSeenNomLatRes[i], 
                             serenM$firstSeenIm[i])[!is.na(c(serenM$firstSeenRes[i]
                             ,serenM$firstSeenNomNevRes[i],
                              serenM$firstSeenNomLatRes[i], serenM$firstSeenIm[i]))]
  }
} 
rm(i)
serenM$datInd <- rep("serenM", nrow(serenM))
save.image("data/serengeti/serenMal05Nov.Rdata")
rm(serenM)

serenF <- read.csv("data/serengeti/serenFem05Nov.csv",
                   header = T, na.strings = "")
serenF[ ,which(names(serenF) == "lsDate")] <- factor(serenF[ ,which(names(serenF) == "lsDate")], levels = c(levels(serenF[ ,which(names(serenF) == "lsDate")]), "2013-08-01"))


serenF$alive <- rep(0, nrow(serenF))
for (i in 1:nrow(serenF)) {
  if (serenF[i,which(names(serenF) == "lsDate")] == "ok") {
    serenF$alive[i] <- 1
    serenF[i,which(names(serenF) == "lsDate")] <- "2013-08-01" 
  }
}
rm(i)
serenF[ ,which(names(serenF) == "lsDate")] <- factor(serenF[ ,which(names(serenF) == "lsDate")], levels = levels(serenF[ ,which(names(serenF) == "lsDate")])[! levels(serenF[ ,which(names(serenF) == "lsDate")]) %in% c("ok", "OK", "Ok", "0k")])

serenF$birthDate <- as.Date(serenF$birthDate)
serenF$firstSeenDate <- as.Date(serenF$firstSeenDate)
serenF$lsDate <- as.Date(serenF$lsDate)
serenF$datInd <- rep("serenF", nrow(serenF))
save.image("data/serengeti/serenFem05Nov.Rdata")
rm(serenF)

rm(list = ls())
serenUX <- read.csv("data/serengeti/serenLostUnseen05Nov.csv",
                   header = T, na.strings = "")
serenUX$birthDate <- as.Date(serenUX$birthDate)
serenUX <- serenUX[which(serenUX$birthDate <= "2013-08-01"), ]
serenUX[ ,4] <- factor(serenUX[ ,4], levels = c(levels(serenUX[ ,4]), "2013-08-01"))

serenUX$alive <- rep(0, nrow(serenUX))
for (i in 1:nrow(serenUX)) {
  if (serenUX[i,4] == "ok") {
    serenUX$alive[i] <- 1
    serenUX[i,4] <- "2013-08-01" 
  }
}
rm(i)
serenUX[ ,4] <- factor(serenUX[ ,4], levels = levels(serenUX[ ,4])[!levels(serenUX[ ,4]) == "ok"])

serenUX$lsDate <- as.Date(serenUX$lsDate)
serenUX$datInd <- rep("serenUX", nrow(serenUX))
save.image("data/serengeti/serenLostUnseen05Nov.Rdata")


rm(list = ls())
serenDead <- read.csv("data/serengeti/serenDeaths05Nov.csv",
                    header = T, na.strings = "")
serenDead$birthDate <- as.Date(serenDead$birthDate)
serenDead$lsDate <- as.Date(serenDead$lsDate)
serenDead$alive <- rep(0, nrow(serenDead))
serenDead$datInd <- rep("serenDead", nrow(serenDead))
save.image("data/serengeti/serenDeaths05Nov.Rdata")
rm(list = ls())

load("data/serengeti/serenDeaths05Nov.Rdata")
load("data/serengeti/serenMal05Nov.Rdata")
load("data/serengeti/serenFem05Nov.Rdata")
load("data/serengeti/serenLostUnseen05Nov.Rdata")

serenM <- serenM[ ,- which(names(serenM) == "ageYrs")]

seren1 <- merge(serenDead, serenM, all = T)
seren2 <- merge(seren1, serenF, all = T)
seren <- merge(seren2, serenUX, all = T)
seren$id <- factor(seren$id, levels = c(levels(seren$id), "sacFem"))
seren$id[which(seren$id == "sac")] [1] <- "sacFem"
seren$id <- factor(seren$id, levels = c(levels(seren$id), paste("lost", 1:length(seren$id[which(seren$id == "lost")]), sep = "")))
seren$id[which(seren$id == "lost")] <- paste("lost", 1:length(seren$id[which(seren$id == "lost")]), sep = "")
seren$id <- factor(seren$id, levels = c(levels(seren$id), paste("wait", 1:length(seren$id[which(seren$id == "wait")]), sep = "")))
seren$id[which(seren$id == "wait")] <- paste("wait", 1:length(seren$id[which(seren$id == "wait")]), sep = "")

double <- seren$id[duplicated(seren$id)]

# transfer the information from "dead" from the observation form the serenDead data to the respektive
# row from the other data set. Turns out the data are really clean and regular, so that I can 
# just copy the information from all the uneven rows to the even rows until the length of double
# is reached
for (i in seq(1, 379, 2)){
  seren$dead[i+1] <- seren$dead[i]
}

seren <- seren[-(seq(1, 379, 2)), ]
# change one entry of sex from m? to u
seren$sex[seren$sex == "m?"] <- "u"

names(seren)
seren$ageLS <- (seren$lsDate - seren$birthDate)/365.25
head(seren$ageLS)
seren <- seren[is.na(seren$nomFirstSeenMoth), ]
seren <- seren[-which(seren$ageLS < 0), ]

write.csv(seren, "data/Serengeti/seren.csv")
rm(list = ls())