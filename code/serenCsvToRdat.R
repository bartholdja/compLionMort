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
save.image("data/serengeti/serenLostUnseen05Nov.Rdata")


rm(list = ls())
serenDead <- read.csv("data/serengeti/serenDeaths05Nov.csv",
                    header = T, na.strings = "")
serenDead$birthDate <- as.Date(serenDead$birthDate)
serenDead$lsDate <- as.Date(serenDead$lsDate)
serenDead$alive <- rep(0, nrow(serenDead))

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

write.csv(seren, "data/Serengeti/seren.csv")