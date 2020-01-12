library(data.table)
library(ggplot2)
library(stringr)

#script----
Metadata[,S1Date:=as.Date.character(S1Date, "%d/%m/%Y")]
Metadata[,Born:=as.Date.character(Born, "%d/%m/%Y")]
Metadata[,Group:=factor(Group, levels = c("1W", "1M", "3M", "AD"))]
#ggplot(Metadata, aes(x=Group, y=HeartWght))+geom_boxplot()

##-technically- more reliable than a dumb group mean but WHO CARES
#b<-lm(HeartWght~Group, data = Metadata)
#summary(b)
#Metadata[is.na(HeartWght) & !is.na(DaysOld), HeartWght:=predict(b, data.frame(Group=Group))]

#c<-lm(BodyWght~DaysOld, data = Metadata); summary(c)

Metadata[,GroupMean:=mean(HeartWght, na.rm = T), Group]
Metadata[is.na(HeartWght), HeartWght:=GroupMean]
#Metadata[,GroupMean:=NULL]

#ggplot(Metadata, aes(x=Group, y=BodyWght))+geom_boxplot()
Metadata[,GroupMean:=mean(BodyWght, na.rm = T), Group]
Metadata[is.na(BodyWght), BodyWght:=GroupMean]

#ggplot(Metadata, aes(x=Group, y=DaysOld))+geom_boxplot()

Metadata[!is.na(Born),ImputationBDay := Born]
Metadata[is.na(Born) & Group == "AD" & !is.na(S1Date),ImputationBDay := as.Date.character(paste0("01/05/", year(S1Date)-1), "%d/%m/%Y")]

Metadata[!is.na(DaysOld), ImputedAgeDays:=DaysOld]
Metadata[is.na(DaysOld) & !is.na(S1Date) & !is.na(ImputationBDay), ImputedAgeDays:= as.integer(S1Date - ImputationBDay)]

#Metadata[Group=="AD" & ImputedAge < 365, ImputationBDay := as.Date.character(paste0("01/05/", year(S1Date)-2), "%d/%m/%Y")]
#Metadata[Group=="AD" & ImputedAge < 365 & !is.na(S1Date), ImputedAge:= (S1Date - ImputationBDay)]

Metadata[,GroupMean:=mean(DaysOld, na.rm = T), Group]
Metadata[is.na(ImputedAgeDays), ImputedAgeDays:=as.integer(GroupMean)]
Metadata[,GroupMean:=NULL]

#CharlotteTTData <- fread("./Metadata/170926_CERS_MeanLeftAtria-TTDensity-Development.csv")
#CharlotteTTData.PerGroup <- CharlotteTTData[,.(meanFA=mean(`F.A.`)), Age]

#cor.test(Metadata[,HeartWght], Metadata[,ImputedAgeDays]) #0.78!! strong positive correlation between HeartWeight and age (duh)
#Therefore, we combine them - allowing us to normalise for heart weight and test for age in a single parameter.
Metadata[, HWAgeRatio:=(HeartWght/log(ImputedAgeDays)) ]
#ggplot(Metadata, aes(x=ImputedAgeDays, y=HWAgeRatio))+geom_point()+geom_text(aes(label=sample), nudge_x=1, nudge_y=1)
#pretty linear,but 161208-06 (3M) is an outlier... We're already down to only 6 3M samples, though.

set(Metadata, j = Metadata[, which((colnames(Metadata) %in% c("S1Date", "ImputedBDay")))], value=NULL)

