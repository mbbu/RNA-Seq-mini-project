#Comparison of tools
#percentage allignment comparison
library(readr)
aln <-read.csv("~/Downloads/percentage_allignment.csv",header = TRUE, sep = ",")
View(aln)
head(aln)
summary(aln)
library(ggplot2)
ggplot(aln,aes(sample,Percentage_allignment,fill=tools))+geom_bar(stat = "identity", position = "dodge")+labs(title = "Percentage allignment comparison")


#percentage reads after trimming
library(readr)
reads <-read.csv("~/Downloads/reads_after_trimming2.csv",header = TRUE, sep = ",")
View(reads)
library(ggplot2)
ggplot(reads,aes(sample,percentage_after_trimming,fill=tools))+geom_bar(stat = "identity", position = "dodge")+labs(title = "Percentage reads after trimming")


