#training
library(randomForest)
df<-read.table("outputdataframe.csv", sep=",", header=T, row.names=1)
featurevector <- substr(row.names(df), 1, 1)
featurematrix <-df
featurematrix$contigID <- NULL
rf <- randomForest(formula = y ~ x, x=featurematrix,y=as.factor(featurevector), ntree=5000, sampsize=(c(1000, 1000)))
save(rf,file = "training.rfo")
rf.classifications <- predict(rf,featurematrix)
rf.votes <- predict(rf,featurematrix, type="vote")

#combine votes and classifcation predict 100 percent
combined1 <- merge(rf.classifications,rf.votes,by="row.names")
row.names(combined1)<-combined1$Row.names
combined1$Row.names <-NULL

#combine votes and classication training 66/33 percent
combined2 <- merge(rf$predicted,rf$votes,by="row.names")
row.names(combined2)<-combined2$Row.names
combined2$Row.names <-NULL

#combine predict and training values
combined3 <- merge(combined2,combined1,by="row.names")
row.names(combined3)<-combined3$Row.names
combined3$Row.names <-NULL

colnames(combined3) <- c("classification training", "votes chromosomal training", "votes plasmid training","classification predict", "votes chromosomal predict", "votes plasmid predict")

output <- merge(combined3,df,by="row.names")
row.names(output)<-output$Row.names
output$Row.names <-NULL
write.csv(combined3, file="classified.csv")
write.csv(output, file="classified_full.csv")
