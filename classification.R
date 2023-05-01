#classifying
Args <- commandArgs(TRUE)
library(randomForest)
model<-paste(Args[2],"/",Args[1],".rfo",sep="")
load(model)
df<-read.table("outputdataframe.csv", sep=",", header=T, row.names=1, comment.char = "")
featurematrix <-df
featurematrix$contigID <- NULL
featurematrix$genome <- NULL
#notice the missing featurevector generation and the use of predict
rf.classifications <- predict(rf,featurematrix)
rf.votes <- predict(rf,featurematrix, type="vote")
combined1 <- merge(rf.classifications,rf.votes,by="row.names")
row.names(combined1)<-combined1$Row.names
combined1$Row.names <-NULL
colnames(combined1) <- c("prediction", "votes chromosomal", "votes plasmid")
output <- merge(combined1,df,by="row.names")
row.names(output)<-output$Row.names
output$Row.names <-NULL
write.csv(output[1:4], file="prediction.csv")
write.csv(output, file="prediction_full.csv")
