bupa = read.arff(file="../../../../liver-disorders_normalized.arff")
fold_level = 10
size = length(bupa[,1])
window = ceiling(size/fold_level)
for (i in 1:(fold_level-1)){
cat("i",i,"\n")
cat("frame:",((i-1)*window+1),"till",(i*window),"\n")
test_indexes = c(((i-1)*window+1):(i*window))
test_set <- bupa[test_indexes,]
cat("test indexes",test_indexes,"\n")
train_indexes = setdiff(c(1:size), test_indexes)
train_set <- bupa[train_indexes,]
write.arff(test_set, file=paste("../../../../liver_disorders_test.fold",i,".arff",sep=""))
write.arff(train_set, file=paste("../../../../liver_disorders_train.fold",i,".arff",sep=""))
}

test_indexes = ((fold_level-1)*window+1):size
test_set <- bupa[test_indexes,]
train_indexes = setdiff(c(1:size), test_indexes)
train_set <- bupa[train_indexes,]
write.arff(test_set, file=paste("../../../../liver_disorders_test.fold",fold_level,".arff",sep=""))
write.arff(train_set, file=paste("../../../../liver_disorders_train.fold",fold_level,".arff",sep=""))
