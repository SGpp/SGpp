data = read.table(file="Documents/data/DR5_MGS_fuat.csv", sep=",", head=T)
library(foreign)
normalize<- function(x) {
x = (x-min(x))/(max(x) - min(x))
return(x)
}
data$dered_u <- normalize(data$dered_u)
data$dered_g <- normalize(data$dered_g)
data$dered_r <- normalize(data$dered_r)
data$dered_i <- normalize(data$dered_i)
data$dered_z <- normalize(data$dered_z)
data$eClass <- normalize(data$eClass)
data_permuted = data[sample(1:length(data[,1])),]
data_permuted$classes = data_permuted$redshift
data_permuted$redshift = NULL
write.arff(data_permuted[1:60000,], file="DR5_MGS_ugrizeC_train.60000.arff")
write.arff(data_permuted[60001:120000,], file="DR5_MGS_ugrizeC_test.60000.arff")
write.arff(data_permuted, file="DR5_MGS_ugrizeC_all.arff")
