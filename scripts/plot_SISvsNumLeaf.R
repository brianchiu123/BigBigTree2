fileName <- commandArgs();
input <- fileName[length(fileName)-1];
output <- fileName[length(fileName)];

r<-read.csv(input,header=TRUE)
jpeg(output)
plot(r)
dev.off()
q(save="no")
