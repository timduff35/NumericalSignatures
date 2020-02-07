library(data.table)
                                        # read data
# standard
D2 <- read.csv("../experiment_data/diff-standard-degree-2",skip=1)
D3 <- read.csv("../experiment_data/diff-standard-degree-3",skip=1)
D4 <- read.csv("../experiment_data/diff-standard-degree-4",skip=1)
D5 <- read.csv("../experiment_data/diff-standard-degree-5",skip=1)
D6 <- read.csv("../experiment_data/diff-standard-degree-6",skip=1)
D=rbindlist(list(D2,D3,D4,D5,D6),idcol=TRUE)
D[,1] = D[,1]+1
# multiprojective
d2 <- read.csv("../experiment_data/diff-1-degree-2",skip=1)
d3 <- read.csv("../experiment_data/diff-1-degree-3",skip=1)
d4 <- read.csv("../experiment_data/diff-1-degree-4",skip=1)
d5 <- read.csv("../experiment_data/diff-1-degree-5",skip=1)
d6 <- read.csv("../experiment_data/diff-1-degree-6",skip=1)
d=rbindlist(list(d2,d3,d4,d5,d6),idcol=TRUE)
d[,1] = as.integer(d[,1]+1)

# set plot parameters
COLS <- c("orange","red","purple","blue","green")
PCH <- 19
TYPE <- "b"

# standard plot
png("../PIX/diff-standard-sensitivity.png")
plot(D2$noise,100*D2$FN,log="x",col="orange",pch=PCH,xlab="Magnitude of noise",ylab="Negative test percentage",ylim=c(0,100),type=TYPE,main="Differential signature sensitivity")
legend(x="bottomright",lty=1,lwd=1,col=COLS,legend=c("2", "3","4", "5", "6"),title="degree")
points(D3$noise,100*D3$FN,col="red",pch=PCH,type=TYPE)
points(D4$noise,100*D4$FN,col="purple",pch=PCH,type=TYPE)
points(D5$noise,100*D5$FN,col="blue",pch=PCH,type=TYPE)
points(D6$noise,100*D6$FN,col="green",pch=PCH,type=TYPE)
dev.off()

# multiprojective plot
png("../PIX/diff-1-sensitivity.png")
plot(d2$noise,100*d2$FN,log="x",col="orange",pch=PCH,xlab="Magnitude of noise",ylab="Negative test percentage",ylim=c(0,100),type=TYPE,main="Differential signature sensitivity (multiprojective)")
legend(x="bottomright",lty=1,lwd=1,col=COLS,legend=c("2", "3","4", "5", "6"),title="degree")
points(d3$noise,100*d3$FN,col="red",pch=PCH,type=TYPE)
points(d4$noise,100*d4$FN,col="purple",pch=PCH,type=TYPE)
points(d5$noise,100*d5$FN,col="blue",pch=PCH,type=TYPE)
points(d6$noise,100*d6$FN,col="green",pch=PCH,type=TYPE)
dev.off()

# make table
library(xtable)
dAgg=aggregate(d[,c(3,4)],list(d$.id), mean)
DAgg=aggregate(D[,c(3,4)],list(D$.id), mean)
colnames(dAgg) <- c("degree", "track time (ms)", "lookup time (ms)")
print(xtable(dAgg),include.rownames=FALSE)
