# read data
D2 <- read.csv("../experiment_data/diff-standard-degree-2",skip=1)
D3 <- read.csv("../experiment_data/diff-standard-degree-3",skip=1)
D4 <- read.csv("../experiment_data/diff-standard-degree-4",skip=1)
D5 <- read.csv("../experiment_data/diff-standard-degree-5",skip=1)
D6 <- read.csv("../experiment_data/diff-standard-degree-6",skip=1)

# standard plot
png("../PIX/diff-standard-sensitivity.png")
plot(D2$noise,100*D2$FN,log="x",col="orange",pch=PCH,xlab="Magnitude of noise",ylab="Negative test percentage",ylim=c(0,100),type=TYPE,main="Differential signature sensitivity")
legend(x="bottomright",lty=1,lwd=1,col=COLS,legend=c("2", "3","4", "5", "6"),title="degree")
points(D3$noise,100*D3$FN,col="red",pch=PCH,type=TYPE)
points(D4$noise,100*D4$FN,col="purple",pch=PCH,type=TYPE)
points(D5$noise,100*D5$FN,col="blue",pch=PCH,type=TYPE)
points(D6$noise,100*D6$FN,col="green",pch=PCH,type=TYPE)
dev.off()

