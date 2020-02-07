png("../PIX/joint-3pan-sensitivity.png")
COLS <- c("orange","red","purple","blue","green")
DEGS <- c(2,3,4,5)
PCH <- 19
TYPE <- "b"
d2 <- read.csv("../experiment_data/joint-3pan-degree-2",skip=1)
plot(d2$noise,100*d2$FN,log="x",col="orange",pch=PCH,xlab="Magnitude of noise",ylab="Negative test percentage",ylim=c(0,100),type=TYPE,main="Joint signature sensitivity (multiprojective)")
legend(x="bottomright",lty=1,lwd=1,col=COLS,legend=c("2", "3","4", "5", "6"),title="degree")
d3 <- read.csv("../experiment_data/joint-3pan-degree-3",skip=1)
points(d3$noise,100*d3$FN,col="red",pch=PCH,type=TYPE)
d4 <- read.csv("../experiment_data/joint-3pan-degree-4",skip=1)
points(d4$noise,100*d4$FN,col="purple",pch=PCH,type=TYPE)
d5 <- read.csv("../experiment_data/joint-3pan-degree-5",skip=1)
points(d5$noise,100*d5$FN,col="blue",pch=PCH,type=TYPE)
d6 <- read.csv("../experiment_data/joint-3pan-degree-6",skip=1)
points(d6$noise,100*d6$FN,col="green",pch=PCH,type=TYPE)
dev.off()

mean(d2$track)
mean(d2$lookup)

mean(d3$track)
mean(d3$lookup)

mean(d4$track)
mean(d4$lookup)

mean(d5$track)
mean(d5$lookup)

mean(d6$track)
mean(d6$lookup)
