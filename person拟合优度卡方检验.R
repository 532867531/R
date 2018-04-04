X=scan()
y<-read.table("clipboard",header=F,sep="\t")
hist(as.numeric(y[,-c(1,2)]))

A=table(cut(X,br=c(0,69,79,89,100)))
p=pnorm(c(70,80,90,100,120),mean(X),sd(X))
