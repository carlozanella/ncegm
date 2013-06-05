#generate raw data
system("../bin/ncegm > ../R/raw")

data <- read.fwf("raw",skip=9,widths=rep(26,5),col.names = c("a","V","c","aprime","d"))

plot(data$a,data$V,type="l")
