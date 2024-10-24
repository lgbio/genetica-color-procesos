
setClass("pop", slots=c(theta="matrix", r="matrix"))

setClass("bipop", slots=c(geno="matrix", info="data.frame", parent1="character", parent2="character"), contains="pop")

setClass("anypop", slots=c(geno="matrix", concordance="numeric", pops="list"), contains="pop")

