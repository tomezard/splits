`remove.terminal.zeros` <-
function(tr) {
num.tip<-length(tr$tip.label)
while (sum((tr$edge[,2]<=num.tip)&(tr$edge.length==0))>0) {
terminal.zeros<-which((tr$edge[,2]<=num.tip)&(tr$edge.length==0))
tr<-drop.tip(tr,tr$edge[terminal.zeros[!duplicated(tr$edge[terminal.zeros,1])],2])
num.tip<-length(tr$tip.label)}
return (tr)
}

