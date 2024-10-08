###Code for Fevrier and Barraclough 2024 analysis of dS and concatenated analysis
##Tim Barraclough tim.barraclough@biology.ox.ac.uk
##Code for looking at gene trees and alignments of orthogroups

setwd("simulations")

library(ape)
library(phangorn)

out.files<-"PierreConcat.lower20.1.out"

i<-1
	
##read in output from site model (same across all branches), i.e. trim.out files
	x<-scan(out.files[i],what=character(),sep="\n")
	
##read in template for evolver control file
	y<-scan("MCcodonNSsites.template.dat",what=character(),sep="\n")
	
##set random number seed using current time
	tmp<-strsplit(date()," ")[[1]][4]
	tmp<-gsub(":","",tmp)
	y<-gsub("randnum",tmp,y)
	
##set number of sequences and codons
	tmp<-unlist(strsplit(x[1]," "))
	tmp<-tmp[!tmp==""]
	y<-gsub("numseqs",tmp[1],y)
	y<-gsub("numcodons",as.integer(tmp[2])/3,y)
	
##load in the right tree - hopefully this always works to find the tree in the output file
	tmp<-grep("Detailed",x)
	y[5]<-x[tmp-1]

##load in class p and ws	
	tmp<-which(substring(x,1,2)=="p:")
	tmp<-unlist(strsplit(x[tmp]," "))
	y[7]<-paste0(tmp[!tmp==""][2:3],collapse="  ")
	
	tmp<-which(substring(x,1,2)=="w:")
	tmp<-unlist(strsplit(x[tmp]," "))
	y[8]<-paste0(tmp[!tmp==""][2:3],collapse="  ")
	
##load kappa
	y<-gsub("kappa.param",unlist(strsplit(x[grep("kappa",x)]," "))[5],y)
	
##codon frequencies
	tmp<-grep("Codon frequencies under model",x)
	y[10:25]<-x[(tmp+1):(tmp+16)]

##write template file	
	write(y,file="MCcodonNSsites.dat")

		line<-grep("tree length",x)
		tr<-read.tree(text=x[line+2])

		maxdist<-max(cophenetic(tr))
		line<-grep("Time used",x)
		split=strsplit(x[line-1]," " )[[1]]
		split<-split[!split==""]
		maxdist/(1+as.numeric(split[5]))

	##observed = 4.813654
	
##run evolve
	system("./evolverNSsites")
	file.copy("mc.txt","lower20sim.2.phy")
	
	
##make labelled tree
library(ape)
library(phangorn)

##set up color-blind palette with 4 colours
palette(c( "#E69F00","#D55E00","#56B4E9","#009E73"))

##extract species names and allocate order names
species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
                 "Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
                 "Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
                 "Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
order.colours<-as.integer(as.factor(order.names))

sp.tr<-tr
sp.tr<-ladderize(sp.tr,right=F)
plot(sp.tr)

  ##extract the names removing the part after "_" and match them with names above
  tip.names<-sp.tr$tip.label
  tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
  tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]

  within.order<-NULL
    for (k in (1:4)) {
      mrca.node<-getMRCA(sp.tr,tip=which(tip.colours==k))
      desc<-Descendants(sp.tr,mrca.node,type="all")
      within.order<-c(within.order,desc[desc>length(sp.tr$tip.label)])}
    sp.tr$node.label<-rep("#1",sp.tr$Nnod)
    sp.tr$node.label[1]<-""	
    sp.tr$node.label[(within.order-length(sp.tr$tip.label))]<-"#0"
    sp.tr$tip.label<-paste0(sp.tr$tip.label, "#0")

	write.tree(sp.tr,file="PierreConcat.lower20.codon.model1.tre")


##Simulations manipulating  

##run scaling whole tree length - not used in paper in end
	
#	for (scale in (c(1.5,2,5))) {
	
#	y<-scan("MCcodonNSsites.dat",what=character(),sep="\n")
#	new.tr<-tr
#	new.tr$edge.length<-new.tr$edge.length*scale
#	y[5]<-write.tree(new.tr)
	
#	for (reps in 1:3) {
	
#	tmp<-strsplit(date()," ")[[1]][4]
#	tmp<-gsub(":","",tmp)
#	y[3]<-paste(c(21,tmp,1),collapse=" ")
	
#	write(y,file="MCcodonNSsites.dat")
#	system("./evolverNSsites")
#	file.copy("mc.txt",paste0("lower20sim.scale.",scale,".",reps,".phy"))

#	}}
	
##run scaling between order branches

	for (scale in (c(1.5,2,5))) {
	
	y<-scan("MCcodonNSsites.dat",what=character(),sep="\n")
	new.tr<-tr
	new.tr$edge.length[new.tr$edge[,2]%in%between.order]<-new.tr$edge.length[new.tr$edge[,2]%in%between.order]*scale
	y[5]<-write.tree(new.tr)
	
	for (reps in 1:3) {
	
	tmp<-strsplit(date()," ")[[1]][4]
	tmp<-gsub(":","",tmp)
	y[3]<-paste(c(21,tmp,1),collapse=" ")
	
	write(y,file="MCcodonNSsites.dat")
	system("./evolverNSsites")
	file.copy("mc.txt",paste0("lower20sim.orderscale.",scale,".",reps,".phy"))

	}}


##Use relative dates from Wiegmann et al http://www.biomedcentral.com/1741-7007/7/34
##Not used in the end
	#calib<-data.frame(cbind(node=between.order,age.min=c(350,169,300,228,280,195,111),age.max=c(350,169,300,228,280,195,111),soft.bounds=rep(NA,7)))
	#timetr<-chronos(tr, lambda = 1, model = "correlated", quiet = FALSE,
    #    calibration = calib)
  ##rescale branch lengths back to total original tree length
  	#timetr$edge.length<-timetr$edge.length*sum(tr$edge.length)/sum(timetr$edge.length)
	#y<-scan("MCcodonNSsites.dat",what=character(),sep="\n")
	#y[5]<-write.tree(timetr)
		
	#tmp<-strsplit(date()," ")[[1]][4]
	#tmp<-gsub(":","",tmp)
	#y[3]<-paste(c(21,tmp,1),collapse=" ")
	
	#write(y,file="MCcodonNSsites.dat")
	#system("./evolverNSsites")
	#file.copy("mc.txt",paste0("lower20sim.wiegmann.phy"))

	
	
##run iqtree

	cd /Users/user/Desktop/split.sims/sim.files
	for f in *.phy
	do
	./iqtree2 -s $f -alrt 0 -quiet -m GTR+I+G -fast
	done	
rm *.log
rm *.gz
rm *.bionj
rm *.iqtree
rm *.mldist


paml.path<-"split.sims"

setwd(paml.path)

file.names<-dir("./sim.files")
tree.files<-file.names[grep(".treefile",file.names)]
seq.files<-gsub(".treefile","",tree.files)

for (i in (1:length(seq.files))) {


##First substitute tmp in the template file with sequence names
##The template is set up to run the model with 2 branch classes
x<-scan("codeml.template.ctl",what=character(),sep="\n")
x<-gsub("tmp",seq.files[i],x)
write(x,file=paste0("sim.files/",gsub(".phy",".template.ctl",seq.files[i])))
#	system("./codeml")
	
##Next change the output file name and the model and rerun
##This is the model with a single dn/ds set for all branches

x<-scan("codeml.template.label.ctl",what=character(),sep="\n")
x<-gsub("tmp",seq.files[i],x)
write(x,file=paste0("./sim.files/",gsub(".phy",".template.label.ctl",seq.files[i])))

	
}



###LOOK AT RESULTS OUTPUT FROM PAML


##load packages
	library(ape)
	library(phangorn)
	library(seqinr)

##set up color-blind palette with 4 colours
	palette(c( "#E69F00","#D55E00","#56B4E9","#009E73"))
	
##extract species names and allocate order names
	species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
					"Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
					"Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
					"Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
	order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
	order.colours<-as.integer(as.factor(order.names))

###MODEL 1

out.path<-"simulation.out"
out.files<-dir(out.path)
out.files<-unique(substring(out.files,1,regexpr("model",out.files)-2))

##set up results tables for the models
	model.results<-matrix(NA,nrow=length(out.files)*2,ncol=11)
	colnames(model.results)<-c("p0","p1","p2a","p2b","w0","w1","w2","LnL","AIC","MaxPairwise","MaxPairwiseDs")
	rownames(model.results)<-paste(rep(out.files,each=2),c(".model1",".model3"),sep="")

##set up tree list
	trlist<-list()
	
for (j in (1:length(out.files))) {
		
##model 1
	x<-scan(file=paste0(out.path,"/",out.files[j],".model1.out"), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[2*j-1,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[2*j-1,9]=(-2*as.numeric(model.results[2*j-1,8])+2*2)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-grep("p:",substring(x,1,2))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[2*j-1,1:2]=as.numeric(split[2:3])
	#corresponding w
		line<-grep("w:   ",x)
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[2*j-1,5:6]=as.numeric(split[2:3])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[2*j-1]]<-read.tree(text=x[line+2])

		model.results[2*j-1,10]=max(cophenetic(trlist[[j]]))
		line<-grep("Time used",x)
		split=strsplit(x[line-1]," " )[[1]]
		split<-split[!split==""]
		model.results[2*j-1,11]=model.results[2*j-1,10]/(1+as.numeric(split[5]))

##model 3
	x<-scan(file=paste0(out.path,"/",out.files[j],".model3.out"), what=character(),sep="\n")
			
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[2*j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[2*j,9]=(-2*as.numeric(model.results[2*j,8])+2*3)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-(grep("proportion",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[2*j,1:4]=as.numeric(split[2:5])
	#w
		line<-(grep("foreground",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[2*j,5:7]=as.numeric(split[3:5])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[2*j]]<-read.tree(text=x[line+2])



}

write.csv(model.results,"Lower20Sims.PamlResults.csv")

par(mfrow=c(2,4),mar=c(2,2,2,1))
##plot figure of trees
for (i in c(1,7,13,19)) {
	p<-plot(trlist[[i]],cex=1.2) 
	add.scale.bar(x=0,y=20)
	dat<-c(model.results[i,1],ifelse(is.na(model.results[i,4]),model.results[i,2],model.results[i,2]+ model.results[i,4]),model.results[i,3])
	rect(rep(p$x.lim[2]*0.9,4),c(0,cumsum(dat[-length(dat)]))*6+1,rep(p$x.lim[2]*0.95,4),cumsum(dat)*6+1,col=c("lightblue","white","darkred"))}
	for (i in c(2,8,14,20)) {
	p<-plot(trlist[[i]],cex=1.2)
	add.scale.bar(x=0,y=20) 
		dat<-c(model.results[i,1],ifelse(is.na(model.results[i,4]),model.results[i,2],model.results[i,2]+ model.results[i,4]),model.results[i,3])
	rect(rep(p$x.lim[2]*0.9,4),c(0,cumsum(dat[-length(dat)]))*6+1,rep(p$x.lim[2]*0.95,4),cumsum(dat)*6+1,col=c("lightblue","white","darkred"))}
	
	
##AIC weights
	exp( -0.5 * (model.results[1:4,9]-min(model.results[1:4,9])))/sum(exp( -0.5 * (model.results[1:4,9]-min(model.results[1:4,9]))))
##exclude informative models
	exp( -0.5 * (model.results[c(1,3),9]-min(model.results[c(1,3),9])))/sum(exp( -0.5 * (model.results[c(1,3),9]-min(model.results[c(1,3),9]))))

model.results[c(1,3),9]-min(model.results[c(1,3),9]

	