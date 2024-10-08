##Code for Fevrier and Barraclough 2024 analysis of dS and concatenated analysis
##Tim Barraclough tim.barraclough@biology.ox.ac.uk
##Code for looking at gene trees and alignments of orthogroups
	
##CONCATENATION APPROACH

##load packages
	library(ape)
	library(phangorn)
	library(seqinr)
	
##set up color-blind palette with 4 colours
	palette(c( "#E69F00","#D55E00","#56B4E9","#009E73"))
	
##paths to tree and sequence alignment inputs 
	input.path<-"trees_trimmed_alignments"

##extract species names and allocate order names
	species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
					"Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
					"Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
					"Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
	order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
	order.colours<-as.integer(as.factor(order.names))

##pull out the list of tree files in the input.path folder
	tree.files<-dir(input.path)[grep(".treefile",dir(input.path))]

	results<-read.csv("PierreGeneStatsConcat.csv")

##focus on the 2142 genes with all species present and 30 or fewer sequences
	focal.genes<-which((results$num_species==21)&(results$num_tips<=30))
	
	par(mfrow=c(2,5))
	
##results file to save statistics of orthogroup properties
concat.results<-NULL

concat.seqs<-NULL

	for (j in (focal.genes)) {
		
	gene.tr<-read.tree(paste0(input.path,"/",tree.files[j]))

##use midpoint rooting
	gene.tr<-midpoint(gene.tr)
	
##extract the names removing the part after "_" and match them with names above
	tip.names<-substring(gene.tr$tip.label ,1,regexpr("_",gene.tr$tip.label)-1)	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]

##nodes with all species descended
	all.species<-which(unlist(lapply(Descendants(gene.tr),function(x) length(unique(tip.names[x]))))==21)	
##if there are several, select the highest node number = most nested node with all species descended
	if (length(all.species)>1) {	
	new.root<-max(all.species)
	gene.tr<-extract.clade(gene.tr,new.root)
	tip.names<-substring(gene.tr$tip.label ,1,regexpr("_",gene.tr$tip.label)-1)	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]}

##if there are duplicates, remove them
	if (any(duplicated(duplicated(tip.names)))) {
	##remove the duplicates - keeps copy that is 'first' in the tree description
	gene.tr<-drop.tip(gene.tr,which(duplicated(tip.names)))	
	tip.names<-substring(gene.tr$tip.label ,1,regexpr("_",gene.tr$tip.label)-1)	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]}
			
#plot(gene.tr,tip.color=tip.colours)

	
	##running loop gives 95%ile of 1.9755736 for length of maximum branch, median is 0.7216353
	
	##Criterion to choose if will use tree, removing top 5% for maximum branch length
	if (max(gene.tr$edge.length)<=1.9755736) {

	##calculate kaks from alignments - need converting to fasta to read into seqinr
	seqs<-read.dna(file=paste0(input.path,"/",gsub(".treefile","",tree.files[j])),format="sequential")
	seqs<-seqs[pmatch(gene.tr$tip.label,rownames(seqs)),]		
	
	write.dna(seqs,file=paste0(input.path,"/",gsub("phy.treefile","fas",tree.files[j])),format="fasta",nbcol=-1,colsep="")
	seqs2<-read.alignment(file=paste0(input.path,"/",gsub("phy.treefile","fas",tree.files[j])),format="fasta")
	file.remove(file=paste0(input.path,"/",gsub("phy.treefile","fas",tree.files[j])))
	
	##ks and ka matrices
		ks.mat<-as.matrix(kaks(seqs2)$ks)
		ka.mat<-as.matrix(kaks(seqs2)$ka)
	##set diagonal to NA	
		diag(ks.mat)<-NA
		diag(ka.mat)<-NA
	##pull out the stats for each order
		within.ka<-array(NA,4)
		within.ks<-array(NA,4)
		GCcont<-array(NA,4)
		inter.mat<-matrix(TRUE,nrow=nrow(ka.mat),ncol=ncol(ka.mat))
		for (i in (1:length(unique(order.names)))) {
			if (any(tip.orders==unique(order.names)[i])) {
				within.ks[i]<-mean(ks.mat[tip.orders==unique(order.names)[i],tip.orders==unique(order.names)[i]],na.rm=T)
				within.ka[i]<-mean(ka.mat[tip.orders==unique(order.names)[i],tip.orders==unique(order.names)[i]],na.rm=T)
				GCcont[i]<-GC.content(seqs[tip.orders==unique(order.names)[i],])
				inter.mat[tip.orders==unique(order.names)[i],tip.orders==unique(order.names)[i]]<-FALSE
			}}
	
	##recheck which ones are monophyletic; 1=yes, 0=no, NA=absent
		check.monophyly<-NULL
	for (i in (1:length(unique(tip.orders)))) {
		which.tips<-which(tip.orders==unique(tip.orders)[i])
		check.monophyly<-c(check.monophyly,is.monophyletic(gene.tr,tips=which.tips))
	}
	names(check.monophyly)<-unique(tip.orders)
		
	##maximum pairwise distance
		max.pairwise<-max(cophenetic(gene.tr))
	##read into results file	
		concat.results<-rbind(concat.results,c(orthogroup=j,check.monophyly[pmatch(unique(order.names),names(check.monophyly))],
		within.ka,mean(ka.mat[!inter.mat],na.rm=T),within.ks,mean(ks.mat[!inter.mat],na.rm=T),
		mean(ka.mat[inter.mat],na.rm=T),
		mean(ks.mat[inter.mat],na.rm=T),GCcont,ncol(seqs),max.edge=max(gene.tr$edge.length),max.pairwise=max.pairwise))

	##concatenate sequences
	seqs<-seqs[pmatch(gene.tr$tip.label,rownames(seqs)),]
	seqs<-seqs[order(rownames(seqs)),]
	rownames(seqs)<-substring(rownames(seqs),1,regexpr("_",rownames(seqs))-1)
	if (j==focal.genes[1]) {concat.seqs<-seqs} else {concat.seqs<-cbind(concat.seqs,seqs)}

	}
	
	}
	
	concat.results<-data.frame(concat.results)
	rownames(concat.results)<-gsub("_dna.trim.phy.treefile","",tree.files[concat.results[,1]])
	colnames(concat.results)<-c("orthogroup",paste0(unique(order.names),"_mono"),
paste0(unique(order.names),"_ka"),"meanka_within",paste0(unique(order.names),"_ks"),"meanks_within",
"ka_between","ks_between",paste0(unique(order.names),"_gc"),"lseq","max.edge","max.pairwise")
	
	write.csv(concat.results,file="PierreConcat2034orthologs.csv")

	write.dna(concat.seqs,file="PierreConcat2034orthologs.phy",format="sequential",nbcol=-1,colsep="")
	
##summarise GC

GC.dat<-array(NA,nrow(concat.seqs))

for (i in (1:nrow(concat.seqs))) {
	GC.dat[i]<-GC.content(concat.seqs[i,])}
	
names(GC.dat)<-rownames(concat.seqs)
GC.order<-order.names[order(species.names)]

m<-lm(GC.dat~GC.order)
anova(m)

barplot(GC.dat[order(GC.order)],las=2,col=as.factor(GC.order[order(GC.order)]))
boxplot(GC.dat~GC.order)

##run IQTree

##./iqtree2 -s PierreConcat2038orthologs.phy -pre PierreConcat2038orthologs -alrt 0 -quiet -m TEST -fast -nt AUTO


##Make labelled tree

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

sp.tr<-read.tree("/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/IQTREE/PierreConcat2034orthologs.tre")
sp.tr<-ladderize(sp.tr,right=F)
root.branch<-which(sp.tr$edge[,1]==length(sp.tr$tip.label)+1)
root.branch.length<-sp.tr$edge.length[root.branch]
sp.tr$edge.length[root.branch]<-c(sum(root.branch.length)-0.08,0.08)
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

	write.tree(sp.tr,file="/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/IQTREE/PierreConcat2034orthologs.labelled.tre")
	
	
###LOOK AT RESULTS OUTPUT FROM PAML

##load packages
	library(ape)
	library(phangorn)
	library(seqinr)

##set up color-blind palette with 4 colours
	palette(c( "#E69F00","#D55E00","#56B4E9","#009E73"))
	
##paths to tree and sequence alignment inputs 
	input.path<-"/Users/tgb/Documents/Projects/Insects/pierre/trees_trimmed_alignments"

##extract species names and allocate order names
	species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
					"Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
					"Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
					"Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
	order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
	order.colours<-as.integer(as.factor(order.names))

out.path<-"/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/outputs"
out.files<-dir(out.path)

empirical.path<-"/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/empiricaloutputs"
empirical.files<-dir(empirical.path)

##set up results tables for the models
	model.results<-matrix(NA,nrow=length(out.files)+length(empirical.files),ncol=11)
	colnames(model.results)<-c("p0","p1","p2a","p2b","w0","w1","w2","LnL","AIC","MaxPairwise","MaxPairwiseDs")
	rownames(model.results)<-c(out.files,empirical.files)

##set up tree list
	trlist<-list()
	
	##Model 1
	j=grep("model1",out.files)
		
	x<-scan(file=paste0(out.path,"/",out.files[j]), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[j,9]=(-2*as.numeric(model.results[j,8])+2*2)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-grep("p:",substring(x,1,2))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[j,1:2]=as.numeric(split[2:3])
	#corresponding w
		line<-grep("w:   ",x)
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[j,5:6]=as.numeric(split[2:3])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[j]]<-read.tree(text=x[line+2])

		model.results[j,10]=max(cophenetic(trlist[[j]]))
		line<-grep("Time used",x)
		split=strsplit(x[line-1]," " )[[1]]
		split<-split[!split==""]
		model.results[j,11]=model.results[j,10]/(1+as.numeric(split[5]))


	
	##Model 2
	j=grep("model2",out.files)
		
	x<-scan(file=paste0(out.path,"/",out.files[j]), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[j,9]=(-2*as.numeric(model.results[j,8])+2*4)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-grep("p:",substring(x,1,2))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[j,1:3]=as.numeric(split[2:4])
	#corresponding w
		line<-(grep("w:   ",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[j,5:7]=as.numeric(split[2:4])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[j]]<-read.tree(text=x[line+2])

	
	
	
		##Model 3
	j=grep("model3",out.files)
		
	x<-scan(file=paste0(out.path,"/",out.files[j]), what=character(),sep="\n")
			
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[j,9]=(-2*as.numeric(model.results[j,8])+2*3)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-(grep("proportion",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[j,1:4]=as.numeric(split[2:5])
	#w
		line<-(grep("foreground",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[j,5:7]=as.numeric(split[3:5])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[j]]<-read.tree(text=x[line+2])




		##Model 4
	j=grep("model4",out.files)
		
	x<-scan(file=paste0(out.path,"/",out.files[j]), what=character(),sep="\n")
		
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[j,9]=(-2*as.numeric(model.results[j,8])+2*4)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-(grep("proportion",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[j,1:4]=as.numeric(split[2:5])
	#w
		line<-(grep("foreground",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[j,5:7]=as.numeric(split[3:5])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[j]]<-read.tree(text=x[line+2])
	
	##empirical

		##Model 2
	
	j=grep("model2",empirical.files)
		
	x<-scan(file=paste0(empirical.path,"/",empirical.files[j]), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[j+4,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[j+4,9]=(-2*as.numeric(model.results[j+4,8])+2*4)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-grep("p:",substring(x,1,2))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[j+4,1:3]=as.numeric(split[2:4])
	#corresponding w
		line<-(grep("w:   ",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[j+4,5:7]=as.numeric(split[2:4])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[j+4]]<-read.tree(text=x[line+2])
		
	j=grep("model4",empirical.files)
		
	x<-scan(file=paste0(empirical.path,"/",empirical.files[j]), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[j+4,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[j+4,9]=(-2*as.numeric(model.results[j+4,8])+2*4)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-(grep("proportion",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[j+4,1:4]=as.numeric(split[2:5])
	#w
		line<-(grep("foreground",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[j+4,5:7]=as.numeric(split[3:5])
	##pull out tree
		line<-grep("tree length",x)
		trlist[[j+4]]<-read.tree(text=x[line+2])

##NB ran model 2 and 4, but these collapse back to model 1 and 3 (w2=1), hence report those values and AIC accordingly
## Didn't run model 1 and 3 separately due to long run times

write.csv(model.results,"ConcatenatedModelResults.csv")

##AIC weights
	exp( -0.5 * (model.results[2:5,9]-min(model.results[2:5,9])))/sum(exp( -0.5 * (model.results[2:5,9]-min(model.results[2:5,9]))))
##exclude informative models
	to.include<-c(1,3)
	exp( -0.5 * (model.results[to.include,9]-min(model.results[to.include,9])))/sum(exp( -0.5 * (model.results[to.include,9]-min(model.results[to.include,9]))))
##and for empirical codon freq
	to.include<-c(5,6)
	exp( -0.5 * (model.results[to.include,9]-min(model.results[to.include,9])))/sum(exp( -0.5 * (model.results[to.include,9]-min(model.results[to.include,9]))))
	

###########
TRY RUNNING MODEL 3 WITH MODEL 1 TREE

##plot trees

par(mfrow=c(1,2))

plot(trlist[[1]])

j<-3
   tip.names<-trlist[[j]]$tip.label
  tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
  tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]

	within.order<-NULL
    for (k in (1:4)) {
      mrca.node<-getMRCA(trlist[[j]],tip=which(tip.colours==k))
      desc<-Descendants(trlist[[j]],mrca.node,type="all")
      within.order<-c(within.order,mrca.node,desc[desc>length(trlist[[j]]$tip.label)])}

edge.cols<-as.integer(trlist[[j]]$edge[,1]%in%within.order)+1

plot(trlist[[j]],edge.width=edge.cols)


##LOOK AT GC
library(ape)

##extract species names and allocate order names
	species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
					"Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
					"Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
					"Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
	order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
	order.colours<-as.integer(as.factor(order.names))
	
	palette(c( "#E69F00","#D55E00","#56B4E9","#009E73"))

seqs<-read.dna(file="/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/PierreConcat2034orthologs.phy",format="sequential")
 
seqs<-seqs[pmatch(species.names ,rownames(seqs)),] 
 
gc<-matrix(NA,nrow(seqs),5) 
 colnames(gc)<-c("First","Second","Third","All","GCaminos")

gc.long

for (i in (1:nrow(seqs))) {
	gc[i,1]<-GC.content(seqs[i,seq(1,ncol(seqs),3)]) 
	gc[i,2]<-GC.content(seqs[i,seq(2,ncol(seqs),3)]) 
	gc[i,3]<-GC.content(seqs[i,seq(3,ncol(seqs),3)]) 
	gc[i,4]<-GC.content(seqs[i,]) 
	GCfirst<-as.character(seqs[i,seq(1,ncol(seqs),3)])%in%c("g","c")
	GCsecond<-as.character(seqs[i,seq(2,ncol(seqs),3)])%in%c("g","c")
	gc[i,5]<-sum(GCfirst&GCsecond)/length(GCfirst)
	}


b<-barplot(gc[,1:4],col=as.factor(order.names),las=1,beside=T,ylab="GC %")
m<-lm(gc[,4]~order.names)
 
 
gc.mean<-apply(gc,2,function(x) by(x,INDICES=as.factor(order.names),FUN=mean))
gc.se<-apply(gc,2,function(x) by(x,INDICES=as.factor(order.names),FUN=sd))/sqrt(apply(gc,2,function(x) by(x,INDICES=as.factor(order.names),FUN=length)))

colMeans(gc)
apply(gc,2,function(x) sd(x)/sqrt(length(x)))


##look at the BEB codons

	x<-scan(file="/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/outputs/PierreConcat2034orthologs.model4.out", skip=33000,what=character(),sep="\n")
	abovey<-grep("Positive sites for foreground",x)
	endy<-grep("The grid",x)
	sig.codons<-grep("\\*",x)
	sig.codons<-sig.codons[sig.codons>abovey]
	sig.codons<-sig.codons[sig.codons<endy]
	tmp<-strsplit(x[sig.codons]," ")
	tmp<-lapply(tmp,function(x) x[!x==""][1])
	sig.codons<-unlist(tmp)
	sig.codons<-as.integer(sig.codons)
	
	seqsbeb<-seqs[,sort(c(3*sig.codons-2,3*sig.codons-1,3*sig.codons))]
	
gcbeb<-matrix(NA,nrow(seqs),5) 
 colnames(gcbeb)<-c("First","Second","Third","All","GCaminos")


for (i in (1:nrow(seqs))) {
	gcbeb[i,1]<-GC.content(seqsbeb[i,seq(1,ncol(seqsbeb),3)]) 
	gcbeb[i,2]<-GC.content(seqsbeb[i,seq(2,ncol(seqsbeb),3)]) 
	gcbeb[i,3]<-GC.content(seqsbeb[i,seq(3,ncol(seqsbeb),3)]) 
	gcbeb[i,4]<-GC.content(seqsbeb[i,]) 
	GCfirst<-as.character(seqsbeb[i,seq(1,ncol(seqsbeb),3)])%in%c("g","c")
	GCsecond<-as.character(seqsbeb[i,seq(2,ncol(seqsbeb),3)])%in%c("g","c")
	gcbeb[i,5]<-sum(GCfirst&GCsecond)/length(GCfirst)
	}	
	
b<-barplot(gcbeb[,1:4],col=as.factor(order.names),las=1,beside=T,ylab="GC %")

gcbeb.mean<-apply(gcbeb,2,function(x) by(x,INDICES=as.factor(order.names),FUN=mean))
gcbeb.se<-apply(gcbeb,2,function(x) by(x,INDICES=as.factor(order.names),FUN=sd))/sqrt(apply(gcbeb,2,function(x) by(x,INDICES=as.factor(order.names),FUN=length)))

colMeans(gcbeb)
apply(gcbeb,2,function(x) sd(x)/sqrt(length(x)))


	
          CodonFreq = 7  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
                         * 4:F1x4MG, 5:F3x4MG, 6:FMutSel0, 7:FMutSel
            estFreq = 0  * 0: use observed freqs; 1: estimate freqs by ML
 
 
 ##EXPORT 2nd and 3rd position sites
 
	GCsecond<-seqs[,seq(2,ncol(seqs),3)]
	GCthird<-seqs[,seq(3,ncol(seqs),3)]
     write.dna(GCsecond,file="PierreConcat2034orthologs.2nd.phy",nbcol=-1,colsep="",format="sequential")
      write.dna(GCthird,file="PierreConcat2034orthologs.3rd.phy",nbcol=-1,colsep="",format="sequential")
    
  cd /Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/second_third_trees
  ./iqtree2 -s PierreConcat2034orthologs.2nd.phy -pre PierreConcat2034orthologs.2nd -alrt 0 -quiet -m TEST -fast -nt AUTO         
  ./iqtree2 -s PierreConcat2034orthologs.3rd.phy -pre PierreConcat2034orthologs.3rd -alrt 0 -quiet -m TEST -fast -nt AUTO         
 
   ./iqtree2 -s PierreConcat2034orthologs.2nd.phy -t PierreConcat2034orthologs.tre -pre PierreConcat2034orthologs.2nd.fix -alrt 0 -quiet -m TEST -fast -nt AUTO         
   ./iqtree2 -s PierreConcat2034orthologs.3rd.phy -t PierreConcat2034orthologs.tre -pre PierreConcat2034orthologs.3rd.fix -alrt 0 -quiet -m TEST -fast -nt AUTO         

   
 ##read in tree fix 2nd and 3rd
 
 	sectr<-read.tree("/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/second_third_trees/2ndTrFix.tre")
 	thitr<-read.tree("/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/second_third_trees/3rdTrFix.tre")
 	secthitr<-sectr
	secthitr$edge.length<-sectr$edge.length/thitr$edge.length

 par(mfrow=c(1,3))
  	plot(sectr,main="A) 2nd positions",cex=1.3)
  	add.scale.bar(x=0,y=20,cex=1.3)
  	plot(thitr,main="B) 3rd positions",cex=1.3)
  	add.scale.bar(x=0,y=20,cex=1.3)
 	plot(secthitr,main="C) Ratio 2nd/3rd",cex=1.3)
  	add.scale.bar(x=0,y=20,cex=1.3)
 
 gc.tree<-gc[pmatch(tip.names ,rownames(seqs)),] 

layout(matrix(c(1,2,3,4,5),nrow=1),widths=c(4,1,4,1,4))
  	par(mar=c(5,4,3,0))
 	plot(sectr,main="A) 2nd positions",cex=1.3)
  	add.scale.bar(x=0,y=20,cex=1.3)
  	par(mar=c(5,0,3,2))
  	dotchart(gc.tree[,2],main="B) GC 2nd",xlab="GC%")
  	par(mar=c(5,4,3,0))
  	plot(thitr,main="C) 3rd positions",cex=1.3)
  	add.scale.bar(x=0,y=20,cex=1.3)
  	par(mar=c(5,0,3,2))
	dotchart(gc.tree[,3],main="D) GC 3rd",,xlab="GC%")
 	  par(mar=c(5,4,3,0))
	plot(secthitr,main="E) Ratio 2nd/3rd",cex=1.3)
	
  		
 	edge.table<-cbind(sectr$edge.length,thitr$edge.length,secthitr$edge.length)
 tip.names<-sectr$tip.label
  tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
  tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]

	within.order<-NULL
    for (k in (1:4)) {
      mrca.node<-getMRCA(sectr,tip=which(tip.colours==k))
      desc<-Descendants(sectr,mrca.node,type="all")
      within.order<-c(within.order,mrca.node,desc[desc>length(sectr$tip.label)])}

edge.cols<-as.integer(sectr$edge[,1]%in%within.order)+1

 	
 	m1<-lm(edge.table[,1]~as.factor(edge.cols))
 	m2<-lm(edge.table[,2]~as.factor(edge.cols))
 	m3<-lm(edge.table[,3]~as.factor(edge.cols))
 	anova(m1)
 	anova(m2)
 	anova(m3)
 	
 	m1<-lm(sectr$edge.length ~thitr$edge.length,subset=as.factor(edge.cols)==1)
 	m2<-lm(sectr$edge.length ~thitr$edge.length,subset=as.factor(edge.cols)==2)
	plot(thitr$edge.length, sectr$edge.length,col=edge.cols)
	abline(m1$coeff,col="red")
	abline(m2$coeff,col="black") 	
            
 ##PLOT IQTREE SPECIES TREE
 library(ape)
 tr<-read.tree("/Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/concatenated/IQTREE/PierreConcat2034orthologs.labelled.tre")

##between nodes
 bet.nodes<-which(tr$node.label=="#1")+ length(tr$tip.label)
 
 bcol<-array("black",nrow(tr$edge))
 bcol[which(tr$edge[,2]%in%bet.nodes)]<-"grey"
plot(tr,show.tip.label=0,edge.color=bcol)
add.scale.bar(y=21,x=0)


##LOOK AT AA COMPOSITION

aacomps<-apply(as.character(aaseqs),1,table)
aacomps<-t(aacomps)
ordercomps<-apply(aacomps,2,function(x) by(x,order.names,mean))
barplot(ordercomps,beside=T)

essential<-c("R","H","I","L","K","M","F","T","W","V")
GC2nd<-c("S","P","T","A","C","W","R","G")

which.essential<-colnames(ordercomps)%in%essential
which.GC<-colnames(ordercomps)%in%GC2nd

pvals<-apply(aacomps,2,function(x) anova(lm(x~order.names))[1,5])

pvals.plot<-pvals[order(!which.GC)]

##plot with GC 2nd position AAs first

bp<-barplot(ordercomps[,order(!which.GC)]/rowSums(ordercomps[,order(!which.GC)]),beside=T,ylim=c(0,0.12),ylab="Proportion",xlab="Amino acid")
mtext(at=colMeans(bp),side=1,line=0.4,text=ifelse(pvals.plot<0.001,"**",""),cex=1.3)
mtext(at=colMeans(bp),side=1,line=0.4,text=ifelse(pvals.plot<0.01,"*",""),cex=1.3)

aaseqsbeb<-trans(seqsbeb)
bebcomps<-apply(as.character(aaseqsbeb),1,table)
bebcomps<-t(bebcomps)
bebordercomps<-apply(bebcomps,2,function(x) by(x,order.names,mean))
barplot(bebordercomps,beside=T)

bebpvals<-apply(bebcomps,2,function(x) anova(lm(x~order.names))[1,5])

##2nd position T
phenylalanine,leucine,isoleucine,methionine, valine
11111

##2nd position C
serine, proline, threonine, alanine
0010

##2ndposition A
tyrosine, histidine, glutamine, asparagine, lysine, aspartate, glutamate
0100100

##2nd position G
cysteine, tryptophan, arginine, serine, glycine
011(0)0

##BEES ESSENTIAL AMINO ACIDS: 
arginine, histidine, isoleucine, leucine, lysine, methionine, phenylalanine, threonine, tryptophan and valine
GATTATTCGT

##ESSENTIAL FOR INSECTS:
arginine, histidine, isoleucine, leucine, lysine, methionine, phenylalanine, threonine, tryptophan and valine

R,H,I,L,K,M,F,T,W,V

##NON ESSENTIAL FOR INSECTS (NOT EXCLUSIVE - https://www.pnas.org/doi/full/10.1073/pnas.072346699
alanine (Ala), glycine (Gly), serine (Ser), proline (Pro), aspartate (Asp), and glutamate (Glu)
CGCCAA

##missing non-essential
tyrosine, glutamine,asparaginine, cysteine
AAAG

##ESSENTIAL FOR INSECTS
threonine (Thr), valine (Val) leucine (Leu), isoleucine (Ile), phenylalanine (Phe), and lysine (Lys) 
CTTTTA