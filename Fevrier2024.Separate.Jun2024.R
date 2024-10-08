##Code for Fevrier and Barraclough 2024 analysis of dS and concatenated analysis
##Tim Barraclough tim.barraclough@biology.ox.ac.uk
##Code for looking at gene trees and alignments of orthogroups

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
	
##function to find nodes defining a paralog set with all orders included

find.focal.nodes<-function() {
##calculate the number of orders descended from each node of the tree
	numb.orders.desc.node<-unlist(lapply(Descendants(gene.tr,type="tips"),function(x) length(unique(tip.orders[x]))))
	numb.orders.mat<-matrix(numb.orders.desc.node[gene.tr$edge],ncol=2)
##go through this and find nodes that have all 4 but both descendents have fewer - these define a potential 'unit' 
##for analysis within orthogroups with divergent paralogs
	focal.node<-NULL
	for (k in (sort(unique(gene.tr$edge[,1])))) {
		if((numb.orders.mat[gene.tr$edge[,1]==k,][1,1]==4)&!any(numb.orders.mat[gene.tr$edge[,1]==k,][,2]==4)) {focal.node<-c(focal.node,k)}}
	return(focal.node)}

##results file to save statistics of orthogroup properties
results<-matrix(NA,nrow=length(tree.files),ncol=25)

colnames(results)<-c(paste0(unique(order.names),"_mono"),"num_tips","num_orders","num_species",
paste0(unique(order.names),"_ka"),"meanka_within",paste0(unique(order.names),"_ks"),"meanks_within",
"ka_between","ks_between",paste0(unique(order.names),"_gc"),"lseq")

##loop through the trees
for (j in (1:length(tree.files))) {

##print progress to screen
if (round(j/100,0)==j/100) print(j)

##load the relevant gene tree
	gene.tr<-read.tree(paste0(input.path,"/",tree.files[j]))

##midpoint root
	gene.tr<-midpoint(gene.tr)
	
##extract the names removing the part after "_" and match them with names above
	tip.names<-substring(gene.tr$tip.label ,1,regexpr("_",gene.tr$tip.label)-1)	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]

##plot the gene tree colouring orders
#plot(gene.tr,tip.color=tip.colours)
	
	##check which ones are monophyletic; 1=yes, 0=no, NA=absent
		check.monophyly<-NULL
	for (i in (1:length(unique(tip.orders)))) {
		which.tips<-which(tip.orders==unique(tip.orders)[i])
		check.monophyly<-c(check.monophyly,is.monophyletic(gene.tr,tips=which.tips))
	}
	names(check.monophyly)<-unique(tip.orders)

	##record monophyly and some other stats	
	results[j,1:4]<-check.monophyly[pmatch(unique(order.names),names(check.monophyly))]	
	results[j,5]<-length(gene.tr$tip.label)	
	results[j,6]<-length(unique(tip.orders))
	results[j,7]<-length(unique(tip.names))
	
	##calculate kaks from alignments - need converting to fasta to read into seqinr
	seqs<-read.dna(file=paste0(input.path,"/",gsub(".treefile","",tree.files[j])),format="sequential")
	write.dna(seqs,file=paste0(input.path,"/",gsub("phy.treefile","fas",tree.files[j])),format="fasta",nbcol=-1,colsep="")
	seqs2<-read.alignment(file=paste0(input.path,"/",gsub("phy.treefile","fas",tree.files[j])),format="fasta")
	file.remove(file=paste0(input.path,"/",gsub("phy.treefile","fas",tree.files[j])))
	
	##ks and ka matrices
		ks.mat<-as.matrix(kaks(seqs2)$ks)
		ka.mat<-as.matrix(kaks(seqs2)$ka)
	##reorder to be same as gene tree tip labels	
		ks.mat<-ks.mat[pmatch(gene.tr$tip.label,rownames(ks.mat)),pmatch(gene.tr$tip.label,rownames(ks.mat))]
		ka.mat<-ka.mat[pmatch(gene.tr$tip.label,rownames(ka.mat)),pmatch(gene.tr$tip.label,rownames(ka.mat))]
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
			
	##read into results file	
		results[j,8:11]<-within.ka
		results[j,12]<-mean(ka.mat[!inter.mat],na.rm=T)
		results[j,13:16]<-within.ks
		results[j,17]<-mean(ks.mat[!inter.mat],na.rm=T)
		results[j,18]<-mean(ka.mat[inter.mat],na.rm=T)
		results[j,19]<-mean(ks.mat[inter.mat],na.rm=T)
		results[j,20:23]<-GCcont
		results[j,24]<-ncol(seqs)
		results[j,25]<-max(gene.tr$edge.length)
			
	}
	
	##convert to dataframe
	results<-data.frame(results)
	
	##write the results file
	write.csv(results,file="PierreGeneStatsConcat.csv")
	
##Numbers of orthologs meeting particular criteria	

##at least one representative of each order
	sum(results[,6]==4)

##all species present
	sum(results[,7]==21)

##just one sequence for each of all species present
	sum((results[,7]==21)&(results[,5]==21))

	
##SEPARATE ANALYSES APPROACH


##load packages
	library(ape)
	library(phangorn)
	library(seqinr)
	
##set up color-blind palette with 4 colours
	palette(c( "#E69F00","#D55E00","#56B4E9","#009E73"))
	
##paths to tree and sequence alignment inputs 
	input.path<-"trees_trimmed_alignments"

##output path for paml
	paml.path<-"paml_revisions"

##extract species names and allocate order names
	species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
					"Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
					"Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
					"Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
	order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
	order.colours<-as.integer(as.factor(order.names))
	orders<-c("Coleoptera","Diptera","Hymenoptera","Lepidoptera")

##pull out the list of tree files in the input.path folder
	tree.files<-dir(input.path)[grep(".treefile",dir(input.path))]
	
##function to find nodes defining a paralog set with all orders included

find.focal.nodes<-function() {
##calculate the number of orders descended from each node of the tree
	numb.orders.desc.node<-unlist(lapply(Descendants(gene.tr,type="tips"),function(x) length(unique(tip.orders[x]))))
	numb.orders.mat<-matrix(numb.orders.desc.node[gene.tr$edge],ncol=2)
##go through this and find nodes that have all 4 but both descendents have fewer - these define a potential 'unit' 
##for analysis within orthogroups with divergent paralogs
	focal.node<-NULL
	for (k in (sort(unique(gene.tr$edge[,1])))) {
		if((numb.orders.mat[gene.tr$edge[,1]==k,][1,1]==4)&!any(numb.orders.mat[gene.tr$edge[,1]==k,][,2]==4)) {focal.node<-c(focal.node,k)}}
	return(focal.node)}
	
results<-read.csv(file="PierreGeneStatsConcat.csv")

##focus on the genes with all 4 orders present
	focal.genes<-which((results$num_order==4))
	
	#par(mfrow=c(2,5))
	
##results file to save statistics of orthogroup properties
sep.results<-NULL

	for (j in (focal.genes)) {
		
	gene.tr<-read.tree(paste0(input.path,"/",tree.files[j]))

##use midpoint rooting
	gene.tr<-read.tree(text=write.tree(midpoint(gene.tr)))
##extract the names removing the part after "_" and match them with names above
	tip.names<-substring(gene.tr$tip.label ,1,regexpr("_",gene.tr$tip.label)-1)	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]

#plot(gene.tr,tip.color=tip.colours)

##nodes with all species descended
	all.orders<-which(unlist(lapply(Descendants(gene.tr),function(x) length(unique(tip.orders[x]))))==4)	
##if there are several, select the highest node number = most nested node with all species descended
	if (length(all.orders)>1) {	
	new.root<-max(all.orders)
	gene.tr<-extract.clade(gene.tr,new.root)
	tip.names<-substring(gene.tr$tip.label ,1,regexpr("_",gene.tr$tip.label)-1)	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]}

#plot(gene.tr,tip.color=tip.colours)

##if there are duplicates, remove them
	if (any(duplicated(duplicated(tip.names)))) {
	##remove the duplicates - keeps copy that is 'first' in the tree description
	gene.tr<-drop.tip(gene.tr,which(duplicated(tip.names)))	
	tip.names<-substring(gene.tr$tip.label ,1,regexpr("_",gene.tr$tip.label)-1)	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]}
			
#plot(gene.tr,tip.color=tip.colours)

	##running loop gives 95%ile of 1.766438 for length of maximum branch, median is 0.7216353
	
	##Criterion to choose if will use tree, removing top 5% for maximum branch length
	##and inner 95% for mean number of samples per order =  2.5% 97.5%, 2.75  8.25 

	if ((max(gene.tr$edge.length)<=1.766438)&(length(tip.orders)>=10)) {

	##calculate kaks from alignments - need converting to fasta to read into seqinr
		seqs<-read.dna(file=paste0(input.path,"/",gsub(".treefile","",tree.files[j])),format="sequential")
		seqs<-seqs[pmatch(gene.tr$tip.label,rownames(seqs)),]		
		
	##trim tip labels
		rownames(seqs)<-substring(rownames(seqs),1,regexpr("_transcript",rownames(seqs))-1)
		gene.tr$tip.label<-substring(gene.tr$tip.label,1,regexpr("_transcript",gene.tr$tip.label)-1)

	##label gene tree
	tr.labels<-unlist(lapply(Descendants(gene.tr,type="tips"),function(x) ifelse(length(unique(tip.orders[x]))==1,"#0","#1")))
	
	new.labels<-array(NA,length(tr.labels))
	new.labels[gene.tr$edge[,2]]<-tr.labels[gene.tr$edge[,1]]
	gene.tr$node.label<-new.labels[(length(gene.tr$tip.label)+1):length(new.labels)]
    gene.tr$tip.label<-paste0(gene.tr$tip.label, new.labels[1:(length(gene.tr$tip.label))])	
		
		rownames(seqs)<-paste0(rownames(seqs)," ")
		write.dna(seqs,file=paste0(paml.path,"/",gsub("phy.treefile","paml.phy",tree.files[j])),format="sequential",nbcol=-1,colsep="")
		write.tree(gene.tr,file=paste0(paml.path,"/",gsub("phy.treefile","paml.tre",tree.files[j])))
	
	##recheck which ones are monophyletic; 1=yes, 0=no, NA=absent
		check.monophyly<-NULL
		GCcont<-array(NA,4)	
	for (i in (1:length(orders))) {
		which.tips<-which(tip.orders==orders[i])
		check.monophyly<-c(check.monophyly,is.monophyletic(gene.tr,tips=which.tips))
		GCcont[i]<-GC.content(seqs[tip.orders==orders[i],])
	}
	names(check.monophyly)<-unique(tip.orders)
		
	##read into results file	
		sep.results<-rbind(sep.results,c(orthogroup=j,check.monophyly[pmatch(sort(unique(order.names)) ,names(check.monophyly))],
			lseq=ncol(seqs),max.edge=max(gene.tr$edge.length),table(tip.orders),GCcont))

	##concatenate sequences
		#seqs<-seqs[pmatch(gene.tr$tip.label,rownames(seqs)),]
		#seqs<-seqs[order(rownames(seqs)),]
		#rownames(seqs)<-substring(rownames(seqs),1,regexpr("_",rownames(seqs))-1)
		#if (j==focal.genes[1]) {concat.seqs<-seqs} else {concat.seqs<-cbind(concat.seqs,seqs)}

	}
	
	}
	
	##leaves 4851 genes
	
	sep.results<-data.frame(sep.results)
	rownames(sep.results)<-gsub("_dna.trim.phy.treefile","",tree.files[sep.results[,1]])
	colnames(sep.results)<-c("orthogroup",paste0(orders,"_mono")
	,"lseq","max.edge",paste0(orders,"_num"),paste0(orders,"_gc"))
	
	write.csv(sep.results,file="PierreConcat2034orthologs.Separate.csv")

	
##copy from ARC
rsync -avz magd5205@arc-login.arc.ox.ac.uk:/data/zool-barralab/magd5205/PierrePaml/paml_revisions/ /Users/tgb/Documents/Projects/Insects/pierre/RevisionAnalyses/separate/paml_revisions_out


#####################
##READ IN OUT FILES##
#####################

library(ape)

out.path<-"paml_revisions_out"

out.files<-dir(out.path)

sep.results<-read.csv(file="PierreConcat2034orthologs.Separate.csv")

count.outfiles<-table(substring(out.files,1, regexpr("_",out.files)-1))
count.outfiles[which(count.outfiles<4)]

##these are missing 1 or more model output files = use them to redo
redo.orthos<-names(count.outfiles[which(count.outfiles<4)])

##these are missing altogether
sep.results[which(!sep.results[,1]%in%unique(substring(out.files,1, regexpr("_",out.files)-1))),1]

paste(sep.results[which(!sep.results[,1]%in%unique(substring(out.files,1, regexpr("_",out.files)-1))),1],collapse=" ")

##check whether any more have missing final bit in model 4 output file

model4.files<-out.files[grep("model4",out.files)]

for (i in (1:length(model4.files))) {
x<-scan(paste0(out.path,"/",model4.files[i]),what=character(),sep="\n")
if (any(grepl("lnL",x))) {} else {tmp<-c(tmp,substring(model4.files[i],1, regexpr("_",model4.files[i])-1))}}

##read in the results

model.results<-list()

tr1<-list()
tr2<-list()
tr3<-list()
tr4<-list()

paml.tree.stats<-matrix(NA,nrow=nrow(sep.results),ncol=2)
colnames(paml.tree.stats)<-c("MaxPairwise","MaxPairwiseDs")

##set up results tables for the four models
for (i in (1:4)) {
	model.results[[i]]<-matrix(NA,nrow=nrow(sep.results),ncol=9)
	colnames(model.results[[i]])<-c("p0","p1","p2a","p2b","w0","w1","w2","LnL","AIC")
	rownames(model.results[[i]])<-sep.results[,1]
	}

##loop through each of the OGs
for (j in (1:nrow(sep.results))) {
	
	##Model 1
	i=1
	
	if (file.exists(paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"))) {
	
	x<-scan(file=paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[[i]][j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[[i]][j,9]=(-2*as.numeric(model.results[[i]][j,8])+2*2)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-grep("p:",substring(x,1,2))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[[i]][j,1:2]=as.numeric(split[2:3])
	#corresponding w
		line<-grep("w:   ",x)
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[[i]][j,5:6]=as.numeric(split[2:3])
	##pull out tree
		line<-grep("tree length",x)
		tr1[[j]]<-read.tree(text=x[line+2])

		paml.tree.stats[j,1]=max(cophenetic(tr1[[j]]))
		line<-grep("Time used",x)
		split=strsplit(x[line-1]," " )[[1]]
		split<-split[!split==""]
		paml.tree.stats[j,2]=paml.tree.stats[j,1]/(1+as.numeric(split[5]))


	}  #end if
	
	##Model 2
	i=2
	
	if (file.exists(paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"))) {
	
	x<-scan(file=paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[[i]][j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[[i]][j,9]=(-2*as.numeric(model.results[[i]][j,8])+2*4)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-grep("p:",substring(x,1,2))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[[i]][j,1:3]=as.numeric(split[2:4])
	#corresponding w
		line<-(grep("w:   ",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""]
		model.results[[i]][j,5:7]=as.numeric(split[2:4])
	##pull out tree
		line<-grep("tree length",x)
		tr2[[j]]<-read.tree(text=x[line+2])

	}  #end if
	
		##Model 3
	i=3
	
	if (file.exists(paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"))) {
	
	x<-scan(file=paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[[i]][j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[[i]][j,9]=(-2*as.numeric(model.results[[i]][j,8])+2*3)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-(grep("proportion",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[[i]][j,1:4]=as.numeric(split[2:5])
	#w
		line<-(grep("foreground",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[[i]][j,5:7]=as.numeric(split[3:5])
	##pull out tree
		line<-grep("tree length",x)
		tr3[[j]]<-read.tree(text=x[line+2])

	}  #end if

		##Model 4
	i=4
	
	if (file.exists(paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"))) {
	
	x<-scan(file=paste0(out.path,"/",sep.results[j,1],"_dna.trim.model",i,".out"), what=character(),sep="\n")
	
	##find the lnL line, split it and extract the lnL value only then put it in the table
		line<-(grep("lnL",x))	
		split=unlist(strsplit(x[line]," " ))
		split<-split[!split==""] 
		model.results[[i]][j,8]=as.numeric(split[grep("-",split)])
	##calculate AIC
		model.results[[i]][j,9]=(-2*as.numeric(model.results[[i]][j,8])+2*4)
	##same thing with site class 0 p, background w and foreground w
	#proportion
		line<-(grep("proportion",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[[i]][j,1:4]=as.numeric(split[2:5])
	#w
		line<-(grep("foreground",x))	
		split=strsplit(x[line]," " )[[1]]
		split<-split[!split==""] 
		model.results[[i]][j,5:7]=as.numeric(split[3:5])
	##pull out tree
		line<-grep("tree length",x)
		tr4[[j]]<-read.tree(text=x[line+2])

	}  #end if


}  ##end j loop


save(model.results,paml.tree.stats,file="SeparateModelResults")

load("SeparateModelResults")


#which model is best?
table(unlist(apply(cbind(model.results[[1]][,9],model.results[[2]][,9],model.results[[3]][,9],model.results[[4]][,9]),1,which.min)))

AICs <- cbind(model.results[[1]][,9],model.results[[2]][,9],model.results[[3]][,9],model.results[[4]][,9])
liks <- cbind(model.results[[1]][,8],model.results[[2]][,8],model.results[[3]][,8],model.results[[4]][,8])
AICs<-na.omit(AICs)
liks<-na.omit(liks)

##remove AICs if model 2 or 4 are uninformative relative to 1 and 3. 
AICs[which(liks[,1]==liks[,2]),2]<-NA
AICs[which(liks[,3]==liks[,4]),4]<-NA

AICmin<-apply(AICs,1,function(x) min(x,na.rm=T))
AICmin<-cbind(AICmin,AICmin,AICmin,AICmin)


AICweights<- exp( -0.5 * (AICs-AICmin))
AICweights<-AICweights/rowSums(AICweights,na.rm=T)

par(mfrow=c(3,2))

hist(model.results[[4]][,5],main="A) w0",xlab="w0")
hist(model.results[[4]][model.results[[4]][,7]<=10,7],main="B) w2",xlab="w2")
hist(model.results[[4]][,1],main="C) p0",xlab="p0")
hist(model.results[[4]][,3],main="D) p2a",xlab="p2a")
hist((paml.tree.stats[paml.tree.stats[,2]<=100,2]),main="E) max pairwise dS",xlab="log max dS",xpd=F)
barplot(table(unlist(apply(cbind(model.results[[1]][,9],model.results[[2]][,9],model.results[[3]][,9],model.results[[4]][,9]),1,which.min))),main="F) # Orthologs supporting each model",xlab="Model")

#B and E truncated to show main part of distribution
 
 
##check covariates

##First p2a as y variable
y<-model.results[[3]][,3]

all.mono<-as.factor(rowSums(sep.results[,3:6])==4)
num.seqs<-rowSums(sep.results[,9:12])
gc.var<-apply(sep.results[,13:16],1,var)
lseq<-sep.results$lseq
dS<-paml.tree.stats[,2]
m<-lm(y~all.mono+num.seqs+lseq+log10(dS))
m1<-step(m)
anova(m1)
summary(m1)

##R2
anova(m1)[,2]/sum(anova(m1)[,2])

which(y<=quantile(dS,0.2,na.rm=T))

##Now try deltaAIC for model 1


concat.results<-read.csv("PierreConcat2034orthologs.csv")
concat.orthos<-sep.results[,1]%in%concat.results[,1]

##20th percentile for max dS pairwise from tree
quantile(paml.tree.stats[which(concat),2],na.rm=T,c(0.2))
#    20% 
#10.20463 

##list of concat genes
lower20.concat.ogs<-sep.results[which((paml.tree.stats[,2]<=quantile(paml.tree.stats[concat,2],na.rm=T,c(0.2)))&(concat)),1]
save(lower20.concat.ogs,"lower20.concat.ogs")


##make supp table

supp.table<-matrix(NA,nrow=nrow(model.results[[1]])*4,ncol=ncol(model.results[[1]]))

supp.table[seq(1,nrow(supp.table),4),]<-model.results[[1]]
supp.table[seq(2,nrow(supp.table),4),]<-model.results[[2]]
supp.table[seq(3,nrow(supp.table),4),]<-model.results[[3]]
supp.table[seq(4,nrow(supp.table),4),]<-model.results[[4]]

supp.table<-data.frame(supp.table)
colnames(supp.table)<-colnames(model.results[[1]])
supp.table<-cbind(OG=rep(rownames(model.results[[1]]),each=4),Model=rep(1:4,nrow(model.results[[1]])),supp.table)
