
##Code for Fevrier and Barraclough 2024 analysis of dS and concatenated analysis
##Tim Barraclough tim.barraclough@biology.ox.ac.uk
##Code for running different levels of splitting into ESUS


library(ape)
library(phangorn)

input.dir<-"concat_levels"
output.dir<-input.dir
setwd(input.dir)


gene.tr<-read.tree(paste0(input.dir,"/ConcatePAMLmodel4.tre"))



	species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
					"Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
					"Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
					"Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
	order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
	order.colours<-as.integer(as.factor(order.names))



##function to output tree and control template file

output.files<-function(tr.tmp,index) {
		tree.file.name<-gsub(".tre",paste0(".",index,".tre"),tree.files[i])
		write.tree(tr.tmp,paste0(output.dir,"/",tree.file.name))
		x<-scan(paste0(prefix,".template.label.ctl"),what=character(),sep="\n")
		x[1]<-paste0("       seqfile = ",prefix,"_dna.trim.phy")
		x[2]<-paste0("     treefile = ",tree.file.name)
		x[3]<-paste0("      outfile = ",tree.file.name,".out")
		x[11]<-"        model = 2"
		x[15]<-"        icode = 0"
		write(x,file=paste0(output.dir,"/",gsub("\\.tre","\\.template.label.ctl",tree.file.name)))
}


	
	##extract the names removing the part after "_" and match them with names above
	tip.names<-gene.tr$tip.label	
	tip.orders<-order.names[pmatch(tip.names ,species.names,duplicates.ok=T)]	
	tip.colours<-order.colours[pmatch(tip.names ,species.names,duplicates.ok=T)]
##plot the tree using tip.colours by Order - comment out when running loop
	plot(gene.tr,tip.color=tip.colours)

##label this tree as per the order model

	within.order<-NULL
	mrca.nodes<-NULL
		for (k in (1:4)) {
			mrca.node<-getMRCA(gene.tr,tip=which(tip.colours==k))
			desc<-Descendants(gene.tr,mrca.node,type="all")
			within.order<-c(within.order,desc[desc>length(gene.tr$tip.label)])
			mrca.nodes<-c(mrca.nodes,mrca.node)}
    gene.tr$node.label<-rep("#1",gene.tr$Nnod)
    gene.tr$node.label[1]<-""	
    gene.tr$node.label[(within.order-length(gene.tr$tip.label))]<-"#0"
    gene.tr$tip.label<-paste0(gene.tr$tip.label, "#0")


 ##make directories to receive the different labelled trees
	for (i in (0:16)) {
		dir.create(paste0(output.dir,"/tree",i))
		file.copy(from="PierreConcat2034orthologs.phy",to=paste0(output.dir,"/tree",i))
		file.copy(from="Model4.template.label.ctl",to=paste0(output.dir,"/tree",i,"/codeml.ctl"))
		}
					
	##a) Coalesce order nodes to compare more aggregated model
		
		##ancestors of the order nodes
		anc.nodes<-unique(gene.tr$edge[,1][gene.tr$edge[,2]%in% mrca.nodes])
		##first coalesce the highest number

par(mfrow=c(1,2))
		tr.tmp<-gene.tr
		##tree 0 = one coalescence
			i=0
			tr.tmp$node.label[tr.tmp$edge[tr.tmp$edge[,1]==max(anc.nodes),2]-length(tr.tmp$tip.label)]<-"#0"
			plot(tr.tmp,show.node.label=1)
			write.tree(tr.tmp,paste0(output.dir,"/tree",i,"/concat.tree",i,".tre"))
			x<-scan(paste0(output.dir,"/tree",i,"/codeml.ctl"),what=character(),sep="\n")
			x[2]<-gsub("orthogroup_dna.trim.paml.tre",paste0("concat.tree",i,".tre"),x[2])
			x[3]<-gsub("orthogroup_dna.trim.model4.out",paste0("concat.tree",i,".out"),x[3])
			write(x,file=paste0(output.dir,"/tree",i,"/codeml.ctl"))

		##tree 1 = two coalescences
			i=1
			tr.tmp$node.label[tr.tmp$edge[tr.tmp$edge[,1]==sort(anc.nodes,decreasing=T)[2],2]-length(tr.tmp$tip.label)]<-"#0"
			plot(tr.tmp,show.node.label=1)
			write.tree(tr.tmp,paste0(output.dir,"/tree",i,"/concat.tree",i,".tre"))
			x<-scan(paste0(output.dir,"/tree",i,"/codeml.ctl"),what=character(),sep="\n")
			x[2]<-gsub("orthogroup_dna.trim.paml.tre",paste0("concat.tree",i,".tre"),x[2])
			x[3]<-gsub("orthogroup_dna.trim.model4.out",paste0("concat.tree",i,".out"),x[3])
			write(x,file=paste0(output.dir,"/tree",i,"/codeml.ctl"))
		
	##b) Split mrca nodes to compare more divided model
	
	par(mfrow=c(1,4))
	
	##Split one order at a time: goes C, D, H, L -ptera, tree 2 to 5
	
		tmp.index=1
		for (j in (1:length(mrca.nodes))) {
		tr.tmp<-gene.tr
		desc<-tr.tmp$edge[,2][tr.tmp$edge[,1]==mrca.nodes[j]]	
		if (any(desc<=length(tr.tmp$tip.label))) {
			tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]]<-gsub("#0","#1",tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]])}
		if (any(desc>length(tr.tmp$tip.label))) {
			tr.tmp$node.label[desc[which(desc>length(tr.tmp$tip.label))]-length(tr.tmp$tip.label)]<-"#1"}
		
		i=j+tmp.index
		plot(tr.tmp,show.node.label=1)
		write.tree(tr.tmp,paste0(output.dir,"/tree",i,"/concat.tree",i,".tre"))
			x<-scan(paste0(output.dir,"/tree",i,"/codeml.ctl"),what=character(),sep="\n")
			x[2]<-gsub("orthogroup_dna.trim.paml.tre",paste0("concat.tree",i,".tre"),x[2])
			x[3]<-gsub("orthogroup_dna.trim.model4.out",paste0("concat.tree",i,".out"),x[3])
			write(x,file=paste0(output.dir,"/tree",i,"/codeml.ctl"))


		}

	par(mfrow=c(2,3))

	##Split two orders at a time: goes C, D, H, L -ptera, 1:2,1:3,1:4,2:3 etc. tree 6 to 11
		tmp.index=6
		for (j in (1:3)) {
			for (k in ((j+1):4)) {
				
		tr.tmp<-gene.tr
		desc<-tr.tmp$edge[,2][tr.tmp$edge[,1]==mrca.nodes[j]]	
		if (any(desc<=length(tr.tmp$tip.label))) {
			tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]]<-gsub("#0","#1",tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]])	}
		if (any(desc>length(tr.tmp$tip.label))) {
			tr.tmp$node.label[desc[which(desc>length(tr.tmp$tip.label))]-length(tr.tmp$tip.label)]<-"#1"}
		desc<-tr.tmp$edge[,2][tr.tmp$edge[,1]==mrca.nodes[k]]	
		if (any(desc<=length(tr.tmp$tip.label))) {
			tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]]<-gsub("#0","#1",tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]])}
		if (any(desc>length(tr.tmp$tip.label))) {
			tr.tmp$node.label[desc[which(desc>length(tr.tmp$tip.label))]-length(tr.tmp$tip.label)]<-"#1"}

		i=tmp.index
		plot(tr.tmp,show.node.label=1,cex=1.5)
		write.tree(tr.tmp,paste0(output.dir,"/tree",i,"/concat.tree",i,".tre"))
			x<-scan(paste0(output.dir,"/tree",i,"/codeml.ctl"),what=character(),sep="\n")
			x[2]<-gsub("orthogroup_dna.trim.paml.tre",paste0("concat.tree",i,".tre"),x[2])
			x[3]<-gsub("orthogroup_dna.trim.model4.out",paste0("concat.tree",i,".out"),x[3])
			write(x,file=paste0(output.dir,"/tree",i,"/codeml.ctl"))
		tmp.index=tmp.index+1
	
				}}
			
	##Split three orders at a time
		tmp.index<-12
		for (k in (1:4)) {
		which.nodes<-combn(1:4,3)[,k]
		tr.tmp<-gene.tr
		for (j in (which.nodes)) {
		desc<-tr.tmp$edge[,2][tr.tmp$edge[,1]==mrca.nodes[j]]	
		if (any(desc<=length(tr.tmp$tip.label))) {
			tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]]<-gsub("#0","#1",tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]])	}
		if (any(desc>length(tr.tmp$tip.label))) {
			tr.tmp$node.label[desc[which(desc>length(tr.tmp$tip.label))]-length(tr.tmp$tip.label)]<-"#1"}
		}
	i=tmp.index
			plot(tr.tmp,show.node.label=1,cex=1.5)
	write.tree(tr.tmp,paste0(output.dir,"/tree",i,"/concat.tree",i,".tre"))
			x<-scan(paste0(output.dir,"/tree",i,"/codeml.ctl"),what=character(),sep="\n")
			x[2]<-gsub("orthogroup_dna.trim.paml.tre",paste0("concat.tree",i,".tre"),x[2])
			x[3]<-gsub("orthogroup_dna.trim.model4.out",paste0("concat.tree",i,".out"),x[3])
			write(x,file=paste0(output.dir,"/tree",i,"/codeml.ctl"))
		tmp.index=tmp.index+1
		}
		
	##And lastly, split all four
		
	j<-1
	desc<-tr.tmp$edge[,2][tr.tmp$edge[,1]==mrca.nodes[j]]	
		if (any(desc<=length(tr.tmp$tip.label))) {
			tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]]<-gsub("#0","#1",tr.tmp$tip.label[desc[which(desc<=length(tr.tmp$tip.label))]])	}
		if (any(desc>length(tr.tmp$tip.label))) {
			tr.tmp$node.label[desc[which(desc>length(tr.tmp$tip.label))]-length(tr.tmp$tip.label)]<-"#1"}
	i=tmp.index
				plot(tr.tmp,show.node.label=1,cex=1.5)
	write.tree(tr.tmp,paste0(output.dir,"/tree",i,"/concat.tree",i,".tre"))
			x<-scan(paste0(output.dir,"/tree",i,"/codeml.ctl"),what=character(),sep="\n")
			x[2]<-gsub("orthogroup_dna.trim.paml.tre",paste0("concat.tree",i,".tre"),x[2])
			x[3]<-gsub("orthogroup_dna.trim.model4.out",paste0("concat.tree",i,".out"),x[3])
			write(x,file=paste0(output.dir,"/tree",i,"/codeml.ctl"))
		tmp.index=tmp.index+1


##copy in sequence files for the selected trees



###LOOK AT RESULTS OUTPUT FROM PAML

##load packages
	library(ape)
	library(phangorn)
	library(seqinr)

##set up color-blind palette with 4 colours
	palette(c( "#E69F00","#D55E00","#56B4E9","#009E73"))
	
##paths to tree and sequence alignment inputs 
	input.path<-"concat_levels_out"

##extract species names and allocate order names
	species.names<-c("Synanthedon","Agriopis","Pieris","Zygaena","Tinea",
					"Coremacera","Sarcophaga","Eupeodes","Bombylius","Bibio",
					"Pyrochroa","Malachius","Coccinella","Rhagonycha","Ocypus","Nebria",
					"Bombus","Vespa","Anoplius","Ichneumon","Tenthredo")
	order.names<-c(rep("Lepidoptera",5),rep("Diptera",5),rep("Coleoptera",6),rep("Hymenoptera",5))
	order.colours<-as.integer(as.factor(order.names))

out.path<-"concat_levels_out"
out.files<-dir(out.path)

##set up results tables for the models
	model.results<-matrix(NA,nrow=length(out.files),ncol=11)
	colnames(model.results)<-c("p0","p1","p2a","p2b","w0","w1","w2","LnL","AIC","MaxPairwise","MaxPairwiseDs")
	rownames(model.results)<-out.files

##set up tree list
	trlist<-list()
		
check.finished<-array(0,length(out.files))

		##Model 4
for (j in (1:length(out.files))) {
			
	x<-scan(file=paste0(out.path,"/",out.files[j]), what=character(),sep="\n")
	
	if (any(grepl("lnL",x))) {
		
		check.finished[j]<-1
	
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
	
	}}
	
	
##load in the order-level results
	d<-read.csv("ConcatenatedModelResults.csv")		

	model.results<-rbind(model.results,concat.tree1.5.out=d[4,-1])
	
	to.order<-as.numeric(substring(rownames(model.results),regexpr("tree", rownames(model.results))+4, regexpr("out", rownames(model.results))-2))
	
	model.results<-model.results[order(to.order),]
	
	##level 0 is actually -1, level 1 is -2
	
	model.results[1:2,]<-model.results[2:1,]
	
	model.names<-c("-2","-1","orders","C","D","H","L","CD","CH","CL","DH","DL","HL","CDH","CDL","CHL","DHL","CDHL")

	if (is.na(model.results[1,8])) {model.results[1,8]= -22447478.657042; model.results[1,9]=(-2*model.results[1,8]+2*4)}

	barplot((model.results[,9]-min(model.results[,9],na.rm=T))/10000,names.arg=model.names,las=2,ylab="delta AIC",xlab="clade structure",col="lightblue")

##AIC weights
	AICw<-exp( -0.5 * (model.results[,9]-min(model.results[,9],na.rm=T)))/sum(exp( -0.5 * (model.results[,9]-min(model.results[,9],na.rm=T))))

##plot AA seq divergence

seqs<-read.dna("PierreConcat2034orthologs.phy")
aa.seqs<-trans(seqs)
aadist<-dist.aa(aa.seqs)/ncol(aa.seqs)
library(MASS)
plot(mm$points[pmatch(species.names ,rownames(mm$points)),],col=as.factor(order.names),pch=substring(order.names,1,1),ylab="MDS2",xlab="MDS1")



##PLOT ALTERNATIVE MODELS

setwd("concat_levels")

my.files<-dir()

par(mfrow=c(5,4))
	model.names<-c("-1","-2","C","D","H","L","CD","CH","CL","DH","DL","HL","CDH","CDL","CHL","DHL","CDHL")

for (i in (c(1,0,2:16))) {
	tr<-read.tree(paste0("tree",i,"/concat.tree",i,".tre"))

##between nodes
 bet.nodes<-which(tr$node.label=="#1")+ length(tr$tip.label)
 
 bcol<-array("black",nrow(tr$edge))
 bcol[which(tr$edge[,2]%in%bet.nodes)]<-"grey"
plot(tr,show.tip.label=1,edge.color=bcol,main=model.names[i+1])
}


