##Code for Fevrier and Barraclough 2024 analysis of dS and concatenated analysis
##Tim Barraclough tim.barraclough@biology.ox.ac.uk
##Code for looking at gene trees and alignments of orthogroups
##Earlier steps produce output files used by later steps

##1) Fevrier2024.Separate.Jun2024.R compiles stats and sets up files for separate orthogroup analyses
##	 After running PAML analyses, it also compiles the output results

##2) Fevrier2024.Concatenated.Jun2024.R constructs a concatenated alignment
##	 After running PAML analyses, it also compiles the output results

##3) Fevrier2024.ConcatLevels.Jun2024.R constructs trees and reads in outputs
##	 from the analyses splitting at different levels of the phylogenetic hierarchy

##4) Fevrier.PAML.EvolverSimulations.R generates command files for simulations
##	 including manipulating branch lengths to explore the effects of saturation

##5) FevrierPamlCommandFile.sh is the bash script for submitting parallel PAML 
##   using the model templates in ModelTemplateFiles.zip

##6) paml_revisions.tar.gz contains the filtered/trimmed individual orthogroup alignments
##   and trees with branch labels

##7) concat.tar.gz contains the concatenated alignment of orthologs with all species represented
##   and tree with branch labels

##8) concat.lower20.tar.gz contains the concatenated alignment of the orthogroups with
##   the lowest max pairwise dS value (among those with all species represented)
