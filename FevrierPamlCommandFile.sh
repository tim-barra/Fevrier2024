#!/bin/bash

#SBATCH --job-name=timpaml
#SBATCH --time=12:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --mem=48g
#SBATCH --partition=short

module load Anaconda3/2022.05
source activate paml

##relevant directories

##a) Contains 4 codeml template files in ModelTemplateFiles.zip
templatedir="codeml_templates"  

homedir="paml_revisions"
cd ${homedir}

##all the tree files
g=(*.tre)

##maximum index
maxindex=$((${#g[@]}-1))

#RUN FOR A BATCH OF 20 BASED ON SLURM INDEX, index 1-243
n1=$(($SLURM_ARRAY_TASK_ID*20-20))
n2=$(($SLURM_ARRAY_TASK_ID*20-1))

if (( $n2 >= ${maxindex} ))
then 
	n2=${maxindex}
fi

##loop through OGs
for i in $(seq $n1 $n2)

do

##extract OG name
og="${g[$i]%_*}"

##make a temporary folder to store results
mkdir ${og}

##copy in the alignment and tree files
cp ${og}_* ${og}/

##copy in the template.ctl files
cp ${templatedir}/* ${og}/

##go to folder
cd ${og}

##replace in the correct file names to the template files
sed -i '' -e "s/orthogroup/${og}/g" *template*

##MODEL 1 - NEUTRAL MODEL
cp *1.template* codeml.ctl
codeml
cp *1.out ${homedir}/

##MODEL 2 - SITE MODEL
cp *2.template* codeml.ctl
codeml
cp *2.out ${homedir}/

##MODEL 3 - CONSTRAINED BRANCH-SITE MODEL
cp *3.template* codeml.ctl
codeml
cp *3.out ${homedir}/

##MODEL 4 - BRANCH-SITE MODEL
cp *4.template* codeml.ctl
codeml
cp *4.out ${homedir}/

##return to home directory and delete the temporary folder
cd ${homedir}
rm -r ${og}

done