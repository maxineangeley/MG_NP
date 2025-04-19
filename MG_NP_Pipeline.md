Assuming the data is downloaded in my server database 
```bash
#log into the server, it will prompt for password
ssh mangeley@fossa.rcs.uvic.ca

#source allows us to use modules
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

#navigate to my directory with my fastq files
cd /project/biol470-grego/MAngeley/final

```

## Install Conda package manager
program page: https://www.anaconda.com/docs/getting-started/miniconda/main

#### Conda allows us to install most of the packages that we will need for the rest of the analysis 
```bash
#download the Miniconda (smaller version of conda) files
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

#install minconda 
bash Miniconda3-latest-*.sh

#restart bash - this just makes sure nothing previously installed will interact poorly with miniconda
source ~/.bashrc

#you need this for a lot of base bash stuff
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

#activate the base conda envrionment
source ~/miniconda3/bin/activate

#check for proper installation (current version = conda 25.3.1)
conda --version

```

## NanoPlot
current version: https://anaconda.org/bioconda/nanoplot
github: https://github.com/wdecoster/NanoPlot
citation: De Coster, W. & Rademakers, R. NanoPack2: population-scale evaluation of long-read sequencing data. _Bioinformatics_ **39**, btad311 (2023).
##### Create and activate virtual environment - new environments for each package avoid errors 
```bash
#create virtual envrionment and download the nanoplot package
conda create -n nanoplot_env nanoplot

#activate the virtual envrionment: allows you to use nanoplot
conda activate nanoplot_env

#check the NanoPlot version: should be 1.44.1
NanoPlot --version
```

#### Run NanoPlot!

#### First we will create a script that can run NanoPlot on both samples at once
```bash
#creat the script to run the NanoPlot
nano NanoPlot.sh
```

#### Paste the following script into that nano
```bash
#Notes for running 
	#input directory should be the location where your raw read fastq files are!
	#THINGS TO EDIT FOR YOUR USE
		#loop:change name of files if different than _minION_reads.fastq, it also works with gzipped files
		#base_name: change extension to the one that your file has
		#out_dir: change path to location where you want your output files to be stored
		#NanoPlot: change any parameters to suit your data and computer, given threads are too many for a normal computer


#path to fastq
INPUT_DIR="$1"

#loop that runs the NanoPlot
for fq_file in "$INPUT_DIR"/*_minION_reads.fastq; do

  #this extracts base sample name without the extension for use in naming files later
  base_name=$(basename "$fq_file" _minION_reads.fastq)
  #this defines where we want to put our file
  out_dir="$INPUT_DIR/nanoplot/NP2_${base_name}"

  #parameters for NanoPlot can be changed as needed, currently:
	  #--fastq tells it what file type
	  #--outdir: output location
	  #--threads: how many cores to use
	  #--plots + legacy: extra plots to make
	  #--verbose: lets you know whats going on
	  
  echo "Running NanoPlot for $base_name"
  NanoPlot \
    --fastq "$fq_file" \
    --outdir "$out_dir" \
    --threads 40 \
    --plots kde dot \
    --verbose
done
```

#### Now we run the script
```bash
#run the script with the directory that houses the fastq files as the input (this assumes youre in the directory with the files)
bash NanoPlot.sh .
```

#### In order to view the results we need to move it to the local computer (must be in a new terminal window, not logged into server)
```bash
scp -r "mangeley@fossa.rcs.uvic.ca:/project/biol470-grego/MAngeley/final/nanoplot/NP2*" "C:\Users\mnjan\OneDrive\Desktop\data\final\final\nanoplot"
```

#### You can open the report.html file on your computer to view the statistics and determine the next course of action 

## Chopper
current version: https://anaconda.org/bioconda/chopper
github: https://github.com/wdecoster/chopper
citation: De Coster, W. & Rademakers, R. NanoPack2: population-scale evaluation of long-read sequencing data. _Bioinformatics_ **39**, btad311 (2023).

#### After viewing the NanoPlot I decided to filter by a quality score of Q>8 and read length >1000
#### Installing chopper, a utility for manipulating fastq files 
```bash
#we will install chopper in the same envrionment because it uses similar packages and doesnt have negative interactions
conda install -c bioconda chopper

#test intall: should be 0.9.2
chopper --version 

#create script
nano chopper.sh
```

#### Paste into the chopper.sh script:
```bash
#Notes for running 
	#input directory should be the location where your read fastq files are!
	#THINGS TO EDIT FOR YOUR USE
		#loop: change name of files if different than _minION_reads.fastq
		#chopper: change quality (-q) and length (--minlength) settings as needed
		#output: change where filtered files should go (currently makes /filtered/ subdir)

INPUT_DIR="$1"

for fq_file in "$INPUT_DIR"/*_minION_reads.fastq; do

  #get sample base name without extension
  base_name=$(basename "$fq_file" _minION_reads.fastq)

  echo "Filtering $base_name"

  #run chopper with quality and length cutoffs
  chopper \
    --input "$fq_file" \
    --quality 8 \
    --minlength 1000 \
    --threads 40 \
    > "$INPUT_DIR/${base_name}_filtered.fastq"

done
```

#### Run chopper.sh with input as directory containing the fastq files (in this case current)
```bash
bash chopper.sh .
```

## Filtered NanoPlot

#### We will redo NanoPlot to visualize the differences that filtering made
```bash
nano filteredNanoPlot.sh
```

#### Paste into the nano
```bash
#Notes for running 
	#input directory should be the location where your read fastq files are!
	#THINGS TO EDIT FOR YOUR USE
		#loop:change name of files if different than _minION_reads.fastq, it also works with gzipped files
		#base_name: change extension to the one that your file has
		#out_dir: change path to location where you want your output files to be stored
		#NanoPlot: change any parameters to suit your data and computer, threads are too many for a normal computer

#path to fastq
INPUT_DIR="$1"

#loop that runs the NanoPlot and scp
for fq_file in "$INPUT_DIR"/*_filtered.fastq; do

  #base sample name without the extension
  base_name=$(basename "$fq_file" _filtered.fastq)
  out_dir="$INPUT_DIR/nanoplot/filtered/NP2_filtered_${base_name}"


  #parameters for NanoPlot can be changed as needed, currently:
	  #--fastq tells it what file type
	  #--outdir: output location
	  #--threads: how many cores to use
	  #--plots + legacy: extra plots to make
	  #--verbose: lets you know whats going on

  #parameters for this can be changed as needed, currently
  echo "Running NanoPlot for $base_name"
  NanoPlot \
    --fastq "$fq_file" \
    --outdir "$out_dir" \
    --minlength 1000 \
    --threads 40 \
    --plots kde dot \
    --verbose

done
```

#### Run the script as before
```bash
bash filteredNanoPlot.sh .
```

#### Transfer to computer (in a new terminal window)
```bash
scp -r "mangeley@fossa.rcs.uvic.ca:/project/biol470-grego/MAngeley/final/nanoplot/filtered/*" "C:\Users\mnjan\OneDrive\Desktop\data\final\nanoplot\"
```

#### Compare key statistics like average read length and quality (we will visualize this later)

### metaFlye assembly
current version: https://anaconda.org/bioconda/flye
github: https://github.com/mikolmogorov/Flye
citation: Kolmogorov, M. _et al._ metaFlye: scalable long-read metagenome assembly using repeat graphs. _Nat. Methods_ **17**, 1103–1110 (2020).

#### Now we make a new environment to avoid conflicts with NanoPlot packages (Flye contains metaFlye)
```bash
#deactivate old conda environment
conda deactivate 

#create a new one and install flye at the same time
conda create -n flye_env flye

#activate the environment
conda activate flye_env

#test install: should be 2.9.5-b1801
flye --version 

```

#### Create script to run the metaflye assembly
```bash
nano metaflye.sh
```

#### Paste into nano
```bash
#Notes for running 
	#input directory should be the location where your read fastq files are!
	#THINGS TO EDIT FOR YOUR USE
		#loop: change name of files if different than _minION_reads.fastq, it also works with gzipped files
		#base_name: change extension to the one that your file has
		#out_dir: change path to location where you want your output files to be stored
		#flye: change any parameters to suit your data and computer, threads are too many for a normal computer

INPUT_DIR="$1"

for fq_file in "$INPUT_DIR"/*_filtered.fastq; do

  #extract sample name without extension using basename
  base_name=$(basename "$fq_file" _filtered.fastq)
  #define output location
  out_dir="$INPUT_DIR/metaflye/MF_${base_name}"

  echo "Running metaFlye for $base_name"


  #run Flye with your chosen parameters, currently:
	  #--nano-raw: tells it that we are using nanopore data
	  #--out-dir: tells it where to put the output files
	  #--meta: tells it to use metafly, loosens requirements for assembly
	  #--threads: how many cores
	    
  flye \
    --nano-raw "$fq_file" \
    --out-dir "$out_dir" \
    --meta \
    --threads 40

done
```

#### Run script (again inputting the directory with the fastq files)
```bash
bash metaFlye.sh .
```

## metaQUAST
current version: https://anaconda.org/bioconda/quast
github:https://github.com/ablab/quast
manual: https://quast.sourceforge.net/

#### Now that we have assembled we can check the quality with metaQUAST (as subset of QUAST)

#### Another new environment
```bash
#deactivate old environment
conda deactivate

#create environment and install quast 
conda create -n quast_env -c bioconda quast

#check if it is working: should be 5.3.0
quast --version

#create script to run quast on our samples
nano quast.sh
```

#### Paste into quast.sh
```bash
#Notes for running 
	#input directory should be the location where your metaflye assembly file folders are!
	#THINGS TO EDIT FOR YOUR USE
		#loop: change path to assembly if necessary
		#base_name: sed removes the prefix that i chose, change if necessary
		#out_dir: change path to location where you want your output files to be stored
		#flye: change any parameters to suit your data and computer, threads are too many for a normal computer


INPUT_DIR="$1"

for sample in "$INPUT_DIR"/metaflye/MF_*/assembly.fasta; do
  # Extract sample name (e.g., 910 or 917)
  base_name=$(basename "$(dirname "$sample")" | sed 's/^MF_//')

  echo "Running MetaQUAST for MF_${base_name}..."

#current parameters: 
	#--threads: numb of cores
	#--min-conting: (this is min for whole program)
	#--gene-finding: this attempts to map your contigs to reference genomes, makes it take longer but if you have no idea whats in your sample it can be nice
	#-o is your output location

  metaquast.py "$sample" \
    --threads 30 \
    --min-contig 1000 \
    --gene-finding \
    -o "metaquast/MQ_${base_name}"
done
```

#### Run the script
```bash
bash quast.sh .
```

#### Transfer files to open in your computer browser
```bash
scp -r "mangeley@fossa.rcs.uvic.ca:/project/biol470-grego/MAngeley/final/metaquast/*" "C:\Users\mnjan\OneDrive\Desktop\data\final\metaquast\"
```


## metaBAT2
current version: https://anaconda.org/bioconda/metabat2
manual: https://bitbucket.org/berkeleylab/metabat/src/master/README.md

#### Now we will bin the assembly contigs to try and separate them by taxonomy if possible

```bash
#deactivate old env
conda deactivate

#create new env and install metabat2
conda create -n metaBAT_env -c bioconda metabat2

#activate environment
conda activate metaBAT_env

#check the install, should be version 2.17
metabat2 -h
```

####  In order to run metaBAT we need to prep the files
- We need a mapped, sorted, and indexed bam file and a depth file
- Because we did metaQUAST we can estimate the lengths of our contigs
- We will filter by length > 2500 because it improves the quality of bins by reducing overlaps
#### We will make a script to do this: 
```bash
nano MB_prep.sh
```

```bash
#load the required modules
#this lests us map OG reads back to the contigs
module load minimap2
#this lets us manipulate the bam files
module load samtools
#this lets us filter the fasta
module load seqtk

INPUT_DIR="$1"

#loop: change location of assembly if needed
for asm in "$INPUT_DIR"/metaflye/MF_*/assembly.fasta; do

  #defines file location and name for reference later
  sample_dir=$(dirname "$asm")
  sample=$(basename "$sample_dir" | sed 's/^MF_//')

  #filter contigs >2500 bp to improve quality
  seqtk seq -L 2500 "$asm" > "$sample_dir/assembly.filtered.fasta"

  #index filtered assembly
  minimap2 -d "$sample_dir/assembly.mmi" "$sample_dir/assembly.filtered.fasta"

  #map reads to my assembled contigs
  minimap2 -t 40 -ax map-ont "$sample_dir/assembly.mmi" $sample_filtered.fastq | \
    samtools view -Sb - > ${sample}_tmp.bam

  #sort and index bam
  samtools sort -@ 40 ${sample}_tmp.bam -o ${sample}_sorted.bam
  samtools index ${sample}_sorted.bam

  #generate depth file
  jgi_summarize_bam_contig_depths \
    --outputDepth ${sample}_depth.txt \
    ${sample}_sorted.bam

  #clean up stuff we dont need for the binning process
  rm ${sample}_tmp.bam ${sample}_sorted.bam ${sample}_sorted.bam.bai

done

```


#### Now we will make a script that will run the metaBAT on those files
```bash
nano metabat.sh

#make directory for file ouputs
mkdir -p metabat_bins
```

#### Paste into the nano
```bash

INPUT_DIR="$1"

#loop that runs the metabat2

for asm in "$INPUT_DIR"/metaflye/MF_*/assembly.filtered.fasta; do
  #defines the location and name of the sample
  sample_dir=$(dirname "$asm")
  sample=$(basename "$sample_dir" | sed 's/^MF_//')

  echo "running metabat2 on $sample"

#runs the metabat2
	#-i: input, -a:depth file, -o: output location, -m: minimum length of conigs, --unbinned: also creates a bin with the leftovers, --verbose: noisy, --t: threads
  metabat2 \
    -i "$asm" \
    -a "${sample}_depth.txt" \
    -o "metabat_bins/MB_${sample}" \
    -m 1500 \
    --unbinned \
    --verbose \
    -t 40

done
```

## Preparation for visualization:

#### This creates a very useful file that combined all the bins into one fasta. They are named by bin and unique contig number. 
* We will use this for visualization as well as creating a BLAST database

```bash
#this loop grabs all the bins and renames each contig, it takes the basename for the bin (base) and the contig name and pastes them as a fasta heading for each file
#it then pastes them all into one big fasta file

for bin in metabat_bins/*.fa; do
  base=$(basename "$bin" .fa)
  awk -v prefix="$base" '/^>/{print ">" prefix "_" substr($0,2); next}1' "$bin"
done > all_bins_renamed.fasta

```


## Optional: BLASTing against a custom database

#### Prep for blasting by making database:
```bash
#this loads the required module for blasting
module load blast+

#this makes a database out of our named contigs.
	#-in: input fasta
	#-dbtype: type of db, nucleotide, protein etc. 
	#-out: output location
	#-parse_seqids: keeps our nicely named contig names
	
makeblastdb -in all_bins_renamed.fasta -dbtype nucl -out all_bins_db -parse_seqids

#create databses for future use
mkdir ref_genomes
mkdir blast
```

#### create script file
```bash
nano url_to_hits.sh
```
#### Paste into the script
```bash
#how to use: bash url_to_blast.sh <input_dir> <name> '<url>'
module load blast+ 

input_dir=$1
name=$2
url=$3

#download and unzip the genome
curl -L -o $name.zip $url

#this unzips the downloaded file
unzip -q $name.zip -d $input_dir/ref_genomes/$name

#get rid of .zip we dont need it anymore
rm $name.zip

#move and rename genome fasta for easier access in the first layer of the folder
mv $input_dir/ref_genomes/$name/ncbi_dataset/data/GC*/GC*.fna $input_dir/ref_genomes/$name/$name.fna

#run blast
	#-query: our reference genome
	#-db: path to the database created earlier and then the database prefix - this is the prefix for all the files in your database. 
	#-out: output location
	#-outfmt: this is important! this is the format that your output is in... 6 is a tsv file that we can use in the following steps!
blastn \
  -query $input_dir/ref_genomes/$name/$name.fna \
  -db bins_db/all_bins_db \
  -out $input_dir/blast/blast_${name}.tsv
  -outfmt 6

#this does a couple things in one
	#it sorts that blast output bu column 12 which is the bit score, this tells you how good your contig hit is
	#then it cuts out only the second column which is our contig identifier, it includes bin and contig number
	#then it gets rid of all the duplicates that mapped multiple times for some reason
	#then it only takes the top 10 hits - change this number if you want more
	#then it prints that as a list of contig names
sort -k12 -nr $input_dir/blast/blast_${name}.tsv | cut -f2 | uniq | head -n 10 > top_hits_${name}.txt

#this loads the required module for subsectioning our fasta
module load seqtk

#this takes the text file previously created and grabs only those contigs in the list and prints them to a new fasta file.
#this fasta file can be used to blast on the NCBI website if you want a better idea of the taxonomy of your hits. 
seqtk subseq all_bins_renamed.fasta top_hits_${name}.txt > $input_dir/blast/${name}_contigs.fasta
```
#### That script file runs the whole process, However if you already have a reference genome downloaded you can modify it to begin at the blastn part


# Visualization:

## GC coverage by bin

#### We will need the all_bins_renamed.fasta file created earlier for this 
```bash
scp -r "mangeley@fossa.rcs.uvic.ca:/project/biol470-grego/MAngeley/final/all_bins_renamed.fasta" "C:\Users\mnjan\OneDrive\Desktop\data\final\bins"
```

## R script for the process
Adapted from this video (where I found out about the package) https://www.youtube.com/watch?v=M9vRfYC3ULc



```r

#this will allow us to manipulate the fasta file easily
install.packages('seqinr')
library('seqinr')

#tidyverse for general r packages
install.packages('tidyverse')
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

#this tells r where our file is without having to set new working directory
bin_file <- "C:/Users/mnjan/OneDrive/Desktop/data/final/bins/all_bins_renamed.fasta"

#this creates a list type dataset with our sequences in the SeqFastdna format.
all_bin_fasta <- read.fasta(bin_file,as.string = TRUE, seqtype = "DNA", seqonly = FALSE)

#creats an empty dataframe to add to later
GCout_df <- data.frame()

#this is adapted from referenced youtube video
#k = individual contig in the list 
for(k in 1:length(all_bin_fasta)){
  seq_to_GC <- getSequence(all_bin_fasta[k]) #this extracts the sequence from the list
  seq_to_GC <- unlist(seq_to_GC)             #this makes it not a list anymore, much easier to manipulate later
  GCout <- GC(seq_to_GC)                     #this calculates the GC content of each contig and prints it to a value
  GCout_df[k,1] <- getName(all_bin_fasta[k]) #this places the name of the contig in the previously creaetd dataframe
  GCout_df[k,2] <- GCout                     #this places the GC content in that same datafrome
}

#this renames the titles of the columns
GCout_df <- GCout_df %>% rename("Name" = "V1", "GC Content" = "V2")

#this separates out important information from the Name of the contig, Sample, bin and contig number
GCout_df <- GCout_df %>%
  separate(Name, c("MF","Sample.bin","contig", "contig_number"), "_", remove=F) %>%
  separate(Sample.bin, c("Sample", "1or5", "Bin"),remove=T)

#this gets rid of columns that we will not use
GCout_df$`1or5` <- NULL
GCout_df$contig <- NULL

#this makes sure that our bins and samples are correctly refered to as characters (some bins are numbers some are letters)
GCout_df <- GCout_df %>%
  mutate(Bin = as.character(Bin)) %>%
  mutate(Sample = as.character(Sample)) 


#this is a nice color palatte
install.packages("wesanderson")
library("wesanderson")


#now we plot! we are making a violin plot to visualize the distribution within bins and difference between bins

#first defines data frame and axes, as well as coloring by sample
GC_plot <- ggplot(GCout_df, aes(x = Bin, y = `GC Content`, fill = Sample)) +
  #type of plot and style
  geom_violin(scale = "width", trim = FALSE, color = "black", size = 0.8, alpha = 0.8) +
  #this facets the graph by sample and labels each one with the Host species name
  facet_wrap(~ Sample, scales = "free_x", labeller = as_labeller(c("910" = "Nearcticorpus canadense", 
                                                                   "917" = "Spelobia bumamma"))) + 
  #this applies labels to the axes
  labs(x = "Bin",
       y = "GC Content (%)") +
  #these apply themes and adjust text size and formatting
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none", strip.text = element_text(size = 16, face = "italic")) +
  #these make sure that the axes are scaled well and are properly labeled
  scale_y_continuous(labels = function(y) format(y * 100)) +
  scale_x_discrete(labels = c("lowDepth" = "Low Depth", "unbinned" = "Unbinned")) +
  #this puts a dot for the mean GC content
  stat_summary(fun = mean, geom = "point", fill = "black", shape = 21, size = 2, position = position_dodge(width = 0.9)) +
  #this colors it!
  scale_fill_manual(values = wes_palette("Chevalier1")) 

#this saves it to my computer
ggsave("C:/Users/mnjan/OneDrive/Desktop/data/final/bins/GC_Plot.jpeg", GC_plot, width = 8, height = 5)

```



## Now we will visualize some key differences caused by filtering the initial reads for quality and length

#### The NanoPlot output includes a summary text file that we can parse for information. And graph more visually that the nanoplot summary. 

```r
#install packages required:
install.packages('seqinr')
library('seqinr')

install.packages('tidyverse')
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library('tidyverse')


#this is a nice color palatte
install.packages("wesanderson")
library("wesanderson")


#locate files 
NP_dirs <- "C:/Users/mnjan/OneDrive/Desktop/data/final/nanoplot"

#this picks out only the files that we want and allows us to reference them in a loop later
summary_files <- Sys.glob(file.path(NP_dirs, "NP*/NanoStats.txt"))

#this large loop allows us to create individual dataframes with each samples summaries being separate for now

for (file in summary_files) {
  #define the name of the sample (file name) 
  sample_id <- basename(dirname(file))
  #and creates a value with all the information from the file in it
  lines <- readLines(file)
  
  #Step 1: Summary section - this contains the info about read lengths and N50 
  #lines 2-9 have the useful information
  summary_lines <- lines[2:9]
  #this makes a tibble dataframe that has 3 columns - Section (summmary), Metric (what is being measures e.g N50 or mean read length) and Value (which is the nuber reported in the summary)
  summary_df <- tibble(
    Section = "Summary",
    Metric = str_trim(str_extract(summary_lines, ".*(?=:)")), #extracts text before the :
    Value = as.numeric(gsub(",", "", sapply(strsplit(summary_lines, ":"), function(x) x[2])))  #extracts stuff after the : (as numeric) - this is so complicated because i was having trouble with commas and the number 50 in N50 
  )
  
  #Step 2: Quality information - After I did this i realized that i dont use this info but still good to have. 

  #instead of picking specified lines this selects ones with the >Q in them which specifies the quality information
  q_lines <- lines[grepl("^>Q", lines)] 

  #this makes a dataframe by treating the information as a table, this is because it is formated similarly with bp total, % and Mbp on one line following the title,
  #we split it into different columns for easier tracking
  q_df <- read.table(text = gsub("[>%():Mb]", "", q_lines),
                     col.names = c("Q", "Count", "Percent", "Mb"),
                     stringsAsFactors = FALSE) %>%
    pivot_longer(cols = -Q, names_to = "Metric", values_to = "Value") %>% #this makes unique rows for each one 
    mutate(
      Section = "Quality Cutoffs", 
      Metric = paste0(Q, "_", Metric), #NOTE: necessary to have a unique metric (the Q part) - messes up downstream if not
      Value = as.numeric(gsub(",", "", Value)) #important comma removal here
    ) %>%
    select(Section, Metric, Value) 
  
  #this combines into one final tidy data frame for each sample
  final_df <- bind_rows(summary_df, q_df) %>%
    mutate(Sample = sample_id)
    
  #this actually names and creates each dataframe for the specific sample
  assign(paste0(sample_id, "_df"), final_df)
  
  #get rid of temporary files
  rm(final_df, q_df, summary_df, sample_id, summary_lines, q_lines, lines, file)
}


#Now that we have our individual samples we can combine them for easier plotting - a long format combine is best for plotting

#this creates a value that is the names of all the df from before (to reference later)
df_names <- ls(pattern = "NP2")

#this makes a dataframe with sections for those names and all the information from the specfied dataframe
all_df <- bind_rows(lapply(df_names, get))


#this makes the names much prettier (i dont end up using these) and also adds a group column (this is important for the graphs)
all_df <- all_df %>%
  mutate(
    Sample = gsub("^NP2_", "", Sample),
    Sample = gsub("^filtered", "Filtered", Sample),
    Group = ifelse(grepl("910", Sample), "910", 
                   ifelse(grepl("917", Sample), "917", NA))
  )


#Now we can create graphs we will do 3 for now these are very similar so i will only annotate the first one (except for the differences)

#this creates a new df that only keeps the info we want for the graph - in this case total bases - and removes a pesky NA row proble,
plot_TB_df <- all_df %>%
  select(Section, Metric, Value, Sample, Group) %>% #selects all the columns
  filter(Section == "Summary", Metric == "Total bases") %>% #filters only the one that we are going to graph 
  filter(!is.na(Sample), !is.na(Value)) #removes NA


#Plotting! - we will do a bar chart this time 
TB_plot <- ggplot(plot_TB_df, aes(x = Sample, y = Value / 1e6, fill = Group)) + #because the numbers are so high divide by 1 million to get Mb
  geom_col(color = "black", size = 0.8, alpha = 0.8) + 
  facet_wrap(~ Group, scales = "free_x",                        #this creates facets labeled by host and sorted by group (that we defined earlier!)
             labeller = as_labeller(c(
               "910" = "Nearcticorpus canadense",
               "917" = "Spelobia bumamma"
             ))) +
  labs(                                                         #this labels axes (and gets rid of x axis label)
    y = "Total Bases (Mb)",
    x = " "
  ) +
  theme_minimal() +                                             #themeing for text, apperance, and labels for filtered and unfiltered
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 16, face = "italic"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = c("910.1" = "Unfiltered", "Filtered_910.1" = "Filtered", "917.5" = "Unfiltered", "Filtered_917.5" = "Filtered")) + 
  scale_fill_manual(values = wesanderson::wes_palette("Chevalier1")) 



#Now we repeat the process for Mean read length 

plot_ARL_df <- all_df %>%
  select(Section, Metric, Value, Sample, Group) %>%
  filter(Section == "Summary", Metric == "Mean read length") %>%
  filter(!is.na(Sample), !is.na(Value))


#CHANGES: 
	#-removed the /1e6 and group labels (we will stack later)
	#-changesd the y-axis title 
	
ARL_plot <- ggplot(plot_ARL_df, aes(x = Sample, y = Value, fill = Group)) +
  geom_col(color = "black", size = 0.8, alpha = 0.8) +
  facet_wrap(~ Group, scales = "free_x",
             labeller = as_labeller(c(
               "910" = "",
               "917" = ""
             ))) +
  labs(
    y = "Average Read Length (bp)",
    x = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 16, face = "italic"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = c("910.1" = "Unfiltered", "Filtered_910.1" = "Filtered", "917.5" = "Unfiltered", "Filtered_917.5" = "Filtered")) +
  scale_fill_manual(values = wesanderson::wes_palette("Chevalier1"))




#Repeat one more time for N50
plot_N50_df <- all_df %>%
  select(Section, Metric, Value, Sample, Group) %>%
  filter(Section == "Summary", Metric == "Read length N50") %>%
  filter(!is.na(Sample), !is.na(Value))


#CHANGES: 
	#-removed the /1e6 and group labels (we will stack later)
	#-changed the y-axis title 
	
N50_plot <- ggplot(plot_N50_df, aes(x = Sample, y = Value, fill = Group)) +
  geom_col(color = "black", size = 0.8, alpha = 0.8) +
  facet_wrap(~ Group, scales = "free_x",
             labeller = as_labeller(c(
               "910" = "",
               "917" = ""
             ))) +
  labs(
    y = "N50 (bp)",
    x = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 16, face = "italic"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = c("910.1" = "Unfiltered", "Filtered_910.1" = "Filtered", "917.5" = "Unfiltered", "Filtered_917.5" = "Filtered")) +
  scale_fill_manual(values = wesanderson::wes_palette("Chevalier1"))


#this is a package for putting graphs together 
install.packages('patchwork')
library(patchwork)

#this stacks on top of each other (| for side to side)
combined_plot <- TB_plot / ARL_plot / N50_plot

#saves to computer
ggsave("C:/Users/mnjan/OneDrive/Desktop/data/final/combined_Plot.jpeg", combined_plot, width = 6, height = 9)

```