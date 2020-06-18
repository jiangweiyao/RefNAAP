# Reference Based Nanopore Amplicon Analysis Pipeline

This pipeline takes basecalled, demultiplexed, and trimmed fastq files from Nanopore sequencing:
1. It QC the fastq files using fastqc and aggregates the reports using multiqc.
2. It maps the reads to a set of reference sequences (max size 4GB) using Minimap2.
3. It finds the reference sequence segments where the fastq files mapped and extract those sequence segments. You specify the top number of reference sequences to include as well as the minimum average coverage required.
4. It uses Medaka to generate the consensus of the Nanopore reads using the reference sequence as the scaffold. 
5. It aligns the the Medaka consensus sequence against the reference sequence to find indels, and remove these from the consensus. Indel errors in homopolyers are common occurences in Nanopore sequencing. These are usually removed by hand by comparing to a reference sequence. This pipeline automates this step. 

This workflow depends on Bioconda to install the environmental dependencies.
And then run the attached script using the environment.

## Summary - Installation 
1. Clone Repository 
2. Install Conda if not already in environment
3. Create conda environment

## Summary - How to run after installation.
1. Activate conda environment - `conda activate RefNAAP`
2. Run the GUI - `~/RefNAAP/RefNAAP_GUI.py`
4. Test data in test folder

## Clone this code using GIT

### Install git for Debian systems using the following command (if necessary)
```
sudo apt update
sudo apt install git
```

##Installation directions 
These instructions install the code into your home path. Change the instructions if appropriate. 

### Clone the code from repository
```
cd ~
git clone https://github.com/jiangweiyao/RefNAAP.git
```

### Install Both Miniconda and Environment (Skip next 2 steps)
You can run the prepackaged script install_all.sh to install both Miniconda and make the Conda environment into your home directory (recommended) by using the following command
```
. ~/RefNAAP/install_all.sh
```


### Install Miniconda (if no Conda is install on system). 
You can run the prepackaged script install_miniconda.sh to install into your home directory (recommended) by using the following command
```
. ~/RefNAAP/install_miniconda.sh
```

Detailed instruction on the the Miniconda website if anything goes wrong:
https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

### Clone the environment. Need to do once.

We use conda to create an environment (that we can activate and deactivate) to install our dependent software and resolve their dependencies. This environment is called "RefNAAP". 

```
conda create -n RefNAAP medaka python gooey tabulate pandas fastqc multiqc bioconductor-decipher -y
```

Conda will automatically figure out and install all dependencies for you. It might take some time for Conda to resolve the dependencies and install everything for you. 

## Run the code.

Activating your environment makes the software you installed in that environment available for use. You will see "(RefNAAP)" in front bash after activation.
```
conda activate RefNAAP
```
Run the code (change path if it is installed else where.)
```
~/RefNAAP/RefNAAP_GUI.py
```

When you are finished running the workflow, exit out of your environment by running `conda deactivate`. Deactivating your environment exits out of your current environment and protects it from being modified by other programs. You can build as many environments as you want and enter and exit out of them. Each environment is separate from each other to prevent version or dependency clashes. The author recommands using Conda/Bioconda to manage your dependencies.

## Author
- Jiangwei Yao

## License 
Apache License 2.0
