#!/usr/bin/env python

import sys
import os
import glob
import re
import pandas as pd
from tabulate import tabulate
from datetime import date
from gooey import Gooey, GooeyParser

@Gooey(program_name='RefMAAP', 
        default_size=(720, 900),
        progress_regex=r"^progress: (?P<current>\d+)/(?P<total>\d+)$",
        progress_expr="current / total * 100")

def main():
    
    local_path = os.path.dirname(os.path.realpath(__file__))
    #print(local_path)
    data_path = f"{local_path}"
    scaffold_helper = f"{local_path}/scaffold_cutter.R"
    gapfixer_helper = f"{local_path}/gapfixer.R"
    now = date.today()

    cli = GooeyParser(description="Reference Based Minion Amplicon Analysis Pipeline")
    required_args = cli.add_argument_group("Input Output Location", gooey_options={'columns': 1, 'show_border': True})
    required_args.add_argument('--InputFolder', help="Folder containing barcoded fastq", required=True, widget='DirChooser')
    required_args.add_argument('--OutputFolder', help="Output Folder", required=False, default=f'~/refnaap_results/output_{now}', widget='DirChooser')
    required_args.add_argument('--RefFile', help="Reference File ", required=False, default=f'{local_path}/Americas2.fasta', widget='FileChooser')

    parser = cli.add_argument_group("Optional Arguments", gooey_options={'columns': 2, 'show_border': True})
    parser.add_argument('--TopN', help="The top N reference sequences with the most depth are analyzed.", type=int, required=False, default=1)
    parser.add_argument('--MinCov', help="Amplicon regions need a minimum of this average coverage number", type=int, required=False, default=10)
    parser.add_argument('--threads', help="Number of threads. More is faster if your computer supports it", type=int, required=False, default=4)
    parser.add_argument('--verbose', help = "Keep Intermediate Files", required=False, widget='BlockCheckbox', action='store_true', gooey_options={ 'checkbox_label': "Yes" })
    parser.add_argument('--model', help="Basecall Model", required=False, type=str, default='r941_min_high_g303')
    args = cli.parse_args()

    files = sorted([f for f in glob.glob(args.InputFolder+"/**", recursive = True) if re.search(r'(.*)\.((fastq|fq)(|\.gz))$', f)])   
    #InputFolder = os.path.expanduser(args.InputFolder)
    #files = sorted(glob.glob(InputFolder+"/*.fastq"))
    print(files)
    OutputFolder = os.path.expanduser(args.OutputFolder)

    qc_dir= f"{OutputFolder}/qc"
    assembly_dir= f"{OutputFolder}/assembly"
    os.system(f"mkdir -p {qc_dir}")
    os.system(f"mkdir -p {assembly_dir}")
    f=open(f"{OutputFolder}/qc.log", 'w+')
    g=open(f"{OutputFolder}/cov_stat.txt", 'w+')

    for i in range(0, len(files)):
        filec = files[i]

        base = os.path.splitext(os.path.basename(filec))[0]
        base = os.path.splitext(base)[0]
        #print(base)

        fastqc_cmd = f"fastqc {filec} -o {qc_dir}"
        f.write(fastqc_cmd+'\n')
        os.system(fastqc_cmd)

        minimap2_cmd = f"minimap2 -ax map-ont {args.RefFile} {filec} -t {args.threads} > {assembly_dir}/{base}.sam "
        f.write(minimap2_cmd+'\n')
        os.system(minimap2_cmd)
        samtools_cmd1 = f"samtools view -S -b {assembly_dir}/{base}.sam > {assembly_dir}/{base}.bam"
        f.write(samtools_cmd1+'\n')
        os.system(samtools_cmd1)
        samtools_cmd2 = f"samtools sort {assembly_dir}/{base}.bam -o {assembly_dir}/{base}.sorted.bam"
        f.write(samtools_cmd2+'\n')
        os.system(samtools_cmd2)
        samtools_cmd3 = f"samtools index {assembly_dir}/{base}.sorted.bam"
        f.write(samtools_cmd3+'\n')
        os.system(samtools_cmd3)
        samtools_cov = f"samtools depth {assembly_dir}/{base}.sorted.bam > {assembly_dir}/{base}.coverage"
        os.system(samtools_cov)
        f.write(samtools_cov+'\n')
        
        scaffold_cmd = f"{scaffold_helper} {args.TopN} {args.MinCov} {args.RefFile} {assembly_dir}/{base}.coverage {assembly_dir}/{base}_scaffold.fasta {assembly_dir}/{base}_cov.txt"
        os.system(scaffold_cmd)
        f.write(scaffold_cmd+'\n')

        medaka_cmd = f"medaka_consensus -i {filec} -d {assembly_dir}/{base}_scaffold.fasta -o {assembly_dir}/{base}_medaka -m {args.model} -t {args.threads}"
        f.write(medaka_cmd+'\n')
        os.system(medaka_cmd)

        gapfixer_cmd = f"{gapfixer_helper} {assembly_dir}/{base}_scaffold.fasta {assembly_dir}/{base}_medaka/consensus.fasta {assembly_dir}/{base}_final.fasta"
        f.write(gapfixer_cmd+'\n')
        os.system(gapfixer_cmd)
        
        cp_cmd = f"cp {assembly_dir}/{base}_final.fasta {OutputFolder}/{base}_final.fasta"
        f.write(cp_cmd+'\n')
        os.system(cp_cmd)

        if(os.path.isfile(f"{assembly_dir}/{base}_cov.txt")):
            df = pd.read_csv(f"{assembly_dir}/{base}_cov.txt")
            g.write(f"{base} Results: \n")
            g.write(tabulate(df, headers="keys", showindex=False, tablefmt="psql", floatfmt=".2f"))
            g.write(f"\n\n")
        else:
            g.write(f"{base} did not have enough aligned reads: \n\n")

        print("progress: {}/{}".format(i+1, len(files)))


    multiqc_cmd = f"multiqc {qc_dir} -o {OutputFolder}"
    f.write(multiqc_cmd+'\n')
    os.system(multiqc_cmd)
    f.close()
    #g.close()
    print(f"verbose status is {args.verbose}")
    if not args.verbose:
        os.system(f"rm -rf {qc_dir}")
        os.system(f"rm -rf {assembly_dir}")

if __name__ == "__main__":
    sys.exit(main())
