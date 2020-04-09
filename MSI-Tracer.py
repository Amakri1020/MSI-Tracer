from Bio import SeqIO
import twobitreader
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import os
from os import listdir
from os.path import isfile, join
from tqdm import tqdm
import mmap
from tqdm import tqdm_notebook
import gc

def reverse(seq : str) -> str:
    """A simple function for generating the opposite strand in a DNA sequence.
    
    Takes a string as input, e.g. 'ACTG', and outputs the reverse strand, 'CAGT'."""
    rev = ""
    for base in reversed(seq):
        if base.upper() == 'T':
            rev += 'A'
        elif base.upper() == 'A':
            rev += 'T'
        elif base.upper() == 'G':
            rev += 'C'
        elif base.upper() == 'C':
            rev += 'G'
        else:
            #print("Invalid Sequence")
            break
    return rev

 def extract_polyA_pure(seq : str, min_length=10) -> str:
    """Extracts a pure polyA of at least min_length bases following the 17-base Alu primer.
    
    Parameters
    ----------
    seq : str
        The DNA sequence beginning with the 17-base Alu tail primer
    min_length : int, optional
        Minimum length required to consider the pure polyA valid, default is 10

    Returns
    -------
    polyA
        The pure polyA, unless one is not found, in which case an empty string is returned
    """
    
    countA = 0
    polyA = ""
    
    for i in range(17, len(seq)):
        if seq[i].upper() == 'A':
            countA += 1
            polyA += seq[i]
        elif seq[i].upper() != 'A' and countA >= min_length:
            return polyA
        else:
            polyA = ""
            return polyA
        
    return polyA

def parse_file_pure_primer(file : str) -> list:
    relevant = []
    with open(file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[0][0] != 'M':
                print('Whoops')
                pass
            else:
                if(row[9][:17] == 'GAGCGAGACTCCGTCTC'):
                    polyA = extract_polyA_pure(row[9], 10)
                    if not polyA == "":
                        relevant.append([row[2], row[3], polyA, len(polyA), row[9]])
                if(row[9][-17:] == 'GAGACGGAGTCTCGCTC'):
                    polyA = extract_polyA_pure(reverse(row[9]), 10)
                    if not polyA == "":
                        relevant.append([row[2], row[3], polyA, len(polyA), row[9]])
    return relevant

# Helper function so TQDM can show progress, possibly not worth for large files, slows things down
def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def info_single_file(fname : str, outdir: str, required_occurrences : int):
    """Generates an info file for the given .CSV file, containing all reads occurring more than 'required_occurrences' times.
    
    Parameters
    ----------
    fname : str
        Filename of the .csv file for which the info file is being generated. This file must be in
        the format of a converted SAM file, e.g. by exporting a SAM file to CSV in Excel,
        or by using the sam_to_csv function provided in this file.
    outdir : str
        Directory for info files to be written into
    required_occurrences : int
        Minimum number of reads from the same genomic position required to qualify a position
        for inclusion in the output info file.

    Returns
    -------
    None
        A .csv file will be written to the current directory.
    """
    
    group = []
    last = 0
    head = 'TGGTCTCGATCTCCTGACCTC'
    head_r = 'GAGGTCAGGAGATCGAGACCA'
    tail = 'GAGCGAGACTCCGTCTCA'
    tail_r = 'TGAGACGGAGTCTCGCTC'
    
    # Get first position number for group identification
    with open(fname) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[0][0] != 'M':
                pass
            else:
                last = row[3]
                break
        
    # Main loop
    with open(fname) as csvfile, open('temp.csv', 'w', newline='') as myfile:
        wr = csv.writer(myfile)
        reader = csv.reader(csvfile, delimiter=',')
        
        # Initialize all info variables
        total_reads = 0
        primer_polyAs = 0
        total_groups = 0
        #unique_groups_primer = 0
        thresh_groups = 0
        primers = 0
        th = 0
        hh = 0
        tt = 0
        h = 0
        t = 0
        
        
        for row in tqdm(reader, total=get_num_lines(fname)):
            total_reads += 1
            if row[0][0] != 'M':
                pass
            else:
                # Individual primer occurrence counts in full sequence
                seq = row[9]
                if head in seq or head_r in seq:
                    h += 1
                if tail in seq or tail_r in seq:
                    t += 1
                if tail in seq and head_r in seq:
                    th += 1
                elif tail in seq and tail_r in seq:
                    tt += 1
                elif head in seq and head_r in seq:
                    hh += 1
                
                current = row[3]
                
                # On position change, write current position group to output if we have enough occurrences
                if current != last:
                    total_groups += 1
                    if len(group) >= required_occurrences:
                        thresh_groups += 1
                        for entry in group:
                            wr.writerow(entry)
                    group = []
                
                # Get pure polyA if we find primer at beginning or end of sequence 
                if(row[9][:17] == 'GAGCGAGACTCCGTCTC'):
                    primers += 1
                    polyA = extract_polyA_pure(row[9], 10)
                    if not polyA == "":
                        primer_polyAs += 1
                        group.append([row[2], row[3], polyA, len(polyA), row[9]])
                        
                # Reverse sequence if primer is at the end so all seq's are in terms of top strand
                elif(row[9][-17:] == 'GAGACGGAGTCTCGCTC'):
                    primers += 1
                    polyA = extract_polyA_pure(reverse(row[9]), 10)
                    if not polyA == "":
                        primer_polyAs += 1
                        group.append([row[2], row[3], polyA, len(polyA), row[9]])
                last = current
                    
                        
    split = fname.split('/')
    outname = outdir + split[-1][:-4] +'_info'+str(required_occurrences)+'.csv'
        
    with open(outname, 'w', newline='') as myfile, open('temp.csv') as csvfile:
        wr = csv.writer(myfile)
        wr.writerow(['Tail-Head: '+str(th)+ ' -- Tail-Tail: ' + str(tt)+' -- Head-Head: '+str(hh)
                     + ' -- All Head Occurrences: '+str(h)+ ' -- All Tail Occurrences: '+str(t)])
        wr.writerow(['Total Reads: ' + str(total_reads) + ' -- Total Groups: ' + str(total_groups) 
                     + ' -- Normalized(100k): ' + str(int((100000 * (float(total_groups) / float(total_reads)))))
                     + ' -- Groups with over ' + str(required_occurrences) + ' occurrences: ' + str(thresh_groups)])
        wr.writerow(['All primers: ' + str(primers) 
                     + '-- Pure polyAs next to tail primer: ' + str(primer_polyAs)]) #+ ' -- Total tail-polyA groups: ' + str(unique_groups_primer)])
        wr.writerow([])
        wr.writerow(['Chromosome', 'Position', 'PolyA', 'Length', 'Original Sequence'])
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            wr.writerow(row)
    os.remove('temp.csv')

def shared_reads(fname, fname_tumor, outname1, outname2, verbose=False):
    """Compares two info files and outputs a .csv file for each, containing only shared reads between the files.
    
    Parameters
    ----------
    fname : str
        Filename of the first info file, typically the one that does not contain any tumor sample although this is
        not required.
    fname_tumor : str
        Filename of the second info file, typically the one containing some trace of the tumor.
    outname1 : str
        Name of the first output file containing the reads from fname that are shared with fname_tumor
        Typical output name would be something like 'MMX_MMY_shared_reads.csv'
    outname2 : str
        Name of the second output file containing the reads from fname_tumor that are shared with fname
        Typical output name would be something like 'MMY_MMX_shared_reads.csv'
    verbose : boolean, optional
        If true, some progress messages will be written to stdout while the function executes.

    Returns
    -------
    None
        Two .csv files will be written to the output paths specified in the function arguments.
    """
    
    IDs = set()
    IDs_tumor = set()
    
    # Iterate through both info files and gather all IDs (format 'Chromosome#-Position') in a set
    with open(fname) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i in range(5):
            next(reader)
        for row in reader:
            IDs.add(str(row[0])+'-'+str(row[1]))
            
    if verbose: print(fname + ' IDs added.') 

    with open(fname_tumor) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i in range(5):
            next(reader)
        for row in reader:
            IDs_tumor.add(str(row[0])+'-'+str(row[1]))
            
    if verbose: print(fname_tumor + ' IDs added.') 

    # Set intersection gives the positions shared by both files
    shared = IDs & IDs_tumor
    count = len(shared)
    
    # Write shared reads files, for each only keeping rows where the position is in the shared set
    with open(fname) as csvfile, open(outname1, 'w', newline='') as myfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i in range(5):
            next(reader)
        wr = csv.writer(myfile)
        wr.writerow(['Shared Groups: ' + str(count)])
        wr.writerow([])
        wr.writerow(['Chromosome', 'Position', 'PolyA', 'Length', 'Original Sequence'])
        for row in reader:
            key = str(row[0])+'-'+str(row[1])
            if key in shared:
                wr.writerow(row)
                
    if verbose: print(outname1 + ' written.')
                
    with open(fname_tumor) as csvfile, open(outname2, 'w', newline='') as myfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i in range(5):
            next(reader)
        wr = csv.writer(myfile)
        wr.writerow(['Shared Groups: ' + str(count)])
        wr.writerow([])
        wr.writerow(['Chromosome', 'Position', 'PolyA', 'Length', 'Original Sequence'])
        for row in reader:
            key = str(row[0])+'-'+str(row[1])
            if key in shared:
                wr.writerow(row)
                
    if verbose: print(outname2 + ' written.')

def generate_graphs(file1 : str, file2 : str, outdir : str, r3 : bool, r6 : bool, X : float, verbose=False):
    """Generates comparison graphs for two input shared_reads files if the given conditions are met.
    
    Parameters
    ----------
    file1 : str
        Name of the first shared_reads file.
    fname_tumor : str
        Name of the second shared_reads file.
    outdir : str
        Directory for all generated graphs to be written to (automatically created if it does not already exist).
    r3 : bool
        Rule 3: The two smallest polyA lengths must be from the tumor sample, and at least one must have
        more than two occurrences. If set to True, this rule must be met for any graphs generated.
    r6 : bool
        Rule 6: X % of tumor reads occur before the first non-tumor read. If set to True, this rule must
        be met for any graphs generated.
    X : float
        Percentage of reads required to occur before the first non-tumor read for rule 6. Range 0 to 1.
    verbose : boolean, optional
        If true, some progress messages will be written to stdout while the function executes.

    Returns
    -------
    None
        Graphs will be written in outdir for all samples meeting the input conditions.
    """
    
    lengths_normal = {}
    lengths_tumor = {}
    
    if X < 0 or X > 1:
        print("X must be a value between 0 and 1")
        return None
    
    ## Initialize keys with empty dicts for every position
    ## Key format 'Chromosome#-Position'
    with open(file1, 'r') as myfile:
        cr = csv.reader(myfile)
        # Skip header rows
        for i in range(4):
            next(cr)
            
        for row in cr:
            key = str(row[0])+'-'+str(row[1])
            if not key in lengths_normal:
                lengths_normal[key] = [] 
                lengths_tumor[key] = []
    
    if verbose : print("Dict keys populated.")
        
    ## For each file, store the lengths of all polyAs associated with every key
    with open(file1,  'r') as myfile:
        cr = csv.reader(myfile)
        for i in range(4):
            next(cr)
            
        for row in cr:
            key = str(row[0]+'-'+str(row[1]))
            lengths_normal[key].append(len(row[2]))
    if verbose : print("Normal Lengths stored.")
    
    with open(file2,  'r') as myfile:
        cr = csv.reader(myfile)
        for i in range(4):
            next(cr)
            
        for row in cr:
            key = str(row[0]+'-'+str(row[1]))
            lengths_tumor[key].append(len(row[2]))
    if verbose : print("Tumor Lengths Stored")
    
    ## For each key, determine whether it qualifies to be graphed, and if so generate and save graphs in outdir.
    for key, value in lengths_normal.items():
        f1 = lengths_normal[key]
        f2 = lengths_tumor[key]
        # Range of lengths for this key
        index = np.arange(min(min(f1),min(f2)), max(max(f1),max(f2)+1))
        # Lists containing identically ordered counts of polyA lengths for each sample
        f1count = []
        f2count = []
        for idx in index:
            f1count.append(f1.count(idx))
            f2count.append(f2.count(idx))

        
        two_lowest = False
        cond1 = False
        cond2 = False
        
        ## Rule 3: The two smallest polyA lengths must be from the tumor sample, and at least one must have
        ## more than two occurrences.
        for idx in index:
            # If both conditions are met, Rule 3 is met and we set it to True and break
            if(cond1 and cond2):
                two_lowest = True
                break
                
            # If we find a non-zero read from the non-tumor sample before meeting our conditions, Rule 3 is broken
            if f1.count(idx) != 0:
                break
            elif f2.count(idx) >= 2: # One of our reads must have 2 or more occurrences (Cond2).
                # If Cond2 is already met and we have a second read with 2+ occurrences, that also satisfies Cond1.
                if cond2:
                    cond1 = True
                else:
                    cond2 = True
            elif f2.count(idx) >= 1:  # One of our reads needs only 1 or more occurrences (Cond1).
                cond1 = True
                
        ## Rule 6: X % of tumor reads occur before the first non-tumor read.
        bluecount = 0
        for idx in index:
            if f1.count(idx) == 0:
                bluecount += f2.count(idx)
            else:
                break
        if (float(float(bluecount) / float(sum(f2count)))) > X:
            rule6 = True
        else:
            rule6 = False
        
        # Labels for the graphs to be created, ex. 'DSN/Shared_reads/MM12_MM5_shared_reads.csv' -> 'MM12'
        # TODO redo this using path module and account for other possible file names.
        split = file1.split('/')
        label1 = split[-1][:4]
        split = file2.split('/')
        label2 = split[-1][:4]
        
        # Boolean formula, applies whichever rules were set to True in function input
        condition = (not r3 or two_lowest) and (not r6 or rule6)
        
        if condition:
            fig, ax = plt.subplots()
            bar_width = 0.35
            opacity = 0.9
            ax.bar(index, f1count, bar_width, alpha=opacity, color='r', label=label1)
            ax.bar(index+bar_width, f2count, bar_width, alpha=opacity, color='b', label=label2)
            plt.xticks(np.arange(min(index), max(index)+1, 1.0))
            ax.set_xlabel('Length of PolyA')
            ax.set_ylabel('Reads')
            split = key.split('-')
            # Separate chromosome# and position for Graph title
            title = (split[0], split[1])
            ax.set_title('Chromosome ' + title[0] + ': ' + title[1])
            ax.legend()
            direc = outdir
            if not os.path.exists(direc):
                os.makedirs(direc)
            #TODO : Fix this in the *extremely* unlikely case of duplicate sequence positions in separate genomes
            plt.savefig(os.path.join(direc, entry[1]))
            plt.close(fig)
            gc.collect()

def shared_graphs(dirs : list, outname : str):
    """Identifies shared graphs in a list of directories containing graphs, outputs results to a csv file.
    
    Parameters
    ----------
    dirs : list
        List of directories
    fname_tumor : str
        Name of the second shared_reads file.
    outdir : str
        Directory for all generated graphs to be written to (automatically created if it does not already exist).
    r3 : bool
        Rule 3: The two smallest polyA lengths must be from the tumor sample, and at least one must have
        more than two occurrences. If set to True, this rule must be met for any graphs generated.
    r6 : bool
        Rule 6: X % of tumor reads occur before the first non-tumor read. If set to True, this rule must
        be met for any graphs generated.
    X : float
        Percentage of reads required to occur before the first non-tumor read for rule 6. Range 0 to 1.
    verbose : boolean, optional
        If true, some progress messages will be written to stdout while the function executes.

    Returns
    -------
    None
        Graphs will be written in outdir for all samples meeting the input conditions.
    """
    allgraphs = []
    for direc in dirs:
        onlyfiles = [f for f in listdir(direc) if isfile(join(direc, f))]
        allgraphs.append(onlyfiles)
    
    res = list(set.intersection(*map(set, allgraphs)))
    
    with open(outname, 'w', newline='') as myfile:
        wr = csv.writer(myfile)
        for x in res:
            wr.writerow([x])
    
    return res

