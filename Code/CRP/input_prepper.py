#I/O
import sys
import getopt
import os
import glob

#data
import pandas as pd

OUTDIR = os.getcwd()


def accum_files(inpdir, search_str):
    inp_folders = glob.glob(os.path.join(inpdir, "*pat*"))
    path_search = os.path.join(inpdir, search_str)
    files = glob.glob(path_search, recursive=True)
    return inp_folders, files

def make_input_files(folders, files):
    for folder in folders:
        subj = folder.split("/")[-1]
        inpdir = os.path.join(OUTDIR,'input_files' ,subj)
        if not os.path.exists(inpdir):
            os.mkdir(inpdir)
        open(os.path.join(OUTDIR, 'input_files',subj, 'inp.csv'), 'w').close() #clears files
    
    files = filter_pulses(files)

    print(f"Writing out input: {len(files)}")
    for i,f in enumerate(files):
        subj, stim, ma = f.split("_")
        with open(os.path.join(OUTDIR, 'input_files', subj,'inp.csv'), 'a') as f:
            line =f'{subj},{stim},{ma}\n'
            f.write(line)


    #TODO address this newline issue
def filter_pulses(files):
    files = [f.split('/')[-1].split("_")[0:3] for f in files ]
    files = set(["_".join(f) for f in files])
    return files


def main(argv):
    inpdir = '/mnt/ernie_main/000_Data/SPES/data/preprocessed/'

    opts, _ = getopt.getopt(argv, "i:", ["inpdir="])
    for opt, arg in opts:
        if opt in ('-i', '--inpdir'):
            inpdir = arg
    
    search_str = os.path.join(inpdir,'*pat*/*pat*_CCEP_single_pulses/*.mat')

    folders, files = accum_files(inpdir, search_str)
    make_input_files(folders, files)

if __name__ == "__main__":
    main(sys.argv[1:])

