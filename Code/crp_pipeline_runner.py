import sys
from pathos.pools import ParallelPool
import getopt
import crp
    
def call_crp(inp_line):
    subj,stim,ma = inp_line.split(",")
    print(f"Running on {subj}, with stim sesh: {stim} at {ma}")
    pathout ='/mnt/ernie_main/Ghassan/ephys/data/'
    crp.run_crp_pipeline(subj, pathout, ma, stim)

def main(argv):
    input_f = ''
    cores = ''
    opts, _ = getopt.getopt(argv,"i:c:",["input_f=","cores="])
    for opt, arg in opts:
        if opt in ("-i", '--input_f'):
            input_f = arg
        elif opt in ("-c", "--cores"):
            cores = int(arg)
    print(f'cores specified: {cores}')
    #parse args into lists to execute
    pool = ParallelPool(nodes=cores)
    with open(input_f, 'r') as f:
        input_lines = f.readlines()
    pool.map(call_crp, input_lines)

if __name__ == "__main__":

    main(sys.argv[1:])