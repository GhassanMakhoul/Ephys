import sys
import os
from pathos.pools import ProcessPool
import getopt
import crp
    
def prep_dirs(inp_line):
    _, _, _, pathout = get_params(inp_line)
    if not os.path.exists(pathout):
        os.mkdir(pathout)
    fig_dir = os.path.join(pathout,"figs")
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)
    deriv_dir = os.path.join(pathout, 'derivatives')
    if not os.path.exists(deriv_dir):
        os.mkdir(deriv_dir)

def get_params(inp_line):
    params = inp_line.split(",")
    subj = params[0].strip("\n")
    stim = params[1].strip("\n")
    ma = params[2].strip("\n")
    print(f"Running on {subj}, with stim sesh: {stim} at {ma}")
    pathout ='/mnt/ernie_main/Ghassan/ephys/data/'
    pathout = os.path.join(pathout, subj,f'{stim}_{ma}')
    return subj, stim, ma, pathout

# def validate_input_file(inp_line, inp_dir):
#     """Returns true if input file exsts

#     Args:
#         inp_line (str): unparsed inputline
#     """
#     subj, stim, ma, pathout = get_params(inp_line)
#     f = f'{subj}_{stim}_{ma}_'
#     return os.path.isfile(inp_f)

def call_crp(inp_line):
    subj, stim, ma, pathout, = get_params(inp_line)
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
    pool = ProcessPool(nodes=cores)
    with open(input_f, 'r') as f:
        input_lines = f.readlines()
    for l in input_lines:
        prep_dirs(l)

    pool.map(call_crp, input_lines)

if __name__ == "__main__":

    main(sys.argv[1:])