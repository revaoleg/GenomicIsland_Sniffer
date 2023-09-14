import sys, os.path, string, array, math, time, shelve
import sys, os
path = os.getcwd()
sys.path.append(os.path.join(path,"lib"))
import cli
import lib, blast

###############################################################################

if __name__ == "__main__":
    options = {"-c":"MGE",     # name of the scenario
               "-u":"yes",              # use BLASTN to filter rrn clusters
               "-n":"no",               # use BLASTN to search tRNA on the borders of GI
               "-l":"8000",             # sliding window length
               "-b":"2000",             # sliding window big step
               "-m":"500",              # sliding window medium step
               "-s":"100",              # sliding window small step
               "-e":"Contrasting/Iteration",    # refinement No | Contrasting | Iteration | Contrasting/Iteration 
               "-i":"input",            # input folder
               "-o":"output",           # output folder
               "-f":"fasta+gbk",            # save GI sequemces: no | fasta | gbk | gbk+fasta
               "-v":"yes"               # save SVG file
            }
    arguments = sys.argv[1:]
    if arguments:
        for i in range(0,len(arguments)-1,2):
            key = arguments[i]
            if key not in options:
                raise IOError("Unknown argument " + key + "!")
            if i <= len(arguments)-2:
                options[key] = arguments[i+1]
        oInterface = cli.Interface(options)
    else:
        response = ""
        while response != "Q":
            print("SeqWord Sniffer 3.0 2023/09/03")
            print()
            print("Copy whole genome sequence in fasta format to folder 'input' and press Y+Enter")
            print("To show menu of run parameters, press M+Enter")
            print("To exit, press Q+Enter")
            print()
            response = input(": ").upper()
            if response=="Q":
	            exit()
            elif response=="M":
	            oInterface = cli.Interface([],True)
            elif response=="Y":
	            oInterface = cli.Interface()
            else:
	            continue
    
