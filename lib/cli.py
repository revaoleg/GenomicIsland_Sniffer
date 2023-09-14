import string, os, re
import lib, seq_io

# Command line interface
class Interface:
    def __init__(self,options=None,show_menu=False):
        self.oValidator = Validator()
        self.IO = seq_io.IO()
        self.scripts = {}
        self.open_scripts()
        self.task_list = {}
        self.options = {"-c":"MGE",# name of the scenario
                   "-u":"yes",              # use BLASTN to filter rrn clusters
                   "-n":"no",               # use BLASTN to search tRNA on the borders of GI
                   "-l":"8000",             # sliding window length
                   "-b":"2000",             # sliding window big step
                   "-m":"500",              # sliding window medium step
                   "-s":"100",              # sliding window small step
                   "-e":"Contrasting/Iteration",    # refinement No | Contrasting | Iteration | Contrasting/Iteration 
                   "-i":"input",            # input folder
                   "-o":"output",           # output folder
                   "-f":"fasta+gbk",            # save GI sequemces: no/fasta/gbk/gbk+fasta
                   "-v":"yes"               # save SVG file
                }
        if options:
            self.options.update(options)
            self.execute()
        elif show_menu:
            self.main_menu()
        else:
            self.execute()

    # Execute selected program
    def execute(self):
        options = self.oValidator.validate(self.options,self.scripts.keys())
        if options:
            self.options = options
        else:
            return
        self.task_list.update(self.scripts[self.options["-c"]])
        oMainWin = lib.Main(self.task_list,self.options)
        del oMainWin
    
    # show command prompt interface
    def edit_tasklist(self):
        e = E = "E"
        a = A = "A"
        r = R = "R"
        q = Q = "Q"
        response = ""
        while response != "Q":
            print()
            print("Edit task list")
            print()
            print("  E   Edit task in the list;")
            print("  A   Add a new task to the list;")
            print("  R   Remove a task from the list;")
            print("  Q   Quit;")
            try:
                response = input("?")
            except:
                continue
            response = string.upper(response)
            if response == "A":
                self.add_task()
                return
            elif response == "E":
                self.edit_task()
                return
            elif response == "R":
                self.remove_task()
                return
            else:
                continue

    def add_task(self,flg_subtrahend=None,flg_divisor=None):
        task_code = "D"
        wlength = 4
        norm = 1
        subtrahend = None
        divisor = None
        task_categories = {"GC":"GC-content",
                            "AT":"AT-content",
                            "GCS":"GC-skew",
                            "ATS":"AT-skew",
                            "D":"pattern deviation",
                            "GD":"generalized pattern deviation",
                            "PS":"pattern skew",
                            "GPS":"generalized pattern skew",
                            "RV":"relative variance",
                            "GRV":"generalized relative variance"}
        response = ""
        while response != "Q":
            print()
            print("Add a task")
            print()
            print("  C   task category:\t\t" + task_code + " (" + task_categories[task_code] + ")")
            print("  W   word length:\t\t" + str(wlength))
            print("  N   normalization:\t\t" + str(norm))
            if not flg_subtrahend and not flg_divisor:
                if subtrahend:
                    print("  S   subtrahend:\t\t" + subtrahend)
                else:
                    print("  S   subtrahend:\t\tNone")
                if divisor:
                    print("  D   divisor:\t\t\t" + divisor)
                else:
                    print("  D   divisor:\t\t\tNone")
            print("  A   add task to the list")
            print("  Q   quit;")
            try:
                response = input("?")
            except:
                continue
            response = string.upper(response)
            if response == "C":
                new_task = self.select_task(task_categories)
                if new_task:
                    task_code = new_task
            elif response == "W":
                val = input("Enter the length of the words in range from 2 to 7 ")
                try:
                    val = int(val)
                    if val > 1 and val < 8:
                        wlength = val
                except:
                    print()
                    print("Word length must be an integer in range from 2 to 7")
                    print()
                    continue
            elif response == "N":
                val = input("Enter the normalization value in range from 0 to word_length - 1: ")
                try:
                    val = int(val)
                    if val >= 0 and val < wlength:
                        norm = val
                except:
                    print()
                    print("Word length must be an integer in range from 0 to word_length - 1")
                    print()
                    continue
            elif response == "S":
                subtrahend = self.add_task(1)
                divisor = None
            elif response == "D":
                divisor = self.add_task(0,1)
                subtrahend = None
            elif response == "A":
                if task_code in ("GC","AT","GCS","ATS"):
                    new_task = task_code
                else:
                    new_task = "n" + str(norm) + "_" + str(wlength) + "mer:" + task_code
                if flg_subtrahend or flg_divisor:
                    return new_task
                if subtrahend:
                    new_task += "-"+subtrahend
                if divisor:
                    new_task += "/"+divisor
                if new_task not in self.task_list:
                    self.task_list[new_task] = {}
                    self.task_list[new_task].update(self.set_conditions(new_task,{"mode":"absolute","condition":"bigger than","val1":0,"val2":0}))
                self.options["-c"] = "User defined"
                print()
                return
            else:
                continue

    def remove_task(self):
        index = None
        while index != 0:
            tasks = self.task_list.keys()
            if not tasks:
                return
            if len(tasks) > 1:
                tasks.sort()
            print()
            print("Remove task from the list")
            print()
            print("  0   escape;")
            for i in range(len(tasks)):
                if self.task_list[tasks[i]]["condition"] in ("bigger than","smaller than"):
                    str_val = str(self.task_list[tasks[i]]["val1"])
                else:
                    str_val = str(self.task_list[tasks[i]]["val1"]) + " and " + str(self.task_list[tasks[i]]["val2"])
                print(("  " + str(i+1) + "   " + tasks[i] +
                               " \t" + self.task_list[tasks[i]]["condition"] +
                               " (" + self.task_list[tasks[i]]["mode"] + ") " + str_val))
            print()
            index = input("Select the task by its index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(tasks):
                    del self.task_list[tasks[index-1]]
                    self.options["-c"] = "User defined"
                elif index != 0:
                    print()
                    print("\tWrong task index")
                    print()
                else:
                    return
                continue
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(tasks)))
                print()
                continue
            
    def edit_task(self):
        index = None
        while index != 0:
            tasks = self.task_list.keys()
            if not tasks:
                return
            if len(tasks) > 1:
                tasks.sort()
            print()
            print("Edit conditions of the selected task")
            print()
            print("  0   escape;")
            for i in range(len(tasks)):
                if self.task_list[tasks[i]]["condition"] in ("bigger than","smaller than"):
                    str_val = str(self.task_list[tasks[i]]["val1"])
                else:
                    str_val = str(self.task_list[tasks[i]]["val1"]) + " and " + str(self.task_list[tasks[i]]["val2"])
                print(("  " + str(i+1) + "   " + tasks[i] +
                               " \t" + self.task_list[tasks[i]]["condition"] +
                               " (" + self.task_list[tasks[i]]["mode"] + ") " + str_val))
            print()
            index = input("Select the task by its index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(tasks):
                    self.task_list[tasks[index-1]].update(self.set_conditions(tasks[index-1],self.task_list[tasks[index-1]]))
                    if self.task_list[tasks[index-1]]["condition"] == "between":
                        val = self.task_list[tasks[index-1]]["val1"]
                        if val > self.task_list[tasks[index-1]]["val2"]:
                            self.task_list[tasks[index-1]]["val1"] = self.task_list[tasks[index-1]]["val2"]
                            self.task_list[tasks[index-1]]["val2"] = val
                    self.options["-c"] = "User defined"
                elif index != 0:
                    print()
                    print("\tWrong task index")
                    print()
                else:
                    return
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(tasks)))
                print()
                continue
            
    def open_scripts(self):
        self.scripts = {}
        path = os.path.join("lib","scripts.txt")
        data,fname = self.IO.open(path)
        if data:
            data = data.replace("\r","")
            data = data.split("END\n")
            for script in data:
                if not script:
                    continue
                script_data = script.split("\n")
                script_name = script_data[0]
                self.scripts[script_name] = {}
                for i in range(1,len(script_data)-1,2):
                    task = script_data[i]
                    self.scripts[script_name][task] = {}
                    settings = script_data[i+1].split(",")
                    for values in settings:
                        name,value = values.split(":")
                        try:
                            value = float(value)
                        except:
                            pass
                        self.scripts[script_name][task][name] = value
            if self.scripts:
                return
        
        self.scripts = {"MGE":
            {
            "n0_4mer:D":{"mode":"sigmas","condition":"bigger than","val1":2.0,"val2":0.0},
            "n1_4mer:GRV/n1_4mer:RV":{"mode":"absolute","condition":"bigger than","val1":2.0,"val2":0.0},
            "n0_4mer:PS":{"mode":"absolute","condition":"smaller than","val1":55.0,"val2":0.0},
            },
            "Fitness genes":
            {
            "n1_4mer:RV":{"val2":0.0,"val1":0.0,"mode":"sigmas","condition":"bigger than"},
            "n0_4mer:D":{"val2":0.0,"val1":0.3,"mode":"fraction","condition":"smaller than"},
            "n0_4mer:PS":{"val2":0.0,"val1":25.0,"mode":"absolute","condition":"smaller than"},
            },
            "Ribosomal RNA":
            {
            "n1_4mer:GRV/n1_4mer:RV":{"val2":0.0,"val1":1.5,"mode":"absolute","condition":"bigger than"},
            "n0_4mer:D":{"val2":0.0,"val1":3.0,"mode":"sigmas","condition":"bigger than"},
            "n0_4mer:PS":{"val2":0.0,"val1":55.0,"mode":"absolute","condition":"bigger than"},
            },
            "Ribosomal proteins":
            {
            "n1_4mer:RV":{"val2":1.0,"val1":0.0,"mode":"sigmas","condition":"between"},
            "n0_4mer:D":{"val2":2.0,"val1":1.0,"mode":"sigmas","condition":"between"},
            "n0_4mer:PS":{"val2":0.0,"val1":55.0,"mode":"absolute","condition":"smaller than"}
            },
            "Giant genes":
            {
            "n1_4mer:RV":{"val2":0.0,"val1":1.2,"mode":"sigmas","condition":"bigger than"},
            "n0_4mer:D":{"val2":0.0,"val1":2.0,"mode":"sigmas","condition":"bigger than"}
            }
        }    
        self.save_scripts(path)

    def save_scripts(self,path=None):
        if path==None:
            path = os.path.join("lib","scripts.txt")
        output = ""
        for script_name in self.scripts:
            output += script_name + "\n"
            for task in self.scripts[script_name]:
                output += task + "\n"
                for name in self.scripts[script_name][task]:
                    output += name + ":" + str(self.scripts[script_name][task][name]) + ","
                output = output[:-1] + "\n"
            output += "END\n"
        self.IO.save(output,path)

    def select_scenario(self):
        options = self.scripts.keys()
        while 5 > 0:
            print()
            print("Select scenario")
            print()
            print("  0   Quit")
            for i in range(len(options)):
                print("  " + str(i+1) + "   " + options[i])
            print()
            index = input("Select scenario by the index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(options):
                    self.options["-c"] = options[index-1]
                    self.task_list = {}
                    self.task_list.update(self.scripts[self.options["-c"]])
                    return
                elif index == 0:
                    return
                else:
                    print()
                    print("\tWrong scenario index")
                    print()
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(options)))
                print()
                continue
        
    def add_scenario(self):
        print()
        script_name = input("Enter name for the current scenario: ")
        print()
        if not script_name:
            return
        self.scripts[script_name] = {}
        self.scripts[script_name].update(self.task_list)
        self.save_scripts()
        
    def remove_scenario(self):
        options = self.scripts.keys()
        while 5 > 0:
            print()
            print("Remove scenario")
            print()
            print("  0   Quit")
            for i in range(len(options)):
                print("  " + str(i+1) + "   " + options[i])
            print()
            index = input("Select scenario by the index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(options):
                    del self.scripts[options[index-1]]
                elif index == 0:
                    self.save_scripts()
                    return
                else:
                    print()
                    print("\tWrong scenario index")
                    print()
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(options)))
                print()
                continue

    def select_task(self,task_categories):
        index = None
        while index != 0:
            codes = []
            count = 1
            print()
            print("Select task category")
            print()
            print("  0   escape;")
            for task_code in task_categories.keys():
                print("  " + str(count) + "   " + task_code + " (" + task_categories[task_code] + ");")
                codes.append(task_code)
                count += 1
            print()
            index = input("Select the task category by its index: ")
            try:
                index = int(index)
                if index > 0 and index <= len(codes):
                    return codes[index-1]
                elif index != 0:
                    print()
                    print("\tWrong category index")
                    print()
                else:
                    break
            except:
                print()
                print("\tEnter an integer from 0 to " + str(len(codes)))
                print()
                continue
        return None

    def set_conditions(self,task,conditions):
        response = ""
        while response != "Q":
            print()
            print("Select condition for the task " + task)
            print()
            print("  G   Bigger than")
            print("  S   Smaller than")
            print("  B   Between")
            print("  M   " + conditions["mode"])
            print("  Q   quit;")
            try:
                response = input("?")
            except:
                continue
            response = string.upper(response)
            if response == "G":
                conditions["condition"] = "bigger than"
            elif response == "S":
                conditions["condition"] = "smaller than"
            elif response == "B":
                conditions["condition"] = "between"
            elif response == "M":
                if conditions["mode"] == "absolute":
                    conditions["mode"] = "sigmas"
                elif conditions["mode"] == "sigmas":
                    conditions["mode"] = "fraction"
                else:
                    conditions["mode"] = "absolute"
                continue
            elif response == "Q":
                return conditions
            else:
                continue
            
            question = conditions["condition"]
            while 5 > 0:
                print()
                try:
                    val = input(question + ": ")
                except:
                    print()
                    if question == "and":
                        print("Enter the upper limit of variation:")
                    else:
                        print("Enter the limit of variation:")
                    print()
                    continue
                try:
                    val = float(val)
                    if question == "and":
                        conditions["val2"] = val
                    else:
                        conditions["val1"] = val
                except:
                    print()
                    print("\tEnter a float point number")
                    print()
                    continue
                if conditions["condition"] != "between" or question == "and":
                    break
                question = "and"
            return conditions

    def main_menu(self):
        # Terminal interface
        if "MGE" in self.options["-c"]:
            self.task_list.update(self.scripts["MGE"])
        else:
            key = self.scripts.keys()[0]
            self.task_list.update(self.scripts[key])
        response = ''
        while response != "Q":
            tasks = list(self.task_list.keys())
            if not tasks:
                self.add_task()
                continue
            if len(tasks) > 1:
                tasks.sort()
            print("SeqWord Sniffer 3.0 2023/09/03")
            print()
            print("Settings for this run:\n")
            print("  C   Scenario?\t: " + self.options["-c"])
            print("  A   Add user defined scenario")
            print("  R   Remove scenario")
            print()
            if self.task_list[tasks[0]]["condition"] in ("bigger than","smaller than"):
                str_val = str(self.task_list[tasks[0]]["val1"])
            else:
                str_val = str(self.task_list[tasks[0]]["val1"]) + " and " + str(self.task_list[tasks[0]]["val2"])
            print(("  T   Tasks to perform?\t: " + tasks[0] +
                           "; " + self.task_list[tasks[0]]["condition"] +
                           " (" + self.task_list[tasks[0]]["mode"] + ") " + str_val))
            if len(tasks) > 1:
                for tn in range(1,len(tasks)):
                    if self.task_list[tasks[tn]]["condition"] in ("bigger than","smaller than"):
                        str_val = str(self.task_list[tasks[tn]]["val1"])
                    else:
                        str_val = str(self.task_list[tasks[tn]]["val1"]) + " and " + str(self.task_list[tasks[tn]]["val2"])
                    print(("\t\t\t: " + tasks[tn] +
                           "; " + self.task_list[tasks[tn]]["condition"] +
                           " (" + self.task_list[tasks[tn]]["mode"] + ") " + str_val))
                print()

            if self.options["-c"].find("MGE") > -1 or self.options["-c"] == "Ribosomal RNA":
                print("  U   Use BLASTn?\t: " + self.options["-u"])
            print("  L   Sliding window?\t: " + str(self.options["-l"]) + " bp.")
            print("  B   Big step?\t\t: " + str(self.options["-b"]) + " bp.")
            print("  M   Medium step?\t: " + str(self.options["-m"]) + " bp.")
            print("  S   Small step?\t: " + str(self.options["-s"]) + " bp.")
            print("  E   Refinement?\t: " + self.options["-e"])
            print("  I   Input folder?\t: " + self.options["-i"])
            print("  O   Output folder?\t: " + self.options["-o"])
            print("  F   Save sequences?\t: " + self.options["-f"])
            print("  V   Save SVG file?\t: " + self.options["-v"])
            print("  Q   to quit")
            print()
            print("Y to accept these settings, type the letter for one to change or Q to quit")
            response = input("?").upper()
            if response not in ("Y","L","T","O","I","S","M","C","A","R","F","E","B","V","U","N"):
                continue
            if response == "Y":
                print()
                self.execute()
            elif response == "C":
                self.select_scenario()
            elif response == "A":
                self.add_scenario()
            elif response == "R":
                self.remove_scenario()
            elif response == "T":
                self.edit_tasklist()
            elif response == "L":
                val = input("Length of sliding window? ")
                try:
                    self.options["-l"] = int(val)
                except:
                    print("Length of the sliding window must be an integer!")
                    print()
            elif response == "B":
                val = input("Big step of the sliding window? ")
                try:
                    self.options["-b"] = int(val)
                except:
                    print("Step of the sliding window must be an integer!")
                    print()
            elif response == "M":
                while 5 > 0:
                    val = input("Medium step of the sliding window? ")
                    try:
                        val = int(val)
                    except:
                        print()
                        print("Step of the sliding window must be an integer!")
                        print()
                        continue
                    if val < self.options["-s"] or val > self.options["-b"]:
                        print()
                        print("Medium step must be an integer between " + str(self.options["-s"]) + " and " + str(self.options["-b"]))
                        print()
                        continue
                    else:
                        self.options["-m"] = val
                        break
            elif response == "S":
                while 5 > 0:
                    val = input("Small step of the sliding window? ")
                    try:
                        val = int(val)
                    except:
                        print()
                        print("Step of the sliding window must be an integer!")
                        print()
                        continue
                    if val < 10 or val > self.options["-m"]:
                        print()
                        print("Small step must be an integer between 10 and " + str(self.options["-m"]))
                        print()
                        continue
                    else:
                        self.options["-s"] = val
                        break
                continue
            elif response == "I":
                folder = input("Input folder? ")
                if not os.path.exists(folder):
                    print("Folder " + folder + " doesn't exist!")
                    print()
                    continue
                self.options["-i"] = folder
            elif response == "O":
                folder = input("Output folder? ")
                for symbol in ("\\","/","|","<",">","*",":","\"","?"):
                    if string.find(folder,symbol)>-1:
                        print("Folder name must not containe symbol " + symbol)
                        print()
                        continue
                self.options["-o"] = folder
            elif response == "F":
                if self.options["-f"] == "no":
                    self.options["-f"] = "fasta"
                elif self.options["-f"] == "fasta":
                    self.options["-f"] = "gbk"
                elif self.options["-f"] == "gbk":
                    self.options["-f"] = "fasta+gbk"
                else:
                    self.options["-f"] = "no"
            elif response == "V":
                if self.options["-v"] == "no":
                    self.options["-v"] = "yes"
                else:
                    self.options["-v"] = "no"
            elif response == "E":
                if self.options["-e"] == "Contrasting/Iteration":
                    self.options["-e"] = "No"
                elif self.options["-e"] == "No":
                    self.options["-e"] = "Contrasting"
                elif self.options["-e"] == "Contrasting":
                    self.options["-e"] = "Iteration"
                elif self.options["-e"] == "Iteration":
                    self.options["-e"] = "Contrasting/Iteration"
            elif response == "U" and (self.options["-c"].find("MGE") > -1 or self.options["-c"] == "Ribosomal RNA"):
                if self.options["-u"] == "no":
                    self.options["-u"] = "yes"
                else:
                    self.options["-u"] = "no"
            elif response == "N" and self.options["-c"].find("MGE") > -1:
                if self.options["-n"] == "no":
                    self.options["-n"] = "yes"
                else:
                    self.options["-n"] = "no"
            else:
                continue

###############################################################################
    
# Validator
class Validator:
    def __init__(self):
        self.IO = seq_io.IO()
        self.options = {"-c":"",# name of the scenario
                   "-u":"",     # use BLASTN to filter rrn clusters
                   "-n":"",     # use BLASTN to search tRNA on the borders of GI
                   "-l":0,      # sliding window length
                   "-b":0,      # sliding window big step
                   "-m":0,      # sliding window medium step
                   "-s":0,      # sliding window small step
                   "-e":"",     # refinement No | Contrasting | Iteration | Contrasting/Iteration 
                   "-i":"",     # input folder
                   "-o":"",     # output folder
                   "-f":"",     # save GI sequemces: no/fasta/gbk/gbk+fasta
                   "-v":""      # save SVG file
                }
        
    def validate(self,options,scenario):
        for key in options:
            if key not in self.options:
                print()
                print("Wrong option "+str(key)+";")
                print()
                return 
            self.options[key] = options[key]
        if self.options["-c"] not in scenario:
            print()
            print("Wrong scenarium name "+str(self.options["-c"]))
            print("Must be in "+"|".join(scenario)+";")
            print()
            return 
        if self.options["-u"] not in ("yes","no"):
            print()
            print("Wrong option -u "+str(self.options["-u"]))
            print("Must be in yes|no;")
            print()
            return
        if self.options["-n"] not in ("yes","no"):
            print()
            print("Wrong option -n "+str(self.options["-n"]))
            print("Must be in yes|no;")
            print()
            return
        try:
            self.options["-l"] = int(self.options["-l"])
        except:
            print()
            print("Sliding window length -l must be an integer;")
            print()
            return
        try:
            self.options["-b"] = int(self.options["-b"])
        except:
            print()
            print("Sliding window step -b must be an integer;")
            print()
            return
        try:
            self.options["-m"] = int(self.options["-m"])
        except:
            print()
            print("Sliding window step -m must be an integer;")
            print()
            return
        try:
            self.options["-s"] = int(self.options["-s"])
        except:
            print()
            print("Sliding window step -s must be an integer;")
            print()
            return
        if self.options["-s"] > self.options["-m"]:
            print()
            print("Small window step -s must be smaller or equal to the medium step -m;")
            print()
            return
        if self.options["-m"] > self.options["-b"]:
            print()
            print("Medium window step -m must be smaller or equal to the big step -b;")
            print()
            return
        if self.options["-b"] > self.options["-l"]:
            print()
            print("Big window step -b must be smaller or equal to the window length -l;")
            print()
            return
        if self.options["-e"] not in ("No","Contrasting","Iteration","Contrasting/Iteration"):
            print()
            print("Wrong option -e "+str(self.options["-e"]))
            print("Must be in "+"|".join(["No","Contrasting","Iteration","Contrasting/Iteration"])+";")
            print()
            return
        if not os.path.exists(self.options["-i"]):
            print()
            print("Folder -i %s does not exists;" % str(self.options["-i"]))
            print()
            return
        if not os.path.exists(self.options["-o"]):
            folder_name = self.IO.new_folder(self.options["-o"])
            if folder_name == None or not os.path.exists(folder_name):
                print()
                print("Folder -o %s does not exists and can not be created;" % str(self.options["-o"]))
                print()
                return
            self.options["-o"] = folder_name
        if self.options["-f"] not in ("no","fasta","gbk","fasta+gbk"):
            print()
            print("Wrong option -f "+str(self.options["-f"]))
            print("Must be in "+"|".join(["no","fasta","gbk","gbk+fasta"])+";")
            print()
            return
        if self.options["-v"] not in ("yes","no"):
            print()
            print("Wrong option -v "+str(self.options["-v"]))
            print("Must be in yes|no;")
            print()
            return
        return self.options
                
###############################################################################

if __name__ == "__main__":
    oInterface = Interface()
    
