import sys, os, string, array, math, time, shelve, subprocess
import blast, seq_io

class Main:
    def __init__(self,task_list,options):
        
        # ATTRIBUTES
        self.IO = seq_io.IO()
        self.task_list = {}
        self.task_list.update(task_list)
        self.frame = options["-l"]
        self.bigstep = options["-b"]
        self.mediumstep = options["-m"]
        self.smallstep = options["-s"]
        self.input_path = options["-i"]
        self.output_path = options["-o"]
        self.flg_saveFastaFile = options["-f"]
        self.flg_saveSVGimage = options["-v"]
        self.flg_doBLAST = options["-u"]
        self.flg_search_tRNA = options["-n"]
        self.scenario = options["-c"]
        self.echo = None
        if options["-e"] == "Contrasting/Iteration":
            self.flg_contrasting = True
            self.flg_iteration = True
        elif options["-e"] == "No":
            self.flg_contrasting = False
            self.flg_iteration = False
        elif options["-e"] == "Contrasting":
            self.flg_contrasting = True
            self.flg_iteration = False
        elif options["-e"] == "Iteration":
            self.flg_contrasting = False
            self.flg_iteration = True
            
        if self.flg_doBLAST == "yes" and self.scenario.find("MGE") > -1:
            self.flg_doBLAST = False
        elif self.flg_doBLAST == "yes" and self.scenario == "Ribosomal RNA":
            self.flg_doBLAST = True
        else:
            self.flg_doBLAST = None
        
        self.StartTime = None
        self.DataSet = None
        self.StandardPatterns = {}
        self.oSVG = None

        # Validator
        self.validator = Validator()
        self.tasks = self.validator.validateTasks(task_list,self.frame)

        # CONSTANTS
        self.symbols = {'GC-content':' %;',
                        'G/C-skew':';',
                        'A/T-skew':';',
                        'Pattern deviation':' %;',
                        'Absolute deviation':' %;',
                        'Pattern skew':' %;',
                        'Variance':';'}

        # IMPLEMENTATION
        self.initiateDataSet()
        self.process()

    # METHODS

    def initiateDataSet(self):
        self.DataSet = {'Sequence name':'',
                        'Sequence link':'',
                        'Sequence description':'',
                        'Total sequence length':0,
                        'Locus length':0,
                        'Left border':0,
                        'Frame':0,
                        'Step':0,
                        'Time':'',
                        'Tasks':{"gc":{}}
                        }
        
    # process source data files
    def process(self):
        if not self.tasks:
            return
        # get file list
        filelist = self.get_fileList(self.input_path)
        if len(filelist) == 0:
            return None
        # run cicle for all files in the list
        outdata = ''
        tl = {}
        for task in self.task_list:
            tl[task] = {}
            tl[task].update(self.task_list[task])
        for fname in filelist:
            #StartTime = time.clock()
            if len(fname) <= 4:
                continue
            if fname[-4:].upper() == ".GBK":
                seqlist = self.getSequenceFromGBK(fname)
            elif fname[-5:].upper() == ".GBFF":
                seqlist = self.getSeqFromGBFF(fname)
            else:
                seqlist = self.getSeqFromFASTA(fname)
            if not seqlist:
                return
            acc = ""
            for seqname in seqlist:
                genes = {}
                if seqlist[seqname]['dataset'] and seqlist[seqname]['dataset']['Gene map']:
                    genes.update(seqlist[seqname]['dataset']['Gene map'])
                seq = seqlist[seqname]["sequence"]
                acc = "GI%s" % self.IO.random_id(4)
                if seqlist[seqname]["dataset"] and seqlist[seqname]["dataset"]["Accession"]:
                    acc = seqlist[seqname]["dataset"]["Accession"]
                if len(seq) == 0:
                    continue
                print(seqname)
                                
                self.get_output(seq,seqname,self.bigstep,genes)
                self.oSVG = None
                if self.flg_saveSVGimage == "yes":
                    self.oSVG = SVG(seqname,len(seq),self.task_list,acc)
                if self.oSVG and "GC" not in self.task_list and "AT" not in self.task_list:
                    task_description = {"statistics":[50.0,5.0],"condition":""}
                    self.oSVG.add_task("GC-content",task_description,self.DataSet["Tasks"]["gc"])
                self.evaluateTasks()

                if self.flg_iteration:
                    borders = self.getPerturbationBorders()
                    if borders:
                        short_seq = seq[:borders[0][0]]
                        for i in range(1,len(borders),1):
                            short_seq += seq[borders[i-1][1]:borders[i][0]]
                        short_seq += seq[borders[-1][1]:]
                        self.updateStdPatterns(short_seq)
                        if self.oSVG:
                            self.oSVG.clear_tasks()
                        self.get_output(seq,seqname,self.bigstep,genes,False,False)
                        if self.oSVG and "GC" not in self.task_list and "AT" not in self.task_list:
                            task_description = {"statistics":[50.0,5.0],"condition":""}
                            self.oSVG.add_task("GC-content",task_description,self.DataSet["Tasks"]["gc"])
                        self.evaluateTasks()
                        
                output = self.getGenomicIslands(seq)
                outfilename = "_".join([self._generic_name(seqname,acc),self.scenario])
                for symbol in ("/","\\",":","*","?","<",">""|"):
                    outfilename = outfilename.replace(symbol,"_")
                mge_counter = 0
                if len(output):
                    print(str(len(output)) + " candidate islands have been predicted!")
                    fasta = []
                    stat = []
                    outStr = ""
                    gi_counter = 1
                    for coord in output:
                        flg_tRNA = False
                        lb, rb, stat_data = coord
                        if self.flg_search_tRNA == "yes":
                            if self.scenario.find("MGE") > -1 and rb-lb > 6000:
                                start, stop, gi_seq = self.get_GISeq(lb-500,rb,seq)
                                oBlast = blast.BLAST("blastn","bin","trna_bact trna_arch",gi_seq)
                                oBlast.execute()
                                result = oBlast.get_close_to_point(500)
                                if result:
                                    flg_tRNA = True
                                    lb -= min(result[0],result[1])
                                start, stop, gi_seq = self.get_GISeq(lb,rb+500,seq)
                                oBlast = blast.BLAST("blastn","bin","trna_bact trna_arch",gi_seq)
                                oBlast.execute()
                                result = oBlast.get_close_to_point(500)
                                if result:
                                    flg_tRNA = True
                                    rb += min(result[0],result[1])
                        
                        stat.append("; ".join(map(lambda key: "%s = %f" % (key,stat_data[key]),stat_data.keys())))

                        if self.flg_doBLAST != None:
                            if int(lb) <= 0:
                                query = seq[int(lb):]+seq[:int(rb)]
                            elif rb >= len(seq):
                                query = seq[int(lb):]+seq[:int(rb)-len(seq)+1]
                            else:
                                query = seq[int(lb)-1:int(rb)-1]
                            
                            oBlast = blast.BLAST("blast","dna","bin","16S",query)
                            oBlast.execute()
                            result = oBlast.get_top_alignment()
                            if self.flg_doBLAST == False and result and result[0] == 0 and result[1] >= 500:
                                if self.oSVG:
                                    self.oSVG.add_gi(lb,rb,"grey","rrn operon %d-%d" % (lb,rb))
                                continue
                            elif self.flg_doBLAST == True and not result and (result[0] != 0 or result[1] < 500):
                                continue
                        mge_counter += 1
                        if self.oSVG:
                            color = ""
                            if flg_tRNA:
                                color = "red"
                            self.oSVG.add_gi(lb,rb,color,"%s:%d [%d-%d];" % (acc,mge_counter,lb,rb))
                        if int(lb) <= 0:
                            lb = 1
                        if seqlist[seqname]['dataset'] and seqlist[seqname]['dataset']['Gene map']:
                            annotation,lb,rb = self.getAnnotation(seqlist[seqname]['dataset']['Gene map'],lb,rb)
                        else:
                            annotation = ""
                        outStr += "<GI> %s:%d <COORDINATES> %d-%d <STAT> %s\n%s\n<END>\n" % (acc,gi_counter,lb,rb,stat[-1],annotation)
                        if self.flg_saveFastaFile != "No":
                            lb, rb, gi_seq = self.get_GISeq(lb,rb,seq)
                            fasta.append(">%s:%d|%s [%s-%s]" % (acc,gi_counter,seqname,lb,rb))
                            fasta.append(gi_seq)
                        gi_counter += 1
                    if self.flg_saveFastaFile == "fasta" or self.flg_saveFastaFile == "fasta+gbk":
                        self.IO.save("\n".join(fasta),os.path.join(self.output_path,outfilename+".fas"))
                    if genes and (self.flg_saveFastaFile == "gbk" or self.flg_saveFastaFile == "fasta+gbk"):
                        for i in range(0,len(fasta),2):
                            seq_id,seqname = fasta[i].split("|")
                            locus_name = seq_id[1:]
                            gbk_fname = os.path.join(self.output_path,locus_name.replace(":","_")+".gbk")
                            gi_seq = fasta[i+1]
                            try:
                                start,stop = list(map(lambda s: int(float(s)),seqname[seqname.find("[")+1:-1].split("-")))
                            except:
                                start,stop = list(map(lambda s: int(float(s)),seqname[seqname.find("[")+1:-1].split("..")))
                            self.IO.saveGBK(gbk_fname,gi_seq,start,stop,locus_name,stat[i//2])
                    
                    self.IO.save(outStr[:-1],os.path.join(self.output_path,outfilename+".out"))
                    print(str(mge_counter) + " genomic islands were identified!")
                else:
                    print("No genomic islands were identified!")
                if self.oSVG:
                    self.save_svg(outfilename)
            #print("has been done in " + str(time.clock() - StartTime) + " sec.")
            print()
            for task in tl:
                self.task_list[task] = {}
                self.task_list[task].update(tl[task])
        return
    
    def get_GISeq(self,lb,rb,seq,flg_circular=True):
        before = after = ""
        start = lb
        end = rb
        if lb < 1 and not flg_circular:
            lb = start = 1
        elif lb < 1 and flg_circular:
            start = 1
            before = seq[len(seq)+lb:]
        elif lb==0:
            lb = start = 1
            
        if rb > len(seq) and not flg_circular:
            rb = end = len(seq)
        elif rb > len(seq) and flg_circular:
            after = seq[:rb-len(seq)]
            end = len(seq)
        GISeq = "".join([before,seq[int(start)-1:int(end)],after])
        return lb,rb,GISeq
    
    def evaluateTasks(self,flg_updateStat=True,flg_generateSVG=True):
        for task in self.task_list:
            if flg_updateStat:
                self.task_list[task]["statistics"] = {}
            values = []
            # append values of pattern deviations
            taskId = self.tasks[task]['ID']
            subtr = self.tasks[task]['subtr']
            divisorId = self.tasks[task]['divisor']
            windows = {}
            for win in self.DataSet["Tasks"][taskId]:
                if subtr:
                    subtrahend = self.DataSet["Tasks"][subtr][win]['value']
                else:
                    subtrahend = 0
                if divisorId:
                    divisor = self.DataSet["Tasks"][divisorId][win]['value']
                else:
                    divisor = 1.0
                if divisor == 0:
                    curr_value = self.DataSet["Tasks"][taskId][win]['value']-subtrahend
                else:
                    curr_value = (self.DataSet["Tasks"][taskId][win]['value']-subtrahend)/divisor
                values.append(curr_value)
                windows[win] = curr_value
            if not values:
                continue
            if self.oSVG and flg_generateSVG:
                task_data = {}
                task_data.update(self.task_list[task])
                task_data["statistics"] = self.getMeanAndStDev(values)
                self.oSVG.add_task(task,task_data,windows)
            # calculate average and std. dev. for pattern deviations
            if self.task_list[task]["mode"] == "absolute":
                continue
            elif self.task_list[task]["mode"] == "sigmas": 
                if not self.task_list[task]["statistics"] and flg_updateStat:
                    self.task_list[task]["statistics"] = self.getMeanAndStDev(values)
            else:
                values.sort()
                if self.task_list[task]['val1']:
                    index = int(self.task_list[task]['val1']*len(values))
                    self.task_list[task]['val1'] = values[index]
                if self.task_list[task]['val2']:
                    index = int(self.task_list[task]['val2']*len(values))
                    self.task_list[task]['val2'] = values[index]
    
    def getAnnotation(self,data,lb,rb):
        annotation = "\t"
        coordinates = list(data.keys())
        if not coordinates:
            return
        try:
            coordinates.sort(key=lambda s: int(s.split("-")[0]),reverse=True)
        except:
            coordinates.sort(key=lambda s: int(s.split("..")[0]),reverse=True)
        #coordinates.sort(self.sort_by_coordinates)
        for coord in coordinates:
            if ((data[coord]['stop'] >= lb and data[coord]['stop'] <= rb) or
                (data[coord]['start'] >= lb and data[coord]['start'] <= rb)):
                description = [data[coord]['name'],data[coord]['description'],data[coord]['remark']]
                for i in range(2,-1,-1):
                    if not description[i]:
                        del description[i]
                if description:
                    description = "; ".join(description)
                else:
                    description = ""
                annotation += ("[" + str(data[coord]['start']) + ":" +
                    str(data[coord]['stop'])+":"+str(data[coord]['direction']) + "] " +
                    description + "\n\t")
                if data[coord]['start'] < lb:
                    lb = data[coord]['start']
                if data[coord]['stop'] > rb:
                    rb = data[coord]['stop']
            elif data[coord]['start'] > rb:
                break
            else:
                pass
        return annotation[:-2],lb,rb

    def getSeqFromFASTA(self,fname,names=[]):
        result = self.IO.openFasta(fname)
        if not result:
            return
        seqlist,path = result
        key = seqlist.keys()[0]
        return {key:{"sequence":seqlist[key],"dataset":None}}

    def getSequenceFromGBK(self,fname):
        result = self.IO.openGBK(fname,"ALL")
        if not result:
            return
        dataset,seq,path = result
        return {dataset['Sequence name']:{"sequence":seq,"dataset":dataset}}

    # check if the proposed fname already exists in the database
    # if so, change the name or return the proposed name    
    def getSeqFromGBFF(self,fname,names=[]):
        try:
            f = open(fname,"r")
            data = f.read()
            f.close()
        except:
            return {}
        data = data.split("\n")
        seqlist = {}
        seqname = ""
        sequence = ""
        counter = 0
        END_COUNT = len(data)
        while counter < END_COUNT:
            if not data[counter]:
                counter += 1
                continue
            if len(data[counter]) >= 12 and data[counter][:12] == "DEFINITION  ":
                seqname = data[counter][12:]
                if names and seqname not in names:
                    continue
                i = 1
                if seqname in seqlist:
                    seqname += " #" + str(i)
                while seqname in seqlist:
                    seqname = seqname[:-len(str(i))] + str(i+1)
                    i += 1
                seqlist[seqname] = ""
            if len(data[counter]) >= 12 and data[counter][:12] == "ORIGIN      ":
                counter += 1
                while data[counter] != "//":
                    line = data[counter][10:]
                    line = line.replace(" ","").upper()
                    seqlist[seqname] += line
                    counter += 1
                seqname = ""
            counter += 1
        return seqlist
        
    def parse_seqname(self,line,seq_names):
        if line > 3 and line[:3]=="gi|":
            seqname_elements = line.split("|")
            if len(seqname_elements)==5:
                seqname = seqname_elements[4]
                pos = seqname.find(", complete sequence")
                if pos > -1:
                    seqname = seqname[:pos]
                if seqname_elements[3].find("NC_")==0:
                    accession = seqname_elements[3]
                    pos = accession.find(".")
                    if pos > -1:
                        accession = accession[:pos]
                    seqname += " [" + accession + "]"
            else:
                seqname = seqname_elements[0]
        else:
            seqname = line
        seqname = self.check_seqname(seqname,seq_names)
        return seqname
        
    def check_seqname(self,seqname,names=[]):
        mark = seqname.rfind(", complete")
        if mark > 0:
            seqname = seqname[:mark]
        for symb in ("\\","/",":","*","?","\"","<",">","|"):
            seqname = seqname.replace(symb," ")
        while seqname and seqname[0] == " ":
            seqname = seqname[1:]
        pos = seqname.find(", complete")
        if pos != -1:
            seqname = seqname[:pos]
        while seqname and (seqname[-1] == "." or seqname[-1]) == " ":
            seqname = seqname[:-1]
        if not seqname:
            seqname = "#1"
        if not names:
            return seqname
        while seqname.upper() in names:
            sep = seqname.rfind("#")+1
            if sep == 0:
                seqname += " #1"
            else:
                try:
                    n = int(seqname[sep:])
                except:
                    return seqname + " #1"
                seqname = seqname[:sep] + str(n+1)
        return seqname
       
    def getGenomicIslands(self,seq):
        print("\tlooking for genomic islands")
        output = []
        borders = self.getPerturbationBorders()
        self.initiateDataSet() 
        for lb,rb,stat_data in borders:
            start,stop = self.setBorders(lb,rb,self.bigstep,len(seq))
            self.get_output(seq[start:stop],"",self.mediumstep,None,False)
            self.evaluateTasks(False,False)
            key = list(self.DataSet['Tasks'].keys())[0]
            coord = list(self.DataSet['Tasks'][key].keys())
            try:
                coord.sort(key=lambda s: int(s.split("-")[0]),reverse=True)
            except:
                coord.sort(key=lambda s: int(s.split("..")[0]),reverse=True)
            #coord.sort(self.sort_by_coordinates)
            innerBorders = self.getPerturbationBorders()
            if not innerBorders:
                innerBorders = [[0,rb-lb,stat_data]]
            left = innerBorders[0][0] + lb
            if left+self.frame <= len(seq)-1:
                right = left+self.frame
            else:
                right = len(seq)-1
            if right-left >= self.frame:
                innerBorders[0][0] = self.getExactValue(left,seq[left:right],0)
            right = innerBorders[0][1] + lb + 3*self.frame
            if right - self.frame >= 0:
                left = right - self.frame
            else:
                left = 0
            if right-left >= self.frame:
                innerBorders[0][1] = self.getExactValue(left,seq[left:right],1)
            start,stop,stat_data = innerBorders[0]
            output.append([start,stop,stat_data])
            self.initiateDataSet()
            # concatenate overlapping islands
            if len(output) > 1:
                for i in range(len(output)-1,0,-1):
                    if output[i][0] <= output[i-1][1]:
                        output[i-1] = self.join_regions(output[i-1],output[i])
                        del output[i]
        return output
    
    def join_regions(self,first,second):
        box = [first[0],second[1],{}]
        for key in first[2]:
            a = float(first[1]-first[0])
            b = float(second[1]-second[0])
            c = float(first[1]-first[0]+second[1]-second[0])
            box[2][key] = (first[2][key]*a + second[2][key]*b)/c
        return box

    def getExactValue(self,left,seq,mode):
        if len(seq) < self.frame:
            return left
        self.initiateDataSet()
        self.get_output(seq,"",self.smallstep,None,False)
        key = list(self.DataSet["Tasks"].keys())[0]
        coordinates = list(self.DataSet["Tasks"][key].keys())
        try:
            coordinates.sort(key=lambda s: int(s.split("-")[0]),reverse=True)
        except:
            coordinates.sort(key=lambda s: int(s.split("..")[0]),reverse=True)
        #coordinates.sort(self.sort_by_coordinates)
        for i in range(1,len(coordinates)):
            isSignal = 1
            for task in self.tasks:
                subtr = self.tasks[task]['subtr']
                if subtr:
                    subtrahend = self.DataSet["Tasks"][subtr][coordinates[i]]['value']
                else:
                    subtrahend = 0
                divisorId = self.tasks[task]['divisor']
                if divisorId:
                    divisor = self.DataSet["Tasks"][divisorId][coordinates[i]]['value']
                else:
                    divisor = 1.0
                taskId = self.tasks[task]['ID']
                value = self.DataSet["Tasks"][taskId][coordinates[i]]['value']
                if divisor==0 or not self.check_condition(self.task_list[task],(value-subtrahend)/divisor):
                    isSignal = 0
                    break
            if (mode and not isSignal) or (not mode and isSignal):
                return left + i*self.smallstep - self.frame/2
            elif mode:
                return left + len(coordinates)*self.smallstep - self.frame/2
            else:
                return left + self.frame/2
                     
    def getPerturbationBorders(self,flg_report=False):
        data = {}
        data.update(self.DataSet)
        borders = []
        inIsland = False
        key = list(data["Tasks"].keys())[0]
        coordinates = list(data["Tasks"][key].keys())
        try:
            coordinates.sort(key=lambda s: int(s.split("-")[0]),reverse=True)
        except:
            coordinates.sort(key=lambda s: int(s.split("..")[0]),reverse=True)
        #coordinates.sort(self.sort_by_coordinates)
        tk = data["Tasks"].keys()
        stat_data = {}
        output = ""
        for i in range(1,len(coordinates)):
            report = str(coordinates[i])+"\n"
            isSignal = True
            for task in self.tasks:
                subtr = self.tasks[task]['subtr']
                if subtr:
                    subtrahend = data["Tasks"][subtr][coordinates[i]]['value']
                else:
                    subtrahend = 0
                divisorId = self.tasks[task]['divisor']
                if divisorId:
                    divisor = data["Tasks"][divisorId][coordinates[i]]['value']
                else:
                    divisor = 1.0
                taskId = self.tasks[task]['ID']
                value = float((data["Tasks"][taskId][coordinates[i]]['value']-subtrahend)/divisor)
                if task not in stat_data:
                    stat_data[task] = []
                stat_data[task].append(value)
                report += task + ": " + str(value)+"; "
                if divisor==0 or not self.check_condition(self.task_list[task],value):
                    isSignal = False
                    report += "%f, %s, %f;" % self.get_condition(self.task_list[task],value)
                    break
            if isSignal and not inIsland:
                inIsland = True
                try:
                    if i:
                        lb = int(coordinates[i-1].split("-")[0])-self.frame
                    else:
                        lb = int(coordinates[i].split("-")[0])-self.frame
                except:
                    if i:
                        lb = int(coordinates[i-1].split("..")[0])-self.frame
                    else:
                        lb = int(coordinates[i].split("..")[0])-self.frame
                output += report + "\nstart\n"
            elif not isSignal and inIsland:
                inIsland = False
                try:
                    if i == len(coordinates)-1:
                        rb = int(coordinates[i].split("-")[1])-self.frame
                    else:
                        rb = int(coordinates[i+1].split("-")[1])-self.frame
                except:
                    if i == len(coordinates)-1:
                        rb = int(coordinates[i].split("..")[1])-self.frame
                    else:
                        rb = int(coordinates[i+1].split("..")[1])-self.frame
                output += report + "\nstop\n\n"
                borders.append([lb,rb,{}])
                for key in stat_data:
                    stat_data[key] = float(sum(stat_data[key]))/len(stat_data[key])
                borders[-1][2].update(stat_data)
                stat_data = {}
            elif not isSignal and not inIsland:
                stat_data = {}
                output += report + "\npass\n\n"
        # Saving the report
        if flg_report:
            self.IO.save(output,"report.txt")
        if len(borders) > 1:
            borders.sort()
        return borders
    
    def check_condition(self,data,val):
        if data["mode"] == "sigmas":
            mean,stdev = data["statistics"]
            if stdev:
                val = (val-mean)/stdev
        if data["condition"] == "bigger than" and val > data["val1"]:
            return 1
        elif data["condition"] == "smaller than" and val < data["val1"]:
            return 1
        elif data["condition"] == "between" and val > data["val1"] and val < data["val2"]:
            return 1
        else:
            return 0
        
    def get_condition(self,data,val):
        if data["mode"] == "sigmas":
            mean,stdev = data["statistics"]
            if stdev:
                val = (val-mean)/stdev
        return val,data["condition"],data["val1"]
        
    # return list of ".fst" files in the home directory
    def get_fileList(self,home_folder):
        extensions = ('FNA','FAS','FST','FASTA',"GBK","GBFF")
        l = os.listdir(home_folder)
        filelist = []
        for fname in l:
            if os.path.isfile(os.path.join(home_folder,fname)):
                if len(fname)>4 and fname[-3:].upper() in extensions:
                    filelist.append(os.path.join(home_folder,fname))
                else:
                    continue
            elif os.path.isdir(os.path.join(home_folder,fname)):
                local_files = self.get_fileList(os.path.join(home_folder,fname))
                for item in local_files:
                    filelist.append(item)
            else:
                pass
        return filelist
    
    def setBorders(self,lb,rb,step,seqlength):
        lb -= step
        if lb < 0:
            lb = 0
        rb += step
        if rb >= seqlength:
            rb = seqlength-1
        return lb,rb
    
    # return requested data table for the sequence in the given file
    def get_output(self, seq, name, step, genes, flg_update_references=True, flg_renew_dataset=True):
        #self.StartTime = time.clock()        
        seqLen = len(seq)
        
        CurrTime = self.getCurrTime()
        # dictionary of standard patterns
        self.DataSet['Time'] = CurrTime
        if flg_renew_dataset:
            self.DataSet['Sequence name'] = name
            self.DataSet['Total sequence length'] = seqLen
            self.DataSet['Locus length'] = seqLen
            self.DataSet['Left border'] = 0
            self.DataSet['Frame'] = self.frame
            self.DataSet['Step'] = step
            for currTask in self.tasks:
                if self.tasks[currTask]['subtr']:
                    self.DataSet["Tasks"][self.tasks[currTask]['subtr']] = {}
                if self.tasks[currTask]['divisor']:
                    self.DataSet["Tasks"][self.tasks[currTask]['divisor']] = {}
                self.DataSet["Tasks"][self.tasks[currTask]['ID']] = {}
        
        if flg_update_references:
            self.StandardPatterns = {}
            for currTask in self.tasks:
                task = {}
                task.update(self.tasks[currTask])
                if self.tasks[currTask]['ID'] not in self.StandardPatterns:
                    self.StandardPatterns[self.tasks[currTask]['ID']] = None
                if self.tasks[currTask]['subtr'] and self.tasks[currTask]['subtr'] not in self.StandardPatterns:
                    self.StandardPatterns[self.tasks[currTask]['subtr']] = None
                if self.tasks[currTask]['divisor'] and self.tasks[currTask]['divisor'] not in self.StandardPatterns:
                    self.StandardPatterns[self.tasks[currTask]['divisor']] = None
            
        # maximal number of digits in coordinates
        maxnumlen = len(str(self.DataSet['Total sequence length']))
        output = ''
        start = 0
        stop = self.frame-1
        while stop <= seqLen+self.frame:
            if self.echo:
                print([start,stop])
            if stop > len(seq):
                locus = seq[start:] + seq[:stop-len(seq)-1]
                output += str(start) + "\t" + str(stop-len(seq))
            else:
                locus = seq[start:stop]
                output += str(start) + "\t" + str(stop+1)
            if not locus:
                break
            key = (maxnumlen - len(str(start)))*" " + str(start)+'-'+str(stop)         
            # GC-content
            self.DataSet["Tasks"]["gc"][key] = 100.0*float(locus.upper().count("G")+locus.upper().count("C"))/len(locus)

            for currTask in self.tasks:
                task = {}
                task.update(self.tasks[currTask])
                if self.echo:
                    TaskStartTime = time.clock()
                    print(task)
                # define a standard pattern for a new task
                if task['task'] in ("D","GD","GPS","GRPS","GRV") and self.StandardPatterns[task['ID']] == None:
                    self.StandardPatterns[task['ID']] = self.getPattern(seq,task['wlength'],task['norm'],task['type'],genes)
                oPattern = None
                if key in self.DataSet["Tasks"][task['ID']] and self.DataSet["Tasks"][task['ID']][key]:
                    oPattern = self.DataSet["Tasks"][task['ID']][key]['oup']
                value,oPattern = self.get_value(locus,
                                       self.validator.tasks[task['task']],
                                       task['norm'],
                                       task['wlength'],
                                       task['type'],
                                       self.StandardPatterns[task['ID']])
                try:
                    value = float(value)
                except:
                    value = float(value[:-1])
                DataSet = {'value':value,
                           'start':start,
                           'stop':stop,
                           'oup':oPattern}
                #key = (maxnumlen - len(str(start)))*" " + str(start)+'-'+str(stop)
                self.DataSet["Tasks"][task['ID']][key] = {}
                self.DataSet["Tasks"][task['ID']][key].update(DataSet)
                
                for id in (task['subtr'],task['divisor']):
                    if not id or (id in self.DataSet["Tasks"] and key in self.DataSet["Tasks"][id]):
                        continue
                    tsk = id[8:]
                    ptype = id[0]
                    norm = int(id[1])
                    wlength = int(id[3])
                    if tsk in ("D","GD","GPS","GRPS","GRV") and not self.StandardPatterns[id]:
                        self.StandardPatterns[id] = self.getPattern(seq,wlength,norm,ptype,genes)
                    oPattern = None
                    if key in self.DataSet["Tasks"][id] and self.DataSet["Tasks"][id][key]:
                        oPattern = self.DataSet["Tasks"][id][key]['oup']
                    value,oPattern = self.get_value(locus,
                                           self.validator.tasks[tsk],
                                           norm,
                                           wlength,
                                           ptype,
                                           self.StandardPatterns[id])
                    try:
                        value = float(value)
                    except:
                        value = float(value[:-1])
                    DataSet = {'value':value,
                               'start':start,
                               'stop':stop,
                               'oup':oPattern}
                    self.DataSet["Tasks"][id][key] = {}
                    self.DataSet["Tasks"][id][key].update(DataSet)
                
            start += step
            stop += step
            if self.echo:
                print("has been done in " + str(time.clock() - TaskStartTime) + " sec.")
                print()
        return output
        
    def get_value(self, locus, mode, normalization, wlength, pattern_type, oStdPattern=None, oCurrPattern=None):
        # create local pattern
        if not oCurrPattern:
            oCurrPattern = self.getPattern(locus,wlength,normalization,pattern_type)
        oSubtrPattern = None
        if mode == 'GC-content':
            return '%s' % oCurrPattern.getPercentage("GC")
        elif mode == 'G/C-skew':
            return '%s' % oCurrPattern.getGCskew()
        elif mode == 'A/T-skew':
            return '%s' % oCurrPattern.getATskew()
        elif mode in ('Generalized distance',
                    'Generalized pattern skew',
                    'Generalized relative pattern skew',
                    'Generalized variance',
                    'Generalized relative variance') and normalization:
            oNormalizationTable = oStdPattern.getNormalizationTable()
            oCurrPattern.setNormalizationTable(oNormalizationTable)
        else:
            pass
            
        if mode in ('Distance', 'Generalized distance'):
            return oCurrPattern-oStdPattern,oCurrPattern
        
        elif mode in ('Pattern skew',
                      'Relative pattern skew',
                      'Generalized pattern skew',
                      'Generalized relative pattern skew'):
            return '%s' % oCurrPattern.getPS(),oCurrPattern
        
        elif mode in ('Variance',
                      'Generalized variance',
                      'Relative variance',
                      'Generalized relative variance'):
            return '%s' % oCurrPattern.getOUV(),oCurrPattern
        
        else:
            print('Error mode ' + mode)
            return None
    
    def getPattern(self,seq,wlength,norm,ptype,genes=None):
        if genes:
            seq = self.getCodingSequence(seq,genes)
        oPattern = Pattern(wlength)
        oPattern.setPattern(seq, norm, ptype)
        return oPattern
    
    def getCodingSequence(self,seq,genes):
        if not self.flg_contrasting or not genes:
            return seq
        cseq = ""
        for gene in genes:
            if (genes[gene]['name']+genes[gene]['description']+genes[gene]['remark']).find('hypothetical') > -1:
                continue
            try:
                lb,rb = map(lambda v: int(v),gene.split(".."))
            except:
                lb,rb = map(lambda v: int(v),gene.split("-"))
            cseq += seq[lb:rb]
        if len(cseq) < 10000:
            print()
            print("\tCoding sequence < 10 000 bp; initial sequence was returned!")
            print()
            return seq
        return cseq
    
    def updateStdPatterns(self,seq):
        for key in self.StandardPatterns:
            ptype = key.split(":")[0]
            norm = int(ptype[1])
            wlength = int(ptype[3])
            ptype = ptype[0]
            self.StandardPatterns[key] = self.getPattern(seq,wlength,norm,ptype)
    
    def getMeanAndStDev(self,DataList):
        sum = 0
        sum_sq = 0
        N = len(DataList)
        for i in range(N):
            sum = sum + DataList[i]
            sum_sq = sum_sq + DataList[i] * DataList[i]
        if N < 2:
            return [0,1.0]
        val = (sum_sq - sum * sum/N)/(N - 1)
        stdev = math.sqrt(val)
        if stdev == 0:
            stdev = 1.0
        mean = sum/N
        return [mean,stdev]
    
    def _generic_name(self,seqname,acc):
        cgr = seqname.rfind(", complete")
        if cgr > -1:
            seqname = seqname[:cgr]
        seqname = seqname.replace(" ","_")
        return "%s_[%s]" % (seqname,acc)
        
    def getCurrTime(self):
        CurrTime = time.localtime()
        CurrTime = (str(CurrTime[3]) + ":" + str(CurrTime[4]) + ":" + str(CurrTime[5]) + " " +
                str(CurrTime[2]) + "." + str(CurrTime[1]) + "." + str(CurrTime[0]))
        return CurrTime

    def save_svg(self,seqname):
        oIO = seq_io.IO()
        oIO.save(self.oSVG.get_svg(),os.path.join(self.output_path,seqname+".svg"))

    # SORTING METHODS
    def sort_by_words(self, a, b):
        if a[0] < b[0]:
            return -1
        elif a[0] > b[0]:
            return 1
        return 0

    def sort_by_order(self, a, b):
        if b[1] < a[1]:
            return -1
        elif b[1] > a[1]:
            return 1
        return 0
    
    def sort_by_values(self, a, b):
        if b[2] < a[2]:
            return -1
        elif b[2] > a[2]:
            return 1
        else:
            if b[0] < a[0]:
                return -1
            elif b[0] > a[0]:
                return 1
            else:
                return 0
    def sort_by_coordinates(self, a, b):
        try:
            valA = int(a.split("..")[0])
        except:
            valA = int(a.split("-")[0])
        try:
            valB = int(b.split("..")[0])
        except:
            valB = int(b.split("-")[0])
        if valA < valB: return -1
        elif valA > valB: return 1
        else: return 0
        
##############################################################################################################
class SVG:
    def __init__(self,seqname,seqlength,task_list,title=""):
        self.seqlength = seqlength
        self.seqname = seqname
        self.task_list = task_list
        self.title = title
        self.flg_wrongLoci = False
        self.added_tasks = 0
        self.colors = ["black","blue","red","green","brown","orange","pink","magenta"]
        self.size = 550
        self.indend = 35
        self.r = (self.size - 2.0*self.indend)/2.0
        self.center = (self.size-2.0*self.indend)/2.0+self.indend
        self.svg_holder = ""
        self.task_lines = {"graphs":{},"circles":""}
        self.mge_boxes = []
        
    def clear_tasks(self):
        self.task_lines = {"graphs":{},"circles":""}
        self.added_tasks = 0
       
    def add_task(self,task,task_description,windows):
        mean,stdev = task_description["statistics"]
        if "GC" in self.task_list or "AT" in self.task_list:
            band = self.r/2.0/len(self.task_list)
        else:
            band = self.r/2.0/(len(self.task_list)+1)
        mid = self.r - self.r/6.0 - band/2.0 - self.added_tasks*band
        if not task_description["condition"]:
            self.task_lines["circles"] = ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"white\" stroke=\"grey\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid))
        elif task_description["condition"] == "bigger than":
            self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid+band/2.0,"aliceblue"))
            self.task_lines["circles"] += self.task_inner_circle(task_description["mode"],"white",band,mid,mean,stdev,float(task_description["val1"]))
        elif task_description["condition"] == "smaller than":
            self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid+band/2.0,"white"))
            self.task_lines["circles"] += self.task_inner_circle(task_description["mode"],"aliceblue",band,mid,mean,stdev,float(task_description["val1"]))
        elif task_description["condition"] == "between":
            self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,mid+band/2.0,"white"))
            self.task_lines["circles"] += self.task_inner_circle(task_description["mode"],"aliceblue",band,mid,mean,stdev,float(task_description["val1"]),float(task_description["val2"]))
        self.task_lines["circles"] += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"none\" stroke=\"black\" stroke-width=\"0.5\" />\n" %
            (self.center,self.center,mid-band/2.0))
        start = 0
        VAL = None
        wins = list(windows.keys())
        try:
            wins.sort(key=lambda s: int(s.split("-")[0]),reverse=True)
        except:
            wins.sort(key=lambda s: int(s.split("..")[0]),reverse=True)
        self.task_lines['graphs'][task] = ""
        for win in wins:
            try:
                first,second = map(lambda v: float(v),win.split("-"))
            except:
                first,second = map(lambda v: float(v),win.split(".."))
            stop = (first+second)/2.0
            deviation = band*(windows[win]-mean)/6.0/stdev
            if VAL == None:
                VAL = deviation
                start = stop
                a = math.pi/2.0 - 2.0*math.pi*start/self.seqlength
                x = self.center + (mid+VAL)*math.cos(a)
                y = self.center - (mid+VAL)*math.sin(a)
                self.task_lines['graphs'][task] += ("<path style=\"stroke:%s;stroke-width:1.0\" d=\"M%f,%f" % (self.colors[self.added_tasks],x,y))
                continue
            b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
            x = self.center + (mid+deviation)*math.cos(b)
            y = self.center - (mid+deviation)*math.sin(b)
            self.task_lines['graphs'][task] += ("L%f,%f" % (x,y))
            VAL = deviation
            start = stop
        self.task_lines['graphs'][task] += "Z\" />\n"
        self.task_lines['graphs'][task] += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"%s\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (50,self.size+25*self.added_tasks,100,self.size+25*self.added_tasks,self.colors[self.added_tasks]))
        self.task_lines['graphs'][task] += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start;fill:black\">%s</text>\n" % 
            (125,self.size+3+25*self.added_tasks,task))
        self.added_tasks += 1
    
    def task_inner_circle(self,mode,fill_color,band,mid,mean,stdev,value1,value2=None):
        if mode=="absolute":
            r = mid+band*(value1-mean)/6.0/stdev
            if r > mid+band/2.0-1:
                r = mid+band/2.0-1
            line = ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"none\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,r,fill_color))
        elif mode=="sigmas":
            r = mid+band*value1/6.0
            if r > mid+band/2.0-1:
                r = mid+band/2.0-1
            line = ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" stroke=\"none\" stroke-width=\"0.5\" />\n" %
                (self.center,self.center,r,fill_color))
        elif mode=="fraction":
            pass
        if value2 != None:
            line += "\n" + self.task_inner_circle(self,mode,"white",band,mid,mean,stdev,value2)
        return line
        
    def add_gi(self,lb,rb,color="",title=""):
        lb = int(lb)
        rb = int(rb)
        r1 = self.r-10
        r2 = self.r-self.r/8.0
        a = math.pi/2.0 - 2.0*math.pi*lb/self.seqlength
        x = self.center + r1*math.cos(a)
        y = self.center - r1*math.sin(a)
        if color:
            if color == "grey":
                self.flg_wrongLoci = True
            color = "fill=\"%s\"" % color
        if title:
            self.mge_boxes.append("<a href=\"\" xlink:title=\"%s\">" % title)
        self.mge_boxes.append("<path %s d=\"M%f,%f" % (color,x,y))
        step = self.seqlength/360.0
        stop = lb + step
        b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
        while stop <= rb:
            b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
            x = self.center + r1*math.cos(b)
            y = self.center - r1*math.sin(b)
            self.mge_boxes[-1] += ("L%f,%f" % (x,y))
            stop += step
        x = self.center + r2*math.cos(b)
        y = self.center - r2*math.sin(b)
        self.mge_boxes[-1] += ("L%f,%f" % (x,y))
        stop -= step
        while stop >= lb:
            b = math.pi/2.0 - 2.0*math.pi*stop/self.seqlength
            x = self.center + r2*math.cos(b)
            y = self.center - r2*math.sin(b)
            self.mge_boxes[-1] += ("L%f,%f" % (x,y))
            stop -= step
        self.mge_boxes[-1] += "Z\" />\n"
        if title:
            self.mge_boxes.append("</a>")
        
    def get_svg(self):
        self.svg_holder = "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewbox=\"0 0 %d %d\">\n" % (self.size,self.size+25*(len(self.task_list)+1))
        self.svg_holder += ("<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"none\" stroke=\"black\" stroke-width=\"1.0\" />\n" %
            ((self.size-2*self.indend)/2+self.indend,(self.size-2*self.indend)/2+self.indend,
            (self.size-2*self.indend)/2))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            ((self.size-2*self.indend)/2+self.indend,self.indend-20,
            (self.size-2*self.indend)/2+self.indend,self.indend+20))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (self.size-self.indend-20,(self.size-2*self.indend)/2+self.indend,
            self.size-self.indend+20,(self.size-2*self.indend)/2+self.indend))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            ((self.size-2*self.indend)/2+self.indend,self.size-self.indend-20,
            (self.size-2*self.indend)/2+self.indend,self.size-self.indend+20))
        self.svg_holder += ("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" fill=\"none\" stroke=\"grey\" stroke-width=\"1.0\" stroke-linejoin=\"round\" stroke-miterlimit=\"10\" />\n" %
            (self.indend-20,(self.size-2*self.indend)/2+self.indend,
            self.indend+20,(self.size-2*self.indend)/2+self.indend))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:middle\" font-size=\"15\">%s</text>" %
            (self.size/2,10,self.seqname))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%s</text>" %
            ((self.size-2*self.indend)/2+self.indend+10,self.indend-10,self.format_numeric_string(self.seqlength)))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%G Mbp</text>" %
            (self.size-self.indend+10,(self.size-2*self.indend)/2+self.indend-10,self.format_number(self.seqlength/4.0,2,-6)))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:start\">%G Mbp</text>" %
            ((self.size-2*self.indend)/2+self.indend+10,self.size-self.indend+15,self.format_number(self.seqlength/2.0,2,-6)))
        self.svg_holder += ("<text x=\"%d\" y=\"%d\" style=\"text-anchor:end\" transform=\"rotate(-90)\">%G Mbp</text>" %
            (-(self.size-2*self.indend)/2-self.indend+75,self.indend-10,self.format_number(3.0*self.seqlength/4.0,2,-6)))
        dot = 100000
        while dot < self.seqlength:
            a = math.pi/2.0 - 2.0*math.pi*dot/self.seqlength
            if dot%1000000 == 0:
                s = 3.0
                fill_color = "black"
            else:
                s = 2.5
                fill_color = "grey"
            self.svg_holder += ("<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"%s\" stroke=\"black\" stroke-width=\"1.0\" />\n" %
                (self.center+self.r*math.cos(a),self.center-self.r*math.sin(a),s,fill_color))
            dot += 100000
        if self.task_lines["circles"]:
            self.svg_holder += "\n"+self.task_lines["circles"]+"\n"
        if self.task_lines['graphs']:
            self.svg_holder += "<g style=\"fill:none;stroke-linejoin:round;stroke-miterlimit:10\">\n"
            for task in self.task_lines['graphs']:
                self.svg_holder += self.task_lines['graphs'][task]
            self.svg_holder += "</g>\n"
        if self.mge_boxes:
            self.svg_holder += "<g style=\"fill:pink;stroke:black\">\n"
            for box in self.mge_boxes:
                self.svg_holder += box
            self.svg_holder += "</g>\n"
        if len(self.title) <= 11:
            self.svg_holder += "<text x=\"275\" y=\"280\" style=\"text-anchor:middle\" font-size=\"25\" font-weight=\"bold\">%s</text>\n" % self.title
        self.svg_holder += "<path  style=\"fill:pink;stroke:black\" d=\"M350,550L375,550L375,565L350,565Z\" />\n"
        self.svg_holder += "<text x=\"405\" y=\"565\" style=\"text-anchor:start;fill:black\">Selected loci;</text>\n"
        if self.flg_wrongLoci:
            self.svg_holder += "<path  style=\"fill:grey;stroke:black\" d=\"M350,600L375,600L375,615L350,615Z\" />\n"
            self.svg_holder += "<text x=\"405\" y=\"615\" style=\"text-anchor:start;fill:black\">Falsely selected rrn operons;</text>\n"
        return self.svg_holder + "</svg>"

    def format_number(self,num,dig,zoom=0):
        return int((10**(dig+zoom))*num)/float(10**dig)
    
    def format_numeric_string(self,num):
        outnum = str(num)
        if len(outnum) <= 3:
            return outnum
        for i in range(len(outnum)-3,0,-3):
            outnum = outnum[:i]+","+outnum[i:]
        return outnum

    def sort_by_coordinates(self, a, b):
        try:
            valA = int(a.split("-")[0])
            valB = int(b.split("-")[0])
        except:
            valA = int(a.split("..")[0])
            valB = int(b.split("..")[0])
        if valA < valB: return -1
        elif valA > valB: return 1
        else: return 0

##############################################################################################################
class Validator:
    def __init__(self):
        pass
        # ATTRIBUTES
        self.MinimalLength = {
            7:295000,
            6:74000,
            5:18500,
            4:4600,
            3:1200,
            2:300
            }
        self.tasks = {'GC':'GC-content',
                      'GCS':'G/C-skew',
                      'ATS':'A/T-skew',
                      'D':'Distance',
                      'GD':'Generalized distance',
                      'PS':'Pattern skew',
                      'RPS':'Relative pattern skew',
                      'GPS':'Generalized pattern skew',
                      'GRPS':'Generalized relative pattern skew',
                      'V':'Variance',
                      'GV':'Generalized variance',
                      'RV':'Relative variance',
                      'GRV':'Generalized relative variance'}

    # METHODS
    def validateTasks(self,task_list,frame):
        tasks = {}
        # check tasks
        LengthLimits = []
        for currTask in task_list:
            SuperTask = currTask.split('-')
            subtrahend = None
            divisor = None
            if len(SuperTask) > 2:
                showError("Wrong task: " + currTask)
                return None
            elif len(SuperTask) == 2:
                subtrahend = SuperTask[1]
                divisor = None
            else:
                SuperTask = SuperTask[0].split('/')
                if len(SuperTask) > 2:
                    showError("Wrong task: " + currTask)
                    return None
                elif len(SuperTask) == 2:
                    divisor = SuperTask[1]
                    subtrahend = None
                else:
                    pass
            TaskId = SuperTask[0]
            
            Values = TaskId.split(":")
            if len(Values) == 1:
                Values = [None,Values[0]]
            elif len(Values) == 2:
                pass
            else:
                showError("Wrong task: '" + str(TaskId) + "'.\nMust be like 'n0_4mer:D'")
                return None
            
            Task = Values[1]
            if Values[0]:
                Type = str(Values[0]).split("_")
                if len(Type) != 2:
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                
                Norm = str(Type[0])
                if len(Norm) != 2 or Norm[0] != "n":
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                try:
                    pattern_type = Norm[0]
                    Norm = int(Norm[1])
                except:
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                
                Wlength = str(Type[1])
                if len(Wlength) != 4 or Wlength[1:] != "mer":
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
                try:
                    Wlength = int(Wlength[0])
                except:
                    showError("Wrong type: " + str(Values[0]) + "\nMust be like 'n0_4mer'")
                    return None
            else:
                Task = Values[1]
                Norm = 0
                Wlength = 4
                pattern_type = "n"
            
            # add task to the task list
            tasks[currTask] = {'task':Task,
                            'norm':Norm,
                            'type':pattern_type,
                            'wlength':Wlength,
                            'subtr':subtrahend,
                            'divisor':divisor,
                            'ID':TaskId}
            subtrahend = None
            divisor = None
            LengthLimits.append(self.get_MinimalLength(Wlength))

            # check frame sizes
            if frame < min(LengthLimits):
                result = showError("The sliding window size must be at least " + str(min(LengthLimits)) + " bp.")
                return None
            else:
                pass
        return tasks

    # ACCESSIONS
    def get_MinimalLength(self, n):
        return self.MinimalLength[n]

    def get_TaskCategory(self, task):
        return self.tasks[task]

###################################################################################################################
class Pattern:
    def __init__(self,wlength):
        # ATTRIBUTES
            # Number of words in collection
        self.TotalWordNumber = 0
            # Word length
        self.wlength = wlength
            # Collection of words: word name like {'AAAA' : [real number, expected number, word spreading]}
            # diviation is measured in values of expected standard diviation calculated by regression
            # equation for given word length and sequence length
        self.words = {}
            # AGCT-box
        self.boxAGCT = {'A':0,'G':0,'C':0,'T':0}
            # Normalization type
        self.normalization = None
            # object NormalizationTable
        self.oNormalizationTable = None
            # reverse complementation
        self.RevComplTable = {'A':'T','T':'A','G':'C','C':'G'}
            # pattern type
        self.PatternType = None
            # word statistics
        self.flg_addWordStatistics = 0
        
        self.initiate()

    # METHODS
    # initiate words
    def initiate(self):
        # Fill in the generic pattern with words
        pmlib = PatternMethod()
        WordList = pmlib.getWordList(self.wlength)
        for i in range(len(WordList)):
            # word data [real,expected,regularity]
            self.words[WordList[i]] = [0,0,None]
            

    # Methods
    def setPattern(self, strSeq, normalization, pattern_type):
        if not len(strSeq):
            print(861)
        self.PatternType = pattern_type
        self.normalization = normalization
        if self.PatternType in ("d","s"):
            flg_addWordStatistics = 1
        else:
            flg_addWordStatistics = 0
        if normalization == 0:
            normalization = 1
        
        # set mononucleotite matrix
        strSeq = strSeq.upper()
        self.setMatrix(strSeq)
        # Calculate expected and real numbers of words: [real number, expected number]
        self.setWordStatistics(strSeq, flg_addWordStatistics)
        if self.PatternType in ("d","s"):
            self.processStatData()
        # set normalization table
        if self.oNormalizationTable == None:
            WordFrequencies = {}
            if normalization == 1:
                WordFrequencies.update(self.boxAGCT)
            else:
                oNormPattern = Pattern(normalization)
                oNormPattern.setPattern(strSeq)
                for word in oNormPattern.words.keys():
                    WordFrequencies[word] = oNormPattern.getWordFrequency(word)
                del oNormPattern
            oNormalizationTable = NormalizationTable(WordFrequencies,normalization)
            self.setNormalizationTable(oNormalizationTable)
            self.setExpectation()
        else:
            self.oNormalizationTable = None

    def setExpectation(self):
        wordlist = self.words.keys()
        pVal = float(len(wordlist)/self.TotalWordNumber)
        for word in wordlist:
            if self.oNormalizationTable:
                pVal = self.oNormalizationTable.getWordLikelihood(word)
            self.words[word][1] = pVal*self.TotalWordNumber
    
    def setMatrix(self, strSeq):
        # Calculate numbers of A, G, C and T in the sequence and assign values to the 4-cell table boxAGCT
        letters = self.boxAGCT.keys()
        for l in letters:
            self.boxAGCT[l] = float(strSeq.count(l))/len(strSeq)
        return 1

    def setWordStatistics(self, strSeq, flg_addWordStatistics, option=None):
        keys = self.words.keys()
        # set real numbers of words and total number of words depending on allowing of word overlapping
        badWords = ''
        # total number of words
        self.TotalWordNumber = len(strSeq)-self.wlength+1
        # real number of words
        WordPosition = {}
        # number of lots for intrinsic PS
        lNum = self.getLotLength(self.TotalWordNumber)
        lot = self.TotalWordNumber/lNum
        k = 1
        for i in range(self.TotalWordNumber):
            if k and i >= k*lot:
                k = k + 1
                if k > lNum:
                    k = None
            word = strSeq[i:i+self.wlength]
            try:
                self.words[word][0] += 1
                if flg_addWordStatistics:
                    # set statistics for intrinsic pattern skew
                    if self.words[word][2]:
                        distance = i - WordPosition[word]
                        WordPosition[word] = i
                        self.words[word][2]['Sum 1'] = self.words[word][2]['Sum 1'] + distance
                        self.words[word][2]['Sum 2'] = self.words[word][2]['Sum 2'] + distance*distance
                        if k:
                            self.words[word][2]['PartialFrequ'][k-1] = self.words[word][2]['PartialFrequ'][k-1] + 1
                    else:
                        WordPosition[word] = i
                        self.words[word][2] = {'Sum 1':0,
                                               'Sum 2':0,
                                               'Start':i,
                                               'PartialFrequ':[]}
                        for j in range(lNum):
                            self.words[word][2]['PartialFrequ'].append(0)
            except:
                if word.find('*') == -1:
                    badWords = badWords + word + '; '
        return badWords

    def processStatData(self):
        tmp_words = []
        for word in self.words:
            try:
                RevComplWord = self.getReverseComplement(word)
            except:
                tmp_words.append([word,None])
            try:
                lNum = len(self.words[word][2]['PartialFrequ'])
            except:
                tmp_words.append([word,None])
            if word == RevComplWord:
                tmp_words.append([word,0])

            words = self.words.keys()
            sum = 0
            
            # distance in ranks
            for i in range(lNum):
                data = []
                for w in words:
                    if self.words[w][2]:
                        data.append([self.words[w][2]['PartialFrequ'][i],w])
                    else:
                        data.append([0,w])
                #data.sort(self.sort_by_ranks)
                data.sort(key=lambda ls: [ls[0],ls[1]],reverse=True)
                r = 1
                for j in range(len(data)):
                    if data[j][1] == word or data[j][1] == RevComplWord:
                        if r:
                            rank = j
                            r = r - 1
                        else:
                            sum = sum + abs(rank - j)
                            break
            tmp_words.append([word,float(sum)/float(lNum)])

        for item in tmp_words:
            if self.words[item[0]][2]:
                self.words[item[0]][2]['PartialFrequ'] = item[1]
            else:
                self.words[item[0]][2] = {'Sum 1':0,
                                               'Sum 2':0,
                                               'Start':0,
                                               'PartialFrequ':0}

    def getLotLength(self, slength):
        MinimalSequencesLength = {
            7:295000,
            6:74000,
            5:18500,
            4:4600,
            3:1200,
            2:300
            }
        minLength = MinimalSequencesLength[self.wlength]
        lLen = int(2.0*math.sqrt(float(slength+self.wlength)/float(minLength)))
        if lLen:
            return lLen
        else:
            return 1

    def getReverseComplement(self, word):
        RevComplWord = ''
        for l in word:
            RevComplWord = self.RevComplTable[l] + RevComplWord
        return RevComplWord
            
    def getWordLength(self):
        return self.wlength
    
    def getPatternType(self):
        return self.PatternType
    
    def getPatternName(self):
        return self.getPatternType() + str(self.getNormalization()) + "_" + str(self.getWordLength()) + "mer"
    
    def getNormalization(self):
        return self.normalization

    def getSeqLength(self):
        return self.TotalWordNumber + self.wlength
    
    def getWordRanks(self):
        oWordList = self.getWordList()
        return oWordList.getWordList()

    # take pattern, standard list of words and mode, return ['word','pattern type'],
    # if take None, return current WordList
    def getWordList(self):
        oWordList = WordList(self.getWordLength())
        #normstdiv = 0.14+(4**self.getWordLength())/self.getSeqLength()
        normal_expect = self.getMeanWordNumber()
        for item in oWordList.wordList:
            value = self.getWordAbundance(item[0])
            # set word frequencies
            item[3] = self.getWordFrequency(item[0])
            if self.PatternType == 'n' and self.normalization == 0:
                item[2] = self.getDeviation(value,normal_expect,normal_expect)
            elif self.PatternType == 's':
                # item[2] = self.getWordSkewness(item[0])/stdivIntrinsicPS
                item[2] = self.getWordSkewness(item[0])/10.0
                if item[2] == None:
                    return None
            elif self.PatternType == 'd':
                IPS = self.getWordSkewness(item[0])/10.0
                if IPS == None:
                    return None
                if self.normalization:
                    # normalized patterns
                    # expected number of word
                    expect = self.getWordExpectation(item[0])
                    # correction for the letter inequivalence
                    if expect == 0:
                        expect = 1
                    if normstdiv == 0:
                        normstdiv = 1
                    item[2] = [self.getDeviation(value,expect,normal_expect),IPS]
                else:
                    item[2] = [self.getDeviation(value,normal_expect,normal_expect),IPS]
            else:
                # normalized patterns
                # expected number of word
                expect = self.getWordExpectation(item[0])
                # correction for the letter inequivalence
                if expect == 0:
                    expect = 1
                '''
                if normstdiv == 0:
                    normstdiv = 1
                '''
                item[2] = self.getDeviation(value,expect,normal_expect)
        return oWordList
    
    def getDeviation(self,value,expect,stdev,flg_lognormal=0):
        if flg_lognormal:
            if value == 0:
                value = 1
            numerator = math.log((value**2)*math.sqrt(stdev**2 + expect**2)/(expect**2)/math.sqrt(stdev**2 + value**2))
            denumerator = math.sqrt(math.log((stdev/expect)**2 + 1))
            return 6.0*numerator/denumerator
        else:
            normstdiv = 0.14+(4**self.getWordLength())/self.getSeqLength()
            return float((value - expect)/stdev)/normstdiv

    def getNormalizationTable(self):
        return self.oNormalizationTable
    
    def getNormalizationTableParametr(self):
        if self.oNormalizationTable:
            return self.oNormalizationTable.getWordLength()
        else:
            return None

    def setNormalizationTable(self, oNormalizationTable):
        del self.oNormalizationTable
        self.oNormalizationTable = oNormalizationTable
        self.setExpectation()

    def getWordAbundance(self, word):
        return self.words[word][0]

    def getWordFrequency(self, word):
        return float(self.words[word][0])/float(self.TotalWordNumber)

    def getWordExpectation(self, word):
        if self.normalization == 0:
            return self.getMeanWordNumber()
        else:
            return self.words[word][1]

    def getWordSkewness(self, word):
        try:
            return float(self.words[word][2]['PartialFrequ'])
        except:
            return None

    def getMeanWordNumber(self):
        return self.TotalWordNumber/len(self.words.keys())
    
    def getOUV(self):
        oWordList1 = self.getWordList()
        wordList = oWordList1.getWordList(self.PatternType)[0]
        sum = 0
        sum_sq = 0
        N = len(wordList)
        for i in range(N):
            sum = sum + wordList[i][2]
            sum_sq = sum_sq + wordList[i][2] * wordList[i][2]
        val = (sum_sq - sum * sum/N)/(N - 1)
        #return math.sqrt(val)
        return val
    
    def getPS(self):
        # set variables
        oWordList1 = self.getWordList()
        wl1 = oWordList1.getWordList(self.PatternType)
        wl2 = self.getComplement(wl1)
        if self.PatternType == "d":
            wl1 = [wl1[0],]
            wl2 = [wl2[0],]
        return self.getDistance(wl1,wl2,self.getMinimalDistance(),0,0)[0]
    
    def getPercentage(self,symbol):
        if len(symbol)==1:
            val = self.getNucleotideFrequency(symbol)
        elif symbol == "GC":
            val = self.getNucleotideFrequency("G") + self.getNucleotideFrequency("C")
        elif symbol == "AT":
            val = self.getNucleotideFrequency("A") + self.getNucleotideFrequency("T")
        else:
            val = 0
        return val*100
    
    def getGCskew(self):
        G = self.getNucleotideFrequency("G")
        C = self.getNucleotideFrequency("C")
        return (float(G-C)/float(G+C))
        
    def getATskew(self):
        A = self.getNucleotideFrequency("A")
        T = self.getNucleotideFrequency("T")
        return (float(A-T)/float(A+T))

    def getMinimalDistance(self):
        if self.wlength%2 == 0:
            return (4**self.wlength - 2**self.wlength)/2.0
        else:
            return (4**self.wlength)/2.0

    def getNucleotideFrequency(self,letter):
        return self.boxAGCT[letter]
    
    def getComplement(self, wl):
        NewWordList = []
        for d in range(len(wl)):
            WordList = []
            for item in wl[d]:
                WordList.append([])
                for val in item:
                    WordList[len(WordList)-1].append(val)
            WordList = sorted(WordList,key=lambda ls: [ls[2],ls[0]],reverse=True)
            WordDic = {}
            WordNumber = len(WordList)
            c = 1
            for i in range(WordNumber):
                WordDic[WordList[i][0]] = [WordList[i][0],WordList[i][1],WordList[i][2],WordList[i][3]]
            for i in range(WordNumber):
                arWord = array.array('u',list(WordList[WordNumber-i-1][0]))
                #arWord.fromstring(WordList[WordNumber-i-1][0])
                arWord.reverse()
                try:
                    WordDic[arWord.tostring()][2] = WordList[i][2]
                    WordDic[arWord.tostring()][3] = None
                    del arWord
                except:
                    c = c+1
            NewWordList.append(WordDic.values())

        del WordList
        del WordDic
        # sort values and assign ranks
        for i in range(len(NewWordList)):
            NewWordList[i] = sorted(NewWordList[i],key=lambda ls: [ls[2],ls[0]],reverse=True)
            for j in range(len(NewWordList[i])):
                NewWordList[i][j][3] = j
                
        return NewWordList

    # return a copy of the current pattern object
    def copy(self):
        oPattern = Pattern(self.wlength)
        oPattern.TotalWordNumber = self.TotalWordNumber
        oPattern.PatternType = self.PatternType
        oPattern.words = self.words
        oPattern.boxAGCT = self.boxAGCT
        oPattern.normalization = self.normalization
        oPattern.flg_addWordStatistics = self.flg_addWordStatistics
        try:
            oPattern.setNormalizationTable(self.getNormalizationTable())
        except:
            oPattern.setNormalizationTable(None)
        return oPattern
    
    def convert(self, pattern_type=None, normalization=None):
        oPattern = self.copy()
        if pattern_type == None:
            pattern_type = self.getPatternType()
        if normalization == None:
            normalization = self.getNormalization()
        if pattern_type == self.getPatternType():
            if normalization == self.getNormalization():
                return oPattern
            if normalization == 0 or normalization == self.getNormalizationTableParametr():
                oPattern.normalization = normalization
            else:
                return None
            oPattern.setExpectation()
            return oPattern
        elif pattern_type == "n" or oPattern.flg_addWordStatistics:
            oPattern.PatternType = pattern_type
        else:
            return None
        
    def compare(self,other):
        if type(other) != type(Pattern(2)):
            return None
        if (self.normalization != other.normalization) or (self.PatternType != other.PatternType):
            return None
        else:
            return 1
    
    # calculate distance between two pattern lists 
    def getDistance(self,first_list,second_list,minDist=0,flg_BestHit=0,flg_Normolize=0):
        wl1 = first_list
        dm = min((len(first_list),len(second_list)))
        total_sum = 0
        wlength = len(first_list[0][0][0])
        if dm > 1:
            WordNumber = 4**wlength/2
        else:
            WordNumber = 4**wlength
        maxDist = float(WordNumber*(WordNumber+1))
        StdDev = 0
        dist = maxDist
        StdDev = 0
        complement = 0
        if flg_BestHit:
            count = 3
        else:
            count = 0
        # calculate distance
        while count >= 0:
            sum = []
            values = []
            for i in range(dm):
                if count==1:
                    wl1 = first_list
                    #wl1[i].sort(self.sort_by_order)
                    wl1[i].sort(key = lambda ls: ls[1])
                    wl2 = self.getComplement(second_list)
                    #wl2[i].sort(self.sort_by_order)
                    wl2[i].sort(key = lambda ls: ls[1])
                elif count==2:
                    wl1 = self.getComplement(first_list)
                    #wl1[i].sort(self.sort_by_order)
                    wl1[i].sort(key = lambda ls: ls[1])
                    wl2 = second_list
                    #wl2[i].sort(self.sort_by_order)
                    wl2[i].sort(key = lambda ls: ls[1])
                elif count==3:
                    wl1 = self.getComplement(first_list)
                    #wl1[i].sort(self.sort_by_order)
                    wl1[i].sort(key = lambda ls: ls[1])
                    wl2 = self.getComplement(second_list)
                    #wl2[i].sort(self.sort_by_order)
                    wl2[i].sort(key = lambda ls: ls[1])
                else:
                    wl1 = first_list
                    #wl1[i].sort(self.sort_by_order)
                    wl1[i].sort(key = lambda ls: ls[1])
                    wl2 = second_list
                    #wl2[i].sort(self.sort_by_order)
                    wl2[i].sort(key = lambda ls: ls[1])
                sum.append(0)
                values.append([])
                for j in range(WordNumber):
                    val = abs(wl1[i][j][3]-wl2[i][j][3])/2
                    sum[i] += val
                    values[i].append(val)
            subtotal = 0
            for val in sum:
                subtotal = subtotal + val*val
            if subtotal:
                subtotal = math.sqrt(subtotal)
            if  subtotal < dist:
                dist = subtotal
                StdDev = self.getMeanAndStDev(values)[1]
            if count == 0:
                break
            else:
                count -= 1
        # normolize distance
        if flg_Normolize:
            OUV1 = getOUV(wl1[0])
            OUV2 = getOUV(wl2[0])
            m = -.25*math.sqrt(OUV1*OUV1+OUV2*OUV2)+4
        else:
            m = 1.0
        if not dist:
            dist = 0.0000001
        dist = 400.0*m*(dist-minDist)/(math.sqrt(dm)*maxDist-minDist)/dm
        return [dist,StdDev,complement]

    # Get list of values and return [mean value,standard deviation]
    def getMeanAndStDev(self,DataList,ind=2,binome=None):
        sum = 0.0
        sum_sq = 0.0
        sum_cube = 0.0
        # if there is a number of lists, recalculate stdDev for all lists and return the average and the skew for the last list
        if type(DataList[0]) == type([1,2]):
            dm = len(DataList)
            N = len(DataList[0])
        else:
            dm = 1
            N = len(DataList)
            DataList = [DataList,]
        StdDev = 0
        for k in range(dm):
            for i in range(N):
                val = float(DataList[k][i])
                sum = sum + val
                if ind > 1:
                    sum_sq = sum_sq + val*val
                if ind > 2:
                    sum_cube = sum_cube + val*val*val
            N = float(N)
            mean = sum/N
            if ind > 1:
                if binome:
                    if mean > 1:
                        bmean = mean/100.0
                        stdev = math.sqrt(bmean*(1-bmean))*100.0
                    else:
                        stdev = math.sqrt(mean*(1-mean))
                else:
                    var = sum_sq/N - mean*mean
                    if var>0:
                        stdev = math.sqrt(var)
                    else:
                        stdev = 0
            else:
                stdev = None
            if ind > 2:
                if binome:
                    skew = 0
                else:
                    if stdev:
                        skew = (sum_cube - 3.0*sum*sum_sq/N + 2.0*sum*sum*sum/N/N)/stdev/stdev/stdev/N
                    else:
                        skew = 0
            else:
                skew = None
            if stdev:
                StdDev = math.sqrt(stdev*stdev + StdDev*StdDev)
            else:
                StdDev = 0
        return [mean,StdDev,skew]

    # Subtracting of patterns
    def __sub__(self,other):
        # set variables
        if self.getPatternName() != other.getPatternName():
            return None
        oWordList1 = self.getWordList()
        oWordList2 = other.getWordList()
        if not oWordList1 or not oWordList2:
            return None
        wl1 = oWordList1.getWordList(self.PatternType)
        wl2 = oWordList2.getWordList(other.PatternType)
        return self.getDistance(wl1,wl2,0,1,0)[0]

    # Summation of patterns
    def __add__(self,other):
        if (self.wlength != other.getWordLength()) or (self.PatternType != other.PatternType):
            return None
        oPattern = Pattern(self.wlength)
        oPattern.PatternType = self.PatternType
        for key in self.boxAGCT.keys():
            oPattern.boxAGCT[key] = (self.boxAGCT[key]*self.getSeqLength() + other.boxAGCT[key]*other.getSeqLength())/(self.getSeqLength() + other.getSeqLength())
        oPattern.TotalWordNumber = self.TotalWordNumber + other.TotalWordNumber
        first_table = self.getNormalizationTable()
        second_table = other.getNormalizationTable()
        if first_table:
            first_words = first_table.getWords()
            second_words = second_table.getWords()
            for key in first_words:
                first_words[key] = (first_words[key]*self.getSeqLength() + second_words[key]*other.getSeqLength())/(self.getSeqLength() + other.getSeqLength())
            new_table = NormalizationTable(first_words)
            oPattern.setNormalizationTable(new_table)
        for key in self.words.keys():
            oPattern.words[key][0] = self.words[key][0] + other.words[key][0]
        oPattern.setExpectation()
        return oPattern

    '''
    def sort_by_values(self, a, b):
        if b[2] < a[2]:
            return -1
        elif b[2] > a[2]:
            return 1
        else:
            if b[0] < a[0]:
                return -1
            elif b[0] > a[0]:
                return 1
            else:
                return 0

    def sort_by_order(self, a, b):
        if b[1] < a[1]:
            return -1
        elif b[1] > a[1]:
            return 1
        return 0

    def sort_by_ranks(self, a, b):
        if a[0] > b[0]:
            return 1
        elif a[0] < b[0]:
            return -1
        else:
            if a[1] > b[1]:
                return 1
            elif a[1] < b[1]:
                return -1
            else:
                return 0
    '''
###################################################################################################################
class WordList:
    def __init__(self,wlength):

        # ATTRIBUTES
        self.wlength = wlength
        self.wordList = []
        # Create library of pattern setting methods
        pmlib = PatternMethod()
        # Create list and dictionary of the words
        for i in range(4**self.wlength):
            word = pmlib.getNextWord(self.wlength)
            CurrWord = ['',0,0,None]
            CurrWord[0] = word
            CurrWord[1] = i
            self.wordList.append(CurrWord)

    # METHODS

    def getRange(self):
        if type(self.wordList[0][2]) == type([0,]):
            return len(self.wordList[0][2])
        else:
            return 1

    # n is int - index of values in the list
    def getWordList(self,pattern_type="n"):
        wl = []
        wl1 = []
        wl2 = []
        for item in self.wordList:
            if pattern_type == "d":
                wl1.append([item[0],item[1],item[2][0],0])
                wl2.append([item[0],item[1],item[2][1],0])
            else:
                wl1.append([item[0],item[1],item[2],0])
        wl.append(wl1)
        if len(wl2):
            wl.append(wl2)
        for l in wl:
            self.assignRanks(l)                
        return wl

    def assignRanks(self, WordList):
        WordList = sorted(WordList,key=lambda ls: [ls[2],ls[0]],reverse=True)
        for j in range(len(WordList)):
            WordList[j][3] = j
        WordList = sorted(WordList,key=lambda ls: ls[1]) 

    def getLength(self):
        return len(self.wordList)

    def getVariance(self):
        sum = 0
        sum_sq = 0
        N = len(self.wordList)
        for i in range(N):
            sum = sum + self.wordList[i][2]
            sum_sq = sum_sq + self.wordList[i][2] * self.wordList[i][2]
        val = (sum_sq - sum * sum/N)/(N - 1)
        return val

##############################################################################################################
class PatternMethod:
    def __init__(self, method="BOXAGCT"):
            # Method
        self.CurrMethod = method
            # Current word
        self.Word = ""
            # The letter of the word in focus
        self.Cursor = None
            # State of the word processing by any function
        self.flg_State = 0
            # Dictionary of methods of oligomer plot arrangment
        self.Methods = {}
            # Dictionary of letter functions 
        self.LetterFunk = {}
        
        self.setDictionary()

    # Methods
        # Block of functions for method BOXAGCT
    def A_BOXAGCT(self):
        self.changeLetters("G","C",self.Cursor)
        self.flg_State = 0       

    def G_BOXAGCT(self):
        self.changeLetters("A","C",self.Cursor)
        if self.Cursor >= 0:
            self.Cursor = self.Cursor - 1
        else:
            self.flg_State = 0
            
    def C_BOXAGCT(self):
        self.changeLetters("T","A",self.Cursor)
        if self.Cursor >= 0:
            self.flg_State = 0
        else:
            self.Cursor = self.Cursor - 1

    def T_BOXAGCT(self):
        self.changeLetters("C","A",self.Cursor)
        self.Cursor = self.Cursor - 1
 
    def changeLetters(self, a, b, position):
        if position >= 0:
            self.Word = self.Word[:position] + a + self.Word[(position+1):]
        else:
            if position == -1:
                self.Word = self.Word[:position] + b + self.Word[:(position+1)]
            else:
                self.Word = self.Word[:position] + b + self.Word[(position+1):]

    def setDictionary(self):
        # initialisation of dictionary of functions for method BOXAGCT
        self.LetterFunk["A"] = self.A_BOXAGCT
        self.LetterFunk["G"] = self.G_BOXAGCT
        self.LetterFunk["C"] = self.C_BOXAGCT
        self.LetterFunk["T"] = self.T_BOXAGCT
        self.Methods["BOXAGCT"] = self.LetterFunk

    def getNextWord(self, wlength):
        if self.Word == '':
            self.Word = "A"*wlength
        else:
            # create new word
            self.Cursor = wlength - 1
            self.flg_State = 1
            while self.flg_State == 1:
                self.Methods[self.CurrMethod][self.Word[self.Cursor]]()
        return self.Word

    def getWordList(self, wlength):
        if wlength != 0:
            self.Word = ''
            wordList = []
            for i in range(4**wlength):
                wordList.append(self.getNextWord(wlength))
            return wordList

###################################################################################################################
# words is a dictionary {'word':frequency}
class NormalizationTable:
    def __init__(self, words, wlength):
        
        # ATTRIBUTES
        self.words = words
        self.wlength = wlength

    # METHODS

    def getWordLikelihood(self, word):
        frame = len(list(self.words.keys())[0])
        pVal = self.words[word[:frame]]
        subword = word[1:frame]
        subset = self.setSubset(subword)
        start = 1
        stop = frame+1
        while stop <= len(word):
            pVal = pVal * subset[word[start:stop]]
            start = start + 1
            subword = word[start:stop]
            subset = self.setSubset(subword)
            stop = stop + 1
        return pVal

    def setSubset(self, subword):
        subset = {}
        wordlist = self.words.keys()
        TotalLikelihood = 0
        for word in wordlist:
            if subword == word[:len(subword)]:
                subset[word] = self.words[word]
                TotalLikelihood = TotalLikelihood + subset[word]
        SubSet = {}
        SubSet.update(subset)
        keys = SubSet.keys()
        for key in keys:
            if TotalLikelihood:
                SubSet[key] = SubSet[key]/TotalLikelihood
            else:
                SubSet[key] = 0
        return SubSet
    
    def getWords(self):
        words = {}
        words.update(self.words)
        return words
    
    def setWords(self,words):
        self.words = {}
        self.words.update(words)

    def getWordLength(self):
        return self.wlength

