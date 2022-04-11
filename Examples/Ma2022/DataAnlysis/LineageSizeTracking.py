import pickle
import CellModeller
import csv
import statistics
import numpy
import os

desiredPath = "/Users/mayinyin/CellModeller/data/RP4SpatialPositioning_finalized/GroupC_130Step/"

def mylistdir(directory):
    """A specialized version of os.listdir() that ignores files that
    start with a leading period."""
    filelist = os.listdir(directory)
    return [x for x in filelist
            if not (x.startswith('.'))]

allfolderlist = mylistdir(desiredPath)
for folder in allfolderlist:
    with open("/Users/mayinyin/CellModeller/output_csv/SpatialPositioning/lineage_size/RawData/" + str(folder) + ".csv","w") as f:
        writer = csv.writer(f)
        writer.writerow(["Step","TotalCell","NumofTrans","NumofLineage","MeanSize","MedianSize"])
        folderpath = os.path.join(desiredPath,folder)
        #allfilelist = sorted(os.listdir(folderpath))
        allfilelist = sorted(os.listdir(os.path.join(desiredPath,folder)))

        for filename in allfilelist:
            if filename.endswith(".pickle"):
                ohnepickle = filename.split(".")[0]
                Step = ohnepickle.split("step-0")[1]
                filepath = os.path.join(desiredPath, folder,filename)
                data = pickle.load(open(filepath,"rb"))
                #HGT = data["HGTevents"] #### count the number of HGT events
                cs = data["cellStates"]
                lineage = data["lineage"]
                ini_cell = 437
                it=iter(cs)
                n = len(cs)
                type2id = []
                cn = [] ## cell number

                i = 0
                for it in cs:
                    ct = cs[it].cellType
                    cn.append(ct)
                    if ct == 2:
                        i = i + 1
                        ci = cs[it].id
                        while ci > ini_cell:
                            ci = lineage[ci]
                        type2id.append(ci)
                NumofTrans = i
                lengthOfLineage = []
                uniqueLineage = list(set(type2id))
                for u in uniqueLineage:
                    lengthOfu = type2id.count(u)
                    lengthOfLineage.append(lengthOfu)

                TotalCell = len(cn)
                lengthOfLineage = sorted(lengthOfLineage) ### Sort the list from small to big
                median = statistics.median(lengthOfLineage)
                ####Calculate the mean/the average lineage size
                NumofLineage = len(set(type2id))
                ## print("Number of Lineage is:" + str(NumofLineage))
                AveLineageSize = round(len(type2id)/NumofLineage)
                ##print("Average size of Lineage is:" + str(AveLineageSize))
                ##print("Median of lineage size is:" + str(median))
                writer.writerow([Step,TotalCell,NumofTrans,NumofLineage,AveLineageSize,median])
