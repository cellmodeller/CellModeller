import pickle
import CellModeller
import csv
import statistics
import numpy as np
import os

desiredPath = "/Users/mayinyin/CellModeller/data/Dissertation_Metabolic_Figure3"

def mylistdir(directory):
    """A specialized version of os.listdir() that ignores files that
    start with a leading period."""
    filelist = os.listdir(directory)
    return [x for x in filelist
            if not (x.startswith('.'))]

allfolderlist = mylistdir(desiredPath)
for folder in allfolderlist:
    folderpath = os.path.join(desiredPath,folder)
        #allfilelist = sorted(os.listdir(folderpath))
    allfilelist = sorted(os.listdir(os.path.join(desiredPath,folder)))

    for filename in allfilelist:
        if filename.endswith(".pickle"):
            ohnepickle = filename.split(".")[0]
            Step = ohnepickle.split("step-0")[1]
            filepath = os.path.join(desiredPath, folder,filename)
            data = pickle.load(open(filepath,"rb"))

            with open("/Users/mayinyin/CellModeller/output_csv/RP4Paper_Figure2/IM_Overtime/" + str(folder) + "_" + str(Step)+ ".csv","w") as f:
                writer = csv.writer(f)
                #writer.writerow(["Step","HGTevents","TotalCell","NumofTrans","NumofLineage"])
                writer.writerow(["FocalCellNumber","FocalCellType","NumofNeighbour","IM_addsup","IM_index"])
                #HGT = data["HGTevents"] #### count the number of HGT events
                cs = data["cellStates"]

                numofnb = [] ## to calculate how many "neighbours" each cell has
                focal_ct = [] ## get the list of focal cell number
                focal_cellType = [] ## get a list of focal cellType
                IM_each = [] ## sum it up for each cell.

                for it in cs:
                    focal_ct.append(it)
                    ct = cs[it].cellType
                    focal_cellType.append(ct)

                    cn = cs[it].neighbours
                    numofnb_percell = len(cn)
                    numofnb.append(numofnb_percell)

                    initial = 0
                    for nb in cn:
                        nbcellType = cs[nb].cellType
                        if ct == 1:
                            if ct == nbcellType:
                                initial = initial - 1
                            else:
                                initial = initial + 1
                        if ct == 0:
                            if ct == nbcellType or ct == 2: ## 0 and 2 are original the same cell type, 2 is transconjugant of 0
                                initial = initial - 1
                            else:
                                initial = initial + 1
                        if ct == 2:
                            if ct == nbcellType or ct == 0:
                                initial = initial - 1
                            else:
                                initial = initial + 1

                    IM_each.append(initial)

                lt = len(focal_ct) + 1
                normalised_IM_each = list(np.array(IM_each)/np.array(numofnb))
                ## write out
                for i in range(1,lt):
                    r_focal_ct = focal_ct[i-1]
                    r_focal_cellType = focal_cellType[i-1]
                    r_numofnb = numofnb[i-1]
                    r_im_quant =IM_each[i-1]
                    r_Ave_im_each = normalised_IM_each[i-1]
                    writer.writerow([r_focal_ct,r_focal_cellType,r_numofnb,r_im_quant,r_Ave_im_each])
