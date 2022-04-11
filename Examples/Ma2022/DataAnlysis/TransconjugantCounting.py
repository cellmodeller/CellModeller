### To count the total number of cells.

import pickle
import CellModeller
import pandas as pd

data = pickle.load(open("/Users/mayinyin/CellModeller/data/CF_model_resize-21-09-08-19-09/step-00000.pickle","rb"))
cs = data["cellStates"]
lineage = data["lineage"]
HGT = data["HGTevents"]
it=iter(cs)

cn = []
type2id = []

for it in cs:
    ct = cs[it].cellType
    cn.append(ct)
    if ct == 2:
        ci = cs[it].id
        while ci > ini_cell:
            ci = lineage[ci]
        type2id.append(ci)

TotalCell = len(cn)
NumofTrans = cn.count(2) ### To count the number of transconjugants

NumofLineage = len(set(type2id))
print("Number of Lineage is:" + str(NumofLineage))
#AveLineageSize = len(type2id)/NumofLineage
#print("Average size of Lineage is:" + str(AveLineageSize))



with open("/Users/mayinyin/CellModeller/output_csv/SpatialPositioning/test.csv","w") as f:
    writer = csv.writer(f)
    writer.writerow(["TotalCell","Transconjugants","HGTevents","NumofLineage"])
    for num in nbisdiff:
        writer.writerow([num])
