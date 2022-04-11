import os
import pickle
import CellModeller
import csv

desiredPath = "/Users/mayinyin/CellModeller/data/FinalSpatialPositioning_ProducerFront/RecipientFront5"
with open("/Users/mayinyin/CellModeller/output_csv/SpatialPositioning/Resized/ProducerDonor/RecipientFront5.csv","w") as f:

    writer = csv.writer(f)
    writer.writerow(["Step","HGTevents","TotalCell","NumofTrans","NumofLineage"])
    allfilelist = sorted(os.listdir(desiredPath))
    ini_cell = 437
    for filename in allfilelist:
        if filename.endswith(".pickle"):
            ohnepickle = filename.split(".")[0]
            Step = ohnepickle.split("step-00")[1]
            filepath = os.path.join(desiredPath, filename)
            data = pickle.load(open(filepath,"rb"))
            HGT = data["HGTevents"] #### count the number of HGT events
            cs = data["cellStates"]
            lineage = data["lineage"]
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

            TotalCell = len(cn) ##### to count the total number of cells in simulations.
            NumofTrans = cn.count(2) ##### to count the number of transconjugants
            NumofLineage = len(set(type2id))

            writer.writerow([Step,HGT,TotalCell,NumofTrans,NumofLineage])
