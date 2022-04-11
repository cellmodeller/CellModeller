import pickle
import CellModeller
import csv
import statistics
import numpy as np
import os

desiredPath = "/Users/mayinyin/CellModeller/data/RP4SpatialPositioning_finalized/Final_SpatialPositioning_ProducerDonor"

def mylistdir(directory):
    """A specialized version of os.listdir() that ignores files that
    start with a leading period."""
    filelist = os.listdir(directory)
    return [x for x in filelist
            if not (x.startswith('.'))]

allfolderlist = mylistdir(desiredPath)
for folder in allfolderlist:
    with open("/Users/mayinyin/CellModeller/output_csv/SpatialPositioning/Resized/IM_Step0/ProducerFront/" + str(folder) +".csv","w") as f:
        writer = csv.writer(f,delimiter=",")
        #writer.writerow(["Step","HGTevents","TotalCell","NumofTrans","NumofLineage"])
        writer.writerow(["FocalCellNumber","FocalCellType","NumofNeighbour","IM_addsup","IM_index"])
        folderpath = os.path.join(desiredPath,folder)
        #allfilelist = sorted(os.listdir(folderpath))
        allfilelist = sorted(os.listdir(os.path.join(desiredPath,folder)))
        for filename in allfilelist:
            if filename.endswith("00000.pickle"):
                filepath = os.path.join(desiredPath, folder,filename)
                data = pickle.load(open(filepath,"rb"))
                #HGT = data["HGTevents"] #### count the number of HGT events
                cs = data["cellStates"]
                lineage = data["lineage"]


                a_step = [-5, +5, 0]
                b_step = [0, +5, 0]
                c_step = [+5, +5, 0]
                d_step = [-5, 0, 0]
                e_step = [+5, 0, 0]
                f_step = [-5, -5, 0]
                g_step = [0, -5, 0]
                h_step = [+5, -5, 0]
                numofnb = [] ## to calculate how many "neighbours" each cell has
                focal_ct = [] ## get the list of focal cell number
                focal_cellType = []
                for it in cs:
                    cp = cs[it].pos
                    focal_ct.append(it)
                    ct = cs[it].cellType
                    focal_cellType.append(ct)
                    p_a = list(map(lambda i, j: i + j, cp, a_step)) # top left
                    p_b = list(map(lambda i, j: i + j, cp, b_step)) # right above
                    p_c = list(map(lambda i, j: i + j, cp, c_step)) # top right
                    p_d = list(map(lambda i, j: i + j, cp, d_step)) # left
                    p_e = list(map(lambda i, j: i + j, cp, e_step)) # right
                    p_f = list(map(lambda i, j: i + j, cp, f_step)) # bottom left
                    p_g = list(map(lambda i, j: i + j, cp, g_step)) # right below
                    p_h = list(map(lambda i, j: i + j, cp, h_step)) # bottom right

                    nb_each = []
                    nb_coord = [p_a,p_b,p_c,p_d,p_e,p_f,p_g,p_h]
                    for it in cs:
                        for each_nb in nb_coord:
                            if cs[it].pos == each_nb:
                                nb_each.append(it)
                    numofnb.append(nb_each)

                All_focal_nb = dict(zip(focal_ct,numofnb)) ### create a dictionary of focal cell and its neighbours.
                ### Next step to calculate the number of neighbours each cell has.
                lengthofnb = []
                for nbs in numofnb:
                    lengthofnb_each = len(nbs)
                    lengthofnb.append(lengthofnb_each)

                ### Next to identify the cellType each neighbour has
                lt = len(focal_ct) + 1
                im_all = [] ### if the neighbour has the same cellType then -1, otherwise +1
                for i in range(1,lt):
                    im_each = 0
                    f = focal_ct[i-1]
                    m = cs[f].cellType
                    nb = numofnb[i-1]
                    for sg_nb in nb:
                        n = cs[sg_nb].cellType
                        if m == n:
                            m_each = im_each - 1
                        else:
                            im_each = im_each + 1
                    im_all.append(im_each)

                Ave_im_each = list(np.array(im_all)/np.array(lengthofnb)) ### average to each focal cell
                Ave_im = sum(Ave_im_each)/len(focal_ct)
                dataset = [focal_ct,focal_cellType,lengthofnb,Ave_im_each]

                for i in range(1,lt):
                    r_focal_ct = focal_ct[i-1]
                    r_focal_cellType = focal_cellType[i-1]
                    r_lengthofnb = lengthofnb[i-1]
                    r_im_quant =im_all[i-1]
                    r_Ave_im_each = Ave_im_each[i-1]
                    writer.writerow([r_focal_ct,r_focal_cellType,r_lengthofnb,r_im_quant,r_Ave_im_each])
