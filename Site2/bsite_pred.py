#!/usr/bin/env python3

vdw_R={"C":1.7,"A":1.7,"O":1.52,"F":1.47,"CL":1.75,"H":1.20,"N":1.55,"I":1.98,"S":1.80,"P":1.8,\
       "BR":1.85,"MG":1.73,"ZN":1.39,"MN":1.7,"CA":2.31,"NA":2.27,"FE":1.65,"NI":1.63,"CU":1.40,\
       "CO":1.50,"RU":2.0,"SE":1.9}
vdw_R_default = 1.8
CUTOFF = 0.5

def predict_bsite(prot, lig):
    min_d, contact = identify_contact_residue(prot, lig)
    wrt = []
    wrt.append("REMARK  Protein-ligand binding site prediction by GalaxySite\n")
    wrt.append("REMARK  BINDING : Ligand binding site residues\n")
    wrt.append("REMARK  DISTANCE: Distance criteria for binding site residue selection\n")
    wrt.append("REMARK            D = (minimum distance between ligand atom and each residue)\n")
    wrt.append("REMARK                - (sum of vdw radii of both atoms + 0.5A)\n")
    for resNo, resName in contact:
        wrt.append("BINDING   %4d  %3s\n"%(resNo, resName))
    resNo_s = prot.get_residue_number()
    resName_s = prot.get_residue_name()
    for i_res, resNo in enumerate(resNo_s):
        wrt.append('DISTANCE  %4d  %3s  %8.3f\n'%(resNo, resName_s[i_res], min_d[resNo]))
    return wrt

def identify_contact_residue(prot, lig):
    from Galaxy.core.vector import distance
    min_d = {}
    contact = []
    for residue in prot.get_residues():
        ds = []
        for ia in range(len(residue)):
            atm_type = residue._atmName[ia][0]
            if atm_type in vdw_R:
                vdw1 = vdw_R[atm_type]
            else:
                vdw1 = vdw_R_default
            for atom in lig.get_atoms():
                atm_type = atom.mol2_type.split('.')[0].upper()
                if atm_type in vdw_R:
                    vdw2 = vdw_R[atm_type]
                else:
                    vdw2 = vdw_R_default
                crit = distance(residue._R[ia], atom.R) - (vdw1+vdw2) - CUTOFF
                ds.append(crit)
                if (crit > 0.0): continue
                if (residue.resNo(), residue.resName()) not in contact:
                    contact.append((residue.resNo(), residue.resName()))
        min_d[residue.resNo()] = min(ds)
    contact.sort()
    return min_d, contact


