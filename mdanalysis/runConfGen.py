#! /home/omert/miniconda3/bin/python

from gmx_md_utils import *

import sys, os, shutil
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina


def mainCalcRMS(mols_path, ref_path):

    ref_mol = Chem.MolFromMolFile(ref_path)
    suppl = Chem.SDMolSupplier(mols_path)
    for i, mol in enumerate(suppl):
        print("RMSD for %d. mol" %i, calcRMS(mol, ref_mol))



def mainGenConf():
    mol_path = "Fabienne_project/Aa12b12BCD/ZnPc_SCH3_Aa1a2b1b2BCD_1.mol2"
    fileBase = mol_path.split("/")[-1].replace(".mol", "")
    print(fileBase)
    numConfs = 500
    maxAttempts = 5000
    pruneRmsThresh = 0.3
    clusterMethod = "RMSD"
    clusterThreshold = 2.0
    minimizeIterations = 0
    #  suppl = Chem.ForwardSDMolSupplier(input_file)
    print(mol_path)
    #  suppl = Chem.MolFromMolFile(mol_path)
    suppl = Chem.MolFromMol2File(mol_path)
    i=0

    xyzDIR = "xyz"
    if os.path.exists(xyzDIR):
        shutil.rmtree(xyzDIR)
    os.mkdir(xyzDIR)
    sdfDIR = "sdf"
    if os.path.exists(sdfDIR):
        shutil.rmtree(sdfDIR)
    os.mkdir(sdfDIR)

    for mol in [suppl]:
        i = i+1
        if mol is None: continue
        m = Chem.AddHs(mol, addCoords=True)
        # generate the confomers
        conformerIds = genGonformers(m, numConfs, maxAttempts, pruneRmsThresh, True, True, True)
        # align conformers
        #  AllChem.AlignMolConformers(m, conformerIds)
        conformerPropsDict = {}
        for j, conformerId in enumerate(conformerIds):
            conf_file_base = fileBase + "_conf_" + str(j)
            writeConf2sdf(m, "%s/%s.sdf" % (sdfDIR, conf_file_base), conformerId)
            mol = read("%s/%s.sdf" % (sdfDIR, conf_file_base))
            write("%s/%s.xyz" % (xyzDIR, conf_file_base), mol)
            # energy minimise (optional) and energy calculation
            #  props = calcEnergy(m, conformerId, minimizeIterations)
            #  conformerPropsDict[conformerId] = props
        #  # cluster the conformers
        #  rmsClusters = getClusterConf(m, clusterMethod, clusterThreshold)
        #  print("Molecule", i, ": generated", len(conformerIds), "conformers and", len(rmsClusters), "clusters")
        #  rmsClustersPerCluster = []
        #  clusterNumber = 0
        #  minEnergy = 9999999999999
        #  for cluster in rmsClusters:
        #      clusterNumber = clusterNumber+1
        #      rmsWithinCluster = alignConfs(m, cluster)
        #      for conformerId in cluster:
        #          e = props["energy_abs"]
        #          if e < minEnergy:
        #              minEnergy = e
        #          props = conformerPropsDict[conformerId]
        #          props["cluster_no"] = clusterNumber
        #          props["cluster_centroid"] = cluster[0] + 1
        #          idx = cluster.index(conformerId)
        #          if idx > 0:
        #              props["rms_to_centroid"] = rmsWithinCluster[idx-1]
        #          else:
        #              props["rms_to_centroid"] = 0.0
        #  writeMinEConf2sdf(m, "target_conf_" + str(i) + ".sdf", rmsClusters, conformerPropsDict, minEnergy)

mainGenConf()
