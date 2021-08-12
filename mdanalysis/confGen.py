import sys, os, shutil
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina

from ase.io import read, write

def genGonformers(mol, numConfs=100,
                   maxAttempts=1000,
                   pruneRmsThresh=0.1,
                   useExpTorsionAnglePrefs=True,
                   useBasicKnowledge=True,
                   enforceChirality=True):
	confs = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs,
                                    maxAttempts=maxAttempts,
                                    pruneRmsThresh=pruneRmsThresh,
                                    useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
                                    useBasicKnowledge=useBasicKnowledge,
                                    enforceChirality=enforceChirality, numThreads=0,
                                   )
	return list(confs)

def writeMinEConf2sdf(mol, filename, rmsClusters, conformerPropsDict, minEnergy):
	w = Chem.SDWriter(filename)
	for cluster in rmsClusters:
		for confId in cluster:
			for name in mol.GetPropNames():
				mol.ClearProp(name)
			conformerProps = conformerPropsDict[confId]
			mol.SetIntProp("conformer_id", confId + 1)
			for key in conformerProps.keys():
				mol.SetProp(key, str(conformerProps[key]))
			e = conformerProps["energy_abs"]
			if e:
				mol.SetDoubleProp("energy_delta", e - minEnergy)
			w.write(mol, confId=confId)
	w.flush()
	w.close()

def writeConf2sdf(mol, filename, confId):
	w = Chem.SDWriter(filename)
	w.write(mol, confId)
	w.close()

def calcEnergy(mol, conformerId, minimizeIts):
	ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId)
	ff.Initialize()
	ff.CalcEnergy()
	results = {}
	if minimizeIts > 0:
		results["converged"] = ff.Minimize(maxIts=minimizeIts)
	results["energy_abs"] = ff.CalcEnergy()
	return results

def getClusterConf(mol, mode="RMSD", threshold=2.0):
	if mode == "TFD":
		dmat = TorsionFingerprints.GetTFDMatrix(mol)
	else:
		dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
	rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
	return rms_clusters

def alignConfs(mol, clust_ids):
	rmslist = []
	AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
	return rmslist

def calcRMS(mol, ref_mol):
    rms = Chem.rdMolAlign.CalcRMS(mol, ref_mol)
    return rms

#  mol = Chem.MolFromMolFile("./mutemel/article_No1_fin.mol")
#  ref_mol = Chem.MolFromMolFile("./mutemel/article_No1_fin_conf_20.sdf")
#  print(calcRMS(mol, ref_mol))

