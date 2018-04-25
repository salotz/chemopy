# -*- coding: utf-8 -*-
"""
This is a test suite for chemopy using 40 structurally diverse molecules.

It takes several minutes to calculate 40 molecules due to the molecular 

optimization using MOPAC. 

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2013.01.25

Email: oriental-cds@163.com
"""

from pychem import pychem


def dotest(smi):
    """
    Test all modules in ChemoPy
    """
    #smi='CC(=O)OC1=CC=CC=C1C(=O)O'
    
    mol=pychem.getmol.ReadMolFromSmile(smi)
    
    print "********************************************************"
    try:
        temp=pychem.constitution.GetConstitutional(mol)
    except:
        print "Constitutional module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "Constitutional module works!!!"
        print "The number of constitutional descriptors is", len(temp)
    print "********************************************************"
    
    
    print "********************************************************"     
    try:
        temp=pychem.connectivity.GetConnectivity(mol)
    except:
        print "Connectivity module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "Connectivity module works!!!"
        print "The number of connectivity descriptors is", len(temp)
    
    print "********************************************************"  
    print "********************************************************" 

    try:
        temp=pychem.basak.Getbasak(mol)
    except:
        print "basak module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "basak module works!!!"
        print "The number of basak descriptors is", len(temp)    
    print "********************************************************" 
    
    
    
    print "********************************************************" 
    try:
        temp=pychem.bcut.GetBurden(mol)
    except:
        print "bcut module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "bcut module works!!!"
        print "The number of Burden descriptors is", len(temp)    
    print "********************************************************" 
    
    
    
    print "********************************************************" 
    try:
        temp=pychem.charge.GetCharge(mol)
    except:
        print "charge module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "charge module works!!!"
        print "The number of charge descriptors is", len(temp)
    print "********************************************************" 
    
    
    print "********************************************************" 
    try:
        temp=pychem.estate.GetEstate(mol)
    except:
        print "estate module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "estate module works!!!"
        print "The number of Estate descriptors is", len(temp)
    print "********************************************************"  
    


    print "********************************************************" 
    try:
        temp=pychem.kappa.GetKappa(mol)
    except:
        print "Kappa module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "Kappa module works!!!"
        print "The number of Kappa descriptors is", len(temp)
    print "********************************************************" 



    print "********************************************************" 
    try:
        temp=pychem.moe.GetMOE(mol)
    except:
        print "MOE module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "MOE module works!!!"
        print "The number of MOE descriptors is", len(temp)
    print "********************************************************" 



    print "********************************************************" 
    try:
        temp=pychem.molproperty.GetMolecularProperty(mol)
    except:
        print "MolecularProperty module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "MolecularProperty module works!!!"
        print "The number of Molecular Properties is", len(temp)
    print "********************************************************" 




    print "********************************************************" 
    try:
        temp=pychem.moran.GetMoranAuto(mol)
    except:
        print "moran module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "moran module works!!!"
        print "The number of Moran autocorrelation descriptors is", len(temp)
    print "********************************************************" 



    print "********************************************************" 
    try:
        temp=pychem.moreaubroto.GetMoreauBrotoAuto(mol)
    except:
        print "moreaubroto module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "moreaubroto module works!!!"
        print "The number of Moreau-Broto autocorrelation descriptors is", len(temp)
    print "********************************************************" 



    print "********************************************************" 
    try:
        temp=pychem.geary.GetGearyAuto(mol)
    except:
        print "geary module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "geary module works!!!"
        print "The number of Geary autocorrelation descriptors is", len(temp)
    print "********************************************************" 



    print "********************************************************" 
    try:
        temp=pychem.topology.GetTopology(mol)
    except:
        print "topology module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "topology module works!!!"
        print "The number of topology descriptors is", len(temp)
    print "********************************************************"   
    
    
    print "********************************************************" 
    
    mol1=pychem.pybel.readstring('smi',smi)
    
    
    try:
        pychem.GetARCFile(mol1)
        temp=pychem.whim.GetWHIM()
    except:
        print "whim module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "whim module works!!!"
        print "The number of whim descriptors is", len(temp)
    print "********************************************************" 


    print "********************************************************" 
    
    mol1=pychem.pybel.readstring('smi',smi)
    
    
    try:
        pychem.GetARCFile(mol1)
        temp=pychem.cpsa.GetCPSA()
    except:
        print "cpsa module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "cpsa module works!!!"
        print "The number of cpsa descriptors is", len(temp)
    print "********************************************************" 
    
    print "********************************************************" 
    
    mol1=pychem.pybel.readstring('smi',smi)
    
    
    try:
        pychem.GetARCFile(mol1)
        temp=pychem.geometric.GetGeometric(mol1)
    except:
        print "geometric module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "cpsa geometric works!!!"
        print "The number of geometric descriptors is", len(temp)
    print "********************************************************" 
    
    
    print "********************************************************" 
    
    mol1=pychem.pybel.readstring('smi',smi)
    
    
    try:
        pychem.GetARCFile(mol1)
        temp=pychem.morse.GetMoRSE(mol1)
    except:
        print "morse module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "morse module works!!!"
        print "The number of morse descriptors is", len(temp)
    print "********************************************************" 
    
    
    print "********************************************************" 
    
    mol1=pychem.pybel.readstring('smi',smi)
    
    
    try:
        pychem.GetARCFile(mol1)
        temp=pychem.rdf.GetRDF(mol1)
    except:
        print "rdf module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "rdf module works!!!"
        print "The number of rdf descriptors is", len(temp)
    print "********************************************************" 
    

    print "********************************************************" 
    
    try:
        temp=pychem.fingerprint.CalculateFP4Fingerprint(smi)
    except:
        print " CalculateFP4Fingerprint module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "CalculateFP4Fingerprint module works!!!"
        print "The number of FP4 fingerprints is", temp[0]
    print "********************************************************" 
    
    
    print "********************************************************" 
    
    try:
        temp=pychem.fingerprint.CalculateMACCSFingerprint(mol)
    except:
        print " CalculateMACCSFingerprint module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "CalculateMACCSFingerprint module works!!!"
        print "The number of MACCS fingerprints is", temp[0]
    print "********************************************************" 
    
    print "********************************************************" 
    
    try:
        temp=pychem.fingerprint.CalculateAtomPairsFingerprint(mol)
    except:
        print " CalculateAtomPairsFingerprint module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "CalculateAtomPairsFingerprint module works!!!"
    print "********************************************************" 
    
    print "********************************************************" 
    
    try:
        temp=pychem.fingerprint.CalculateDaylightFingerprint(mol)
    except:
        print " CalculateDaylightFingerprint module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "CalculateDaylightFingerprint module works!!!"
        print "The number of Daylight fingerprints is", temp[0]
    print "********************************************************" 

    print "********************************************************" 
    
    try:
        temp=pychem.fingerprint.CalculateEstateFingerprint(mol)
    except:
        print " CalculateEstateFingerprint module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "CalculateEstateFingerprint module works!!!"
        print "The number of Estate fingerprints is", temp[0]
    print "********************************************************" 
    
    print "********************************************************" 
    
    try:
        temp=pychem.fingerprint.CalculateMorganFingerprint(mol)
    except:
        print " CalculateMorganFingerprint module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "CalculateMorganFingerprint module works!!!"
    print "********************************************************" 
    
    
    print "********************************************************" 
    try:
        temp=pychem.fingerprint.CalculateTopologicalTorsionFingerprint(mol)
    except:
        print " CalculateTopologicalTorsionFingerprint module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "CalculateTopologicalTorsionFingerprint module works!!!"
    print "********************************************************"



    print "********************************************************" 

    try:
        from pychem.pychem import PyChem2d
        des1 = PyChem2d()
        des1.ReadMolFromSmile(smi)
        temp=des1.GetAllDescriptor()
    except:
        print " PyChem2d module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "PyChem2d module works!!!"
        print "The number of all 2D descriptors is", len(temp)
    print "********************************************************"
    
    
    print "********************************************************" 

    try:
        from pychem.pychem import PyChem3d
        des1 = PyChem3d()
        des1.ReadMol(smi)
        temp=des1.GetAllDescriptor()
    except:
        print " PyChem3d module failed!!!"
        temp=[]
        
    if temp!=[]:
        #print temp
        print "PyChem3d module works!!!"
        print "The number of all 3D descriptors is", len(temp)
    print "********************************************************"
    return 'OK'
    
    
    
if "__main__"==__name__:
    
   smi=['CC(=O)OC1=CC=CC=C1C(=O)O',
        'N[C@@H](CC1=CC=CC=C1)C(O)=O',
        'CSCC[C@H](N)C(O)=O',
        'CC[C@H](C)[C@H](N)C(O)=O',
        'OC(=O)CCCCC1CCSS1',
        'CC1=NC=C(CO)C(CO)=C1O',
        'NC1=C2NC=NC2=NC=N1',
        'OC1N=C(C2=CC=CC=C2Cl)C2=C(NC1=O)C=CC(Cl)=C2',
        'CC(C)(N)CC1=CC=CC=C1',
        'CN(C)CCCC1(OCC2=C1C=CC(=C2)C#N)C1=CC=C(F)C=C1',
        'NC1=CC(C(O)=O)=C(O)C=C1',
        'ClCCNC(=O)N(CCCl)N=O',
        'COC1=CC2=C(C=C1)C=C(CCC(C)=O)C=C2',
        'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)NC(=O)C(CC(=O)NC)N)C(=O)O)C',
        'CCCCOC(=O)C1=CC=CC=C1C(=O)OCC2=CC=CC=C2',
        'CCN(CC)CCCC(C)NC1=C2C=C(C=CC2=NC3=C1C=CC(=C3)Cl)OC',
        'C(CCC(=O)O)CCN',
        'OC(=O)CNC(=O)NC1CCCCC1',
        'CN1N=C(C2=C(N)N=CN=C12)C1=CC2=CC=CC=C2C=C1',
        'OC1=CC=CC=C1C1=CC=CC=C1[S@@](O)=O',
        'CCOCCOCCOCC',
        'N[C@@H](CSCC1=CC=C(Br)C=C1)C(O)=O',
        'OC(=O)N1CCC(CC2=CC(OC3=NC=C(C=C3)C(F)(F)F)=CC=C2)CC1',
        'C[C@H](O)[C@@H](N)CC1=CC=CC=C1',
        'OC1=CC(O)=C(C(=O)N2C=CC=C2)C(O)=C1',
        'CN(C)CC1=CNC2=C1C=CC=C2CCCO',
        'COC1=CC(OC)=CC(NC(=N)NC(N)=N)=C1',
        'CCOC(=O)NCCCNC(=O)OCC',
        'FC1=CC=C(C=C1)C1=C(N(CC2CC2)C=N1)C1=CC=NC=C1',
        'OC(=O)CCC(=O)NC1=CC=C(Br)C=N1',
        'O=C1CN(C2=CC=CC=C2)S(=O)(=O)N1',
        'CCCCCCCCSCC(=O)C(F)(F)F',
        'CN(C)C1=CC=C(C=C1)C1=[N+](C)C2=CC=C(O)C=C2S1',
        'CCO[P@](O)(=O)N(CC)CC',
        'CCC1=C(SC=C1)C1=CC(O)=C(O)C=C1',
        'NC[C@H]1CN(CCO1)C1=C(Br)C=NC2=C1C=NN2',
        'CN(C)CCN(CC1=CC=C(Cl)C=C1)C1=CC=CC=N1',
        'OC(=O)C(Cl)Cl',
        'CC(C)[C@H]1SC(NC2=C(F)C=CC=C2)=NC1=O',
        'CN(C)C1=CC=C(C=C1)C(O)=O']
   for index,i in enumerate(smi):
       print "The", index, "molecule is running ......"
       dotest(i)
