# -*- coding: utf-8 -*-
"""
##############################################################################

The calculation of 3D RDF descriptors. You can get 180 molecular

decriptors. You can freely use and distribute it. If you hava  

any problem, you could contact with us timely!

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.11.13

Email: oriental-cds@163.com

##############################################################################
"""

from GeoOpt import _ReadCoordinates
from AtomProperty import GetRelativeAtomicProperty
import scipy
import math


Version=1.0
#########################################################################
###set the parameters in RDF equation
_beta=100
#########################################################################

def _GetR(n=30):
    """
    #################################################################
    Obtain the parameter R in RDF equation.
    #################################################################
    """
    R=[]
    for i in range(2,n+2):
        R.append(float(i*0.5))
    return R


def GetAtomDistance(x,y):
    """
    #################################################################
    Obtain the Elucidian distance based on the 
    
    coordinates of two atoms
    #################################################################
    """

    temp=[math.pow(x[0]-y[0],2),math.pow(x[1]-y[1],2),math.pow(x[2]-y[2],2)]
    res=math.sqrt(sum(temp))
    return res
    

def GetGementricalDistanceMatrix(CoordinateList):
    """
    #################################################################
    Obtain the distance matrix of a molecule based 
    
    on coordinate list
    #################################################################
    """
    NAtom=len(CoordinateList)
    DistanceMatrix=scipy.zeros((NAtom,NAtom))
    for i in range(NAtom-1):
        for j in range(i+1,NAtom):
            DistanceMatrix[i,j]=GetAtomDistance(CoordinateList[i],CoordinateList[j])
            DistanceMatrix[j,i]=DistanceMatrix[i,j]
    return DistanceMatrix

    
def CalculateUnweightRDF(ChargeCoordinates):
    """
    #################################################################
    The calculation of unweighted radial distribution 
    
    function (RDF) descriptors.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        
    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'U'+str(kkk+1)]=round(res,3)
        
    return RDFresult

def CalculateChargeRDF(ChargeCoordinates):
    
    """
    #################################################################
    The calculation of  radial distribution function 
    
    (RDF) descriptors based on atomic charge.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    Charge=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        Charge.append(float(i[4]))
        
    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+Charge[j]*Charge[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'C'+str(kkk+1)]=round(res,3)
        
    return RDFresult


def CalculateMassRDF(mol,ChargeCoordinates):
    """
    #################################################################
    The calculation of radial distribution function (RDF) 
    
    descriptors based on atomic mass.
    #################################################################
    """
    mol.addh()
    mass=[i.atomicmass for i in mol.atoms]
    R=_GetR(n=30)
    temp=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        
    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+mass[j]*mass[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'M'+str(kkk+1)]=round(res/144,3)
        
    return RDFresult
 
        

def CalculatePolarizabilityRDF(ChargeCoordinates):
    """
    #################################################################
    The calculation of  radial distribution function 
    
    (RDF) descriptors based on atomic polarizability.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    polarizability=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        polarizability.append(GetRelativeAtomicProperty(i[0],'alapha'))
        
    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+polarizability[j]*polarizability[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'P'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def CalculateSandersonElectronegativityRDF(ChargeCoordinates):
    """
    #################################################################
    The calculation of  radial distribution function 
    
    (RDF) descriptors based on atomic electronegativity.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    EN=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        EN.append(GetRelativeAtomicProperty(i[0],'En'))
        
    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+EN[j]*EN[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'E'+str(kkk+1)]=round(res,3)
        
    return RDFresult


def CalculateVDWVolumeRDF(ChargeCoordinates):
    """
    #################################################################
    The calculation of  radial distribution function 
    
    (RDF) descriptors based on atomic van der Waals volume.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    VDW=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        VDW.append(GetRelativeAtomicProperty(i[0],'V'))
        
    DM=GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+VDW[j]*VDW[k]*math.exp(-_beta*math.pow(Ri-DM[j,k],2))
        RDFresult['RDF'+'V'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def GetRDFUnweighed(mol):
    
    """
    #################################################################
    Obtain all Unweighed radial distribution function descriptors.
    #################################################################
    """

    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateUnweightRDF(ChargeCoordinates)
     
    return result



def GetRDFCharge(mol):
    
    """
    #################################################################
    Obtain all radial distribution function descriptors based 
    
    on Charge schems.
    #################################################################
    """

    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateChargeRDF(ChargeCoordinates)
     
    return result
    
    
def GetRDFMass(mol):
    
    """
    #################################################################
    Obtain all radial distribution function descriptors based 
    
    on Mass schems.
    #################################################################
    """

    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateMassRDF(mol,ChargeCoordinates)
     
    return result
    
    
def GetRDFPolarizability(mol):
    
    """
    #################################################################
    Obtain all radial distribution function descriptors based 
    
    on Polarizability schems.
    #################################################################
    """

    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculatePolarizabilityRDF(ChargeCoordinates)
     
    return result



def GetRDFSandersonElectronegativity(mol):
    
    """
    #################################################################
    Obtain all radial distribution function descriptors based 
    
    onSanderson Electronegativity schems.
    #################################################################
    """

    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateSandersonElectronegativityRDF(ChargeCoordinates)
     
    return result


def GetRDFVDWVolume(mol):
    
    """
    #################################################################
    Obtain all radial distribution function descriptors based 
    
    on VDW Volume schems.
    #################################################################
    """

    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateVDWVolumeRDF(ChargeCoordinates)
     
    return result



def GetRDF(mol):
    
    """
    #################################################################
    Obtain all radial distribution function descriptors based 
    
    on different weighted schems.
    #################################################################
    """
    result={}
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result.update(CalculateUnweightRDF(ChargeCoordinates))
    result.update(CalculateChargeRDF(ChargeCoordinates))
    result.update(CalculateMassRDF(mol,ChargeCoordinates))
    result.update(CalculatePolarizabilityRDF(ChargeCoordinates))
    result.update(CalculateSandersonElectronegativityRDF(ChargeCoordinates))
    result.update(CalculateVDWVolumeRDF(ChargeCoordinates))
     
    return result


def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc
    pydoc.writedoc('rdf')
############################################################################
if __name__=="__main__":

    import pybel
    from GeoOpt import GetARCFile
    mol='C1C=CCCS1'
    inputmol=pybel.readstring('smi',mol)  
    GetARCFile(inputmol)
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename)
    res=CalculateVDWVolumeRDF(ChargeCoordinates)
    print res
    print len(GetRDF(inputmol))
    
