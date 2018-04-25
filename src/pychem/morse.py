# -*- coding: utf-8 -*-
"""
##############################################################################

The calculation of 3D morse descriptors. You can get 210 molecular

decriptors. You can freely use and distribute it. If you hava  

any problem, you could contact with us timely!

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.11.13

Email: oriental-cds@163.com

##############################################################################
"""

from GeoOpt import _ReadCoordinates
from AtomProperty import GetRelativeAtomicProperty

import pybel
import scipy

import math


Version=1.0
###set the parameters in RDF equation
_beta=100

def _GetR(n=32):
    """
    #################################################################
    *Internal Use Only*
    
    Obtain the parameter R in RDF equation.
    #################################################################
    """
    R=[]
    for i in range(1,n+1):
        R.append(float(i*1))
    return R


def _GetAtomDistance(x,y):
    """
    #################################################################
    *Internal Use Only*
    
    Obtain the Elucidian distance based on the coordinates of two atoms
    #################################################################
    """

    temp=[math.pow(x[0]-y[0],2),math.pow(x[1]-y[1],2),math.pow(x[2]-y[2],2)]
    res=math.sqrt(sum(temp))
    return res
    

def _GetGementricalDistanceMatrix(CoordinateList):
    """
    #################################################################
    *Internal Use Only*
    
    Obtain the distance matrix of a molecule based on coordinate list
    #################################################################
    """
    NAtom=len(CoordinateList)
    DistanceMatrix=scipy.zeros((NAtom,NAtom))
    for i in range(NAtom-1):
        for j in range(i+1,NAtom):
            DistanceMatrix[i,j]=_GetAtomDistance(CoordinateList[i],CoordinateList[j])
            DistanceMatrix[j,i]=DistanceMatrix[i,j]
    return DistanceMatrix

    
def CalculateUnweightMoRSE(ChargeCoordinates):
    """
    #################################################################
    The calculation of  unweighted 3-D MoRse descriptors 
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
#    ChargeCoordinates=_ReadCoordinates('temp.arc')
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'U'+str(kkk+1)]=round(res,3)
        
    return RDFresult

def CalculateChargeMoRSE(ChargeCoordinates):
    
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic charge.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    charge=[]
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        charge.append(float(i[4]))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+charge[j]*charge[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'C'+str(kkk+1)]=round(res,3)
        
    return RDFresult


def CalculateMassMoRSE(mol,ChargeCoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic mass.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    mass=[i.atomicmass for i in mol.atoms]
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+mass[j]*mass[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'M'+str(kkk+1)]=round(res/144,3)
        
    return RDFresult    


def CalculateAtomicNumberMoRSE(mol,ChargeCoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic number.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    mass=[i.atomicnum for i in mol.atoms]
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+mass[j]*mass[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'N'+str(kkk+1)]=round(res/144,3)
        
    return RDFresult       




def CalculatePolarizabilityMoRSE(ChargeCoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic polarizablity.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    polarizability=[]
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        polarizability.append(GetRelativeAtomicProperty(i[0],'alapha'))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+polarizability[j]*polarizability[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'P'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def CalculateSandersonElectronegativityMoRSE(ChargeCoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic sanderson electronegativity.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    En=[]
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        En.append(GetRelativeAtomicProperty(i[0],'En'))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+En[j]*En[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'E'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def CalculateVDWVolumeMoRSE(ChargeCoordinates):
    """
    #################################################################
    The calculation of  3-D MoRse descriptors 
    
    based on atomic van der Waals volume.
    #################################################################
    """
    R=_GetR(n=30)
    temp=[]
    VDW=[]
    for i in ChargeCoordinates:
        #if i[0]!='H':
        temp.append([float(i[1]),float(i[2]),float(i[3])])
        VDW.append(GetRelativeAtomicProperty(i[0],'V'))
    DM=_GetGementricalDistanceMatrix(temp)
    nAT=len(temp)
    RDFresult={}
    
    for kkk,Ri in enumerate(R):        
        res=0.0
        for j in range(nAT-1):
            for k in range(j+1,nAT):
                res=res+VDW[j]*VDW[k]*math.sin(Ri*DM[j,k])/(Ri*DM[j,k])
        RDFresult['MoRSE'+'V'+str(kkk+1)]=round(res,3)
        
    return RDFresult



def GetMoRSEUnweighted(mol):
    
    """
    #################################################################
    Obtain all unweighted 3D-Morse descriptors .
    #################################################################
    """
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateUnweightMoRSE(ChargeCoordinates)
    
    return result


def GetMoRSECharge(mol):
    
    """
    #################################################################
    Obtain all 3D-Morse descriptors baed on charge schems.
    #################################################################
    """

    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateChargeMoRSE(ChargeCoordinates)
    
    return result

    
    
def GetMoRSEMass(mol):
    
    """
    #################################################################
    Obtain all 3D-Morse descriptors baed on mass schems.
    #################################################################
    """
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateMassMoRSE(mol,ChargeCoordinates)
    
    return result



    
def GetMoRSEAtomicNumber(mol):
    
    """
    #################################################################
    Obtain all 3D-Morse descriptors baed on atomic number schems.
    #################################################################
    """
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateAtomicNumberMoRSE(mol,ChargeCoordinates)
    
    return result
    
    
    
    
def GetMoRSEPolarizability(mol):
    
    """
    #################################################################
    Obtain all 3D-Morse descriptors baed on polarizability schems.
    #################################################################
    """
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculatePolarizabilityMoRSE(ChargeCoordinates)
    
    return result
    
    
    
def GetMoRSESandersonElectronegativity(mol):
    
    """
    #################################################################
    Obtain all 3D-Morse descriptors baed on Sanderson Electronegativity schems.
    #################################################################
    """
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateSandersonElectronegativityMoRSE(ChargeCoordinates)
    
    return result



def GetMoRSEVDWVolume(mol):
    
    """
    #################################################################
    Obtain all 3D-Morse descriptors baed on VDW Volume schems.
    #################################################################
    """
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result=CalculateVDWVolumeMoRSE(ChargeCoordinates)
    
    return result    
    
def GetMoRSE(mol):
    
    """
    #################################################################
    Obtain all 3D-Morse descriptors baed on different weighted schems.
    #################################################################
    """
    result={}
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename) 
    result.update(CalculateUnweightMoRSE(ChargeCoordinates))
    result.update(CalculateChargeMoRSE(ChargeCoordinates))
    result.update(CalculateMassMoRSE(mol,ChargeCoordinates))
    result.update(CalculateAtomicNumberMoRSE(mol,ChargeCoordinates))
    result.update(CalculatePolarizabilityMoRSE(ChargeCoordinates))
    result.update(CalculateSandersonElectronegativityMoRSE(ChargeCoordinates))
    result.update(CalculateVDWVolumeMoRSE(ChargeCoordinates))
     
    return result

def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc
    pydoc.writedoc('morse')
############################################################################
if __name__=="__main__":
    

    from GeoOpt import GetARCFile
    mol='C1C=CCCS1'
    inputmol=pybel.readstring('smi',mol)  
    GetARCFile(inputmol)
    #filename='temp'
    ChargeCoordinates=_ReadCoordinates()
    print CalculateVDWVolumeMoRSE(ChargeCoordinates)
    print len(GetMoRSE(inputmol))