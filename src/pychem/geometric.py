# -*- coding: utf-8 -*-
"""
##############################################################################
This module is to compute geometrical descriptors based on the optimized

molecular structure by MOPAC.If you have any question please contact me

via email.

Created on Tue Apr 19 11:37:14 2011

@author: Dongsheng Cao
##############################################################################
"""

import pybel
import scipy
#import scipy.linalg
from GeoOpt import _ReadCoordinates

import math

Version=1.0
############################################################################
            
    
def GetAtomDistance(x,y):
    """
    #################################################################
    Obtain the Elucidian distance based on the coordinates of two atoms
    #################################################################
    """

    temp=[math.pow(x[0]-y[0],2),math.pow(x[1]-y[1],2),math.pow(x[2]-y[2],2)]
    res=scipy.sqrt(sum(temp))
    return res


def GetGementricalDistanceMatrix(CoordinateList):
    """
    #################################################################
    Obtain the distance matrix of a molecule based on coordinate list
    #################################################################
    """
    NAtom=len(CoordinateList)
    DistanceMatrix=scipy.zeros((NAtom,NAtom))
    for i in range(NAtom-1):
        for j in range(i+1,NAtom):
            DistanceMatrix[i,j]=GetAtomDistance(CoordinateList[i],CoordinateList[j])
            DistanceMatrix[j,i]=DistanceMatrix[i,j]
    return DistanceMatrix
 

           
def _GetMassCenter(MassCoordinates):
    
    """
    #################################################################
    Get the center of mass.
    INPUT: MassCoordinates is [[atommass,[x,y,z]],......].
    #################################################################
    """
    
    res1=0.0
    res2=0.0
    res3=0.0
    temp=[]
    for i in MassCoordinates:
        res1=res1+i[0]*i[1][0]
        res2=res2+i[0]*i[1][1]
        res3=res3+i[0]*i[1][2]
        temp.append(i[0])
    result=[res1/sum(temp),res2/sum(temp),res3/sum(temp)]
    return result
    

def _GetGeometricalCenter(ChargeCoordinates):
    """
    #################################################################
    Get the geometrical center
    #################################################################
    """
    res1=[]
    res2=[]
    res3=[]

    for i in ChargeCoordinates:
        res1.append(float(i[1]))
        res2.append(float(i[2]))
        res3.append(float(i[3]))
    
    result=[scipy.mean(res1),scipy.mean(res2),scipy.mean(res3)]
    
    return result    


def Calculate3DWienerWithH(ChargeCoordinates):
    """
    #################################################################
    The calculation of 3-D Wiener index based gemetrical distance matrix optimized
    by MOPAC(including Hs)
    -->W3DH
    #################################################################
    """
    temp=[]
    for i in ChargeCoordinates:        
        temp.append([float(i[1]),float(i[2]),float(i[3])])
    
    DistanceMatrix=GetGementricalDistanceMatrix(temp)
   
    return round(scipy.sum(DistanceMatrix)/2.0,3)


def Calculate3DWienerWithoutH(ChargeCoordinates):
    """
    #################################################################
    The calculation of 3-D Wiener index based 
    gemetrical distance matrix optimized
    by MOPAC(Not including Hs) 
    -->W3D
    #################################################################
    """
    temp=[]
    for i in ChargeCoordinates:
        if i[0]!='H':
            temp.append([float(i[1]),float(i[2]),float(i[3])])
    
    DistanceMatrix=GetGementricalDistanceMatrix(temp)
   
    return round(scipy.sum(DistanceMatrix)/2.0,3)


        
def CalculatePetitjean3DIndex(ChargeCoordinates):
    """
    #################################################################
    CalculatePetitjean Index based on molecular gemetrical distance matrix
    -->Petitj3D
    The 3D Petitjean shape index (PJI3) is calculated 
    dividing the difference between geometric diameter and 
    radius by the geometric radius [P.A. Bath, A.R. Poirrette, 
    P. Willett, F.H. Allen, J.Chem.Inf.Comput.Sci. 1995, 35, 714-716]. 
    The geometric radius of a molecule is defined as the minimum 
    geometric eccentricity and the diameter is defined as the 
    maximum geometric eccentricity in the molecule, the atom 
    geometric eccentricity being the longest geometric distance 
    from the considered atom to any other atom in the molecule. 
    #################################################################
    """
    temp=[]
    for i in ChargeCoordinates:        
        temp.append([float(i[1]),float(i[2]),float(i[3])])
    
    DistanceMatrix=GetGementricalDistanceMatrix(temp)    
    temp1=scipy.amax(DistanceMatrix,axis=0)
    
    return round(max(temp1)/min(temp1)-1.0,3)


def CalculateGemetricalDiameter(ChargeCoordinates):
    """
    #################################################################
    The longest distance between two atoms (gemetrical diameter)
    -->GeDi
    #################################################################
    """
    temp=[]
    for i in ChargeCoordinates:        
        temp.append([float(i[1]),float(i[2]),float(i[3])])

    DistanceMatrix=GetGementricalDistanceMatrix(temp)    
    temp1=scipy.amax(DistanceMatrix,axis=0)
    
    return round(max(temp1),3)


def CalculateTopoElectronic(ChargeCoordinates):
    """
    #################################################################
    #################################################################
    """
    pass



def CalculateGravitational3D1(mol,ChargeCoordinates):
    
    """
    #################################################################
    Calculation of Gravitational 3D index.
    --->grav1
    #################################################################
    """
    mol.removeh()
    mol.addh()
    temp=[]  
    for i,j in enumerate(ChargeCoordinates):        
        temp.append([mol.atoms[i].atomicmass,[float(j[1]),float(j[2]),float(j[3])]])
    
    
    nAT=len(temp)
    result=0.0
    for i in range(nAT-1):
        for j in range(i+1,nAT):
            dis=GetAtomDistance(temp[i][1],temp[j][1])
            
            result=result+temp[i][0]*temp[j][0]/scipy.power(dis,p=2)
    
    return round(float(result)/100,3)
            

def CalculateGravitational3D2((mol,ChargeCoordinates)):
    """
    #################################################################
    Gravitational indices are molecular descriptors 
    reflecting the mass distribution in a molecule, 
    defined as [A.R. Katritzky, L. Mu, V.S. Lobanov,
    M. Karelson, J.Phys.Chem. 1996, 100, 10400-10407]
    ---->grav2
    #################################################################
    """
    pass
        


def CalculateRadiusofGyration(mol,ChargeCoordinates):
    
    """
    #################################################################
    Calculation of Radius of gyration. 
    --->rygr
    #################################################################
    """
    mol.addh()
    temp=[]  
    for i,j in enumerate(ChargeCoordinates):        
        temp.append([mol.atoms[i].atomicmass,[float(j[1]),float(j[2]),float(j[3])]])
    nAT=len(temp)
    
    
    masscenter=_GetMassCenter(temp)
    result=0.0
    for i in range(nAT):
        dis=GetAtomDistance(temp[i][1],masscenter)
        result=result+temp[i][0]*scipy.power(dis,p=2)
    
    
    return round(scipy.sqrt(float(result/mol.molwt)),3)
    



def GetInertiaMatrix(mol,ChargeCoordinates):
    """
    #################################################################
    Get Inertia matrix based on atomic mass and optimized coordinates.
    #################################################################
    """
    mol.removeh()
    mol.addh()
    
    temp=[]  
    for i,j in enumerate(ChargeCoordinates):        
        temp.append([mol.atoms[i].atomicmass,[float(j[1]),float(j[2]),float(j[3])]])   
#     
#    masscenter=_GetMassCenter(temp)  
#    
#    for i,j in enumerate(temp):
#        temp[i][1]=[d-masscenter[k] for k,d in enumerate(j[1])]     
     

    nAT=len(temp)    
    
    InertiaMatrix=scipy.zeros((3,3))
    res11=0.0
    res22=0.0
    res33=0.0
    res12=0.0
    res23=0.0
    res13=0.0
    for i in range(nAT):
        res11=res11+temp[i][0]*(math.pow(temp[i][1][1],2)+math.pow(temp[i][1][2],2))
        res22=res22+temp[i][0]*(math.pow(temp[i][1][0],2)+math.pow(temp[i][1][2],2))
        res33=res33+temp[i][0]*(math.pow(temp[i][1][0],2)+math.pow(temp[i][1][1],2))
        res12=res12+temp[i][0]*(temp[i][1][0]*temp[i][1][1])
        res13=res13+temp[i][0]*(temp[i][1][0]*temp[i][1][2])
        res23=res23+temp[i][0]*(temp[i][1][1]*temp[i][1][2])
    InertiaMatrix[0,0]=res11
    InertiaMatrix[1,1]=res22
    InertiaMatrix[2,2]=res33
    InertiaMatrix[0,1]=res12
    InertiaMatrix[0,2]=res13
    InertiaMatrix[1,2]=res23
    InertiaMatrix[1,0]=res12
    InertiaMatrix[2,0]=res13
    InertiaMatrix[2,1]=res23
    
    return InertiaMatrix
    
        


def CalculatePrincipalMomentofInertia(mol,ChargeCoordinates):
    """
    #################################################################
    X,Y and Z-principal geometric moment.   
    drived from ADAPT developed by Jurs.
    #################################################################
    """
    InertiaMatrix=GetInertiaMatrix(mol,ChargeCoordinates)
    ma=scipy.mean(InertiaMatrix,axis=1)
    ms=scipy.std(InertiaMatrix,axis=1,ddof=1)
    bb=scipy.ones((3,1))
    InertiaMatrix=(InertiaMatrix-bb*ma.T)/(bb*ms.T)  
    u,s,v=scipy.linalg.svd(InertiaMatrix)

    res={}
    res['IA']=round(s[2],3)
    res['IB']=round(s[1],3)
    res['IC']=round(s[0],3)
    
    return res
    
        
def CalculateRatioPMI(mol,ChargeCoordinates):
    """
    #################################################################
    The ratio of X/Y, Y/Z and X/Z (principal moment of inertia)
    drived from ADAPT developed by Jurs.
    #################################################################
    """
    temp=CalculatePrincipalMomentofInertia(mol,ChargeCoordinates)
    res={}
        
    res['IA/B']=round(temp['IA']/temp['IB'],3)
    res['IA/C']=round(temp['IA']/temp['IC'],3)
    res['IB/C']=round(temp['IB']/temp['IC'],3)
    return res






def CalculateHarary3D(ChargeCoordinates):
    """
    #################################################################
    The 3D-Harary index (H3D) is calculated as 
    the sum of all the reciprocal geometric distances 
    in a molecule. 
    --->Harary3D
    #################################################################
    """
    temp=[]
    for i in ChargeCoordinates:       
        temp.append([float(i[1]),float(i[2]),float(i[3])])    
    DistanceMatrix=GetGementricalDistanceMatrix(temp) 
    nAT=len(temp)
    res=0.0
    for i in range(nAT-1):
        for j in range(i+1,nAT):
            if DistanceMatrix[i,j]==0:
                cds=0.0
            else:
                cds=1./DistanceMatrix[i,j]
            res=res+cds
    return round(res,3)
            
            


def CalculateAverageGeometricalDistanceDegree(ChargeCoordinates):
    """
    #################################################################
    The average geometric distance degree (AGDD) is 
    calculated dividing the sum of all geometric distance 
    degrees by the total number of molecule atoms (nAT). 
    ---->AGDD
    #################################################################
    """
    temp=[]
    for i in ChargeCoordinates:       
        temp.append([float(i[1]),float(i[2]),float(i[3])])    
    DistanceMatrix=GetGementricalDistanceMatrix(temp) 
    nAT=len(temp)
    
    res=sum(sum(DistanceMatrix))/nAT
    
    return round(res,3)
    


def CalculateAbsEigenvalueSumOnGeometricMatrix(ChargeCoordinates):
    """
    #################################################################
    The absolute eigenvalue sum on geometry matrix (SEig) 
    is the sum of the absolute eigenvalues of the geometry matrix. 
    --->SEig
    #################################################################
    """
    temp=[]
    for i in ChargeCoordinates:       
        temp.append([float(i[1]),float(i[2]),float(i[3])])    
    DistanceMatrix=GetGementricalDistanceMatrix(temp) 
    
    u,s,vt=scipy.linalg.svd(DistanceMatrix)
    
    return round(sum(abs(s)),3)


def CalculateSPANR(mol,ChargeCoordinates):
    """
    #################################################################
    The span R (SPAN) is a size descriptor defined as 
    the radius of the smallest sphere, centred on the centre 
    of mass, completely enclosing all atoms of a molecule 
    [G.A. Arteca, Molecular Shape Descriptors in Reviews in 
    Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.), 
    VCH Publishers, New York (NY), pp. 191-253, 1991]
    --->SPAN
    #################################################################
    """
    mol.removeh()
    mol.addh()
    temp=[]  
    for i,j in enumerate(ChargeCoordinates):        
        temp.append([mol.atoms[i].atomicmass,[float(j[1]),float(j[2]),float(j[3])]])   
     
    masscenter=_GetMassCenter(temp)  
    
    res=[]
    for i in temp:
        res.append(GetAtomDistance(i[1],masscenter)) 

    return round(float(max(res)),3)


def CalculateAverageSPANR(mol,ChargeCoordinates):
    """
    #################################################################
    The average span R (SPAM) is the root square of 
    the ratio of SPAN over the number of atoms.
    --->ASPAN
    #################################################################
    """
    mol.removeh()
    mol.addh()
    temp=[]  
    for i,j in enumerate(ChargeCoordinates):        
        temp.append([mol.atoms[i].atomicmass,[float(j[1]),float(j[2]),float(j[3])]])   
    
    nAT=len(temp)
    masscenter=_GetMassCenter(temp)      
    res=[]
    for i in temp:
        res.append(GetAtomDistance(i[1],masscenter)) 

    return round(math.pow(float(max(res))/nAT,0.5),3)


def CalculateMolecularEccentricity(mol,ChargeCoordinates):
    """
    #################################################################
    The molecular eccentricity (MEcc) is a shape descriptor 
    calculated from the eigenvalues l of the molecular inertia matrix 
    [G.A. Arteca, Molecular Shape Descriptors in Reviews 
    in Computational Chemistry - Vol. 9, K.B. Lipkowitz, D. Boyd (Eds.), 
    VCH Publishers, New York (NY), pp. 191-253, 1991].
    --->MEcc
    #################################################################
    """
    InertiaMatrix=GetInertiaMatrix(mol,ChargeCoordinates)
    u,s,v=scipy.linalg.svd(InertiaMatrix)
    
    res1=s[0]
    res3=s[2]
    
    res=math.pow(res1*res1-res3*res3,1./2)/res1
    return round(res,3)
    



#############################################################################    

def GetGeometric(mol):
    """
    #################################################################
    Wrapper for Geometrical descriptors
    #################################################################
    """
    filename='temp'
    ChargeCoordinates=_ReadCoordinates(filename)
    res={}
    res['W3DH']=Calculate3DWienerWithH(ChargeCoordinates)
    res['W3D']=Calculate3DWienerWithoutH(ChargeCoordinates)
    res['Petitj3D']=CalculatePetitjean3DIndex(ChargeCoordinates)
    res['GeDi']=CalculateGemetricalDiameter(ChargeCoordinates)
    res['grav']=CalculateGravitational3D1(mol,ChargeCoordinates)
    res['rygr']=CalculateRadiusofGyration(mol,ChargeCoordinates)
    res['Harary3D']=CalculateHarary3D(ChargeCoordinates)
    res['AGDD']=CalculateAverageGeometricalDistanceDegree(ChargeCoordinates)
    res['SEig']=CalculateAbsEigenvalueSumOnGeometricMatrix(ChargeCoordinates)
    res['SPAN']=CalculateSPANR(mol,ChargeCoordinates)
    res['ASPAN']=CalculateAverageSPANR(mol,ChargeCoordinates)
    res['MEcc']=CalculateMolecularEccentricity(mol,ChargeCoordinates)
    #res.update(CalculatePrincipalMomentofInertia(mol,ChargeCoordinates))
    #res.update(CalculateRatioPMI(mol,ChargeCoordinates))
    
    return res

def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc
    pydoc.writedoc('geometric')
    
#############################################################################
if __name__=="__main__":
    
    
    from GeoOpt import GetARCFile
    mol='C1C=CCCS1'
    mol='ClC(Cl)(Cl)Cl'


    inputmol=pybel.readstring('smi',mol)  
    GetARCFile(inputmol)
    result=GetGeometric(inputmol)
    print result
    print len(result)
