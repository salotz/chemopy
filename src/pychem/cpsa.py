# -*- coding: utf-8 -*-
"""
##############################################################################
This module is used for calculating charged partial surface area descriptors

(CPSA) proposed by Peter. Jurs. If you have any question please contact me

via email.

Created on Tue Apr 19 10:40:45 2011

@author: Dongsheng Cao, Yizeng Liang, Qingsong Xu
##############################################################################
"""

import asa
import pybel
from GeoOpt import GetAtomClassList, _ReadCoordinates


Version=1.0
############################################################################
filename='temp'

def GetChargeSA(RadiusProbe=1.5,n_sphere_point=960):
    """
    #################################################################
    Get the list form for all atoms in a molecule.
    
    It includes the atom symbol, charge and partial solvent-accessible 
    
    surface areas. 
    
    Note that this is list form whose element is still list form of each atom.
    #################################################################
    """
    ChargeCoordinates=_ReadCoordinates(filename)
    atoms=GetAtomClassList(ChargeCoordinates)
    FASA=asa.calculate_asa(atoms, RadiusProbe, n_sphere_point)
    
    res=[]
    for i in range(len(FASA)):
        res.append([ChargeCoordinates[i][0],ChargeCoordinates[i][4],FASA[i]])
        
    return res
        

def CalculateASA(ChargeSA):
    """
    #################################################################
    The calculation of solvent-accessible surface areas
    
    -->ASA
    #################################################################
    """
    res=0.0
    for i in ChargeSA:
        res=res+i[2]
    return res
    
def CalculateMSA():
    """
    #################################################################
    The calculation of molecular surface areas
    
    -->MSA
    #################################################################
    """
    ChargeSA=GetChargeSA(RadiusProbe=0,n_sphere_point=960)
    res=0.0
    for i in ChargeSA:
        res=res+i[2]
    return res
    
    
def CalculatePNSA1(ChargeSA):
    """
    #################################################################
    The calculation of partial negative area
    
    It is the sum of the solvent-accessible surface areas of all 
    
    negatively charged atoms.
    
    -->PNSA1
    #################################################################
    """
    res=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            res=res+i[2]
            
    return res
            

def CalculatePPSA1(ChargeSA):
    """
    #################################################################
    The calculation of partial negative area
    
    It is the sum of the solvent-accessible surface areas of 
    
    all positively charged atoms.
    
    -->PPSA1
    #################################################################
    """
    res=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            res=res+i[2]
            
    return res

def CalculatePNSA2(ChargeSA):
    """
    #################################################################
    The calculation of total charge wighted negative surface area
    
    It is the partial negative solvent-accessible surface area 
    
    multiplied by the total negative charge.
    
    -->PNSA2
    #################################################################
    """
    temp1=temp2=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            temp1=temp1+float(i[1])
            temp2=temp2+i[2]
    res=temp1*temp2
    
    return res
            


def CalculatePPSA2(ChargeSA):
    """
    #################################################################
    The calculation of total charge wighted negative surface area
    
    It is the partial negative solvent-accessible surface area 
    
    multiplied by the total positive charge.
    
    -->PPSA2
    #################################################################
    """
    temp1=temp2=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            temp1=temp1+float(i[1])
            temp2=temp2+i[2]
    res=temp1*temp2
    
    return res

def CalculatePNSA3(ChargeSA):
    """
    #################################################################
    The calculation of atom charge weighted negative surface ares
    
    It is the sum of the products of atomic solvent-accessible 
    
    surface area and partial charges over all negatively charges atoms.
    
    -->PNSA3
    #################################################################
    """
    res=0.0
    for i in ChargeSA:
        if float(i[1])<0:
            res=res+float(i[1])*i[2]
    return res


def CalculatePPSA3(ChargeSA):
    """
    #################################################################
    The calculation of atom charge weighted positive surface ares
    
    It is the sum of the products of atomic solvent-accessible
    
    surface area and partial charges over all positively charges atoms.
    
    -->PPSA3
    #################################################################
    """
    res=0.0
    for i in ChargeSA:
        if float(i[1])>0:
            res=res+float(i[1])*i[2]
    return res


def CalculateDPSA1(ChargeSA):
    """
    #################################################################
    The calculation of difference in charged partial surface area
    -->DPSA1
    #################################################################
    """
    return CalculatePPSA1(ChargeSA)-CalculatePNSA1(ChargeSA)


def CalculateDPSA2(ChargeSA):
    """
    #################################################################
    The calculation of difference in total charge weighted partial
    
    surface area
    
    -->DPSA2
    #################################################################
    """
    return CalculatePPSA2(ChargeSA)-CalculatePNSA2(ChargeSA)

def CalculateDPSA3(ChargeSA):
    """
    #################################################################
    The calculation of difference in atomic charge weighted surface area
    
    -->DPSA3
    #################################################################
    """
    return CalculatePPSA3(ChargeSA)-CalculatePNSA3(ChargeSA)

def CalculateFNSA1(ChargeSA):
    """
    #################################################################
    The calculation of fractional charged partial negative surface areas
    
    -->FNSA1
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePNSA1(ChargeSA)/temp

def CalculateFNSA2(ChargeSA):
    """
    #################################################################
    The calculation of fractional charged partial negative surface areas
    
    -->FNSA2
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePNSA2(ChargeSA)/temp


def CalculateFNSA3(ChargeSA):
    """
    #################################################################
    The calculation of fractional charged partial negative surface areas
    
    -->FNSA3
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePNSA3(ChargeSA)/temp

def CalculateFPSA1(ChargeSA):
    """
    #################################################################
    The calculation of fractional charged partial negative surface areas
    
    -->FPSA1
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePPSA1(ChargeSA)/temp

def CalculateFPSA2(ChargeSA):
    """
    #################################################################
    The calculation of fractional charged partial negative surface areas
    
    -->FPSA2
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePPSA2(ChargeSA)/temp

def CalculateFPSA3(ChargeSA):
    """
    #################################################################
    The calculation of fractional charged partial negative surface
    
    areas
    
    -->FPSA3
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePPSA3(ChargeSA)/temp

def CalculateWNSA1(ChargeSA):
    """
    #################################################################
    The calculation of surface weighted charged partial negative 
    
    surface areas
    
    -->WNSA1
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePNSA1(ChargeSA)*temp/1000

def CalculateWNSA2(ChargeSA):
    """
    #################################################################
    The calculation of surface weighted charged partial negative 
    
    surface areas
    
    -->WNSA2
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePNSA2(ChargeSA)*temp/1000

def CalculateWNSA3(ChargeSA):
    """
    #################################################################
    The calculation of surface weighted charged partial negative 
    
    surface areas
    
    -->WNSA3
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePNSA3(ChargeSA)*temp/1000

def CalculateWPSA1(ChargeSA):
    """
    #################################################################
    The calculation of surface weighted charged partial negative 
    
    surface areas
    
    -->WPSA1
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePPSA1(ChargeSA)*temp/1000

def CalculateWPSA2(ChargeSA):
    """
    #################################################################
    The calculation of surface weighted charged partial negative 
    
    surface areas
    
    -->WPSA2
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePPSA2(ChargeSA)*temp/1000

def CalculateWPSA3(ChargeSA):
    """
    #################################################################
    The calculation of surface weighted charged partial negative
    
    surface areas
    
    -->WPSA3
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
        
    return CalculatePPSA3(ChargeSA)*temp/1000


def CalculateTASA(ChargeSA):
    """
    #################################################################
    The calculation of total hydrophobic surface area
    
    -->TASA
    #################################################################
    """
    res=0.0
    for i in ChargeSA:
        if abs(float(i[1]))<0.2:
            res=res+i[2]
    return res

def CalculateTPSA(ChargeSA):
    """
    #################################################################
    The calculation of total polar surface area
    
    -->PSA
    #################################################################
    """
    res=0.0
    for i in ChargeSA:
        if abs(float(i[1]))>=0.2:
            res=res+i[2]
    return res
    
    
def CalculateFractionTATP(ChargeSA):
    """
    #################################################################
    The fraction between TASA and TPSA
    
    --->FrTATP
    #################################################################
    """
    res=0.0
    if CalculateTPSA(ChargeSA)==0:
        return res
    else:
        return CalculateTASA(ChargeSA)/CalculateTPSA(ChargeSA)

def CalculateRASA(ChargeSA):
    """
    #################################################################
    The calculation of relative hydrophobic surface area
    
    -->RASA
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    return CalculateTASA(ChargeSA)/temp

def CalculateRPSA(ChargeSA):
    """
    #################################################################
    The calculation of relative polar surface area
    
    -->RPSA
    #################################################################
    """
    temp=0.0
    for i in ChargeSA:
        temp=temp+i[2]
    return CalculateTPSA(ChargeSA)/temp


def CalculateRNCS(ChargeSA):
    """
    #################################################################
    The calculation of relative negative charge surface area
    
    -->RNCS
    #################################################################
    """
    charge=[]
    for i in ChargeSA:
        charge.append(float(i[1]))
    
    temp=[]
    for i in ChargeSA:
        temp.append(i[2])
     
    RNCG=min(charge)/sum([i for i in charge if i<0])
    
    
    return  temp[charge.index(min(charge))]/RNCG
    
          

def CalculateRPCS(ChargeSA):
    """
    #################################################################
    The calculation of relative positive charge surface area
    
    -->RPCS
    #################################################################
    """
    charge=[]
    for i in ChargeSA:
        charge.append(float(i[1]))
    
    temp=[]
    for i in ChargeSA:
        temp.append(i[2])
     
    RPCG=max(charge)/sum([i for i in charge if i>0])
      
    return  temp[charge.index(min(charge))]/RPCG
    
    
############################################################################
def GetCPSA():
    """
    #################################################################
    Wrapper for the CPSA descriptors
    #################################################################
    """
    ChargeSA=GetChargeSA(RadiusProbe=1.5,n_sphere_point=5000)
    
    res={}
    res['ASA']=round(CalculateASA(ChargeSA),3)
    res['MSA']=round(CalculateMSA(),3)
    res['PNSA1']=round(CalculatePNSA1(ChargeSA),3)
    res['PNSA2']=round(CalculatePNSA2(ChargeSA),3)
    res['PNSA3']=round(CalculatePNSA3(ChargeSA),3)
    res['PPSA1']=round(CalculatePPSA1(ChargeSA),3)
    res['PPSA2']=round(CalculatePPSA2(ChargeSA),3)
    res['PPSA3']=round(CalculatePPSA3(ChargeSA),3)
    res['DPSA1']=round(CalculateDPSA1(ChargeSA),3)
    res['DPSA2']=round(CalculateDPSA2(ChargeSA),3)
    res['DPSA3']=round(CalculateDPSA3(ChargeSA),3)
    res['FNSA1']=round(CalculateFNSA1(ChargeSA),3)
    res['FNSA2']=round(CalculateFNSA2(ChargeSA),3)
    res['FNSA3']=round(CalculateFNSA3(ChargeSA),3)
    res['FPSA1']=round(CalculateFPSA1(ChargeSA),3)
    res['FPSA2']=round(CalculateFPSA2(ChargeSA),3)
    res['FPSA3']=round(CalculateFPSA3(ChargeSA),3)
    res['WNSA1']=round(CalculateWNSA1(ChargeSA),3)
    res['WNSA2']=round(CalculateWNSA2(ChargeSA),3)
    res['WNSA3']=round(CalculateWNSA3(ChargeSA),3)
    res['WPSA1']=round(CalculateWPSA1(ChargeSA),3)
    res['WPSA2']=round(CalculateWPSA2(ChargeSA),3)
    res['WPSA3']=round(CalculateWPSA3(ChargeSA),3)
    res['TASA']=round(CalculateTASA(ChargeSA),3)
    res['PSA']=round(CalculateTPSA(ChargeSA),3)
    res['RASA']=round(CalculateRASA(ChargeSA),3)
    res['RPSA']=round(CalculateRPSA(ChargeSA),3)
    res['RNCS']=round(CalculateRNCS(ChargeSA),3)
    res['RPCS']=round(CalculateRPCS(ChargeSA),3)
    res['FrTATP']=round(CalculateFractionTATP(ChargeSA),3)


    return res


def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc
    pydoc.writedoc('cpsa')
#############################################################################    
if __name__=="__main__":
    
    from GeoOpt import GetARCFile
    mol='C1C=CCS1'

    inputmol=pybel.readstring('smi',mol)  
    GetARCFile(inputmol)
    result=GetCPSA()
    print result
    print len(result)
