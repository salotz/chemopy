# -*- coding: utf-8 -*-
"""
##############################################################################

The calculation of whole holistic invariant molecular descriptors (WHIM). 

You can get 70 molecular decriptors. You can freely use and distribute it.

If you hava any problem, you could contact with us timely!

Authors: Dongsheng Cao and Yizeng Liang, Qingsong Xu

Date: 2011.04.19

Email: oriental-cds@163.com

##############################################################################
"""
import pybel
import scipy
import scipy.linalg
from GeoOpt import _ReadCoordinates
from AtomProperty import GetRelativeAtomicProperty



Version=1.0
#############################################################################
_filename='temp'

def GetAtomCoordinateMatrix():
    
    """
    #################################################################
    Get the atom coordinate matrix
    #################################################################
    """
    ChargeCoordinates=_ReadCoordinates(_filename)
    nAtom=len(ChargeCoordinates)
    CoordinateMatrix=scipy.zeros([nAtom,3])
    
    AtomLabel=[]
    
    for i,j in enumerate(ChargeCoordinates):

        CoordinateMatrix[i,:]=[j[1],j[2],j[3]]
        AtomLabel.append(j[0])
    
    return scipy.matrix(CoordinateMatrix),AtomLabel


def XPreCenter(X):
    """
    #################################################################
    Center the data matrix X
    #################################################################
    """
    Xdim=scipy.size(X,axis=0)
    Xmean=scipy.mean(X,axis=0)
    Xmean=scipy.matrix(Xmean)
    Xp=X-scipy.ones([Xdim,1])*Xmean
   
    return Xp


def GetPropertyMatrix(AtomLabel,proname='m'):
    """
    #################################################################
    #################################################################
    """
    res=[]
    for i in AtomLabel:
        res.append(GetRelativeAtomicProperty(i,proname))
    
    return scipy.matrix(scipy.diag(res))
        


def GetSVDEig(CoordinateMatrix,AtomLabel,proname='u'):
    """
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    
    S=XPreCenter(CoordinateMatrix)

    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
    
    return s




def GetWHIM1(CoordinateMatrix,AtomLabel,proname='u'):
    """
    #################################################################
    WHIM descriptors
    --->L1u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[0],3)


def GetWHIM2(CoordinateMatrix,AtomLabel,proname='u'):
    """
    #################################################################
    WHIM descriptors
    
    --->L2u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[1],3)


def GetWHIM3(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->L3u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[2],3)


def GetWHIM4(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Tu
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    T=round(sum(s),3)
    
    return T

def GetWHIM5(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Au
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    A=s[0]*s[1]+s[0]*s[2]+s[1]*s[2]
    
    return round(A,3)


def GetWHIM6(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Vu
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    A=s[0]*s[1]+s[0]*s[2]+s[1]*s[2]
    T=sum(s)
    V=A+T+s[0]*s[1]*s[2]
    
    return round(V,3)

def GetWHIM7(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->P1u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[0]/(s[0]+s[1]+s[2]),3)


def GetWHIM8(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->P2u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[1]/(s[0]+s[1]+s[2]),3)




def GetWHIM9(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Ku
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    res=0.0
    for i in s:
        res=res+abs(i/sum(s)-1/3.0)
        
    Ku=3.0/4*res
    
    return round(Ku,3) 


def GetWHIM10(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E1u
    #################################################################
    """
    
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    res=scipy.power(s[0],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,0]).T,4))
    
    return round(float(res.real),3)
    
    
def GetWHIM11(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E2u
    #################################################################
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    res=scipy.power(s[1],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,1]).T,4))
    
    return round(float(res.real),3)
    

def GetWHIM12(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E3u
    #################################################################
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    res=scipy.power(s[2],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,2]).T,4))
    
    return round(float(res.real),3)
    
def GetWHIM13(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Du
    #################################################################
    """
    c1=GetWHIM10(CoordinateMatrix,AtomLabel,proname)
    c2=GetWHIM11(CoordinateMatrix,AtomLabel,proname)
    c3=GetWHIM12(CoordinateMatrix,AtomLabel,proname)
    Du=c1+c2+c3
    
    return round(float(Du),3)




def GetWHIM14(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->P3u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[2]/(s[0]+s[1]+s[2]),3)


###########################################################################
def GetWHIMUnweighted():
    
    """
    #################################################################
    Wrapper for the unweighted WHIM descriptors. 
    #################################################################
    """
    
    res={}
    CoordinateMatrix,AtomLabel=GetAtomCoordinateMatrix()  
    res['L1u']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='u')
    res['L2u']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='u')
    res['L3u']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='u')
    res['Tu']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='u')
    res['Au']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='u')
    res['Vu']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='u')
    res['P1u']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='u')
    res['P2u']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='u')
    res['Ku']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='u')
    res['E1u']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='u')
    res['E2u']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='u')
    res['E3u']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='u')
    res['Du']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='u')
    res['P3u']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='u')
    
    return res


def GetWHIMMass():
    
    """
    #################################################################
    Wrapper for the WHIM descriptors based on atomic mass. 
    #################################################################
    """
    
    res={}
    CoordinateMatrix,AtomLabel=GetAtomCoordinateMatrix()  
    res['L1m']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='m')
    res['L2m']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='m')
    res['L3m']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='m')
    res['Tm']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='m')
    res['Am']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='m')
    res['Vm']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='m')
    res['P1m']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='m')
    res['P2m']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='m')
    res['Km']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='m')
    res['E1m']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='m')
    res['E2m']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='m')
    res['E3m']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='m')
    res['Dm']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='m')
    res['P3m']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='m')
    
    return res


def GetWHIMSandersonElectronegativity():
    
    """
    #################################################################
    Wrapper for the WHIM descriptors based on Sanderson Electronegativity. 
    #################################################################
    """
    
    res={}
    CoordinateMatrix,AtomLabel=GetAtomCoordinateMatrix()  

    res['L1e']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='En')
    res['L2e']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='En')
    res['L3e']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='En')
    res['Te']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='En')
    res['Ae']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='En')
    res['Ve']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='En')
    res['P1e']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='En')
    res['P2e']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='En')
    res['Ke']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='En')
    res['E1e']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='En')
    res['E2e']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='En')
    res['E3e']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='En')
    res['De']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='En')
    res['P3e']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='En')
    
    return res
    
    

def GetWHIMVDWVolume():
    
    """
    #################################################################
    Wrapper for the WHIM descriptors based on VDW Volume. 
    #################################################################
    """
    
    res={}
    CoordinateMatrix,AtomLabel=GetAtomCoordinateMatrix()  

    res['L1v']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='V')
    res['L2v']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='V')
    res['L3v']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='V')
    res['Tv']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='V')
    res['Av']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='V')
    res['Vv']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='V')
    res['P1v']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='V')
    res['P2v']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='V')
    res['Kv']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='V')
    res['E1v']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='V')
    res['E2v']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='V')
    res['E3v']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='V')
    res['Dv']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='V')
    res['P3v']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='V')
    
    return res



def GetWHIMPolarizability():
    
    """
    #################################################################
    Wrapper for the WHIM descriptors based on Polarizability. 
    #################################################################
    """
    
    res={}
    CoordinateMatrix,AtomLabel=GetAtomCoordinateMatrix()  
    res['L1p']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='alapha')
    res['L2p']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='alapha')
    res['L3p']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Tp']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Ap']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Vp']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='alapha')
    res['P1p']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='alapha')
    res['P2p']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Kp']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='alapha')
    res['E1p']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='alapha')
    res['E2p']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='alapha')
    res['E3p']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Dp']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='alapha')
    res['P3p']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='alapha')
    
    return res



def GetWHIM():
    
    """
    #################################################################
    Wrapper for the WHIM descriptors. 
    #################################################################
    """
    
    res={}
    CoordinateMatrix,AtomLabel=GetAtomCoordinateMatrix()  
    res['L1u']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='u')
    res['L2u']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='u')
    res['L3u']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='u')
    res['Tu']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='u')
    res['Au']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='u')
    res['Vu']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='u')
    res['P1u']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='u')
    res['P2u']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='u')
    res['Ku']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='u')
    res['E1u']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='u')
    res['E2u']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='u')
    res['E3u']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='u')
    res['Du']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='u')
    res['L1m']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='m')
    res['L2m']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='m')
    res['L3m']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='m')
    res['Tm']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='m')
    res['Am']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='m')
    res['Vm']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='m')
    res['P1m']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='m')
    res['P2m']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='m')
    res['Km']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='m')
    res['E1m']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='m')
    res['E2m']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='m')
    res['E3m']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='m')
    res['Dm']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='m')
    res['L1e']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='En')
    res['L2e']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='En')
    res['L3e']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='En')
    res['Te']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='En')
    res['Ae']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='En')
    res['Ve']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='En')
    res['P1e']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='En')
    res['P2e']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='En')
    res['Ke']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='En')
    res['E1e']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='En')
    res['E2e']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='En')
    res['E3e']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='En')
    res['De']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='En')
    res['L1v']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='V')
    res['L2v']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='V')
    res['L3v']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='V')
    res['Tv']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='V')
    res['Av']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='V')
    res['Vv']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='V')
    res['P1v']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='V')
    res['P2v']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='V')
    res['Kv']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='V')
    res['E1v']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='V')
    res['E2v']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='V')
    res['E3v']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='V')
    res['Dv']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='V')
    res['L1p']=GetWHIM1(CoordinateMatrix,AtomLabel,proname='alapha')
    res['L2p']=GetWHIM2(CoordinateMatrix,AtomLabel,proname='alapha')
    res['L3p']=GetWHIM3(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Tp']=GetWHIM4(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Ap']=GetWHIM5(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Vp']=GetWHIM6(CoordinateMatrix,AtomLabel,proname='alapha')
    res['P1p']=GetWHIM7(CoordinateMatrix,AtomLabel,proname='alapha')
    res['P2p']=GetWHIM8(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Kp']=GetWHIM9(CoordinateMatrix,AtomLabel,proname='alapha')
    res['E1p']=GetWHIM10(CoordinateMatrix,AtomLabel,proname='alapha')
    res['E2p']=GetWHIM11(CoordinateMatrix,AtomLabel,proname='alapha')
    res['E3p']=GetWHIM12(CoordinateMatrix,AtomLabel,proname='alapha')
    res['Dp']=GetWHIM13(CoordinateMatrix,AtomLabel,proname='alapha')
    res['P3p']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='alapha')
    res['P3u']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='u')
    res['P3m']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='m')
    res['P3e']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='En')
    res['P3v']=GetWHIM14(CoordinateMatrix,AtomLabel,proname='V')
    
    return res



def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc
    pydoc.writedoc('whim')    
#############################################################################
if __name__=="__main__":

    from GeoOpt import GetARCFile
    mol='c1ccccc1N'
    inputmol=pybel.readstring('smi',mol)  
    GetARCFile(inputmol)
    result=GetWHIMSandersonElectronegativity()
    print result
    print len(result)