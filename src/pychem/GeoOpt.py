"""
This module is used for optimizing molecular structures and getting ARC file
 
based on pybel and MOPAC!

Written by Dongsheng Cao

Date:2011.3.23
"""


import pybel
import vector3d

import os
import string


Version=1.1
################################################################################
class Atom:
    """
    #################################################################
    A atom class used for wrapping some properties of atoms.
    
    Note that Coordinates is the output of the function 
    
    (_ReadCoordinates).
    #################################################################
    """

    def __init__(self,Coordinates):
        
        self.pos = vector3d.Vector3d()
        self.radius=0.0
        self.Coordinates=Coordinates
        self.Element=''
        
        
    def SetCoordinates(self):
        
        temp=self.Coordinates
        self.pos.x=float(temp[1])
        self.pos.y=float(temp[2])
        self.pos.z=float(temp[3])

    def GetCoordinates(self):
        
        self.SetCoordinates()
        
        return self.pos
        
    def SetElement(self):
        
        temp=self.Coordinates
        
        self.Element=temp[0]
        
    def GetElement(self):
        
        self.SetElement()
        
        return self.Element
        
    def SetRadius(self):
        
        radii={ 'H': 1.20,'N': 1.55,'Na': 2.27,'Cu': 1.40,'Cl': 1.75,'C': 1.70,
        'O': 1.52,'I': 1.98,'P': 1.80,'B': 1.85,'Br': 1.85,'S': 1.80,'Se': 1.90,
        'F': 1.47,'Fe': 1.80,'K': 2.75,'Mn': 1.73,'Mg': 1.73,'Zn': 1.39,'Hg': 1.8,
        'Li': 1.8,'.': 1.8}
        
        temp=self.GetElement()
        
        if temp in radii.keys():
            self.radius=radii[temp]
        else: 
            self.radius=radii['.']
            
    def GetRadius(self):
        
        self.SetRadius()
        
        return self.radius
###########################################################################

def GetAtomClassList(Coordinates):
    """
    #################################################################
    Combine all atoms in a molecule into a list form.
    
    Note that Coordinates is the output of the function (_ReadCoordinates).
    #################################################################
    """
    Atoms=[]
    for i in Coordinates:
        atom=Atom(i)
        atom.SetCoordinates()
        atom.SetElement()
        atom.SetRadius()
        Atoms.append(atom)
    return Atoms
    
###########################################################################    

def _ReadCoordinates(filename="temp"):
    """
    #################################################################
    Read the coordinates and charge of each atom in molecule from .arc file.
    #################################################################
    """
    res=[]
    
    f=file(filename,'r')
    templine=f.readlines()
    f.close()
    
    for line in range(len(templine)):
        if templine[line][-7:-1]=="CHARGE":
            k=line
            break
        
    for i in templine[k+4:len(templine)-1]:
        
        temp=i.split()
        ElementCoordinate=[string.strip(temp[0]),string.strip(temp[1]),
                           string.strip(temp[3]),string.strip(temp[5]),
                           string.strip(temp[10])]
        res.append(ElementCoordinate)

    return res

#############################################################################


def FormatConversion(inputmol):
    """
    #################################################################
    Using Pybel to convert the smi/sdf formats to mop format!
    #################################################################
    """
    #inputmol.removeh()
    inputmol.addh()
    inputmol.make3D(forcefield='mmff94',steps=50)    ##Gemetrical optimization
    ##forcefields = ['uff', 'mmff94', 'ghemical']
    #make3D(self, forcefield = "mmff94", steps = 50)
    ##inputmol.localopt(forcefield='mmff94',steps=50)
    outputmol=pybel.Outputfile('mop',"temp.dat",overwrite=True)
    outputmol.write(inputmol)
    outputmol.close()
    f=file('temp.dat','r+')
    f.write('AM1              ')
    f.close()


def RunMOPAC(filename):
    """
    #################################################################
    Run the MOPAC using os.system
    #################################################################
    """
    
    itest=os.system("run_mopac7"+" "+filename)
    #time.sleep(1)
    return itest

############################################################################ 
def GetARCFile(inputmol):  
    """
    #################################################################
    Get ARC file for each molecule
    #################################################################
    """
    
    FormatConversion(inputmol)
    
    itest=RunMOPAC('temp')

    if not itest:
        print itest,'\t', 'Finshed successfully!'
    else:
        print itest,'\t', 'Failed Finished........'
 
    os.remove('temp.dat')
    os.remove('temp.log')
    os.remove('temp.OUT')    
    #os.remove('temp.arc')
    oldpath=os.getcwd()+'/temp.arc'
    newpath=os.getcwd()+'/temp'
    os.rename(oldpath,newpath)
    
##############################################################################   
if __name__=="__main__":
    
    mol='C1C=CCS1'
    mol='SCCC(=O)N1[C@@H](CCC1)C(=O)OCC'
    inputmol=pybel.readstring('smi',mol)  
    GetARCFile(inputmol)
    res=_ReadCoordinates('temp')
    print res









