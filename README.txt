Copyright (c) 2012CBDD03, Computational Biology and Drug Design Group, 
Central South University, China.
All rights reserved.

INSTRUCTION

The molecular descriptors from molecular structures have been widely
applied in the QSAR/QSPR studies, database searching, clustering, ranking,
drug ADME/T predcition, chemoinformatics, and current chemogeonomics. 
These descriptors capture and magnify distinct aspects of molecular structures.
The usefulness of these descriptors covered by PyChem for representing structure 
features of small molecules have been sufficiently demonstrated by a number of 
published studies of the development of machine learning classification systems
in chemoinformatics fields. We developed a powerful python package to calculate 
a large number of molecular descriptors in the field of chemoinforamtics.

If you have any problem, contact Dongsheng Cao (oriental-cds@163.com)
#######################FEATURES###########################
The drug descriptors calculated by PyChem:

(1) Constitutional descriptors (30)
(2) Topologcial descriptors (35)
(3) Molecular connectivity descriptors (44)
(4) Kappa descriptors(7)
(5) E-state descriptors (245)
(6) Autocorrelation descriptors (96) 
including Moreau-Broto,Moran and Geary descriptors
(7) Charge descriptors (25)
(8) Molecular property descriptors (6)
(9) MOE-type descriptors (60)
(10) Burden descriptors (64)
(11) Basak information descriptors (21)
(12) 3-D MoRSE descriptors (210)
(13) 3-D RDF descriptors (180)
(14) Geometric descriptors (12)
(15) WHIM descriptors (70)
(16) CPSA descriptors (30)
(17) Daylight fingerprint (2048)
(18) MACCS keys (166)
(19) FP4 fingerprints (307) 
(20) E-state fingerprints (79)
(21) Atom Paris fingerprints and Morgan fingerprints

For detailed information on these descriptors, please refer to 
the PyChem manual.
##########################################################
INSTALL

PyChem has been successfully tested on Linux and Windows systems.
The author could download the PyChem package from
http://code.google.com/p/pychem/downloads/list (.zip and .tar.gz). 
The install process of PyChem is very easy:
***************************************************************************
*You first need to install RDkit, Openbabel, MOPAC and pybel successfully.*
***************************************************************************
Openbabel and pybel can be downloaded via: http://openbabel.org/wiki/Main_Page
RDkit can be downloaded via: http://code.google.com/p/rdkit/
MOPAC can be downloaded via: http://openmopac.net/

Note: pychem was tested in MOPAC 7.

On Windows:
(1): download the pychem package (.zip)
(2): extract or uncompress the .zip file 
(3): cd pychem-1.0
(4): python setup.py install

On Linux:
(1): download the pychem package (.tar.gz) 
(2): tar -zxf pychem-1.0.tar.gz 
(3): cd pychem-1.0 
(4): python setup.py install or sudo python setup.py install 
############################################################
EXAMPLE:
For more examples, please see the user guide.
############################################################

