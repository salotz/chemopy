# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 00:27:37 2012

@author: orient
"""

from distutils.core import setup 


packagedata={'pychem': ['html/*','data/*','manual/*']}


setup(name = 'chemopy', 

	version = '1.0', 
	
	description ="A powerful tool for chemoinformatics study",
	
	author = "Dongsheng Cao",
	
	author_email = "oriental-cds@163.com",
	
	url ="http://cbdd.csu.edu.cn/index",
	
	license = "GPL",
	
	packages = ['pychem'],
	
	package_data=packagedata,
	
#	data_files = datafiles,
	
	package_dir={'pychem':'src/pychem'},
	
	scripts = [],
	
	py_modules = []

	)

