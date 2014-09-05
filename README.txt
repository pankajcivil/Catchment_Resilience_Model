    Auxiliary Material Submission for Paper 2012WR013003
    Multiple hydrological attractors under stochastic daily forcing: 1. Can multiple attractors exist?
    Tim J. Peterson and Andrew William Western
    (Department of Infrastructure Engineering, The University of Melbourne, Australia)
    WRR., XXX (YYY), doi:10.1029/013003, 2012

INTRODUCTION
----------------------------------

The supplementary material accompanying this file comprise of MatLab code and data to build and run the hillslope Boussinesq model presented within the paper. Using this code, the model can be show to have two steady states for the same parameter set when under daily stochastic climate forcing.

The supplementary material is contained within the file "suplementary material.wrr". It is a zipped file with the extension changed from ".zip" to ".wrr". The extension change was undertaken because the AGU paper submission process unzipped the file to many hundreds of files when the ".zip" extension was used. To access the file, simply chnage the file extension to ".zip" and extract the contents.

CONTENTS OF SUPPLEMENTARY MATERIAL
----------------------------------

Below is a summary of the supplementary material:

* CatchmentResilienceModel.m	- MatLab class definition for creating the model.

* ContinuationAnalysis.m	- MatLab class definition for continuation analysis (i.e. for identifying attractors and repellors)

* Folder 'bsxops'		- Folder containing MatLab functions for efficient matrix manipulation. Copyright (c) 2009, Bruno Luong.

* Folder 'continuation'		- Folder containing the functions to undertake the continuation analysis. Te folder contains: (1) modified MatCont functions and (2) non-modified MatCont functions (see Dhoogle et al 2003) 

* ClimateData.mat		- MatLab data file containing the daily precipitation and ETo data used within Peterson and Western 2012. It also contains the climate replicate data investigated within Peterson, Western and Argent 2012.

USING THE SUPPLEMENTARY MATERIAL
---------------------------------

To begin building the vadose zone-Boussinesq model, please first read the model documentation. It can be accessed by opening MatLab and changing the current path (within MatLab) to the location of 'CatchmentResilienceModel.m'. Once there, enter the following command within the MatLab command window (ignore the quotation marks): "doc CatchmentResilienceModel". The documentation that should appear containing a step by step example of how to build and run the model. 

COPYRIGHT AND LICENSE
---------------------
The files listed below are Copyright (C) 2012  Dr. Tim J. Peterson, Prof. Andrew W. Western and Dr. Robert M. Argent. The files listed below are also licensed under the GNU General Public License 3. For the other files contained within this software that were not entirely written by Dr. Tim J. Peterson, Prof. Andrew W. Western and Dr. Robert M. Argent, their copyright and licensed are as detailed within their respective folders.

* CatchmentResilienceModel.m
* ContinuationAnalysis.m

The GNU General Public License 3 is detailed within "GNU GPL3.txt". A summary of the license is given below:

"This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details."


REFERENCES
----------
Dhooge, A., W. Govaerts, and Y. A. Kuznetsov (2003), MATCONT: a MATLABpackage for numerical bifurcation analysis of ODEs, ACM transactionson mathematical software, 29 (2), 141–164

