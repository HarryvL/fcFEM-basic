# fcFEM-basic
Finite Element macro for [FreeCAD](https://freecad.org)

<img src="https://github.com/HarryvL/fcFEM-basic/blob/main/Pictures/Slope_Clay_Fine_PV.png" height="200"/> <img src="https://github.com/HarryvL/fcFEM-basic/blob/main/Pictures/Slope_Clay_Fine.png" height="220" raw=true/> 

### Background
fcFEM is a finite element solver for performing collapse analysis of structures and soil bodies. It is based on the theory of elastoplasticity and gives insight in ductility and reserve strength beyond first yield. The theory underpinning fcFEM can be found in the document included in this repository.

<img src="https://github.com/HarryvL/fcFEM-basic/blob/main/Pictures/punch_indentation_VTK.png" height="200"/>

### Prerequisites
* FreeCAD >= v0.20

### Installation
1. Copy the following files `fcFEM.FCMacro`, `femTools.py` and the `"name".inp`files in a single directory on your machine, where "name" comes from the "name".FCStd` file of the freeCAD problem.  
1. Within the [FreeCAD macro dialog](https://wiki.freecad.org/Macros) set the macro path to the above directory.  
1. Run `fcFEM.FCMacro` with the FreeCAD macro editor.

<img src="https://github.com/HarryvL/fcFEM-basic/blob/main/Pictures/Plate_with_Hole_3.png" height="200"/>

### Dependencies
fcFEM imports (from) the following packages: numpy, scipy, numba, matplotlib and scikit-sparse. On Linux and Mac-OS systems these can be installed with the usual package managers. On Windows this requires more effort:
1. download Miniforge: https://github.com/conda-forge/miniforg ... x86_64.exe
1. run the installer: Miniforge-pypy3-Windows-x86_64.exe
1. find and run Miniforge on your system - this opens a Miniforge Prompt: (base) C:\"path">
1. create a new virtual environment: (base) C:\"path"> mamba create --name fcFEM (or any other name of your choice)
1. change into the new environment: (base) C:\"path"> mamba activate fcFEM (or the other name you chose)
1. instal freeCAD and dependencies: (fcFEM) C:\"path"> mamba install freeCAD scipy numba matplotlib scikit-sparse (with spaces and no commas)
1. check with python if the dependencies can be imported: (fcFEM) C:\"path"> python
1. >>> import scipy.sparse
1. >>> import sksparse.cholmod
1. if no problems quit python and start freecad: (fcFEM) C:\"path"> freecad

### Documentation
Please refer to source code for in-line comments and to the FreeCAD forum (https://forum.freecadweb.org/viewforum.php?f=18)

### TODO

- [ ] Frictional material for the analysis of soils and concrete.
- [ ] Tension cut-off and reinforcement capabilities.
- [ ] Large deformation analysis capability.
- [ ] Addition of beam and shell elements.
- [ ] Linear buckling and initial imperfections for non-linear buckling.
- [ ] Loading stages.
- [ ] Advanced material modelling.

<img src="https://github.com/HarryvL/fcFEM-basic/blob/main/Pictures/Beam_Collapse%20(2).jpg" width="450"/>

### Licence information

Copyright (c) 2019 - Harry van Langen <hvlanalysis@icloud.com>  


This program is free software; you can redistribute it and/or modify  
it under the terms of the GNU Lesser General Public License (LGPL)    
as published by the Free Software Foundation; either version 2 of     
the License, or (at your option) any later version.                   
for detail see the LICENCE text file.                                 
                                                                         
This program is distributed in the hope that it will be useful,       
but WITHOUT ANY WARRANTY; without even the implied warranty of        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
GNU Library General Public License for more details.                  
                                                                         
You should have received a copy of the GNU Library General Public     
License along with this program; if not, write to the Free Software   
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
USA                                                                   
