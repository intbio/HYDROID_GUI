# HYDROID_GUI
GUI for HYDROID - Python package for analyzing hydroxyl-radical footprinting experiments of protein-DNA complexes (https://github.com/ncbi/HYDROID)
HYROID_GUI is a sister package that wraps some basic gel lane quantification functionality into a more user friendly graphical interface.

## Video tutorial:
[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/dJVoKrpH4f4/0.jpg)](http://www.youtube.com/watch?v=dJVoKrpH4f4)

## Installation
HYDROID_GUI gui is tested to work on python 2.7 through conda environment manager, it runs on Linux (ubuntu 16.04, 14.04 tested), Windows and MacOS
### install hydroid
- conda install -c hydroid hydroid
- conda install pyqt
- conda install -c conda-forge billiard
### clone this repo
- git clone https://github.com/intbio/HYDROID_GUI
- cd hydroid gui
### run 
- python hydroid_gui
### usage
HYROID_GUI wraps basic gel lane quantification functionality of HYDROID into a more user friendly graphical interface.

Hydroid gui allows user to 
- open lane profile files [extracted from gel images via ImageJ](https://github.com/ncbi/HYDROID/blob/master/examples/example1/exp_s1_extract_lp.md) [video](https://youtu.be/7UCb0IkXL2g)
- perform peak detection, naming and quantification for multiple lanes simultaneously, [see HYDROID workflow for details](https://github.com/ncbi/HYDROID/blob/master/docs/INDEX.md)
- export results directly to a spreadsheet editor or scv files
- save config files for later use with HYDROID Python scripting interface.

