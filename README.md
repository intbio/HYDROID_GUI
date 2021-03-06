[![](https://anaconda.org/hydroid/hydroid_gui/badges/version.svg)](https://anaconda.org/hydroid/hydroid_gui)
[![](https://anaconda.org/hydroid/hydroid_gui/badges/platforms.svg)](https://anaconda.org/hydroid/hydroid_gui)
[![](https://anaconda.org/hydroid/hydroid_gui/badges/installer/conda.svg)](https://anaconda.org/hydroid/hydroid_gui)

# HYDROID_GUI
Graphical User Interface for [HYDROID](https://github.com/ncbi/HYDROID) - A Python package for analyzing hydroxyl-radical footprinting experiments of protein-DNA complexes.

HYROID_GUI is a sister add-on package that wraps some basic gel lane quantification functionality into a more user friendly graphical interface.

Please, refer to the HYDROID package web-page [https://github.com/ncbi/HYDROID](https://github.com/ncbi/HYDROID) for the detailed description of the HYDROID approach.

The manual for HYDROID_GUI is provided as a video tutorial (see below).

## Video tutorial:
[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/yQGMtqKzzbA/0.jpg)](https://youtu.be/yQGMtqKzzbA)

## Features
HYDROID_GUI allows user to 
- open lane profile files [extracted from gel images via ImageJ](https://github.com/ncbi/HYDROID/blob/master/examples/example1/exp_s1_extract_lp.md) ([see video tutorial](https://youtu.be/7UCb0IkXL2g))
- perform peak detection, naming and quantification for multiple lanes simultaneously, [see HYDROID workflow for details](https://github.com/ncbi/HYDROID/blob/master/docs/INDEX.md)
- export results directly to a spreadsheet editor or scv files
- save config files for later use with HYDROID Python scripting interface.

## Installation
HYDROID_GUI gui is tested to work on python 2.7 through conda environment manager, it runs on Linux (ubuntu 16.04, 14.04 tested), Windows and MacOS
### Install hydroid through conda
- conda install -c hydroid hydroid_gui -c conda-forge
### run 
- HYDROID_GUI


