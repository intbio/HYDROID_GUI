## Building and uploading to Anaconda cloud
```
# all should be done in hydroid_gui_conda
conda install -y conda-build
conda install -y anaconda-client
# billiard is in conda-forge
conda-build -c hydroid -c conda-forge hydroid_gui

anaconda login
anaconda upload path_to_package
conda convert --platform all path_to_package -o output/

## Installing HYDROID from Anaconda cloud

```
conda install -c hydroid -c conda-forge hydroid_gui
```
