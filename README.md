# FLASH-AVM
A Matlab-based analysis and visualization toolkit for the FLASH simulation code
(http://flash.uchicago.edu/site/flashcode/). While many visualization software
are compatible with FLASH data, few are well-suited for making
publication-quality figures. The intent of FLASH-AVM is to provide a
Matlab-based alternative to the yt python package
(https://yt-project.org/docs/dev/index.html), capable of producing high quality
visualizations of volumetric data. The toolkit allows the user to read from and
write to multi-dimensional, multi-block, adaptive mesh refinement (AMR) FLASH
HDF5-formatted checkpoint or plot files, and to design their own custom
analysis and visualization scripts. Reader files allow the user to read in AMR
data and plot or post-process the data. Capabilities include: prolonging coarse
data to the finest available level, plotting mesh block outlines, reading from
integral quantities files, and more. One example is provided of the
two-dimensional Sedov problem with multiple levels of adaptive mesh refinement
active. The code is currently in its infancy, and collaboration with interested
parties is welcome.

## Quick-start:

Set the following environment variable:
FLASHAVM="path-to-flashavm/Release/src"

Add the following command to any Matlab script:
addpath(genpath(getenv('FLASHAVM')));
