# LWimager-SASIR
# This version of LWimager is used by the docker image at SASIR-docker (github.com/JulienNGirard/SASIR-docker

List of patched files:

U = updated file
A = added file

Requirements for the compressed sensing mode to work (check out github.com/JulienNGirard/SASIR-library for building):
* Header :  CEA_comp_sens.h  (put in /usr/include)
* CEA library:  libCEA_comp_sens.a  (put in /usr/lib)

In order to insert the compressed sensing mode as a new "operation" of lwimager, some files need to be updated/added (which has already be done in this repository)

LWimager file alteration:

1) Cmake files:

U CMakeLists.txt		
U synthesis/CMakeLists.txt
A cmake/FindCEACompSens.cmake

2) LWimager:

U- synthesis/apps/lwimager.cc  # adding "compsens" input parameters & call

A- synthesis/MeasurementComponents/CoSeImageSkyModel.h  # CS method header
A- synthesis/MeasurementComponents/CoSeImageSkyModel.cc # CS method source

U- synthesis/MeasurementEquations/Imager.h # declaring CS method
U- synthesis/MeasurementEquations/Imager.cc # adding "compsens"
U- synthesis/MeasurementEquations/Imager2.cc # to implement the final convolution option of the image

U- synthesis/MeasurementEquations/ClarkCleanLatModel.h  # (getPSFpatch made public)


Questions:
julien.girard@cea.fr or 4m0nr3@gmail.com


