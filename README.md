# DRAGON-G4-LaBr3
Development of GEANT4 code for the LaBr3 Gamma Ray Detector
# Overview
This is the development code for a LaBr3 gamma ray detector for DRAGON at TRIUMF. The intention of the simulation development is to test various properties and characteristics of LaBr3 to be used in a possible future replacement for the DRAGON BGO gamma array.The code is being built to be generalized for various geometries of LaBr3 crystal that the user can set at run-time. The code will include the internal radiation present in all LaBr3 scintillation detectors. The GPS module of G4 is being used to expose the crystal to various energies of gamma rays. The source characteristics can be controlled at run-time. 
# How to Compile
The user should create a GEANT4 working directory such as:
```
> mkdir G4WorkDir
```
In your work directory create a source directory (e.g. LaBr3_v1)
```
> mkdir LaBr3_v1
```
Pull the entire contents of this repository to that directory
```
> git pull origin master
```
Go back to your work directory and create a build directory (e.g LaBr3_v1-build)
```
> cd ..
> mdkir LaBr3_v1-build
````
Switch to your build directory and run cmake using two arguments, the -D option to point to the directory holding the Geant4Config.cmake file that Geant4 installs to help CMake find and use Geant4, and the second argument is the source directory for the application you want to build. For example, if you were using the version 10.3 of G4 and it was installed in the conventional way then you would do:
```
> cd LaBr3_v1-build
> cmake -DGeant4_DIR=/home/username/geant4-install/lib64/Geant4-10.3.0 /home/username/G4WorkDir/LaBr3_v1
```
Then from you build directory make/compile the code
```
> make
```
# How to Run the LaBr3 Simulation
The executable will be created in your build directory by default but you may change the location to somewhere else if you wish. The name of the executable can be changed by editting the CMakeLists.txt file and the filename of the executable. If it is unchanged the executable can be run by typing:

In Batch Mode:
```
> ./LaBr3_v4
```
In Interactive Mode:
```
> ./LaBr3_v4 vis.mac
```
