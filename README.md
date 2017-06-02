# DRAGON-G4-LaBr3
Development of GEANT4 code for the LaBr3 Gamma Ray Detector
# Overview
This is the development code for a LaBr3 gamma ray detector for DRAGON at TRIUMF. The intention of the simulation development is to test various properties and characteristics of LaBr3 to be used in a possible future replacement for the DRAGON BGO gamma array.The code is being built to be generalized for various geometries of LaBr3 crystal that the user can set at run-time. The code will include the internal radiation present in all LaBr3 scintillation detectors. The GPS module of G4 is being used to expose the crystal to various energies of gamma rays. The source characteristics can be controlled at run-time. 
# How to Compile
The user should create a GEANT4 working directory such as:
```
> mkdir G4WorkDir
```
In your work directory create a source directory (e.g. yoursourcedir)
```
> mkdir yoursourcedir
```
Initialize the source directory for use with Git
```
> git init
```
Pull the entire contents of this repository to that directory
```
> git pull origin master
```
Go back to your work directory and create a build directory (e.g yoursourcedir-build)
```
> cd ..
> mkdir yoursourcedir-build
````
Switch to your build directory and run cmake using two arguments, the -D option to point to the directory holding the Geant4Config.cmake file that Geant4 installs to help CMake find and use Geant4, and the second argument is the source directory for the application you want to build. For example, if you were using version 10.3 of G4 and it was installed in the conventional way, do:
```
> cd yoursourcedir-build
> cmake -DGeant4_DIR=/home/username/geant4-install/lib64/Geant4-10.3.0 /home/username/G4WorkDir/yoursourcedir
```
Then, from your build directory make/compile the code
```
> make
```
# How to Run the LaBr3 Simulation
The recommended location for the executable is you build directory and it will be created there by default. The name of the executable can be changed to anything you wish by editting the CMakeLists.txt file and the filename of the executable. If it is unchanged the default name of the executable is "LaBr3_v4" (the files uploaded to Github were the 4th iteration of the code)

The repository includes a folder called "macros" which contains various examples for input files that can be used and modified to run the code under different conditions. When you build with CMake the macros will be copied into a "macros" directory in your build folder. 

If you wish to run an input macro in *batch mode* then do:
```
> ./LaBr3_v4 macros/yourinputfile
```
You can also run an input file in *interactive mode* directly by doing:
```
> ./LaBr3_v4
```
Then on the **G4 "session" line** type:
```
> /control/macroPath macros
> /control/execute vis.mac
> /control/execute yourinputfile
```
You can continue to run multiple macros in the same interactive screen my repeating the /control/execute line as needed.
