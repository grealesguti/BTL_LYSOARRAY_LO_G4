# README

## BTL LYSO Scintillation Simulation

This code creates simulations of the impact of Na gamma rays and Muons onto a scintillator.

The output provides information of the Light Output and timing resolution of the attached SiPMs.

Furthermore the code is expected to work with a NSGAII genetic algorithm to study the shape optimization of the crystals.

## Installation and running

Download the latest code using:

`git clone --single-branch --branch g/f/G4asfunction/src/TierII-gmsh https://github.com/grealesguti/BTL_LYSOARRAY_LO_G4.git`

The code is expected to work in conjonction with Gmsh for the shape of the crystals.

### Apptainer/Singularity installation

Apptainer, previously called singularity, provides a conteinerized environment that can run this code without any need for local installation. This provides the advantage to work in any operating system independently of any future updates at the cost of slighlty lower efficiency.

The singularity file ( https://drive.google.com/file/d/1Yjz-a1oVOpPOJ554Z07GvCtMKcX5Dai7/view?usp=sharing , request access if needed) provides the required libraries to run the code.

An example installed and run using:

`
singularity exec <path-to-the-singularity-file> make -f Makefile.localgmsh
singularity exec <path-to-the-singularity-file> ./sim -GeomConfig 1 -runevt 10
`

or with a more complex example and a given path:

`singularity exec ../../SingDir/sandG4Gmsh ./sim -GeomConfig 11 -runevt 10 -Volume -nDetected -incrSiPM 125 -LYSO_L 20 -ResinMach -rnd 0 -Znode 1 -Ypos {1-1} -Acte -BC400`


### Local installation (2023.10)

To install locally we need to install GMSH, GEANT4 and their dependencies. 

The compilation and running is similar as in the previous case without the need of singularity/apptainer.


#### GEANT4 installation


#### GMSH installation



### How to run under a singularity file (Singularity has changed its name!! Apptainer, same fucntionality)

Download the singularity file from:
`singularity build <singularity_file_name>.sif docker-daemon://<docker_registry>/<repository>:<tag>`


The singularity file needs to be in your node as a sandbox. Once downloaded the file remember to run:
`singularity build --sandbox <sandboxname> <path/singularity_file_name.sif>`

Now introduce the path to the singularity file within your submission file
`+Singularity=<path>`

To run it locally remember to install singularity:
https://www.linuxwave.info/2022/02/installing-singularity-in-ubuntu-2004.html

### Running in TierII

Login into TierII: `username@login-1.hep.caltech.edu`

Set your respective proxy: `voms-proxy-init -voms cms -rfc -valid 192:0 --bits 2048`

Run using one of the example submission and job files in their respective folders (SubFiles & JobFiles folders).

### Running numerical optimization (Runs In TierII)

First you will need to update the submodule git project. If it's **the first time** you check-out a repo you need to use `--init` first:

```
git submodule update --init --recursive
```

For **git 1.8.2** or above, the option `--remote` was added to support updating to latest tips of remote branches:

```
git submodule update --recursive --remote
```

The new submodule will be downloaded in NSGAII/PyNSGA/NSGA-II, from where the optimization needs to be run, `cd NSGAII/PyNSGA/NSGA-II`.

We need to have root in our python environment
```
export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases/
eval "` $LCGENV_PATH/lcgenv/latest/lcgenv -p LCG_102rc1 x86_64-centos7-gcc11-opt ROOT `"
```
Navigate to the folder NSGAII/PyNSGA/NSGAII and run:

`python example/G4example.py`

This file contains the python scripts that executes the NSGAII algorithm.

It is useful to run under tmux to avoid the code from stopping if you log out. An example to how to open a new tmux window:
`tmux new -s sessionname`

The root file that contain the resulting individuals can be found in NSGAII/PyNSGA/NSGAII/ROOT

Standard use is with 100 individuals. Notice that the algorithm can break if we use Short jobs that can not run for longer than 2h. In these cases try to reduce the duration of the simulation with simpler geometries.
- TODO : penalize with large objectives if the files are not found or can not be opened. remove all submited jobs before continuing.

## Optional arguments

The commands work as follow `-command argument`, in `[]` we set the number of options. The <> defines the type of input, if no <> is given no further values are required for this command.

-GeomConfig <int> : options [1 2 3 11 13] This argument changes the configuration of the detector. Configuration 1-3 provide use the default geometry construction of G4 while 11 and 13 use Gmsh and allow shape changes through splines and the -YPos argument. Config. 3 and 13 use 16 crystals while the rest use a single one.
-Zelem <int> : number of mesh nodes along half of Z of the LYSO crystal
-Znode <int> : number of sections of the LYSO mesh, number of <int> within the -Ypos argument minus 1
-Ypos {-<int>-<int>} : value that multiplies the 3mm, by default height of the LYSO, of each spline control node of the mesh.
-incrV <int> : range [1-199] value between 1 and 199 that provides the 1D study of 2 half tetrahedra for the crystal shape decreasing and increasing the LYSO thickness in the middle and edges of the crystal

-o <string> :  name identifies the name.root file created in the Results folder with the information of the simulation.
-runevt <int> : argument that provide how many events will be run. If not given and if you use G4vis.cc you can use the graphical interface.

-rnd <int> <int> <int> : taking 3 arguments in the form of 1 or 0. The second argument sets no random position to the particle impacts, the second sets no random geometry tolerances during the detector construction.
-rndangle : the particle impacts in random directions and not only along the Z direction
-gunmesh <int> <int> : provides a uniform impact pattern in a quarter of the crystal along X and Z. With the given defaults it would run 60 events.
-Muon : sets the particle impact to 2MeV Muons.

-ESRbackpainted : changes the surface coating model so that there is an air gap between the ESR and the LYSO
-ESRdefaultmodel: changes the default G4 model with a polished LYSO, air gap and ESR (LUT DAVIS)
-noESR : no ESR coating is applied to the crystal

-LYSO_L <int> : default 28.5 where the number given is half the lenght of the LYSO crystal along Z
-Acte : command that maintains the XxZ area constant increasint X if LYSO_L is changed
-matchSiPM : command that sets the SiPM width X equal to the LYSO crystal width
-incrSiPM  <int> :range [25-200], increases in percentage the height of the SiPM from its default 3mm 
-SiPM_Z <int>
-RESIN_Z <int>
-RESIN_Z1000 <int>
-Glue_Z <int>
-Glue_Z1000 <int>
-RESIN_W <int>

-GmshView: opens Gmsh before running G4. Right now does not allow to continue the simulation.
-SaveSTL: saves an .stl file of the mesh for further postprocessing

-BC408: changes the crystal material from LYSO to BC408 (polymer).

-TileV0 : changes configuration from bar to tile. (TODO: right now the design is assymmetrical, the triangulation needs to be changed in Gmsh)
-ForceBottomLine: To always be used with -TileV0. Forces a given flat bottom surface.
-SZloc <double> : range [0-1] changes the SiPM location along Z at the bottom of the tile as a percentage of LYSO_L.

## How to start the graphical interface

In the sim.cc file change 

//#include "src_G4/G4simTierII.hh"
#include "src_G4/G4sim.hh"

and 

//G4simulationNOVIS *sim1 = new G4simulationNOVIS(runManager,argc, argv, Onode, Znode, radinit);
G4simulation *sim = new G4simulation(argc, argv);

now remake and compile. If no -runevt is given you should be running the gui.

THIS CAN NOT WORK IN TIERII or HPCs in general!!!

