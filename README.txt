CUDA code for developmental biology using Subcellular Element Method

Hardware requirement: 
Nvidia video card that supports SM 2.0+ and CUDA 4.0 

Software environment requirement: 
CMAKE ----- Build system.
CUDA  ----- Provide runtime support for parallel GPU computation.
CGAL  ----- Computational geometry library.
Thrust ---- Build-in library of cuda, similar to STL of C++
Paraview -- (Optional) Visualization software for animation purpose. 

************************
To run simulation on slurm cluster (acms-gpu is powered by slurm) 
 (1) In project root folder, cd ./scripts
 (2) sbatch *.sh, for example, sbatch discN01G02.sh means take 
     the first configuration file and then submit it to gpu02 compute node 
     so. The actual GPU device number that is going to run the program is 
     controled by slurm, and specified GPUDevice in the config file is ignored 
     if running on cluster.

Location of configuration files:
 ./resources
*******************************************
To run simulation on UCR HPCC cluster: 
   After uploading all the files under this repository, follow the command below.
   (1) cd [The folder where all your files are placed in]
   (2) module load singluarity
   (3) module load centos
   (4) centos.sh
   (5) module load extra; module load GCC; module load cuda;
   (6) (IF it is the first time compiling) module load CGAL; module load cmake;
   (7) (IF it is the first time compiling) cmake .
   (8) make
   (9) After compilation, enter: exit
   (10) sbatch -p gpu --gres=gpu:1 --time=X:00:00 EpiScale_run.sh; (X here is the number of hours you want to keep the simulation running)
