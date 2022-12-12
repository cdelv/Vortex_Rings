# Basilisk

## Instalation

first, you need to clone the repository. You can do it here. In the .gitignore I have excluded the basilisk source code directory. For details go to http://basilisk.dalembert.upmc.fr/src/INSTALL.  

First get the source code

darcs get http://basilisk.fr/basilisk 

or 

wget http://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz

then enter to the src directory and add it to the path

cd basilisk/src
export BASILISK=$PWD
export PATH=$PATH:$PWD

Now, install de dependancies

sudo apt install flex bison make gawk gfortran
sudo apt-get install bison libglu1-mesa-dev libglew-dev libgl1-mesa-dev

Now, we are going to compile basilisk with MPI. you can use the config file provided in the repository.

cp ../../config.custom .
ln -s config.custom config

Before doing Make, its important to create the configuration file and to compile the OpenGL dependant libraries for the view.h visualizations. 

cd $BASILISK/gl
make 

cd $BASILISK/ppr
make

At last, compile Basilisk.

cd $BASILISK
make

Install some optional but usefull dependencies for visualization.

sudo apt install gnuplot imagemagick ffmpeg graphviz valgrind gifsicle pstoedit



