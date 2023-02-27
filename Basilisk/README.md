# Basilisk

## Compile a Basilisk Program

We have set up a makefile to help with the compilation process. First, copy the `basilisk_config.mk` template.
```
cp basilisk_config.mk local_config.mk
```
 `local_config.mk` is ignored from the repo syncronization. Then, set up the necesary fields on  `local_config.mk`.
 
 After that, go to a directory and type
 ```
make
```

## Install Basilisk

First, you need to clone the repository. You can do it here. In the .gitignore I have excluded the Basilisk source code directory. For more details, go to http://basilisk.dalembert.upmc.fr/src/INSTALL.  

- Get the source code:

```
darcs get http://basilisk.fr/basilisk 
```

or 

```
wget http://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz
```

- Enter to `basilisk/src` directory and add it to the path:

```
cd basilisk/src
export BASILISK=$PWD
export PATH=$PATH:$PWD
```

- Install de dependencies:

```
sudo apt install flex bison make gawk gfortran
sudo apt-get install bison libglu1-mesa-dev libglew-dev libgl1-mesa-dev
sudo apt install gnuplot imagemagick ffmpeg graphviz valgrind gifsicle pstoedit bison libglu1-mesa-dev libosmesa6-dev
```

The last ones are for visualization and are optional. We are going to compile Basilisk with MPI. For that, you can use the config file provided in the repository.

- Create configuration file:

```
cp ../../config.custom .
ln -s config.custom config
```

Before doing Make, compile the OpenGL libraries for basilisk `view.h`. You can also compile the optional `ppr` libraries.

- Compile visualization and optional libraries:

```
cd $BASILISK/gl
make 

cd $BASILISK/ppr
make
```

- At last, compile Basilisk:

```
cd $BASILISK
make -k
make
```
