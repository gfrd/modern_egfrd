
### Prerequisites

1. an ISO C++14 compiler, like GNU GCC v5.0 or Microsoft Visual C++ 2017*
2. GNU Scientific Library (GSL) 1.16*
3. [optional] FreeGlut3 for visualization

*or newer version



### Building on Linux

1. Get prerequisites

If you use Fedora or any other RedHat lines of distributions: 

`sudo yum install gcc-c++ gsl gsl-devel`

or on Ubuntu Linux and some other Debian families, install the following:

`sudo apt-get install gcc-c++ gsl-bin libgsl0-dev`


2. Get source codes: 

```
git clone https://github.com/gfrd/modern_egfrd
	or
download and unzip
```
   
3. Run make to build:
	
`make samples`


4. Execute the sample:

```
cd bin
export LD_LIBRARY_PATH=$(pwd)
./RunGfrd Equilbrium
```


5. Execute UnitTest

```
make tests
cd bin
export LD_LIBRARY_PATH=$(pwd)
./TestGreensFunctions
./TestGFRD
```


6. Optional 3D Visualizer

`sudo yum install freeglut freeglut-devel`
or
`sudo apt-get install freeglut3-dev`

```
make visualize
cd bin
export LD_LIBRARY_PATH=$(pwd)
./gfrdVisualizer
```



This package has been tested on:

| **Distribution** |  **GNU GCC/G++**  | **GSL** |
|--------------|:-----:|:----:|
| CentOS 7     | v5.3.1 | v1.15 |
| Fedora 25    | v6.3.1 | v2.1 |
| Ubuntu 16    | v5.4.0 | v2.1 |
| Debian 8.8   | v4.9.2 | v1.16 |
| Ubuntu 14    | v4.9.2 | v1.16 |
| CentOS 6.6   | v4.9.2 | v1.13 |


For all distributions we have checked so far everthing compiles with zero errors and zero warnings!

If you have compilation problems please let us know!



### Building on Microsoft Windows


1. Get [Visual Studio Community 2017](https://www.visualstudio.com/)
   install "Desktop development with C++" (check if 'WindowsSDK 10.0.15063.0 for Desktop' is selected)

2. Get source codes: 

```
git clone https://github.com/gfrd/modern_egfrd
	or
download and unzip
```

3. Get windows versions of the dependancies: FreeGlut3 / GSL1.16 / LibCCD.
   For convineance we bundeled these in a [support package](http://egfrd.org/includes/packages/WinSupport.zip)
   extract it in the project root, so you get the Libs-folder.
   
	FreeGlut3 for MSVC from [Transmission Zero](http://www.transmissionzero.co.uk/software/freeglut-devel/), 
	GSL for MSVC from Brian Gladman.

3. Open solution newGfrd.sln

4. Select build configuration (e.g. x64 Release) and startup project RunGfrd 

5. Build and Run






