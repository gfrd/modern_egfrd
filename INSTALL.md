
### Prerequisites

1. an ISO C++14 compiler, like GNU GCC v5.0 or Microsoft Visual C++ 2017
2. GNU Scientific Library (GSL) 1.16.


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


has been tested on:
* Ubuntu 14.04 LTS
* Ubuntu 10.04 LTS
* Debian 7.7.0
* Fedora 20.1
* CentOS 7.0.14.06


### Building on Windows

1. Get Visual Studio 2017 (...)

2. Get FreeGlut / GSL1.16 / libccd  (link to win support package)
   extract under folder Libs

3. Open solution newGfrd.sln

4. Select startup project (RunGfrd)

5. Build and Run






