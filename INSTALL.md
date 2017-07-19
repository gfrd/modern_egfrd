
### Prerequisites

1. an ISO C++14 compiler, like GNU GCC v5.0 or Microsoft Visual C++ 2017*
2. GNU Scientific Library (GSL) 1.16*
3. [optional] FreeGlut3 for visualization

*or newer version



### Building on Linux

1. Get prerequisites

If you use Fedora or any other RedHat lines of distributions: 

```
> sudo yum install gcc-c++ gsl gsl-devel
```

or on Ubuntu Linux and some other Debian families, install the following:

```
> sudo apt-get install gcc-c++ gsl-bin libgsl0-dev
```


2. Get source codes: 

```
> git clone https://github.com/gfrd/modern_egfrd
```
   
3. Run make to build:

```
> make samples
```


4. Execute the sample:

```
> cd bin
> export LD_LIBRARY_PATH=$(pwd)
> ./RunGfrd Equilibrium -ka 1e-19 -kd 2e-2 -p 100 -e 200 > data.out
```


5. Build and execute UnitTests

```
> make tests
> cd bin
> export LD_LIBRARY_PATH=$(pwd)
> ./TestGreensFunctions
> ./TestGFRD
```


6. Optional run 3D Visualizer

Get prerequisites for you plaform:
```
fedora> sudo yum install freeglut freeglut-devel
ubuntu> sudo apt-get install freeglut3-dev
```

build the project:
```
> make visualize
> cd bin
> export LD_LIBRARY_PATH=$(pwd)
> ./gfrdVisualizer
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


For all distributions we have checked so far everthing compiles with zero errors and (almost) zero warnings!

If you have compilation problems please let us know!



### Building on Microsoft Windows


1. Get [Visual Studio Community 2017](https://www.visualstudio.com/)


   install "Desktop development with C++" (check if 'WindowsSDK 10.0.15063.0 for Desktop' is selected)

   
2. Get source codes: 

```
> git clone https://github.com/gfrd/modern_egfrd
```

3. Get the dependencies: 

**Method A**: (vcpkg)

    Install and setup Microsoft Vcpkg (https://github.com/Microsoft/vcpkg)

    Get and build the packages:
   
   
    ```
    > vcpkg install gsl:x64-windows freeglut:x64-windows
    ```
    
	If you wish to do so, you can also install the x86-windows (32bit) versions of these packages.

    Vcpkg automatically handles include directorties, library directories, linkage and runtime-dll copying!


**Method B**: (manual)
   
    Download the [WinSupportPackage](http://egfrd.org/includes/packages/WinSupport.zip)
    and extract it to some `<path>` you prefer.
   
    For project(s): genBesselTables, gfrdVisualizer, RunGfrd, TestGFRD, TestGreensFunctions, eGFRD and greensfunctions:
        *   Add '`<path>`\gsl-1.16' to the VC++ include directories.
        *   Add '`<path>`\gsl-1.16\lib\[win32/x64]' to the VC++ library directories.
        *   Add 'gsl.lib' and 'cblas.lib' to the linker additional dependencies.
        *   Copy 'gsl.dll' and 'cblas.dll' from '[Win32\X64]\[Debug/Release]' to '_bin'-target folder.

    For project(s): gfrdVisualizer:
        *   Add '`<path>`\freeglut\include' to the VC++ include directories.
        *   Add '`<path>`\freeglut\lib\[win32/x64]' to the VC++ library directories.
        *   Add 'freeglut.lib' to the linker additional dependencies.
        *   Copy 'freeglut.dll' from '[Win32\X64]\[Debug/Release]' to '_bin'-target folder.

    The WinSupportPackage contains prebuild versions of the packages, created by:
        * FreeGlut3 for MSVC from [Transmission Zero](http://www.transmissionzero.co.uk/software/freeglut-devel/), 
        * GSL for MSVC from Brian Gladman.

		
3. Open solution modern_egfrd.sln

4. Select build configuration (e.g. x64 Release) and startup project RunGfrd 

5. Build and Run






