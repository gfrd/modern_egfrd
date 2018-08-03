
### Prerequisites

1. an ISO standard C++11 compiler, like GNU GCC v5.0* or Microsoft Visual C++ 2017*
2. CMake v3.2* 
3. GNU Scientific Library (GSL) 1.16*
4. [optional] FreeGlut3 for visualization

   *or newer version

### Building on Linux
1. Get prerequisites

   If you use Fedora or any other RedHat lines of distributions: 
   ```
   > sudo yum install gcc-c++ cmake gsl gsl-devel
   ```
   or on Ubuntu Linux and some other Debian families, install the following:
   ```
   > sudo apt-get install gcc-c++ cmake gsl-bin libgsl0-dev
   ```
2. Get source codes: 

   You can git clone the repository
   ```
   > git clone https://github.com/gfrd/modern_egfrd
   > cd modern_egfrd
   ```
   or download and unzip:
   ```
   > wget https://codeload.github.com/gfrd/modern_egfrd/zip/master
   > unzip master
   > cd modern_egfrd-master
   ```
   
 
3. Configure and build:
   
   In general 'out-of-source' builds are preferred, so first make a sub-folder. Then generate and build a Release configuration.
   ```
   > mkdir build
   > cd build
   > cmake -DCMAKE_BUILD_TYPE=Release ..
   > make RunGfrd
   ```

4. Execute the sample:

   ```
   > cd bin
   > ./RunGfrd ../../samples/simple/simple.gfrd
   ```

5. Build and execute tests

   ```
   > make TestGreensFunctions TestGFRD test
   ```

6. Build and execute 3D Visualizer [Optional]

   Get prerequisites for you platform:
   ```
   fedora> sudo yum install freeglut-devel
   ubuntu> sudo apt-get install freeglut3-dev
   ```
   build the project:
   ```
   > make gfrdVisualizer
   > cd bin
   > ./gfrdVisualizer
   ```
   When running press keys f, d, i, s and a for some action, press h for help.

   



### Building on Microsoft Windows



1. Get a [Visual Studio IDE 2017](https://www.visualstudio.com/) (any flavor)

   Run the installer and select "Desktop development with C++" and under 'Individual Components' check 'Visual C++ tools for CMake' and 'Git for Windows'.

   
2. Get the dependencies:

   Install and setup Microsoft Vcpkg (https://github.com/Microsoft/vcpkg)
   
   Get and build the required packages:
   ```
   > vcpkg install gsl:x64-windows freeglut:x64-windows
   ```
   
   Setup a VCPKG_ROOT environment variable:
   ```
   > SET VCPKG_ROOT=<your vcpkg path>
   ```
   or for a more permanent setting, open 'Advanced System Settings' and set the _VCPKG_ROOT_ environment here. Then also check that 'vcpkg' directory and 'CMake' are in the PATH environment variable.
   

3. Get source codes:

   You can clone the repository, either with _git_ or _TortoiseSVN_
   ```
   > git clone https://github.com/gfrd/modern_egfrd
   > cd modern_egfrd
   ```
   or [download](https://codeload.github.com/gfrd/modern_egfrd/zip/master) and unzip.
      

**Workflow I: Command Line**:
   
   
4. Configure and build:
   
   In general 'out-of-source' builds are preferred, so first make a sub-folder. Supply cmake with the desired 'x64' generator and the path to the VCPKG toolchain. Then build 'RunGfrd' in Release configuration.
   ```
   > mkdir build
   > cd build
   > cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake ..
   > cmake --build . --config Release --target RunGfrd
   > 	
   ```

5. Execute the sample:

   ```
   > cd bin\Release
   > RunGfrd.exe ..\..\..\samples\simple\simple.gfrd
   ```

6. Build and execute tests

   ```
   > cmake --build . --config Release
   > ctest -C Release
   ```

7. Build and execute 3D Visualizer [Optional]

   build the project:
   ```
   > cmake --build . --config Release --target gfrdVisualizer
   > cd bin/Release
   > gfrdVisualizer
   ```
   When running press keys f, d, i, s and a for some action, press h for help.
   


**Workflow II: Visual Studio IDE**:
   
		
4. Configure and build:

   Start Visual Studio 2017 and choose 'File', 'Open', 'Folder', select 'modern_egfrd' root folder.

   Visual Studio will process the CMakeList.txt and CMakeSettings.json configuration files and generate the projects.
   
   When ready click 'CMake' from the menubar and choose 'Build All'.

   
5. Execute the sample:   
   
   To run the sample, the command line arguments  need to be specified. To do this click 'CMake' from the menubar and choose 'Debug and Launch Settings' for 'RunGfrd'.
   In the 'launch.vs.json' file add a comma and the following line after "name":
   
   ```
   "args": ["${workspaceRoot}\\samples\\simple\\simple.gfrd"]
   ```
 
   Then click 'CMake' from the menubar and choose 'Debug', 'RunGfrd'.
   or use the 'Select Startup Item'-toolbar to choose 'RunGfrd' and click it to start.
   

6. Build and execute tests
   
   Click 'CMake' from the menubar and choose 'Run Test' for 'CMakeLists.txt (Project modern_egfrd).

   
7. Build and execute 3D Visualizer [Optional]

   Then click 'CMake' from the menubar and choose 'Debug', 'gfrdVisualizer'.
   or use the 'Select Startup Item'-toolbar to choose 'gfrdVisualizer' and click it to start.
   
   When running press keys f, d, i, s and a for some action, press h for help.
   

   
### Building on Apple MacOS

1. Todo



### Compliance


   This package has been tested on:

   | **Distribution** | **CMake** | **C++** | **GSL** |
   |------------------|:---------:|:-------:|:-------:|
   | CentOS 7         | v3.6.3    | v5.3.1  | v1.15   |
   | Fedora 25        |           | v6.3.1  | v2.1    |
   | Ubuntu 16        | v3.5.1    | v5.4.0  | v2.1    |
   | Debian 8.8       |           | v4.9.2  | v1.16   |
   | Ubuntu 14        | v3.2.2    | v4.9.2  | v1.16   |
   | CentOS 6.6       |           | v4.9.2  | v1.13   |
   | Win7 / VS2017    | v3.7      | v19.10  | v1.16   |
   | Win10 / VS2017   | v3.7      | v19.10  | v2.3    |

   
   For all distributions checked so far everything builds with zero errors and (almost) zero warnings!
   
   If you have compilation problems please let us know by filing an [issue](https://github.com/gfrd/modern_egfrd/issues)

