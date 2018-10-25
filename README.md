# PP

The code is an implementation of the siggraph 2018 paper [Progressive Parameterizations](http://staff.ustc.edu.cn/~fuxm/projects/ProgressivePara/index.html) .

1.	Download and install Mingw-w64 from [here](http://www.mingw-w64.org/doku.php/download).
2.	Download and install Cmake from [here](https://cmake.org/download/).
3.	Apply for download following the official instructions from https://pardiso-project.org/, it will send you an email which contains the download link and the license number of Pardiso v6. 

	- Create a new file named as pardiso.lic and copy the license number into it.  
	- Download the pardiso library from the download link provided by the email.  
	- Put the pardiso.lic and pardiso library(libpardiso600-WIN-X86-64.dll) into the same folder.  
	- Set the environment variable PARDISO_LIC_PATH to the path of that folder(use "/" instead of "\\" on Windows when filling the path).  
	- The environment variable OMP_NUM_THREADS should be set according to your CPU cores number.  
4.	Download the OpenMesh from [here](http://www.openmesh.org/media/Releases/7.1/OpenMesh-7.1.zip) and compile it with Mingw-w64, copy the libOpenMeshCore.a and libOpenMeshTools.a to the folder "Library/".
5.	Other dependences such as libgfortran-3.dll, libgfortran-4.dll, libopenblas.dll are provided in the folder “build/”.
6.	Cmake with the generator options="MinGW Makefiles", mingw32-make, and then it will generate the executable file named PP.exe in the "build/".
7.	Decompress the data file("data/*.zip") into the "data/" folder directly, run the batch processing file "build/run.bat". The unfolded meshes of all the models appear in the paper will be produced and saved into "build/result/".

