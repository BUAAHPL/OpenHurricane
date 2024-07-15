# OpenHurricane

*Open parts of Hurricane project*

## OpenHurricane introduction

OpenHurricane is an open-source library written in C++ 17 and CUDA for computing supersonic reacting flows.
In the current version, the OpenHurricane is able to compute chemistry source terms and gas physical properties on the CPU and GPU platforms alternatively.
OpenHurricane is in its development phase, and there would be more features in it.

OpenHurricane is built using a typical ```camke/make/make install``` process. When make install is complete, the built library will be installed in the "lib" folder.
There also exists an example in the "examples" folder which is about how the functions of the library are called.


## Build

Building OpenHurricane requires C++ compiler and [CMake 3.18](https://cmake.org/ "CMake") or later.

### Dependencies


* [CGNS](https://github.com/CGNS/CGNS "CGNS"): a general, portable, and extensible standard for the storage and retrieval of CFD analysis data.
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/ "hdf5"): a set of file formats (HDF4, HDF5) designed to store and organize large amounts of data, needed by CGNS and OpenHurricane.
* [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/download "Metis"): a set of serial programs for partitioning graphs, partitioning finite element meshes.
* [MPI](https://computing.llnl.gov/tutorials/mpi/ "MPI"): Massage Passing Interface.A standardized and portable message-passing standard for parallel computing.
  Users can use [MS-MPI](https://github.com/Microsoft/Microsoft-MPI "MS-MPI"), Intel MPI from [Intel oneAPI toolkits](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html, "Intel oneAPI"), [MPICH](https://github.com/pmodels/mpich "MPICH") and [OpenMPI](https://github.com/open-mpi/ompi "OpenMPI") packages.
* [Eigen](http://eigen.tuxfamily.org/ "Eigen"): Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
* [CUDA](https://developer.nvidia.com/cuda-downloads "CUDA"): CUDA is a parallel computing platform and programming model that makes using a GPU for general purpose computing simple and elegant.

### Windows

For Windows platform, the [Visual Studio Community](https://visualstudio.microsoft.com/ "Visual Studio IDE") is recommended.
If the C++ Cmake tools for Windows is installed with the Visual Studio, then users can use Visual Studio directly to run Cmake configuration and build applications for OpenHurricane.

Users can also use the [CMake](https://cmake.org/, "CMake") tool to configure the project of OpenHurricane.
However, in-source builds are not allowed.

#### Using CMake GUI

To do this, run the cmake-gui.exe, which should be in your Start menu or on your desktop,
and browse to the location of OpenHurricane's
source with the "Browse Source" button. You can also change the binary
directory, where the Visual Studio project files will be placed, with the "Browse Build" button. 
Click "Generate" to select the correct visual studio version and build the project files.
Then the Visual Studio project will be generated.

#### Using command line

Open the command prompt and cd to the OpenHurricane source directory.
Make a directory (e.g., *build*) where the resulting binaries should be placed:

    > mkdir build
    > cd build

Run

    > cmake --help

and look at the list of generators for the Visual Studio you
want to build for. For example, the generator for Visual Studio 2022
is called "Visual Studio 17".

After you have found the appropriate generator, run

    > cmake .. -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release

to generate the project files. The project files will be placed in directory *build*.
The Visual Studio project will be called OpenHurricane.sln. Open it in Visual
Studio. 

#### Using Visual Studio IDE

To do this, run the "Visual Studio 2022.exe", which should be in your Start menu or on your desktop,
and click "Open a local folder" to select the location of OpenHurricane's source directory. 
Then the project will be automatically configured by the Visual Studio IDE with the CMake.

### Linux

To compile OpenHurricane in Linux systems, users should use the [CMake](https://cmake.org/ "cmake") tool 
to configure the project and then use the [GNU Make](https://ftp.gnu.org/gnu/make/ "make") tool 
to control the generations of the applications. 
The following additional requirements must be met in your Linux systems:

#### Require in Linux
* [gcc](https://ftp.gnu.org/gnu/gcc/ "gcc"):  OpenHurricane has been only tested by using the GCC compiler in Linux systems. Therefore, the GCC compiler is recommended for Linux platforms to build OpenHurricane.
        Require version 9.0.0 or later.
        Check with command ```gcc --version```
* [cmake](https://cmake.org/ "cmake"): Require version 3.18 or later.
        Check with command ```cmake --version```
* [make](https://ftp.gnu.org/gnu/make/ "make"): GNU make. Require version 3.82 or later. 
        Check with command ```make --version```

The source codes are already together with prebuilt third-party libraries.
You can also build them by yourself.
 If users want to use the prebuilt third-party libraries in Linux systems,
 then the below requirements should be met.
* [GNU binutils](https://ftp.gnu.org/gnu/binutils/ "GNU binutils"): GNU binutils (including ```as``` and ```ld``` ), require version 2.27 or later. 
        Check with command ```ld --version```
* [ldd](https://ftp.gnu.org/gnu/glibc/ "ldd"): Require version 2.17 or later. 
        Check with command ```ldd --version```

Then users can build OpenHurricane with CMake or CCMake tools.
In-source builds are not allowed.
#### Using CMake from the command line

The example of using the ```cmake``` command for building OpenHurricane object is

    > mkdir build
    > cd build
    > cmake .. -DCMAKE_BUILD_TYPE=Release
    > make && make install


#### Using CCMake from the command line

If the CMake curses interface is supported in your Linux platform, then the ```ccmake``` command can be used to configure this project.
The example of using the ```ccmake``` command for building OpenHurricane object is

    > mkdir build
    > cd build
    > ccmake ..
    > make && make install

Thanks for building, and happy optimizing!

## OpenHurricane DEVELOPERS

The current OpenHurricane project was developed by **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")'s** research group at ***Beihang University***.

If you have any problems in building or running the code or any suggestions, please do not hesitate to contact.

 * **Ph.D Rao Sihang**. Email:<raosihang@hotmail.com>
 * **[A.P. Chen Bing](http://shi.buaa.edu.cn/chenbing/zh_CN/index.htm "A.P. Chen Bing")**. Email:<MarkChien@buaa.edu.cn>
 * **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")**. Email:<xuxu_1521@126.com>

## Current Maintainer

The OpenHurricane project is maintained by  

* **Ph.D Rao Sihang** of **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")'s** group at ***Beihang University***. Email:raosihang@hotmail.com
* **Ph.D Candidates Peng Jian** of **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")'s** group at ***Beihang University***. Email:p308224999@outlook.com
* **Ph.D Candidates Yang Hongzhen** of **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")'s** group at ***Beihang University***. Email:826134437@qq.com
* **Ph.D Candidates Zhang Jiebo** of **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")'s** group at ***Beihang University***.
* **Zheng Putian** of **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")'s** group at ***Beihang University***.
* **Qin Xu** of **[Prof. Xu Xu](https://shi.buaa.edu.cn/xuxu/en/index.htm "Prof. Xu Xu")'s** group at ***Beihang University***.


## License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html "GNU General Public License version 3")

Copyright (C) 2019-2024, Prof. Xu Xu's group at Beihang University   
