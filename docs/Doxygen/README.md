# To build code guides for OpenHurricane with Doxygen tool 

*Open parts of Hurricane project*

## Install Doxygen

The code guides document for OpenHurricane can be built by using the [Doxygen](https://www.doxygen.nl/index.html "Doxygen software") software.
Images in the document are generated using the [Graphviz](https://graphviz.org/ "Graphviz software") software.
To build the documentation either on the Ubuntu GNU/Linux system or on the Windows system, the user should install both [Doxygen](https://www.doxygen.nl/index.html "Doxygen software")
and [Graphviz](https://graphviz.org/ "Graphviz software") package.

## Configure Doxygen file

The Doygen configuration file, Doxyfile, in the ```OpenHurricane/docs/Doxygen```
directory is configured to work with Doxygen versions 1.10.0.

## Run Doxygen

In the ```OpenHurricane/docs/Doxygen``` directory run
> doxygen
