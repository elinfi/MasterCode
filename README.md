# An assessment of methods for determining and visualiazing differential Hi-C contact maps

Code related to my mester's thesis, titled "An assessment of methods for determining and visualiazing differential Hi-C contact maps".

To run the code, Hi-C data is needed. This Hi-C data needs to be stored in a `.cool/.mcool` file. The entire code is written in `Python` with the use of packages supporting this file format, including `Cooler`and `cooltools`.

A simple clustering procedure is included using `FasterPAM`, provided by the python package `kmedoids`.
