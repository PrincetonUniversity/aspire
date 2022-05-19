GPU kernel implmenetation of cryo_clmatrix.
The folder contains a series of GPU kernels used both to learn how to write kernels, and to implement cryo_clmatrix step by step.

Kernels are compiled using
	nvcc -ptx <kernel file>.cu
