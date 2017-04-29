# r3dm

ADAPTIVE GRID REGISTRATION ALGORITHM WITH MULTITHREAD (R3D- Multi) MANUAL
====
By `Gustavo` & `Eason`

Introduction
----
  This document describes the usage of the Adaptive Grid Algorithm registration 
  software improved by EasonYu in the summer of 2017 with QT. Given two input images,
  the software finds the nonrigid geometric transformation that best matches the spatial
  coordinates of one image with those of another. The Normalized Mutual Informayion 
  is used as the measure of similarity during optimization. A detailed description of
  the Adaptive Grid Algorithm used here can be found in Gustavo’s masters’ thesis 
  entitled “The Adaptive Grid Algorithm: A New Spline Modeling Approach for Nonrigid 
  Registration of Medical Images.” Note that nonrigid alignment is a computationally
  complex problem. Thus, expect long execution times when running on CPUs of low clock
  speed. Also note that for best results it is advantageous to align the images through
  affine transformations before attempting nonrigid registration. 
  
  
Usage
----
	To start the program 
      
    simply type ./r3dmulti inputfile at the command prompt.

	Here “inputfile” refers to an ASCII text file that provides the program with the 
  location of the images as well as the registration parameters to be used.
  
Input File 
----
	The following is an example input file that can be used for nonrigid registration. Note
  that each line of this file is expected to contain specific information. The lines in 
  the example file are numbered and the input expected for each line is described below.

  1.	SOURCE:
  2.	s01_r.im
  3.	128
  4.	TARGET:
  5.	ref.im
  6.	0
  7.	
  8.	256
  9.	256
  10.	124
  11.	PARAMETERS:
  12.	
  13.	7
  14.	2
  15.	2
  16.	2 3 5 9 17 32 60
  17.	2 3 5 9 15 28 44
  18.	9 17
  19.	255
  20.	0
  21.	32
  22.	32
  23.	0.05
  24.	
  25.	CONFIGURATION:
  26.	2
  27.	0
  28.	1.7
  29.	0.0001
  30.	0.0005
  31.	0.002
  32.	0.1
  33.	0.0001
  34.	1 1 1
  35.	1
  36.	
  37.	OUTPUT:
  38.	s01_nr.im

Difination of inputfile:
----

1. “SOURCE:” in ASCII characters.
2. The file name of the source image. The source image must be a 3D image stored in a rowXcolumnXslice manner. Each pixel must be described by a two byte integer. Negative intensity values are not allowed. Note that the file extension must be “.im” for the software to work correctly. Lastly, the file name must not contain the character “.”, with the only exception being the extension “.im”.
3. The number of bytes used as a header in the source image.
4. “TARGET:” in ASCII characters
5. The file name of the target image. The target image must be a 3D image stored in a
rowXcolumnXslice manner. Each pixel must be described by a two byte integer. Negative intensity values are not allowed. The image dimensions must be the same as those of the source image. Note that the file extension must be “.im” for the software to work correctly. Lastly, the file name must not contain the character “.”, with the only exception being the extension “.im”.
6. The number of bytes used as a header in the target image.
7. Blank line.
8. The number of pixel in the x dimension of the source and target images.
9. The number of pixel in the y dimension of the source and target images.
10. The number of pixel in the z dimension of the source and target images.
11. Blank line.
12. “PARAMETERS:” in ASCII characters.
13. The number of levels to be used by the algorithm.
14. The number of resolutions to be used. 0 specifies that multiresolutions will not be used and the algorithm will use the images in their original state to obtain solution. Note that multiresolution is done by a factor of two. Thus each dimension of the images being used must be exactly divisible by 2 n , with n referring to the number of resolutions to be used.
15. Maximum image resolution with which the algorithm will work. If zero, for example, the algorithm will only use one image resolution to achieve solution. The resolution used will that in which the input image dimensions are divided by 2 n , with n referring to the number of resolutions to be used specified in line 14. If the this number is equal to the number in line 14, the algorithm will be allowed to work at the full image resolution.
16. Basis function grid resolution to be used in each level in the x and y directions. Note that the size of the sequence of numbers in this line should equal the number of levels specified in line 13.
17. Basis function grid resolution to be used in each level in the z. Note that the size of the sequence of numbers in this line should equal the number of levels specified in line 13.
18. A sequence of numbers determining when the images should be upsampled during
optimization. The images are upsampled when the basis function grid resolution (in the x and y directions) reaches a value specified in this line. In this example file, the images will be upsampled after the basis function grid of size 9X9X9 has been optimized. The images will be upsampled again when the optimization of the grid of size 17X17X15 has been optimized.
19. Maximum pixel value to allow in the cost function computation
20. Minimum pixel value to allow in the cost function computation
21. Number of bins to use in the computation of the source image histogram.
22. Number of bins to use in the computation of the source image histogram.
23. Threshold to be enforced on the Jacobian of the transformation at each level.
24. Blank line
25. “CONFIGURATION:” in ASCII characters.
26. Values allowed: 0, 1, 2. If zero, use the regular grid approach. If one, use adaptive grid approach. If two, use adaptive grid multi-thread approach.
27. If “0” normal approach. If “1” place basis functions only within the support of the foreground of the source image.
28. Parameter specifying the size of the basis functions to be used in each level. More precisely, two times the amount (measured by the radius of the basis functions being used) by which one basis function overlaps with another.
29. Epsilon to be used in the procedure to bracket one optimization interval. Function values must increase by more than epsilon from the minimum for the interval to be declared.
30. Epsilon used to declare the function minimum. If two consecutive steps fail to decrease the cost function by more than epsilon optimization is halted.
31. Size of initial step used during interval bracketing. More precisely, the size of initial step used during interval bracketing will be the value in this line times the magnitude of the gradient of the cost function with respect to the basis functions being used.32. Size of step in used in finite difference gradient evaluation.
33. Threshold based on which points will be chosen as hot spots for optimization.
34. Sequence of three numbers, each being either zero or one. In this example, optimize in each x, y, and z directions.
35. Zero or one. If zero, will not save the deformation fields that register the source to the target image. If one will save them. The deformation fields will be saved with the same file name given in line
38. The x direction deformation field will be have the “.dx” extension, and so on.
36. Blank line.
37. “OUTPUT:” in ASCII characters.
38. The file name where the output image will be saved. The output image will be the source image resampled to match the given target image.

Thanks for downloading
----
