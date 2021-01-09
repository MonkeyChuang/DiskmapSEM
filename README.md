# DiskmapSEM
DiskmapSEM implements the disk-shaped SEM algorithm in [1] using Matlab code. In brief, the disk-shaped SEM algorithm produces the equiareal parametrization between closed simply connected surface and unit disk.
##  Getting Started
Run the `demo_SEM.m` file will give a quick and thorough navigation to the disk-shaped SEM algorithm.

The demo code includes several parts:

1. Load 3D mesh data in folder `data`:  
	```matlab
	M = load(filename);
	```
1. Preprocess the 3D data:
	this step is to be encapsulated in one function in the future release.
1. Calculate the equiareal map:
	```matlab
	[uv,C,D] = DiskmapSEM(F,V);
	```
1. Compute the evaluation metric(評價指標) of SEM algorithm
1. Display the resulting equiareal map.

## Release History
- 1.0
	- Add: The most basic version of DiskmapSEM, PlotMesh, and others. An introductory demo code: demo_SEM.
- 1.1
	- Change: Comments/Documentation of demo_SEM and DiskmapSEM
## Contact
Pei-Sheng Chuang - pschuang96@gmail.com
## Reference
[1] M.-H. Yueh, W.-W. Lin, C.-T. Wu, and S.-T. Yau, A Novel Stretch Energy Minimization Algorithm for Equiareal Parameterizations, J. Sci. Comput. 78(3): 1353–1386, 2019.