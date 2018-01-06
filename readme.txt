ASPIRE – Algorithms for Single Particle Reconstruction

Current version: 0.3
Date: 10/2017

Installation
------------
The following steps should be executed only once.
1. Download and extract the package.
   We will assume that the package has been extracted to a folder named ASPIRE.
2. Change into the directory ASPIRE and run “install”.
3. Change into the directory ASPIRE/projections/class_average/ and run “gen_simulation_data”.


Getting started
---------------
1. Enter the directory ASPIRE.
2. Run “initpath”
    This script should be executed each time you start your Matlab session.
3. Change into the directory ASPIRE/examples
4. Enjoy the example scripts


Revisions
---------
*Changes from version 0.1
1.	Updated class averaging code.
2.	Revised the function cryo_project to generate projections whose size is different from the size of the projected volume.


*Changes from version 0.11
1. 	Added memory efficient FIRM reconstruction routine reconstruction/FIRM/recon3d_firm_ctf_large.m 
2.	Added example file for using cryo_project examples/simulated_projections.m
3.	Update documentation of projections/simulation/cryo_project.m
4. 	Remove obsolete benchmark results from comment at the end of  ./sinograms/test_commonlines_gaussian.m


*Changes from version 0.12
1. Added improved denoising for images
2. Added volume covariance estimation
3. Added improved noise estimation
4. Added Fast steerable PCA
5. New example scripts
6. Improved ab-initio reconstruction for non-symmetric molecules
7. Added routines for reconstruction workflow

