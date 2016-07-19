This is a MATLAB package for denoising CTF -affected cryo-EM images based on the following manuscripts:
1) Denoising and Covariance Estimation of Single Particle Cryo-EM Images
Tejal Bhamre, Teng Zhang, Amit Singer
http://arxiv.org/abs/1602.06632

2) Fast Steerable Principal Component Analysis
Zhizhen Zhao, Yoel Shkolnisky, Amit Singer
http://arxiv.org/abs/1412.0781

The folder kn_rankest includes code for rank estimation by S. Kritchman and B. Nadler.

DEPENDENCIES
----------------

This package should be used in conjunction with the cryo-EM tool box ASPIRE (http://spr.math.princeton.edu/) and will be included in the latest version of ASPIRE. It also requires the NUFFT package (to be included in the latest version of ASPIRE) available
at http://www.cims.nyu.edu/cmcl/nufft/nufft.html.

INSTRUCTIONS
----------------

1) Download and install ASPIRE from http://spr.math.princeton.edu/ following the instructions for installation.
2) Add ASPIRE files in your MATLAB path using initpath.m in ASPIRE.
3) If this package is in a separate location than ASPIRE, add the package to your MATLAB path using cwf_paths.m  
4) Enjoy the example simulation scripts in cwf_scripts to denoise images.

INSTRUCTIONS FOR NUFFT
---------------------------

wget http://www.fftw.org/fftw-3.3.4.tar.gz
tar -zxf fftw-3.3.4.tar.gz
cd fftw-3.3.4
./configure --prefix=$HOME/local --enable-openmp --enable-shared
make
make install

Then we do the same thing for NFFT:

wget https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-3.2.3.tar.gz
tar -zxf nfft-3.2.3.tar.gz
cd nfft-3.2.3
./configure --prefix=$HOME/local --enable-openmp --with-matlab=/usr/local/Matlab/R2013a --with-fftw3=$HOME/local
make
make install

To run everything in MATLAB, we type

addpath([getenv('HOME') '/local/share/nfft/matlab/nfft']);
addpath([getenv('HOME') '/local/lib']);

and check by running

simple_test;

In case of issues or questions, please email Tejal (tbhamre@princeton.edu) and Jane (jzhao@cims.nyu.edu).
