%NFFT_INIT_GURU Initialise plans, no error handling
%   Matlab might run into a segmentation violation for wrong parameters
%
%   nfft_init_guru(d,N1,...,Nd,M,n1,...,nd,m,nfft_flags,fftw_flags)
%
%   d            spatial dimension
%   N1,...,Nd    bandwidths
%   M            number of nodes
%   n1,...,nd    fft lengths
%   m            cut-off parameter
%   nfft_flags   PRE_PHI_HUT | {FG_PSI, PRE_LIN_PSI, PRE_FG_PSI, PRE_PSI,
%	             PRE_FULL_PSI} | FFT_OUT_OF_PLACE
%   fftw_flags   {FFTW_ESTIMATE, FFTW_MEASURE}
%
%   Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts

% Copyright (c) 2002, 2012 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%
% $Id: nfft_init_guru.m 3776 2012-06-03 13:29:25Z keiner $
function p = nfft_init_guru(varargin)

p = nfftmex('init_guru',varargin);
