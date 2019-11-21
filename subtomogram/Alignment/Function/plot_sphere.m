function plot_sphere(f, thetas, phis)
% ssht_plot_sphere - Plot function on sphere
%
% Plots a functions defined on the sphere for various sampling schemes.
%
% Default usage is given by
%
%   ssht_plot_sphere(f, L, <options>)
%
% where f is the sampled function and L is the harmonic band-limit.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%  'Method'          = { 'MW'         [McEwen & Wiaux sampling (default)],
%                        'MWSS'       [McEwen & Wiaux symmetric sampling],
%                        'DH'         [Driscoll & Healy sampling],
%                        'GL'         [Gauss-Legendre sampling] }
%  'Type'            = { 'colour'     [generate colour plot (default)],
%                        'parametric' [generate parametric plot] }
%  'Close'           = { true         [close plot in phi (default)],
%                        false        [do not close plot in phi] }
%  'PlotSamples'     = { false        [do not plot sample positions (default)],
%                        true         [plot sample positions] }
%  'ParametricScale' = scale          [scaling for parametric plot (default=0.5)]
%  'ParametricMin'   = { false        [do not clip min to zero for parametric plot (default)],
%                        true         [clip min to zero for parametric plot] }
%  'Lighting'        = { false        [do not light plot (default)],
%                        true         [light plot] }
%  'ColourBar'       = { false        [do not add colour bar (default)],
%                        true         [add colour bar] }
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
%p = inputParser;
%p.addRequired('f', @isnumeric);     
% p.addRequired('L', @isnumeric);          
% p.addParamValue('Type', 'colour', ...
%    @(x)strcmpi(x,'colour') | strcmpi(x,'parametric'));
% p.addParamValue('Method', 'MW', @ischar);
%p.addParamValue('Close', true, @islogical);
% p.addParamValue('PlotSamples', false, @islogical);
% p.addParamValue('ParametricScale', 0.5, @isnumeric);
% p.addParamValue('ParametricMin', false, @islogical);
% p.addParamValue('Lighting', false, @islogical);
% p.addParamValue('ColourBar', false, @islogical);
%p.parse(f, thetas, phis);
%args = p.Results;

% Define parameters.
TOL = 1e-10;
%PARAMETRIC_SCALE = args.ParametricScale;

% Compute grids.
minf = min(f(:));
maxf = max(f(:));
% if args.ParametricMin
%    minf = 0;
% end
   
%[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Method', ...
%                                                args.Method, 'Grid', true);
if (size(thetas) ~= size(f)) 
  error('Size of f does not match sampling');
end

% Compute position scaling for parametric plot.
% if strcmpi(args.Type, 'parametric')
% 
%    if abs(maxf - minf) < TOL
%       f_normalised = f;
%    else
%       f_normalised = (f - minf)./(maxf - minf).*PARAMETRIC_SCALE; % + 0.1;
%       f_normalised(find(f_normalised < 0.0)) = 0.0;
%    end
% else
f_normalised = ones(size(f));
% end

% Close plot.
% if args.Close
%    close = @(x) [x, x(:,1)];
%    f = close(f);
%    f_normalised = close(f_normalised);
%    thetas = close(thetas);
%    phis = close(phis);
% end

% Compute location of vertices.
[x, y, z] = ssht_s2c(thetas, phis, 1.0);
x = x .* f_normalised;
y = y .* f_normalised;
z = z .* f_normalised;

% Plot.
h = surf(x,y,z,f);
if (abs(minf-maxf)>TOL) caxis([minf, maxf]); end
% if args.ColourBar 
%    colorbar('vert'); 
% end
set(h, 'LineStyle', 'none')
% if args.PlotSamples
%     hold on;
%     plot3(x, y, z, 'k.');
%     hold off;
% end
axis([-1,1,-1,1,-1,1]);
view([-1,1,1]);
axis equal;
% if args.Lighting
%   axis off; 
%   light; 
%   lighting phong; 
%   camzoom(1.3);
% end


