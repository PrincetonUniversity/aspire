% Test the function recenter
%
% Same as test_recenter_3 but with a real volume.
%
% Yoel Shkolnisky, November 2014

clear;
load cleanrib;

% Mask the volme. That really improves the results.
n=size(volref,1);
mask = fuzzymask(size(volref),3,floor(n*0.45),floor(n*0.1));
volref=volref.*mask;

[vol1,cm1]=recenter(volref);
[vol2,cm2]=recenter(vol1);

assert(norm(cm2)<1.0e-2);
fprintf('Alignment error %6.4e (pixels)\n',norm(cm2));