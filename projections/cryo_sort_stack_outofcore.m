function cryo_sort_stack_outofcore(instackname,idx,outstackname)
%
% CRYO_SORT_STACK_OUTOFCORE Sort images in a stack
%
% cryo_sort_stack_outofcore(instackname,vals,outstackname)
%   Sort the images in the MRCS stack instackname in descending order of
%   vals. The sorted stack is written to outstackname. The length of vals
%   must equal the number of images in instackname.
%
%   Input parameters:
%       instackname   Name of MRCS file containing input images.
%       idx           New ordering of the images.
%       outstackname  Name of MRCS file into which sorted image are
%                     written. 
%
%   Examples:     
%       cryo_sort_stack_outofcore('instack.mrcs',1000:-1:1,'outstack.mrcs');
%
% Yoel Shkolnisky, March 2018.

instack=imagestackReader(instackname,100);
outstack=imagestackWriter(outstackname,instack.dim(3),1,100);
for k=1:numel(idx)
    im=instack.getImage(idx(k));
    outstack.append(im);
end
outstack.close;

