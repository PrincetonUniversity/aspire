function poolreopen(n)
%
% POOLREOPEN    Open parallerl pool
%
% poolreopen(n)
%   Open parallel pool with n workers. Pool is opened only if no already
%   open.
%
% Yoel Shkolsnisky, July 2015.

if nargin<1
    n=8;
end

[mjv,mnv]=matlabversion; % Get MATLAB version
% Pool commands were changes at version 8.5

if mjv>8 || (mjv==8 && mnv >=5)
    % Use parpool
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers;
    end
    if poolsize==0 || poolsize~=n
        if poolsize~=0
            delete(poolobj);
        else
            parpool(n);
        end
    end
else
    % Use matlabpool        
    poolsize=matlabpool('size');
    if poolsize==0 || poolsize~=n
        if poolsize~=0
            matlabpool close;
        else
            matlabpool('open',n)
        end
    end
end