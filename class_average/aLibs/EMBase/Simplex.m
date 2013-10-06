function [newPars, t]=Simplex(y, oldPars, steps, mask)
% function [newPars, t]=Simplex(y, oldPars, steps, mask)
% Nelder-Mead Simplex minimization, implemented as a state machine.
% Usage:
% [newPars, t]=Simplex('init', startPars, steps);  % Initialization (in this
%   case the returned newPars are equal to startPars)
% [newPars, t]=Simplex(y);                      % Iteration
% ok=Simplex('converged', epsilon);             % Test for convergence
% finalPars=Simplex('centroid');                % Obtain final value
% 
% - startPars is the vector of starting parameter values (or alternatively a
%    struct, see below.
% - steps is a vector of initial step sizes, a struct, or a scalar which
%    multiplies StartPars.
% - mask is a boolean array telling which parameters are "active", that is
%    will be varied.
% - newPars is the returned array of parameters to try next, a row vector.
% - t is a structure containing the internal state.
% 
% The oldPars argument can alternatively be a struct.  In that case steps and
% mask, if given, can be structs with the same field names; e.g.
%   pars.a=3.14;  steps.a =.1;  mask.a=1;
%   for any missing field mask is taken to be zero.
%   steps may be a scalar for simplicity.
%   the returned newPar is a struct with the same fields as oldPars.
% 
% The test for convergence tests if the maximum and minimum y values in the
% simplex differ by less than epsilon.
% 
% --- Example of use:  minimize (p-q).^6 ---
% p=[2 2 1.5 2];  % inital guess
% q=[1 1 1 1];       % true values
% p=Simplex('init',p);
% for i=1:200
%     y=sum((p-q).^6);
%     p=Simplex(y);
%     if ~mod(i, 10)  % every 10 iterations print out the fitted value
%         p'
%     end;
% end;
% p=Simplex('centroid'); % obtain the final value.
% 
% ---Example of use with structs---
% par.kX=1;
% par.kY=2;
% msk.kX=1;
% msk.kY=0;  % don't vary par.kY (default is also zero if not assigned)
% stp.kX=.1;
% stp.kY=.2;
% x=(1:10)'; % set up the problem
% y=(1:10)';
% z=exp(-.5*x)+exp(-.5*y)+.5*randn(10,1);  % function to fit
% 
% par=Simplex('init',par,stp,msk);
% for i=1:100
%     y=exp(-par.kX*x)+exp(-par.kY*y);
%     err=(z-y)'*(z-y);
%     par=Simplex(err);
% end;
% finalPars=Simplex('centroid')  % returns a new struct with same fields
% 
% % Alternatively, the struct problem could be initialized these ways:
% par=Simplex('init',par);     % all parameters varied, steps = .1 * par values
% par=Simplex('init',par,.1);  % same as above
% par=Simplex('init',par,.1,msk);  % vary only where msk.field=1.
% 
% F. Sigworth, 15 March 2003
%  Mask added 11 March 2010
%  Pars as a struct added 20 March 2011
% Based on a Modula-2 implementation by S. H. Heinemann, 1987 which
% in turn was based on M. Caceci and W. Cacheris, Byte, p. 340, May 1984.

persistent PState;
% We make a temporary copy of our state variables for local use.
% This will allow the function to be interrupted without corrupting the
% state variables.
t=PState;   
% The structure elements are the following:
%   t.n  number of active elements in the parameter vector
%   t.active indices of active elements in the parameter vector
%   t.allpars copy of the StartPars vector
%   t.prow row vector of parameters presently being tested
%   t.simp the simplex matrix, (n+1) row vectors
%   t.vals the (n+1) element column vector of function values
%   t.high index of the highest t.vals element
%   t.low index of the lowest t.vals element
%   t.centr centroid row vector of the best n vertices
%   t.index counter for loops
%   t.state the state variable of the machine
%   t.structFieldnames  -- cell array of field names if we're using a
%   structure for parameter i/o

% default arguments for steps
defaultStep=.1;
zeroStep = 1e-3;

% Interpret the first argument.
if ischar(y)  % an option?
    switch lower(y)
        
        case 'init'  % Initialize the machine, with StartPars being the anchor vertex.
            if nargin<4
                mask=ones(1,numel(oldPars));
            end;
            if isstruct(oldPars)  % if the input is a structure, decode it.
                p=oldPars;
                t.structFieldnames=fieldnames(p);
                np=numel(t.structFieldnames);
                oldPars=zeros(np,1);
                for i=1:numel(t.structFieldnames)
                    field=t.structFieldnames{i};
                    val=p.(field);
                    oldPars(i)=val(1);
                end;
                if isstruct(steps)  % if Steps is given as a struct, use it
                    sSteps=steps;
                    steps=ones(np,1)*defaultStep;
                    for i=1:numel(t.structFieldnames)
                        field=t.structFieldnames{i};
                        if isfield(sSteps,field)
                            steps(i)=sSteps.(field);
                        end;
                    end;
                end;
                if isstruct(mask)  % if Mask is given as a struct, use it
                    sMask=mask;
                    mask=zeros(np,1);  % default is _not_ to vary
                    for i=1:numel(t.structFieldnames)
                        field=t.structFieldnames{i};
                        if isfield(sMask,field)
                            mask(i)=sMask.(field);
                        end;
                    end;
                end;
            else
                t.structFieldnames={};  % Mark that we're not using a struct.
            end;
                
            % Pick up the active variables (Mask entries > 0)
            t.allpars=oldPars(:)';  % row vector
            t.active=find(mask);
            t.n=numel(t.active);
            t.prow=reshape(oldPars(t.active),1,t.n);
            
            % Handle defaults for the step size.
            if nargin <3  % No step size given
                steps = defaultStep;
            end;
            if numel(steps)<numel(oldPars)  % not enough elements; assume a scalar
                steps=steps(1)*oldPars+zeroStep*(oldPars==0);
            end;

            % Pick up the relevant Steps elements
            steps=steps(t.active);
            % The simplex is (n+1) row vectors.
            t.simp=repmat(t.prow,t.n+1,1);
            for i=1:t.n
                t.simp(i+1,i)=t.simp(i+1,i)+steps(i);
            end;

            % vals is a column vector of function values
            t.vals=zeros(t.n+1,1);
            
            % Initialize the other variables
            t.index=1;
            t.state=1;
            t.high=0;
            t.low=0;
            t.centr=t.prow;
            newPars=t.allpars;
            newPars(t.active)=t.prow;
            PState=t;
            
        case 'centroid'  % Return the centroid of the present simplex
            active=(sum(t.simp)/(t.n+1))';
            newPars=t.allpars;
            newPars(t.active)=sum(t.simp)/(t.n+1);
        case 'converged'  % Do a convergence test on the vals array.
            err=max(t.vals)-min(t.vals);
            newPars=(t.state==3) && all(err < oldPars);  % return a boolean
            
        otherwise
            error('Simplex: unrecognized option');
    end; % switch
    
    
else  % y has a numeric value: this is running mode
    switch t.state
        
        case 1 % Start-up
            t.vals(t.index)=y;  % pick up the function value from last time.
            t.index=t.index+1;
            if t.index <= t.n+1  % continue to fill up the simplex
                t.prow=t.simp(t.index,:);
            else  % Simplex is full, make the first move
                t.state=3;
            end;
            
        case 3  % Test a new vertex
            i=t.high;
            if y < t.vals(i)  % The new vertex is better than some.
                t.simp(i,:)=t.prow;  % replace the worst one.
                t.vals(i)=y;
                if y < t.vals(t.low)  % The new vertex is better than the best,
                   t.prow=t.simp(i,:)+1.1*(t.simp(i,:)-t.centr); % so, expand in the new direction.
                   t.prevy=y;
                   t.state=4;
                else
                    t.state=3;
                end;
            else  % the new vertex is worse than the worst: contract.
                t.prow=0.5*(t.simp(t.high,:)+t.centr);
                t.state=5;
            end;
            
        case 4 % Test an expansion
            if y < t.prevy %t.vals(t.low)  % Accept the expansion
                t.simp(t.high,:)=t.prow;
                t.vals(t.high)=y;
            end;
            t.state=3;
            
        case 5 % Test a contraction
            if y<t.vals(t.high) % Accept the contraction
                t.simp(t.high,:)=t.prow;
                t.vals(t.high)=y;
                t.state=3;
            else %  contract the whole simplex toward the best vertex.
                t.index=1;
                t.simp(1,:)=0.5*(t.simp(1,:)+t.simp(t.low,:));
                prow=t.simp(1,:);
                t.state=6;
            end;
            
        case 6
            t.vals(t.index)=y;  % pick up the function value.
            t.index=t.index+1;
            i=t.index;
            if i <= t.n+1
                t.simp(i,:)=0.5*(t.simp(i,:)+t.simp(t.low,:));
                t.prow=t.simp(i,:);
                t.state=6;  % continue evaluating the vertices
            else
                t.state=3;  % 
            end;
    end; % switch
        
    if t.state==3  % Normal exit mode: sort the vertices and try a reflection.
        % assign min and max
        [z,ind]=sort(t.vals);
        t.low=ind(1);
        t.high=ind(t.n+1);
        % find the excluded centroid
        t.centr=(sum(t.simp)-t.simp(t.high,:))/(t.n);
        % reflect about the centroid from the highest vertex
        t.prow=2*t.centr-t.simp(t.high,:);
    end;
    
    
    % Copy the output and persistent variables.
    newPars=t.allpars;
    newPars(t.active)=t.prow;
    PState=t;
end;

% If we're using structs, return one.
np=numel(t.structFieldnames);
if np>0  % need to return a struct
    vals=newPars;
    newPars=struct;
    for i=1:np
        newPars.(t.structFieldnames{i})=vals(i);
    end;
end

