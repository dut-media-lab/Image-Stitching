%-------------------------------------------------------------------------
% Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
%-------------------------------------------------------------------------
% The demo code in this package implements the guided-sampling method for
% multi-structure robust fitting proposed in:
% T.-J. Chin, J. Yu and D. Suter
% Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
% In Proc. European Conf. on Computer Vision, Crete, Greece, 2010.
%
% Copyright (c) 2010 Tat-Jun Chin and Jin Yu
% School of Computer Science, The University of Adelaide, South Australia
% http://www.cs.adelaide.edu.au/~{tjchin,jinyu}
%
% The program is free for non-commercial academic use. Any commercial use
% is strictly prohibited without the authors' consent. Please acknowledge
% the authors by citing the above paper in any academic publications that
% have made use of this package or part of it.
%
% If you encounter any problems or questions please email to 
% tjchin@cs.adelaide.edu.au.
% 
% This program makes use of Peter Kovesi and Andrew Zisserman's MATLAB
% functions for multi-view geometry
% (http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/
%  http://www.robots.ox.ac.uk/~vgg/hzbook/code/).


function [ par, res, inx, tim, err ] = multigsSampling(lim,data,M,blksiz)
% input: 
% lim (1x1) = Maximum CPU seconds allowed.
% data (dxn) = Input data of dimensionality d.
% M (1x1) = Maximum number of hypotheses to be generated.
% blksize (1x1) = Block size of Multi-GS.
%
% output:
% par (dxM) = Parameters of the putative models.
% res (nxM) = Residuals as measured to the putative models.
% inx (pxM) = Indices of p-subsets
% tim (1xM) = CPU time for generating each model.


%---------------------------
% Model specific parameters.
%---------------------------
% Now defined globally at file main.m.
global fitfn resfn degenfn psize numpar

%------------------------
% Check other parameters.
%------------------------
degenmax = 10;  % Max number of inner loop to avoid degeneracy.

if (mod(M,blksiz)~=0)
    error('Bad block size!');
end

%-----------------
% Prepare storage.
%-----------------
n = size(data,2);
par = zeros(numpar,M);
res = zeros(n,M);
inx = zeros(psize,M);
tim = zeros(1,M);
err = 0;
%-----------------
% Guided sampling.
%-----------------
fprintf('> Multi-GS (RANSAC) sampling for %.2f seconds...',lim);
t0 = cputime;
for m=1:M

    degencnt = 0;
    isdegen = 1;
    while (isdegen==1)&&(degencnt<=degenmax)        
        % Increment degeneracy count.
        degencnt = degencnt + 1;
        
        if m<=blksiz            
            % Uniform sampling for the first block.
            [ pinx ] = randsample(n,psize);
        else            
            % Weighted sampling
            [ pinx ] = weightedSampling(n,psize,resinx,win);
        end
        psub = data(:,pinx);        
        
        % Check for degeneracy.
        isdegen = feval(degenfn,psub);
    end
    
    if (isdegen==1)
        fprintf('\n  ');
        warning('Cannot find a valid p-subset!');
        err = 1;
        return;
    end    
    
    % Fit the model on the p-subset
    st = feval(fitfn,psub);
    
    % Compute residuals.
    ds = feval(resfn,st,data);
       
    % Store.
    par(:,m) = st;
    res(:,m) = ds;
    inx(:,m) = pinx;
    tim(1,m) = cputime-t0;
    
    if tim(1,m)>=lim
        par = par(:,1:m);
        res = res(:,1:m);
        inx = inx(:,1:m);
        tim = tim(:,1:m);
        break;
    end

    % Update indices of sorted hypotheses.
    if (m>=blksiz)&&(mod(m,blksiz)==0)
        
        % Intersection width.
        win = round(0.1*m);
        
        % Compute indices of sorted hypotheses.     
        [ ~, resinx ] = sort(res(:,1:m),2);
    end

end
fprintf('done (%fs)\n',tim(end));

end

%-------------------
% Weighted sampling.
%-------------------
function [ pinx ] = weightedSampling(n,psize,resinx,win)
% n (1x1)      = Size of data.
% psize (1x1)  = Size of p-subset.
% resinx (nxm) = Indices of sorted hypotheses for each datum.
% win (1x1)    = Intersection width.

% Storage.
pinx = zeros(psize,1);

% First index is found by uniform sampling.
seedinx = randsample(n,1);
pinx(1,:) = seedinx;

% Weighted sampling.
w = ones(1,n);
for i=2:psize
    new_w = computeIntersection(resinx(seedinx,:)',resinx',win);
    
    % Ensures sampling without replacement.
    new_w(seedinx) = 0;
    
    % Aggregate weights
    w = w.*new_w;
    
    if sum(w) > 0
        othinx = randsample(n,1,true,w);
    else
        % If weights are zero, do random sampling.
        othinx = randsample(n,1);
    end
    pinx(i,:) = othinx;
    seedinx = othinx;
end

end
