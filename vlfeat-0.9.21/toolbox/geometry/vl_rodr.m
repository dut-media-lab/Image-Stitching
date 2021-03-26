% VL_RODR  Rodrigues' formula
%   R = VL_RODR(OM) where OM a 3-dimensional column vector computes the
%   Rodrigues' formula of OM, returning the rotation matrix R =
%   expm(vl_hat(OM)).
%
%   [R,DR] = VL_RODR(OM) computes also the derivative of the Rodrigues
%   formula. In matrix notation this is the expression
%
%           d(vec expm(vl_hat(OM)) )
%     dR = ----------------------.
%                  d om^T
%
%   [R,DR]=VL_RODR(OM) when OM is a 3xK matrix repeats the operation for
%   each column (or equivalently matrix with 3*K elements). In this
%   case R and DR are arrays with K slices, one per rotation.
%
%   See also: VL_IRODR(), VL_HELP().

% Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
% All rights reserved.
%
% This file is part of the VLFeat library and is made available under
% the terms of the BSD license (see the COPYING file).
