function build_graphCutMex
% build_graphCutMex builds package graphCutMex
%
% Anton Osokin (firstname.lastname@gmail.com),  19.05.2013

maxFlowPath = 'maxflow-v3.03.src';

mexCmd = ['mex graphCutMex.cpp -output graphCutMex -largeArrayDims ', '-I', maxFlowPath];
eval(mexCmd);
