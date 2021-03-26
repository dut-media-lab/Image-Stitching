% graphCutMex - Matlab wrapper to the implementation of min-cut algorithm by Yuri Boykov and Vladimir Kolmogorov:
% 	http://pub.ist.ac.at/~vnk/software/maxflow-v3.03.src.zip
% This version can automatically perform reparametrization on all submodular edges.
% 
% The algorithm is described in:
% Yuri Boykov and Vladimir Kolmogorov, 'An Experimental Comparison of Min-Cut/Max-Flow Algorithms for
% Energy Minimization in Vision', IEEE Transactions on Pattern Analysis and Machine Intelligence, vol.
% 26, no. 9, pp. 1124-1137, Sept. 2004.
% 
% Usage:
% [cut] = graphCutMex(termWeights, edgeWeights);
% [cut, labels] = graphCutMex(termWeights, edgeWeights);
% 
% Inputs:
% termWeights	-	the edges connecting the source and the sink with the regular nodes (array of type double, size : [numNodes, 2])
% 				termWeights(i, 1) is the weight of the edge connecting the source with node #i
% 				termWeights(i, 2) is the weight of the edge connecting node #i with the sink
% 				numNodes is determined from the size of termWeights.
% edgeWeights	-	the edges connecting regular nodes with each other (array of type double, array size [numEdges, 4])
% 				edgeWeights(i, 3) connects node #edgeWeights(i, 1) to node #edgeWeights(i, 2)
% 				edgeWeights(i, 4) connects node #edgeWeights(i, 2) to node #edgeWeights(i, 1)
%				The only requirement on edge weights is submodularity: edgeWeights(i, 3) + edgeWeights(i, 4) >= 0
%
% Outputs:
% cut           -	the minimum cut value (type double)
% labels		-	a vector of length numNodes, where labels(i) is 0 or 1 if node #i belongs to S (source) or T (sink) respectively.
% 
% To build the code in Matlab choose reasonable compiler and run build_graphCutMex.m
% Run example_graphCutMex.m to test the code
% 
% 	Anton Osokin (firstname.lastname@gmail.com),  19.05.2013
