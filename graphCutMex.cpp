
#include "graphCutMex.h"
#include "mex.h"

#include <limits>
#include <cmath>

#define INFTY INT_MAX

//define types
typedef double EnergyType;
mxClassID MATLAB_ENERGYTERM_TYPE = mxDOUBLE_CLASS;

typedef double EnergyTermType;
mxClassID MATLAB_ENERGY_TYPE = mxDOUBLE_CLASS;

typedef double LabelType;
mxClassID MATLAB_LABEL_TYPE = mxDOUBLE_CLASS;
/*
typedef int EnergyType;
mxClassID MATLAB_ENERGYTERM_TYPE = mxINT32_CLASS;

typedef int EnergyTermType;
mxClassID MATLAB_ENERGY_TYPE = mxINT32_CLASS;

typedef int LabelType;
mxClassID MATLAB_LABEL_TYPE = mxINT32_CLASS;
*/

typedef Graph<EnergyTermType,EnergyTermType,EnergyType> GraphType; 

double round(double a);
int isInteger(double a);

#define MATLAB_ASSERT(expr,msg) if (!(expr)) { mexErrMsgTxt(msg);}

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef int mwSize;
typedef int mwIndex;
#endif



void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
	MATLAB_ASSERT( nrhs == 2, "graphCutMex: Wrong number of input parameters: expected 2");
    MATLAB_ASSERT( nlhs <= 2, "graphCutMex: Too many output arguments: expected 2 or less");
	
	//Fix input parameter order:
	const mxArray *uInPtr = (nrhs >= 1) ? prhs[0] : NULL; //unary
	const mxArray *pInPtr = (nrhs >= 2) ? prhs[1] : NULL; //pairwise
	
	//Fix output parameter order:
	mxArray **cOutPtr = (nlhs >= 1) ? &plhs[0] : NULL; //cut
	mxArray **lOutPtr = (nlhs >= 2) ? &plhs[1] : NULL; //labels

	 //node number
	int numNodes;
    
	// get unary potentials
	MATLAB_ASSERT(mxGetNumberOfDimensions(uInPtr) == 2, "graphCutMex: The first paramater is not 2-dimensional");
	MATLAB_ASSERT(mxGetClassID(uInPtr) == MATLAB_ENERGYTERM_TYPE, "graphCutMex: Unary potentials are of wrong type");
	MATLAB_ASSERT(mxGetPi(uInPtr) == NULL, "graphCutMex: Unary potentials should not be complex");
	
	numNodes = mxGetM(uInPtr);
	
	MATLAB_ASSERT(numNodes >= 1, "graphCutMex: The number of nodes is not positive");
	MATLAB_ASSERT(mxGetN(uInPtr) == 2, "graphCutMex: The first paramater is not of size #nodes x 2");
	
	EnergyTermType* termW = (EnergyTermType*)mxGetData(uInPtr);

	//get pairwise potentials
	MATLAB_ASSERT(mxGetNumberOfDimensions(pInPtr) == 2, "graphCutMex: The second paramater is not 2-dimensional");
	
	mwSize numEdges = mxGetM(pInPtr);

	MATLAB_ASSERT( mxGetN(pInPtr) == 4, "graphCutMex: The second paramater is not of size #edges x 4");
	MATLAB_ASSERT(mxGetClassID(pInPtr) == MATLAB_ENERGYTERM_TYPE, "graphCutMex: Pairwise potentials are of wrong type");

	EnergyTermType* edges = (EnergyTermType*)mxGetData(pInPtr);
	for(int i = 0; i < numEdges; i++)
	{
		MATLAB_ASSERT(1 <= round(edges[i]) && round(edges[i]) <= numNodes, "graphCutMex: error in pairwise terms array: wrong vertex index");
		MATLAB_ASSERT(isInteger(edges[i]), "graphCutMex: error in pairwise terms array: wrong vertex index");
		MATLAB_ASSERT(1 <= round(edges[i + numEdges]) && round(edges[i + numEdges]) <= numNodes, "graphCutMex: error in pairwise terms array: wrong vertex index");
		MATLAB_ASSERT(isInteger(edges[i + numEdges]), "graphCutMex: error in pairwise terms array: wrong vertex index");
	    MATLAB_ASSERT(edges[i + 2 * numEdges] + edges[i + 3 * numEdges] >= 0, "graphCutMex: error in pairwise terms array: nonsubmodular edge");
	}


	// start computing
	if (nlhs == 0){
		return;
	}

	//prepare graph
	GraphType *g = new GraphType( numNodes, numEdges); 
	
	for(int i = 0; i < numNodes; i++)
	{
		g -> add_node(); 
		g -> add_tweights( i, termW[i], termW[numNodes + i]); 
	}
	
	for(int i = 0; i < numEdges; i++)
		if(edges[i] < 1 || edges[i] > numNodes || edges[numEdges + i] < 1 || edges[numEdges + i] > numNodes || edges[i] == edges[numEdges + i] || !isInteger(edges[i]) || !isInteger(edges[numEdges + i])){
			mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Some edge has invalid vertex numbers and therefore it is ignored");
		}
		else
			if(edges[2 * numEdges + i] + edges[3 * numEdges + i] < 0){
				mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Some edge is non-submodular and therefore it is ignored");
			}
			else
			{
				if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] >= 0)
					g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), edges[2 * numEdges + i], edges[3 * numEdges + i]);
				else
					if (edges[2 * numEdges + i] <= 0 && edges[3 * numEdges + i] >= 0)
					{
						g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i] + edges[2 * numEdges + i]);
						g -> add_tweights((GraphType::node_id)round(edges[i] - 1), 0, edges[2 * numEdges + i]); 
						g -> add_tweights((GraphType::node_id)round(edges[numEdges + i] - 1),0 , -edges[2 * numEdges + i]); 
					}
					else
						if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] <= 0)
						{
							g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), edges[3 * numEdges + i] + edges[2 * numEdges + i], 0);
							g -> add_tweights((GraphType::node_id)round(edges[i] - 1),0 , -edges[3 * numEdges + i]); 
							g -> add_tweights((GraphType::node_id)round(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i]); 
						}
						else
							mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Something strange with an edge and therefore it is ignored");
			}

	//compute flow
	EnergyType flow = g -> maxflow();

	//output minimum value
	if (cOutPtr != NULL){
		*cOutPtr = mxCreateNumericMatrix(1, 1, MATLAB_ENERGY_TYPE, mxREAL);
		*(EnergyType*)mxGetData(*cOutPtr) = (EnergyType)flow;
	}

	//output minimum cut
	if (lOutPtr != NULL){
		*lOutPtr = mxCreateNumericMatrix(numNodes, 1, MATLAB_LABEL_TYPE, mxREAL);
		LabelType* segment = (LabelType*)mxGetData(*lOutPtr);
		for(int i = 0; i < numNodes; i++)
			segment[i] = g -> what_segment(i);
	}
    
    delete g;
}

double round(double a)
{
	return floor(a + 0.5);
}


int isInteger(double a)
{
	return (abs(a - round(a)) < 1e-6);
}
