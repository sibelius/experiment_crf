#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Variables */
    int i,n,s,e,n1,n2,Vind,maxIter,
    nNodes, nEdges, maxState, sizeS[2],
    *edgeEnds, *nStates, *V, *E,*y,*S;
    
    double *pot,z,U,u,
    *nodePot, *edgePot;
    
   /* Input */
    
    nodePot = mxGetPr(prhs[0]);
    edgePot = mxGetPr(prhs[1]);
    edgeEnds = mxGetPr(prhs[2]);
    nStates = mxGetPr(prhs[3]);
    V = mxGetPr(prhs[4]);
    E = mxGetPr(prhs[5]);
    maxIter = mxGetPr(prhs[6])[0];
    
   /* Compute Sizes */
    
    nNodes = mxGetDimensions(prhs[0])[0];
    maxState = mxGetDimensions(prhs[0])[1];
    nEdges = mxGetDimensions(prhs[2])[0];
    decrementEdgeEnds(edgeEnds,nEdges);
    decrementVector(V,nNodes+1);
    decrementVector(E,nEdges*2);

    
   /* Output */
    pot = mxCalloc(maxState,sizeof(double));
    y = mxCalloc(nNodes,sizeof(int));
    sizeS[0] = nNodes;
    sizeS[1] = maxIter;
    plhs[0] = mxCreateNumericArray(2,sizeS,mxINT32_CLASS,mxREAL);
    S = mxGetPr(plhs[0]);
    
   /* Initialize to States with highest node potentials*/
    for(n = 0; n < nNodes; n++)
    {
        u = -1;
        U = 0;
        for(s = 0; s < nStates[n]; s++)
        {
            if(nodePot[n+nNodes*s] > u)
            {
                u = nodePot[n+nNodes*s];
                U = s;
            }
        }
        y[n] = U;
    }
    
    for(i = 0; i < maxIter; i++)
    {
        for(n = 0; n < nNodes; n++)
        {
            /* Compute Node Potential */
            for(s = 0; s < nStates[n]; s++)
            {
                pot[s] = nodePot[n + nNodes*s];
            }
            
           /* Iterate over Neighbors */
            for(Vind = V[n]; Vind < V[n+1]; Vind++)
            {
                e = E[Vind];
                n1 = edgeEnds[e];
                n2 = edgeEnds[e+nEdges];
                 
                /* Multiply Edge Potentials */
                if(n == n1)
                {
                   for(s = 0; s < nStates[n]; s++)
                   {
                        pot[s] *= edgePot[s+maxState*(y[n2] + maxState*e)];
                   }
                    
                }
                else
                {
                    for(s = 0; s < nStates[n]; s++)
                    {
                        pot[s] *= edgePot[y[n1]+maxState*(s + maxState*e)];
                    }
                }
                
            }
            
            /* Normalize */
            z = 0;
            for(s = 0; s < nStates[n]; s++)
                z = z + pot[s];
            for(s = 0; s < nStates[n]; s++)
                pot[s] /= z;
            
            /* Display */
            for(s = 0; s < nStates[n]; s++)
            
            /* Sample Discrete State */
            U = rand()/((double)RAND_MAX + 1);
            u = 0;
            for(s = 0; s < nStates[n]; s++)
            {
                u += pot[s];
                if(u > U)
                {
                    break;
                }
            }
            y[n] = s;
        }
        /* Record Sample */
        for(n = 0; n < nNodes; n++)
        {
            S[n + nNodes*i] = y[n]+1;
        }
    }
    
   /* Free memory */
    mxFree(pot);
    mxFree(y);
}
