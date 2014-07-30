function [samples] = UGM_Sample_Gibbs(nodePot,edgePot,edgeStruct,infoStruct)
% Single Site Gibbs Sampling

if infoStruct.useMex
    samples = double(UGM_Sample_GibbsC(nodePot,edgePot,int32(edgeStruct.edgeEnds),int32(infoStruct.nStates),int32(edgeStruct.V),int32(edgeStruct.E),infoStruct.maxIter));
else
    samples = Sample_Gibbs(nodePot,edgePot,edgeStruct,infoStruct);
end

end

function [samples] = Sample_Gibbs(nodePot,edgePot,edgeStruct,infoStruct)
[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = infoStruct.nStates;
maxIter = infoStruct.maxIter;

% Initialize
[junk y] = max(nodePot,[],2);

samples = zeros(nNodes,0);

for i = 1:maxIter

    for n = 1:nNodes

        % Compute Node Potential
        pot = nodePot(n,1:nStates(n));

        % Find Neighbors
        edges = E(V(n):V(n+1)-1);

        % Multiply Edge Potentials
        for e = edges(:)'
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);

            if n == edgeEnds(e,1)
                ep = edgePot(1:nStates(n1),y(n2),e)';
            else
                ep = edgePot(y(n1),1:nStates(n2),e);
            end
            pot = pot .* ep;
        end

        % Sample State;
        y(n) = sampleDiscrete(pot./sum(pot));
    end
    samples(:,i) = y;
end
end