function [infoStruct] = UGM_makeInfoStruct(X,Xedge,edgeStruct,nStates,ising,tied,useMex,maxIter)

if nargin < 8
    maxIter = 50;
end

[nInstances nNodeFeatures nNodes] = size(X);
[nEdgeFeatures] = size(Xedge,2);
maxState = max(nStates);
edgeEnds = edgeStruct.edgeEnds;
nEdges = size(edgeEnds,1);

if isscalar(tied)
    tieNodes = tied;
    tieEdges = tied;
else
    tieNodes = tied(1);
    tieEdges = tied(2);
end

% Check that nStates is scalar for completely tied models
if (tieNodes == 1 || tieEdges == 1) && ~isscalar(nStates)
    if max(nStates) ~= min(nStates)
        assert(1==0,'Untied number of states not supported for tied models');
    end
end

if tieNodes
    wSize = [nNodeFeatures maxState-1];
else
    wSize = [nNodeFeatures maxState-1 nNodes];
end

if tieEdges
    if ising
        vSize = [nEdgeFeatures 1];
    else
        vSize = [nEdgeFeatures maxState^2-1];
    end
else
    if ising
        vSize = [nEdgeFeatures 1 nEdges];
    else
        vSize = [nEdgeFeatures maxState^2-1 nEdges];
    end
end
   

% Make Linear Parameter Index
if isscalar(nStates)
   wLinInd = find(ones(wSize));
   vLinInd = find(ones(vSize));
   nStates = repmat(nStates,[nNodes 1]);
else

   % Remove variables for states that do not occur

   % Irrelevant Node Weights
   wInd = ones(wSize);
   for n = 1:nNodes
      for s = nStates(n):maxState-1
         wInd(:,s,n) = 0;
      end
   end
   wLinInd = find(wInd);

   % Irrelevant Edge Weights
   if ising
      vLinInd = find(ones(vSize));
   else
      vInd = ones(vSize);
      for e = 1:nEdges
         n1 = edgeEnds(e,1);
         n2 = edgeEnds(e,2);

         for s1 = 1:maxState
            for s2 = 1:maxState
               if s1 > nStates(n1) || s2 > nStates(n2)
                  s = sub2ind([maxState maxState],s1,s2);
                  if s < maxState^2
                     vInd(:,s,e) = 0;
                  end
               end
            end
         end

         if nStates(n1) < maxState || nStates(n2) < maxState
            s = sub2ind([maxState maxState],nStates(n1),nStates(n2));
            vInd(:,s,e) = 0;
         end
      end
      vLinInd = find(vInd);
   end
   
   if size(nStates,1) == 1
       nStates = nStates';
   end
end


infoStruct.wSize = wSize;
infoStruct.vSize = vSize;
infoStruct.wLinInd = wLinInd;
infoStruct.vLinInd = vLinInd;
infoStruct.nStates = nStates;
infoStruct.maxIter = maxIter;
infoStruct.ising = ising;
infoStruct.tieNodes = tieNodes;
infoStruct.tieEdges = tieEdges;
infoStruct.useMex = useMex;