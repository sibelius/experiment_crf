natops = load('NATOPS24.mat');

X = natops.seqs;
Y = natops.labels;

% Put into UGM format
fprintf('Putting NATOPS data into UGM format\n');
nFeatures = size(X{1},1);
nStates = max(cell2mat(Y));

ising = 0;
tied = 1;
paramLastState = 1;

nExamples = length(X);
examples = cell(nExamples,1);

for i = 1:nExamples,
    nNodes = size(X{i},2);
    examples{i}.edgeStruct = UGM_makeEdgeStruct(chainAdjMatrix(nNodes),nStates);
    examples{i}.Y = int32(Y{i});
    for n = 1:nNodes,
        examples{i}.Xnode(1,:,n) = [1 X{i}(:,n)'];
    end
    examples{i}.Xedge = ones(1,1,nNodes-1);
    [examples{i}.nodeMap examples{i}.edgeMap w] = UGM_makeCRFmaps(examples{i}.Xnode,examples{i}.Xedge,examples{i}.edgeStruct, ...
        ising,tied,paramLastState);
end

% Split training/testing data
fprintf('Splitting into training and testing data\n');
[trainfolds, testfolds] = Kfold(nExamples, 2, 0);
trainExamples = examples(trainfolds{1});
testExamples = examples(testfolds{1});

% Train L2-regularized CRF
fprintf('Finding maximum likelihood parameters\n');
UGM_CRFcell_NLL(w,trainExamples,@UGM_Infer_Chain);
funObj = @(w)UGM_CRFcell_NLL(w,trainExamples,@UGM_Infer_Chain);
lambda = ones(size(w));
penalizedFunObj = @(w)penalizedL2(w,funObj,lambda);
options.Display = 'full';
w = minFunc(penalizedFunObj,w,options);

% Compute test error
testErrs = 0;
Z = 0;
for i = 1:length(testExamples)
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,testExamples{i}.Xnode,testExamples{i}.Xedge,testExamples{i}.nodeMap,testExamples{i}.edgeMap,  testExamples{i}.edgeStruct);
    yMAP = UGM_Decode_Chain(nodePot,edgePot,testExamples{i}.edgeStruct);

    testErrs = testErrs + sum(yMAP'~=testExamples{i}.Y);
    Z = Z + length(testExamples{i}.Y);
end
testErrorRate = testErrs/Z

