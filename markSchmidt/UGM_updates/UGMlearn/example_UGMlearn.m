function [] = example_UGMlearn(display)

%% Parameters of Experiment
nTrain = 500; % number of examples to use for training
nTest = 1000; % number of examples to use for test
nFeatures = 10; % number of features for each node
nNodes = 10; % number of nodes
nStates = 2; % number of states that each node can take
edgeProb = .5; % probability of each edge being included in the graph
edgeType = 1; % set to 0 to make the edge features normally distributed
ising = 0; % set to 1 to use ising potentials
trainType = 'pseudo'; % set to 'pseudo' for pseudo-likelihood, 'loopy' for loopy belief propagation, 'exact' for 'exact inference
testType = 'exact'; % set to 'loopy' to test with loopy belief propagation, 'exact' for exact inference
structureSeed = 0; % change this to generate different structures
trainSeed = 0; % vary seed from 0:9 to get paper results
useMex = 1; % use mex files in UGM to speed things up

% Regularization Parameters
lambdaNode = 10;
lambdaEdge = 10;
% In the paper, we picked these two values by two-fold
% cross-validation, testing the values 2.^[7:-1:-5] for each model
% In the paper, we also used warm-starting to speed up the optimization
%   for this sequence of values

%% Generate data
fprintf('Generating Synthetic CRF Data...\n');
rand('state',structureSeed);
randn('state',structureSeed);
nInstances = nTrain+nTest;
tied = 0;
[y,adjTrue,X] = UGM_generate(nInstances,nFeatures,nNodes,edgeProb,nStates,ising,tied,useMex,edgeType);

if display
    f = 1;
    figure(f);clf;hold on;
    drawGraph(adjTrue);
    title('True Graph');
end

%% Form Training/testing premutation
rand('state',trainSeed);
randn('state',trainSeed);
perm = randperm(nInstances);
trainNdx = perm(1:nTrain);
testNdx = perm(nTrain+1:end);

if isscalar(nStates)
    nStates = repmat(nStates,[nNodes 1]);
end

%% Fixed Structure Methods

edgePenaltyType = 'L2';
subDisplay = display;

% Empty Structure
type = 'Fixed: Empty';
adjInit = zeros(nNodes);
example_UGMlearnSub;

% Chain Structure
type = 'Fixed: Chain';
adjInit = fixed_Chain(nNodes);
example_UGMlearnSub;

% Full Structure
type = 'Fixed: Full';
adjInit = fixed_Full(nNodes);
example_UGMlearnSub;

% True Structure
type = 'Fixed: True';
adjInit = adjTrue;
example_UGMlearnSub;

%% Generative Non-L1 Methods

edgePenaltyType = 'L2';

% Optimal Tree
type = 'Generative Non-L1: Tree';
[junk1,junk2,junk3,adjTree] = chowliu(y',nStates,0);
adjInit = adjTree+adjTree';
example_UGMlearnSub

% Greedy DAG-Search
type = 'Generative Non-L1: DAG';
if all(nStates == 2) && ising % binary ising
    CPD = 1;
    yDS = sign(y-1.5);
else
    for s = 1:nNodes
        CPD.type(s) = 'M';
        CPD.nStates(s) = nStates(s);
    end
    yDS = y;
end
adjDS = DAGsearch(yDS,inf,0,log(nTrain)/2,CPD,zeros(size(yDS)),ones(nNodes),display,zeros(nNodes));
for i = 1:nNodes
    parents = find(adjDS(:,i));
    for p1 = parents
        for p2 = parents
            if p1 ~= p2
                adjDS(p1,p2) = 1; % Moralize
            end
        end
    end
end
adjInit = sign(adjDS+adjDS'); % Make symmetric
example_UGMlearnSub;

%% Generative L1 Methods

% Regular L1
type = 'Generative L1: L1-L1';
edgePenaltyType = 'L1'; % Train with L1 regularization on the edges
subDisplay = -1;
Xold = X;
X = zeros(nInstances,0,nNodes); % Train while ignoring features
adjInit = fixed_Full(nNodes);
example_UGMlearnSub

adjInit = adjFinal;
X = Xold;
edgePenaltyType = 'L2';
subDisplay = display;
example_UGMlearnSub

% Group-L1, using 2-norm of groups
type = 'Generative L1: L1-L2';
edgePenaltyType = 'L1-L2'; % Train with L1 regularization on the edges
subDisplay = -1;
Xold = X;
X = zeros(nInstances,0,nNodes); % Train while ignoring features
adjInit = fixed_Full(nNodes);
example_UGMlearnSub

adjInit = adjFinal;
X = Xold;
edgePenaltyType = 'L2';
subDisplay = display;
example_UGMlearnSub

% Group-L1, using inf-norm of groups
type = 'Generative L1: L1-Linf';
edgePenaltyType = 'L1-Linf'; % Train with L1 regularization on the edges
subDisplay = -1;
Xold = X;
X = zeros(nInstances,0,nNodes); % Train while ignoring features
adjInit = fixed_Full(nNodes);
example_UGMlearnSub

adjInit = adjFinal;
X = Xold;
edgePenaltyType = 'L2';
subDisplay = display;
example_UGMlearnSub

%% Discriminative L1 Methods

if nFeatures > 0
    type = 'Discriminative L1: L1-L1';
    edgePenaltyType = 'L1'; % Train with L1 regularization on the edges
    adjInit = fixed_Full(nNodes);
    example_UGMlearnSub

    type = 'Discriminative L1: L1-L2';
    edgePenaltyType = 'L1-L2'; % Train with L1 regularization on the edges
    adjInit = fixed_Full(nNodes);
    example_UGMlearnSub

    type = 'Discriminative L1: L1-Linf';
    edgePenaltyType = 'L1-Linf'; % Train with L1 regularization on the edges
    adjInit = fixed_Full(nNodes);
end
example_UGMlearnSub