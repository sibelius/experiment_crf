% Plots Regularization Path for L1L1/L1Linf/L1L2
clear all
close all

nInstances = 100;
p = 15;
lambdaValues = 4000:-100:0;

% Divide data into 3 groups
groups = zeros(p,1);
groups(1:3:p) = [1:p/3]';
groups(2:3:p) = [1:p/3]';
groups(3:3:p) = [1:p/3]';
nGroups = length(unique(groups(groups>0)));


%% Make Data
X = randn(nInstances,p);
w = randn(p,1).*groups; % Higher groups have higher weight multiplier
y = X*w + randn(nInstances,1);

%% Minimize Gaussian Loss function with the different regularizers
funObj_sub = @(w)GaussianLoss(w,X,y);

%% L1L1

i = 1;
w_init = zeros(p,1);
clear w
for lambda = lambdaValues
    lambdaVect = lambda*ones(p,1);
    
    funObj = @(w)nonNegGrad(w,lambdaVect,funObj_sub);
    wPosNeg = minConF_BC(funObj,[w_init;-w_init],zeros(2*p,1),inf(2*p,1));
    w(:,i) = wPosNeg(1:p)-wPosNeg(p+1:end);
    
    w_init = w(:,i);
    i = i + 1;
end

% Plot
figure(1);clf; hold on;
hTitle  = title ('(p=1) path');
plotGroupL1
print -dpng 'p1Path.png'
fprintf('(paused)\n');
pause;

%% L1L2 (parameterized in terms of lambda)

i = 1;
w_init = zeros(p,1);
clear w
for lambda = lambdaValues

    % Make Initial Value
    w_init = [w_init;zeros(nGroups,1)];
    for g = 1:nGroups
        w_init(p+g) = norm(w_init(groups==g));
    end

    % Make Objective and Projection Function
    funObj = @(w)auxGroupLoss(w,groups,lambda,funObj_sub);
    funProj = @(w)auxGroupL2Proj(w,groups);

    % Solve
    wAlpha = minConF_SPG(funObj,w_init,funProj);

    w(:,i) = wAlpha(1:p);
    w_init = w(:,i);
    i = i + 1;
end

% Plot
figure(2);clf; hold on;
hTitle  = title ('(p=2) path');
plotGroupL1
print -dpng 'p2Path.png'
fprintf('(paused)\n');
pause;

%% L1Linf
w_init = zeros(p,1);
clear w

i = 1;
for lambda = lambdaValues

    % Make Initial Value
    w_init = [w_init;zeros(nGroups,1)];
    for g = 1:nGroups
        w_init(p+g) = max(abs(w_init(groups==g)));
    end

    % Make Objective and Projection Function
    funObj_sub = @(w)GaussianLoss(w,X,y);
    funObj = @(w)auxGroupLoss(w,groups,lambda,funObj_sub);
    funProj = @(w)auxGroupLinfProj(w,groups);

    % Solve
    wAlpha = minConF_SPG(funObj,w_init,funProj);

    w(:,i) = wAlpha(1:p);
    w_init = w(:,i);
    i = i + 1;
end

% Plot
figure(3);clf; hold on;
hTitle  = title ('(p=inf) path');
plotGroupL1
print -dpng 'pinfPath.png'
fprintf('(paused)\n');
pause;

