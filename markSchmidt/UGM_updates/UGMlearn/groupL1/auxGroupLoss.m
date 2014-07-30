function [f,g] = groupLinfLoss(w,groups,lambda,funObj)
p = length(groups);
nGroups = length(p+1:length(w));
[f,g] = funObj(w(1:p));

f = f + sum(lambda.*w(p+1:end));
g = [g;lambda.*ones(nGroups,1)];