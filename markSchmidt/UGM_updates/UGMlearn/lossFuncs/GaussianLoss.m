function [nll,g,H] = GaussianLoss(w,X,y)
% w(feature,1)
% X(instance,feature)
% y(instance,1)
    
XX = X.'*X;
XXw = XX*w;
Xy = X.'*y;

nll = w.'*XXw - 2*w.'*Xy + y.'*y;

if nargout > 1
    g = 2*XXw - 2*Xy;
end

if nargout > 2
    H = 2*XX;
end