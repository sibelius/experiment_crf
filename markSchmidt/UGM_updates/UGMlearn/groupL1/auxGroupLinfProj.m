function w = groupLinfProj(w,groups)
p = length(groups);
alpha = w(p+1:end);
w = w(1:p);
u = unique(groups(groups>0));

for i = 1:length(u)
    % Sorting approach
    [w(groups==u(i)) alpha(i)] = projectAuxSort(w(groups==u(i)),alpha(i));
end
w = [w;alpha];

end

%% Function to solve the projection for a single group
function [w,alpha] = projectAuxSort(w,alpha)
if ~all(abs(w) <= alpha)
    sorted = [sort(abs(w),'descend');0];
    s = 0;
    for k = 1:length(sorted)

        % Compute Projection with k largest elements
        s = s + sorted(k);
        projPoint = (s+alpha)/(k+1);
       
        if projPoint > 0 && projPoint > sorted(k+1)
            w(abs(w) >= sorted(k)) = sign(w(abs(w) >= sorted(k)))*projPoint;
            alpha = projPoint;
            break;
        end

        if k == length(sorted)
            % alpha is too negative, optimal answer is 0
            w = zeros(size(w));
            alpha = 0;
        end
    end
end
end