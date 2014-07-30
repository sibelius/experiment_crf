function w = groupLinfProj(w,groups)
p = length(groups);
alpha = w(p+1:end);
w = w(1:p);
u = unique(groups(groups>0));

for i = 1:length(u)
    [w(groups==u(i)) alpha(i)] = projectAux(w(groups==u(i)),alpha(i));
end
w = [w;alpha];

end

%% Function to solve the projection for a single group
function [w,alpha] = projectAux(w,alpha)
p = length(w);
nw = norm(w);
    if nw > alpha
       avg = (nw+alpha)/2;
       if avg < 0
           w(:) = 0;
           alpha = 0;
       else
           w = w*avg/nw;
           alpha = avg;
       end 
    end
end
