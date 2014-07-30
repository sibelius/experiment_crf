colors = getColorsRGB;

% Plot Groups
for g = 1:nGroups
    h(:,g) = plot(lambdaValues,w(groups==g,:));
    set(h(:,g),'LineWidth',3,'Color',colors(g,:));
    groupName{g} = sprintf('Group %d\n',g);
end
legend(h(1,:),groupName{:});

% Set X-axis
[i j] = find(w);
xlim([0 lambdaValues(max(min(j)-1,1))]);

set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

% Plot Vertical Lines
for i = 1:length(lambdaValues)-1
    for g = nGroups:-1:1
        if any(abs(w(groups==g,i+1)) < 1e-4 ~= abs(w(groups==g,i)) < 1e-4)
           v = vline(lambdaValues(i));
           set(v,'Color',colors(g,:),'LineWidth',1+(g/nGroups));
        end
    end
end