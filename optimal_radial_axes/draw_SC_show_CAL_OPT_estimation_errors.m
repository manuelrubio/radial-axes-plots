function V_star = draw_SC_show_CAL_OPT_estimation_errors(V,X,variable_names,variable_names_opt,colors)

[N,n] = size(X);


%% SC plot

P = X*V; % Two-dimensional embedded points


% initialize plot
hf = figure;
set(hf, 'Renderer', 'painters');
ha = axes();
hold
set(ha,'FontSize',12);



% Draw plotted points
for j=1:N
    plot(P(j,1),P(j,2),'k.','Color',colors.points,'MarkerSize',10); 
end

% Draw vectors
arrow_length = 10;
width = 1;
for i=1:n
    arrow([0 0],[V(i,1), V(i,2)],'Length',arrow_length,'BaseAngle',60,'Width',width,'FaceColor',colors.SC_vectors,'EdgeColor',colors.SC_vectors);
end

% Draw variable names
for i=1:n
    text(V(i,1), V(i,2), variable_names(i));
end

box on
axis equal


% Draw lines and circle for reference
t = 0:0.01:2*pi;
plot(cos(t),sin(t),'k:','Color',colors.reference_lines);

limitsX = get(ha,'XLim');
limitsY = get(ha,'YLim');
plot(limitsX*1.05,[0,0],'k:','Color',colors.reference_lines);
plot([0,0],limitsY*1.05,'k:','Color',colors.reference_lines);



% Total squared estimation error without calibration:
e = norm(P*V'-X,'fro')^2;
fprintf('Total squared estimation error without calibration: %.2f\n',e);



%% CAL

alphas = zeros(n,1);
betas = zeros(n,1);
for i=1:n

    sumDotProducts = 0;
    for j=1:N
        sumDotProducts = sumDotProducts + P(j,:)*V(i,:)';
    end

    sumWeightedDotProducts = 0;
    for j=1:N
        sumWeightedDotProducts = sumWeightedDotProducts + X(j,i)*P(j,:)*V(i,:)';
    end

    sumSquaredDotProducts = 0;
    for j=1:N
        sumSquaredDotProducts = sumSquaredDotProducts + (P(j,:)*V(i,:)')^2;
    end

    meanXi = mean(X(:,i));


    alphas(i) = (sumWeightedDotProducts - meanXi*sumDotProducts) / (sumSquaredDotProducts - sumDotProducts^2/N);
    betas(i) = meanXi - alphas(i)*sumDotProducts/N;

end

% Total squared estimation error with calibration (CAL):
e_cal = norm( P*V'*diag(alphas) - X + ones(N,1)*betas','fro')^2;
fprintf('Total squared estimation error with calibration (CAL): %.2f\n',e_cal);



%% OPT

% Compute centered plotted points
Pc = P;
mP = mean(P);
for j=1:N
    Pc(j,:) = Pc(j,:) - mP;
end 

V_star = (pinv(Pc)*X)';


% Draw OPT vectors
arrow_length = 10;
width = 1;
for i=1:n
    arrow([0 0],[V_star(i,1), V_star(i,2)],'Length',arrow_length,'BaseAngle',60,'Width',width,'FaceColor',colors.OPT_vectors,'EdgeColor',colors.OPT_vectors);
end

% Draw variable names
for i=1:n
    text(V_star(i,1), V_star(i,2), variable_names_opt(i));
end



% Compute shifts (gammas) for OPT
gammas = zeros(n,1);
for i=1:n
    
    sumDotProducts = 0;
    for j=1:N
        sumDotProducts = sumDotProducts + P(j,:)*V_star(i,:)';
    end
    
    meanXi = mean(X(:,i));
    
    gammas(i) = meanXi - sumDotProducts/N;
end



% Total squared estimation error for OPT:
e_opt = norm( P*V_star' - X + ones(N,1)*gammas','fro')^2;
fprintf('Total squared estimation error for OPT: %.2f\n',e_opt);


str = sprintf('Est errors. No cal: %.1f  CAL: %.1f  OPT: %.1f',e,e_cal,e_opt);
title(str);