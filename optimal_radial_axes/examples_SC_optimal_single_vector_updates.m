
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2021 Manuel Rubio-Sánchez

% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required additional functions and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Olives data set
% http://www.interactivegraphics.org/Datasets_files/olives.txt

% arrow.m
% Available at https://www.mathworks.com/matlabcentral/fileexchange/278-arrow



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assistance and feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please write to manuel.rubio@urjc.es if encountering any issues




% Close figures, clear variables, clear console
close all;
clear;
clc;

% set initial random seed
seed = 1;
s = RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% constants

colors.reference_lines = [0.75 0.75 0.75];
colors.points = [0.7 0.8 1];
colors.SC_vectors = [0.5 0.5 0.5];
colors.OPT_vectors = [0 0 0];
colors.best_vector21 = [0 0 0.3];


% Select the number of variables 3 <= n <= 7
n = 4;




%% Load data

load olives_data.txt; 
X = olives_data;
Labels = X(:,1); % Class (This is a combination of the Region and Area attributes, which we coded into a single variable. There are 9 different classes)
X = X(:,2:end); % 8 numerical variables
original_variable_names = {'palmitic','palmitoleic','stearic','oleic','linoleic','linolenic','arachidic','eicosenoic'};
original_variable_names_best21 = {'palmitic+','palmitoleic+','stearic+','oleic+','linoleic+','linolenic+','arachidic+','eicosenoic+'};

[N,n_vars] = size(X);



% Standardize data 
Xz = zscore(X);





% choose n variables at random
perm = randperm(n_vars);
% Choose new variable
sel_var = perm(n+1);


Xdata = Xz(:,perm(1:n));
xe = Xz(:,sel_var);



variable_names = original_variable_names(perm(1:n));
variable_names_best21 = original_variable_names_best21(perm(1:n));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First example

% Generate random axis vectors
V = randn(n,2); % Random matrix of axis vectors with components drawn from a standard normal distribution


% Find optimal axis vector for new variable, through gradient descent with
% 4 initial points in W


step = 1/n/N/100; 
EPSILON = 0.00001;
  
W = [sqrt(2)/2 sqrt(2)/2; -sqrt(2)/2 sqrt(2)/2; sqrt(2)/2 -sqrt(2)/2; -sqrt(2)/2 -sqrt(2)/2];
mx = mean(Xdata)';
my = mean(xe);

I = eye(size(V,1));
    
f_min = Inf;
for i=1:4
    w = W(i,:)';

    E = Xdata'*Xdata;
    D = V'*V;
    C = V'*E*V;
    a = V'*Xdata'*xe;
    b = D*a;
    alpha = xe'*xe;


    gradient = compute_gradient(w,a,b,C,D,alpha);
    while norm(gradient)>EPSILON
        w_old = w;
        w = w - step*gradient;
        gradient = compute_gradient(w,a,b,C,D,alpha);
    end    
    
    Ve = [V; w'];   
    
    f = norm((Ve*Ve'-eye(n+1))*[Xdata xe]','fro')^2; % objective function
    
    if f<f_min
        f_min = f;
        v_best21 = w;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% SC plot

P = [Xdata xe]*[V; v_best21']; % Two-dimensional embedded points


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


% Draw best vector according to (21)
arrow([0 0],[v_best21(1), v_best21(2)],'Length',arrow_length,'BaseAngle',60,'Width',width,'FaceColor',colors.best_vector21,'EdgeColor',colors.best_vector21);
text(v_best21(1), v_best21(2), original_variable_names_best21(sel_var));

box on
axis equal

% Draw lines and circle for reference
t = 0:0.01:2*pi;
plot(cos(t),sin(t),'k:','Color',colors.reference_lines);

limitsX = get(ha,'XLim');
limitsY = get(ha,'YLim');
plot(limitsX*1.05,[0,0],'k:','Color',colors.reference_lines);
plot([0,0],limitsY*1.05,'k:','Color',colors.reference_lines);






%% Compute colored curve through the procedure related to (20)

n = n+1;
Xdata = [Xdata xe];

% Centere data
Xc = Xdata;
mX = mean(Xdata);
for i=1:N
    Xc(i,:) = Xc(i,:) - mX;
end 
Xdata = Xc;


I = eye(n-1);
range_a = 0:pi/180:pi;
s = size(range_a,2);
curve = zeros(s,2);
lambdas = zeros(s,1);
fs = zeros(s,1); % objective function values
truelambdas = zeros(s,1);

minfs = Inf;
ia = 0;
for a=range_a
    ia = ia+1;
    
    w = [cos(a); sin(a)];
    

    fp = zeros(4,1);
    for j=1:N
        x = Xdata(j,1:n-1)';
        y = Xdata(j,n);

        fp(4) = fp(4) + 2*x'*(V*V'-2*I)*V*w*y; % constant
        fp(3) = fp(3) + 2*(x'*V*(w*w')*V'*x + (w'*(V'*V)*w-2*(w'*w))*y^2);
        fp(2) = fp(2) + 6*x'*V*w*(w'*w)*y;
        fp(1) = fp(1) + 4*((w'*w)*y)^2;
    end       

    r = roots(fp);
    
    
    real_roots = zeros(3,1);
    n_real = 0;
    for j=1:3
        if imag(r(j))==0
            n_real = n_real+1;
            real_roots(n_real) = r(j);
        end
    end
    real_roots = real_roots(1:n_real);

    
    min_f_val = Inf;
    for k=1:n_real
        alfa = real_roots(k);
        
        A = [V*V'-I, alfa*V*w; alfa*w'*V', alfa^2*(w'*w)-1];
        gamma = - A*[mx;my];
        
        f = 0;
        for j=1:N
            x = Xdata(j,1:n-1)';
            y = Xdata(j,n);            
            
            f = f + alfa^4*((w'*w)*y)^2 + 2*alfa^3*x'*V*w*(w'*w)*y + alfa^2*(x'*V*(w*w')*V'*x + (w'*(V'*V)*w-2*(w'*w))*y^2) + 2*alfa*x'*(V*V'-2*I)*V*w*y + x'*(V*V'-I)^2*x + y^2;
        end
        
        
        
        % Test: (f and f_alt should be the same)
%         Ve = [V; alfa*w'];
%         f_alt = norm((Ve*Ve'-eye(n))*Xdata','fro')^2;
        
        if f<min_f_val
            min_f_val = f;
            lambda = alfa;
        end
    end
    
    truelambdas(ia) = lambda;
    
    if lambda<0
        lambda = -lambda;
        w = -w;
    end
    
    fs(ia) = min_f_val;
    lambdas(ia) = lambda;
    curve(ia,:) = lambda*w';
    
    if min_f_val<minfs
        minfs = min_f_val;
    end
end

C = colormap();

maxval = max(fs);
minval = min(fs);

for i=2:s

    val = 63*(fs(i) - minval)/(maxval - minval) + 1;
    color = floor(val);

    %color = 65-color;

    if color==0
        color = 1;
    end

    plot([curve(i-1,1),curve(i,1)],[curve(i-1,2),curve(i,2)],'k-','LineWidth',5,'Color',C(color,:));
    
end

val = 63*(fs(1) - minval)/(maxval - minval) + 1;
color = floor(val);
color = 65-color;
if color==0
    color = 1;
end
plot([curve(1,1),curve(s,1)],[curve(1,2),curve(s,2)],'k-','LineWidth',5,'Color',C(color,:));





% Draw color bar
hc = colorbar;
labels = cell(6,1);
i=0;
for val=minval:(maxval-minval)/5:maxval
    i=i+1;
    labels{i} = sprintf('%.2f',val); 
end
set(hc,'TickLabels',labels);



