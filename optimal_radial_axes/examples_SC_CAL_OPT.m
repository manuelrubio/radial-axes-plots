
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
colors.SC_vectors = [0.25 0.25 0.75];
colors.OPT_vectors = [0 0 0];


% Select the number of variables1 3 <= n <= 8
n = 5;




%% Load data

load olives_data.txt; 
X = olives_data;
Labels = X(:,1); % Class (This is a combination of the Region and Area attributes, which we coded into a single variable. There are 9 different classes)
X = X(:,2:end); % 8 numerical variables
original_variable_names = {'palmitic','palmitoleic','stearic','oleic','linoleic','linolenic','arachidic','eicosenoic'};
original_variable_names_opt = {'palmitic*','palmitoleic*','stearic*','oleic*','linoleic*','linolenic*','arachidic*','eicosenoic*'};

[N,n_vars] = size(X);
perm = randperm(n_vars);
% choose n variables at random
X = X(:,perm(1:n));
variable_names = original_variable_names(perm(1:n));
variable_names_opt = original_variable_names_opt(perm(1:n));



X01 = zeros(N,n); % Data normalized to lie in the [0,1] interval
min_data = min(X);
max_data = max(X);
for i=1:n
    if max_data(i) ~= min_data(i)
        for j=1:N
            X01(j,i) = (X(j,i) - min_data(i))/(max_data(i) - min_data(i));
        end
    else
        X01(:,i) = 0.5*ones(N,1);
    end
end

% For standardized data use Xz
%Xz = zscore(X); % Standardized data


Xdata = X01;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First example

% Generate random axis vectors
V = randn(n,2); % Random matrix of axis vectors with components drawn from a standard normal distribution

% draw example
V_star = draw_SC_show_CAL_OPT_estimation_errors(V,Xdata,variable_names,variable_names_opt,colors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same example, but scaling V as described in (18) and (19)

% Scale V
theta = sqrt(norm(V_star,'fro')/norm(V,'fro')); % (18)
V = theta*V; % (19) 

% draw example
draw_SC_show_CAL_OPT_estimation_errors(V,Xdata,variable_names,variable_names_opt,colors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same example, but scaling V as described in (18) and (19)

% Scale factor
theta = sqrt(norm(V_star,'fro')/norm(V,'fro')); % (18)

% draw example
draw_SC_show_CAL_OPT_estimation_errors(theta*V,Xdata,variable_names,variable_names_opt,colors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Replace one vector by its optimal counterpart

% Update V
sel_var = ceil(rand(1)*n); % select variable at random
V(sel_var,:) = V_star(sel_var,:); % update axis vector

% draw example
draw_SC_show_CAL_OPT_estimation_errors(V,Xdata,variable_names,variable_names_opt,colors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Replace all vectors by the ones obtained through OPT

% Update V
V = V_star; % update axis vector

% draw example
draw_SC_show_CAL_OPT_estimation_errors(V,Xdata,variable_names,variable_names_opt,colors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


