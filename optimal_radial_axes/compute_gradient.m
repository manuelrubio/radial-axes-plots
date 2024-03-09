function gradient = compute_gradient(v_n,a,C,D,x_n)

% Equation (22) of
% M. Rubio-Sánchez, Dirk Lehmann, Alberto Sanchez, José Luis Rojo-Álvarez. 
% Optimal Axes for Data Value Estimation in Star Coordinates and Radial 
% Axes Plots. Computer Graphics Forum 40 (3), pp. 483-494. June, 2021.
% DOI: 10.1111/cgf.14323.

gradient = 2*C*v_n + 2*(v_n'*v_n)*a + (2*D - 4*eye(2) + 4*(v_n*v_n'))*(a + x_n'*x_n*v_n);
