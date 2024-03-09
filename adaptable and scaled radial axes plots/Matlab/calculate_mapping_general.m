function P = calculate_mapping_general(algorithm,X,V,W,vector_norm,normalize_axes,chosen_variable)
%CALCULATE_MAPPING_GENERAL
%   
% P = CALCULATE_MAPPING_GENERAL(X,V,W,algorithm,normalize) computes the mappings of several
% radial axes methods
% algorithm defines the radial method ("SC": Star Coordinates, "RadViz", "SRA": Scaled Radial Axes, "Adaptable": Adaptable Radial Axes Plots
% X is an N by n matrix whose rows contain the n dimensional data samples.
% V is an n by m matrix whose rows define the method's axis vectors. 
% W is an n by n diagonal matrix defining nonnegative weights for each variable
% vector_norm is the vector norm associated with adaptable radial axes plots
% normalize_axes indicates whether to divide the vectors by their squared norm (for Scaled Radial Axes)
% chosen_variable is the selected attribute for constrained adaptable radial axes plots
% The low-dimensional embeddings are stored in the N by m matrix P.

% Adaptable Radial Axes Plots:
% M. Rubio-Sánchez, A. Sanchez, Dirk J. Lehmann: Adaptable Radial Axes Plots for Improved Multivariate Data Visualization. Computer Graphics Forum. Vol 36, no. 3, pp. 389-399, Jun. 2017. DOI: 10.1111/cgf.13196. (Q1, CORE-A, GGS Class 3)

% Scaled Radial Axes:
% A. Sanchez, C. Soguero-Ruiz, I. Mora-Jiménez, F. J. Rivas-Flores, Dirk J. Lehmann, M. Rubio-Sánchez: Scaled Radial Axes for Interactive Visual Feature Selection: A Case Study for Analyzing Chronic Conditions. Expert Systems with Applications. Vol 100, pp. 182-196, Jun. 2018. DOI: 10.1016/j.eswa.2018.01.054. (Q1)


[N,n] = size(X);
m = size(V,2);
    
if strcmp(algorithm,'SC')
    
    P = X*V;
	
else if strcmp(algorithm,'SRA')	% Scaled Radial Axes

	for i=1:n
		norm_squared_of_ith_row_of_V = V(i,:)*V(i,:)';
		if norm_squared_of_ith_row_of_V ~= 0
			V(i,:) = V(i,:)/norm_squared_of_ith_row_of_V;
		end
	end 
	
	P = X*V;
    
else if strcmp(algorithm,'RadViz')

        minimum = min(X);
        maximum = max(X);
        X_Radviz = (X - repmat(minimum,N,1))./(repmat(maximum,N,1) - repmat(minimum,N,1));
      
        for i=1:N
            sum_row = sum(X_Radviz(i,:));
            if sum_row==0
                X_Radviz(i,:) = ones(1,n)/n;
            else
                X_Radviz(i,:) = X_Radviz(i,:)/sum_row;
            end
        end
        
        P = X_Radviz * V;
        
    else
        
        if normalize_axes==1 % (For Scaled Radial Axes)
            for i=1:n
                norm_squared_of_ith_row_of_V = V(i,:)*V(i,:)';
                if norm_squared_of_ith_row_of_V ~= 0
                    V(i,:) = V(i,:)/norm_squared_of_ith_row_of_V;
                end
            end 
        end

        if strcmp(algorithm,'Adaptable')
            
            switch vector_norm
                case 1,
                    
                    options = optimset('Algorithm', 'interior-point', 'Display', 'off');
                    A = [-eye(n), -W * V ; -eye(n), W * V ];
                    f = [ones(1,n),zeros(1,m)];

                    P = zeros(N,m);            
                    for i=1:N
                        b = [ -W * X(i,:)' ; W * X(i,:)'];

                        p_star = linprog(f,A,b,[],[],[],[],[],options);

                        P(i,:) = p_star(end-m+1:end)';
                    end 
            
                case 2,    
                    P = X*(pinv(W * V)*W)';                
                    
                case Inf,

                    options = optimset('LargeScale', 'off', 'Algorithm', 'active-set', 'Display', 'off');

                    A = [-ones(n,1), - W * V ; -ones(n,1), + W * V];

                    f = [1,zeros(1,m)];

                    P = zeros(N,m);               
                    for i=1:N
                        b = [ -W * X(i,:)' ; W * X(i,:)'];

                        p_star = linprog(f,A,b,[],[],[],[],[],options);

                        P(i,:) = p_star(end-m+1:end)';
                    end                      
            end
            
        else if strcmp(algorithm,'Adaptable exact')
 
                v = V(chosen_variable,:);
                x = X(:,chosen_variable);        

                switch vector_norm
                    case 1,

                        options = optimset('Display', 'off');
                        
                        A = [-eye(n), -V; -eye(n), V];
                        f = [ones(1,n),0,0];
                        Aeq = [zeros(1,n), v];
                        
                        P = zeros(N,2);
                        for i=1:N
                            b = [ -X(i,:)'; X(i,:)'];

                            p_star = linprog(f,A,b,Aeq,x(i),[],[],[],options);

                            if ~isempty(p_star)
                                P(i,:) = p_star(end-1:end)';
                            end
                        end 
                        
                    case 2,    
                        
                        options = optimoptions('lsqlin','Algorithm','active-set','Display','off');
                        
                        P = zeros(N,2);
                        for i=1:N
                            p = lsqlin(V,X(i,:)',[],[],v,x(i),[],[],[],options);
                            if ~isempty(p)
                                P(i,:) = p';
                            end
                        end 

                    case Inf,

                        options = optimset('Display', 'off');

                        A = [-ones(n,1), -V ; -ones(n,1), V];
                        f = [1,zeros(1,2)];
                        Aeq = [0, v];

                        P = zeros(N,2);                         
                        for i=1:N
                            b = [ -X(i,:)' ; X(i,:)'];

                            p_star = linprog(f,A,b,Aeq,x(i),[],[],[],options);

                            if ~isempty(p_star)
                                P(i,:) = p_star(end-1:end)';
                            end
                        end  
        
                end        
            
            else if strcmp(algorithm,'Adaptable ordered') % solved through CVX

                    k = chosen_variable
                    [s, I] = sort(X(:,k));
        
                    switch vector_norm
                        case 1,

                        cvx_begin
                            cvx_quiet(true);    
                            variable P(N,2);
                            z = 0;
                            for i=1:N
                                z = z + norm(V*P(i,:)' - X(i,:)',1);
                            end

                            minimize( z );

                            P(I(1:N-1),:)*V(k,:)' <= P(I(2:N),:)*V(k,:)'

                        cvx_end  
        
                        case 2,    
                            
                            cvx_begin
                                cvx_quiet(true);    
                                variable P(N,2);

                                minimize( norm(P*V' - X,'fro') );

                                P(I(1:N-1),:)*V(k,:)' <= P(I(2:N),:)*V(k,:)'
                            cvx_end

                        case Inf,

                            cvx_begin
                                cvx_quiet(true);    
                                variable P(N,2);
                                variable x;
                                expression z(n);
                                for j=1:n
                                    z(j) = norm(P*V(j,:)' - X(:,j),Inf);
                                end

                                minimize( x );

                                for j=1:n
                                    z(j) <= x
                                end

                                P(I(1:N-1),:)*V(k,:)' <= P(I(2:N),:)*V(k,:)'

                            cvx_end 
                            
                            A = P*V'-X;
                            val = 0;
                            for i=1:n
                                for j=1:N
                                    if abs(A(j,i))>val
                                        val = abs(A(j,i));
                                    end
                                end
                            end
                            
                            stop = 1;
                    end        

                end            
            end            
        end
    end
end

		


            




 