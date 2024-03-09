import numpy as np
import numpy.matlib
from scipy.optimize import linprog  
# CVX
import cvxpy as cvx

# MAPPING
#   
# P = MAPPING(algorithm, X, V, W, vector_norm, chosen_variable) 
# computes the mappings of several radial axes methods algorithm defines the radial method
# X is an N by n matrix wose rows contain the n dimensional data samples.
# V is an n by m matrix whose rows define the method's axis vectors. 
# W is an n by n diagonal matrix defining nonnegative weights for each variable
# vector_norm is the vector norm associated with adaptable radial axes plots
# chosen_variable is the selected attribute for constrained adaptable radial axes plots
# The low-dimensional embeddings are stored in the N by m matrix P.
def mapping(algorithm, X, V, W=None, vector_norm=None, chosen_variable=None):
    X = np.matrix(X)
    V = np.matrix(V)
    W = np.matrix(W)

    [N,n] = X.shape
    m = V.shape[1]

    if algorithm == 'SC':
        P = X * V

    elif algorithm == 'RadViz':
        minimum = np.amin(X, axis = 0)
        maximum = np.amax(X ,axis = 0)        
        X_Radviz = (X - numpy.matlib.repmat(minimum, N, 1)) / (numpy.matlib.repmat(maximum, N, 1) - numpy.matlib.repmat(minimum, N, 1))     
        for i in range(N):
            sum_row = np.sum(X_Radviz[i,])
            if sum_row == 0:
                X_Radviz[i,] = np.ones((n)) / n
            else:
                X_Radviz[i,] = X_Radviz[i,] / sum_row        
        P = X_Radviz * V    
    
    elif algorithm == 'Adaptable':
        if vector_norm == 1:
            A = np.vstack((np.hstack((-np.identity(n) , -W * V)), np.hstack((-np.identity(n) , W * V))))  
            f = np.hstack((np.ones((n)) , np.zeros((m))))
            P = np.zeros((N, m))
            for i in range(N):
                b = np.vstack((-W * X[i,].transpose(), W * X[i,].transpose()))
                p_star = linprog(f, A_ub=A, b_ub=b, bounds=(None, None), method='interior-point')     
                P[i,] = p_star.x[-m:].transpose()

        elif vector_norm == 2:
            P = X * (np.linalg.pinv(W * V) * W).transpose()

        elif vector_norm == 'Inf':
            A = np.vstack((np.hstack((-np.ones((n,1)), -W * V)), np.hstack((-np.ones((n,1)), W * V))))	
            f = np.hstack((1, np.zeros((m))))
            P = np.zeros((N, m))  
            for i in range(N):    
                b = np.vstack((-W * X[i,].transpose() , W * X[i,].transpose()))
                p_star = linprog(f, A_ub=A, b_ub=b, bounds=(float('-inf'), None), method='interior-point')     
                P[i,] = p_star.x[-m:].transpose()       

    elif algorithm == 'Adaptable exact':
        v = np.squeeze(np.asarray(V[chosen_variable,]))
        x = np.squeeze(np.asarray(X[:,chosen_variable]))
        if vector_norm == 1:
            A = np.vstack((np.hstack((-np.identity(n), -V)),np.hstack((-np.identity(n), V)))) 
            f = np.hstack((np.ones((n)), 0, 0))
            Aeq = np.matrix([np.hstack((np.zeros((n)), v))])
            P = np.zeros((N, m))
            for i in range(N):
                b = np.vstack((-X[i,].transpose(),X[i,].transpose()))
                p_star = linprog(f, A, b, Aeq, np.array([x[i]]), bounds=(float('-inf'), None), method='interior-point') 
                if p_star:
                    P[i,] = p_star.x[-m:].transpose()
                    
        elif vector_norm == 2:
            P = cvx.Variable(N, m)
            obj = cvx.Minimize(cvx.norm(P * V.T - X, "fro"))
            constraints = [P * v.T == x]
            prob = cvx.Problem(obj, constraints)
            prob.solve()
            P = P.value

        elif vector_norm == 'Inf':
            A = np.vstack((np.hstack((-np.ones((n, 1)), -V)), np.hstack((-np.ones((n, 1)), V))))
            f = np.hstack((1, np.zeros((2))))
            Aeq = np.matrix([np.hstack((0, v))])
            P = np.zeros((N, m))
            for i in range (N):
                b = np.vstack((-X[i,].transpose(), X[i,].transpose()))
                p_star = linprog(f, A, b, Aeq, x[i], bounds=(float('-inf'), None), method='interior-point')
                if p_star:
                    P[i,] = p_star.x[-m:].transpose()
            
    elif algorithm == 'Adaptable ordered':
        k = chosen_variable
        I = np.argsort(np.asarray(X)[:,k])
        constraints_list = []
        if vector_norm == 1:
            P = cvx.Variable(N, m)
            z = 0
            for i in range(N):
                z = z + cvx.norm(V * P[i,:].T - X[i,].T, 1)
            obj = cvx.Minimize(z)
            for i in range (N-1):
                constraints_list.append(P[I[i],:] * V[k,].T <= P[I[i+1],:] * V[k,].T)
            prob = cvx.Problem(obj, constraints_list)
            prob.solve()
            P = P.value

        elif vector_norm == 2:
            P = cvx.Variable(N, m)
            obj = cvx.Minimize(cvx.norm(P * V.T - X, "fro"))
            for i in range (N-1):
                constraints_list.append(P[I[i],:] * V[k,].T <= P[I[i+1],:] * V[k,].T)
            prob = cvx.Problem(obj, constraints_list)    
            prob.solve()
            P = P.value

        elif vector_norm == 'Inf':
            P = cvx.Variable(N, m)
            x = cvx.Variable()
            obj = cvx.Minimize(x)
            for j in range(n):
                z = cvx.norm(P * V[j,].T - X[:,j], "inf")
                constraints_list.append(z <= x)                
            for i in range(N-1):
                constraints_list.append(P[I[i],:] * V[k,].T <= P[I[i+1],:] * V[k,].T)
            prob = cvx.Problem(obj, constraints_list)
            prob.solve()
            P = P.value   
            
    return(P)     