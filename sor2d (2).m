% Nagaprasad Rudrapatna

function [sol, sampled_comp, sampled_convergence, it] = sor2d(beta, n_x, n_y, h, tol, S, b1, b2, b3, b4, L, max_iter)
% Description: sor2d implements the 2D SOR method to solve the discretized form of Poisson's Equation
% Inputs:
%    beta: relaxation parameter
%    n_x, n_y: number of x-points / y-points in a uniform grid
%    h: scalar step size
%    tol: the user-specified error tolerance level 
%    S: RHS of Poisson's Equation
%    b1, b2, b3, b4: four boundary conditions
%    L: dimension of box
%    max_iter: maximum number of iterations
% Returns:
%    sol: the solution matrix
%    sampled_comp: the matrix of sampled components used to determine the rate of convergence of SOR
%    sampled_convergence: the sampled component used to determine order of convergence 
%    it: number of iterations to achieve tol

sol = zeros(n_x, n_y);
x = linspace(-L, L, n_x); % (-L,L): gravity
y = linspace(-L, L, n_y);

% applying BC's
for i = 1:n_x
    sol(i, n_y) = b1(x(i)); % top
    sol(i, 1) = b2(x(i)); % bottom
end

for j = 1:n_y
    sol(n_x, j) = b3(y(j)); % right
    sol(1, j) = b4(y(j)); % left
end

sol_old = sol;
sampled_comp = zeros(1, 1);
it = 0; 
error = 10; 

while (error > tol && it < max_iter)
        
    for i = 2:n_x-1
        for j = 2:n_y-1
           sol(i,j) = (1/4) * (sol_old(i+1, j) + sol(i-1, j) + sol_old(i, j+1) + sol(i, j-1) - (h^2) * S(x(i), y(j)));
        end
    end
    sol = beta*sol + (1 - beta)*sol_old;    
    error = max(max(abs(sol_old - sol))); % infinity vector norm
    sol_old = sol;     
    sampled_comp(it+1) = sol((n_x + 1) / 2, (n_y + 1) / 2);
    it = it + 1;
end

sampled_convergence = sol((n_x + 1) / 2, (n_y + 1) / 2);

end