% Nagaprasad Rudrapatna

N = 8;
SAMPLE = zeros(N-1, 1);
abs_error = zeros(N-1, 1);
h = zeros(N-1, 1);
beta = 1.7; 
tol = 1e-10; 
L = 1;
max_iter = 150000;

for m = 2:N   
    
% setting up uniform grid
n_x = 1 + 2^m; 
n_y = 1 + 2^m;

x = linspace(0, L, n_x);
y = linspace(0, L, n_y);

% setting up increments
delta_x = x(2) - x(1);
h(m-1) = delta_x;

% RHS & BC's
exact = @(x,y) exp(x^2 + y^2);
s = @(x,y) exp(x^2)*(4*x^2 * exp(y^2) + (4*y^2 + 4) * exp(y^2)); % Laplacian of exact
B1 = @(x) exp(x^2 + L^2); 
B2 = @(x) exp(x^2); 
B3 = @(y) exp(L^2 + y^2); 
B4 = @(y) exp(y^2); 

ref = zeros(n_x, n_y);

for k = 1:n_x
    for j = 1:n_y
        ref(k, j) = exact(x(k), y(j));
    end
end

[phi, sample_SOR, sample_conv, iter] = sor2d(beta, n_x, n_y, h(m-1), tol, s, B1, B2, B3, B4, L, max_iter);

fprintf("\nSOR: %d iterations \n", iter) 

SAMPLE(m-1) = sample_conv;
abs_error(m-1) = max(max(abs(phi - ref))); 

end

figure(1)
loglog(h, abs_error, '--hm');
xlabel('Step Size');
ylabel('Absolute Error');

set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[4 3]);
set(gcf,'PaperPosition',[0 0 4 3]);

rho1 = zeros(length(SAMPLE)-1, 1); 

for k = 2:length(SAMPLE)-1
     rho1(k) = (SAMPLE(k+1) - SAMPLE(k)) / (SAMPLE(k) - SAMPLE(k-1));
     fprintf("\n%.5f\n", rho1(k));
end

rho2 = zeros(length(sample_SOR)-1, 1);

for k = 2:length(sample_SOR)-1
    rho2(k) = (sample_SOR(k+1) - sample_SOR(k)) / (sample_SOR(k) - sample_SOR(k-1));
  %  fprintf("%.5f\n", rho2(k));
end
