% Nagaprasad Rudrapatna

N = 8;
SAMPLE = zeros(N-1, 1);
h = zeros(N-1, 1);
beta = 1;
tol = 1e-9; 
L = 196; 
max_iter = 150000;

for m = 2:N   
    
% setting up uniform grid
n_x = 1 + 2^m; 
n_y = 1 + 2^m;

x = linspace(-L, L, n_x); % centered at the origin (make appropriate change in .m)
y = linspace(-L, L, n_y);

% setting up increments
delta_x = x(2) - x(1);
h(m-1) = delta_x; 

% RHS & BC's
G = 1; % scaling to get phi / G
s = @(x,y) 4 * pi * G * exp(-(x^2 + y^2)); 
B1 = @(x) 0; 
B2 = @(x) 0; 
B3 = @(y) 0; 
B4 = @(y) 0; 

[phi, sample_SOR, sample_conv, iter] = sor2d(beta, n_x, n_y, h(m-1), tol, s, B1, B2, B3, B4, L, max_iter);
fprintf("\nSOR: %d iterations \n", iter) 

SAMPLE(m-1) = sample_conv;

end

figure(2)
surf(y, x, phi)
grid off
xlabel('X');
ylabel('Y');
zlabel('Gravitational Potential')
colormap spring

set(gcf,'PaperUnits','inches'); 
set(gcf,'PaperSize',[4 3]);
set(gcf,'PaperPosition',[0 0 4 3]);

rho1 = zeros(length(SAMPLE)-1, 1);

for k = 2:length(SAMPLE)-1
     rho1(k) = (SAMPLE(k+1) - SAMPLE(k)) / (SAMPLE(k) - SAMPLE(k-1));
     % fprintf("%.5f\n", rho1(k));
end

rho2 = zeros(length(sample_SOR)-1, 1);

for k = 2:length(sample_SOR)-1
    rho2(k) = (sample_SOR(k+1) - sample_SOR(k)) / (sample_SOR(k) - sample_SOR(k-1));
    % fprintf("%.5f\n", rho2(k));
end
