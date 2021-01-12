clearvars; clc; clf
figure(1); axis image; axis off; colormap(gray);
% set(gcf, 'Position', get(0, 'Screensize'));
global I0 S0 M N
global a b r theta
global S0DM xstar ystar C1 C2
global eps alpha lambda1 lambda2 nu dt
eps = 1e-6; % For Heaviside and delta functions
%% PARAMETERS
Image_filename = "Drawing4.png";
ShapePrior_filename = "Shape1.png";
alpha = 4.0;
lambda1 = 0.2;
lambda2 = 0.2;
nu = 0.5;
dt = 5e-2;
%% Initialize
[I0, S0, M, N] = ReadInputs(Image_filename, ShapePrior_filename);
[Phi0, L0] = Initialize(50, 1, 1);
a = 10; b = 10; r = 1.6; theta = -pi/4;
[Psi0, S0DM, xstar, ystar] = DomainExtension(a, b, r, theta);
Psi_old = Psi0; L_old = L0; Phi_old = Phi0; 
%% SOLVE
Iter = 0;
while(Iter < 100)
    [C1, C2] = CalculateC1C2(Phi_old, Psi_old);
    Psi = UpdateShapePrior(Phi_old, L_old, Psi_old);
    L = UpdateLabel(Phi_old, L_old, Psi);
    [Phi, E] = UpdateSegmentation(Phi_old, L, Psi);
    errPsi = norm(HeaviSide(Psi) - HeaviSide(Psi_old), inf);
    errL = norm(HeaviSide(L) - HeaviSide(L_old), inf); 
    errPhi = norm(HeaviSide(Phi) - HeaviSide(Phi_old), inf); 
    Psi_old = Psi;
    L_old = L;
    Phi_old = Phi;
    subplot(1, 3, 1);
    imshow(I0); title(['\phi, ' 'Error = ' num2str(errPhi)]); hold on;
    contour(Phi, [0 0], 'r', 'LineWidth', 2); hold off;
    drawnow;
    subplot(1, 3, 2);
    imshow(I0); title(['L, ' 'Error = ' num2str(errPhi)]); hold on;
    contour(L, [0 0], 'r', 'LineWidth', 2); hold off;
    drawnow;
    subplot(1, 3, 3);
    imshow(I0); title(['\psi, ' 'Error = ' num2str(errPsi)]); hold on;
    contour(Psi, [0 0], 'r', 'LineWidth', 2); hold off;
    drawnow;  
    Iter = Iter + 1;
%     break;
end

%% FUNCTIONS
function Y = HeaviSide(Phi)
global eps;
Y = 0.5 + atan(Phi/eps)/pi;
% Y = (Phi >= 0)*1;
end

function Y = DiracDelta(Phi)
global eps;
Y = (eps./(eps^2 + Phi.^2))/pi*10;
% Y = (abs(Phi) <= 2);
end

function [I0, S0, M, N] = ReadInputs(Image_filename, ShapePrior_filename)
I0 = imread(Image_filename);
I0 = rgb2gray(I0);
I0 = im2double(I0);
% I0 = imnoise(I0, 'salt & pepper', 0.02);

S0 = imread(ShapePrior_filename);
S0 = rgb2gray(S0);
S0 = im2double(S0);
S0 = S0 > 0;

[M, N] = size(I0);
end

function [Phi0, L0] = Initialize(r, nx, ny)
global M N
r = min([r, N/2/nx, M/2/ny]);
cx = 0+floor(N/2/nx):round(N/nx):N-floor(N/2/nx);
cy = 0+floor(M/2/ny):round(M/ny):M-floor(M/2/ny);
n = length(cx)*length(cy);
[cx, cy] = meshgrid(cx, cy);
Phi0 = zeros(M, N);
for i = 1:n
    Cx = cx(i); Cy = cy(i);
    for x = 1:N
        for y = 1:M
            Phi0(y, x) = Phi0(y, x) + ((x-Cx).^2 + (y-Cy).^2 < r^2).*1;
        end
    end
end
Phi0 = bwdist(1-Phi0) - bwdist(Phi0) - 0.5;
L0 = Phi0;
end

function [Psi, S0DM, xstar, ystar] = DomainExtension(a, b, r, theta)
global S0 M N
CornerPoints = [1-round(N/2)-a, 1-round(N/2)-a, N-round(N/2)-a, N-round(N/2)-a; ...
                1-round(M/2)-b, M-round(M/2)-b, 1-round(M/2)-b, M-round(M/2)-b];
CornerPoints = [cos(theta) sin(theta); -sin(theta) cos(theta)]/r*CornerPoints;
CornerPoints = round(CornerPoints);
[x_min, x_max] = bounds([CornerPoints(1, :), 1-round(N/2), N-round(N/2)]);
[y_min, y_max] = bounds([CornerPoints(2, :), 1-round(M/2), M-round(M/2)]);
NN = x_max - x_min + 1; MM = y_max - y_min + 1;
S0DM = zeros(MM, NN);
S0DM(y_max-M+round(M/2)+1:y_max+round(M/2), 2-x_min-round(N/2):1+N-round(N/2)-x_min) = S0;
S0DM = imfill(S0DM, 'holes');
S0DM = bwdist(1-S0DM) - bwdist(S0DM) - 0.5;
xstar = zeros(M, N); ystar = zeros(M, N);
Psi = zeros(M, N);
for i = 1:M
    for j = 1:N
        temp_star = [cos(theta) sin(theta); -sin(theta) cos(theta)]/r * [j-round(N/2)-a; M-round(M/2)+1-i-b];
        xstar(i, j) = round(temp_star(1)); ystar(i, j) = round(temp_star(2));
        ii = y_max - round(temp_star(2)) + 1; jj = round(temp_star(1)) - x_min + 1;
        Psi(i, j) = r*S0DM(ii, jj);
    end
end
end

function [C1, C2] = CalculateC1C2(Phi_old, Psi_old)
global I0
global nu
C1 = sum(I0.*(HeaviSide(Phi_old) + nu*HeaviSide(Psi_old)))/sum(HeaviSide(Phi_old) + nu*HeaviSide(Psi_old));
C2 = sum(I0.*(1-HeaviSide(Phi_old) + nu*(1-HeaviSide(Psi_old))))/sum(1-HeaviSide(Phi_old) + nu*(1-HeaviSide(Psi_old)));
end

function Psi = UpdateShapePrior(Phi_old, L_old, Psi_old)
global I0 M N
global a b r theta
global S0DM xstar ystar C1 C2
global alpha nu dt
S0DMx = S0DM - circshift(S0DM, 1, 2);
S0DMy = S0DM - circshift(S0DM, -1, 1);
temp = DiracDelta(Psi_old).*(2*alpha*(HeaviSide(Psi_old) - HeaviSide(Phi_old).*HeaviSide(L_old)) + nu*((I0 - C1).^2 - (I0 - C2).^2));
temp_a = 0; temp_b = 0; temp_r = 0; temp_theta = 0;
y_max = max(max(ystar)); x_min = min(min(xstar));
for i = 2:M
    for j = 2:N
        ii = y_max - ystar(i, j) + 1; jj = xstar(i, j) - x_min + 1;        
        temp_a = temp(i, j)*(S0DMx(ii, jj)*cos(theta) - S0DMy(ii, jj)*sin(theta)) + temp_a;
        temp_b = temp(i, j)*(S0DMx(ii, jj)*sin(theta) + S0DMy(ii, jj)*cos(theta)) + temp_b;
        temp_r = temp(i, j)*(-S0DM(ii, jj) + S0DMx(ii, jj)*xstar(i, j) + S0DMy(ii, jj)*ystar(i, j)) + temp_r;
        temp_theta = temp(i, j)*(-r*S0DMx(ii, jj)*ystar(i, j) + r*S0DMy(ii, jj)*xstar(i, j)) + temp_theta;
    end
end
a = dt*temp_a + a; b = dt*temp_b + b; r = dt*temp_r + r; theta = dt*temp_theta + theta;
[Psi, S0DM, xstar, ystar] = DomainExtension(a, b, r, theta);
end

function L = UpdateLabel(Phi_old, L_old, Psi)
global M N
global alpha lambda1 lambda2
dt = 1;
L = L_old;
Iter = 0; IterMax = 3;
while(Iter < IterMax)   
    delta = DiracDelta(L)/eps*100;
    for i = 2:M-1
        for j = 2:N-1
            L_x = L(i+1, j) - L(i, j);
            L_y = (L(i, j+1) - L(i, j-1))/2;
            co1 = 1/sqrt(0.00001^2 + L_x^2 + L_y^2);
            
            L_x = L(i, j) - L(i-1, j);
            L_y = (L(i-1, j+1) - L(i-1, j-1))/2;
            co2 = 1/sqrt(0.00001^2 + L_x^2 + L_y^2);
            
            L_x = (L(i+1, j) - L(i-1, j))/2;
            L_y = L(i, j+1) - L(i, j);
            co3 = 1/sqrt(0.00001^2 + L_x^2 + L_y^2);
            
            L_x = (L(i+1, j-1) - L(i-1, j-1))/2;
            L_y = L(i, j) - L(i, j-1);
            co4 = 1/sqrt(0.00001^2 + L_x^2 + L_y^2);
            
            co = co1 + co2 + co3 + co4;
%             L_x = (L(i+1, j-1) - L(i-1, j-1))/2;
%             L_y = (L(i-1, j+1) - L(i-1, j-1))/2;
%             grad = sqrt(L_x^2 + L_y^2);
            grad = 1;
            seg = 2*HeaviSide(Phi_old(i, j))*(HeaviSide(Phi_old(i, j))*HeaviSide(L(i, j)) - HeaviSide(Psi(i, j)));
            div = co1*L(i+1, j) + co2*L(i-1, j) + co3*L(i, j+1) + co4*L(i, j-1);
            L(i, j) = (L(i, j) + dt*delta(i, j)*grad*(lambda2*div + lambda1 - alpha*seg))/(1+dt*delta(i, j)*lambda2*co*grad);
        end
    end
    L(:, 1) = L(:, 2); L(:, N) = L(:, N-1);
    L(1, :) = L(2, :); L(M, :) = L(M-1, :);
    L(1, 1) = L(2, 2); L(1, N) = L(2, N-1);
    L(M, 1) = L(M-1, 2); L(M, N) = L(M-1, N-1);
    
    Iter = Iter + 1; 
end
end

function [Phi, E] = UpdateSegmentation(Phi_old, L, Psi)
global I0
global alpha lambda1 lambda2 nu
global C1 C2
temp = -((I0 - C1).^2 - (I0 - C2).^2 + 2*alpha*HeaviSide(L).*(HeaviSide(Phi_old).*HeaviSide(L) - HeaviSide(Psi)));
Phi = (temp > 0).*1 + (temp < 0).*-1;
E = (I0 - C1).^2.*HeaviSide(Phi) + (I0 - C2).^2.*(1-HeaviSide(Phi)) ...
    + alpha*(HeaviSide(Phi).*HeaviSide(L) - HeaviSide(Psi)).^2 ...
    + lambda1*(1 - HeaviSide(L)) + lambda2*DiracDelta(L) ...
    + nu*((I0 - C1).^2.*HeaviSide(Psi) + (I0 - C2).^2.*(1-HeaviSide(Psi)));
E = sum(E, 'all');
end