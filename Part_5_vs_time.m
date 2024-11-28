%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap using fast fourier transport (fft) method
%   In matlab, there is a internal function fft(...) which can do the job
%   For more detials on fft, type "help fft" in matlab command window and
%   press enter button to check the matlab help file
%   Only need some basics on split-operator method and Fourier transform:
%   All quantities are in dimensionless unit 
%% 
%   Unit of energy: hbar*omega, where h_bar is the Planck constant and
%   omega is the frequency of the trap
%   Unit of length: l=sqrt(h_bar/(m*omega)), where sqrt(...) is the square
%   root function and m is the mass of the particle
%   Unit of momentum: hbar/l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% capital or small letters are different!
clear all; clc; tic
%--------------------------------------------------------------------------
a = -80;                        % Left end point of the trap
b = +80;                        % Right end point of the trap
L = b-a;                        % Width of the trap
N = 2048;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T = 1200;                       % Time duration of the evolution, can adjust this to see the difference between fast and slow move
dt = 0.02;
M = T/dt
% T = 2000;                         % Time duration of the evolution, can adjust this to see the difference between fast and slow move
% M = 12*10^4;                     % Total No. of steps in the evolution
% dt = T/M;                       % Time step
Binsize=25;                     % this is to control every Binsize steps, we take snapshots
A = 0.001;
w = 0.499;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    % One-step propagator in position space
UT = exp(-1i*(P.^2/2)*dt);      % One-setp propagator in momentum space

UV = @(m) exp(-1i*(X.^2/2+ A*sin(X)*cos(m*w*dt))*dt/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket with wave vector K0, centered at X0 and width DEL0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K0 = 0;  % Wavevector of the Gaussian  %to demonstrate FFT-p, use K0=0, for dynamics, K0=3
X0 = 0;     % Center of the Gaussian 
DEL0 = 1;  % Width of the Gaussian
%--------------------------------------------------------------------------
%   Un-normalized initial state
%try=hermiteH(0,1)
Poly_g = hermiteH(0,X);
Poly_e=hermiteH(1,X);

VE_INI_temp_g = Poly_g.*exp(-(X-X0).^2/(2*DEL0^2));%1i means "i"
%   Normalized initial state as ground state or excited of harmonic
%   oscillators
VE_INI_temp_e = Poly_e.*exp(-(X-X0).^2/(2*DEL0^2));
VE_INI_g = VE_INI_temp_g/sqrt(VE_INI_temp_g*VE_INI_temp_g'); %normalization
VE_INI_e = VE_INI_temp_e/sqrt(VE_INI_temp_e*VE_INI_temp_e');

psi_0 = VE_INI_g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the nth eigenfunction
herm_poly = @(n) hermiteH(n,X);
eigenfxn_temp = @(n) herm_poly(n).*exp(-(X-X0).^2/(2*DEL0^2));

%normalize
eigenfxn = @(n) eigenfxn_temp(n)/norm(eigenfxn_temp(n));

%define the potential
pot = diag(A/2*sin(X - X0));

%define the matrix elements Vab
V = @(a,b) eigenfxn(a)*pot*transpose(eigenfxn(b));

%define the energy difference between two energy levels
omega = @(a,b) (a-b);

%define the summation. m=0 and n=1.
k_max = 15;
coeff = 0;
for k = 0:k_max
    coeff = coeff + V(1,k)*V(k,0)/(omega(k,0)-w);
end
% Plot theoretical 2nd order perturbation results
% t_grid_t = linspace(0,(1-2*w)/(2*pi)*T,M);
% theoreticalP = abs(coeff)^2*4*sin(t_grid_t*pi).^2/(1-2*w)^2;
% plot(t_grid_t, theoreticalP, 'b')
% hold on

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot theoretical 1st order perturbation results with RWA
t_grid = linspace(0,(1-w)/(2*pi)*T,M);
t_gridx = linspace(0,T,M);

sin_arr = diag(A/2*sin(X - X0));
temp = sin_arr*VE_INI_g';
V_10 = dot(VE_INI_e,temp);
% 
theoreticalP = abs(V_10)^2*4*sin(t_grid*pi).^2/(1-w)^2;
% plot(t_grid, theoreticalP, 'b')
% hold on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot theoretical 1st order perturbation results without RWA
no_RWA_prob = zeros(1,M);
for t_ind = 1:M
    t = t_gridx(t_ind);
    term_1 = (exp(1i*(1+w)*t)-1)/(1+ w);
    term_2 = (exp(1i*(1-w)*t)-1)/(1-w);
    no_RWA_prob(t_ind) = abs(V_10*(term_1 + term_2))^2;
end
plot(t_grid, no_RWA_prob, 'g')
hold on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Plot numerical
trans_prob = zeros(1,M);
for m = 1:M
    psi_1 = UV(m).*psi_0;
    phi_2 = fft(psi_1);     %wavefunction in momentum space
    phi_3 = UT.*phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV(m).*psi_3;
    psi_0 = psi_4;          %prepare a new cycle
    temp = abs(dot(VE_INI_e,psi_0))^2;
    trans_prob(m) = temp;
end

trans_prob;
plot(t_grid,trans_prob,'r')
xlabel('$\frac{\omega_{10}-\omega}{2\pi} t_{i}$','Interpreter','latex','FontSize',20)
% xlabel('$\frac{\omega_{10}-2\omega}{2\pi} t_{i}$','Interpreter','latex','FontSize',20)
ylabel('${P_{1\leftarrow 0}}$','Interpreter','latex','FontSize',20)
% legend('Theory (RWA)', 'Theory (no RWA)','Numerical')
legend('Theory (no RWA)','Numerical')
% legend('Theory')
title(['$\omega  = $ ', num2str(w)],'Interpreter','latex','FontSize',20)