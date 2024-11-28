%attempted changes: 
% reduce L to make it faster: convergence check is by comparing two similar
% values of L around that scale and seeing if the wavefunction changes by
% much

%increase T significantly: for RWA to be valid, we need to be close to
%resonance. for the tr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% capital or small letters are different!
clear all; clc; tic
%--------------------------------------------------------------------------
a = -80;                        % Left end point of the trap
b = +80;                        % Right end point of the trap
L = b-a;                        % Width of the trap
N = 1024;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T = 600;                       % Time duration of the evolution, can adjust this to see the difference between fast and slow move
dt = 0.02;
M = T/dt
% T = 2000;                         % Time duration of the evolution, can adjust this to see the difference between fast and slow move
% M = 12*10^4;                     % Total No. of steps in the evolution
% dt = T/M; 
Binsize=25;                     % this is to control every Binsize steps, we take snapshots
A = 0.001;
wwidth = 4*pi/T;
warray = linspace(0.5-wwidth,0.5+wwidth,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    % One-step propagator in position space
UT = exp(-1i*(P.^2/2)*dt);      % One-setp propagator in momentum space

UV = @(m,w) exp(-1i*(X.^2/2+ A*sin(X)*cos(m*w*dt))*dt/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket with wave vector K0, centered at X0 and width DEL0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K0 = 0;  % Wavevector of the Gaussian  %to demonstrate FFT-p, use K0=0, for dynamics, K0=3
X0 = 0;     % Center of the Gaussian 
DEL0 = 1;  % Width of the Gaussian

Poly_g = hermiteH(0,X);
Poly_e=hermiteH(1,X);

VE_INI_temp_g = Poly_g.*exp(-(X-X0).^2/(2*DEL0^2));%1i means "i"
%   Normalized initial state as ground state or excited of harmonic
%   oscillators
VE_INI_temp_e = Poly_e.*exp(-(X-X0).^2/(2*DEL0^2));
VE_INI_g = VE_INI_temp_g/sqrt(VE_INI_temp_g*VE_INI_temp_g'); %normalization
VE_INI_e = VE_INI_temp_e/sqrt(VE_INI_temp_e*VE_INI_temp_e');

psi_0 = VE_INI_g;

% w_grid = linspace(0,(1-w)/(2*pi)*2000,M);
% Plot theoretical 1st order perturbation results

sin_arr = diag(A/2*sin(X - X0));
temp = sin_arr * VE_INI_g';
V_10 = dot(VE_INI_e,temp);

theoreticalP = zeros(1,length(warray));
for w_ind = 1:length(warray)
    w = warray(w_ind);
    theoreticalP(w_ind) = abs(V_10)^2*4*sin((1-w)/(2/T))^2/(1-w)^2;
end
plot((1-warray)/(2*pi/T), theoreticalP, 'b')
hold on

% Plot numerical 1st order perturbation
trans_prob = zeros(1,length(warray));
for w_ind = 1:length(warray)
    w = warray(w_ind);
    psi_0 = VE_INI_g;
    for m = 1:M
        psi_1 = UV(m,w).*psi_0;
        phi_2 = fft(psi_1);   %wavefunction in momentum space
        phi_3 = UT.*phi_2;
        psi_3 = ifft(phi_3);
        psi_4 = UV(m,w).*psi_3;
        psi_0 = psi_4; %prepare a new cycle 
        temp = abs(dot(VE_INI_e,psi_0))^2;
    end
    trans_prob(w_ind) = temp;
end

%x axis is w-wmn, divided by 2pi/T
plot((1-1*warray)/(2*pi/T),trans_prob,'r')
xlabel('$\frac{\omega_{10}-\omega_{i}}{2\pi} T$','Interpreter','latex','FontSize',20)
ylabel('${P_{1\leftarrow 0}}$','Interpreter','latex','FontSize',20)
hold on
Amp = abs(V_10)^2*T^2;
Amp
yline(Amp);
legend('Theory','Numerical','$|\langle\psi_{1}^{0}|\frac{A}{2}sin(x)|\psi_{0}^{0}\rangle|^{2} T^{2}$','Interpreter','latex','FontSize',12)
