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
warray = linspace(2-wwidth,2+wwidth,40);
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
Poly_2e = hermiteH(2,X);
Poly_4e = hermiteH(4,X);

VE_INI_temp_g = Poly_g.*exp(-(X-X0).^2/(2*DEL0^2));%1i means "i"
%   Normalized initial state as ground state or excited of harmonic
%   oscillators
VE_INI_temp_2e = Poly_2e.*exp(-(X-X0).^2/(2*DEL0^2));
VE_INI_g = VE_INI_temp_g/sqrt(VE_INI_temp_g*VE_INI_temp_g'); %normalization
VE_INI_2e = VE_INI_temp_2e/sqrt(VE_INI_temp_2e*VE_INI_temp_2e');

VE_INI_temp_4e = Poly_4e.*exp(-(X-X0).^2/(2*DEL0^2));
VE_INI_4e = VE_INI_temp_4e/sqrt(VE_INI_temp_4e*VE_INI_temp_4e');

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%define the summation. m=0 and n=4.
k_max = 10;
theoreticalcoeff = zeros(1,length(warray));

% Plot theoretical & numerical 2nd order perturbation
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
        temp = abs(dot(VE_INI_4e,psi_0))^2;
    end
    trans_prob(w_ind) = temp;
    coeff = 0;
    for k = 1:2:k_max
        coeff = coeff + V(4,k)*V(k,0)/(omega(k,0)-w);
    end
    theoreticalcoeff(w_ind) = abs(coeff)^2*4*sin((4-2*w)/2*T).^2/(4-2*w)^2;
end
tgrid = (4-2*warray)/(2*pi/T);
plot(tgrid, theoreticalcoeff, 'b')
hold on

%x axis is w-wmn, divided by 2pi/T
plot(tgrid,trans_prob,'r')
xlabel('$\frac{\omega_{40}-2\omega_{i}}{2\pi} T$','Interpreter','latex','FontSize',20)
ylabel('${P_{4\leftarrow 0}}$','Interpreter','latex','FontSize',20)
coeff = 0;
for k = 1:2:k_max
    coeff = coeff + V(4,k)*V(k,0)/(omega(k,0)-2);
end
% coeff
Amp = abs(coeff)^2*T^2;
Amp;
yline(Amp);
% legend('Theory','Numerical','$|\Sigma_{k} \frac{\langle\psi_{4}^{0}|\frac{A}{2}sin(x)|\psi_{k}^{0}\rangle \langle\psi_{k}^{0}|\frac{A}{2}sin(x)|\psi_{0}^{0}\rangle}{\omega_{k0}-\omega_{40}/2}|^{2} T^{2}$','Interpreter','latex')
legend('Theory','Numerical','$|\Sigma_{k} \frac{V_{4k}V_{k0}}{\omega_{k0}-\omega_{40}/2}|^{2} T^{2}$','Interpreter','latex','FontSize',11)
% title(['$\omega  = $ ', num2str(w)],'Interpreter','latex')