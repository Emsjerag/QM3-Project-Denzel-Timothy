%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap using fast fourier transport (fft) method
%   In matlab, there is a internal function fft(...) which can do the job
%   For more detials on fft, type "help fft" in matlab command window and
%   press enter button to check the matlab help file  
% Unit of energy: hbar*omega, where h_bar is the Planck constant and
%   omega is the frequency of the trap 99
%   Unit of length: l=sqrt(h_bar/(m*omega)), where sqrt(...) is the square
%  mo
%   root function and m is the mass of the particle
%   Unit of momentum: hbar/l
%    energy unit: hbar\omega,  Hamiltonian --> dimensionless
%%   time dimensionless: omega*t    i d/dt | >= dimension H |>
%    dimensionless time = 2pi. one classical period
%--------------------------------------------------------------------------
a = -20;                       % Left end point 
b = +20;                       % Right end point 
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;       % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1];
A = 1; %amplitude paramneter
w = 0.5;
T= 0.269*pi;                         % Time duration of the evolution
M = 10^3;                     % Total No. of steps in the evolution
dt = T/M;                       % Time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UV = @(m) exp(-1i*(X.^2/2+ A*sin(X)*cos(m*w*dt))*dt/2);  

% One-step propagator in position space, only taking diagonal form
%UV = exp(-1i*(X.^2/2+0.1*X.^4)*dt/2);
UT = exp(-1i*(P.^2/2)*dt);       % One-setp propagator in momentum space
% note, hbar=1 in our dimensionless units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket located at X0
X0=4.0;
sigma=1.0;  % sigma is the width of the initial wavepacket
%psiprep=exp(-(X(1:N-1)-X0).^2/0.5)  squeezed
psiprep=exp(-(X(1:N)-X0).^2/(2*sigma^2));  %Gaussian state
psi=psiprep/sqrt(sum(abs(psiprep).^2));%normalized state

plot (X(1:N),abs(psi(1:N)).^2, 'g');   % plotting initial state

hold on 
%plot (P(1:N),abs(phi(1:N)).^2) 
psi_0=psi;

for m = 1:M
    psi_1 = UV(m).*psi_0;
    phi_2 = fft(psi_1);   %wavefunction in momentum space
    phi_3 = UT.*phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV(m).*psi_3;
    psi_0 = psi_4; %prepare a new cycle 
   
end
psi=psi_0; %final state updated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot (X(1:N),abs(psi(1:N)).^2,'r')  %plotting the final state profile


%Notes: When the initial width of the wave packet is sigma = 1, 
%then the shape of the wave packet will not change over time. This is known
%as the coherent state. it is the most "classical like" state