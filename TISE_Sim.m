% 208643N - Kaushan G. Ranasinghe - 25 MAY 2021
%=========================================================================
% TISE Simulation using FDM with Central Difference - 1D
%=========================================================================

clear all
close all

% Simulation of a electron in various potential functions

%% ------Initiating the constants------------------------------------------

hbar = 1.055e-34;
m = 9.11e-31; % Rest mass of electron
eV = 1.602e-19;

%% ------Initiating the stepped lattice points-----------------------------

L = 10e-10; % Potential well width
z_min = -L; % -ve Z range
z_max = L; % +ve Z range
n = 1000; % Number of lattice points
z = linspace(z_min,z_max,n); % Vector of 'z' lattice points
d = z(2) - z(1); % Successive lattice spacing

%% ------Generating the potential energy matrix----------------------------

Wh = 10; % Well Height
Vo = Wh*eV; % Max finite potential
B1 = -L/2; % -ve well boundary
B2 = L/2; % +ve well boundary

% -------Generation of potential function----------------------------------

%Select potential well type

%Symmetric Wells

%V = ones(1,n)*Vo; V(z>= B1 & z<=B2) = 0; % Square Well
%V = Vo/(L^2)*(z.^2); % Harmonic Well
%V = ones(1,n)*Vo; V(z>= B1 & z<=B2) = Vo/(B2^2)*(z((z>= B1 & z<=B2)).^2); % Truncated Parabolic Well
%V((z<=B1)) = Vo/(B2-L)^2*((z(z<=B1)) + L).^2; V(z>= B1 & z<=B2) = 0;V((z>=B2)) = Vo/(L-B2)^2*((z(z>=B2)) - L).^2; % Quadratic Well

%Asymmetric Wells

%V = ones(1,n)*Vo; V(z>= B1 & z<=0) = 0; % Offset Square Well
%V = zeros(1,n)*Vo; V(z>= -0.5e-9 & z<=-0.35e-9) = Vo; % Offset Tunnel
%V = ones(1,n)*Vo; V(z>= B1 & z<=0) = 0; V(z>=0 & z<=B2) = 0.5*Vo; % Stepped Well

%One-Electron Atom
% 
% Wh = 0; % Well Height
% epsi=8.854e-12; % Permittivity of free space
% z_min = 1e-10; % Minimum electron radius
% z = linspace(z_min,z_max,n);
% 
% V = -((eV^2)./(4*pi*epsi*z)); % Potential function

% Series of Square Potential Wells 

Wd = 0.2e-9; % Well width
Ws = 0.1e-9; % Well spacing
V=ones(1,n)*Vo;
q = 7; % Number of Series Wells
z0 = z(:,1);
for ser = 1:q
   
   z1 = z0;
   z2 = z1+Wd;
   V(z1<=z & z<=z2) = 0;
   z0 = z2+Ws
    
end

%--------------------------------------------------------------------------

Vn = eye(n,n);
Vp = V'.*Vn;

%% ------Generating the kinetic energy matrix using FDM--------------------

K = eye(n,n); % Initiating the kinetic energy matrix
K = K*(-2); % Setting the main diagonal
for t = 1:n-1 % Setting the two flanking diagonals
   K(t,t+1) = 1;
   K(t+1,t) = 1;
end

%% ------Generating the Hamiltonian----------------------------------------

H = ((-(hbar^2)/(2*m*(d^2)))*K) + Vp ;

%% ------Generating the energy eigen values/states-------------------------

[U,EV] = eig(H); % Obtaining the wavefunctions and the eigen energies

P = U.^2; % Probability density of wavefunctions

E = diag(EV); % Accumulating the energy states to a column vector

Ec = E./(eV);

Ve = (1/eV).*V;

disp(E(1))
disp(E(2))

%% -------Plotting the Functions-------------------------------------------

%---------Plotting the allowable energy states-----------------------------

Ea = Ec(Ec<Wh)

figure
plot(z',Ve,'r','Linewidth',2);
hold on
title('Allowed Energy States','Fontsize',16)
ylabel('Potential, V(\psi)eV','Fontsize',13,'FontWeight','bold')
xlabel('Displacement / z','FontWeight','bold')

for r = 1:length(Ea)
   
    yline(Ea(r),'b--',sprintf('%.3f eV',Ea(r)),'Linewidth',1.5);
   
end
grid on

%---------Plotting the Wavefunction / Probabilities------------------------

for w = 1:length(Ea)
    
    count = w-1;
    figure
    title("Wavefunction \psi_n, for n = " +count,'Fontsize',16)
    xlabel('Displacement / z','FontWeight','bold')
    yyaxis right % Plot the potential function
    
    plot(z',Ve,'--','Linewidth',2) ;
    ylabel('Potential, V(\psi)eV','Fontsize',13,'FontWeight','bold')
    hold on
    grid on
    yyaxis left
    plot(z',(U(:,w)),'Linewidth',2);
    
end

for b = 1:length(Ea)
    
    count = b-1;
    figure
    title("Probability Density |\psi_n|^2, for n = " +count,'Fontsize',16)
    xlabel('Displacement / z','FontWeight','bold')
    yyaxis right % Plot the potential function
    
    plot(z',Ve,'--','Linewidth',2) ;
    ylabel('Potential, V(\psi)eV','Fontsize',13,'FontWeight','bold')
    hold on
    grid on
    yyaxis left
    plot(z',(P(:,b)),'Linewidth',2);
    
end
