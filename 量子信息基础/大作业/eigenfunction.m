% eigenfunction.m

clear all

NN = 800;                                           % Number of points in the problem space.
hbar = 6.63e-34 / 2 / pi;                           % Reduced Plank's constant
melec = 9.1e-31;                                    % Mass of an electron
eV2J = 1.6e-19;                                     % Energy conversion factors
J2eV = 1 / eV2J;

del_x = 0.05e-9;                                    % The step size
dt = 0.5e-17;                                       % Time steps 
DX = del_x * 1e9;                                   % Step size in nm.
XX = (DX : DX : DX * NN);                           % Length in nm for plotting

chi0 = hbar ^ 2 / (2 * melec * del_x ^ 2);          % This is for the eigen calculation.

%-----Define the potential-------------------------------------------------
V = zeros(1,NN);
% V shaped potential
for n = 1 : NN
    %V(n) = eV2J*(.0005)*(abs(NN / 2 - n));         % Triangle potential
    V(n) = eV2J*(.0125)*(abs(NN / 2 - n)/20)^2;    % Parabolic potential
    %V(n) = eV2J*(.002)*(NN / 2 - n);               % With electric field
end
%--------------------------------------------------------------------------

figure
subplot(3,2,1);                                     % Plot the potential function
plot(XX,J2eV * V * 1000,'k');
set(gca,'fontsize',12);
ylabel('V (meV)');
xlabel('nm');
Umax = max(J2eV*V);
title('Potential');
%axis([ DX (DX*NN) 0 Umax ]);

% ----- Eigenvalue calculation --------------------------------------------
% Specify the Hamiltonian
H = zeros(NN,NN);
H(1,1) = 2 * chi0 + V(1);
H(1,2) = -1 * chi0;
for n = 2 : NN-1
    H(n,n-1) = -1 * chi0;
    H(n,n)  =  2 * chi0 + V(n);
    H(n,n+1) = -1 * chi0;
end
H(NN,NN-1) = -1 * chi0;
H(NN,NN) =  2 * chi0 + V(NN);
  
[phi,D] = eig(H);                                   % Put the eigenenergies in D and functions in phi
%--------------------------------------------------------------------------

%Plot the eigenfunctions
subplot(3,2,3);

plot(XX,phi(:,1),'k');                              % Plot the ground state
TT = ylabel('f_1','FontName','Arial','fontsize',12);
TT = text(5,0.04,sprintf('%7.4f meV',J2eV * D(1,1) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
title('Eigen');
subplot(3,2,4);

plot(XX,phi(:,2),'k');                              % Plot the first excited state
TT = ylabel('f_2','FontName','Arial','fontsize',12);
TT = text(5,.03,sprintf('%7.4f meV',J2eV * D(2,2) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
subplot(3,2,5);

plot(XX,phi(:,3),'k');                              % Plot the second excited state
TT = ylabel('f_3','FontName','Arial','fontsize',12);
TT = text(5,.03,sprintf('%7.4f meV',J2eV * D(3,3) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
xlabel('nm');
subplot(3,2,6);

plot(XX,phi(:,4),'k');                              % Plot the third excited state
TT = ylabel('f_4','FontName','Arial','fontsize',12);
TT = text(5,.03,sprintf('%7.4f meV',J2eV * D(4,4) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
xlabel('nm');

saveas(gcf,'eigen.png');