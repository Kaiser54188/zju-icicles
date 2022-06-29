NN = 90;        
hbar = 6.63e-34 / 2 / pi;                          
melec = 9.1e-31;                                    
eV2J = 1.6e-19;                                     
J2eV = 1 / eV2J;
a=10;
del_x = 0.05e-9;                                 
dt = 0.5e-17;                                       
DX = del_x * 1e9;                             
XX = (DX : DX : DX * NN);                          
chi0 = hbar ^ 2 / (2 * melec * del_x ^ 2);         
V = zeros(1,NN);
for n = 41 : 50
    V(n) = eV2J*(.0005)*(abs(NN / 2 ))*0.04;  
end
figure
subplot(5,1,1);                                    
plot(XX,J2eV * V * 1000,'k');
set(gca,'fontsize',12);
ylabel('V (meV)');
xlabel('nm');
Umax = max(J2eV*V);
title('Potential');
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
  
[phi,D] = eig(H);                                
subplot(5,1,2);

plot(XX,phi(:,1),'k');                             
TT = ylabel('f_1','FontName','Arial','fontsize',12);
TT = text(5,0.04,sprintf('%7.4f meV',J2eV * D(1,1) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
title('Eigen');
subplot(5,1,3);

plot(XX,phi(:,2),'k');                         
TT = ylabel('f_2','FontName','Arial','fontsize',12);
TT = text(5,.03,sprintf('%7.4f meV',J2eV * D(2,2) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
subplot(5,1,4);

plot(XX,phi(:,3),'k');                         
TT = ylabel('f_3','FontName','Arial','fontsize',12);
TT = text(5,.03,sprintf('%7.4f meV',J2eV * D(3,3) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
xlabel('nm');
subplot(5,1,5);

plot(XX,phi(:,4),'k');                            
TT = ylabel('f_4','FontName','Arial','fontsize',12);
TT = text(5,.03,sprintf('%7.4f meV',J2eV * D(4,4) * 1000));
set(TT,'fontsize',12);
set(gca,'fontsize',12);
xlabel('nm');