%% Figure 6 - Dispersion relations for Borek synthetic gene circuit
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

% Auxiliary function
DispersionRelation=@(q,J,D)max(real(eig(J-q^2*D)));

% Parameters
bR=0.0156;uR=2;uA=2;kA=2;bA=0.0117;aR=0.5;kdim=2;
a1 = 2142;a3 = 1190;kAHL = 2;kH2O2 = 0.057;b1 = 0.156;b3 = 0.03;b4 = 0.25;
bF = 2;dN = 2;dI = 2;dR = 2;g2 = 2;g3 = 2;k1 = 0.0156;
k2 = 2;kb = 2;kcat = 2;kdim = 2;kf = 0.0117;kluxR = 0.5;kun = 2;nF = 2;
nP = 2;u1 = 2;u3 = 2;u4 = 2;uF = 2;
DAHL=1; DH2O2=100;
qspace=logspace(-5,5,1000);

% Equilibria
AHLSS=1313.75;luxRAHLSS=2.5618;AHLluxRSS=2.5618;H2O2SS=3.03947;AiiaSS=188.41;AiiaAHLSS=725.241;...
luxRSS=0.25;P1SS=0.677162;P3SS=0.17924;P4SS=0.316655;luxISS=725.241;NdhSS=106.648;
FisSS=1.50489;DimerSS=6.56284;

% Full system (35)
JF=[-AiiaSS*kf-k1*luxRSS,k2,0,-AHLSS*kf,kb,-AHLSS*k1,0,0,0,kAHL,0,0,0;
    k1*luxRSS,-k2-4*kdim*luxRAHLSS,0,0,0,AHLSS*k1,0,0,0,0,0,0,2*kun;
    0,0,-g2-bF*(-FisSS+nF),0,0,0,0,0,0,0,kH2O2,bF*H2O2SS+uF,0;
    -AiiaSS*kf,0,0,-g3-AHLSS*kf,kb+kcat,0,0,0,a3,0,0,0,0;
    AiiaSS*kf,0,0,AHLSS*kf,-kb-kcat,0,0,0,0,0,0,0,0;
    0,0,0,0,0,-dR,0,0,0,0,0,0,0;
    0,0,0,0,0,0,-b1*DimerSS-u1,0,0,0,0,0,b1*(nP - P1SS);
    0,0,0,0,0,0,0,-b3*DimerSS-u3,0,0,0,0,b3*(nP-P3SS);
    0,0,0,0,0,0,0,0,-b4*FisSS-u4,0,0,b4*(nP-P4SS),0;
    0,0,0,0,0,0,a1,0,0,-dI,0,0,0;
    0,0,0,0,0,0,0,a3,0,0,-dN,0,0;
    0,0,bF*(-FisSS+nF),0,0,0,0,0,b4*FisSS+u4,0,0,-bF*H2O2SS-b4*(nP-P4SS)-uF,0;
    0,2*kdim*luxRAHLSS,0,0,0,0,b1*DimerSS+u1,b3*DimerSS+u3,0,0,0,0,-kun-b1*(nP-P1SS)-b3*(nP-P3SS)];

DF=zeros(13);DF(1,1)=DAHL;DF(3,3)=DH2O2;

% Intermediate system (36)
J4=[-(bR*aR/dR+bA*AiiaSS*kA/(uA+kA)),0,uR-2*a1*kAHL*b1*kdim*nP*b1*dI*kdim*AHLluxRSS^3/(dI*kun*u1+b1*dI*kdim*AHLluxRSS^2)^2+2*a1*kAHL*b1*kdim*nP*AHLluxRSS/(dI*kun*u1+b1*dI*kdim*AHLluxRSS^2),-bA*kA*AHLSS/(uA+kA);
    0,-g2,-2*b3*kdim*a3*kH2O2*nP*b3*dN*kdim*AHLluxRSS^3/(dN*kun*u3+b3*dN*kdim*AHLluxRSS^2)^2+2*b3*kdim*a3*kH2O2*nP*AHLluxRSS/(dN*kun*u3+b3*dN*kdim*AHLluxRSS^2),0;
    bR*aR/dR,0,-uR,0;
    0,-a3*b4*bF*nF*nP*(b4*bF*nF+u4*bF)*H2O2SS/(uF*u4+(b4*bF*nF+u4*bF)*H2O2SS)^2+a3*b4*bF*nF*nP/(uF*u4+(b4*bF*nF+u4*bF)*H2O2SS),0,-g3];
D4=diag([DAHL DH2O2 0 0]);

% Reduced system (37)
J2=[-2*a1*kAHL*b1*bR^2*kdim*aR^2*nP*b1*dI*bR^2*kdim*aR^2*AHLSS^3/(dI*dR^2*uR^2*kun*u1+b1*dI*bR^2*kdim*aR^2*AHLSS^2)^2+2*a1*kAHL*b1*bR^2*kdim*aR^2*nP*AHLSS/(dI*dR^2*uR^2*kun*u1+b1*dI*bR^2*kdim*aR^2*AHLSS^2)-nP*a3*b4*bF*H2O2SS*kA*bA*nF/(g3*(uA+kA)*(b4*bF*H2O2SS*nF+u4*(bF*H2O2SS+uF))),...
    -(-nP*a3*b4*bF*kA*bA*nF*AHLSS)*g3*(uA+kA)*(b4*bF*nF+u4*bF)*H2O2SS/(g3*(uA+kA)*u4*uF+g3*(uA+kA)*(b4*bF*nF+u4*bF)*H2O2SS)^2+(-nP*a3*b4*bF*kA*bA*nF*AHLSS)/(g3*(uA+kA)*u4*uF+g3*(uA+kA)*(b4*bF*nF+u4*bF)*H2O2SS);...
    -2*b3*bR^2*kdim*aR^2*a3*kH2O2*nP*dN*b3*bR^2*kdim*aR^2*AHLSS^3/(dN*dR^2*uR^2*kun*u3+dN*b3*bR^2*kdim*aR^2*AHLSS^2)^2+2*b3*bR^2*kdim*aR^2*a3*kH2O2*nP*AHLSS/(dN*dR^2*uR^2*kun*u3+dN*b3*bR^2*kdim*aR^2*AHLSS^2),...
    -g2];
D2=diag([DAHL DH2O2]);

%% Create plot
lw = 1;
left = 0.12;
bottom = 0.15;
width = 0.8;
height = 0.8;
xlims = [-4 4]; ylims = [-0.3 0.2];

f1 = figure(1);
f1.Position = [100 300 500 350];
clf;

semilogx(qspace,arrayfun(@(q)DispersionRelation(q,JF,DF),qspace),'LineWidth',lw);
hold on;
semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J4,D4),qspace),'LineWidth',lw);
semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J2,D2),qspace),'LineWidth',lw);
semilogx(qspace,0*qspace,'--k','LineWidth',lw);
hold off
axis([10.^xlims ylims]);
set(gca,'XTick',[1e-4,1e-2,1,1e2,1e4],'LineWidth',lw,'position',[left bottom width height],'color','none','Xscale','log')
xlabel('k');
ylabel('max \{Re \lambda(k)\}')
box off
legend({'Full system (35)','Intermediate system (36)','Reduced system (37)'},'Location','NorthEast','box','off');

%save2pdf('Figures/BorekDispersion',f1,300)

return