%% Figure 2 - Dispersion relations for Brusselator
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

% Auxiliary function
DispersionRelation=@(q,J,D)max(real(eig(J-q^2*D)));

%Parameters
a=1;b=1.88;Dx=1;Dy=10;qspace=logspace(-5,5,1000);

%Equilibria
Xss=1;Yss=b/a;Zss=a^2; 

% Bimolecular Brusselator
J3 = [-(1+b)-4*Xss,Zss,2+Yss;
    b,-Zss,-Yss;
    2*Xss,0,-1];
D3 = diag([Dx Dy 0]);

% Classical Brusselator
J2 = [-(1+b)+2*Xss*Yss,Xss^2;
    b-2*Xss*Yss,-Xss^2];
D2 = diag([Dx Dy]);

%% Create plot
lw = 1;
left = 0.12;
bottom = 0.15;
width = 0.8;
height = 0.8;
xlims = [-4 4]; ylims = [-1.1 0.2];

f1 = figure(1);
f1.Position = [100 300 500 350];
semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J3,D3),qspace),'LineWidth',lw);
hold on;
semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J2,D2),qspace),qspace,0*qspace,'--k','LineWidth',lw);
hold off
axis([10.^xlims ylims]);
set(gca,'XTick',[1e-4,1e-2,1,1e2,1e4],'LineWidth',lw,'position',[left bottom width height],'color','none')
xlabel('k');
ylabel('max \{Re \lambda(k)\}')
box off
legend({'Bimolecular','Classical'},'location','east','box','off');

%save2pdf('Figures/BrusselatorDispersion',f1,300)

return