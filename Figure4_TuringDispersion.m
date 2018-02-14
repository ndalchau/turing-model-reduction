%% Figure 4 - Dispersion relations for Turing's model
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

% Parameters
k1=2;k2=0.2;k3=0.01;k4=0.08;k5=0.04;k6=3.37;k7=2;
nc=6;
Dx=1;Dy=0.04;
qspace=logspace(-4,4,800);

% Auxiliary function
DispersionRelation=@(q,J,D)max(real(eig(J-q^2*D)));

% Equilibria
Xss=(-k1*k4^2*k6-k1*k4*k5*k7+sqrt(k4)*sqrt(k4*k6+k5*k7)*sqrt(k1^2*k4^2*k6+...
    8*k3*k4*k5^2*k6+k1^2*k4*k5*k7+8*k3*k5^3*k7+8*k3*k5^2*k6*k7*nc))/(4*(k3*...
    k4*k5*k6+k3*k5^2*k7));
Yss=k4/k5;Wss=k3*Xss^2/k2+k1*Xss*Yss/k2;Css=k7*nc/(k7+k6*Yss);

% 4 Species
J4=[-k1*Yss-4*k3*Xss,-k1*Xss,0,-k7;
    -k1*Yss,-k1*Xss-k5-k6*Css,2*k2,-k6*Yss;
    k1*Yss+2*k3*Xss,k1*Xss,-k2,0;
    0,-k6*Css,0,-k6*Yss-k7];
D4=diag([Dx Dy 0 0]);

% 3 Species
J3=[-k1*Yss-4*k3*Xss,-k1*Xss-k6^2*k7*nc*Yss/(k7+k6*Yss)^2+k6*k7*nc/(k7+k6*Yss),0;
    -k1*Yss,-k1*Xss-k5+k6^2*k7*nc*Yss/(k7+k6*Yss)^2-k6*k7*nc/(k7+k6*Yss),2*k2;
    k1*Yss+2*k3*Xss,k1*Xss,-k2];
D3=diag([Dx Dy 0]);

% 2 Species
J2=[-k1*Yss-4*k3*Xss,-k1*Xss-k6^2*k7*nc*Yss/(k7+k6*Yss)^2+k6*k7*nc/(k7+k6*Yss);
    k1*Yss+4*k3*Xss,k1*Xss-k5+k6^2*k7*nc*Yss/(k7+k6*Yss)^2-k6*k7*nc/(k7+k6*Yss)];
D2=diag([Dx Dy]);

%% Create plot
lw = 1;
boxcolor = 0.6*[1 1 1];
left = 0.12;
bottom = 0.15;
width = 0.8;
height = 0.8;
xlims = [-4 4]; ylims = [-1 2.5];

f1 = figure(1);
f1.Position = [100 300 500 350];
clf;

% Start by placing polygon connecting inset to box annotation
pl = -1.5; pr = -0.25;
il = 0.18; ib = 0.58; iw = 0.2; ih = 0.2;
patch(10.^[pl pr pr pl pl],[0.01 0.01 -0.01 -0.01 0.01],'w','EdgeColor',boxcolor)
patch(10.^[xlims(1)+(il-left)/width*diff(xlims) pl pr xlims(1)+(il+iw-left)/width*diff(xlims) xlims(1)+(il-left)/width*diff(xlims)],...
  [ylims(1)+(ib-bottom)/height*diff(ylims) 0.01 0.01 ylims(1)+(ib+ih-bottom)/height*diff(ylims) ylims(1)+(ib-bottom)/height*diff(ylims)],...
  0.9*[1 1 1],'EdgeColor','none');

% Now add the plot
hold on;
ph(1) = semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J4,D4),qspace),'LineWidth',lw);
ph(2) = semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J3,D3),qspace),'LineWidth',lw);
ph(3) = semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J2,D2),qspace),'LineWidth',lw);
ph2 = semilogx(qspace,arrayfun(@(q)DispersionRelation(q,J3,D3),qspace),'--','LineWidth',lw);
semilogx(qspace,0*qspace,'--k','LineWidth',lw);
hold off
axis([10.^xlims ylims]);
set(gca,'XTick',[1e-4,1e-2,1,1e2,1e4],'LineWidth',lw,'position',[left bottom width height],'color','none','Xscale','log')
xlabel('k');
ylabel('max \{Re \lambda(k)\}')
box off
legend(ph,{'4 Species','3 Species','2 Species'},'Location','NorthEast','box','off');

% Now plot the inset
qspace2 = logspace(-1,-0.5,50);
inset = axes('position',[il ib iw ih]);
semilogx(qspace2,arrayfun(@(q)DispersionRelation(q,J4,D4),qspace2),'LineWidth',lw);
hold on;
semilogx(qspace2,arrayfun(@(q)DispersionRelation(q,J3,D3),qspace2),'--','LineWidth',lw);
semilogx(qspace2,arrayfun(@(q)DispersionRelation(q,J2,D2),qspace2),'LineWidth',lw)
semilogx(qspace2,0*qspace2,'--k','LineWidth',lw);
hold off
axis([min(qspace2) max(qspace2) -0.01 0.01])
set(inset,'FontSize',8,'Xtick',[],'Ytick',[],'LineWidth',lw)

%save2pdf('Figures/TuringDispersion',f1,300)

return