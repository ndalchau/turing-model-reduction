%% Figure 5 - Bifurcation analysis for Borek synthetic gene circuit
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

plotting = 1;   % Plot PDE simulation progress
compute = 1;    % Recompute everything or load from file

%% Run Bifurcation calculations
% The equilibria are computed using Mathematica. 
% Run the Borek sections of the analysis notebook to regenerate the files.
npoints = 201;
kdimspace = logspace(-1,1,npoints);
DiffH2O2space = logspace(0,3,npoints);
qspace = logspace(-5,5,101);
fname = sprintf('Matfiles/BorekBifurcations%d.mat',npoints);
if (compute || ~exist(fname,'file'))
  [borekF,eqF] = BorekBifurcations('Full',kdimspace,DiffH2O2space,qspace,1);
  [borek4,eq4] = BorekBifurcations('Four',kdimspace,DiffH2O2space,qspace);
  [borek2,eq2] = BorekBifurcations('Two',kdimspace,DiffH2O2space,qspace,1);
  mkdir Matfiles
  save(fname,'borek*','*space','eq*')
else
  load(fname)
end

%% Run PDE simulations
T = 2000;
fnamePDE = sprintf('Matfiles/BorekPDE%d.mat',T);
if compute || ~exist(fnamePDE,'file')
  [x,pdeF] = BorekFullPDE(T,plotting);
  [~,pde4] = Borek4PDE(T,plotting);
  [~,pde2] = Borek2PDE(T,plotting);
  
  mkdir Matfiles
  save(fnamePDE,'pde*','x')
else
  load(fnamePDE)
end

%% Create plot
Bifs = {borekF,borek4,borek2};
pdes = {pdeF,pde4,pde2};
lw = 1;
greyscale = [0 0 0;1 1 1;0.75 0.75 0.75];
titles = {'Full system','Intermediate system','Reduced system'};

left = 0.075;
bottom = 0.1;
width = 0.24;
height = 0.36;
dy = 0.49;
xlims = [1e-1 1e1];
ylims = [1e0 1e3];
clims = [0 2510];

dim = 225;
fwidth = 3*dim+100;
fheight = 2*dim+50;
dx = dim/fwidth;

f1 = figure;
f1.Position = [100 100 fwidth fheight];

for i = 1:3
  ax(i) = subplot('position',[left+(i-1)*dx bottom+dy width height]);
  surf(kdimspace,DiffH2O2space,Bifs{i})
  set(gca,'Yscale','log','Ydir','reverse','Xscale','log','layer','top','LineWidth',lw,'tickdir','out');
  colormap(ax(i),greyscale)
  axis([xlims ylims])
  caxis([0 3])
  shading flat
  grid off
  view([0 90])
  box off
  title(titles{i})
  xlabel('k_{dim}')
  if i==1
    ylabel('D_{H2O2}')
    hold on
    hs1 = area([5 6],[1e4 1e5],1e4,'FaceColor',greyscale(2,:));
    hs2 = area([5 6],[1e4 1e5],1e4,'FaceColor',greyscale(1,:));
    hold off
    legend([hs1 hs2],{'Stable','Patterns'},'location','northwest','box','off')
  end
end

for i = 1:3
  subplot('position',[left+(i-1)*dx bottom width height]);
  surf(x,x,pdes{i});
  set(gca,'layer','top','LineWidth',lw,'tickdir','out')
  caxis(clims)
  shading flat
  grid off
  view([0 90])
  xlabel('x coordinate')
  if i==1,ylabel('y coordinate'),end
end

labels = {'a','b','c','d','e','f'};
for i=1:3
  label(['{\bf' labels{i} '}'],[0.015+(i-1)*dx 2*dy-0.02 0.05 0.05]);
  label(['{\bf' labels{i+3} '}'],[0.015+(i-1)*dx dy-0.02 0.05 0.05]);
end

colorbar('position',[3*dx+0.04 bottom 0.02 height]);
%print('Figures/BorekExample.png','-dpng','-r300')
print('-depsc2','-r600','Figures/BorekExample')

return