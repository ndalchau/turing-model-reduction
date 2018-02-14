%% Figure 1 - Bifurcation analysis for Brusselator
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

plotting = 1;   % Plot PDE simulation progress
compute = 1;    % Recompute everything or load from file

%% Run Bifurcation calculations
npoints = 300;
fnameBif = ['Matfiles/BrusselatorBifurcations' num2str(npoints) '.mat'];
if compute || ~exist(fnameBif,'file')
  a = 1;
  Dx = 1;
  qspace=logspace(-5,5,1000);
  bspace=logspace(log10(0.8),log10(10),npoints);
  Dyspace=logspace(-1,4,npoints);
  
  U=meshgrid(bspace,bspace);
  V=meshgrid(Dyspace,Dyspace)';
  
  disp('Running bifurcation analysis for Bimolecular Brusselator')
  tic
  bimol_raw = arrayfun(@(b,Dy)BimolecularBrusselator(a,b,Dx,Dy,qspace),U,V);
  bimol = reshape(bimol_raw,[npoints,npoints]);
  toc
  
  disp('Running bifurcation analysis for Classical Brusselator')
  tic
  class_raw = arrayfun(@(b,Dy)ClassicalBrusselator(a,b,Dx,Dy,qspace),U,V);
  class = reshape(class_raw,[npoints,npoints]);
  toc
  
  save(fnameBif,'qspace','bimol','class','bspace','Dyspace')
else
  load(fnameBif,'qspace','bimol','class','bspace','Dyspace')
end

%% Run PDE simulations
T=1000;
fnamePDE = ['Matfiles/BrusselatorPDE' num2str(T) '.mat'];
if compute || ~exist(fnamePDE,'file')
  disp('Simulating Bimolecular Brusselator')
  tic
  [x,pde3] = BrussPDE(T,3,plotting);
  toc
  disp('Simulating Bimolecular Brusselator')
  tic
  [~,pde2] = BrussPDE(T,2,plotting);
  toc
  
  save(fnamePDE,'pde3','pde2','x')
else
  load(fnamePDE,'pde3','pde2','x')
end

%% Create plot
Bifs = {bimol,class};
pdes = {pde3,pde2};
lw = 1;
greyscale = [0 0 0;1 1 1;0.75 0.75 0.75];
titles = {'Bimolecular Brusselator','Classical Brusselator'};

left = 0.1;
bottom = 0.1;
width = 0.36;
height = 0.36;
dy = 0.49;

dim = 225;
fwidth = 2*dim+100;
fheight = 2*dim+50;
dx = dim/fwidth + 0.025;

f1 = figure(1);
f1.Position = [100 100 fwidth fheight];

for i = 1:2
  ax(i) = subplot('position',[left+(i-1)*dx bottom+dy width height]);
  surf(bspace,Dyspace,Bifs{i})
  set(gca,'Xtick',[1 2 3 5 10],'Yscale','log','Ydir','reverse','Xscale','log','layer','top','LineWidth',lw,'tickdir','out');
  colormap(ax(i),greyscale)
  axis([0.8 10 1e-1 1e4])
  shading flat
  grid off
  view([0 90])
  title(titles{i})
  xlabel('b')
  if i==1
    ylabel('D_y')
  end
end
hold on
hs(1) = area([5 6],[1e5 1e6],1e5,'FaceColor',greyscale(2,:));
hs(2) = area([5 6],[1e5 1e6],1e5,'FaceColor',greyscale(1,:));
hs(3) = area([5 6],[1e5 1e6],1e5,'FaceColor',greyscale(3,:));
hold off
legend(hs,{'Stable','Patterns','Unstable'},'position',[0.85 0.75 0.12 0.05])

for i = 1:2
  subplot('position',[left+(i-1)*dx bottom width height]);
  surf(x,x,pdes{i});
  set(gca,'layer','top','LineWidth',lw,'tickdir','out')
  caxis([0 4])
  shading flat
  grid off
  view([0 90])
  xlabel('x coordinate')
  if i==1,ylabel('y coordinate'),end
end

labels = {'a','b','c','d','e','f'};
for i=1:2
  label(['{\bf' labels{i} '}'],[0.02+(i-1)*dx 2*dy-0.02 0.05 0.05]);
  label(['{\bf' labels{i+2} '}'],[0.02+(i-1)*dx dy-0.02 0.05 0.05]);
end

colorbar('position',[2*dx+0.06 bottom 0.03 height]);
%print('Figures/BrusselatorExample.png','-dpng','-r300')
print('-depsc2','-r600','Figures/BrusselatorExample')

return