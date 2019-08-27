%% Figure 3 - Bifurcation analysis for Turing's model 
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

plotting = 1;   % Plot PDE simulation progress
compute = 1;    % Recompute everything or load from file

%% Run Bifurcation calculations
npoints = 200;
fnameBif = ['Matfiles/TuringBifurcations' num2str(npoints) '.mat'];
if compute || ~exist(fnameBif,'file')
  k1=2; k2=0.2; k3=0.01; k4=0.08; k5=0.04;
  k6space=logspace(-2,4,npoints);
  k7=2; nc=6; Dx=1;
  Dyspace=logspace(-4,1,npoints);
  qspace=logspace(-5,5,200);
  U = meshgrid(k6space,k6space);
  V = meshgrid(Dyspace,Dyspace)';
  
  disp('Running 4 species')
  raw4 = arrayfun(@(k6,Dy)Bifurcation_Turing4(k1,k2,k3,k4,k5,k6,k7,nc,Dx,Dy,qspace),U,V);
  Bif4 = reshape(raw4,npoints,npoints);  
  disp('Running 3 species')
  raw3 = arrayfun(@(k6,Dy)Bifurcation_Turing3(k1,k2,k3,k4,k5,k6,k7,nc,Dx,Dy,qspace),U,V);
  Bif3 = reshape(raw3,npoints,npoints);  
  disp('Running 2 species')
  raw2 = arrayfun(@(k6,Dy)Bifurcation_Turing2(k1,k2,k3,k4,k5,k6,k7,nc,Dx,Dy,qspace),U,V);
  Bif2 = reshape(raw2,npoints,npoints);
  
  mkdir Matfiles
  save(fnameBif,'qspace','k6space','Dyspace','Bif*')
else
  load(fnameBif)
end

%% Run PDE simulations
T=5000;
fnamePDE = ['Matfiles/TuringPDE' num2str(T) '.mat'];
if compute || ~exist(fnamePDE,'file')
  [x,pde4] = TuringPDE(T,4,plotting);
  [~,pde3] = TuringPDE(T,3,plotting);
  [~,pde2] = TuringPDE(T,2,plotting);
  mkdir Matfiles
  save(fnamePDE,'pde*','x')
else
  load(fnamePDE)
end

%% Create plot
Bifs = {Bif4,Bif3,Bif2};
pdes = {pde4,pde3,pde2};
lw = 1;
greyscale = [0 0 0;1 1 1;0.75 0.75 0.75];
titles = {'4 Species','3 Species','2 Species'};
xlims = [1e-2 1e4];
ylims = [1e-4 1e1];
clims = [0 1.5];

left = 0.075;
bottom = 0.1;
width = 0.24;
height = 0.36;
dy = 0.49;

dim = 225;
fwidth = 3*dim+100;
fheight = 2*dim+50;
dx = dim/fwidth;

f1 = figure(1);
f1.Position = [100 100 fwidth fheight];

for i = 1:3
  ax(i) = subplot('position',[left+(i-1)*dx bottom+dy width height]);
  area([5 6],[100 200],100,'FaceColor','w')
  hold on
  surf(k6space,Dyspace,Bifs{i})
  set(gca,'Xtick',10.^(-2:2:4),'Yscale','log','Ydir','reverse','Xscale','log','layer','top','LineWidth',lw,'tickdir','out');
  colormap(ax(i),greyscale)
  axis([xlims ylims])
  caxis([0 3])
  shading flat
  grid off
  box off
  view([0 90])
  title(titles{i})
  xlabel('k_6')
  if i==1
    ylabel('D_y')
  end
  if i==3
    hold on
    hs(1) = area([5 6],[100 200],100,'FaceColor',greyscale(2,:));
    hs(2) = area([5 6],[100 200],100,'FaceColor',greyscale(1,:));
    hs(3) = area([5 6],[100 200],100,'FaceColor',greyscale(3,:));
    hold off
    legend(hs,{'Stable','Patterns','Unstable'},'position',[0.87 0.75 0.12 0.05])
  end
end

for i = 1:3
  subplot('position',[left+(i-1)*dx bottom width height]);
  surf(x,x,pdes{i});
  set(gca,'layer','top','LineWidth',lw,'tickdir','out')
  axis([min(x) max(x) min(x) max(x)])
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
%print('Figures/TuringExample.png','-dpng','-r300')
print('-depsc2','-r600','Figures/TuringExample')

return