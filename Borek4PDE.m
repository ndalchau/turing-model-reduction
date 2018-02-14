function [x,AHL,H2O2]=Borek4PDE(T,plotting)
% Borek4PDE
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

a1	= 2142.00;a3	= 1190.00; aAHL	= 2.00;aH2O2	= 0.05714;b1	= 0.1563;
b3	= 0.03125; b4	= 0.25; bFis	= 2.00;decayNdh	= 2.00;decayluxI	= 2.00;decayluxR	= 2.00;
g2	= 2.00;g3	= 2.00;k1	= 0.01563;k2	= 2.00;kb	= 2.00;kcat	= 2.00; kdim	= 2.00;
kf	= 0.01172;kluxR	= 0.50;kundim	= 2.00; nFis	= 2.00;nP	= 2.00;u1	= 2.00;u3	= 2.00;u4	= 2.00;uFis	= 2.00;
DiffAHL=1; DiffH2O2=100;

% Solver and domain setup
dx=2;
dt=0.2*dx^2/max(DiffAHL,DiffH2O2);
L=200;
s=0.01;
x=0:dx:L;
n=length(x);

% Initial conditions
AHL0=1060.37;
luxRAHL0=2.07171;
H2O20=2.1367;
Aiia0=173.165;
AHL=AHL0*ones(n)+s*(rand(n,n)-1/2); %1
luxRAHL=luxRAHL0*ones(n)+s*(rand(n,n)-1/2); %2
H2O2=H2O20*ones(n)+s*(rand(n,n)-1/2);%3
Aiia=Aiia0*ones(n)+s*(rand(n,n)-1/2);%4

i=0;
t=0;

if plotting
  ts = [];
  Xmin = [];
  Xmax = [];
  fh = figure;
  fh.Position = [200 300 1000 400];
end

while t<T
  %t/T
  d2AHL=4*del2(AHL,dx);
  d2H2O2=4*del2(H2O2,dx);  
  dAHL=dt*(-AHL.*(Aiia*decayluxR*kcat*kf+k1*(kb+kcat)*kluxR)/(decayluxR*(kb+kcat))+luxRAHL.*(k2+a1*aAHL*b1*kdim*luxRAHL*nP./(b1*decayluxI*kdim*luxRAHL.^2+decayluxI*kundim*u1))+DiffAHL*d2AHL);
  dH2O2=dt*(b3*kdim*luxRAHL.^2*a3*aH2O2*nP./(decayNdh*(b3*kdim*luxRAHL.^2+kundim*u3))-g2*H2O2+DiffH2O2*d2H2O2);
  dluxRAHL=dt*(AHL*k1*kluxR/decayluxR-k2*luxRAHL);
  dAiia=dt*(a3*b4*bFis*H2O2*nFis*nP./(b4*bFis*H2O2*nFis+u4*(bFis*H2O2+uFis))-g3*Aiia);
    
  AHL=AHL+dAHL;
  luxRAHL=luxRAHL+dluxRAHL;
  H2O2=H2O2+dH2O2;
  Aiia=Aiia+dAiia;  
  
  AHL(1,:)=AHL(2,:);AHL(n,:)=AHL(n-1,:);AHL(:,1)=AHL(:,2);AHL(:,n)=AHL(:,n-1);
  luxRAHL(1,:)=luxRAHL(2,:);luxRAHL(n,:)=luxRAHL(n-1,:);luxRAHL(:,1)=luxRAHL(:,2);luxRAHL(:,n)=luxRAHL(:,n-1);
  H2O2(1,:)=H2O2(2,:);H2O2(n,:)=H2O2(n-1,:);H2O2(:,1)=H2O2(:,2);H2O2(:,n)=H2O2(:,n-1);
  Aiia(1,:)=Aiia(2,:);Aiia(n,:)=Aiia(n-1,:);Aiia(:,1)=Aiia(:,2);Aiia(:,n)=Aiia(:,n-1);
  
  t=t+dt;
  i=i+1;
  
  if plotting 
    if mod(i,1000)==0
      %t/T
      subplot(1,2,1)
      surf(x,x,AHL);
      set(gca,'layer','top','tickdir','out')
      shading flat
      grid off
      view([0 90])
      %imagesc(X)
      colorbar
      title(sprintf('t = %1.3f',t))
      drawnow;
      ts = [ts t];
      Xmin = [Xmin min(min(AHL))];
      Xmax = [Xmax max(max(AHL))];
      subplot(1,2,2)
      plot(ts,Xmax,ts,Xmin)
      ylabel('max(X) - min(X)')
    end
  else
    if mod(i,1000)==0
      fprintf('.')
    end
  end  
end

end