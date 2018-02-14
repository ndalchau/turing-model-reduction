function [x,AHL,H2O2] = BorekFullPDE(T,plotting)
% BorekFullPDE
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
AiiaAHL0=538.007;
luxR0=0.25;
P10=0.502341;
P30=0.125695;
P40=0.291034;
luxI0=538.007;
Ndh0=74.7883;
Fis0=1.36239;
Dimer0=4.29197;

AHL=AHL0*ones(n)+s*(rand(n,n)-1/2); %1
luxRAHL=luxRAHL0*ones(n)+s*(rand(n,n)-1/2); %2
H2O2=H2O20*ones(n)+s*(rand(n,n)-1/2);%3
Aiia=Aiia0*ones(n)+s*(rand(n,n)-1/2);%4
AiiaAHL=AiiaAHL0*ones(n)+s*(rand(n,n)-1/2);%5
luxR=luxR0*ones(n)+s*(rand(n,n)-1/2);%6
P4=P40*ones(n)+s*(rand(n,n)-1/2);%7
P3=P30*ones(n)+s*(rand(n,n)-1/2);%8
P1=P10*ones(n)+s*(rand(n,n)-1/2);%9
luxI=luxI0*ones(n)+s*(rand(n,n)-1/2);%10
Ndh=Ndh0*ones(n)+s*(rand(n,n)-1/2);%11
Fis=Fis0*ones(n)+s*(rand(n,n)-1/2);%12
Dimer=Dimer0*ones(n)+s*(rand(n,n)-1/2);%13

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
  
  dAHL=dt*(aAHL*luxI-kf*Aiia.*AHL+kb*AiiaAHL-k1*luxR.*AHL+k2*luxRAHL+DiffAHL*d2AHL);
  dluxRAHL=dt*(k1*luxR.*AHL-k2*luxRAHL-2*kdim*luxRAHL.^2+2*kundim*Dimer);
  dH2O2=dt*(aH2O2*Ndh-g2*H2O2-bFis*(nFis-Fis).*H2O2+uFis*Fis+DiffH2O2*d2H2O2);
  dAiia=dt*(a3*P4-g3*Aiia-kf*Aiia.*AHL+kb*AiiaAHL+kcat*AiiaAHL);
  dAiiaAHL=dt*(kf*Aiia.*AHL-kb*AiiaAHL-kcat*AiiaAHL);
  dluxR=dt*(kluxR-decayluxR*luxR);
  dP4=dt*(b4*Fis.*(nP-P4)-u4*P4);
  dP3=dt*(b3*Dimer.*(nP-P3)-u3*P3);
  dP1=dt*(b1*Dimer.*(nP-P1)-u1*P1);
  dluxI=dt*(a1*P1-decayluxI*luxI);
  dNdh=dt*(a3*P3-decayNdh*Ndh);
  dFis=dt*(bFis*(nFis-Fis).*H2O2-uFis*Fis-b4*Fis.*(nP-P4)+u4*P4);
  dDimer=dt*(kdim*luxRAHL.^2-kundim*Dimer-b3*Dimer.*(nP-P3)+u3*P3-b1*Dimer.*(nP-P1)+u1*P1);
  
  AHL=AHL+dAHL;
  luxRAHL=luxRAHL+dluxRAHL;
  H2O2=H2O2+dH2O2;
  Aiia=Aiia+dAiia;
  AiiaAHL=AiiaAHL+dAiiaAHL;
  luxR=luxR+dluxR;
  P4=P4+dP4;
  P3=P3+dP3;
  P1=P1+dP1;
  luxI=luxI+dluxI;
  Ndh=Ndh+dNdh;
  Fis=Fis+dFis;
  Dimer=Dimer+dDimer;  
  
  AHL(1,:)=AHL(2,:);AHL(n,:)=AHL(n-1,:);AHL(:,1)=AHL(:,2);AHL(:,n)=AHL(:,n-1);
  luxRAHL(1,:)=luxRAHL(2,:);luxRAHL(n,:)=luxRAHL(n-1,:);luxRAHL(:,1)=luxRAHL(:,2);luxRAHL(:,n)=luxRAHL(:,n-1);
  H2O2(1,:)=H2O2(2,:);H2O2(n,:)=H2O2(n-1,:);H2O2(:,1)=H2O2(:,2);H2O2(:,n)=H2O2(:,n-1);
  Aiia(1,:)=Aiia(2,:);Aiia(n,:)=Aiia(n-1,:);Aiia(:,1)=Aiia(:,2);Aiia(:,n)=Aiia(:,n-1);
  AiiaAHL(1,:)=AiiaAHL(2,:);AiiaAHL(n,:)=AiiaAHL(n-1,:);AiiaAHL(:,1)=AiiaAHL(:,2);AiiaAHL(:,n)=AiiaAHL(:,n-1);
  luxR(1,:)=luxR(2,:);luxR(n,:)=luxR(n-1,:);luxR(:,1)=luxR(:,2);luxR(:,n)=luxR(:,n-1);
  P4(1,:)=P4(2,:);P4(n,:)=P4(n-1,:);P4(:,1)=P4(:,2);P4(:,n)=P4(:,n-1);
  P3(1,:)=P3(2,:);P3(n,:)=P3(n-1,:);P3(:,1)=P3(:,2);P3(:,n)=P3(:,n-1);
  P1(1,:)=P1(2,:);P1(n,:)=P1(n-1,:);P1(:,1)=P1(:,2);P1(:,n)=P1(:,n-1);
  luxI(1,:)=luxI(2,:);luxI(n,:)=luxI(n-1,:);luxI(:,1)=luxI(:,2);luxI(:,n)=luxI(:,n-1);
  Ndh(1,:)=Ndh(2,:);Ndh(n,:)=Ndh(n-1,:);Ndh(:,1)=Ndh(:,2);Ndh(:,n)=Ndh(:,n-1);
  Fis(1,:)=Fis(2,:);Fis(n,:)=Fis(n-1,:);Fis(:,1)=Fis(:,2);Fis(:,n)=Fis(:,n-1);
  Dimer(1,:)=Dimer(2,:);Dimer(n,:)=Dimer(n-1,:);Dimer(:,1)=Dimer(:,2);Dimer(:,n)=Dimer(:,n-1);
  
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
      ylabel('Range of X')
    end
  else
    if mod(i,1000)==0
      fprintf('.')
    end
  end  
end

end