function [reshaped,eq] = BorekBifurcations(sys,kdimspace,DiffH2O2space,qspace,loadEq)
% BorekBifurcations
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

if nargin<5
  loadEq=0;
end

% Parameters
k1=0.01563;k2=2;kb=2;kcat=2;kf=0.01172;kluxR=0.5;
a1=2142;a3=1190;kAHL=2;kH2O2=0.057;b1=0.156;b3=0.03;b4=0.25;bF=2;dN=2;dI=2;
dR=2;g2=2;g3=2;bR=0.0156;uR=2;uA=2;kA=2;bA=0.0117;aR=0.5;kun=2;nF=2;nP=2;
u1=2;u3=2;u4=2;uF=2;DiffAHL=1;

  function x = BorekFullEq(x0,kdim)
    %log10(kdim)
    fprintf('.')
    
    fun=@(x) [kAHL*x(10)^2-kf*x(4)^2*x(1)^2+kb*x(5)^2-k1*x(6)^2*x(1)^2+k2*x(2)^2;...
      k1*x(6)^2*x(1)^2-k2*x(2)^2-2*kdim*x(2)^4+2*kun*x(13)^2;...
      kH2O2*x(11)^2-g2*x(3)^2-bF*(nF-x(12)^2)*x(3)^2+uF*x(12)^2;...
      a3*x(7)^2-g3*x(4)^2-kf*x(4)^2*x(1)^2+kb*x(5)^2+kcat*x(5)^2;...
      kf*x(4)^2*x(1)^2-kb*x(5)^2-kcat*x(5)^2;...
      kluxR-dR*x(6)^2;...
      b4*x(12)^2*(nP-x(7)^2)-u4*x(7)^2;...
      b3*x(13)^2*(nP-x(8)^2)-u3*x(8)^2;...
      b1*x(13)^2*(nP-x(9)^2)-u1*x(9)^2;...
      a1*x(9)^2-dI*x(10)^2;...
      a3*x(8)^2-dN*x(11)^2;...
      bF*(nF-x(12)^2)*x(3)^2-uF*x(12)^2-b4*x(12)^2*(nP-x(7)^2)+u4*x(7)^2;...
      kdim*x(2)^4-kun*x(13)^2-b3*x(13)^2*(nP-x(8)^2)+u3*x(8)^2-b1*x(13)^2*(nP-x(9)^2)+u1*x(9)^2];
    
    x=abs(fsolve(fun,x0,optimset('MaxFunEvals',1e5,'MaxIter',1e5,'Display','none')));
    x=x.^2;
  end

  function out = BorekFullClassify(x,kdim,DiffH2O2)
    AHL=x(1);
    luxRAHL=x(2);
    H2O2=x(3);
    Aiia=x(4);
    %AiikAHL=x(5);
    luxR=x(6);
    P4=x(7);
    P3=x(8);
    P1=x(9);
    %luxI=x(10);
    %Ndh=x(11);
    Fis=x(12);
    Dimer=x(13);
    
    J=[-kf*Aiia-k1*luxR,k2,0,-kf*AHL,kb,-k1*AHL,0,0,0,kAHL,0,0,0;
      k1*luxR,-k2-4*kdim*luxRAHL,0,0,0,k1*AHL,0,0,0,0,0,0,2*kun;
      0,0,-g2-bF*(nF-Fis),0,0,0,0,0,0,0,kH2O2,-bF*(-1)*H2O2+uF,0;
      -kf*Aiia,0,0,-g3-kf*AHL,kb+kcat,0,a3,0,0,0,0,0,0;
      kf*Aiia,0,0,kf*AHL,-kb-kcat,0,0,0,0,0,0,0,0;
      0,0,0,0,0,-dR,0,0,0,0,0,0,0;
      0,0,0,0,0,0,b4*Fis*(-1)-u4,0,0,0,0,b4*(nP-P4),0;
      0,0,0,0,0,0,0,b3*Dimer*(-1)-u3,0,0,0,0,b3*(nP-P3);
      0,0,0,0,0,0,0,0,b1*Dimer*(-1)-u1,0,0,0,b1*(nP-P1);
      0,0,0,0,0,0,0,0,a1,-dI,0,0,0;
      0,0,0,0,0,0,0,a3,0,0,-dN,0,0;
      0,0,bF*(nF-Fis),0,0,0,-b4*Fis*(-1)+u4,0,0,0,0,bF*(-1)*H2O2-uF-b4*(nP-P4),0;
      0,2*kdim*luxRAHL,0,0,0,0,0,-b3*Dimer*(-1)+u3,-b1*Dimer*(-1)+u1,0,0,0,-kun-b3*(nP-P3)-b1*(nP-P4)];

    D=zeros(13,13);D(1,1)=DiffAHL;D(3,3)=DiffH2O2;
    DispersionRelation=arrayfun(@(q)max(real(eig(J-q^2*D))),qspace);
    
    swd=(max(real(eig(J)))<0); %StableWithoutDiffusion
    na=(max(real(eig(J(3:end,3:end))))>0); %NoiseAmplifying
    pfsk=(max(DispersionRelation)>0);%PositiveForSomeK
    if swd&&(~na)&&pfsk
      out=0;%patterns
    elseif swd&&(~na)&&(~pfsk)
      out=1;%always stable
    elseif swd&&na
      out=2; %noise amplifying
    elseif ~swd
      out=3; %unstable
    end
  end

  function x = BorekFourEq(x0,kdim)
    fprintf('.')
    
    fun=@(x) [(uR+a1*kAHL*b1*kdim*x(3)^2*nP/(b1*dI*kdim*x(3)^4+dI*kun*u1))*x(3)^2-(bR*aR/dR+bA*x(4)^2*kA/(uA+kA))*x(1)^2;
      b3*kdim*x(3)^4*a3*kH2O2*nP/(dN*(b3*kdim*x(3)^4+kun*u3))-g2*x(2)^2;
      bR*aR*x(1)^2/dR-uR*x(3)^2;
      a3*b4*bF*x(2)^2*nF*nP/(b4*bF*x(2)^2*nF+u4*(bF*x(2)^2+uF))-g3*x(4)^2];
    
    x=abs(fsolve(fun,x0,optimset('MaxFunEvals',1e5,'MaxIter',1e5,'Display','none')));
    x=x.^2;
  end

  function out = BorekFourClassify(x,kdim,DiffH2O2)
    
    AHL=x(1);
    H2O2=x(2);
    AHLluxR=x(3);
    Aiia=x(4);
    
    J=[-(bR*aR/dR+bA*Aiia*kA/(uA+kA)),0,uR-2*a1*kAHL*b1*kdim*nP*b1*dI*kdim*AHLluxR^3/(dI*kun*u1+b1*dI*kdim*AHLluxR^2)^2+2*a1*kAHL*b1*kdim*nP*AHLluxR/(dI*kun*u1+b1*dI*kdim*AHLluxR^2),-bA*kA*AHL/(uA+kA);
      0,-g2,-2*b3*kdim*a3*kH2O2*nP*b3*dN*kdim*AHLluxR^3/(dN*kun*u3+b3*dN*kdim*AHLluxR^2)^2+2*b3*kdim*a3*kH2O2*nP*AHLluxR/(dN*kun*u3+b3*dN*kdim*AHLluxR^2),0;
      bR*aR/dR,0,-uR,0;
      0,-a3*b4*bF*nF*nP*(b4*bF*nF+u4*bF)*H2O2/(uF*u4+(b4*bF*nF+u4*bF)*H2O2)^2+a3*b4*bF*nF*nP/(uF*u4+(b4*bF*nF+u4*bF)*H2O2),0,-g3];
    
    D=zeros(4,4);D(1,1)=DiffAHL;D(2,2)=DiffH2O2;
    DispersionRelation=arrayfun(@(q)max(real(eig(J-q^2*D))),qspace);
    
    swd=(max(real(eig(J)))<0); %StableWithoutDiffusion
    na=(max(real(eig(J(3:end,3:end))))>0); %NoiseAmplifying
    pfsk=(max(DispersionRelation)>0);%PositiveForSomeK
    if swd&&(~na)&&pfsk
      out=0;%patterns
    elseif swd&&(~na)&&(~pfsk)
      out=1;%always stable
    elseif swd&&na
      out=2; %noise amplifying
    elseif ~swd
      out=3; %unstable
    end
  end

  function x = BorekTwoEq(kdim)
    
    fprintf('.')    
    fun=@(x) [x(1)^2*nP*(a1*kAHL*x(1)^2*b1*k1^2*kdim*kluxR^2/(x(1)^4*b1*dI*k1^2*kdim*kluxR^2+dI*dR^2*k2^2*kun*u1)-a3*b4*bF*x(2)^2*kcat*kf*nF/(g3*(kb+kcat)*(b4*bF*x(2)^2*nF+u4*(bF*x(2)^2+uF))));...
      (x(1)^4*b3*k1^2*kdim*kluxR^2*a3*kH2O2*nP)/(dN*(x(1)^4*b3*k1^2*kdim*kluxR^2+dR^2*k2^2*kun*u3))-g2*x(2)^2];
    
    x0=sqrt([1024;2]);    
    x=abs(fsolve(fun,x0,optimset('MaxFunEvals',1e08,'MaxIter',1e08,'Display','none')));
    x=x.^2;
  end
    
  function out = BorekTwoClassify(x,kdim,DiffH2O2)
    
    J=[2*nP*a1*kAHL*x(1)*b1*k1^2*kdim*kluxR^2*dI*dR^2*k2^2*kun*u1/(x(1)^2*b1*dI*k1^2*kdim*kluxR^2+dI*dR^2*k2^2*kun*u1)^2-nP*a3*b4*bF*x(2)*kcat*kf*nF/(g3*(kb+kcat)*(b4*bF*x(2)*nF+u4*(bF*x(2)+uF))),-nP*x(1)*a3*b4*bF*kcat*kf*nF*g3*(kb+kcat)*uF*u4/(g3*(kb+kcat)*(b4*bF*x(2)*nF+u4*(bF*x(2)+uF)))^2;
      2*(x(1)*b3*k1^2*kdim*kluxR^2*(a3*kH2O2*nP))*dR^2*k2^2*kun*u3*dN/(dN*(x(1)^2*b3*k1^2*kdim*kluxR^2+dR^2*k2^2*kun*u3))^2,-g2];
    
    D=zeros(2,2);D(1,1)=DiffAHL;D(2,2)=DiffH2O2;
    DispersionRelation=arrayfun(@(q)max(real(eig(J-q^2*D))),qspace);
    
    swd=(max(real(eig(J)))<0); %StableWithoutDiffusion
    na=false; %NoiseAmplifying
    pfsk=(max(DispersionRelation)>0);%PositiveForSomeK
    if swd&&(~na)&&pfsk
      out=0;%patterns
    elseif swd&&(~na)&&(~pfsk)
      out=1;%always stable
    elseif swd&&na
      out=2; %noise amplifying
    elseif ~swd
      out=3; %unstable
    end
  end

%% Load or compute the equilibria
if loadEq
  fprintf(['Loading equilibria for Borek ' sys])
  dat = csvread(['../Borek' sys 'Eq.csv']);
  kdimspace=10.^dat(:,1);
  npoints = length(kdimspace);
  eq=dat(:,2:end)';
else
  fprintf(['Computing equilibria for Borek ' sys])
  npoints = length(kdimspace);
  switch sys
    case 'Full'
      x0F = sqrt([1060;2.07;2.14;173;538;0.25;0.5;0.13;0.29;538;74.8;1.36;4.29]);
      for i = 1:npoints
        eq(:,i) = BorekFullEq(x0F,kdimspace(i));
        x0F = sqrt(eq(:,i));
      end
    case 'Four'
      x04=sqrt([508;0.319;0.991;67.9]);
      for i = 1:npoints
        eq(:,i) = BorekFourEq(x04,kdimspace(i));
        x04 = sqrt(eq(:,i));
      end
    case 'Two'
      for i = 1:npoints
        eq(:,i) = BorekTwoEq(kdimspace(i));
      end
  end
end
fprintf('\n')

%% Now classify each point, and as a function of DiffH2O2 as well
fprintf(['Classifying Borek ' sys])
switch sys
  case 'Full'
    func = @BorekFullClassify;
  case 'Four'
    func = @BorekFourClassify;
  case 'Two'
    func = @BorekTwoClassify;
end

tic
for i = 1:npoints
  fprintf('.')
  reshaped(:,i) = arrayfun(@(D)func(eq(:,i),kdimspace(i),D),DiffH2O2space);
end
fprintf('\n')
toc

end