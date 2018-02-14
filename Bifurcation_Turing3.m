function out = Bifurcation_Turing3(k1,k2,k3,k4,k5,k6,k7,nc,Dx,Dy,qspace)
% Bifurcation_Turing3
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

Xss=max([(-k1*k4^2*k6-k1*k4*k5*k7+sqrt(k4)*sqrt(k4*k6+k5*k7)*sqrt(k1^2*k4^2*k6+...
  8*k3*k4*k5^2*k6+k1^2*k4*k5*k7+8*k3*k5^3*k7+8*k3*k5^2*k6*k7*nc))/(4*(k3*...
  k4*k5*k6+k3*k5^2*k7)),(-k1*k4^2*k6-k1*k4*k5*k7-sqrt(k4)*sqrt(k4*k6+k5*k7)*sqrt(k1^2*k4^2*k6+...
  8*k3*k4*k5^2*k6+k1^2*k4*k5*k7+8*k3*k5^3*k7+8*k3*k5^2*k6*k7*nc))/(4*(k3*...
  k4*k5*k6+k3*k5^2*k7))]);
Yss=k4/k5;Wss=k3*Xss^2/k2+k1*Xss*Yss/k2;Css=k7*nc/(k7+k6*Yss);


J=[-k1*Yss-4*k3*Xss,-k1*Xss-k6^2*k7*nc*Yss/(k7+k6*Yss)^2+k6*k7*nc/(k7+k6*Yss),0;
  -k1*Yss,-k1*Xss-k5+k6^2*k7*nc*Yss/(k7+k6*Yss)^2-k6*k7*nc/(k7+k6*Yss),2*k2;
  k1*Yss+2*k3*Xss,k1*Xss,-k2];

D=[Dx,0,0;0,Dy,0;0,0,0];


DispersionRelation=arrayfun(@(q)max(real(eig(J-q^2*D))),qspace);
swd=(max(real(eig(J)))<0); %StableWithoutDiffusion
na=(max(real(eig(J(3,3))))>0); %NoiseAmplifying
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