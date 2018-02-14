function[out]=BimolecularBrusselator(a,b,Dx,Dy,qspace)
% BimolecularBrusselator
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

Xss=1;Yss=b/a;Zss=a^2;
J=[-(1+b)-4*Xss,Zss,2+Yss;
  b,-Zss,-Yss;
  2*Xss,0,-1];
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