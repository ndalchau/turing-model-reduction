function [x,X] = BrussPDE(T,ns,plotting)
% BrussPDE
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

if nargin<3
  plotting=0;
end

%filename = 'testAnimated-';
a=1;
b=1.88;
DX=1;
DY=10;
%T=100;
X0=a;
Y0=b/a;
Z0=a^2;
dx=0.4;
dt=0.2*dx^2/DY;
L=40;
s=0.01;
x=0:dx:L;
n=length(x);
X=X0*ones(n)+s*rand(n,n);
Y=Y0*ones(n)+s*rand(n,n);
Z=Z0*ones(n)+s*rand(n,n);
t=0;
i=0;

if plotting
  ts = [];
  Xmin = [];
  Xmax = [];
  fh = figure;
  fh.Position = [200 300 1000 400];
end

while t<T
  d2X=4*del2(X,dx);
  d2Y=4*del2(Y,dx);
  if ns==2
    dX=dt*(a-X-b*X+X.*X.*Y+DX*d2X);
    dY=dt*(b*X-X.*X.*Y+DY*d2Y);
    X=X+dX;
    Y=Y+dY;
    X(1,:)=X(2,:);X(n,:)=X(n-1,:);X(:,1)=X(:,2);X(:,n)=X(:,n-1);
    Y(1,:)=Y(2,:);Y(n,:)=Y(n-1,:);Y(:,1)=Y(:,2);Y(:,n)=Y(:,n-1);
  elseif ns==3
    dX=dt*(a-X-b*X+Z.*Y+DX*d2X);
    dY=dt*(b*X-Z.*Y+DY*d2Y);
    dZ=dt*(X.^2-Z);
    X=X+dX;
    Y=Y+dY;
    Z=Z+dZ;
    X(1,:)=X(2,:);X(n,:)=X(n-1,:);X(:,1)=X(:,2);X(:,n)=X(:,n-1);
    Y(1,:)=Y(2,:);Y(n,:)=Y(n-1,:);Y(:,1)=Y(:,2);Y(:,n)=Y(:,n-1);
    Z(1,:)=Z(2,:);Z(n,:)=Z(n-1,:);Z(:,1)=Z(:,2);Z(:,n)=Z(:,n-1);
  end
  
  
  t=t+dt;
  i=i+1;
  
  if plotting 
    if mod(i,1000)==0
      %t/T
      subplot(1,2,1)
      surf(x,x,X);
      set(gca,'layer','top','tickdir','out')
      shading flat
      grid off
      view([0 90])
      %imagesc(X)
      colorbar
      title(sprintf('t = %1.3f',t))
      drawnow;
      ts = [ts t];
      Xmin = [Xmin min(min(X))];
      Xmax = [Xmax max(max(X))];
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