function[x,X,Y,W,C]=TuringPDE(T,ns,plotting)

fprintf('Running Turing %d for %1.1f time units\n',ns,T)

if nargin<3
  plotting=0;
end

k1=2;k2=0.2;k3=0.01;k4=0.08;k5=0.04;k6=3.37;k7=2;nc=6;DX=1;DY=0.04;
%T=1000;
X0=2.3069;
Y0=2;
W0=46.4039;
C0=1.3730;
%dx=0.05;L=6;
dx=0.1;L=12;
dt=0.05*dx^2/DX;
s=0.01;
x=0:dx:L;
n=length(x);
X=X0*ones(n)+s*rand(n,n);
Y=Y0*ones(n)+s*rand(n,n);
W=W0*ones(n)+s*rand(n,n);
C=C0*ones(n)+s*rand(n,n);
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
  switch ns
    case 4
      dX=dt*(-k1*X.*Y-2*k3*X.^2+k4+k7*(nc-C)+DX*d2X);
      dY=dt*(-k1*X.*Y+2*k2*W-k5*Y-k6*Y.*C+DY*d2Y);
      dW=dt*(k1*X.*Y-k2*W+k3*X.^2);
      dC=dt*(-k6*Y.*C+k7*(nc-C));
      X=X+dX;
      Y=Y+dY;
      W=W+dW;
      C=C+dC;
      X(1,:)=X(2,:);X(n,:)=X(n-1,:);X(:,1)=X(:,2);X(:,n)=X(:,n-1);
      Y(1,:)=Y(2,:);Y(n,:)=Y(n-1,:);Y(:,1)=Y(:,2);Y(:,n)=Y(:,n-1);
      W(1,:)=W(2,:);W(n,:)=W(n-1,:);W(:,1)=W(:,2);W(:,n)=W(:,n-1);
      C(1,:)=C(2,:);C(n,:)=C(n-1,:);C(:,1)=C(:,2);C(:,n)=C(:,n-1);
    case 3
      dX=dt*(-k1*X.*Y-2*k3*X.^2+k4+k6*k7*nc*Y./(k7+k6*Y)+DX*d2X);
      dY=dt*(-k1*X.*Y+2*k2*W-k5*Y-k6*k7*nc*Y./(k7+k6*Y)+DY*d2Y);
      dW=dt*(k1*X.*Y-k2*W+k3*X.^2);
      X=X+dX;
      Y=Y+dY;
      W=W+dW;
      X(1,:)=X(2,:);X(n,:)=X(n-1,:);X(:,1)=X(:,2);X(:,n)=X(:,n-1);
      Y(1,:)=Y(2,:);Y(n,:)=Y(n-1,:);Y(:,1)=Y(:,2);Y(:,n)=Y(:,n-1);
      W(1,:)=W(2,:);W(n,:)=W(n-1,:);W(:,1)=W(:,2);W(:,n)=W(:,n-1);
    case 2
      dX=dt*(-k1*X.*Y-2*k3*X.^2+k4+k6*k7*nc*Y./(k7+k6*Y)+DX*d2X);
      dY=dt*(k1*X.*Y+2*k3*X.^2-k5*Y-k6*k7*nc*Y./(k7+k6*Y)+DY*d2Y);
      X=X+dX;
      Y=Y+dY;
      X(1,:)=X(2,:);X(n,:)=X(n-1,:);X(:,1)=X(:,2);X(:,n)=X(:,n-1);
      Y(1,:)=Y(2,:);Y(n,:)=Y(n-1,:);Y(:,1)=Y(:,2);Y(:,n)=Y(:,n-1);
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
      plot(ts,Xmax-Xmin)
      ylabel('max(X) - min(X)')
    end
  else
    if mod(i,1000)==0
      fprintf('.')
    end
  end  
end

end