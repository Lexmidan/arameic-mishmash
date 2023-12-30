function so=voda
% resi PDR u_t + (u^2)_x = 0
global T t dt dx J x u u2 un B C a F ul up g

B = 20;
J = 200;
T = 2;
C = 0.5;

nlwlf = 4;

t=0;
nstep=1;
init;
du = u;
while t < T
  if(mod(nstep,nlwlf)==0)
    lfstep;
  else
    lwstep;
  end
  plot(x,un(:,1)/g,'-','Marker','x');
  hold on;
  plot(x,un(:,2)./un(:,1),'-r','Marker','x');
  hold off;

  %  axis([-B B 10 100]);
  drawnow;
  if nstep==1
    pause
  end
  nstep=nstep+1;
end

function ii=init
% inicializace, pocatecni podminky
global T dt dx J x u u2 un B C a F ul up g

dx = 2*B/J;
x = zeros(J,1);
for j=1:J
  x(j) = -B + j*dx;
end
F = zeros(J,2);
u = zeros(J,2);
u2 = zeros(J,2);
un = zeros(J,2);
g = 10;
%vell = -1; velp = 1; hl = 1; hp = 1; T =4;
%vell = 1; velp = -1; hl = 0.1; hp = 0.1; T =8;
%vell = 0; velp = 0; hl = 2; hp = 0.1; T=4;
vell = 0; velp = 0; hl = 0.1; hp = 2; T=4;

% jen 1 razova vlna
%vell = sqrt(6)/2; velp = 0; hl = 4; hp = 2; g=1; T=6;
% 1 vlna zredeni
%vell = -sqrt(6)/2; velp = 0; hl = 4; hp = 2; g=1; T=2;

ul = [g*hl  g*hl*vell];
up = [g*hp  g*hp*velp];

for j=1:J
  if x(j) > 0
      u(j,:) = up(:);
  else
      u(j,:) = ul(:);
  end
end
lambda1 = max(abs(vell + sqrt(g*hl)),abs(velp + sqrt(g*hp)));
lambda2 = max(abs(vell - sqrt(g*hl)),abs(velp - sqrt(g*hp)));
dt = C*dx/abs(max(lambda1,lambda2));

function lf=lfstep
% jeden casovy krok Lax-Friendrichsovo schema
global t T dt dx J x u u2 un B C a F
t = t + dt;
% prvni pulkrok
flux(u,J);
for j=1:J-1
  u2(j,:) = (u(j,:) + u(j+1,:))/2 - dt*(F(j+1,:) - F(j,:))/(2*dx);
end
% druhy pulkrok
flux(u2,J-1);
for j=2:J-1
  un(j,:) = (u2(j-1,:) + u2(j,:))/2 - dt*(F(j,:) - F(j-1,:))/(2*dx);
end

% okrajove podminky
un(1,:) = u(2,:);
un(J,:) = u(J-1,:);
u = un;

function lw=lwstep
% jeden casovy krok Lax-Wendroff schema
global t T dt dx J x u u2 un B C a F
t = t + dt;
% prvni pulkrok
flux(u,J);
for j=1:J-1
  u2(j,:) = (u(j,:) + u(j+1,:))/2 - dt*(F(j+1,:) - F(j,:))/(2*dx);
end
% druhy pulkrok
flux(u2,J-1);
for j=2:J-1
  un(j,:) = u(j,:) - dt*(F(j,:) - F(j-1,:))/dx;
end

% okrajove podminky
un(1,:) = u(2,:);
un(J,:) = u(J-1,:);
u = un;

function ff = flux(uu,JJ)
global F
for j=1:JJ
%  F(j) = uu(j)^2/2;
  F(j,1) = uu(j,2);
  F(j,2) = uu(j,2)^2/uu(j,1) + uu(j,1)^2/2;
end


function lf=lfstep1
% jeden casovy krok Lax-Friendrichsovo schema
global t T dt dx J x u u2 un B C a F
t = t + dt;
% jedenkrok
flux(u,J);
for j=2:J-1
  un(j,:) = (u(j-1,:) + u(j+1,:))/2 - dt*(F(j+1,:) - F(j-1,:))/(2*dx);
end
% okrajove podminky
un(1,:) = u(2,:);
un(J,:) = u(J-1,:);
u = un;


