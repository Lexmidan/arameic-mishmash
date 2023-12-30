function so=voda
% resi PDR u_t + (u^2)_x = 0
global T t dt dx J x u u2 un B C a F ul up g

B       = 1;
J       = 1600;
T       = 2;
C       = 0.5;
g       = 1.4; % gamma pro 1 atomový plyn
nlwlf   = 4;

t     = 0;
nstep =1;
init;
du    = u;
while t < T
  if(mod(nstep,nlwlf)==0)
    lfstep;
  else
    lwstep;
  end
  v = u(:,2)./un(:,1);
  p = (un(:,3) - un(:,1).*v.^2/2)*(g-1);
  plot(x,un(:,1),'-','Marker','x');
  
  hold on;
  plot(x,un(:,2)./un(:,1),'-r','Marker','x');
  plot(x,p, '-g','Marker', 'x')
  hold off;

  %  axis([-B B 10 100]);
  drawnow;
  if nstep==1
    pause
  end
  nstep=nstep+1;
end
lfstep;
lfstep;
  v = u(:,2)./un(:,1);
  p = (un(:,3) - un(:,1).*v.^2/2)*(g-1);
  plot(x,un(:,1),'-','Marker','x');
  
  hold on;
  plot(x,un(:,2)./un(:,1),'-r','Marker','x');
  plot(x,p, '-g','Marker', 'x')
  hold off;
function ii=init
% inicializace, pocatecni podminky
global T dt dx J x u u2 un B C a F ul up g

dx = B/J;
x  = zeros(J,1);
for j=1:J
  x(j) = -0 + j*dx;
end
F   = zeros(J,3);
u   = zeros(J,3);
u2  = zeros(J,3);
un  = zeros(J,3);
%g = 10;
%vell = -1; velp = 1; hl = 1; hp = 1; T =4;
%vell = 1; velp = -1; hl = 0.1; hp = 0.1; T =8;
%vell = 0; velp = 0; hl = 2; hp = 0.1; T=4;
%vell = 0; velp = 0; hl = 0.1; hp = 2; T=4;

% jen 1 razova vlna
%vell = sqrt(6)/2; velp = 0; hl = 4; hp = 2; g=1; T=6;
% 1 vlna zredeni
%vell = -sqrt(6)/2; velp = 0; hl = 4; hp = 2; g=1; T=2;


rho_l   = 1;
u_l     = 0.75;
p_l     = 1;

rho_r   = 0.125;
u_r     = 0;
p_r     = 0.1;

x0  = 0.3;
T   = 0.2

e_l = p_l/(g-1) /rho_l;
e_r = p_r/(g-1) /rho_r;

ul  = [rho_l rho_l*u_l rho_l*e_l+1/2 *rho_l *u_l^2 ];
up  = [rho_r rho_r*u_r rho_r*e_r+1/2 *rho_r *u_r^2 ];

for j=1:J
  if x(j) > x0
      u(j,:) = up(:);
  else
      u(j,:) = ul(:);
  end
end
c_l = sqrt(g*p_l/rho_l);
c_r = sqrt(g*p_r/rho_r);

lambda1 = max(abs(u_l + (c_l)),abs(u_r + c_r));
lambda2 = max(abs(u_l - (c_l)),abs(u_r - c_r));
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
%{
function ff = flux(uu,JJ)
global F
for j=1:JJ
%  F(j) = uu(j)^2/2;
  F(j,1) = uu(j,2);
  F(j,2) = uu(j,2)^2/uu(j,1) + uu(j,1)^2/2;
end
%}

function ff = flux(uu,JJ)
global F g
for j=1:JJ
  rho   = uu(j,1);
  rychl = uu(j,2)/uu(j,1);
  e     = uu(j,3)/uu(j,1) - 0.5* rychl^2 *rho;
  p     = e *(g-1) *rho;  
  
  F(j,1) = uu(j,2);
  F(j,2) = rychl^2 *rho + p ;
  F(j,3) = rychl*(uu(j,3) + p);

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


