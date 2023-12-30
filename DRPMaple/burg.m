function so=burg
% resi PDR u_t + (u^2)_x = 0
global T t dt dx J x u u2 un B C a F ul up

B = 20;
J = 200;
T = 8;
C = 0.8;

ul=1;
up=2;
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
  plot(x,un,'-','Marker','x');
  drawnow;
  %pause;
  %  axis([-B B 0 1.4]);
  if nstep==1
    pause
  end
  nstep=nstep+1;
end
lfstep;
plot(x,un,'-','Marker','x');
drawnow;
  
function ii=init
% inicializace, pocatecni podminky
global T dt dx J x u u2 un B C a F ul up 

dx = 2*B/J;
x = zeros(J,1);
for j=1:J
  x(j) = -B + j*dx;
end
F = zeros(J,1);
u = zeros(J,1);
u2 = zeros(J,1);
un = zeros(J,1);
for j=1:J
    if x(j) > 0
      u(j) = up;
    else
      u(j) = ul;
    end
%  if abs(x(j)) < 2*pi
%      u(j) = (1 + cos(x(j)/2))/2;
%  end
end
dt = C*dx/abs(max(u));

function lf=lfstep
% jeden casovy krok Lax-Friendrichsovo schema
global t T dt dx J x u u2 un B C a F
t = t + dt;
% prvni pulkrok
flux(u,J);
for j=1:J-1
  u2(j) = (u(j) + u(j+1))/2 - dt*(F(j+1) - F(j))/(2*dx);
end
% druhy pulkrok
flux(u2,J-1);
for j=2:J-1
  un(j) = (u2(j-1) + u2(j))/2 - dt*(F(j) - F(j-1))/(2*dx);
end

% okrajove podminky
un(1) = u(2);
un(J) = u(J-1);
u = un;

function lw=lwstep
% jeden casovy krok Lax-Wendroff schema
global t T dt dx J x u u2 un B C a F
t = t + dt;
% prvni pulkrok
flux(u,J);
for j=1:J-1
  u2(j) = (u(j) + u(j+1))/2 - dt*(F(j+1) - F(j))/(2*dx);
end
% druhy pulkrok
flux(u2,J-1);
for j=2:J-1
  un(j) = u(j) - dt*(F(j) - F(j-1))/dx;
end

% okrajove podminky
un(1) = u(2);
un(J) = u(J-1);
u = un;

function ff = flux(uu,JJ)
global F
for j=1:JJ
  F(j) = uu(j)^2/2;
end


function lf=lfstep1
% jeden casovy krok Lax-Friendrichsovo schema
global t T dt dx J x u u2 un B C a F
t = t + dt;
% jedenkrok
flux(u,J);
for j=2:J-1
  un(j) = (u(j-1) + u(j+1))/2 - dt*(F(j+1) - F(j-1))/(2*dx);
end
% okrajove podminky
un(1) = u(2);
un(J) = u(J-1);
u = un;


