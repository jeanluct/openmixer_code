%Plots Mixer

% evolution of vk, eigenvectors, support of eigenvectors, expected
% residence times, follow particles (without transition matrix),
% streamlines

%% DG
% test points and box discretization from "MixingDGparameterstudy.m"
uwx=0.5;
f= @(x)myrk4_end(@gyresMixer,0,1,0.01,x,0.5,0.4,2*pi,0.5,0);
P=tpmatrix(tree,f, X, -1,1);
P=P';

%% LDC
% test points and box discretization from "MixingLDC.m"
%uwx=1;
% f= @(x)myrk4_end5(@liddrivencavityMixer,0,1,h,x,1,0,1,9,8,6,1);
% P=tpmatrix(tree,f, X, -1,1);
% P=P';

%% invariante Massenverteilung
Vinv=transpose((speye(tc,tc)-P))\mu;
% 
%% Spektrum
[R1, S1] = eigs(P,30,'LR'); 
[L1, ~] = eigs(P',30,'LR'); 

%% Evolution of vk
 
k=8;

vk=sparse(1,tc);
Pok=speye(tc);

for ts=1:k
       
 vk=vk+mu'*Pok;
 Pok=Pok*P;
 
figure;
boxplot2(tree,'depth',tree.depth,'density',vk);
axis tight; axis equal; shading flat; caxis([-1,1]); %colorbar;% colormap('parula(7)')

end   


%% plot eigenvectors, support of eigenvectors, residence times

%% eigenvectors

for i=1:2
figure;

boxplot2(tree,'depth',tree.depth,'density',L1(:,i));
axis tight; axis equal; %axis([2,2.5,0,1])
shading flat; colorbar;

axis off;
end

for i=1:2
figure;

boxplot2(tree,'depth',tree.depth,'density',R1(:,i));
axis tight; axis equal; %axis([2,2.5,0,1])
shading flat; colorbar;

axis off;
end

%% support of the eigenvectors

for i=1:2
id1=(abs(L1(:,i))>1*10^(-12)) | (abs(R1(:,i))>1*10^(-12));
id2=(abs(L1(:,i))>1*10^(-12)) & (abs(R1(:,i))>1*10^(-12));

% id1=(abs(L1(:,i))>1*10^(-3)) | (abs(R1(:,i))>1*10^(-3));
% id2=(abs(L1(:,i))>1*10^(-3)) & (abs(R1(:,i))>1*10^(-3));
 
  C1=b(1:4,id1);
  C2=b(1:4,id2);
  
 figure; hold on
  h=boxplot2(C1', 'color', 'green');
   %set(h,'FaceColor',[0.87,0.49,0])
   set(h,'EdgeColor','none')
  h=boxplot2(C2', 'color', 'blue');
   set(h,'EdgeColor','none')
   set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
  axis equal; axis tight; axis([-0.5,2.5,0,1])
  %axis([-1,7,-1,1])
box on;  
end

%% expected residence times

%etM=(speye(tc,tc)-P')\ones(tc,1);  %forward
etM=(speye(tc,tc)-P')\ones(tc,1);  %backward
figure;
boxplot2(tree,'depth',tree.depth,'density',etM);
axis tight; axis equal;
shading flat; colorbar;
axis off;



%% DG: follow particles (without transition matrix)

%% Partikel, die hineingeschickt werden
c = [1 0.5]; r = [1.5 0.5];               %region with particles that we want to follow
%c = [-0.25 0.5]; r = [0.25 0.25]; 
%c = [-0.13 0.13]; r = [0.01 0.01]; %small 
tree = Tree(c, r); 
% the box collection
steps=16;
for i = 1:steps
    tree.set_flags('all', 8);           
    tree.subdivide;                      
end
b=tree.boxes(-1);
c = b(1:2,:); r = b(3:4,1);  
n = size(c,2); E = ones(n,1);    
BB = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1)); %particles

% parameters of the DG mixer
Eps=0.4;
h=0.01;
A=0.5; %
omega=2*pi;
uwx=0.5;% [0.25,0.5,1];
uwy=0;

T=10;

%
for i=1:1

f= @(x)myrk4_end(@gyresMixer,0,T,h,x,A,Eps,omega, uwx, uwy);
Bild=f(BB);

figure;hold on;

Index1=(BB(:,2)<1 & BB(:,2)>0  );  
plot(BB(Index1,1),BB(Index1,2),'.','MarkerSize',2); %Partikeln in t=0
plot(Bild(Index1,1),Bild(Index1,2),'.','MarkerSize',2); %Position dieser in t=T

Index2=( Bild(:,1)>0 & Bild(:,1)<2.5  ); 
plot(Bild(Index2,1),Bild(Index2,2),'.','MarkerSize',2); %Partikel die noch im Mixer sind in t=T
plot(BB(Index2,1),BB(Index2,2),'.','MarkerSize',2); %Anfangsposition von diesen

axis equal;axis tight;
axis([-0.5,6,0,1])

end


%% Streamlines

%% Streamlines Double Gyre Flow    
T=0;
omega=2*pi;
A=0.5;
eps=0.4;

g = @(t,z) eps*sin(omega*t)*z.^2+(1-2*eps*sin(omega*t))*z;
dg = @(t,z) 2*z.*eps*sin(omega*t)+(1-2*eps*sin(omega*t));


[x,y] = meshgrid(0:0.001:2,0:0.001:1);

dx1 = -pi*A*sin(pi*g(T,x)).*cos(pi*y);%+0.5;
dy1 = pi*A*cos(pi*g(T,x)).*sin(pi*y).*dg(T,x);

figure
startx = [0.1:0.1:0.4  1.9:-0.1:1.6  1 1.9999 1.00001 2 0.0001 0];
starty = [ones(1,length(startx)-6)*0.5  0 1 0 0.0001 1 0.0001];
h=streamline(x,y,dx1,dy1,startx,starty, [0.05, 50000]);
axis equal, axis tight; axis([0,2,0,1]), 

 xticks([0 1 2])
 yticks([0 1])

set(gca,'FontSize', 14);  
set(h, 'Color' , 'blue' , 'LineWidth',1 ) 


%%

T=0.25;
omega=2*pi;
A=0.5;
eps=0.4;

g = @(t,z) eps*sin(omega*t)*z.^2+(1-2*eps*sin(omega*t))*z;
dg = @(t,z) 2*z.*eps*sin(omega*t)+(1-2*eps*sin(omega*t));

[x,y] = meshgrid(0:0.001:2,0:0.001:1);

dx1 = -pi*A*sin(pi*g(T,x)).*cos(pi*y);%+0.5;
dy1 = pi*A*cos(pi*g(T,x)).*sin(pi*y).*dg(T,x);

figure
startx = [0.2:0.2:0.7  1.95:-0.05:1.75 1.35 1.9999 1.351 2 0.0001 0];
starty = [ones(1,length(startx)-6)*0.5 0 1 0 0.0001 1 0.0001];
h=streamline(x,y,dx1,dy1,startx,starty, [0.05, 65000]);
axis equal, axis tight; axis([0,2,0,1])

 xticks([0 1 2])
 yticks([0 1])
set(gca,'FontSize', 14);  
set(h, 'Color' , 'blue', 'LineWidth',1 ) 


%%
T=0.75;
omega=2*pi;
A=0.5;
eps=0.4;

g = @(t,z) eps*sin(omega*t)*z.^2+(1-2*eps*sin(omega*t))*z;
dg = @(t,z) 2*z.*eps*sin(omega*t)+(1-2*eps*sin(omega*t));

[x,y] = meshgrid(0:0.001:2,0:0.001:1);

dx1 = -pi*A*sin(pi*g(T,x)).*cos(pi*y);%+0.5;
dy1 = pi*A*cos(pi*g(T,x)).*sin(pi*y).*dg(T,x);

figure
startx = [startx*(-1)+2 2];
starty = [ones(1,length(startx)-7)*0.5 0 1 0 0.0001 1 0.0001 0.001];
h=streamline(x,y,dx1,dy1,startx,starty, [0.05, 65000]);
axis equal, axis tight; axis([0,2,0,1])

 xticks([0 1 2])
 yticks([0 1])
set(gca,'FontSize', 14);  
set(h, 'Color' , 'blue', 'LineWidth',1 ) 


%% Streamlines Lid-Driven Cavity Flow

T=0;
U1=9;
U2=8;
tf=1;
a=6;
b=1;

[x,y] = meshgrid(0:0.001:6,-1:0.001:1);

f = @(k,y)(2*pi*y/a)*cosh(k*pi*b/a).*sinh(k*pi*y/a)-(2*pi*b/a)*sinh(k*pi*b/a)*cosh(k*pi*y/a);
df = @(k,y) cosh(k*pi*b/a)*((2*pi/a)*sinh(k*pi*y/a)+(2*pi*y/a).*cosh(k*pi*y/a)*(k*pi/a))-(2*pi*b/a)*sinh(k*pi*b/a)*sinh(k*pi*y/a)*(k*pi/a);

C1=(a^2/(2*pi^2*b))/((a/(2*pi*b))*sinh(2*pi*b/a)+1);
C2=(a^2/(4*pi^2*b))/((a/(4*pi*b))*sinh(4*pi*b/a)+1);

dx1 = (((mod(T,tf)<tf/2)*(-1)+(mod(T,tf)>=tf/2)*(1))*(U1*C1*df(1,y).*sin(pi*x/a))+U2*C2*df(2,y).*sin(2*pi*x/a));
dy1 = (((mod(T,tf)<tf/2)*(-1)+(mod(T,tf)>=tf/2)*(1))*(-U1*C1*f(1,y).*cos(pi*x/a)*(pi/a))-U2*C2*f(2,y).*cos(2*pi*x/a)*2*pi/a);


%%

figure
startx = [0.2:0.2: 0.8, 1.8:0.4:4  1 1.8595];
starty = [ones(1,length(startx)-2)*0  1 -0.9999];

h=streamline(x,y,dx1,dy1,startx,starty, [0.05, 215000]);
axis equal, axis tight; axis([0,6,-1,1])

 xticks([0 6])
 yticks([-1 1])
set(gca,'FontSize', 14);  
set(h, 'Color' , 'blue', 'LineWidth',1 ) 

%%

T=0.5;
U1=9;
U2=8;
tf=1;
a=6;
b=1;

[x,y] = meshgrid(0:0.001:6,-1:0.001:1);

f = @(k,y)(2*pi*y/a)*cosh(k*pi*b/a).*sinh(k*pi*y/a)-(2*pi*b/a)*sinh(k*pi*b/a)*cosh(k*pi*y/a);
df = @(k,y) cosh(k*pi*b/a)*((2*pi/a)*sinh(k*pi*y/a)+(2*pi*y/a).*cosh(k*pi*y/a)*(k*pi/a))-(2*pi*b/a)*sinh(k*pi*b/a)*sinh(k*pi*y/a)*(k*pi/a);

C1=(a^2/(2*pi^2*b))/((a/(2*pi*b))*sinh(2*pi*b/a)+1);
C2=(a^2/(4*pi^2*b))/((a/(4*pi*b))*sinh(4*pi*b/a)+1);

dx1 = (((mod(T,tf)<tf/2)*(-1)+(mod(T,tf)>=tf/2)*(1))*(U1*C1*df(1,y).*sin(pi*x/a))+U2*C2*df(2,y).*sin(2*pi*x/a));
dy1 = (((mod(T,tf)<tf/2)*(-1)+(mod(T,tf)>=tf/2)*(1))*(-U1*C1*f(1,y).*cos(pi*x/a)*(pi/a))-U2*C2*f(2,y).*cos(2*pi*x/a)*2*pi/a);

%%
figure
sx = [0.2:0.2: 0.8, 1.8:0.4:4  ];

startx = [sx*(-1)+6 1.5 4.14 ];  
starty = [ones(1,length(startx)-2)*0 1 -0.9999 ]; 
h=streamline(x,y,dx1,dy1,startx,starty, [0.05, 215000]);
axis equal, axis tight; axis([0,6,-1,1])

 xticks([0 6])
 yticks([-1 1])
set(gca,'FontSize', 14);  
set(h, 'Color' , 'blue', 'LineWidth',1 ) 


