% DG Mixer, parameter study eps, A

%clear all;

%% boxes and test points

n = 10; x = linspace(-1,1,n+2)'; % n^2 sample points for each box (uniform inner grid)
x=x(2:(length(x)-1));
[XX,YY] = meshgrid(x,x);
X = [ XX(:) YY(:) ];     
c = [1 0.5]; r = [2 0.5]; 
tree = Tree(c, r); 

% the box collection
sd=tree.sd;
sd(1:4)=[0,1,0,0];
tree.sd=sd;

steps=16;
for i = 1:steps
    tree.set_flags('all', 8);           
    tree.subdivide;                      
end

b=tree.boxes(-1);
c = b(1:2,:); r = b(3:4,1);  
w=c(1,:)>-0.5 & c(1,:)<2.5;
w=reshape(char(w+'0'),2^steps,1);
tree.set_flags(w,9,-1);
tree.remove(9);              % delete boxes
tc=tree.count(-1);
b=tree.boxes(-1);
c = b(1:2,:); r = b(3:4,1);  % center and radii of the boxes



%% DG mixer parameters
A=0.5;
uwx=0.5;
uwy=0;
h=0.01;
omega=2*pi;

%eps=0.25
%A=0.1:0.1:12
Eps=0:0.025:2.5;


%% Source-distribution
mu=zeros(tc,1);
Idex=c(1,:)<0 & c(2,:)>0.5;
mu(Idex)=1;
Idex2=c(1,:)<0 & c(2,:)<0.5;
mu(Idex2)=-1;

%%
par=Eps;  
%par=A;

Eigenvalues=zeros(10,length(par));
EigenvecR1=zeros(tc,length(par));
EigenvecL1=zeros(tc,length(par));
EigenvecR2=zeros(tc,length(par));
EigenvecL2=zeros(tc,length(par));
ETM=zeros(tc,length(par));
Vinv=zeros(length(par),tc);


for i = 1:length(par)

f= @(x)myrk4_end(@gyresMixer,0,1, h,x,A,par(i),omega, uwx, uwy);
%f= @(x)myrk4_end(@gyresMixer,0,1, h,x,A(i),eps,omega, uwx, uwy);

% transition matrix
P=tpmatrix(tree,f, X, -1,1);
P=P';  

% % invariante Massenverteilung v_inv
Vinv(i,:)=transpose((speye(tc,tc)-P))\mu;
% 
% % Spektrum
[R1, S1] = eigs(P,10,'LR'); 
Eigenvalues(:,i)=diag(S1);
EigenvecR1(:,i)=R1(:,1);                           

EigenvecR2(:,i)=R1(:,2);
[L1, ~] = eigs(P',10,'LR'); 
EigenvecL1(:,i)=L1(:,1);
EigenvecL2(:,i)=L1(:,2);

% erwartete Zeiten im Mixer bei Start in Boxen 
etM=(speye(tc,tc)-P)\ones(tc,1);
ETM(:,i)=etM;

end

%P=0
%save('NeuMixingDGATP100')
%save('NeuMixingDGepsTP100')

%% v_inv
%pI=[1,4,9,15,21,23,28,32,36,41,43,54,56,60,64,72,81]+4;
%pI=1:4:101
Eps(pI)

ki=1;

for i=pI
figure;
boxplot2(tree,'depth',tree.depth,'density',Vinv(i,:));
axis tight; axis equal; axis([2,2.5,0,1])
shading flat; caxis([-1,1]);
axis off;
ki=ki+1;
end


