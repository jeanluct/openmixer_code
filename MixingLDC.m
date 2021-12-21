%LDC Mixer
%parameter study U2



%% Markov chain

% tree
n = 10; x = linspace(-1,1,n+2)'; % n^2 sample points for each box (uniform inner grid)
x=x(2:(length(x)-1));
[XX,YY] = meshgrid(x,x);
X = [ XX(:) YY(:) ];     
c = [3 0]; r = [4 1]; 
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
tc=tree.count(-1);
b=tree.boxes(-1);

c = b(1:2,:); r = b(3:4,1);  % center and radii of the boxes
n = size(c,2); E = ones(n,1);    % sample points in all boxes
BB = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1));


%% source
mu=zeros(1,tc);
Idex=c(1,:)<0 & c(2,:)>0;
mu(Idex)=1;
Idex2=c(1,:)<0 & c(2,:)<0;
mu(Idex2)=-1;


%parameter
%a=6;
%b=1;
%tf=1;

U1=ones(1,25)*9;
U2=6:0.25:12;

h=0.01;
uwx=1;
uwy=0;


%%

par=U1;
%
Eigenvalues=zeros(10,length(par));
EigenvecR1=zeros(tc,length(par));
EigenvecL1=zeros(tc,length(par));
EigenvecR2=zeros(tc,length(par));
EigenvecL2=zeros(tc,length(par));
ETM=zeros(tc,length(par));

Vinv=zeros(length(par),tc);
%%

for i = 1:length(par)

f= @(x)myrk4_end5(@liddrivencavityMixer,0,1, h,x, uwx,uwy,1,U1(i),U2(i),6,1);
P=tpmatrix(tree,f, X, -1,1);
P=P';

% % invariante Massenverteilung
Vinv(i,:)=transpose((speye(tc,tc)-P))\mu';
% 
% % Spektrum
[R1, S1] = eigs(P,50,'LR'); 
S1=diag(S1);
[s1, is]=sort(S1, 'descend');
R1=R1(:,is);

Eigenvalues(:,i)=s1(1:10);
EigenvecR1(:,i)=R1(:,1);
EigenvecR2(:,i)=R1(:,2);

[L1, S1] = eigs(P',50,'LR');   
S1=diag(S1);
[s1, is]=sort(S1, 'descend');
L1=L1(:,is);
EigenvecL1(:,i)=L1(:,1);
EigenvecL2(:,i)=L1(:,2);

% erwartete Zeiten im Mixer bei Start in Boxen 
etM=(speye(tc,tc)-P)\ones(tc,1);
ETM(:,i)=etM;

end

%%
save('lidmixerQuadraticBoxes', 'Eigenvalues','EigenvecR1', 'EigenvecR2', 'EigenvecL1', 'EigenvecL2', 'ETM', 'Vinv')


%%
for i=1:25
figure; %subplot(5,5,i)
boxplot2(tree,'depth',tree.depth,'density',Vinv(i,:));
axis tight; axis equal;
shading flat; caxis([-1,1]);
axis([6,7,-1,1])
set(gca,'XTick',[])
set(gca,'YTick',[])
%t=sprintf('U_1=%g, U_2=%g', U1(i),U2(i));
%title(t)
end








