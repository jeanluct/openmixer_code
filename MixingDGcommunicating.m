
% Example: (Self-)communicating classes DG mixer

%% test points and boxes

n = 10; x = linspace(-1,1,n+2)'; % n^2 sample points for each box (uniform inner grid)
x=x(2:(length(x)-1));
[XX,YY] = meshgrid(x,x);
X = [ XX(:) YY(:) ];     
c = [1 0.5]; r = [2 0.5]; 
tree = Tree(c, r); 

% 
sd=tree.sd;
sd(1:4)=[0,1,0,0];
tree.sd=sd;

steps=16;
for i = 1:steps
    tree.set_flags('all', 8);           
    tree.subdivide;                      
end

b=tree.boxes(-1);
c = b(1:2,:); %r = b(3:4,1);  % center and radii of the boxes


w=c(1,:)>-0.5 & c(1,:)<2.5;
w=reshape(char(w+'0'),2^steps,1);


 tree.set_flags(w,9,-1);
 tree.remove(9);
tc=tree.count(-1);
b=tree.boxes(-1);
c = b(1:2,:); r = b(3:4,1);  % center and radii of the boxes



%% DG mixer parameters
A=0.5;
uwx=0.5;
uwy=0;
h=0.01;
omega=2*pi;

%%
Jmixer=tree.search(b(1:2, c(1,:)<2 & c(1,:)>0));

                           
%%

f= @(x)myrk4_end(@gyresMixer,0,1,h,x,A,0.4,omega,uwx,uwy);
P=tpmatrix(tree,f, X, -1,1);
P=P';

X2=P(Jmixer, Jmixer);

Z=double(X2>0);
 

Rr=(speye(32768)+Z);
%clearvars -except Rr 

%% Reachability matrix
 tic
 R=Rr^(32767);
 toc

%save('MixingDGcommunicatingbeispiel','R')
 %%
 
 Jn=32768;
 
 C = R & R';
 U = unique(C,'rows');
%  figure;
%  spy(U);
nCC=size(U,1); %number of communicatiog classes
AnzC=sum(U,2); %number of boxes in classes

Cv=zeros(1,Jn); %class id
for j=1:Jn
Cv(j)=find(U(:,j)>0);
end

CV=Cv;
iC=find(AnzC>1);
nCC2=length(iC);


%%
[n,bin] = hist(CV,unique(CV));
[~,idx] = sort(-n);
nB=n(idx); % count instances
NB=nB(1:2);
idB=bin(idx); % corresponding values
IDB=idB(1:2);

inC1= CV==IDB(1)';
inC2= CV==IDB(2)';

save('MixingDGcommunicatingbeispiel','R', 'inC1' ,'inC2' ,'CV', 'nCC','nCC2', 'NB') 
 
 
 
%%

BB=b(1:4, c(1,:)<2& c(1,:)>0);  %Boxen im Mixer

%plot of the communication classes

  figure;    
 boxplot2(BB','density',CV);axis tight;shading flat;
axis equal; axis tight; 

 
% plot of the 2 self-communicating classes

       
  C1=BB(1:4,inC1);
  C2=BB(1:4,inC2);
  
 figure; hold on
  h=boxplot2(C1', 'color', 'red');
   set(h,'FaceColor',[0.87,0.49,0])
   set(h,'EdgeColor','none')
  h=boxplot2(C2', 'color', 'blue');
   set(h,'EdgeColor','none')
   set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
  axis equal; axis tight; axis([-0.5,2.5,0,1])
box on;  





