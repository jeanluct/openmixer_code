%% Mixing measures

% load v_inv, b and c (box-discretization), and par, uwx

%% Mixnorm

[X,Y]=meshgrid(-0.5:0.001:0, 0:0.001:1);  %% unmixed region of the Double-Gyre Mixer
%[X,Y]=meshgrid(-1:0.001:0, -1:0.001:1);  %% unmixed region of the Lid Driven Cavity Mixer

[Nx, Ny]=size(X);
b=tree.boxes(-1);
v=Vinv(1,c(1,:)<0);

[ndX,ndY]=ndgrid(unique(b(1, c(1,:)<0)'),unique(b(2, c(1,:)<0)'));
V=zeros(size(ndX));

for li=1:size(ndX,1)
    for lj=1:size(ndY,2)
V(li,lj)=v(tree.search([ndX(li,lj),ndY(li,lj)]',-1));
    end
end

F = griddedInterpolant(ndX,ndY,V,'nearest');
vint=F(X',Y'); %%% interpoliert v bzgl. des Gitters X,Y

figure; subplot(7,2,1)
contourf(X,Y, vint'); 
mixnorm0=computeMixNorm(vint', Nx, Ny);
axis tight; axis equal;
t=sprintf('rm=%g',mixnorm0/mixnorm0);
title(t);


RM=zeros(8,1);

[X,Y]=meshgrid(2:0.001:2.5, 0:0.001:1);  %% domain of the mixed region Double Gyre
%[X,Y]=meshgrid(6:0.001:7, -1:0.001:1);  %% domain of the mixed region Lid
%Driven Cavity
[Nx, Ny]=size(X);

for i=1:length(par)
v=Vinv(i,:);
[ndX,ndY]=ndgrid(unique(b(1, :)'),unique(b(2, :)'));
V=zeros(size(ndX));
for li=1:size(ndX,1)
    for lj=1:size(ndY,2)
V(li,lj)=v(tree.search([ndX(li,lj),ndY(li,lj)]',-1));
    end
end

F = griddedInterpolant(ndX,ndY,V,'nearest');
vint=F(X',Y'); %%% interpoliert v bzgl. des Gitters X,Y

%figure(53); subplot(7,2,i+1)
%contourf(X,Y, vint'); 
mixnorm=computeMixNorm(vint', Nx, Ny);
RM(i)=mixnorm;
%axis tight; axis equal;
end


%% Scale of segregation

%% box-variogram: kategoriale Daten
VinvK=double(Vinv(:,c(1,:)>2)<0);

% variogram
np=55;                                  %number of bins
minx = min(b(1:2,b(1,:)>2)',[],1);
maxx = max(b(1:2,b(1,:)>2)',[],1);
maxd = sqrt(sum((maxx-minx).^2))/2;
delta = maxd/((np-1)*2); %bin width
periodx=uwx;
periody=0;

% figure; hold on;
slK=0;
k=1;
for i=1:length(par)
[hv,vario,nh] = vario2(b(1:2,b(1,:)>2)',VinvK(i,:)',np,delta,periodx,periody);
%plot(hv,vario,'.-');
% magnitude of the slope near the origin 
slK(k)=(vario(2)-vario(1))/hv(2);
k=k+1;
end
% title('Variogram');
% axis([0 0.6 0 1.3]);


%% Plot mixing measures

figure; hold on;
plot(par, RM/mixnorm0,'.-')
plot(par,var(Vinv(:,c(1,:)>2),1,2) ,'.-')
plot(par,0.5*1./slK, '.-');


