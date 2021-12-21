

function [h,vario,nh]=vario2(X,Z,n,delta,periodx,periody)

minx = min(X,[],1);
maxx = max(X,[],1);
maxd = sqrt(sum((maxx-minx).^2))/2;
%delta = maxd/((n-1)*2);

h=linspace(0,maxd,n);
vario=zeros(1,length(h));
nh=zeros(1,length(h));



%distance matrix
%D1=(X(:,1) - X(:,1).');
%D2=(X(:,2) - X(:,2).');

dm=sqrt( (X(:,1) - X(:,1).').^2 + (X(:,2) - X(:,2).').^2 ) ;
%dm=sqrt(D1.^2+D2.^2);


if periodx>0
%D1=(X(:,1)+periodx - X(:,1).');
 
dm2=sqrt( (X(:,1)+periodx - X(:,1).').^2 + (X(:,2) - X(:,2).').^2 ) ;
%dm2=sqrt(D1.^2+D2.^2);

dm=min(dm,dm2);
end

if periody>0
 %D1=(X(:,1) - X(:,1).');
 %D2=(X(:,2)+periody - X(:,2).');  
    
dm2=sqrt( (X(:,1) - X(:,1).').^2 + (X(:,2)+periody - X(:,2).').^2 ) ;
dm3=sqrt( (X(:,1)+periodx - X(:,1).').^2 + (X(:,2)+periody - X(:,2).').^2 ) ;

%dm2=sqrt(D1.*D1+D2.*D2);
%D1=(X(:,1)+periodx - X(:,1).');
%dm3=sqrt(D1.^2+D2.^2);


dm=min(dm,dm2);
dm=min(dm,dm3);
end


%variogram


D=(Z-Z.');
%Zm=(Z-Z.').^2;
Zm=D.*D;


for i=1:length(h)

Ir=(h(i)-delta<=dm & dm<=h(i)+delta);
vario(i)=sum(sum(Ir.*Zm))/(2*nnz(Ir));

nh(i)=nnz(Ir)/2;
end