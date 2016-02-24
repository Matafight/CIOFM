function [B,wtr,wte,Ytr,Yte]=loadData();

data=load('diabetes.data.ord');
[datam,datan]=size(data);

newindex=randperm(datam);
data=data(newindex,:);

trlen=datam*0.7;
telen=datam-trlen;

Xtr=data(1:trlen,1:datan-1);
Ytr=data(1:trlen,datan);
Xte=data(trlen+1:datam,1:datan-1);
Yte=data(trlen+1:datam,datan);

%generate full matrix of covariates
wtr=[];
wte=[];
[m,n]=size(Xtr);
X1=[ones(m,1),Xtr];

[tem,ten]=size(Xte);
X2=[ones(tem,1),Xte];
B=zeros(n+1,n+1);   

for  i=1:n+1
    for j =1:n+1
   wtr=[wtr,X1(:,j).*X1(:,i)];
   wte=[wte,X2(:,j).*X2(:,i)];
    end
end

