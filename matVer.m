function  matVer();
%生成数据,X and Z
Xtr=randn(100,10);
Ztr=randn(100,15);

Xte=randn(100,10);
Zte=randn(100,15);

%需要标准化?

%generate full matrix of covariates
wtr=[];
wte=[];
[m,n]=size(Xtr);
X1=[ones(m,1),Xtr];
Z1=[ones(m,1),Xtr];

X2=[ones(m,1),Xte];
Y2=[ones(m,1),Yte];

for  i=1:16 :
    for j =1:11:
   wtr=[wtr,X1(:,i).*Z1(:,j)]
   wte=[wte,X2(:,i).*Z1(:,j)]
    end
end

%generate response variables with signal from
%first 5 x features and 5 z features
B=zeros(16,11);
B_high_SH=B;
B_high_SH(1:6,1:6)=1;

[m,n]=size(B_high_SH);
Y_high_SHte =wtr*reshape(B_high_SH,m*n,1)+randn(100,1);
Y_high_SHte=wte*reshpae(B_high_SH,m*n,1)+randn(100,1);



