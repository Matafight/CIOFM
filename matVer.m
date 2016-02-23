function  [B,wtr,wte,Y_high_SHtr,Y_high_SHte]=matVer();
%this is the function that generate both training and testing data
%ç”Ÿæˆæ•°æ®,X and Z
Xtr=randn(100,10);
Ztr=randn(100,15);

Xte=randn(100,10);
Zte=randn(100,15);

%éœ?¦æ ‡å‡†åŒ?

%generate full matrix of covariates
wtr=[];
wte=[];
[m,n]=size(Xtr);
X1=[ones(m,1),Xtr];
Z1=[ones(m,1),Ztr];

X2=[ones(m,1),Xte];
Z2=[ones(m,1),Zte];

for  i=1:16
    for j =1:11
   wtr=[wtr,X1(:,j).*Z1(:,i)];
   wte=[wte,X2(:,j).*Z2(:,i)];
    end
end

%generate response variables with signal from
%first 5 x features and 5 z features
B=zeros(16,11);
B_high_SH=B;
B_high_SH(1:6,1:6)=1;

[m,n]=size(B_high_SH);
Y_high_SHtr =wtr*reshape(B_high_SH,m*n,1)+randn(100,1);
Y_high_SHte =wte*reshape(B_high_SH,m*n,1)+randn(100,1);



