function  [ratotest,rato,B,matD,matE,matF,matG,matH]  =ADMMmat(edata);
%test set
%x=[1,2;3,4];
%z=[5,6;7,8];
%产生的w中，main effects ???interactions 的顺序是混合???%w=generateW(x,z)
%可以用matVer这个代码来生成数据，暂时不用genereateW这个函数
lambda.fir=0.01;
lambda.sec=0.5;
lambda.thi=0.99;
lambda.four=0.2;%random
%初始化的B是全0

%synthesis data 
%[B,wtr,wte,Ytr,Yte,dataname]=matVer();

%load real world data
[B,wtr,wte,Ytr,Yte,dataname]=loadData(edata);
matD=B;
matE=B;
matF=B;
matG=B;
matH=B;
%g.gamma 拉格朗日乘子,这里都是.与B大小相同的矩阵形???

g.gamma1=B;
g.gamma2=B;
g.gamma3=B;
g.gamma4=B;
g.gamma5=B;
    
rho=1;  



iter=100;
e.abs=1e-4;
e.rel=1e-4;
[B,matD,matE,matF,matG,matH]=estimate(wtr,Ytr,lambda,rho,B,matD,matE,matF,matG,matH,g,iter,e);
rato=mypredict(wtr,B,Ytr);
ratotest=mypredict(wte,B,Yte);
saveResult(rato,ratotest,dataname);


function saveResult(rato,ratotest,dataname)
    fs=fopen('result.txt','at');
    fprintf(fs,'%s\n',dataname)
    fprintf(fs,'%f\n',rato);
    fprintf(fs,'%f\n',ratotest);
    fclose(fs);



%参数w是wtr的这种矩阵形???
function [B,matD,matE,matF,matG,matH]=estimate(w,y,lambda,rho,B,matD,matE,matF,matG,matH,g,iter,e);

for i=1:iter
    %uupdate each variable
    newB=update_b(w,y,B,matD,matE,matF,matG,matH,rho,g);
    newmatD=update_D(newB,w,rho,g,lambda);
    newmatE=update_E(newB,w,rho,g,lambda);
    newmatF=update_F(newB,w,rho,g,lambda);
    newmatG=update_G(y,newB,g,rho,lambda);
    newmatH=update_H(w,y,newB,g,rho,lambda);
    newg=updategammas(w,newB,newmatD,newmatE,newmatF,newmatG,newmatH,g,rho);
    

    %dual resudual s
    d.diff=sum(sum(rho^2*(newmatD-matD).^2));
    e.diff=sum(sum(rho^2*(newmatE-matE).^2));
    f.diff=sum(sum(rho^2*(newmatF-matF).^2));
    g.diff=sum(sum(rho^2*(newmatG-matG).^2));
    h.diff=sum(sum(rho^2*(newmatH-matH).^2));
    s.norm=sqrt(d.diff+e.diff+f.diff+g.diff+h.diff);

    %primal residual r
    d.dif2=sum(sum((newB-newmatD).^2));
    e.dif2=sum(sum((newB-newmatE).^2));
    f.dif2=sum(sum((newB-newmatF).^2));
    g.dif2=sum(sum((newB-newmatG).^2));
    h.dif2=sum(sum((newB-newmatH).^2));
    r.norm=sqrt(d.dif2+e.dif2+f.dif2+g.dif2+h.dif2);
    
    allmat=[newmatD,newmatE,newmatF,newmatG,newmatH];
    
    [brow,bcol]=size(newB);
    crit1=max(sqrt(sum(sum(newB.^2))),sqrt(sum(sum(allmat.^2))))*e.rel+sqrt(brow*bcol)*e.abs;
    vg=[];
    vg=[vg,newg.gamma1(:)];
    vg=[vg,newg.gamma2(:)];
    vg=[vg,newg.gamma3(:)];
    vg=[vg,newg.gamma4(:)];
    vg=[vg,newg.gamma5(:)];
    crit2=norm(vg)*e.rel+sqrt(brow*bcol)*e.abs;

    if(r.norm<crit1 & s.norm<crit2)
        B=newB;
        matD=newmatD;
        matE=newmatE;
        matF=newmatF;
        matG=newmatG;
        matH=newmatH;
        g=newg;
        return 
    else
        B=newB;
        matD=newmatD;
        matE=newmatE;
        matF=newmatF;
        matG=newmatG;
        matH=newmatH;
        g=newg;
        %update rho
        if(r.norm>10*s.norm)
            rho=2*rho;
        else if(r.norm*10<s.norm)
            rho=rho/2;
              end
         end
    end
end
%after iter k times ,just jump out the iteration

%return





function  update_b2(w,y,B,matD,matE,matF,matG,matH,rho,g);

[u,s,v]=svd(w);
v=v';
grank=sum(sum(s~=0,2));
s=s(1:grank,1:grank);
u=u(:,1:grank);
v=v(1:grank,:);





% g.list 可以用结构体
%update d e g 与原代码是相同的???
function X=update_b(w,y,B,matD,matE,matF,matG,matH,rho,g);
%w is a n*(p1+1*p2+1) matrix ,n is the number of examples
%convert D,E,F,G,H to p1+1 * p2+1 vector 
[n,m]=size(w);
[brow,bcol]=size(matD);
vmatD=matD(:);
vmatE=matE(:);
vmatF=matF(:);
vmatG=matG(:);
vmatH=matH(:);
vg.gamma1=g.gamma1(:);
vg.gamma2=g.gamma2(:);
vg.gamma3=g.gamma3(:);
vg.gamma4=g.gamma4(:);
vg.gamma5=g.gamma5(:);
%solve B according to ||y-wx||^2 this formula x=(w^t*w)^{-1}w^t*y,so next construct y and w
Y=y/sqrt(n);
Y_down=(rho*(vmatD+vmatE+vmatF+vmatG+vmatH)-(vg.gamma1+vg.gamma2+vg.gamma3+vg.gamma4+vg.gamma5))/sqrt(5*rho);
Y=[Y;Y_down];

W=w/sqrt(n);
W_down=sqrt(5*n)*ones(m,m);
W=[W;W_down];
%this inverse in time consuming
X=pinv(W'*W)*W'*Y
%X=gradientDescentB(Y,W);

%reshape in column order
X=reshape(X,brow,bcol);




%g.gama1,g.gama2,etc.
function Dnew=update_D(B,w,rho,g,lambda);
D.new=B;
brow=size(B,1);
D.new(1,:)=B(1,:)+g.gamma1(1,:)/rho;
newmat=B(2:brow,:)+g.gamma1(2:brow,:)/rho;
%求newmat 每一行的二范???
normmat=[];
for i =1:size(newmat,1)
    normmat=[normmat;norm(newmat(i,:))];
end
%first para is a vector,so dot divide is necessary , and num subtract
%vector is ok
coef=pmax(1-(lambda.fir/rho)./normmat,0,true);
%newmat 的每?????，乘上每?????的coef
newmat2=bsxfun(@times,newmat,coef);
D.new(2:brow,:)=newmat2;
Dnew=D.new;


function Enew=update_E(B,w,rho,g,lambda);
E.new=B;
bcol=size(B,2);
E.new(:,1)=B(:,1)+g.gamma2(:,1)/rho;
newmat=B(:,2:bcol)+g.gamma2(:,2:bcol)/rho;

normat=[];
for i =1:size(newmat,2)
    normat=[normat,norm(newmat(:,i))];
end
%coef 行
coef=pmax(1-(lambda.sec/rho)./normat,0,false);
newmat2=bsxfun(@times,newmat,coef);
E.new(:,2:bcol)=newmat2;
Enew=E.new;

function Fnew=update_F(B,w,rho,g,lambda);
F.new=B+g.gamma3/rho;
% now we take the part F_{-0,-0}
[frow,fcol]=size(F.new);
tempf=F.new(2:frow,2:fcol);
%matrix - number is ok
tempf2=sign(tempf).*fmax(abs(tempf)-lambda.thi/rho,0);
F.new(2:frow,2:fcol)=tempf2;
Fnew=F.new;

%solving nuclear norm regularized problem using soft-shrinkage
function Gnew=update_G(y,B,g,rho,lambda);
G.new=B+g.gamma4/rho;
%now update data the part G_{-0,-0}
[grow,gcol]=size(G.new);
tempg=G.new(2:grow,2:gcol);
%apply svd to tempg,
%tempg=u*s*v'
[u,s,v]=svd(tempg);
v=v';
grank=sum(sum(s~=0,2));
s=s(1:grank,1:grank);
u=u(:,1:grank);
v=v(1:grank,:);

sdiag=diag(s);%vector
sdiag=fmax(sdiag,lambda.four/rho);
sdiag=diag(sdiag);%matrix

G_0=u*sdiag*v;
G.new(2:grow,2:gcol)=G_0;
Gnew=G.new;


%using gradient descent method
function Hnew=update_H(w,y,B,g,rho,lambda);
[meanfeat,uniqueY]=separateClasses(w,y);
[numcls,numf]=size(meanfeat);
diffmean=zeros(numcls-1,numf);
k=1;
for i =2:numcls
    diffmean(k,:)=meanfeat(i,:)-meanfeat(i-1,:);
    k=k+1;
end
Hnew=gradientDescentH(diffmean,y,B,g,rho,lambda);



%help fun for update H using gradient descent
function  H=gradientDescentH(diffmean,y,B,g,rho,lambda)
%convert matrix into vector for convenience
[orirow,oricol]=size(B);
vB=B(:);
vg.gamma5=g.gamma5(:);
[hrow,hcol]=size(vB);
H_old=zeros(hrow,1);
epsilon=0.0003;
gamma=0.0001;
iterk=0;
H_best=H_old;
%may can't jump out this loop,better set a iteration time
iterset=1000;
for it=1:iterset
    firderi=rho*(H_old-(vB+vg.gamma5/rho));
    fullderi=firderi+gradientHingeLoss(diffmean,H_old);
    H_new=H_old-gamma*fullderi;
    if(norm(H_new-H_old)<epsilon)
        H_best=H_new;
        break;
    end
    H_old=H_new;
    H_best=H_old;%if don't have this line then after iteration H_best would not have been assigned
    iterk=iterk+1;
end
H=H_best;
%convert to matrix again
H=reshape(H,orirow,oricol);


%calculate the gradient of the hinge loss,H is a  column vector
function deri=gradientHingeLoss(diffmean,H)
len=size(diffmean,1);
lenH=size(H,1);
deri=zeros(lenH,1);

for i =1:len
    z=diffmean(i,:)*H;
    if(z<=0)
        deri=deri+(-1*diffmean(i,:)');   
    else if(z>0 & z<1 )
        deri=deri+(diffmean(i,:)'*(diffmean(i,:)*H-1));
     end
 end
end

%update g.gamma
function gnew=updategammas(w,B,matD,matE,matF,matG,matH,g,rho);
 g.new=g;
 g.new.gamma1=g.gamma1+rho*(B-matD);
 g.new.gamma2=g.gamma2+rho*(B-matE);
 g.new.gamma3=g.gamma3+rho*(B-matF);
 g.new.gamma4=g.gamma4+rho*(B-matG);
 g.new.gamma5=g.gamma5+rho*(B-matH);
 gnew=g.new;


%help function to generate W 
function xz=outdot(x,z)
[xrow,xcol]=size(x);
[zrow,zcol]=size(z);
xz=[];
mn=xcol*zcol;
for i =1 : xrow
    temp=x(i,:)'*z(i,:);
    xz=[xz;reshape(temp,1,mn)];
end


function w=generateW(x,z)
xrow=size(x,1);
zrow=size(z,1);
x=[ones(xrow,1),x];
z=[ones(zrow,1),z];
w=outdot(x,z);

%help fun for update d,x is a vector and y is a scalar
function coef=pmax(x,y,isD)
if (isD==true)
xrow=size(x,1);
y=repmat(y,xrow,1);
resu=x>y;
ret=x;
ret(resu,:)=x(resu,:);
ret(~resu,:)=y(~resu,:);
coef=ret;

else 
  xcol=size(x,2);
  y=repmat(y,1,xcol);
  resu=x>y;
  ret=x;
  ret(:,resu)=x(:,resu);
  ret(:,~resu)=y(:,~resu);
  coef=ret;
end

%help fun for update f and g
function ret= fmax(x,y);
%x is a matrix ,y is a scalar
resu=x<y;
x(resu)=0;
ret=x;

%help function for update H,about the hinge loss , separate data for different class
function [meanfeat,uniqueY] = separateClasses(w,y);
%y is a vector,sort in ascending order,Y=y(I)
[Y,I]=sort(y);
W=w(I,:);
[numEx,numfeat]=size(W);
uniqueY=unique(Y);
numrealCls=length(uniqueY);
%calculate the mean features of different classes
%numclass is the number of different classes,cla is the class in ascending order 
%note that numclasses can be zero due to the continuous number in cla
%[numclass,cla]=hist(Y);      hist may not work here 
%lencla=length(numclass);

%这个代码写的很好！！,不过没有考虑样本数为0的类
d=diff([Y;max(Y)+1]);
count=diff(find([1;d]));
numclass=count;
cla=Y(find(d));

lencla=length(numclass);
meanfeat=zeros(numrealCls,numfeat);
startpos=1;     
realClspos=1;
for i =1:lencla;
    if(numclass(i,:)~=0);
       meanfeat(realClspos,:)=mean(W(startpos:startpos+numclass(i,:)-1,:));
       realClspos=realClspos+1;
       startpos=startpos+numclass(i,:);
    end
end


%gradient descent method for update B,in case the inverse opeation in b is very time-consuming
function B=gradientDescentB(Y,W);
epsilon = 0.0003;
gamma=0.0001;
x_old=zeros(size(W,2),1);
iterk=0;
siter=100;
for it=1:100
    x_new=x_old-gamma*(W'*W*x_old-W'*Y);
    if(norm(x_new-x_old)<epsilon)
        x_best=x_new;
        break;
     end
     x_old=x_new;
     x_best=x_old;
     iterk=iterk+1;
 end
 B=x_best;

            
