function  w=ADMMmat();
%test set
x=[1,2;3,4];
z=[5,6;7,8];
%产生的w中，main effects 和 interactions 的顺序是混合的
w=generateW(x,z)



% g.list 可以用结构体
%update d e g 与原代码是相同的。
function update_b();


%g.gama1,g.gama2,etc.
function D.new=update_D(B,w,rho=0,g,lambda);
D.new=B;
brow=size(B,1);
D.new(1,:)=B(1,:)+g.gamma1(1,:)/rho;
newmat=B(2:brow,:)+g.gamma1(2:brow,:)/rho
%求newmat 每一行的二范数

normmat=[];
for i =1:size(newmat,1)
    normmat=[normmat;norm(newmat(i,:))];
end
%first para is a vector,so dot divide is necessary , and num subtract
%vector is ok
coef=pmax(1-(lambda.fir/rho)./normmat,0);
%newmat 的每一行，乘上每一行的coef
newmat2=bsxfun(@times,newmat,coef);
D.new(2:brow)=newmat2;



function E.new=update_E(B,w,rho=0,g,lambda);
E.new=B;
bcol=size(B,2);
E.new(:,1)=B(:,1)+g.gamma2(1,:)/rho;
newmat=B(:,2:bcol)+g.gamma2(:,2:bcol)/rho;

normat=[];
for i =1:size(newmat,2)
    normat=[normat,norm(newmat(:,i))];
end
%coef 行向量
coef=pmax(1-(lambda.sec/rho)./newmat,0,false);
newmat2=bsxfun(@times,newmat,coef);
E.new(:,2:bcol)=newmat2;


function F.new=update_F(w,B,g,rho=0,lambda);
F.new=B+g.gamma3/rho
% now we take the part F_{-0,-0}
[frow,fcol]=size(F.new);
tempf=F.new(2:frow,2:fcol);
%matrix - number is ok
tempf2=sign(tempf).*fmax(abs(tempf)-lambda.thi/rho,0);
F.new(2:frow,2:fcol)=tempf2;


%solving nuclear norm regularized problem using soft-shrinkage
function G.new=update_G(y,B,g,rho=0,lambda);
G.new=B+g.gamma3/rho;
%now update data the part G_{-0,-0}
[grow,gcol]=size(G.new);
tempg=G.new(2:grow,2:gcol);
%apply svd to tempg,
[u,s,v]=svd(tempg);
v=v';
grank=sum(s~=0,2);
s=s(1:grank,1:grank);
u=u(:,1:grank);
v=v(1:grank,:);

sdiag=diag(s);%vector
sdiag=fmax(sdiag,lambda.four/rho);
sdiag=diag(sdiag);%matrix

G_0=u*sdiag*v;
G.new(2:grow,2:gcol)=G_0;
return G.new


function H.new=update_H();



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
function coef=pmax(x,y,isD=true)
if (isD==true)
xrow=size(x,1);
y=repmat(y,xrow,1);
resu=x>y;
ret=x;
ret(resu,:)=x(resu,:);
ret(~resu,:)=y(~resu,:)
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



