function  w=ADMMmat();
%test set
x=[1,2;3,4];
z=[5,6;7,8];
%产生的w中，main effects 和 interactions 的顺序是混合的
w=generateW(x,z)



% g.list 可以用结构体
%update d e g 与原代码是相同的。
function update_b();

function update_D(B,w,rho=0,g.list);




function update_E();

function update_F();




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





