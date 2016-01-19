%help fun for update d and e,x is a vector and y is a scalar
function coef=pmax(x,y)
%可以用bsxfun
coef =bsxfun(@mycompare,x,y);

%a,b 是大小相同的向量，而不是标量
function ret=mycompare(a,b)
resu=a>b;
resu
ret=a;
a
ret(resu,:)=a(resu,:);
ret
b
ret(~resu,:)=b(~resu,:);

