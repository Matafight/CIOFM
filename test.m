function ret=test();
%
a=[1,2;,3,4;5,6];
b=3
ret=fmax(a,b)


%help fun for update f and g
function ret= fmax(x,y);
%x is a matrix ,y is a scalar
resu=x<y;
x(resu)=0;
ret=x;
