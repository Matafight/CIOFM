
function rato=mypredict(w,B,y);

%B is a matrix .
[m,n]=size(w);
preval=w*reshape(B,n,1);
diffval=round(preval-y);
ind=find(diffval);
rato=length(ind)/m;



