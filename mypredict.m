function rato=mypredict(w,B,y);

%B is a matrix .
[m,n]=size(w);
preval=w*reshape(B,n,1);
diffval=round(preval-y);
 
%return accuracy
ind=find(diffval == 0);
rato=length(ind)/m;



