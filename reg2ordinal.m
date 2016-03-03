function data=reg2ordinal(data);
[m,n]=size(data);

minv=min(data(:,n));
maxv=max(data(:,n));

interval=(maxv-minv)/5;

data=sortrows(data,n);

for i =1:m;
    dival=round((data(i,n)-minv)/interval);
    data(i,n)=dival+1;
end



