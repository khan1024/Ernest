function y=fasttoeplitz(x,index)


if nargin>1;
    y=x(index);
else
    x=x(:);
    n=floor((length(x)+1)/2);
    y=zeros(n,n);
    for ii=1:n;
        y(:,ii)=x(n-ii+1:2*n-ii);
    end;
    %y=retsparse(y);%  
end;