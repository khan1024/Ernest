clear all

N_X=5;

m=-N_X:1:N_X;
n=m;

diff_orders=zeros((2*N_X+1)^2,3);
for i=1:1:length(m)
    for j=1:1:length(n)
        r_index=1+(m(j)+N_X)+(2*N_X+1)*(n(i)+N_X);
        diff_orders(r_index,:)=[r_index,m(j),n(i)];
    end
end