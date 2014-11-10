function results=function_sum_results(M,results,fourier_1D_u,fourier_1D_v)
N=M;

for m=-M:1:M    
    for n=-N:1:N
        results(m+M+1,n+N+1)=results(m+M+1,n+N+1)+fourier_1D_u(m+M+1)*fourier_1D_v(n+N+1);
    end
end