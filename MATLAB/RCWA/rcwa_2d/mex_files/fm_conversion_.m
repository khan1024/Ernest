function fM=fm_conversion_(mMAX,aMAX,fMpp,fMmp,fMpm,fMmm)

fM=zeros(aMAX,aMAX);

for a1=1:aMAX
    for a2=1:aMAX
        
        n1 = floor((a1-1)/(2*mMAX+1)) - mMAX;
        m1 = rem(a1-1, 2*mMAX+1) - mMAX;
        
        n2 = floor((a2-1)/(2*mMAX+1)) - mMAX;
        m2 = rem(a2-1, 2*mMAX+1) - mMAX; 
        
        m = (m2-m1); 
        n = (n2-n1);
        
        if (m>=0) && (n>=0)
            fM(a1,a2) = fMpp(m+1,n+1);
        elseif (m<0) && (n>=0)
            fM(a1,a2) = fMmp(-m+1,n+1);
        elseif (m>=0) && (n<0)
            fM(a1,a2) = fMpm(m+1,-n+1);
        else
            fM(a1,a2) = fMmm(-m+1,-n+1);
        end
    end
end