function fM = fft2D(FM, maxx)
% fM = Fourierova transformace matice FM

aMAX = maxx(3);
Nf1=size(FM,1);
Nf2=size(FM,2);

fMpp=fft2(FM/(Nf1*Nf2)); fMmp=fft2(FM(Nf1:-1:1,:)/(Nf1*Nf2)); fMpm=fft2(FM(:,Nf2:-1:1)/(Nf1*Nf2)); fMmm=fft2(FM(Nf1:-1:1,Nf2:-1:1)/(Nf1*Nf2));

fM=fm_conversion(maxx(1),aMAX,fMpp,fMmp,fMpm,fMmm);


%{
fM=zeros(aMAX,aMAX);

for a1=1:aMAX
    for a2=1:aMAX
        [m1,n1] = mnAlpha(a1,maxx);
        [m2,n2] = mnAlpha(a2,maxx);
        m = (m2-m1); % puvodne zde neni -()
        n = (n2-n1); % puvodne zde neni -()
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
%}


