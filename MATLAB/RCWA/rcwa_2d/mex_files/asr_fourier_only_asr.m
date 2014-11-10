% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function [pred_asr]=asr_fourier_only_asr(disc_x,disc_u,M,d,G)

length_disc_x=length(disc_x);

% nastaveni x a x_m1
x_m1=zeros(1,length_disc_x-1);
x=zeros(1,length_disc_x-1);

% nastaveni u a u_m1
u_m1=zeros(1,length_disc_x-1);
u=zeros(1,length_disc_x-1);

for s=2:1:length_disc_x
    x_m1(1,s-1)=disc_x(1,s-1);
    x(1,s-1)=disc_x(1,s);    
    u_m1(1,s-1)=disc_u(1,s-1);
    u(1,s-1)=disc_u(1,s);
end

% nastaveni asr parametru a_1, a_2, a_3
a_1=zeros(1,length_disc_x-1);
a_2=zeros(1,length_disc_x-1);
a_3=zeros(1,length_disc_x-1);

for s=1:1:length_disc_x-1
    a_1(1,s)=(u(1,s)*x_m1(1,s)-u_m1(1,s)*x(1,s))/(u(1,s)-u_m1(1,s));
    a_2(1,s)=(x(1,s)-x_m1(1,s))/(u(1,s)-u_m1(1,s));    
    a_3(1,s)=G*(u(1,s)-u_m1(1,s))-(x(1,s)-x_m1(1,s));
end

pred_asr=zeros(1,2*M+1,length(x));

for i=1:1:length(x)
    a_dif=u(i)+u_m1(i);
    pred_asr(:,:,i)=(-1./2).*(-2.*u(i)+a_dif)./d./(u(i)-u_m1(i)).*(2.*exp(-1i.*pi.*(-M:M).*a_dif./d).*a_2(i).*sinc(-(-M:M).*(-2.*u(i)+...
        a_dif)./d).*u(i)-2.*exp(-1i.*pi.*(-M:M).*a_dif./d).*a_2(i).*sinc(-(-M:M).*(-2.*u(i)+a_dif)./d).*u_m1(i));
end