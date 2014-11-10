% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

thickness_total=2.24;
d=thickness_total;
Lambda=1;
number_of_layers=10;

rozdeleni=zeros(1,number_of_layers+1);
pomocna=zeros(1,number_of_layers+1);
duty_cycle=zeros(1,number_of_layers);

x=-Lambda/2:Lambda/100:Lambda/2;

rozdeleni=linspace(-d/2,d/2,number_of_layers+1);
pomocna=Lambda/(2*pi)*acos(2*rozdeleni/d);

for i=1:1:number_of_layers %udělám průměr dvou soesedních x
    nova(i)=(pomocna(i)+pomocna(i+1))/2;
end
duty_cycle=2*nova/Lambda;

hold on
y1=(d/2)*cos(2*pi/Lambda*x);
%y2=Lambda/(2*pi)*acos(2*x/d);
plot(x,y1,'r')
xlabel('\Lambda','FontSize',18);
ylabel('d','FontSize',18);
for k=1:1:number_of_layers
    z=-nova(k):.01:nova(k);
    plot(z,0*z+rozdeleni(k))
    plot(z,0*z+rozdeleni(k)+d/number_of_layers)
    %zz=-pomocna(k):.01:pomocna(k);
    %plot(zz,0*zz+rozdeleni(k),'r')
    line([nova(k),nova(k)],[rozdeleni(k),rozdeleni(k)+d/number_of_layers]);
    line([-nova(k),-nova(k)],[rozdeleni(k),rozdeleni(k)+d/number_of_layers]);
end
