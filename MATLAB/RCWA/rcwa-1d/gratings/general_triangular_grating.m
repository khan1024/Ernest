% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%d=6;
%Lambda=6;
%number_of_layers=20;
%nakreslit_graf=1;
%a=.3;
%b=1;

d=thickness_total;
sirka=thickness_total*1E-6;

% definice vrcholů trojúhelníku
y_1=0;
x_1=(Lambda-b*Lambda)/2;

y_2=0;
x_2=(Lambda+b*Lambda)/2;

y_3=d;
x_3=a*Lambda;

% alokace matic
rozdeleni=zeros(1,number_of_layers+1);
nova_y=zeros(1,number_of_layers);
pomocna_1=zeros(1,number_of_layers);
pomocna_2=zeros(1,number_of_layers);
duty_cycle=zeros(1,number_of_layers);
shift=zeros(1,number_of_layers);

rozdeleni=linspace(0,d,number_of_layers+1);%rozdělím zadaný interval
for i=1:1:number_of_layers %udělám průměr dvou sousedních y
    nova_y(i)=(rozdeleni(i)+rozdeleni(i+1))/2;
end

% případ blejzované mřížky je potřeba ošetřit (dělení 0) 
if x_3-x_1==0
    pomocna_1=x_1*ones(1,number_of_layers);
else
    pomocna_1=(nova_y-(y_1*x_3-y_3*x_1)/(x_3-x_1))*((x_1-x_3)/(y_1-y_3));
end

if x_3-x_2==0
    pomocna_2=x_2*ones(1,number_of_layers);
else
    pomocna_2=(nova_y-(y_2*x_3-y_3*x_2)/(x_3-x_2))*((x_2-x_3)/(y_2-y_3));
end

duty_cycle=fliplr((pomocna_2-pomocna_1)/Lambda);
layer_thickness=sirka/number_of_layers*ones(1,number_of_layers);%případně zakomentovat
shift=fliplr((pomocna_1+pomocna_2)/(2*Lambda));

switch plot_grating
    case 1
        figure;
        line([x_1,x_3],[y_1,y_3],'Color','r','LineWidth',2);
        line([x_2,x_3],[y_2,y_3],'Color','r','LineWidth',2);
        line([x_1,x_2],[y_1,y_2],'Color','r','LineWidth',2);
        axis([0 Lambda 0 d])
        hold on
        xlabel('\Lambda [\mu{}m]','FontSize',18);
        ylabel('d [\mu{}m]','FontSize',18);
        for k=1:1:number_of_layers
            z=pomocna_1(k):Lambda/1000:pomocna_2(k);
            plot(z,0*z+rozdeleni(k))
            plot(z,0*z+rozdeleni(k)+d/number_of_layers)
            line([pomocna_1(k);pomocna_1(k)],[nova_y(k)-d/(2*number_of_layers);nova_y(k)+d/(2*number_of_layers)]);
            line([pomocna_2(k);pomocna_2(k)],[nova_y(k)-d/(2*number_of_layers);nova_y(k)+d/(2*number_of_layers)]);
        end
    case 2
end