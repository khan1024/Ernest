% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%nakreslit_graf=1;
%Lambda=3;
%d=5;
%number_of_layers=300;

Lambda=Lambda/2;
d=thickness_total;
sirka=thickness_total*1E-6;

rozdeleni=zeros(1,number_of_layers+1);
nova_y=zeros(1,number_of_layers);
pomocna=zeros(1,number_of_layers);
duty_cycle=zeros(1,number_of_layers);
shift=zeros(1,number_of_layers);


rozdeleni=linspace(0,d,number_of_layers+1);%rozdělím zadaný interval
for i=1:1:number_of_layers %udělám průměr dvou soesedních y
    nova_y(i)=(rozdeleni(i)+rozdeleni(i+1))/2;
end
pomocna=sqrt(-2*(Lambda^2/pi^2)*log(nova_y/d));%inverzní fuknkce


for j=1:1:number_of_layers %omezeni velikosti vrstviček
    if pomocna(j)>Lambda
        pomocna(j)=Lambda;
    else
        pomocna(j)=pomocna(j);
    end
end

duty_cycle=fliplr(2*pomocna/Lambda);
layer_thickness=sirka/number_of_layers*ones(1,number_of_layers);%případně zakomentovat
shift=.5*ones(1,number_of_layers);

switch plot_grating
    case 1
        figure;
        x=-Lambda:Lambda/100:Lambda;
        y1=d*exp(-pi*pi*x.*x/(2*Lambda*Lambda));
        hold on
        plot(x,y1,'r','LineWidth',2);
        line([-Lambda,Lambda],[0,0],'Color','r','LineWidth',2);
        xlabel('\Lambda [\mu{}m]','FontSize',18);
        ylabel('d [\mu{}m]','FontSize',18);
        for k=1:1:number_of_layers
            z=-pomocna(k):Lambda/100:pomocna(k);
            plot(z,0*z+rozdeleni(k))
            plot(z,0*z+rozdeleni(k)+d/number_of_layers)
            line([pomocna(k);pomocna(k)],[nova_y(k)-d/(2*number_of_layers);nova_y(k)+d/(2*number_of_layers)])
            line([-pomocna(k);-pomocna(k)],[nova_y(k)-d/(2*number_of_layers);nova_y(k)+d/(2*number_of_layers)])
        end
    case 2
end