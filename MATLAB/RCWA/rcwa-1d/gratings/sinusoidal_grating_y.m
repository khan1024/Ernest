% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%d=4;
%Lambda=4;
%number_of_layers=5;

d=thickness_total;
sirka=thickness_total*1E-6;

rozdeleni=zeros(1,number_of_layers+1);
nova_y=zeros(1,number_of_layers);
pomocna=zeros(1,number_of_layers);
duty_cycle=zeros(1,number_of_layers);
shift=zeros(1,number_of_layers);

rozdeleni=linspace(-d/2,d/2,number_of_layers+1);%rozdělím zadaný interval
for i=1:1:number_of_layers %udělám průměr dvou soesedních y
    nova_y(i)=(rozdeleni(i)+rozdeleni(i+1))/2;
end
pomocna=Lambda/(2*pi)*acos(2*nova_y/d);%inverzní fuknkce

duty_cycle=fliplr(2*pomocna/Lambda);
layer_thickness=sirka/number_of_layers*ones(1,number_of_layers);%případně zakomentovat
shift=.5*ones(1,number_of_layers);

switch plot_grating
    case 1
        figure;
        x=-Lambda/2:Lambda/100:Lambda/2;
        y1=(d/2)*cos(2*pi/Lambda*x);
        hold on
        plot(x,y1,'r','LineWidth',2);
        line([-Lambda/2,Lambda/2],[-d/2,-d/2],'Color','r','LineWidth',2);
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

