% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%d=4;
%Lambda=7;
%number_of_layers=20;


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
pomocna=Lambda/(d)*nova_y;%inverzní fuknkce

duty_cycle=fliplr((Lambda-pomocna)/Lambda);
layer_thickness=sirka/number_of_layers*ones(1,number_of_layers);%případně zakomentovat
% z jedne strany shift=fliplr(pomocna/(2*Lambda));
shift=fliplr((pomocna/2+Lambda/2)/Lambda);

switch plot_grating
    case 1
        figure;
        x=0:Lambda/100:Lambda;
        y1=d/Lambda*x;
        hold on
        plot(x,y1,'r','LineWidth',2);
        line([0,Lambda],[0,0],'Color','r','LineWidth',2);
        line([Lambda,Lambda],[0,d],'Color','r','LineWidth',2);
        xlabel('\Lambda [\mu{}m]','FontSize',18);
        ylabel('d [\mu{}m]','FontSize',18);
        for k=1:1:number_of_layers
            z=Lambda:-Lambda/100:pomocna(k);
            plot(z,0*z+rozdeleni(k))
            plot(z,0*z+rozdeleni(k)+d/number_of_layers)
            line([pomocna(k);pomocna(k)],[nova_y(k)-d/(2*number_of_layers);nova_y(k)+d/(2*number_of_layers)])
            %line([-pomocna(k);-pomocna(k)],[nova_y(k)-d/(2*number_of_layers);nova_y(k)+d/(2*number_of_layers)])
        end
    case 2
end