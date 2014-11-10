% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function [D_mat_1]=diffraction_layout(lambda,n1,n3,Lambda,theta0)
figure;

lambda=lambda*1E-6;
thetai=theta0;
perioda=Lambda*1E-6;

p=1;

q=ceil(perioda*(n1*sin(thetai*(pi/180))-n3)/lambda):1:floor(perioda*(n1*sin(thetai*(pi/180))+n3)/lambda);

thetad=-asin((-n1*sin(thetai*(pi/180))*perioda+q*lambda)/n3/perioda)*180/pi;
thetar=-asin(n3*sin(thetad*pi/180)/n1)*180/pi;
D_mat=[thetad' q' thetar'];
D_mat_1=rot90(rot90([-thetad' -q' -thetar']));

%plot the figure
if p==1
	kt=2*pi*n3/lambda;
	ki=2*pi*n1/lambda;
    
    k_0 = 2*pi/lambda;
    k_1 = k_0*n1;
    k_2 = k_0*n3;
    
    for ord = -80:80
      kx=k_1*sin(pi*thetai/180)+ord*2*pi/perioda;
      ky=sqrt(k_1.^2-kx^2);
      if (imag(ky)==0 )
          arrow([0 0],[kx ky],.07,[1,0,0]);
          text(1.1*kx,1.1*ky,num2str(ord))
      end
    end
    kx = 1.2*k_1*sin(-pi*thetai/180);
    ky = 1.2*k_1*cos(-pi*thetai/180);
    arrow([kx ky],[0 0],.07,[0,0,1]);
  
  for ord=-80:80
      kx = k_1*sin(pi*thetai/180)+ord*2*pi/perioda;
      ky = sqrt(k_2^2-kx^2);
      if (imag(ky)==0 )
          arrow([0 0],[kx -ky],.07,[1,0,0]);
      end
  end
  
    
    for v=1:length(q)
		pl=line([-(q(v))*2*pi/perioda+kt*sin(pi/180*thetad(find(q==0)))...
			-(q(v))*2*pi/perioda+kt*sin(pi/180*thetad(find(q==0)))],...
			[-kt*cos(pi/180*thetad(find(q==0))) ki*cos(pi/180*thetai) ]);
		set(pl,'Linestyle',':','Color','m')
		
		text(1.1*kt*sin(thetad(v)*pi/180),-1.1*kt*cos(thetad(v)*pi/180)...
			,num2str(-q(v)),'HorizontalAlignment','center')
	end
	hold on;
	plot(linspace(-ki,ki,1000),sqrt(ki^2-(linspace(-ki,ki,1000)).^2),'g')% indexova plocha n_1
	plot(linspace(-kt,kt,1000),-sqrt(kt^2-(linspace(-kt,kt,1000)).^2),'g')% indexova plocha n_2
	plot([-max(ki,kt) max(ki,kt)],[0 0],'g'); %osa x
	axis off;
	title(['Vector diagram for \lambda = ' num2str(lambda) 'm and \Lambda =' num2str(perioda) 'm'],'FontSize',14)
end