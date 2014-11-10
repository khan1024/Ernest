clear all
lambda_min=0.000125;
lambda_max=9.919;
lambda_krok=0.001;

% vyber model
model=6; % 1-Drude, 2-Herzberger, 3-Sellmeier, 4-Schott, 5-Polynom, 6-Tabulka

% nastav parametry

parametry=[216,72500];

% v případě výpočtu z tabulky je nutné zadat zdrojová data
material='gold_pcgrate.txt';

% výpočet
lambda=lambda_min:lambda_krok:lambda_max;
velikost=size(lambda);
n_re=zeros(1,velikost(2));
n_im=zeros(1,velikost(2));

for i=1:1:velikost(2)
    switch model
        case 1
            n_re(i)=real(drude(lambda(i),parametry));
            n_im(i)=-imag(drude(lambda(i),parametry));
            titulek_re=['Drudeho model Re(n)'];
            titulek_im=['Drudeho model Im(n)'];
        case 2
            n_re(i)=real(herzberger(lambda(i),parametry));
            n_im(i)=-imag(herzberger(lambda(i),parametry));
            titulek_re=['Herzbergeruv model Re(n)'];
            titulek_im=['Herzbergeruv model Im(n)'];
        case 3
            n_re(i)=real(sellmeier(lambda(i),parametry));
            n_im(i)=-imag(sellmeier(lambda(i),parametry));
            titulek_re=['Sellmeieruv model Re(n)'];
            titulek_im=['Sellmeieruv model Im(n)'];
        case 4
            n_re(i)=real(schott(lambda(i),parametry));
            n_im(i)=-imag(schott(lambda(i),parametry));
            titulek_re=['Schottuv model Re(n)'];
            titulek_im=['Schottuv model Im(n)'];
        case 5
            n_re(i)=real(polynomial(lambda(i),parametry));
            n_im(i)=-imag(polynomial(lambda(i),parametry));
            titulek_re=['Polynomialni model Re(n)'];
            titulek_im=['Polynomialni Im(n)'];
        case 6
            n_re(i)=real(tabulka(lambda(i),material));
            n_im(i)=-imag(tabulka(lambda(i),material));
            titulek_re=['Data z tabulky Re(n)'];
            titulek_im=['Data z tabulky Im(n)'];
    end
end

switch model
    case 1
        titulek_re=['Drudeho model Re(n)'];
        titulek_im=['Drudeho model Im(n)'];
    case 2
        titulek_re=['Herzbergeruv model Re(n)'];
        titulek_im=['Herzbergeruv model Im(n)'];
    case 3
        titulek_re=['Sellmeieruv model Re(n)'];
        titulek_im=['Sellmeieruv model Im(n)'];
    case 4
        titulek_re=['Schottuv model Re(n)'];
        titulek_im=['Schottuv model Im(n)'];
    case 5
        titulek_re=['Polynomialni model Re(n)'];
        titulek_im=['Polynomialni Im(n)'];
    case 6
        titulek_re=['Data z tabulky Re(n)'];
        titulek_im=['Data z tabulky Im(n)'];
end

plot(lambda,n_re);
xlabel('\lambda [\mu m]','FontSize',14);
ylabel('Re(n)','FontSize',14);
title(titulek_re,'FontSize',15);

figure;
plot(lambda,n_im);
xlabel('\lambda [\mu m]','FontSize',14);
ylabel('Im(n)','FontSize',14);
title(titulek_im,'FontSize',15);