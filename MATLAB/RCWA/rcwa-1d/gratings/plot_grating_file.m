
velikost_vysledna_matice_obr=size(data);
vysledna_matice_obr=data;

y_max=0;
for i=1:1:velikost_vysledna_matice_obr(1)
    y_max=y_max+vysledna_matice_obr(i,1);
end

y_predesla=0;
for i=1:1:velikost_vysledna_matice_obr(1)
    for j=2:2:velikost_vysledna_matice_obr(2)-1
        if isnan(vysledna_matice_obr(i,j+2))==1
            
        else
            x_pocatek=vysledna_matice_obr(i,j);
            x_konec=vysledna_matice_obr(i,j+2);
            x_width=x_konec-x_pocatek;
            y_pocatek=y_max-vysledna_matice_obr(i,1)-y_predesla;
            y_height=vysledna_matice_obr(i,1);
            if y_height==0
                y_height=1e-30;
            end
            barva=[abs(cos(3.5*real(vysledna_matice_obr(i,j+1)))),abs(cos(5.26*real(vysledna_matice_obr(i,j+1)))),abs(sin(35*real(vysledna_matice_obr(i,j+1))))];
            if x_width>0
                rectangle('Position',[x_pocatek,y_pocatek,x_width,y_height],'FaceColor',barva,'LineStyle','none')
                %rectangle('Position',[x_pocatek,y_pocatek,x_width,y_height],'LineWidth',2,'FaceColor','none')
            end
        end
    end
    y_predesla=y_predesla+vysledna_matice_obr(i,1);
end

%axis off
%plot2svg('sensor_mesh.svg', 1);

xlabel('\Lambda [\mu{}m]','FontSize',17);
ylabel('d [\mu{}m]','FontSize',17);