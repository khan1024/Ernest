function em_pole_procedura=zobrazeni_pole(em_pole,procedura_pole)

if procedura_pole==1
    em_pole_procedura=abs(em_pole)./376.76;
elseif procedura_pole==2
    em_pole_procedura=angle(em_pole);
elseif procedura_pole==3
    em_pole_procedura=real(em_pole);
elseif procedura_pole==4
    em_pole_procedura=imag(em_pole);
elseif procedura_pole==5
    em_pole_procedura=em_pole;
end