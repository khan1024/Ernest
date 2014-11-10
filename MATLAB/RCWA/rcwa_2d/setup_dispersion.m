if use_dispersion==1
    %p_dr=[216,72500]; %gold
    %p_dr=[145,72500]; %stribo
    %c_he=[1,-0.00309946,-9.61396e-006,1.72005,0,0.00561194,0,0,0,0,0.028,0,0,0,0,-1.09862e-005];
    %c_se=[0.5,1,1.03961,0.231792,1.01147,0,0,0.0060007,0.0200179,103.561,0,0];
    %c_sc=[2.27189,-0.0101081,0.0105925,0.00020817,-7.64725e-006,4.24099e-007];
    %c_po=[1.5,0.2,0,0,0,0,0,0,0,0,-0.1,0.095,-0.1,0,0,0,0,0,0,0;];

    %gold=drude(lambda,p_dr);
    %herzberger(lambda,c_he);
    %sellmeier(lambda,c_se);
    %schott(lambda,c_sc);
    %polynomial(lambda,c_po);
    %n3=tabulka(lambda,'gold.txt');
    %nr=n3;
    %gold=0.21332132-3.0131975*j;

    %a=1.147;
    %b=2.578;
    %c=0.1472;

    %b_0=2.979864;
    %b_1=8.777808*1E-3;
    %b_2=1.0609*1E-2;
    %b_3=84.06224;
    %b_4=-96;

    %n_g=sqrt(a^2+b*lambda.*lambda./(lambda.*lambda-c^2));
    %n_s=sqrt(b_0+b_1./(lambda.*lambda+b_2)+b_3./(lambda.*lambda+b_4));

    %index_kov=3.5;
    %index_kov=0.537000000-9.580000000*j;
    %gold=tabulka(lambda,'lorentz_drude_au.txt');
    %gold=tabulka(lambda,'clanek2.txt');
    %gold=cepau(lambda);
    %au_cepau=cepau(lambda);
    %p_dr=[216,72500];
    %gold=drude(lambda,p_dr);
    %stribro=tabulka(lambda,'stribro_pcgrate.txt');
    %index_kov=tabulka(lambda,'gold_pcgrate.txt');
    %gold=tabulka(lambda,'gold.txt');
    
    
    
    %       material ==>    'Ag'  = silver
    %                       'Al'  = aluminum
    %                       'Au'  = gold
    %                       'Cu'  = copper
    %                       'Cr'  = chromium
    %                       'Ni'  = nickel
    %                       'W'   = tungsten
    %                       'Ti'  = titanium
    %                       'Be'  = beryllium
    %                       'Pd'  = palladium
    %                       'Pt'  = platinum
    %                       'H2O' = pure water (triply distilled)
    weiss_drude_eps=Drude_model_weiss(lambda*1e-6);
    n_weiss_drude=conj(sqrt(weiss_drude_eps));

    %al_h=rix_spline(lambda*1000,'al_pcgrate.txt');
    gold=rix_spline(lambda,'gold_j_a_ch.txt');
    %ag_h=Drude_model_fun(lambda*1e-6);
    %{
    number_of_poles=1; % 1-Drude, 2-5(6)-Lorentz-Drude
    damping=1; % 1-ano, 2-ne
    imag_sign=2; % 1-imag=+, 2-imag=-
    [real_epsilon,imag_epsilon,gold]=new_LD(lambda,'Au',number_of_poles,damping,imag_sign);
    %}
end