% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%p_dr=[216,72500];
%c_he=[1,-0.00309946,-9.61396e-006,1.72005,0,0.00561194,0,0,0,0,0.028,0,0,0,0,-1.09862e-005];
%c_se=[0.5,1,1.03961,0.231792,1.01147,0,0,0.0060007,0.0200179,103.561,0,0];
%c_sc=[2.27189,-0.0101081,0.0105925,0.00020817,-7.64725e-006,4.24099e-007];
%c_po=[1.5,0.2,0,0,0,0,0,0,0,0,-0.1,0.095,-0.1,0,0,0,0,0,0,0;];

%zlato=drude(lambda,p_dr);
%herzberger(lambda,c_he);
%sellmeier(lambda,c_se);
%schott(lambda,c_sc);
%polynomial(lambda,c_po);
silica=cauchy(lambda,1.491,0.00686,-0.0007648);

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

gold=rix_spline(lambda,'gold_palik.txt');
