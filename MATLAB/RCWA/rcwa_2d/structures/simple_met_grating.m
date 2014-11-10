number_of_layers=1;


l=1; % layer 1
thickness(l)=0.25;

g_Lambda=0.5;
g_d=0.286;
g_thickness=0.25;

thickness(1)=g_thickness;

%coordinate_x=[0,0.75,1.3,2.1,2.65,3.4];
px_1=[0,(g_Lambda-g_d)/2,(g_Lambda+g_d)/2,g_Lambda];
py_1=[0,(g_Lambda-g_d)/2,(g_Lambda+g_d)/2,g_Lambda];


n_mat_1=[al_h, al_h, al_h;
    al_h, 1, al_h;
    al_h, al_h, al_h];

%{
l=2; % layer 2
thickness(l)=0.5;

px_2=1.2.*[0,0.25,0.75,1];
py_2=1.2.*[0,0.25,0.75,1];

n_mat_2=[1, 1, 1;
    1, 1.5, 1;
    1, 1, 1];

l=3; % layer 2
thickness(l)=0.5;

px_3=1.2.*[0,1];
py_3=1.2.*[0,1];

n_mat_3=[1.5];
%}
