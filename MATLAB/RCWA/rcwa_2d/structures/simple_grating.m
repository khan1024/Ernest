number_of_layers=2;

%l=1;
%thickness(l)=0.2;

%px(:,:,l)=1.2*[0, 0.25, 0.75, 1];
%py(:,:,l)=1.2*[0, 0.25, 0.75, 1];

%px(:,:,l)=[0, 0.005, 0.015, 0.02];
%py(:,:,l)=[0, 0.005, 0.055, 0.06];
%n_mat(:,:,l)=[1,1,1; 1,1.5,1; 1,1,1];

%l=1;
%thickness(l)=0.2;
%zl=0.23-6.71*1i;
%px(:,:,l)=[0, 0.45-0.125, 0.45+0.125, 0.9];
%py(:,:,l)=[0, 0.45-0.125, 0.45+0.125, 0.9];
%n_mat(:,:,l)=[zl,zl,zl; zl,1.0,zl; zl,zl,zl];

%
l=1;
thickness(l)=1.0;

co=3.5;
px(:,:,l)=[0,0.5,1.5,1.9,2.9,3.4];
py(:,:,l)=[0,0.5,1.5,1.9,2.9,3.4];
n_mat(:,:,l)=[co, co, co, co, co;
    co, co, co, co, co;
    co, co, co, co, co;
    co, co, co, co, co;
    co, co, co, co, co];

l=2;
thickness(l)=1.0;

coo=2.5;
px(:,:,l)=[0,0.5,1.5,1.9,2.9,3.4];
py(:,:,l)=[0,0.5,1.5,1.9,2.9,3.4];
n_mat(:,:,l)=[coo, coo, coo, coo, coo;
    coo, coo, coo, coo, coo;
    coo, coo, coo, coo, coo;
    coo, coo, coo, coo, coo;
    coo, coo, coo, coo, coo];

%{
l=1;
thickness(l)=.2;

px(:,:,l)=0.5*[0, 0.2, 0.8, 1];
py(:,:,l)=0.5*[0, 0.2, 0.8, 1];
n_mat(:,:,l)=[1,1,1; 1,sqrt(3.18-4.41*1i),1; 1,1,1];
%}

%{
l=1;
thickness(l)=1;

px(:,:,l)=0.1*[0, 0.45, 0.55, 1];
py(:,:,l)=0.1*[0, 0.1, 0.9, 1];
n_mat(:,:,l)=[1,1,1; 1,2,1; 1,1,1];
%}