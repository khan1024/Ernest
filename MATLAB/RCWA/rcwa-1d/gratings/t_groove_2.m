h_1=0.3;
h_2=0.35;
W_1=0.2;
W_2=0.05;
P=1.1;
d=0.03;

gap=(P-4*W_2)/4;
b_1=W_2/2;
b_2=b_1+gap;
b_3=b_2+W_2;
b_4=b_3+gap;
b_5=b_4+W_2;
b_6=b_5+gap;
b_7=b_6+W_2;
b_8=b_7+gap;

data=[
h_1, 0, silver, P/2-W_1/2, 1, P/2+W_1/2, silver, P, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
d, 0, 1, P, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
h_2, 0, 1, b_1, silver, b_2, 1, b_3, silver, b_4, 1, b_5, silver, b_6, 1, b_7, silver, b_8, 1, P;
];

n1=1;
n3=silver;