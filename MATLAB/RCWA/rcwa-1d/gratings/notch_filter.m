notch_thickness_g=0.01;
notch_thickness_s=4;
notch_duty=0.3;
notch_L=0.5;

data=[
notch_thickness_s, 0, silica, notch_L, NaN, NaN;
notch_thickness_g, 0, gold, notch_L*notch_duty, silica, notch_L;
];
n1=1.9;
n3=silica;