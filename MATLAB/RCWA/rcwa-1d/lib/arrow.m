function l = arrow(p1,p2,headlength,color)
%function l = arrow(p1,p2,headlength)
%
% function plotting a black arrow from point p1 to point p2
% headlength specifies the relative length of the head-triangle
% compared to the arrow-length
% the axis should have equal aspect ratio, otherwise the arrows
% look slant

% Bernard Haasdonk 1.8.2005

ang = 150*2*pi/360;
rm = [cos(ang), sin(ang); -sin(ang), cos(ang)];

l1 = line([p1(1),p2(1)],[p1(2),p2(2)]);
set(l1,'Color',color,'Linewidth',2);

p3 = rm *  (p2(:)-p1(:));
p4 = rm' * (p2(:)-p1(:));
p3 = p3 *headlength + p2(:);
p4 = p4 *headlength + p2(:);
p5 = p2(:)/4+p3*3/8+p4*3/8;

l2 = patch([p2(1),p3(1),p5(1),p4(1)],[p2(2),p3(2),p5(2),p4(2)],color);
l = [l1,l2];



