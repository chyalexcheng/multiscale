Point(1) = {0, 0, 0, 0.12};
Point(2) = {1.2, 0, 0, 0.12};
Point(3) = {1.2, 1.2, 0, 0.12};
Point(4) = {0, 1.2, 0, 0.12};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Line {1, 2, 3, 4} = 7;
Transfinite Surface {1} Left;
Recombine Surface {1};

Physical Line("bound") = {1,2,3,4};
Physical Surface("dom") = {1};
