function [pnt, dhk] = tetrahedron()

% TETRAHEDRON creates a tetrahedron
%
% [pnt, dhk] = tetrahedron
% creates a tetrahedron with 4 vertices and 4 triangles

dhk = [
   1   2   3
   1   3   4
   1   4   2
   2   3   4
];

pnt = zeros(4, 3);

pnt(1, :) = [sqrt(8/9), 0, -1/3];
pnt(2, :) = [-sqrt(2/9), sqrt(2/3), -1/3];
pnt(3, :) = [-sqrt(2/9), -sqrt(2/3), -1/3];
pnt(4, :) = [0, 0, 1];

end
