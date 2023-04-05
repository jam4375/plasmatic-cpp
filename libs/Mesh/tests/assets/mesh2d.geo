// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 1
//
//  Geometry basics, elementary entities, physical groups
//
// -----------------------------------------------------------------------------

// The simplest construction in Gmsh's scripting language is the
// `affectation'. The following command defines a new variable `lc':

lc = 1e-1;

Point(1) = {0, 0, 0, lc};
Point(2) = {.1, 0,  0, lc};
Point(3) = {.1, .3, 0, lc};
Point(4) = {0,  .3, 0, lc};

Line(1) = {1, 2};
Line(2) = {3, 2};
Line(3) = {3, 4};
Line(4) = {4, 1};


Curve Loop(1) = {4, 1, -2, 3};


Plane Surface(1) = {1};

Physical Curve("My curve 1") = {1, 2, 4};
Physical Curve("My curve 2") = {3};
Physical Surface("My surface") = {1};

