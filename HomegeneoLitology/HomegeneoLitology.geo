//--------------------------------------------------------------------
//---------------- UNIVERSIDADE FEDERAL DE PERNAMBUCO ----------------
//---------------- CENTRO DE TECNOLOGIA E GEOCIENCIAS ----------------
//--------------------------------------------------------------------

//Work developed by: Marcio Souza, Luiz E. Queiroz e Fernando Contreras
//Adviser Professors: Paulo Lyra & Darlan Carvalho
//Create date: 13/12/2011
//Create date: 2022/3/24;	hour: 20:24h

//--------------------------------------------------------------------
//This file has CAD parameters. It is related to building of domain

cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Physical Point(101) = {1, 4};
Physical Point(102) = {2, 3};
Physical Line(201) = {1, 3};
Physical Line(101) = {4};
Physical Line(102) = {2};
Physical Surface(1) = {6};


Transfinite Line {1,3} = 3 Using Progression 1.000000;
Transfinite Line {2,4} = 2 Using Progression 1.000000;
Transfinite Surface {6} = {1,2,3,4};
