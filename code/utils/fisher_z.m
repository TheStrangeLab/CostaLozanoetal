%%%%%%%%%%% do fisher z transformation
%%%%%%%%%%% Stephan Moratti, UCM, Departamento de Psicología Básica I,
%%%%%%%%%%% Madrid Spain

function [z] = fisher_z(data);


z = (0.5*log((1+data)./(1-data)));

