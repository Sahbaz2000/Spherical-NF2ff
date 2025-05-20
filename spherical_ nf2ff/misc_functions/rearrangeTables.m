function [data_nf] = rearrangeTables(data_nf)
% Diese Funktion strukturiert die Tabelle data_nf um, um den Output für die
% Funktion rotateSphericalData vorzuberaiten.
% 
% Input Arguments:
%
%       data_nf          Tabelle der Daten der Nahfeldaten mit 11 Spalten:
%                        X, Y, Z, ExReal, ExImg, EyReal, EyImg, EzReal, 
%                        EzImg, EabsReal, EabsImg
% 
%                        siehe Funktion NearfieldDatagenerator für 
%                        mehr  details des Aufbaus der Tabelle
%     
%
%
% Output Arguments:
%
%       "neue" data_nf     Tabelle mit 5 Spalten:
%                           X, Y, Z, E, Eabs
%
% Die Ausgabetabelle enthält die kartesischen Koordinaten X, Y und Z, die 
% unverändert aus der Eingabetabelle übernommen werden. Die elektrischen 
% Feldkomponenten in x-, y- und z-Richtung (Real- und Imaginärteil) wurden 
% zusammengeführt, sodass E(:,1) für Ex, E(:,2) für Ey und E(:,3) für Ez steht.
% Schließlich wird der Gesamtbetrag des elektrischen Feldes als Eabs ausgegeben.

Ex = data_nf.ExReal+ 1j*data_nf.ExImg;
Ey = data_nf.EyReal+ 1j*data_nf.EyImg;
Ez = data_nf.EzReal+ 1j*data_nf.EzImg;

% EabsImg sollte immer 0 sein 
Eabs = data_nf.EabsReal;

data_nf = table(data_nf.X,data_nf.Y,data_nf.Z,[Ex,Ey,Ez],Eabs);
data_nf.Properties.VariableNames = {'x' 'y' 'z' 'E' 'Eabs'};

% ALter Code hat xyz durch 1000 geteilt vermutlich da seine Messergebnisse
% in mm waren und nun m nötig sind
% data_nf.x = data_nf.x/1000;
% data_nf.y = data_nf.y/1000;
% data_nf.z = data_nf.z/1000;
end

