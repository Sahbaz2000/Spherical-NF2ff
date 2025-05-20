function [data_nf_raw] = NearfieldDatagenerator(ant,freq,r_nf,theta_nf,phi_nf)

% Diese Funktion simuliert Nahfelddaten für eine vorgegebene Antenne (ant)
% bei einer bestimmten Frequenz (freq). Die Messpunkte liegen auf einer 
% kugelförmigen Oberfläche mit dem Radius r_nf um die Antenne. 
% Zur Abdeckung der Kugel wird ein Gitter aus Theta- und Phi-Winkeln
% (theta_nf & phi_nf) 
% verwendet, sodass die Nahfelddaten flächendeckend über die gesamte 
% Kugeloberfläche erfasst werden können.
% 
% Input Arguments:
%
%       ant          Antenne als Objekt (Antenna Toolbox)
%       freq         Frequenz in Hz
%       r_nf         Radius in m
%       theta_nf     Theta angle in Grad
%       phi_nf       Phi angle in Grad
%
%
% Output Arguments:
%
%       data_nf_raw      Tabelle Data_nf_raw mit 11 Spalten:
%                        X, Y, Z, ExReal, ExImg, EyReal, EyImg, EzReal, EzImg, EabsReal, EabsImg
%
% Die Output Tabelle enthält die kartesichen Koordianten X,Y & Z aus
% sphärischen Koordinaten theta_nf, phi_nf und r_nf. Die elektrischen 
% Feldkomponenten in x,y und z Richtung unterteilt in Real und Imaginärteil.
% Zum Schluss noch den Gesamtbetrag des elektrischen Feldes in EabsReal 

figure (Name='Antenne in der Antenna Toolbox'); show(ant);
%% Nahfeldpunkte definieren
% Für eine sphärische NF2FF-Transformation brauchen wir Nahfelddaten,
% die wir auf einem Kugeloberfläche um die Antenne messen.

% Erzeuge alle (theta, phi)-Kombinationen.
[ThetaGrid_nf, PhiGrid_nf] = meshgrid(theta_nf, phi_nf);

% In Vektoren bringen
theta_vec = ThetaGrid_nf(:);
phi_vec   = PhiGrid_nf(:);

% Nun nach Kartesisch umrechnen:
x_nf = r_nf * sind(theta_vec) .* cosd(phi_vec);
y_nf = r_nf * sind(theta_vec) .* sind(phi_vec);
z_nf = r_nf * cosd(theta_vec);

% Zusammenstellen aller Nahfeldkoordinaten
nfPoints = [x_nf, y_nf, z_nf].';
figure (Name='Messpunkte auf der Kugeloberfläche für die Nahfeldmessung')
plot3(nfPoints(1,:), nfPoints(2,:), nfPoints(3,:), 'x');
title('Messpunkte auf der Kugeloberfläche für die Nahfeldmessung')
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
grid on;
axis equal;

%% Nahfelddaten simulieren (E und H) 
% EHfields(Antenne, Frequenz, [x,y,z]) gibt die komplexen E/H-Felder zurück.

[E_nf, H_nf] = EHfields(ant, freq, nfPoints); %#ok

% EHfields liefert E_nf(1,i), E_nf(2,i), E_nf(3,i) als komplexe Werte (Ex, Ey, Ez).
% Wir benötigen Real- und Imaginärteile getrennt 

Ex_nf = E_nf(1,:);
Ey_nf = E_nf(2,:);
Ez_nf = E_nf(3,:);

% Betrag des E-Feldes
Eabs_nf = sqrt(abs(Ex_nf).^2 + abs(Ey_nf).^2 + abs(Ez_nf).^2);

%% Nahfelddaten in Tabellenformat bringen
% Damit "rearrangeTables" funktioniert, erwartet es üblicherweise Spalten:
% X, Y, Z, ExReal, ExImg, EyReal, EyImg, EzReal, EzImg, EabsReal, EabsImg

% Wir bauen uns eine Tabelle "data_nf_raw", so wie sie
% "rearrangeTables.m"  erwarten würde.

data_nf_raw = table();
data_nf_raw.X = x_nf;
data_nf_raw.Y = y_nf;
data_nf_raw.Z = z_nf;

data_nf_raw.ExReal = real(Ex_nf)';
data_nf_raw.ExImg  = imag(Ex_nf)';
data_nf_raw.EyReal = real(Ey_nf)';
data_nf_raw.EyImg  = imag(Ey_nf)';
data_nf_raw.EzReal = real(Ez_nf)';
data_nf_raw.EzImg  = imag(Ez_nf)';

data_nf_raw.EabsReal = real(Eabs_nf)';

%% Test ob Schrittweite klein genug dass zwischen 2 Punkten max Abstand von Lambda/2 ist Beginn
% Um das Abtasttheorem nicht zu verletzen darf der Abstand maximal zwischen
% 2 Punkten Lambda/2 sein
lambda = 3e8/ freq;
d_max = 0.45 * lambda;  % Sicherheitsfaktor

% Dmax der maximale Abstand zwischen 2 Punkten in Theta und Phi Richtung
delta_theta = theta_nf(3) - theta_nf(2);
Dmax_theta = 2.*r_nf.* sind(delta_theta./2);
delta_phi = phi_nf(3) - phi_nf(2);
Dmax_phi = 2.*r_nf.* sind(delta_phi./2);

if Dmax_phi < Dmax_theta
    Dmax = Dmax_theta;
else
    Dmax = Dmax_phi;
end

% Prüfen ob der größte Abstand dennoch kleiner ist als lamda/2 wenn nicht
% muss die Winkel Schrittweite kleiner gewählt werden
if Dmax < d_max
     disp('✅ Theta & Phi Abstand klein genug gewählt / Max Abstand zwischen 2 punkten darf Lambda/2 sein');
else
  disp('❌ Theta & Phi Abstand zu groß gewählt / Max Abstand zwischen 2 punkten darf Lambda/2 sein');
end
%% Test ob Schrittweite klein genug dass zwischen 2 Punkten max Abstand von Lambda/2 ist Ende
end