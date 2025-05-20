function [data_nf_rotated] = rotateSphericalNFData(data_nf)
% Diese Funktion wandelt die kartesischen Koordinaten (x, y, z) in 
% sphärische Koordinaten um und konvertiert gleichzeitig die elektrischen 
% Feldkomponenten (Ex, Ey, Ez) in die radialen Größen Eradial, Etheta und Ephi.
% 
% Input Arguments:
%
%       data_nf             Tabelle  mit 5 Spalten:
%                           X, Y, Z, E, Eabs
%     
% Output Arguments:
%
%       data_nf_rotated     Tabelle Data_nf_raw mit 8 Spalten:
%                            X, Y, Z, r, theta, phi, E, Eabs
%
% Die Ausgabetabelle enthält die kartesischen Koordinaten X, Y und Z, die 
% unverändert aus der Eingabetabelle übernommen werden, sowie die daraus 
% umgewandelten Kugelkoordinaten r, θ (theta) und φ (phi). Die elektrischen 
% Feldkomponenten Ex, Ey und Ez werden ebenfalls in Er, Eθ und Eφ
% umgewandelt in eine Spalte E also [Er = E(:,1), Etheta = E(:,2), Ephi = E(:,3)
% Schließlich wird der Gesamtbetrag des elektrischen Feldes als Eabs ausgegeben.

% Erstelle Output Tabelle mit Nullen
p = length(data_nf.E);
data_nf_rotated = table(zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,1),zeros(p,3),zeros(p,1));
data_nf_rotated.Properties.VariableNames = {'x','y','z','r','theta','phi','E','Eabs'};

% Übergabe der Kartesischen Koordinaten
data_nf_rotated.x = data_nf.x;
data_nf_rotated.y = data_nf.y;
data_nf_rotated.z = data_nf.z;

% Umformung in Kugelkoordinaten
r = sqrt(data_nf.x.^2 + data_nf.y.^2 + data_nf.z.^2);
theta = acos(data_nf.z./r);
phi = atan2(data_nf.y,data_nf.x);

%% Test Kartesisch in Kugelkoordinaten Beginn
% Test ob die Umformung korrekt erfolgt ist
[az_test, el_test, r_test] = cart2sph(data_nf_rotated.x,data_nf_rotated.y,data_nf_rotated.z);
tol = 1e-9;
% Jetzt theta und phi im "mathematischen" Sinn:
phi_test   = az_test;             % Azimut
theta_test = (pi/2) - el_test;    % Polarwinkel
if all(abs(r - r_test) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz)/ Test Kartesisch in Kugelkoordinaten / Variable r');
else
    disp('❌ Es gibt Unterschiede / Test Kartesisch in Kugelkoordinaten / Variable r');
end

if all(abs(theta - theta_test) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz) / Test Kartesisch in Kugelkoordinaten / Variable theta');
else
    disp('❌ Es gibt Unterschiede / Test Kartesisch in Kugelkoordinaten / Variable theta');
end

if all(abs(phi - phi_test) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz) / Test Kartesisch in Kugelkoordinaten / Variable phi');
else
    disp('❌ Es gibt Unterschiede / Test Kartesisch in Kugelkoordinaten / Variable phi');
end
%% Test Kartesisch in Kugelkoordinaten Ende

% phi um 2pi (360°) verschieben da phi sonst 0-pi dann von 2pi-pi läuft
phi_new = mod(phi_test, 2*pi); % dadurch läuft phi von 0 bis 2pi

for i =1:length(data_nf.E)

    % R ist die Rotationsmatrix für die Efeld komponenten für die
    % Umwandlung von Ex, Ey & Ez zu Er, Etheta & Ephi
    R = [sin(theta(i))*cos(phi(i)) sin(theta(i))*sin(phi(i)) cos(theta(i));...
         cos(theta(i))*cos(phi(i)) cos(theta(i))*sin(phi(i)) -sin(theta(i));...
         -sin(phi(i)) cos(phi(i)) 0];
     
    %E(1)-> Radial,E(2)->Theta,E(3)->Phi
    data_nf_rotated.E(i,:) = (R*data_nf.E(i,:)')';
    data_nf_rotated.r(i) = r(i);
    data_nf_rotated.theta(i) = theta(i);
    data_nf_rotated.phi(i)= phi_new(i);
    data_nf_rotated.Eabs(i) = data_nf.Eabs(i);
    end
    
    %% Test Rotation Korrekt Beginn
    % Getestet ob die Rotationsmatrix korrekt ist
    % Es wird rotiert mithilfe der Formeln und mit den Ergebnissen der Rotationsmatrix verglichen 
    Ex_test = data_nf.E(:,1);
    Ey_test = data_nf.E(:,2);
    Ez_test = data_nf.E(:,3);
 
    E_r_test = Ex_test .* sin(theta_test) .* cos(phi_test) ...
    + Ey_test .* sin(theta_test) .* sin(phi_test) ...
    + Ez_test .* cos(theta_test);

    E_theta_test = Ex_test .* cos(theta_test) .* cos(phi_test) ...
        + Ey_test .* cos(theta_test) .* sin(phi_test) ...
        - Ez_test .* sin(theta_test);

    E_phi_test = -Ex_test .* sin(phi_test) + Ey_test .* cos(phi_test);

    E_abs_test = sqrt(abs(Ex_test).^2 + abs(Ey_test).^2 + abs(Ez_test).^2);

    if all(abs(E_r_test - data_nf_rotated.E(:,1)) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz)/ Test Rotation / Variable E_r');
    else
    disp('❌ Es gibt Unterschiede / Test Rotation / Variable E_r');
    end

    if all(abs(E_theta_test - data_nf_rotated.E(:,2) ) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz)/ Test Rotation / Variable E_theta');
    else
    disp('❌ Es gibt Unterschiede / Test Rotation / Variable E_theta');
    end

    if all(abs(E_phi_test - data_nf_rotated.E(:,3)) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz)/ Test Rotation / Variable E_phi');
    else
    disp('❌ Es gibt Unterschiede / Test Rotation / Variable E_phi');
    end

    if all(abs(E_abs_test - data_nf_rotated.Eabs) < tol)
    disp('✅ Alle Werte sind gleich (mit Toleranz)/ Test Rotation / Variable E_abs');
    else
    disp('❌ Es gibt Unterschiede / Test Rotation / Variable E_abs');
    end
    %% Test Rotation Korrekt Ende
end
  
