function [] = plotNF2FF_thetaCut(data_nf2ff,theta_cut,normalized,logarithmic)
% Diese Funktion erzeugt einen Theta-Cut (bei einem gewünschten theta in Grad (theta_cut)), wobei 
% alle zugehörigen Phi-Werte auf der x-Achse und die entsprechenden Eabs-Werte
% auf der y-Achse dargestellt werden. Die Darstellung von Eabs kann wahlweise 
% logarithmisch in dB (logarithmic = true) oder linear in V/m (logarithmic = false) 
% erfolgen. Zusätzlich kann der Wert entweder auf 0 dB normiert (normalize = true) 
% oder ungenormt (normalize = false) dargestellt werden.
% Wird hauptsächlich dafür genutzt die transformierten Daten zu plotten.
%
% Input Arguments:
%
%       data_nf2ff             Tabelle  mit mindestens 3 Spalten nötig:
%                              phi, theta, Eabs
%
%       theta_cut              gewünschter Winkel Theta(°) bei dem ein Schnitt
%                              gemacht werden soll
%
%       normalized            true für normierung auf 0 Db / false ohne normierung
%
%       logarithmic           true für logarithmisch in Db / false in V/m                    
%
%
% Output Arguments:
%
%       Keine                Es wird ein Plot enstehen ohne Rückgabe wert
%
deg = theta_cut; % für den Titel benötigt
theta_cut = deg2rad(theta_cut); % von Grad in Radian

% Zeilen suchen wo die gesuchten thetas sind
index = find(min(abs(data_nf2ff.theta - theta_cut))==abs(data_nf2ff.theta - theta_cut));

% Alle zugehörigen Phis und Eabs zuordnen
phi_thetaCut_nf2ff = data_nf2ff.phi(index);
eabs_thetaCut_nf2ff = data_nf2ff.Eabs(index);

% Neu ordnen anstannt von 0-360 geht es von 180-360 dann 0-179 
idx1 = (phi_thetaCut_nf2ff >= pi);    % ab 180 bis 360
idx2 = (phi_thetaCut_nf2ff < pi);     % von 0 bis <180
%Neu zusammensetzen: erst idx1, dann idx2
phi_thetaCut_nf2ff = [phi_thetaCut_nf2ff(idx1)-2*pi;  phi_thetaCut_nf2ff(idx2)];
eabs_thetaCut_nf2ff  = [eabs_thetaCut_nf2ff(idx1);    eabs_thetaCut_nf2ff(idx2)];

% Normierung auf 0 Db falls normalized true
if normalized == true
   eabs_thetaCut_nf2ff = eabs_thetaCut_nf2ff ./ max(eabs_thetaCut_nf2ff);
end

% Logarithmische Darstellung falls logarithmic true
if logarithmic == true
    eabs_thetaCut_nf2ff =  20 * log10(eabs_thetaCut_nf2ff);
end

% Plot erzeugen 
figure (Name='Transformierte FF Daten, Theta Schnitt')
plot(rad2deg(phi_thetaCut_nf2ff), eabs_thetaCut_nf2ff, 'LineWidth', 1.5);
xlabel('Phi [°]');
if logarithmic == true
ylabel('|E| [dB]');
else 
ylabel('|E| [V/m]');  
end
title(['Theta = ', num2str(deg),' Schnitt'])
end