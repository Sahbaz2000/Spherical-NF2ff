function [] = Fernfelddaten_Simulation(ant, freq,theta_cut,phi_cut)
% Diese Funktion visualisiert Fernfelddaten einer definierten Antenne (ant)
% bei einer bestimmten Frequenz (freq) mithilfe der MATLAB-Antenna Toolbox.
% Dabei werden verschiedene Darstellungen erzeugt, um die Richtcharakteristik
% der Antenne im Fernfeld grafisch auszuwerten.
%
% Es werden sowohl 3D-Diagramme als auch Phi- und Theta-Schnitte in verschiedenen
% Darstellungsformen (linear, logarithmisch, normiert) generiert. Zusätzlich
% ist eine Darstellung in Polarform möglich.
% 
% Input Arguments:
%
%       ant             Antennenobjekt (z. B. aus Antenna Toolbox)
%       freq            Frequenz in Hz
%       theta_cut       gewünschter Theta-Winkel in Grad (für Theta-Schnitt)
%       phi_cut         gewünschter Phi-Winkel in Grad (für Phi-Schnitt)
%
%
% Output Arguments:
%
%       Keine – die Funktion erzeugt mehrere Plots, liefert aber keine Rückgabewerte
%
% Es werden folgende Diagramme erzeugt:
%   - 3D-Fernfelddiagramm der Antenne
%   - Phi-Schnitt in verschiedenen Darstellungsarten (Gain, dB, normiert, polar)
%   - Theta-Schnitt in verschiedenen Darstellungsarten (Gain, dB, normiert, polar)

el= 90 - theta_cut; % es muss der Elevationswinkel übergeben werden
figure (Name='Matlab simulierte FF Daten, 3D Antennendiagramm')
pattern(ant,freq,Type="powerdb",normalize=true)
xAchse = -180:0.5:180;


%% Phi Cut
figure(Name='Matlab Simulierte FF Daten, Type Gain')
pattern(ant,freq,phi_cut,xAchse,CoordinateSystem="rectangular",Type="gain");
figure(Name='Matlab Simulierte FF Daten, in power Db')
pattern(ant,freq,phi_cut,xAchse,CoordinateSystem="rectangular",Type="powerdb");
% Normieren auf 0 dB Maximum
figure(Name='Matlab Simulierte FF Daten,in power Db Normiert')
pattern(ant,freq,phi_cut,xAchse,CoordinateSystem="rectangular",Type="powerdb",Normalize=true);
% Polar form
figure(Name='Matlab Simulierte FF Daten, in Dbi Polarform')
patternElevation(ant, freq, phi_cut);


%% Theta Cut
figure(Name='Matlab Simulierte FF Daten, Type Gain')
pattern(ant,freq,xAchse,el,CoordinateSystem="rectangular",Type="Gain");
figure(Name='Matlab Simulierte FF Daten, in power Db')
pattern(ant,freq,xAchse,el,CoordinateSystem="rectangular",Type="powerdb");
% Normieren auf 0 dB Maximum
figure(Name='Matlab Simulierte FF Daten,in power Db Normiert')
pattern(ant,freq,xAchse,el,CoordinateSystem="rectangular",Type="powerdb",Normalize=true);
% Polar form
figure(Name='Matlab Simulierte FF Daten, in Dbi Polarform')
patternAzimuth(ant, freq, el);
end
