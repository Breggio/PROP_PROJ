function [Rb, Peff, tburn1, tburn2] = calulate_Rb(p)

web = 0.03;  %propellant thickness, 3 cm

maxitems = max( size( p ));

%Pressure for the definition of the action time
[Pmax, imax] = max( p);
Pact = 0.05 .* Pmax; % controllare dalla registrazione questo valore

i = imax;

%For the ignition
while ( p( i ) > Pact )
  i = i - 1;
  if i == 0  
    break;
  end;
  
end
tact1 = i;

i = imax;

%For the extinction
while ( p( i ) > Pact )
  i = i + 1;

  if i == maxitems;  
    break;
  end;
  
end
tact2 = i;

% Reference pressure for burning time
Pref = trapz( p( tact1:tact2) ) / (tact2 - tact1) / 2.0 ;

i = imax;
% For the ignition
while ( p( i ) > Pref )
  i = i - 1;
  if i == 0  
    break;
  end;
  
end
tburn1 = i;

i = imax;
% For the extinction
while ( p( i ) > Pref )
  i = i + 1;

  if i == maxitems;  
    break;
  end;
  
end
tburn2 = i;

% Burning time, s
Deltat = ( tburn2 - tburn1) / 1000;

% Effective burning pressure
Peff = trapz( p( tburn1:tburn2) ) / (tburn2 - tburn1) ;

% Burning rate, m/s
Rb = web / Deltat;

end