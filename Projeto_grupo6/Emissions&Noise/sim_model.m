function f = sim_model(x)
  f = [];

  %% Inputs, Design Choices and Design Variables
  % x(1) - lift to drag ratio (-)
  % x(2) - disk loading (lb/ft2)
  % x(3) - wing's aspect ratio (-)
  % x(4) - cruise speed (kts)
  % x(5) - 1 trip range (nm)
  % x(6) - payload (lb)
  % x(7) - Mach number at the tip of the blade (-)
  % x(8) - solidity (-)
  % x(9) - maximum thickness to chord ratio (-)
  % x(10) - hybridization factor in cruise (-)

  % Configuration type (LC - Lift+Cruise; TR - Tilt-rotor; TW - Tilt-wing)
  config = 'TR';
  % Configuration Fixed Inputs
  if strcmpi(config,'LC')
    eta_p = 0.8;	% Cruise propeller's efficiency
    fws = 0.28;		% Structural factor
    NR = 8;       % Number of rotors
  elseif strcmpi(config,'TR')
	  eta_p = 0.76;	% Cruise propeller's efficiency
    fws = 0.3;		% Structural factor
    NR = 12;      % Number of rotors
  else
    eta_p = 0.76;	% Cruise propeller's efficiency
    fws = 0.3;		% Structural factor
    NR = 8;       % Number of rotors
  end

  % Propulsive System type (aC - All-combustion; aE - All-electric; hybrid - Hybrid-electric)
  type = 'hybrid';
  % Propulsive System - Efficiencies
  if strcmpi(type,'aE')
    n_es = 0.908;		% Efficiency of the electric system
  elseif strcmpi(type,'aC')
    n_cs = 0.429;		% Efficiency of the combustion system
  else
    HC = x(10);			% Hybridization factor in cruise (0 - all-combustion; 1 - all-electric)
    n_cs = 0.429;		% Efficiency of the combustion system
    n_es = 0.862;		% Efficiency of the electric system
  end

  % Fuel Selection (Fuel's specific energy density and emission in production phase)
  fuel = 'HEFA_jatropha';
  if strcmpi(fuel,'FT_SPK_camelina')
    Efuel = 44.2; % Specific Energy Density (MJ/kg))
    Fuel_production = 30.7; % (gCO2eq/MJ)
  elseif strcmpi(fuel,'FT_SPK_jatropha')
    Efuel = 44.2; % Specific Energy Density (MJ/kg))
    Fuel_production = 37.9; % (gCO2eq/MJ)
  elseif strcmpi(fuel,'FT_SPK_microalgae')
    Efuel = 44.2; % Specific Energy Density (MJ/kg))
    Fuel_production = 39.7; % (gCO2eq/MJ)
  elseif strcmpi(fuel,'HEFA_palm_oil')
    Efuel = 43.7; % Specific Energy Density (MJ/kg))
    Fuel_production = 52.0; % (gCO2eq/MJ)
  elseif strcmpi(fuel,'HEFA_jatropha')
    Efuel = 43.7; % Specific Energy Density (MJ/kg))
    Fuel_production = 73.5; % (gCO2eq/MJ)
  elseif strcmpi(fuel,'ATJ_wheat_straw')
    Efuel  = 43.2; % Specific Energy Density (MJ/kg))
    Fuel_production = 30.4; % (gCO2eq/MJ)
  elseif strcmpi(fuel,'ATJ_wheat_grain')
    Efuel = 43.2; % Specific Energy Density (MJ/kg))
    Fuel_production = 70.3; % (gCO2eq/MJ)
  else
    % Jet Fuel A /Kerosene
    Efuel = 43.28; % Specific Energy Density (MJ/kg))
    Fuel_production = 87.5; % (gCO2eq/MJ)
  end

  % Batteries related inputs
  Ebat = 500;     % energy specific density (W.h/kg)
  Ebat_r = 0.3;   % batteries reserve (-)
  N_cycles = 500; % limit number of battery cycles (-)


  % Maximum allowable values
  Pmax = 670;           % Maximum installed power (hp)
  SPLmax =65;           % Maximum SPL value (dB)
  bmax = 50;            % Maximum wingspan (ft)

  % General Inputs
  g = 9.80665;          % Gravity acceleration (m/s2)

  % Mission related inputs
  Nm = 4;                     	% Number of mission
  Reserve = 20/60;            	% Reserve time (h)
  HoverTime = 120*Nm/3600;    	% Total hover time (h)
  Vc = x(4)*0.5144;             % Cruise speed (m/s)
  Vy = 500*.3048/60;          	% Climb speed in VTOL mode (m/s)
  Range = x(5)*1.852*Nm;        % Range (km)
  Payload = x(6)*0.45359237;  	% Payload (kg)
  Endurance = Range/(Vc*3.6);   % Endurance (h)

  % Air properties
  rho = 1.14549;        % Air Density (kg/m3) - SL, ISA + 20ºC conditions
  rho_cr = 0.984762;    % Air Density (kg/m3) - 5000 ft + 20ºC conditions
  a = 351.906;          % Speed of sound (m/s) - SL, ISA + 20ºC Sea Level Conditions

  % Electric Motors related inputs
  n_ES = 0.862;		% Electric system efficiency (-)
  PWe = 5;        % Power to weight ratio (kW/kg)
  PWesc = 20;     % Power to weight ratio (kW/kg)
  MI = 0.55;      % Electric motor integration factor (-)

  % Combustion Engine related inputs
  PWg = 1/(0.5*0.608277388);  % Power to weight ratio (kW/kg)
  n_CS = 0.429;               % Power train efficiency of the combustion system (-)

  % Rotors related items
  s = x(8);       % solidity (-)
  cd0 = 0.01;     % airfoil base drag (-)
  t_c_max = x(9); % maximum thickness (-)
  ki = 1.2;       % induced power coefficient (-)
  Mtip = x(7);    % Mach number at the tip of the blade (-)
  Vtip = Mtip*a;  % Speed at the blade's tip (m/s)

  % Mass related inputs
  fwo = 0.22;     % fraction of the empty weight (-)


  % Other Inputs
  L_D = x(1);             % cruise lift-to-drag ratio (-)
  DL = x(2)*4.88242764*g; % disk loading (N/m2)
  AR = x(3);              % wing's aspect-ratio (-)


  %% Iterative Block
  error = 1;
  iter = 0;
  maxIter = 50;   % maximum number of iterations (-)
  MTOM = 3000;    % initial maximum take-off mass (kg)
  MTOM0 = MTOM;
  while (error > 0.001 && iter <maxIter)

    %% Power Estimation
    Th = MTOM0*g;                   % Thrust in hover (N)
    r = sqrt(Th/(pi()*NR*DL));      % Rotor's radius (m) - Assumption: all rotors are the same
    % Hover (Out of Ground Effect)
    Thr = Th/NR;                    % Thrust in hover per rotor (N) - Assumption: all rotors contribute equally to the total thrust
    Phr = Thr*(ki*sqrt(DL/(2*rho)) + (rho*Vtip*Vtip*Vtip/DL)*(s*cd0/8)); % Hover power per rotor (W)
    Ph = NR*Phr;                    % Hover power (W)
    FoM = Thr*sqrt(DL/(2*rho))/Phr; % Figure of merit (-)
    % Vertical Climb
    Tclr = Thr;                     % Thrust in vertical climb per rotor (N)- Assumption: the same thrust as in hover
    Pclr = Tclr*(Vy - ki*Vy/2 + ki*0.5*sqrt(Vy*Vy + 2*DL/rho) + ((rho*Vtip*Vtip*Vtip)/DL)*(s*cd0/8)); % Power in vertical climb per rotor (W)
    Pcl = NR*Pclr;                  % Power in vertical climb (W)
    % Cruise
    Pc = (MTOM0*g*Vc/(L_D*eta_p));  % Power in cruise (W)

    if strcmpi(type,'aE')
      %% Mass Estimation
      % Propulsive System
      Mem = NR*((Pclr/1000)/PWe + (Pclr/1000)/PWesc)*(1+MI);  % mass of the electric system (kg) - motor, ESC and integration
      Mce = 0;                                                % mass of the combustion engine (kg)
      Mg = 0;                                                 % mass of the generator (kg)
      Mprop = Mem + Mce + Mg;
      % Structural
      Mstr = fws*MTOM0;                                       % structural mass (kg) - wing, rotors, empennage and fuselage
      % Empty mass
      Mempty = (Mstr + Mprop)/(1-fwo);                        % empty mass (kg)
      % Other masses
      Mow = fwo*Mempty;                                       % other masses (kg) - avionics, flight control, instrumentation, enviromental control, furnishing

      %% Energy estimation
      Mfuel = 0;                                                                          % fuel mass (kg)
      Mbat = ((Ph/n_ES)*HoverTime + (Pc/n_ES)*(Endurance + Reserve))/(Ebat*(1 - Ebat_r)); % battery mass (kg)
      Menergy = Mfuel + Mbat;                                                             % energy mass (kg)

      Energy_fuel = 0;                                                                    % fuel energy (W.h)
      Energy_bat = ((Ph/n_ES)*HoverTime + (Pc/n_ES)*(Endurance + Reserve))/(1 - Ebat_r);  % battery energy (W.h)

    elseif strcmpi(type,'aC')
      %% Mass Estimation
      % Propulsive System
      Mem = NR*((Pclr/1000)/PWe + (Pclr/1000)/PWesc)*(1+MI);  % mass of the electric system (kg) - motor, ESC and integration
      Mce = Pcl/(1000*PWg);                                   % mass of the combustion engine (kg)
      Mg = Pcl/(1000*PWe);                                    % mass of the generator (kg)
      Mprop = Mem + Mce + Mg;
      % Structural
      Mstr = fws*MTOM0;                                       % structural mass (kg) - wing, rotors, empennage and fuselage
      % Empty mass
      Mempty = (Mstr + Mprop)/(1-fwo);                        % empty mass (kg)
      % Other masses
      Mow = fwo*Mempty;                                       % other masses (kg) - avionics, flight control, instrumentation, enviromental control, furnishing

      %% Energy estimation
      Mfuel = ((Pc/n_CS)*(Endurance + Reserve) + (Ph/n_CS)*HoverTime)*1.05*(3.6/1000)/Efuel;  % fuel mass (kg)
      Mbat = 0;                                                                               % battery mass (kg)
      Menergy = Mfuel + Mbat;                                                                 % energy mass (kg)

      Energy_fuel = ((Pc/n_CS)*(Endurance + Reserve) + (Ph/n_CS)*HoverTime)*1.05;             % fuel energy (W.h)
      Energy_bat = 0;                                                                         % battery energy (W.h)

    else
      %% Mass Estimation
      % Propulsive System
      Mem = NR*((Pclr/1000)/PWe + (Pclr/1000)/PWesc)*(1+MI);  % mass of the electric system (kg) - motor, ESC and integration
      Mce = Pc/(1000*PWg);                                    % mass of the combustion engine (kg)
      Mg = Pc/(1000*PWe);                                     % mass of the generator (kg)
      Mprop = Mem + Mce + Mg;
      % Structural
      Mstr = fws*MTOM0;                                       % structural mass (kg) - wing, rotors, empennage and fuselage
      % Empty mass
      Mempty = (Mstr + Mprop)/(1-fwo);                        % empty mass (kg)
      % Other masses
      Mow = fwo*Mempty;                                       % other masses (kg) - avionics, flight control, instrumentation, enviromental control, furnishing

      %% Energy estimation
      Mfuel = (1-HC)*((Pc/n_CS)*(Endurance + Reserve)*1.05)*(3.6/1000)/Efuel;                 % fuel mass (kg)
      Mbat = ((Ph/n_ES)*HoverTime + HC*(Pc/n_ES)*(Endurance + Reserve))/(Ebat*(1 - Ebat_r));  % battery mass (kg)
      Menergy = Mfuel + Mbat;                                                                 % energy mass (kg)

      Energy_fuel = (1-HC)*((Pc/n_CS)*(Endurance + Reserve)*1.05);                            % fuel energy (W.h)
      Energy_bat = ((Ph/n_ES)*HoverTime + HC*(Pc/n_ES)*(Endurance + Reserve))/(1 - Ebat_r);   % battery energy (W.h)
    end

    %% MTOM new estimation
    MTOM = Mempty + Payload + Menergy;
    error = abs((MTOM - MTOM0)/MTOM);
    MTOM0 = MTOM;
    iter = iter + 1;
  end

  %% Wing area and span
  % Stall speed or range requirements to determine the wing area and thus the wingspan
  Vstall = 65;                                                      % stall speed (kts)
  CLmax = 1.2;                                                      % maximum lift coefficient (-)
  S_stall = MTOM*g/(0.5*1.225*(Vstall*.5144)*(Vstall*.5144)*CLmax); % wing area (m2)
  CL_max_range = pi()*AR*0.8/(2*L_D);                               % lift coefficient for maximum range (-)
  S_range = MTOM*g/(0.5*1.225*Vc*Vc*CL_max_range);                  % wing area (m2)
  S = max(S_stall,S_range);                                         % wing area (m2)
  b = sqrt(AR*S)/0.3048;                                            % wingspan (ft)


  %% CO2-eq Emissions
  Electric_mix = 294;												% CO2-eq emissions per kW.h of electric energy (gCO2eq/(kW.h))
  Fuel_emissions = 8.887/142.2;                                     % CO2-eq emissions per MJ of fuel (gCO2eq/MJ)
  Emissions_energy = ((Energy_bat/1000)*Electric_mix + (Energy_fuel*3.6/1000)*Fuel_emissions)*N_cycles/1000; % CO2-eq emissions due to energy consumption (kgCO2eq)
  Battery_production = 80.8*1000;									% CO2-eq emissions per kW.h of battery produced (gCO2eq/(kW.h))
  Emissions_production = ((Energy_bat/1000)*Battery_production + (Energy_fuel*3.6/1000)*Fuel_production*N_cycles)/1000;  % CO2-eq emissions due to energy production (kgCO2eq)
  Emissions_total = Emissions_energy + Emissions_production;        % Total value of CO2-eq emissions (kgCO2eq)

  Emissions_fuel = ((Energy_fuel*3.6/1000)*Fuel_emissions + (Energy_fuel*3.6/1000)*Fuel_production)*N_cycles/1000;  % CO2-eq emissions due to fuel energy consumption and production (kgCO2eq)
  Emissions_battery = ((Energy_bat/1000)*Electric_mix*N_cycles + (Energy_bat/1000)*Battery_production)/1000;        % CO2-eq emissions due to battery energy production and consumption (kgCO2eq)


  %% Noise estimation
  m = 1;                              % (-) Harmonic number
  NB = 3;                             % (-) Number of blades [2-5]
  c = pi()*r*s/NB;                    % (m) Blade's mean chord
  t = t_c_max*c;                      % (m) Blade's maximum thickness
  Omega = Vtip/r;                     % (rad/s)
  Torque = Phr*r/Vtip;                % (N.m)
  Cl = 6*Thr/(rho*Vtip*Vtip*NB*c*r);  % (-) mean lift coefficient
  alt = 500*0.3048;                   % (m) - vertical distance of the observer from the source
  pref = 2.0e-5;                      % (Pa) Reference pressure
  St = 0.28;                          % (-) Strouhal number
  nsteps = 50;
  for i=1:(nsteps+1)
    aa = 89 - (i-1)*80/nsteps;
    ydist(i) = alt/sin(aa*pi()/180);
    DS = ydist(i);
    tt(i) = (pi()/2) + aa*(pi()/180);
    theta = tt(i);

    %% Rotational Noise
    Re = 0.8*r;
    JmB = besselj(m*NB,(m*NB*Omega*Re*sin(theta)/a));
    pmL = (m*NB*Omega/(2*sqrt(2)*pi()*a*DS))*(Thr*cos(theta) - Torque*a/(Omega*Re*Re))*JmB;
    pmT = -(rho*((m*NB*Omega)^2)*NB/(3*sqrt(2)*pi()*DS))*c*t*Re*JmB;
    pd_rot = NR*((pmL*pmL + pmT*pmT)/(pref*pref));

    %% Vortex Noise
    AA = pi()*r*r; % Disk area
    % 1st Step - Find Peak Frequency and the Overall SPL
    alpha = Cl/(2*pi());
    hh = t*cos(alpha) + c*sin(alpha);
    fpeak = (0.7*Vtip)*St/hh;
    k2 = 1.206*10^-2;   % (s^3/ft^3) Calibrated for helicopter
    k2 = k2/(.3048^3);  % (s^3/m^3)
    SPLo = 20*log10(k2*(Vtip/(rho*DS))*sqrt(NR*Thr*(Thr/AA)/s));
    % 2nd Step - Calculate SPL for each of spectrum frequencies
    freq = [fpeak/2; fpeak; 2*fpeak; 4*fpeak; 8*fpeak; 16*fpeak];
    SPL_w = [7.92; 4.17; 8.33; 8.75; 12.92; 13.33];
    SPL_ctr = 0;
    for ii = 1:5
      fr1 = freq(ii)/fpeak;
      SPL1 = SPLo - SPL_w(ii);
      fr2 = freq(ii+1)/fpeak;
      SPL2 = SPLo - SPL_w(ii+1);
      C1 = (SPL2 - SPL1)/(log10(fr2) - log10(fr1));
      C2 = SPL2 - C1*log10(fr2);
      int_val = (10^(C2/10))*((fr2^(C1/10 + 1))/(C1/10 + 1)) - (10^(C2/10))*((fr1^(C1/10 + 1))/(C1/10 + 1));
      SPL_ctr = SPL_ctr + int_val;
    end
    pd_vortex = SPL_ctr;

    %% Total Noise
    SPL(i) = 10*log10(pd_rot + pd_vortex);
  end
  SPL_max = max(SPL);



  %% Test if the concept is feasible for the maximum allowable values
  penalty = 0;
  Pcl = Pcl/745.7;
  Ph = Ph/745.7;
  Pc = Pc/745.7;
  if (b > bmax || Pcl > Pmax || Ph > Pmax || Pc > Pmax || SPL_max > SPLmax)
    penalty = 1e9;
  else
    penalty = 0;
  end

  %% Outputs
  f = zeros(1,26);
  f(1) = (Energy_fuel + Energy_bat);  % (W.h)
  f(2) = MTOM;                        % (kg)
  f(3) = Range;                       % (km)
  f(4) = Payload;                     % (kg)
  f(5) = Mempty;                      % (kg)
  f(6) = Mstr;                        % (kg)
  f(7) = Mprop;                       % (kg)
  f(8) = Mow;                         % (kg)
  f(9) = Mbat;                        % (kg)
  f(10) = Mfuel;                      % (kg)
  f(11) = Ph*745.7;                   % (hp)
  f(12) = Pcl*745.7;                  % (hp)
  f(13) = Pc*745.7;                   % (hp)
  f(14) = r;                          % (m)
  f(15) = FoM;                        % (-)
  f(16) = S;                          % (m2)
  f(17) = b*0.3048;                   % ()
  f(18) = Energy_fuel;                % (W.h)
  f(19) = Energy_bat;                 % (W.h)
  f(20) = SPL_max;                    % (dB)
  f(21) = Emissions_total/1000;       % (tons CO2eq)
  f(22) = Emissions_energy/1000;      % (tons CO2eq)
  f(23) = Emissions_production/1000;  % (tons CO2eq)
  f(24) = Emissions_fuel/1000;        % (tons CO2eq)
  f(25) = Emissions_battery/1000;     % (tons CO2eq)
  f(26) = penalty;                    % (-)
end
