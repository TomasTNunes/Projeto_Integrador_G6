clear all
close all
clc
addpath('../aircraft-design-tool-main/')

% Get data from json
global constants;
constants.g = 9.81; % m/s^2
data = load_project('rescue.json');
data.mission = build_mission(data.mission);
data.vehicle = build_vehicle(data.mission, data.vehicle);
data.vehicle = aero_analysis(data.mission, data.vehicle);
[data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);

mission = data.mission;
vehicle = data.vehicle;
energy = data.energy;

% Fuel characteristics (Jet Fuel A-1)
M_fuel = vehicle.components{11, 1}.mass * (1 + vehicle.components{11, 1}.reserve); % fuel mass w/ reserve (kg)
fuel_spec_ene_dens = 43.2; % Fuel Specific Energy Density (MJ/kg)
fuel_dens = 808; % Fuel density (kg/m^3)
fuel_vol_ene_dens = fuel_spec_ene_dens * fuel_dens / 1000; % Fuel Volumetric Energy Density (MJ/L)
fuel_production = 87.5; % (gCO2eq/MJ)
fuel_emissions = 0.0635; % CO2-eq emissions per MJ of fuel (gCO2eq/MJ)
Energy_fuel = fuel_spec_ene_dens * M_fuel * 1000 / 3.6; % fuel energy (W.h)

% Battery characteristics (Li-Ion Battery +/-)
M_bat = vehicle.components{10, 1}.mass * (1 + vehicle.components{10, 1}.reserve); % Battery mass w/ reserve (kg)
bat_spec_ene_dens = vehicle.components{10, 1}.specific_energy / 1000000; % Battery Specific Energy Density (MJ/kg)
bat_vol_ene_dens = 450; % Fuel Volumetric Energy Density (W.h/L)
battery_production = 147.7 * 1000 ; % CO2-eq emissions per kW.h of battery produced (gCO2eq/(kW.h))
Electric_mix = 324,7; % CO2-eq emissions per kW.h of electric energy (gCO2eq/(kW.h)) [Portugal Eletric Mix]
Energy_bat = bat_spec_ene_dens * M_bat * 1000 / 3.6; % Battery energy (W.h)
N_cycles = 500;

% CO2-eq Emissions p/ mission
Emissions_bat_prod = Energy_bat * battery_production / 1000000 / N_cycles % (kgCO2eq)
Emissions_bat_energy = Energy_bat * Electric_mix / 1000000 % (kgCO2eq)
Emissions_bat = Emissions_bat_prod + Emissions_bat_energy % CO2-eq emissions due to battery energy production and consumption (kgCO2eq)

Emissions_fuel_prod = Energy_fuel * 3.6 / 1000 * fuel_production / 1000 % (kgCO2eq)
Emissions_fuel_energy = Energy_fuel * 3.6 / 1000 * fuel_emissions / 1000 % (kgCO2eq)
Emissions_fuel = Emissions_fuel_prod + Emissions_fuel_energy % CO2-eq emissions due to fuel energy consumption and production (kgCO2eq)

Emissions_production = Emissions_bat_prod + Emissions_fuel_prod % CO2-eq emissions due to energy production (kgCO2eq)
Emissions_energy = Emissions_bat_energy + Emissions_fuel_energy % CO2-eq emissions due to energy consumption (kgCO2eq)
Emissions_total = Emissions_production + Emissions_energy % Total value of CO2-eq emissions (kgCO2eq)

Emissons_use = Emissions_fuel + Emissions_bat_energy % CO2-eq emissions due to use phase(KgCO2eq)
