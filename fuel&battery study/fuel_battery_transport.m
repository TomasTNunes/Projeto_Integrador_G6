clear all
close all
clc
addpath('../aircraft-design-tool-main/')

% Get data from json and algorithm
global constants;
constants.g = 9.81; % m/s^2
data = load_project('transport.json');
data.mission = build_mission(data.mission);
data.vehicle = build_vehicle(data.mission, data.vehicle);
data.vehicle = aero_analysis(data.mission, data.vehicle);
[data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);

mission = data.mission;
vehicle = data.vehicle;
energy = data.energy;

vehicle.components{10, 1}.mass = 0;
vehicle.components{11, 1}.mass = 0;

m_fuel = 0;
m_battery = 0;

% Fuel or Batt for each segment
for j = 1 : length(mission.segments)
    if strcmp(mission.segments{j}.type, 'taxi') % Taxi segment
        [vehicle m_fuel_seg(j)] = taxi(mission.segments{j}, vehicle, energy);
        % programa so faz contas se air taxi fosse a fuel
        m_battery_seg(j) = nan;
        m_fuel_seg(j) = nan;
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;
         
    elseif strcmp(mission.segments{j}.type, 'hover') % Hover segment
        [vehicle mf_batt_hover] = hover(mission.segments{j}, vehicle, energy);
        m_fuel_seg(j) = nan;
        m_battery_seg(j) = mf_batt_hover * vehicle.mass;
        m_battery = m_battery + m_battery_seg(j);
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;

    elseif strcmp(mission.segments{j}.type, 'climb') % Climb segment
        [vehicle m_fuel_seg(j)] = climb(mission.segments{j}, vehicle, energy);
        % nao usa brake specific fuel consumption
        m_battery_seg(j) = nan;
        m_fuel = m_fuel + m_fuel_seg(j);
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;

    elseif strcmp(mission.segments{j}.type, 'vertical_climb') % Vertical climb segment
        [vehicle mf_batt_vclimb] = vertical_climb(mission.segments{j}, vehicle, energy);
        m_fuel_seg(j) = nan;
        m_battery_seg(j) = mf_batt_vclimb * vehicle.mass;
        m_battery = m_battery + m_battery_seg(j);
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;

    elseif strcmp(mission.segments{j}.type, 'cruise') % Cruise segment
        [vehicle m_fuel_seg(j)] = cruise(mission.segments{j}, vehicle, energy);
        m_battery_seg(j) = nan;
        m_fuel = m_fuel + m_fuel_seg(j);
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;

    elseif strcmp(mission.segments{j}.type, 'descent') % Descent segment
        m_fuel_seg(j) = nan;
        m_battery_seg(j) = nan;
        % não gasta nada supostamente
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;

    elseif strcmp(mission.segments{j}.type, 'vertical_descent') % Vertical descent segment
        [vehicle mf_batt_vdescent] = vertical_descent(mission.segments{j}, vehicle, energy);
        m_fuel_seg(j) = nan;
        m_battery_seg(j) = mf_batt_vdescent * vehicle.mass;
        m_battery = m_battery + m_battery_seg(j);
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;

    elseif strcmp(mission.segments{j}.type, 'transition') % transition
        m_fuel_seg(j) = nan;
        m_battery_seg(j) = nan;
        seg{j} = strcat(mission.segments{j}.type,sprintf('(%s)',mission.segments{j}.name));
        seg_time(j) = mission.segments{j}.time;
    end
end

% Table
seg{end+1} = 'Total';
seg_time(end+1) = sum(seg_time);
m_fuel_seg(end+1) = m_fuel;
m_battery_seg(end+1) = m_battery;
row_n = {'time (s)' 'mass fuel (kg)' 'mass battery (kg)'};
T = table(seg_time(:),m_fuel_seg(:),m_battery_seg(:),'VariableNames',row_n,'RowNames',seg)

% Determine Fuel Volume
fuel_dens = 808; %Kg/m^3
fuel_volume = m_fuel/fuel_dens;
fprintf('Volume for fuel = %f m^3 = %f dm^3 = %f cm^3 \n',fuel_volume,fuel_volume*1000,fuel_volume*1000000)
fprintf('Caso seja um cubo é de dimensão %f x %f x %f cm \n\n',(fuel_volume*1000000)^(1/3),(fuel_volume*1000000)^(1/3),(fuel_volume*1000000)^(1/3))

% Determine Battery Volume
spec_ener_batt = data.vehicle.components{10, 1}.specific_energy / 3600; %Wh/Kg
energy_dens_batt = 450; %Wh/L
batt_dens = spec_ener_batt / energy_dens_batt; %L/Kg
batt_volume = batt_dens * m_battery; %L
fprintf('Volume for battery = %f m^3 = %f dm^3 = %f cm^3 \n',batt_volume/1000,batt_volume,batt_volume*1000)
fprintf('Caso seja um cubo é de dimensão %f x %f x %f cm \n',(batt_volume*1000)^(1/3),(batt_volume*1000)^(1/3),(batt_volume*1000)^(1/3))


%%%%%%%%%%%%%%%%%%%%%%-------Functions-------%%%%%%%%%%%%%%%%%%%%%%%

function [vehicle m_fuel] = taxi(segment, vehicle, energy)
[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');

if is_type(source, 'energy.fuel')
    mf_fuel = 1 - 0.9725
    m_fuel = mf_fuel * vehicle.mass;
    vehicle.components{network_ids(source_id)}.mass = source.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - vehicle.components{network_ids(source_id)}.mass;
else
    mf_fuel = 0;
    m_fuel = mf_fuel * vehicle.mass;
end
end


function [vehicle mf_batt] = hover(segment, vehicle, energy)
global constants;

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');

if is_type(source, 'energy.fuel')
    errordlg('Hover not available for fuel energy sources'); % NOT AVAILABLE
    return;
elseif is_type(source, 'energy.electric')
    rotor = find_by_type(network, 'driver.rotor');

    dl = vehicle.mass * constants.g / rotor_area(rotor);
    pl = rotor.efficiency * sqrt(2 * segment.density / dl);
    mf_batt = segment.time * constants.g / source.specific_energy / network_efficiency(network) / pl;
    vehicle.components{network_ids(source_id)}.mass = source.mass + mf_batt * vehicle.mass;
end
end


function [vehicle m_fuel] = climb(segment, vehicle, energy)
[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');

if is_type(source, 'energy.fuel')
    mach = segment.velocity / segment.speed_sound(1);
    if mach < 1
        mf_fuel = 1 - (1 - 0.04 * mach);
    else
        mf_fuel = 1 - (0.96 - 0.03 * (mach - 1));
    end
    m_fuel = mf_fuel * vehicle.mass;
    vehicle.components{network_ids(source_id)}.mass = source.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - vehicle.components{network_ids(source_id)}.mass;
elseif is_type(source, 'energy.electric')
    errordlg('Climb not available for electric energy sources'); % NOT AVAILABLE
    return;
end
end


function [vehicle mf_batt] = vertical_climb(segment, vehicle, energy)
global constants;

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');

if is_type(source, 'energy.fuel')
    errordlg('Vertical climb not available for fuel energy sources'); % NOT AVAILABLE
    return;
elseif is_type(source, 'energy.electric')
    rotor = find_by_type(network, 'driver.rotor');

    altitude_range = segment.altitude(2) - segment.altitude(1);
    dl = vehicle.mass * constants.g / rotor_area(rotor);
    pl = 1 / (segment.velocity - rotor.induced_power_factor / 2 * segment.velocity + rotor.induced_power_factor / 2 * sqrt(segment.velocity^2 + 2 * dl / segment.density(1)) + segment.density(1) * rotor.tip_velocity^3 / dl * rotor.rotor_solidity * rotor.base_drag_coefficient / 8); % Power loading
    mf_batt = altitude_range * constants.g / source.specific_energy / network_efficiency(network) / pl / segment.velocity; % Mass fraction for this segment
    vehicle.components{network_ids(source_id)}.mass = source.mass + mf_batt * vehicle.mass;
end
end


function [vehicle m_fuel] = cruise(segment, vehicle, energy)
global constants;

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');

engine = find_by_type(network, 'engine');
ld = get_ld(vehicle, segment, engine);

if is_type(source, 'energy.fuel')
    if is_type(engine, 'engine.jet')
        mf_fuel = breguet(segment.range, segment.velocity, engine.specific_fuel_consumption, ld);
    elseif is_type(engine, 'engine.prop')
        prop = find_by_type(network, 'driver.rotor');
        equivalent_sfc = engine.brake_specific_fuel_consumption * segment.velocity / prop.efficiency;
        mf_fuel = breguet(segment.range, segment.velocity, equivalent_sfc, ld);
    end
    m_fuel = mf_fuel * vehicle.mass;
    vehicle.components{network_ids(source_id)}.mass = source.mass + mf_fuel * vehicle.mass;
    vehicle.mass = vehicle.mass - vehicle.components{network_ids(source_id)}.mass;
elseif is_type(source, 'energy.electric')
    errordlg('Fiz isto para cruise com energy fuel e nao enegry eletric'); % NOT AVAILABLE
end
end


function [vehicle mf_batt] = vertical_descent(segment, vehicle, energy)
global constants;

[network, network_ids] = find_network_components(vehicle, find_by_name(energy.networks, segment.energy_network));
[source, source_id] = find_by_type(network, 'energy');
rotor = find_by_type(network, 'driver.rotor');

altitude_range = abs(segment.altitude(2) - segment.altitude(1));

if is_type(source, 'energy.fuel')
    errordlg('Vertical descent not available for fuel energy sources'); % NOT AVAILABLE
    return;
elseif is_type(source, 'energy.electric')
    dl = vehicle.mass * constants.g / rotor_area(rotor);
    v_i = sqrt(dl / 2 / segment.density(2)); % Induced velocity in hover
    if segment.velocity / v_i <= -2 % If this condition is met, the vertical climb equation is used for descent, else, an empirical equation is employed
        pl = 1 / (segment.velocity - rotor.induced_power_factor / 2 * (segment.velocity + sqrt(segment.velocity^2 - 2 * dl / segment.density(2))) + segment.density(2) * rotor.tip_velocity^3 / dl * rotor.rotor_solidity * rotor.base_drag_coefficient / 8);
    else
        v_d = v_i * (rotor.induced_power_factor - 1.125 * segment.velocity / v_i - 1.372 * (segment.velocity / v_i)^2 - 1.718 * (segment.velocity / v_i)^3 - 0.655 * (segment.velocity / v_i)^4); % Induced velocity in descent according to an empirical relation (see lecture slides)
        pl = 1 / (segment.velocity + rotor.induced_power_factor * v_d + segment.density(2) * rotor.tip_velocity^3 / dl * rotor.rotor_solidity * rotor.base_drag_coefficient / 8);
    end

    if pl > 0
        mf_batt = altitude_range * constants.g / source.specific_energy / network_efficiency(network) / pl / abs(segment.velocity);
    else
        mf_batt = 0;
    end
    vehicle.components{network_ids(source_id)}.mass = source.mass + mf_batt * vehicle.mass;
end
end


function e = network_efficiency(network)
e = 1.0;
for i = 1 : length(network)
    e = e * network{i}.efficiency;
end
end


function ld = get_ld(vehicle, segment, engine)
ld_max = estimate_ld_max(vehicle);
if strcmp(segment.type, 'cruise')
    if is_type(engine, 'engine.jet')
        ld = 0.886 * ld_max;
    elseif is_type(engine, 'engine.prop')
        ld = ld_max;
    end
elseif strcmp(segment.type, 'hold')
    if is_type(engine, 'engine.jet')
        ld = ld_max;
    elseif is_type(engine, 'engine.prop')
        ld = 0.886 * ld_max;
    end
end
end


function ld = estimate_ld_max(vehicle)
c = find_by_type(vehicle.components, 'wing.main');
ld = c.aspect_ratio + 10;
end


function mf_fuel = breguet(range, velocity, sfc, ld)
mf_fuel = 1 - exp(-range * sfc / velocity / ld);
end
