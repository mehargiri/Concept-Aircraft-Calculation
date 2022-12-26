% Aero 3002 Design Course Calculations

% Redoing Sensitivity analysis and performance matching after finidng new
% L/D max value using drag polar analysis (L/D max = 12.815)
%% Wetted Area Calculations
clear; clc;
% Main Wing
Sexp_wing = 1911; % mm2
Swet_wing = Sexp_wing*(1.977 + (0.52*0.12)); % mm2

% Tail
Sexp_tail = 1068; % mm2
Swet_tail = Sexp_tail*(1.977 + (0.52*0.12)); % mm2

% Canard
Sexp_can = 435; % mm2
Swet_can = Sexp_can*(1.977 + (0.52*0.12)); % mm2

% Fuselage
Atop = 1994; Aside = 1993; % mm2
Swet_fuse = 3.4*((Atop + Aside)/2); % mm2

Swet = Swet_wing + Swet_tail + Swet_can + Swet_fuse; %mm2
fprintf('\nThe wetted area of the concept aircraft is %.2f mm^2\n', Swet);

% Wetted Aspect Ratio (Awetted)

% The wing span (b) of the concept aircraft is 150 mm.
b = 150;
Awetted = b^2/Swet;

fprintf('\nThe wetted aspect ratio of the concept aircraft is %.2f\n',Awetted);

% Using Fig 3.5 in Aircraft Design book, the calculated Awetted value gives
% L/D max of 12.185

%% Initial Sizing

clc; clear;
% Design Specifications

Range = 200; % nmi 
LDmax = 12.2815; Cruise_vel = 250; % knots;
Passenger = 6; Pass_weight = 80; % kg
Payload_weight = 15; % kg
End = 10; % min

% Range Conversion to ft

Rangefts = Range*6076.12; % ft

% Cruise Velocity conversion to ft/s

Cruise_velfts = Cruise_vel*1.68781; % ft/s

% Total Passenger and Paylod Weight in lb

Wtotal = (Passenger*Pass_weight*2.205) + (Passenger*Payload_weight*2.205); % lb;

% Endurance conversion to seconds

Endsec = End * 60; % sec

% Using cruise velocity of 250 knots and Fig 3.3 in Aircraft Design
% Textbook, the SFC (Specific Fuel Consumption) of the concept aircraft is
% 0.7 lb/(lb*hr) during crusie and 0.6 lb/(lb*hr) during loiter.

SFC_cruise = 0.7/3600; % sec
SFC_loit = 0.6/3600; % sec

% Mission Segment Weight Fractions

% Historical Weight Fractions for different mission segment

% Warmup and takeoff = 0.97   Climb  = 0.985    Landing = 0.995
% Info taken form Table 3.2 in Aircraft Design Textbook

W1W0 = 0.97; % takeoff
W2W1 = 0.985; %climb

RC = Rangefts * SFC_cruise;
VLDVal = Cruise_velfts * LDmax;

W3W2 = exp(-RC/VLDVal);

EC = Endsec * SFC_loit;

W4W3 = exp(-EC/LDmax);

W5W4 = W3W2; % Second cruise and loiter is same as first one
W6W5 = W4W3;

W7W6 = 0.995; % Landing

W7W0 = W1W0 * W2W1 * W3W2 * W4W3 * W5W4 * W6W5 * W7W6;

WfW0 = 1.06 * (1- W7W0);

W0gVal = []; WeW0Val = []; WeVal  = []; W0CVal = [];

for W0guess = [10000 5800 5700 5771 5772]
    W0gVal = [W0gVal, W0guess];

    WeW0 = -6*10^-6 * W0guess + 0.6619; % equation derived from research aircrafts
    WeW0Val = [WeW0Val, WeW0];

    W0Cal = Wtotal/(1 - (WfW0) - (WeW0));
    W0CVal = [W0CVal, W0Cal];

    We = WeW0 * W0Cal;
    WeVal = [WeVal, We];
end
 
fprintf('\n W0, guess    We/W0        We         W0, calculated\n')
formatspec = '%7.0f %12.4f %12.3f %15.3f \n';
fprintf(formatspec,[W0gVal; WeW0Val; WeVal; W0CVal]);

W0 = W0CVal(5);
Wf = W0 * WfW0;
We = W0 * ((-6*10^-6 * W0) + 0.6619);
fprintf('\nThe max take-off weight of the concept aircraft is %.3f lb or %.3f kg\n', W0, W0/2.205);

fprintf('\nThe max empty weight is %.3f lb or %.3f kg\n', We, We/2.205);
fprintf('\nThe max fuel weight is %.3f lb or %.3f kg\n', Wf, Wf/2.205);

fprintf('\nThe empty weight fraction is %.7f\n', We/W0);

%% Sensitivity Studies

% values of 25% below and above the max value is taken into consideration
% for these sensitivity studies

% Max take off weight is calculated for each changes in value

%% L/D changes
clear; clc;
for LDVal = 25.625
    Range = 200; % nmi 
    Cruise_vel = 250; % knots;
    Passenger = 6; Pass_weight = 80; % kg
    Payload_weight = 15; % kg
    End = 10; % min

    % Range Conversion to ft

    Rangefts = Range*6076.12; % ft

    % Cruise Velocity conversion to ft/s

    Cruise_velfts = Cruise_vel*1.68781; % ft/s

    % Total Passenger and Paylod Weight in lb

    Wtotal = (Passenger*Pass_weight*2.205) + (Passenger*Payload_weight*2.205); % lb;

    % Endurance conversion to seconds

    Endsec = End * 60; % sec

    % Using cruise velocity of 250 knots and Fig 3.3 in Aircraft Design
    % Textbook, the SFC (Specific Fuel Consumption) of the concept aircraft is
    % 0.7 lb/(lb*hr) during crusie and 0.6 lb/(lb*hr) during loiter.

    SFC_cruise = 0.7/3600; % sec
    SFC_loit = 0.6/3600; % sec

    % Mission Segment Weight Fractions

    % Historical Weight Fractions for different mission segment

    % Warmup and takeoff = 0.97   Climb  = 0.985    Landing = 0.995
    % Info taken form Table 3.2 in Aircraft Design Textbook

    W1W0 = 0.97; % takeoff
    W2W1 = 0.985; %climb

    RC = Rangefts * SFC_cruise;
    VLDVal = Cruise_velfts * LDVal;

    W3W2 = exp(-RC/VLDVal);

    EC = Endsec * SFC_loit;

    W4W3 = exp(-EC/LDVal);

    W5W4 = W3W2; % Second cruise and loiter is same as first one
    W6W5 = W4W3;

    W7W6 = 0.995; % Landing

    W7W0 = W1W0 * W2W1 * W3W2 * W4W3 * W5W4 * W6W5 * W7W6;

    WfW0 = 1.06 * (1- W7W0);

    W0gVal = []; WeW0Val = []; WeVal = []; W0CVal = [];
    for W0guess = 4764
        W0gVal = [W0gVal, W0guess];

        WeW0 = -6*10^-6 * W0guess + 0.6619; % equation derived from research aircrafts
        WeW0Val = [WeW0Val, WeW0];

        W0Cal = Wtotal/(1 - (WfW0) - (WeW0));
        W0CVal = [W0CVal, W0Cal];

        We = WeW0 * W0Cal;
        WeVal = [WeVal, We];
    end
    
    fprintf('\n W0, guess    We/W0        We         W0, calculated\n')
    formatspec = '%7.0f %12.4f %12.3f %15.3f \n';
    fprintf(formatspec,[W0gVal; WeW0Val; WeVal; W0CVal]);
end

%% Range changes
clc; clear;
% Design Specifications

Range = 250; % nmi 
LDmax = 20.5; Cruise_vel = 250; % knots;
Passenger = 6; Pass_weight = 80; % kg
Payload_weight = 15; % kg
End = 10; % min

% Range Conversion to ft

Rangefts = Range*6076.12; % ft

% Cruise Velocity conversion to ft/s

Cruise_velfts = Cruise_vel*1.68781; % ft/s

% Total Passenger and Paylod Weight in lb

Wtotal = (Passenger*Pass_weight*2.205) + (Passenger*Payload_weight*2.205); % lb;

% Endurance conversion to seconds

Endsec = End * 60; % sec

% Using cruise velocity of 250 knots and Fig 3.3 in Aircraft Design
% Textbook, the SFC (Specific Fuel Consumption) of the concept aircraft is
% 0.7 lb/(lb*hr) during crusie and 0.6 lb/(lb*hr) during loiter.

SFC_cruise = 0.7/3600; % sec
SFC_loit = 0.6/3600; % sec

% Mission Segment Weight Fractions

% Historical Weight Fractions for different mission segment

% Warmup and takeoff = 0.97   Climb  = 0.985    Landing = 0.995
% Info taken form Table 3.2 in Aircraft Design Textbook

W1W0 = 0.97; % takeoff
W2W1 = 0.985; %climb

RC = Rangefts * SFC_cruise;
VLDVal = Cruise_velfts * LDmax;

W3W2 = exp(-RC/VLDVal);

EC = Endsec * SFC_loit;

W4W3 = exp(-EC/LDmax);

W5W4 = W3W2; % Second cruise and loiter is same as first one
W6W5 = W4W3;

W7W6 = 0.995; % Landing

W7W0 = W1W0 * W2W1 * W3W2 * W4W3 * W5W4 * W6W5 * W7W6;

WfW0 = 1.06 * (1- W7W0);

W0gVal = []; WeW0Val = []; WeVal = []; W0CVal = [];

for W0guess = 5206
    W0gVal = [W0gVal, W0guess];

    WeW0 = -6*10^-6 * W0guess + 0.6619; % equation derived from research aircrafts
    WeW0Val = [WeW0Val, WeW0];

    W0Cal = Wtotal/(1 - (WfW0) - (WeW0));
    W0CVal = [W0CVal, W0Cal];

    We = WeW0 * W0Cal;
    WeVal = [WeVal, We];
end
 
fprintf('\n W0, guess    We/W0        We         W0, calculated\n')
formatspec = '%7.0f %12.4f %12.3f %15.3f \n';
fprintf(formatspec,[W0gVal; WeW0Val; WeVal; W0CVal]);

%% SFC Cruise changes
clc; clear;
% Design Specifications

Range = 200; % nmi 
LDmax = 20.5; Cruise_vel = 250; % knots;
Passenger = 6; Pass_weight = 80; % kg
Payload_weight = 15; % kg
End = 10; % min

% Range Conversion to ft

Rangefts = Range*6076.12; % ft

% Cruise Velocity conversion to ft/s

Cruise_velfts = Cruise_vel*1.68781; % ft/s

% Total Passenger and Paylod Weight in lb

Wtotal = (Passenger*Pass_weight*2.205) + (Passenger*Payload_weight*2.205); % lb;

% Endurance conversion to seconds

Endsec = End * 60; % sec

% Using cruise velocity of 250 knots and Fig 3.3 in Aircraft Design
% Textbook, the SFC (Specific Fuel Consumption) of the concept aircraft is
% 0.7 lb/(lb*hr) during crusie and 0.6 lb/(lb*hr) during loiter.

for SFC = 0.875
    SFC_cruise = SFC/3600; % sec
    SFC_loit = 0.6/3600; % sec

    % Mission Segment Weight Fractions

    % Historical Weight Fractions for different mission segment

    % Warmup and takeoff = 0.97   Climb  = 0.985    Landing = 0.995
    % Info taken form Table 3.2 in Aircraft Design Textbook

    W1W0 = 0.97; % takeoff
    W2W1 = 0.985; %climb

    RC = Rangefts * SFC_cruise;
    VLDVal = Cruise_velfts * LDmax;

    W3W2 = exp(-RC/VLDVal);

    EC = Endsec * SFC_loit;

    W4W3 = exp(-EC/LDmax);

    W5W4 = W3W2; % Second cruise and loiter is same as first one
    W6W5 = W4W3;

    W7W6 = 0.995; % Landing

    W7W0 = W1W0 * W2W1 * W3W2 * W4W3 * W5W4 * W6W5 * W7W6;

    WfW0 = 1.06 * (1- W7W0);

    W0gVal = []; WeW0Val = []; WeVal = []; W0CVal = [];

    for W0guess = 5206
        W0gVal = [W0gVal, W0guess];

        WeW0 = -6*10^-6 * W0guess + 0.6619; % equation derived from research aircrafts
        WeW0Val = [WeW0Val, WeW0];

        W0Cal = Wtotal/(1 - (WfW0) - (WeW0));
        W0CVal = [W0CVal, W0Cal];

        We = WeW0 * W0Cal;
        WeVal = [WeVal, We];
    end

    fprintf('\n W0, guess    We/W0        We         W0, calculated\n')
    formatspec = '%7.0f %12.4f %12.3f %15.3f \n';
    fprintf(formatspec,[W0gVal; WeW0Val; WeVal; W0CVal]);
end

%% We/Wo Changes
clc; clear;
% Design Specifications

Range = 200; % nmi 
LDmax = 20.5; Cruise_vel = 250; % knots;
Passenger = 6; Pass_weight = 80; % kg
Payload_weight = 15; % kg
End = 10; % min

% Range Conversion to ft

Rangefts = Range*6076.12; % ft

% Cruise Velocity conversion to ft/s

Cruise_velfts = Cruise_vel*1.68781; % ft/s

% Total Passenger and Paylod Weight in lb

Wtotal = (Passenger*Pass_weight*2.205) + (Passenger*Payload_weight*2.205); % lb;

% Endurance conversion to seconds

Endsec = End * 60; % sec

% Using cruise velocity of 250 knots and Fig 3.3 in Aircraft Design
% Textbook, the SFC (Specific Fuel Consumption) of the concept aircraft is
% 0.7 lb/(lb*hr) during crusie and 0.6 lb/(lb*hr) during loiter.

for SFC = 0.7
    SFC_cruise = SFC/3600; % sec
    SFC_loit = 0.6/3600; % sec

    % Mission Segment Weight Fractions

    % Historical Weight Fractions for different mission segment

    % Warmup and takeoff = 0.97   Climb  = 0.985    Landing = 0.995
    % Info taken form Table 3.2 in Aircraft Design Textbook

    W1W0 = 0.97; % takeoff
    W2W1 = 0.985; %climb

    RC = Rangefts * SFC_cruise;
    VLDVal = Cruise_velfts * LDmax;

    W3W2 = exp(-RC/VLDVal);

    EC = Endsec * SFC_loit;

    W4W3 = exp(-EC/LDmax);

    W5W4 = W3W2; % Second cruise and loiter is same as first one
    W6W5 = W4W3;

    W7W6 = 0.995; % Landing

    W7W0 = W1W0 * W2W1 * W3W2 * W4W3 * W5W4 * W6W5 * W7W6;

    WfW0 = 1.06 * (1- W7W0);

    W0gVal = []; WeW0Val = []; WeVal = []; W0CVal = [];

    for W0guess = 8000
        W0gVal = [W0gVal, W0guess];

        WeW0 = 0.790087; % equation derived from research aircrafts
        WeW0Val = [WeW0Val, WeW0];

        W0Cal = Wtotal/(1 - (WfW0) - (WeW0));
        W0CVal = [W0CVal, W0Cal];

        We = WeW0 * W0Cal;
        WeVal = [WeVal, We];
    end

    fprintf('\n W0, guess    We/W0        We         W0, calculated\n')
    formatspec = '%7.0f %12.4f %12.3f %15.3f \n';
    fprintf(formatspec,[W0gVal; WeW0Val; WeVal; W0CVal]);
end

%% Payload Changes

clc; clear;
% Design Specifications

Range = 200; % nmi 
LDmax = 20.5; Cruise_vel = 250; % knots;
Passenger = 6; Pass_weight = 80; % kg
Payload_weight = 18.75; % kg
End = 10; % min

% Range Conversion to ft

Rangefts = Range*6076.12; % ft

% Cruise Velocity conversion to ft/s

Cruise_velfts = Cruise_vel*1.68781; % ft/s

% Total Passenger and Paylod Weight in lb

Wtotal = (Passenger*Pass_weight*2.205) + (Passenger*Payload_weight*2.205); % lb;

% Endurance conversion to seconds

Endsec = End * 60; % sec

% Using cruise velocity of 250 knots and Fig 3.3 in Aircraft Design
% Textbook, the SFC (Specific Fuel Consumption) of the concept aircraft is
% 0.7 lb/(lb*hr) during crusie and 0.6 lb/(lb*hr) during loiter.

SFC_cruise = 0.7/3600; % sec
SFC_loit = 0.6/3600; % sec

% Mission Segment Weight Fractions

% Historical Weight Fractions for different mission segment

% Warmup and takeoff = 0.97   Climb  = 0.985    Landing = 0.995
% Info taken form Table 3.2 in Aircraft Design Textbook

W1W0 = 0.97; % takeoff
W2W1 = 0.985; %climb

RC = Rangefts * SFC_cruise;
VLDVal = Cruise_velfts * LDmax;

W3W2 = exp(-RC/VLDVal);

EC = Endsec * SFC_loit;

W4W3 = exp(-EC/LDmax);

W5W4 = W3W2; % Second cruise and loiter is same as first one
W6W5 = W4W3;

W7W6 = 0.995; % Landing

W7W0 = W1W0 * W2W1 * W3W2 * W4W3 * W5W4 * W6W5 * W7W6;

WfW0 = 1.06 * (1- W7W0);

W0gVal = []; WeW0Val = []; WeVal  = []; W0CVal = [];

for W0guess = 5146
    W0gVal = [W0gVal, W0guess];

    WeW0 = -6*10^-6 * W0guess + 0.6619; % equation derived from research aircrafts
    WeW0Val = [WeW0Val, WeW0];

    W0Cal = Wtotal/(1 - (WfW0) - (WeW0));
    W0CVal = [W0CVal, W0Cal];

    We = WeW0 * W0Cal;
    WeVal = [WeVal, We];
end
 
fprintf('\n W0, guess    We/W0        We         W0, calculated\n')
formatspec = '%7.0f %12.4f %12.3f %15.3f \n';
fprintf(formatspec,[W0gVal; WeW0Val; WeVal; W0CVal]);

%% Performance Matching of various flight condition

% Stall
clear; clc;
row = 1.23; % kg/m3
Vstall = 61; % knots
Clmax_stall = 1.80645;

WS_stall = 0.5*row*Clmax_stall*(Vstall/1.944)^2; % N/m2

% Take-off

TOP = 66.66; % lbf/ft2;
den_ratio = 1;
Cl_TO = Clmax_stall/1.1;

WS_take = @(TW) TOP*47.88*den_ratio*Cl_TO*TW;

% Landing
Sl = 750; Sa = 183;%m
Cl_Land = Clmax_stall; % Cl of landing is same as Cl max at stall
WS_landing = (1/5)*(Sl-Sa)*den_ratio*Cl_Land*9.81;

% Cruise
 % Necessary relationship between T/W (Thrust to Weight ratio) and L/D
 % (Lift to drag ratio) during cruise
 
 % (T/W)cruise = 1/(L/D) cruise
 
 % At cruise, L/D = 1/(qCdo/(W/S) + (W/S*q*pi*A*e))
 
 CD0 = 0.015; % Zero lift drag coefficent for a jet aircraft
 row_25Kft = 0.550196; % kg/m3
 V = 250/1.944; % m/s
 A_wing = 9.9075; % Aspect ratio of wing
 
 q = 0.5*row_25Kft*V^2; % kg/ms2 dynamic pressure at 7620 m 
 
 e = 1.78*(1-0.045*A_wing^0.68) - 0.64; % Oswald efficiency factor for a straight- wing aircraft
 
 TW_cruise = @(WS) ((q*CD0)/WS) + (WS/(q*pi()*A_wing*e));
 
 % Climb
 
 % Thrust to weight equation for climb can be derived through climb
 % equations of motion
 
 % Rate of climb is related to thrust to weight and wing loading using the
 % following equation
 
 % Rc = ((2*WS)/(row*CL_LDmax))^0.5 * (TW - (1/LDmax))
 
 % The rate of climb equation can be isolated for TW (Thrust to weight)
 
 Rc = 0.508; % m/s rate of climb for service
 row_sea = 1.24; % kg/m3 density at sea level
 LDmax = 12.2815;
 
 CL_LDmax = (CD0 * pi()*A_wing*e)^0.5;
 
 TW_climb = @(WS) Rc/(((2*WS)/(row_sea*CL_LDmax))^0.5) + 1/(LDmax);
 
 % Plotting takeoff performance plot
 TakeVal = []; CrusVal = []; ClimbVal = [];
 for TW = linspace(0,0.5)
     Take = WS_take(TW);
     TakeVal = [TakeVal,Take];
 end
 
 for WS = linspace(100,2500)
     Cruise = TW_cruise(WS);
     CrusVal = [CrusVal,Cruise];
 end
 
 for WS = linspace(0,2500)
     Climb = TW_climb(WS);
     ClimbVal = [ClimbVal, Climb];
 end
 figure(1);
 TW = linspace(0,0.6);
 WS_cruise = linspace(100,2500);
 WS_climb = linspace(0,2500);
 plot(TakeVal,TW,'r');
 hold on; grid on;
 plot(WS_cruise,CrusVal,'b');
 plot(WS_climb,ClimbVal,'m');
 plot([WS_stall,WS_stall],[0,0.7],'k');
 plot([WS_landing,WS_landing],[0,0.7],'c');
 title('Performance Matching Analysis for various conditions');
 xlabel('Wo/Sref (N/m^2)');
 ylabel('To/Wo');
 
 % Intersection Points (Design Points)
plot(557.475,0.12763025,'k*');
plot(WS_stall,0.250435992,'g*');
plot(WS_stall, 0.3734, 'r*');
legend('Take-off','Cruise','Climb','Stall','Landing','Design Point 1','Design Point 2', 'Design Point 3');
hold off;

fprintf('\nDesign Point 3 is selected\n');
fprintf('\nThe selected T/W is 0.3734 and the selected W/S is %.4f N/m^2\n',WS_stall);
Thrust = 0.3734*2254.773*9.81; % N
MTOW = 2617.624; % kg Max take-off weight
Sref = (MTOW*9.81)/WS_stall;

fprintf('The selected Sref for the aircraft is %.5f m^2\n',Sref);
fprintf('\nThe selected thrust for the aircraft is %.4f N for two engines\n',Thrust);
fprintf('\nThe selected thrust for the aircraft is %.4f N for one engine\n',Thrust/2);

% Wingspan
b = (A_wing*Sref)^0.5; %m
fprintf('\nThe wing span of the aircraft is %.4f m\n',b);
%% Drag Polar Analysis
clear; clc;
A = 9.90753;
e = 0.759032;

K = 1/(pi()*A*e);
Sref = 23.47513; % m2
Sref_Imp = 19.69207*10.764; % ft2

Cd0 = 8.3/Sref_Imp;

Drag_Coeff = @(Cl) Cd0 + K*Cl.^2;

s = 1/(2*Cd0^(1/2)*K^(1/2));

CdVal = []; a = 1; b = [];
for Cl = 0:0.01:2.5
    Cd = Drag_Coeff(Cl);
    CdVal = [CdVal, Cd];
    b = [b,a];
    a = a+1;
    
end

figure(2);
Cl = 0:0.01:2.5;
Clmax = 2.5;
plot(CdVal,Cl); hold on;
plot([0, Clmax/s],[0,Clmax]);
plot(CdVal(97), Cl(97), 'r*');
title('Drag Polar Analysis');
xlabel('Cd');
ylabel('Cl');
legend('Aircraft Drag','Tangent Line','L/D max');
hold off;
grid on;


% syms Cd Cd0 K s
% solve(Cd == Cd0 + K*(s*Cd)^2,Cd)
% 
% solve(1 - 4*Cd0*K*s^2 == 0, s)
Cl_coeff = 0:0.01:2.5;
y1 = Drag_Coeff(Cl_coeff);
slope = Cl_coeff./y1;

fprintf('\n Iterations      Cl         Cd        Slope\n')
formatspec = '%7.0f %13.4f %10.4f  %10.4f\n';
fprintf(formatspec,[b; Cl_coeff; y1; slope]);

%% V-N Diagram
clear; clc;
% Need to assess the different loads experienced by the aircraft at
% different velocities
% Vdive Vstall different Vgust

Vcruise = 250; % knots
MTOW = 2255.159; %kg
Clmax = 1.80645;
WS = 1093.8763; %N/m^2
row = 1.23; %kg/m3 density at sea level

%Dive Speed
Vdive = 1.4*Vcruise;  % (units: knots) design dive speed from CAR 523.335c

% Stall
% Load Equation for stall : n = (row*Vstall^2*Cl*Sref)/(2*Wo);
 
%Stall Speed
Vstall = ((1*2*WS)/(row*Clmax))^0.5; % m/s
Vstall_knots = Vstall*1.944; % knots

%Positive limit on manoeuvring load factor

n_pos = 2.1 + 24000/(MTOW+10000);

%positive value of n must be 3.8 or lower

% Negative limit on maoeuvring load factor

n_neg = 0.4*n_pos;

% plotting V-N diagram (with maneuver loads)
Vval_pos = []; Vval_neg = [];
Speed = @(n) ((n*2*WS)/(row*Clmax))^0.5;
for n = linspace(0,n_pos,100)
    velo = Speed(n);
    Vval_pos = [Vval_pos,velo];
end

for n = linspace(0,n_neg,100)
    velo = Speed(n);
    Vval_neg = [Vval_neg,velo];
end
Vval_pos_knots = Vval_pos.*1.944;
Vval_neg_knots = Vval_neg.*1.944;
n = linspace(0,n_pos,100);
n2 = linspace(0,n_neg,100);

figure(1);
plot(Vval_pos_knots,n,'b');
axis([0 400 -3 5]);
title('V-n diagram (maneuver)');
xlabel('Knots Equivalent Airspeed');
ylabel('n');
hold on;
plot(Vval_neg_knots,n2.*-1,'b');
plot([Vval_pos_knots(100), Vdive], [n(100), n(100)],'b');
plot([Vval_neg_knots(100), Vcruise], [n2(100)*-1 , n2(100)*-1],'b');
plot([Vdive, Vdive], [n(100) ,0] ,'b');
plot([Vcruise, Vdive], [n2(100)*-1, 0],'b');
p1 = plot([Vstall_knots, Vstall_knots], [-2, 1], 'r');
plot([Vcruise, Vcruise], [-2, -1.7], 'g');
p2 = plot([Vcruise, Vcruise], [-1.5, 0], 'g');
p3 = plot([Vdive, Vdive], [-0.1, -1], 'k');
plot([0, 400], [0, 0],'k');
legend([p1, p2, p3],{'Vstall','Vcruise','Vdive'});
grid on;
hold off;
    
%Gust Loads
% aircraft can experience load due to strong wind gusts

% Change in load factor due to wind gust must be calculated

row_cruise = 0.5502; % kg/m3 density at 25000ft
Clalpha = 4.6; % Slope of the lift coefficient vs angle of attack graph for NACA 4415 (selected airfoil)
taper = 0.4; % taper ratio
b = 14.1542; %m
g = 9.81 ; %m/s^2 gravity
S = 20.221; %m2

croot = 2*S/(b*(1+taper));
Ude_cruise_imp = 50 ; % feets per second
Ude_cruise = Ude_cruise_imp/3.281; % m/s  gust speed at cruise given in CARs 523.333c

Ude_dive_imp = 25; % feets per second
Ude_dive = Ude_dive_imp/3.281; % m/s gust speed at dive


cbar = (2/3)*croot*((1+taper+taper^2)/(1+taper));
mass_ratio = 2*WS/(row_cruise*g*cbar*Clalpha); % mass ratio


K = (0.88*mass_ratio)/(5.3 + mass_ratio); % gust alleviation factor

U_cruise = K*Ude_cruise;
U_dive = K*Ude_dive;


deltan_cruise = @(V) (row_cruise*U_cruise*V*Clalpha)/(2*WS);
deltan_dive = @(V) (row_cruise*U_dive*V*Clalpha)/(2*WS);


% Veocity from 0 to Vcruise
deln = [];
for V = linspace(0,Vcruise,100)
    load = deltan_cruise(V);
    deln = [deln,load];
end

figure(2);
v3 = linspace(0,Vcruise,100);
hold on; grid on;
plot(v3,deln+1,'--k');
plot(v3,(deln.*-1)+1,'--k');

% maneuver v-n diagram
plot(Vval_pos_knots,n,'--k');
plot(Vval_neg_knots,n2.*-1,'--k');
plot([Vval_pos_knots(100), Vdive], [n_pos, n_pos],'--k');
plot([Vval_neg_knots(100), Vcruise], [n_neg*-1 , n_neg*-1],'--k');
plot([Vdive, Vdive], [n_pos ,0] ,'--k');
plot([Vcruise, Vdive], [n_neg*-1, 0],'--k');


% Velocity from 0 to Vdive
deln2 = [];
for V = linspace(0,Vdive,100)
    load = deltan_dive(V);
    deln2 = [deln2,load];
end

v4 = linspace(0,Vdive,100);
plot(v4,deln2+1,'--k');
plot(v4,(deln2.*-1)+1,'--k');
plot([v3(100), v4(100)], [deln(100)+1, deln2(100)+1],'k');
plot([v3(100), v4(100)], [(deln(100)*-1)+1, (deln2(100)*-1)+1],'k');
plot([Vdive, Vdive], [deln2(100)+1, (deln2(100)*-1)+1],'k');
plot([215.111, Vcruise], [n_pos, deln(100)+1], 'k');
plot([215.111, 92.939], [n_pos, 2.321367], 'k');
plot([Vcruise, 40.029], [(deln(100)*-1)+1, 0.4308], 'k');
pl1 = plot(Vcruise,0,'r*');
pl2  = plot(Vdive, 0, 'ko');
plot([0, 400], [0, 0],'b');
legend([pl1 pl2], {'Vcruise', 'Vdive'});
axis([0 400 -3 5]);
title('V-n diagram (gust)');
xlabel('Knots Equivalent Airspeed');
ylabel('n');
hold off;


% Combined V-N diagram
figure(3);
hold on; grid on;
plot(Vval_pos_knots,n,'b');
plot(Vval_neg_knots,n2.*-1,'b');
plot([Vval_pos_knots(100), 215.111], [n_pos, n_pos],'b');
plot([215.111, Vcruise], [n_pos, deln(100)+1], 'b');
plot([Vcruise, 296.519], [deln(100)+1, n_pos],'b');
plot([296.519, Vdive], [n_pos, n_pos], 'b');
plot([Vdive, Vdive], [n_pos, (deln2(100)*-1)+1], 'b');
plot([Vdive, Vcruise], [(deln2(100)*-1)+1, deln(100)*-1], 'b');
plot([Vcruise, 184.514], [deln(100)*-1, n_neg*-1],'b');
plot([184.514, Vval_neg_knots(100)], [n_neg*-1, n_neg*-1], 'b');
pl1 = plot(Vcruise,0,'r*');
pl2  = plot(Vdive, 0, 'ko');
plot([0, 400], [0, 0],'k');
legend([pl1 pl2], {'Vcruise', 'Vdive'});
axis([0 400 -4 5]);
title('Combined V-n diagram');
xlabel('Knots Equivalent Airspeed');
ylabel('n');
hold off;