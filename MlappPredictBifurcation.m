function[BifurcationHeight,tA]=MlappPredictBifurcation(ShockFilePath,VacuumFilePath, NumHeaderLines, timeColumn, PressureColumn, PitchColumn, CatchColumn, Mis,Mrs, gamma, MW)
% "MlappPredictBifurcation" - Mitchell D. Hageman October 2024
% PURPOSE:
%   *Determine liklihood and extent of shock bifurcation in a shock tube
% PREREQUISITES
%   * Experimental Pressure Trace
%   * Experimental laser absorbtion trace with strong schlieren spikes
%   * MlappPredictBifurcationHelp.pdf is needed in the PWD (present working directory) if you want the help button to work.
% REFERENCES
% Ref 1: THE INTERACTION OF A REFLECTED SHOCK WAVE WITH THE BOUNDARY LAYER IN A SHOCK TUBE
%       Herman Mark (NACA TM 1418) 1958
% Ref 2: Measurement of reflected-shock bifurcation over a wide range of gas composition and pressure
%       E. L. Petersen · R. K. Hanson, Shock Waves (2006) 15:333–340DOI 10.1007/s00193-006-0032-3
% Ref 3: Influence of Reflected Shock and Boundary‐Layer Interaction on Shock‐Tube Flows
%       L. Davies; J. L. Wilson Phys. Fluids 12, I-37–I-43 (1969)
% INPUTS:
%   * VacuumFilePath - full file path and name of csv or excel file with your Vacuum sample data
%       -example:  'C:\Users\mitchell.hageman\Desktop\Data\20240923_001_VacuumData.csv'
%       -assumes that any dark (zero) signal offset correction has already been applied to the recorded voltages.
%   * ShockFilePath - full file path and name of csv or excel file with your sample data
%       -example:  'C:\Users\mitchell.hageman\Desktop\Data\20240923_001_ShockData.csv'
%       -assumes that any dark (zero) signal offset correction has already been applied to the recorded voltages.
%       -Assumes the sample data file and vacuum data file are organized identically.
%   * NumHeaderLines - number of header lines in csv file before data begins.
%       -example: My data has two header lines.  Row 1 is data labels, and Row 2 is the offset voltage. So my data  begins in row 3.
%   *timeColumn - column of csv where time trace is.  Mine is Column A so Timecolumn=1.
%   *PressureColumn - column of csv where detector trace is.  Mine is column B so PressureColumn=2.
%   *PitchColumn - column of csv where LAS reference (Pitch) detector trace is.  Mine is column C so PitchColumn=3.
%   *CatchColumn - column of csv where LAS signal (Catch) detector trace is.  Mine is column D so CatchColumn=4.
%   **Note: The following would likely come from a normal shock equation solver, such as the FROzen chemistry SHock solver (FROSH)**
%   *Mis [-] - Incident shock Mach number
%   *Mrs [-] - Reflected shock Mach number
%%   *gamma [-] - specific heat ratio of the test gas.  NEED TO CHECK WHETHER TO USE STATE 1 or 2 FOR THESE CALCS
%   *MW_Mix [kg/kmol] - Molecular weight of the test gas.
%
% OUTPUTS:
%   *PointData
%       -BifurcationHeight [m]                - Value - See Fig1 & Eq(5) Ref 1 (BifurcationHeight=l)
%       -Reflected Shock arrival time, tA [s] - Value - Determined from Schlieren spike
%   * PointLables - Data lables for PointData - Character Matrix
% VERSION NUMBER:
%   * 1.0: January 2025 - initial release, Mitchell D. Hageman

%% Read in Experimental Data
ShockData = readmatrix(ShockFilePath,'NumHeaderLines',NumHeaderLines);
time = ShockData(:,timeColumn); %[ms]
PressureVolts=(ShockData(:,PressureColumn)); %[V] - raw trace.  Record for troubleshooting.
SampleReferenceVoltage=ShockData(:,PitchColumn); %[V]
SampleSignalVoltage=ShockData(:,CatchColumn); %[V]

VacuumData= readmatrix(VacuumFilePath,'NumHeaderLines',NumHeaderLines);
VacuumReferenceVoltage=VacuumData(:,PitchColumn); %[v]
VacuumSignalVoltage=VacuumData(:,CatchColumn); %[v]

It=SampleSignalVoltage./SampleReferenceVoltage; %It -with CMR.(Assumes Dark signal offset is already applied to loaded data)
Io=VacuumSignalVoltage./VacuumReferenceVoltage;%Io -with CMR.(Assumes Dark signal offset is already applied to loaded data)
IoAvg=mean(Io);
IoStdev=std(Io);
absorbance = real(-log((It)./IoAvg));% With Common mode rejection (Assumes Dark signal offset is already applied to loaded data)

%% Create Interactive Figure and ID tA using laser absorbance
fig = uifigure('Name', 'Interactive Plot', 'Position', [500 300 800 600]); % Create the main GUI window (figure)

% Create a button that opens a help file
helpButton = uibutton(fig, 'push', 'Text', 'Help!', ...
    'Position', [600 10 150 25], ...
    'ButtonPushedFcn', @(btn, event) openHelpFile());

% Create axes for plotting
ax = uiaxes(fig, 'Position', [50 100 700 475]);
plot(ax,time,absorbance,'DisplayName',"absorbance [-]"); %ax=axes you created, time=x, Pressure_atm=y, DisplayName=legend command, "Pressure [atm]" = display name.
legend(ax, 'show');
ax.YLabel.String ="Absorbance [-]";
ax.XLabel.String = 'time [ms]';

tOLabel = uilabel(fig, 'Text', 'First Deflection of Laser Beam, t_O [ms]', 'Position', [20 75 250 25]); %Create Label
tOInput = uieditfield(fig, 'numeric', 'Position', [235 75 50 25], 'Value', 0); %Create input box
tALabel = uilabel(fig, 'Text', 'Reflected Shock Arrival Time, t_A [ms]', 'Position', [325 75 250 25]); %Create Label
tAInput = uieditfield(fig, 'numeric', 'Position', [530 75 50 25], 'Value', 0); %Create input box
proceedButton = uibutton(fig, 'push', 'Text', 'Proceed to Pressure Trace', 'Position', [200 10 250 25], ...
    'ButtonPushedFcn', @(btn, event) resumeCalculation());
% Pause execution here until the user presses the Proceed button
uiwait(fig); % Pauses the code execution
tA=tAInput.Value/1000; %[s]
tO=tOInput.Value/1000; %[s]
%% Repeat Interactive Figure to ID tI and State 2/5 boundaries using PressureVolts
cla(ax);
plot(ax,time,PressureVolts,'DisplayName',"Pressure Signal [V]"); %ax=axes you created, time=x, Pressure_atm=y, DisplayName=legend command, "Pressure [atm]" = display name.
legend(ax, 'show');
ax.YLabel.String ="Pressure Transducer Signal [V]";
ax.XLabel.String = 'time [ms]';

%State 2
State2StartLabel = uilabel(fig, 'Text', 'State 2 Start [ms]:', 'Position', [175 50 100 25]); %Create Label
State2StartInput = uieditfield(fig, 'numeric', 'Position', [275 50 50 25], 'Value', 0); %Create input box
State2EndLabel = uilabel(fig, 'Text', '.State 2 End [ms]:', 'Position', [330 50 100 25]); %Create Label
State2EndInput = uieditfield(fig, 'numeric', 'Position', [430 50 50 25], 'Value', 0); %Create input box
%State 5
State5StartLabel = uilabel(fig, 'Text', 'State 5 Start [ms]', 'Position', [480 50 100 25]); %Create Label
State5StartInput = uieditfield(fig, 'numeric', 'Position', [580 50 50 25], 'Value', 0); %Create input box
State5EndLabel = uilabel(fig, 'Text', 'State 5 End [ms]', 'Position', [630 50 100 25]); %Create Label
State5EndInput = uieditfield(fig, 'numeric', 'Position', [730 50 50 25], 'Value', 0);

proceedButton = uibutton(fig, 'push', 'Text', 'Calculate Bifurcation Liklihood and Height', 'Position', [200 10 250 25], ...
    'ButtonPushedFcn', @(btn, event) resumeCalculation());
% Pause execution here until the user presses the Proceed button
uiwait(fig); % Pauses the code execution

[~, State2StartIndex]=min(abs(time-State2StartInput.Value)); %[index]
[~, State2EndIndex]=min(abs(time-State2EndInput.Value)); %[index]
[~, State5StartIndex]=min(abs(time-State5StartInput.Value)); %[index]
[~, State5EndIndex]=min(abs(time-State5EndInput.Value)); %[index]
P5overP2Transducer=mean(PressureVolts(State5StartIndex:State5EndIndex))/mean(PressureVolts(State2StartIndex:State2EndIndex));
%% Step 1 (Ref 2): Obtain delta_ta0_tAO from the available laser schlieren measurements
try
    V_R=evalin('base','Urs'); %[m/s] If you can get Urs directly from FROSH output, do it.
catch
    C=sqrt((gamma*8314*300)/MW); %Speed of sound =sqrt([J/kmol-K]*[K]/[kg/kmol])= sqrt([J/kg])=sqrt([kg*m^2]/[kg*s^2])=[m/s] ASSUMES 300K
    V_R=Mrs/C; %[m/s]
end
delta_ta0=tA-tO; %[s] -- See Fig 2, Ref 2 --    t_A=shock arrival time (Determined from Schlieren spike)

%% Step 2: calculate x1 from Eq. (4):
x_1 = delta_ta0*V_R; %[m] -- See Fig 1, Eq(4), Ref 2 --

%% Step 3: determine the boundary-layer pressure ratio from Eqs. (1), (2), and/or (3)
M_3=sqrt((2*gamma*Mis^2-(gamma-1))/((gamma-1)*Mis^2+2));%[-] Mach number at State 3 (Ref 1, Eq II-23)
M_BL=(2*(gamma-1)*Mis^2+(3-gamma))/((gamma+1)*Mis);%[-] Mach number in the LAMINAR BOUNDARY LAYER (Eq(1), Ref2 AND Eq III-7, Ref1)
%Note: M_BL for TURBULENT BOUNDARY LAYER-->Mark (Ref1) suggests that the transition Reynolds number from laminar to turbulent BL is 1.47E6-1.56E6...
%      and these Re occurred for an M1=2.15 in his experiments. However, he doesn't provide an equation for M_BL in this turbulent boundary layer, and...
%      He calls his desccription of transition to turbulence "crude." So, for now, we'll assume all BL follow the M_BL for the laminar BL.
if M_BL<1 %If it's subsonic
    PBLoverP2=(1+((gamma-1)/2)*M_BL^2)^(gamma/(gamma-1));%[-] ---Eq(2), Ref 2 AND Eq III-15, Ref1 ----
else %If it's supersonic
    PBLoverP2=((((gamma+1)/2)*M_BL^2)^(gamma/(gamma-1)))*((((2*gamma/(gamma+1))*M_BL^2)-((gamma-1)/(gamma+1)))^(1/(1-gamma)));%[-] --Eq(3), Ref2 AND Eq III-15, Ref1--NOTE The Petersen formula has a typo!!! starts with (gamma-1), but should be (gamma+1)
end
P5overP2=((2*gamma/(gamma+1))*M_3^2)-((gamma-1)/(gamma+1)); %[-] --Eq III-14 Ref1--

%% Aside:  make "bifurcation liklihood" plot for this specific case.
% Define Mach range for theoretical lines
startMach = 1;
endMach = 7;
Machstep = 0.1;
% Pre-allocate a vectors
numValues = (endMach - startMach) / Machstep + 1;
MisTheory = zeros(1, numValues);
M3Theory = zeros(1, numValues);
M_BLTheory = zeros(1, numValues);
PBLoverP2SupersonicTheory = zeros(1, numValues);
PBLoverP2SubsonicTheory = zeros(1, numValues);
P5overP2Theory = zeros(1, numValues);
% Initialize index
index = 1;
% Loop through values and store them in the vector
for MisTheorystep = startMach:Machstep:endMach%For a conceptual sweep of Mis
    MisTheory(index)=MisTheorystep;
    M3Theory(index)=sqrt((2*gamma*MisTheory(index)^2-(gamma-1))/((gamma-1)*MisTheory(index)^2+2));% eq II-23, Ref1
    M_BLTheory(index)=(2*(gamma-1)*MisTheory(index)^2+(3-gamma))/((gamma+1)*MisTheory(index));% eq III-7 Ref1
    PBLoverP2SupersonicTheory(index)=((((gamma+1)/2)*M_BLTheory(index)^2)^(gamma/(gamma-1)))*((((2*gamma/(gamma+1))*M_BLTheory(index)^2)-((gamma-1)/(gamma+1)))^(1/(1-gamma))); %Corrected from eq. 3 Ref2
    PBLoverP2SubsonicTheory(index)=(1+((gamma-1)/2)*M_BLTheory(index)^2)^(gamma/(gamma-1));%[-] eq. 3 Ref2
    P5overP2Theory(index)=((2*gamma/(gamma+1))*M3Theory(index)^2)-((gamma-1)/(gamma+1));% eq III-14, Ref1
    index=index+1;
end
SupersonicBLindex = find(M_BLTheory>1); %Find the index where the boundary layer becomes supersonic

figure
plot(MisTheory(1:SupersonicBLindex), PBLoverP2SubsonicTheory(1:SupersonicBLindex),'b--', 'DisplayName', 'P_{BL}/P_{2} Subsonic BL'); %Only plot subsonic PBL/P2 UP TO where M_BL goes supersonic
hold on
plot(MisTheory(SupersonicBLindex:end), PBLoverP2SupersonicTheory(SupersonicBLindex:end),'b-', 'DisplayName', 'P_{BL}/P_{2} Supersonic BL'); %Only plot supersonic PBL/P2 AFTER M_BL goes supersonic
plot(MisTheory, P5overP2Theory, 'k-', 'DisplayName', 'P_{5}/P_{2}'); 
plot([Mis Mis],[PBLoverP2 P5overP2], 'DisplayName','P_{5}/P_{2} to P_{BL}/P_{2} (Mark Calc)')
plot(Mis,P5overP2Transducer, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName','P_{5}/P_{2} Transducer') %Plot P5/P2 determined from the transducer trace, as a reference for the Marks calculation
try %If you have a 3rd measure of P5/P2 already saved in the base workspace, plot that as well for reference. I use FROSH to calculate P5 and P2, so that's the name I gave it. 
    P5overP2FROSH=((evalin('base','P5_atm'))*101325)/(evalin('base','P2'));
    plot(Mis,P5overP2FROSH, 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'DisplayName','P_{5}/P_{2} FROSH') %
catch
    %Otherwise do nothing and move on.
end
xlabel('Incident Shock Mach # [-]')
ylabel('Normalized Pressure [P/P_{2}]')
legend('show');
hold off

%% Step 4) calculate θ1 from Eq.(7) using γ for that mixture;
%Trig functions in MATLAB use radians. See https://www.mathworks.com/matlabcentral/answers/792987-solving-equations-with-trigonometric-functions

%!!!Symbolic Notation - Depends how you use the fuction!!!.
%syms theta1Sym %<-------------------------------This is simpler if you're using the function on it's own.
evalin('base', 'syms theta1Sym');%<--------------These two lines are required if you're using as a subfunction.
theta1Sym = evalin('base', 'theta1Sym'); %<------BUT, these two lines do modify the variable in the base workspace, so more careful attention is required.

eqn1=(sin(theta1Sym))^2==(((gamma+1)*(PBLoverP2)+gamma-1)/(2*gamma*Mrs^2));
theta1sol=double(solve(eqn1,theta1Sym)); %Solve eqn1 for theta1, convert to doubles. There will be a positive and negative result
theta1=theta1sol(theta1sol>=0); %[rad] Keep only the positive result
if PBLoverP2<=P5overP2
    %% Step 5 -calculate l from Eq(5)
    BifurcationHeight = x_1*tan(theta1); %[m] --See Fig1 & Eq(5) Ref 1 (BifurcationHeight=l) --
    BifurcationHeight_Petersenfit=7.5*(Mis^1.07)/(gamma^-2.66*MW^-0.37); % [mm] Ref2 eq 10 - specific to 5cm Diam CST, D=5mm PCB113A PXD, 20mm from endwall
    % Note: "The thickness of the boundary layer...was neglected... this extra height should be added"... (Ref 2)
    %...But I haven't found any way to approximate the thickness of the boundary layer, so I'm also ignoring it. 
    %% STEP5(a) (MDH) Calculate remaining Fig1, Ref2 geometry
    AngleOAB=pi()/2; %[rad]=90deg --MDH Assumption. This is notional.
    AngleCOA =theta1; %[rad] --Ref3, Fig 1

    %!!!Symbolic Notation - Depends how you use the fuction!!!.
    %SymAngleCOB %<--------------------------------------This is simpler if you're using the function on it's own.
    evalin('base', 'syms SymAngleCOB');%<----------------These two lines are required if you're using as a subfunction.
    SymAngleCOB = evalin('base', 'SymAngleCOB'); %<------BUT, these two lines do modify the variable in base workspace, so more careful attention is required.

    MBLsin2COA=(M_BL^2)*(sin(AngleCOA))^2;
    eqn2=tan(AngleCOA-SymAngleCOB)/tan(AngleCOA)==((gamma-1)*MBLsin2COA)/((gamma+1)*MBLsin2COA);% Ref3, top left of page I-39
    AngleCOB=double(solve(eqn2,SymAngleCOB)); %[rad]
    AngleBOA=AngleCOA-AngleCOB; %[rad] - Trigonometry
    LengthAO=BifurcationHeight/sin(AngleCOA);%[m] (Trigonometry) --See Fig1 Ref3 and Fig1 Ref2--
    x_1_check=(BifurcationHeight)*cos(AngleCOA);%[m] (Trig) --See Fig1, Ref2--should match x_1 from delta_ta0*V_R;
    LengthOB=LengthAO/cos(AngleBOA); %[m] - Trig, but relies on MDH assumption that AngleOAB=90deg.
    LengthAB=LengthAO*tan(AngleBOA); %[rad] - Trig, but relies on MDH assumption that AngleOAB=90deg.
    x_2=cos(AngleCOB)*LengthOB;%[m] - Trig, --See Eq(6), Ref2 for a different way, but angles aren't defined there--
    LengthCB=x_2*tan(AngleCOB);
    i=0;
    for semicircleAngle=0:pi()/200:pi()/2
        i=i+1;
        semicircleYaxis(i)=cos(semicircleAngle)*LengthCB;
        semicircleXaxis(i)=x_2+sin(semicircleAngle)*LengthCB;
    end
    %%Plotlines
    LineOA=[0,0;x_1,BifurcationHeight];
    LineOC=[0,0;x_2,x_2*tan(AngleCOB)];
    LineAB=[x_1,BifurcationHeight;x_2,x_2*tan(AngleCOB)];
try
    figure
    plot(LineOA(:,1),LineOA(:,2), 'DisplayName', 'Line OA')
    hold on
    plot(LineOC(:,1),LineOC(:,2), 'DisplayName', 'Line OC')
    plot(LineAB(:,1),LineAB(:,2), 'DisplayName', 'Line AB')
    plot(semicircleXaxis,semicircleYaxis, 'DisplayName', 'Notional Line CD')
    xlabel('Distance from Bifurcation Foot [m]')
    ylabel('Height from Wall [m]')
    if max(semicircleXaxis)>BifurcationHeight
        xlim([0 max(semicircleXaxis)]);
        ylim([0 max(semicircleXaxis)]);
    else
        xlim([0 BifurcationHeight]);
        ylim([0 BifurcationHeight]);
    end
    legend('show');
    hold off
catch
    message=sprintf('Check your command window for t_A, t_0, and Bifurcation Height. \n The Bifurcation Height Plot failed.\n This typically happens when t_0 > t_A, which causes bifurcation height to be incorrect.');
    msgbox(message, 'Information');
end
else
    BifurcationHeight='NA';
    message=sprintf('PBL/P2 > P5/P2.\n No Bifurcation is anticipated. \n Calculated Bifurcation height = NA');
    msgbox(message, 'Information');
end

%% Callback Functions
%Callback function to open the help file
    function openHelpFile()
        helpFile = 'MlappLASTraceHelp.pdf';  % Path to the help file
        if isfile(helpFile) % Use MATLAB's open function to open the file
            open(helpFile);
        else
            % Display a message if the file doesn't exist
            uialert(fig, 'Help file not found!', 'Error');
        end
    end

%Callback function for "Proceed..." button
    function resumeCalculation()
        % Resume program execution
        uiresume(fig);
    end

hold(ax, 'off'); %allow next plug-in to move on to the next figure.
end