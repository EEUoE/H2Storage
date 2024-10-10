clear; clc; clf; close all
mrstModule add co2lab solvent ad-core ad-props mrst-gui compositional diagnostics incomp 
nx = 15;
ny = 15;
nz = 15;
casename = 'eec';
dims = [nx, ny, nz];
pdims = [1000, 1000, 1000]*meter;
G = cartGrid(dims, pdims);
G = computeGeometry(G);
%% Generate permeability using logNormLayers
rng(1);
[K,L] = logNormLayers([nx, ny, nz], [50, 80, 50], 'sigma', 1); 
target_kh = 30* milli * darcy; % Target mean horizontal permeability
mean_kh = mean(K(:)) * milli * darcy;
scale_factor_kh = target_kh / mean_kh;
K1_kh = K * scale_factor_kh;
K1_kh_mapped = K1_kh(G.cells.indexMap);
K1_kv_mapped = 0.1* K1_kh(G.cells.indexMap);
rng(1); 
phi = 0.2 + 0.02*randn(G.cells.num, 1); % Random distribution (with Mean as 0.2 and SD of 0.02)
rock = makeRock(G, [K1_kh_mapped*milli*darcy, K1_kh_mapped*milli*darcy, K1_kv_mapped*milli*darcy], phi);

figure;
subplot(1, 2, 1);
RKh = convertTo(rock.perm(:,1), milli*darcy);
plotCellData(G, RKh); axis on; view(3);
colorbarHist(RKh, [min(RKh) max(RKh)],'West', 100); title('Horizontal Permeability'); shading faceted;

subplot(1, 2, 2);
RKv = convertTo(rock.perm(:,3), milli*darcy);
plotCellData(G, RKv); axis on; view(3);
colorbarHist(RKv, [min(RKv) max(RKv)],'West', 100); title('Vertical Permeability'); shading faceted;

RP = rock.poro;
figure; plotCellData(G, RP); axis on; view(3);
colorbarHist(RP, [min(RP) max(RP)],'West', 100); title('Porosity'); shading faceted;
%% Set up fluid properties
rho_water_stc = 1000*kilogram/meter^3;
rho_oil_stc = 800*kilogram/meter^3;
rho_gas_stc = 0.083*kilogram/meter^3;
rho_co2_stc = 1.977*kilogram/meter^3;
mu_water_stc = 1*centi*poise;
mu_oil_stc = 5*centi*poise;
mu_gas_stc = 0.009*centi*poise;
mu_co2_stc = 0.0148*centi*poise;
%% Corey parameters
S_wr = 0.15;
S_or = 0.2; 
S_gr = 0.05; 
S_cr = 0.2;
n_w = 2.0; 
n_o = 2.5;  
n_g = 2.5;  
%% Define contacts, set up model and initial state
gravity reset on
%% Define GOC and WOC
c = 0.001/barsa;
depthD = 100*meter;
GOC = 500; % Depth of Gas-Oil Contact

[compfluid, info] = getBenchmarkMixtureOG(casename);
p_ref = info.pressure; 
% Initialize state with GOC and WOC
regions.pressure = p_ref + c * (G.cells.centroids(:,3) - depthD);
regions.temp = info.temp;
regions.s = zeros(G.cells.num, 3);
for i = 1:G.cells.num
    z = G.cells.centroids(i, 3);
    if z <= GOC
        regions.s(i, :) = [0, 0.1, 0.9]; % Gas
    else
        regions.s(i, :) = [0, 0.9, 0.1]; % OIL
    end
end
flowfluid = initSimpleADIFluid('phases', 'WOG', ...
    'mu'  , [mu_water_stc, mu_oil_stc, mu_gas_stc], ...
    'rho' , [rho_water_stc, rho_oil_stc, rho_gas_stc], ...
    'n'   , [1.5, 1, 2]);
sOr_i = 0.30;
sOr_m = 0.10;
sGc_i = 0.15;
sGc_m = 0.05;
flowfluid.bG = @(p, vargin) exp((p - p_ref)*c);
model = NaturalVariablesCompositionalModel(G, rock, flowfluid, compfluid);
state0 = initCompositionalState(G, regions.pressure, regions.temp, regions.s, info.initial, model.EOSModel);
state0.s = regions.s;
%% Define relative permeability functions
flowfluid.krW = @(sw) (max(sw - S_wr, 0) / (1 - S_wr - S_or)).^n_w;
flowfluid.krO = @(so) (max(so - S_or, 0) / (1 - S_or - S_gr)).^n_o;
flowfluid.krG = @(sg) (max(sg - S_gr, 0) / (1 - S_gr - S_wr)).^n_g;
%% Define wells
inj_rate = 0.85 * mega * 1e3 / year; % Final injection rate (400% of initial)
W = [];
W = verticalWell(W, G, rock, 1, 1, 1:2, 'comp_i', [0, 0, 1], 'name', 'Inj', ...
    'Val', inj_rate, 'sign', 1, 'type', 'rate');
for i = 1:numel(W)
    if strcmp(W(i).name, 'Inj')
        W(i).components = info.injection;
    end
end
plotGrid(G)
plotWell(G, W,'height', 500, 'color', 'k')
view(30,50)
[i, j, k] = ind2sub(G.cartDims, 1:G.cells.num);
state0.flux = fluxside([], G, 'Top',    0);
state0.flux = fluxside([], G, 'Bottom', 0);
state0.flux = fluxside([], G, 'Front',  0);
state0.flux = fluxside([], G, 'Back',   0);
state0.flux = fluxside([], G, 'Left',   0);
state0.flux = fluxside([], G, 'Right',  0);
figure;
% plot kH
subplot(1, 2, 1);
plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
plotCellData(G, convertTo(rock.perm(:, 1), milli*darcy), j == round(G.cartDims(2)/2))
title('Horizontal Permeability (kh)');
colorbar; view(3)
% plot kv
subplot(1, 2, 2);
plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
plotCellData(G, convertTo(rock.perm(:, 3), milli*darcy), j == round(G.cartDims(2)/2))
title('Vertical Permeability (kv)');
colorbar; view(3)
%% Set up schedule and simulate the problem
time = 365*day*1;
nn = 5*day;
dt = repmat(nn, time./nn, 1);
tt = cumsum(dt);
schedule = simpleSchedule(dt, 'W', W);
schedule.step.control(tt < 0.30*year) = 1;
schedule.step.control(tt >= 0.30*year & tt < 0.65*year) = 2;
schedule.step.control(tt >= 0.65*year) = 3;
schedule.control = repmat(schedule.control, 3, 1);
schedule.control(1).W(1).val = 0.25*inj_rate;
schedule.control(2).W(1).val = 1*inj_rate;
schedule.control(3).W(1).val = 1.75*inj_rate;
nls = NonLinearSolver('useRelaxation', true);
%% Set up schedule and simulate the problem
[ws, state, rep] = simulateScheduleAD(state0, model, schedule,'nonlinearsolver', nls);
%% Plot the results
inspectFluidModel(model)
plotWellSols(ws, dt);
figure
plotToolbar(G, state);
c2 = colorbar;
colorbar, view(30,50)
plotWell(G, W, 'height', 500, 'color', 'k')
figure
plotGridVolumes(G, state{end}.components(:,1), 'basealpha', 2)
c4 = colorbar;
plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
plotWell(G, W, 'height', 500, 'color', 'k');
colorbar, view(30,50);
colormap parula
pause(.5)
%% extract the number of cells that meet a certain criterion
SH = state{end}.components(:, 1);
count = 0;
for i = 1:length(SH)
    if SH(i) < 0.9 && SH(i) > 0.100
        count = count + 1;
    end
end
disp(['numbers of cells=', num2str(count)]);
disp(['percentages of cells=', num2str(count/length(SH)*100), '%']);
%%
figure
plotCellData(G, state0.s(:,1)); colorbar; view (-50, 50)
%% Simulate H2 production
t_ref = 323.15*Kelvin;
maxP = 50*barsa;
Wp = [];
Wp = verticalWell(Wp, G, rock, 1, 1, 1:2, 'comp_i', [0, 0, 1], 'Name', 'Prod', 'Val', maxP, 'sign', -1, 'Type', 'bhp');
Wp(1).components = [1, 0, 0, 0, 0,  0];
modelp = GenericOverallCompositionModel(G, rock, flowfluid, compfluid, 'water', true);
state0p = initCompositionalState(G, p_ref, t_ref, state{end,1}.s, state{end,1}.components, modelp.EOSModel);
schedulep = simpleSchedule(dt, 'W', Wp);
schedulep.step.control(tt < 0.30*year) = 1;
schedulep.step.control(tt >= 0.30*year & tt < 0.65*year) = 2;
schedulep.step.control(tt >= 0.65*year) = 3;

schedulep.control = repmat(schedulep.control, 3, 1);
schedulep.control(1).W(1).val = 0.9*maxP;
schedulep.control(2).W(1).val = 0.8*maxP;
schedulep.control(3).W(1).val = 0.60*maxP;
nls = NonLinearSolver('useRelaxation', true);
[wsp, statep, repp] = simulateScheduleAD(state0p, modelp, schedulep, 'nonlinearsolver', nls);
%% Plot H2 production solution
plotWellSols(wsp, dt);
figure; plotToolbar(G, statep); c6 = colorbar; colorbar, view(30,50)
plotWell(G, Wp, 'height', 500, 'color', 'k')
%% Calculate and plot H2 recovery
figure
T_mrst = convertTo(cumsum(schedulep.step.val), day);
mrstG = convertTo(abs(getWellOutput(wsp, 'qGs', 'Prod')), meter^3/day);
mrstO = convertTo(abs(getWellOutput(wsp, 'qOs', 'Prod')), meter^3/day);
mrstW = convertTo(abs(getWellOutput(wsp, 'qWs', 'Prod')), meter^3/day);
mrstplotG = @(T_mrst, data) plot(T_mrst, data, '-', 'linewidth', 2, 'color', 'k');
mrstplotO = @(T_mrst, data) plot(T_mrst, data, '--r', 'linewidth', 2, 'color', 'y');
mrstplotW = @(T_mrst, data) plot(T_mrst, data, '-.', 'linewidth', 2, 'color', 'm');
title('Fluid production profile (m^3/d)') 
xlabel('Time (days)')
yyaxis left;
mrstplotG(T_mrst, mrstG);
ylabel('Gas production rate (m^3/d)') 
hold on
yyaxis right; 
mrstplotO(T_mrst, mrstO); 
ylabel('Oil production rate (m^3/d)')
hold on
yyaxis right
mrstplotW(T_mrst, mrstW); 
legend({'Gas','Oil','Water'},'Location','northeast','Orientation','vertical')
%Total fluxes of the phases from the well are shown below in kg/s
wspH = zeros(1, numel(state)); wspM = zeros(1, numel(state)); wspN = zeros(1, numel(state)); wspC = zeros(1, numel(state)); 
wspuD = zeros(1, numel(state)); wspdD = zeros(1, numel(state));
for i = 1:numel(state)
    wspH(i) = wsp{i}.Hydrogen;
    wspM(i) = wsp{i}.Methane;
    wspN(i) = wsp{i}.Nitrogen;
    wspC(i) = wsp{i}.Carbondioxide;
    wspuD(i) = wsp{i}.nUndecane;
    wspdD(i) = wsp{i}.nDodecane;
end
% Conversion to volumetric flow rate (m^3/d) using rho in kg/m^3 at surface conditions (20 oC and 1 bara)
rhoH = 0.0830; rhoM = 0.7170; rhoN = 1.2506; rhoC = 1.860; rhouD = 740; rhodD = 750; 
wspH1 = (convertTo(abs(wspH), kilogram/day))./rhoH;
wspM1 = (convertTo(abs(wspM), kilogram/day))./rhoM;
wspN1 = (convertTo(abs(wspN), kilogram/day))./rhoN;
wspC1 = (convertTo(abs(wspC), kilogram/day))./rhoC;
wspuD1 = (convertTo(abs(wspuD), kilogram/day))./rhouD;
wspdD1 = (convertTo(abs(wspdD), kilogram/day))./rhodD;

figure
plot(T_mrst, abs(wspH1), '-', 'linewidth', 2, 'color', 'k');
hold on
plot(T_mrst, abs(wspM1), '--', 'linewidth', 2, 'color', 'm');
hold on
plot(T_mrst, abs(wspN1), '-', 'linewidth', 2, 'color', 'r');
hold on
plot(T_mrst, abs(wspC1), '-', 'linewidth', 2, 'color', 'b');
hold on
plot(T_mrst, abs(wspuD1), '-', 'linewidth', 2, 'color', 'c');
hold on
plot(T_mrst, abs(wspdD1), '-', 'linewidth', 2, 'color', 'g');

legend({'Hydrogen','Methane','Nitrogen','Carbondioxide', 'unDecane', 'doDecane'},'Location','northeast','Orientation','vertical')
title('Species fluxes (m^3/d)'); xlabel('Time (days)'); ylabel('Total flux (m^3/d)')
H2FL= wspH1;
CH4FL=wspM1;
CO2FL=wspC1;
N2FL=wspN1;
% Plot volume fractions
wspHv = zeros(1, numel(state)); wspMv = zeros(1, numel(state)); wspNv = zeros(1, numel(state)); wspCv = zeros(1, numel(state));
wspuDv = zeros(1, numel(state)); wspdDv = zeros(1, numel(state)); 
for i = 1:numel(state)
    wspHv(i) = wspH1(i)/(wspH1(i) + wspM1(i) + wspN1(i) + wspC1(i) + wspuD1(i) + wspdD1(i));
    wspMv(i) = wspM1(i)/(wspH1(i) + wspM1(i) + wspN1(i) + wspC1(i) + wspuD1(i) + wspdD1(i));
    wspNv(i) = wspN1(i)/(wspH1(i) + wspM1(i) + wspN1(i) + wspC1(i) + wspuD1(i) + wspdD1(i));
    wspCv(i) = wspC1(i)/(wspH1(i) + wspM1(i) + wspN1(i) + wspC1(i) + wspuD1(i) + wspdD1(i));
    wspuDv(i) = wspuD1(i)/(wspH1(i) + wspM1(i) + wspN1(i) + wspC1(i) + wspuD1(i) + wspdD1(i));
    wspdDv(i) = wspdD1(i)/(wspH1(i) + wspM1(i) + wspN1(i) + wspC1(i) + wspuD1(i) + wspdD1(i));
   
end
figure
plot(T_mrst, abs(wspHv), '-', 'linewidth', 2, 'color', 'k');
hold on
plot(T_mrst, abs(wspMv), '--', 'linewidth', 2, 'color', 'm');
hold on
plot(T_mrst, abs(wspNv), '-', 'linewidth', 2, 'color', 'r');
hold on
plot(T_mrst, abs(wspCv), '-', 'linewidth', 2, 'color', 'b');
hold on
plot(T_mrst, abs(wspuDv), '-', 'linewidth', 2, 'color', 'c');
hold on
plot(T_mrst, abs(wspdDv), '-', 'linewidth', 2, 'color', 'g');

legend({'Hydrogen','Methane','Nitrogen','Carbondioxide', 'unDecane', 'doDecane'},'Location','northeast','Orientation','vertical')
title('Species volume fractions'); xlabel('Time (days)'); ylabel('Volume fraction')
H2P = wspHv*100;
%% Material balance
Total1 = abs(wspH1(1) + wspM1(1) + wspN1(1) + wspC1(1));
Total2 = abs(mrstG(1));
Diff1 = Total2/Total1; display(Diff1)
Total3 = abs(wspuD1(1) + wspdD1(1));
Total4 = abs(mrstO(1));
Diff2 = Total3/Total4; display(Diff2)

%% Plot recovery as a fxn of time
qH2Inj  = convertTo(abs(getWellOutput(ws, 'qGs', 'Inj')), meter^3/day); %only H2 is injected; so qGs = qH2
TotalH2Inj = sum(qH2Inj); %Total H2 in place before production
qH2Prod = wspH1;
CummH2Prod = cumsum(qH2Prod); %Cummulative sum of produced H2.
H2Recovery = (CummH2Prod/TotalH2Inj)*100;
figure; plot(T_mrst, H2Recovery, '-', 'linewidth', 2, 'color', 'b');
title('Hydrogen recovery'); xlabel('Time (days)'); ylabel('Recovery (%)')
zpp = zeros(73, 1);
for i=1:73
PP = statep{i}.pressure(:, 1);
mPP=mean(PP);
disp(mPP);
zpp(i)=mPP;
end