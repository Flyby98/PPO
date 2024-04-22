%% load modules
mrstModule add yby;
mrstModule add shale;
mrstModule add compositional;
mrstModule add hfm;
mrstModule add ad-core;
mrstModule add ad-props;
mrstModule add mrst-gui;

start_string = y_datetime('second');


%% Grid and rock
G = cartGrid(the_input.G.cartdim, the_input.G.physdim);
G = computeGeometry(G);
G.rock = makeShaleRock(G, the_input.rock.perm, the_input.rock.poro);

if strcmp(the_input.rock.poro_mode, 'file')
    load_poro = load(the_input.rock.poro_file);
    G.rock.poro = load_poro.poro;
end

if strcmp(the_input.rock.perm_mode, 'file')
    load_perm = load(the_input.rock.perm_file);
    G.rock.perm = load_perm.perm;
end

G_plot = G;
if checkPlot(the_input, the_input.plot.grid)
    figure;
    plotGrid(G, 'facealpha', 0.2, 'edgealpha', 0.5);view(3);
    title('Grid');
end

if checkPlot(the_input, the_input.plot.rock)
    figure;
    % there is something wrong if G.cell.num change from 28000 to 28714, but the data if right, just plot it independently.
    plotCellData(G_plot, G_plot.rock.perm);view(3);shading flat;
    slice(reshape(G.rock.perm, G.cartDims(2), G.cartDims(1), G.cartDims(3)), [G.cartDims(1)/2, G.cartDims(1)], [G.cartDims(2)/2, G.cartDims(2)], [G.cartDims(3)/2, G.cartDims(3)]);
    colorbar;set(gca, 'ZDir', 'reverse');axis tight equal;
    title('Rock Permeability')
end


%% fracture planes
tem_index = [5,9,13,17,21];
tem_index_ver = fliplr(tem_index);
for i = 1:length(the_input.frac.index);
    index = the_input.frac.index(i) + tem_index_ver(i);
    centroid = [index*20+10, 100, 50];
    points = [
        [index*20+10, 100-the_input.frac.length(i)/2, 50-the_input.frac.height(i)/2]
        [index*20+10, 100+the_input.frac.length(i)/2, 50-the_input.frac.height(i)/2]
        [index*20+10, 100+the_input.frac.length(i)/2, 50+the_input.frac.height(i)/2]
        [index*20+10, 100-the_input.frac.length(i)/2, 50+the_input.frac.height(i)/2],
        ];
    fracplanes(i).points = points;
    fracplanes(i).aperture = the_input.frac.width(i);
    fracplanes(i).poro = the_input.frac.poro(i);
    fracplanes(i).perm = the_input.frac.perm(i);
end

                        
if checkPlot(the_input, the_input.plot.fracplanes)
    figure;
    plotfracongrid(G,fracplanes);
    view(30,30);
    title('Fracplanes');
end


%% fracture assem
[G,fracplanes]=EDFMshalegrid(G,fracplanes,'Tolerance',the_input.frac.tol,'fracturelist',1:length(the_input.frac.index));

if checkPlot(the_input, the_input.plot.fracgrid)
    figure;
    plotGrid(cartGrid([1 1 1],the_input.G.physdim),'facealpha',0);
    hold on;
    plotGrid(G,G.Matrix.cells.num+1:G.cells.num);
    axis tight equal
    title('Fracture Grid')
    view(30,30);
end


%% wells points
wells = struct;

for i = 1:the_input.well.num
    well(i).radius = the_input.well.radius;
end

% index:[5,9,13,17,21]
% ind:[102, 602, [605,609,613,617,621], 623]
% sub:[[5,4,2], [9,4,2], [13,4,2], [17,4,2], [21,4,2]]

tem_cellprod = [102, 602];
tem_frac = [605,609,613,617,621];
tem_frac_ver = fliplr(tem_frac);
tem_length = length(the_input.frac.index);
tem_value = the_input.frac.index + tem_frac_ver(1:tem_length);
tem_value_ver = fliplr(tem_value);
tem_cellprod = [tem_cellprod, tem_value_ver];
tem_cellprod = [tem_cellprod, 623];

% [nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
% cellinj = 1:nx*ny:(1+(nz-2)*nx*ny);
cellprod = tem_cellprod
% wells(1).points = y_cellind2cellcent(G, cellinj);
wells(1).points = y_cellind2cellcent(G, cellprod);


%% NNC
G = fracturematrixShaleNNC3D(G, the_input.frac.tol);

if checkPlot(the_input, the_input.plot.NNC1)
    figure;
    plotfracongrid(cartGrid([1 1 1],the_input.G.physdim),fracplanes);
    hold on;
    %plotGrid(G,G.nnc.cells(:,1),'facealpha',0,'edgealpha',1);
    axis tight equal
    title('NNC1')
    view(30,45);
end

[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,the_input.frac.tol,'Verbose', true);

if checkPlot(the_input, the_input.plot.NNC2)
    fraccells = G.nnc.cells(strcmp(G.nnc.type,'fracfrac'),:);
    figure;
    plotGrid(cartGrid([1 1 1],the_input.G.physdim),'facealpha',0);
    hold on;
    plotGrid(G,G.Matrix.cells.num+1:G.cells.num,'facealpha',0,'edgealpha',0.5);
    plotGrid(G,fraccells);
    axis equal tight;
    title('NNC2')
    view(30,45);
end

if strcmp(the_input.mechanism.NNC3, 'on')
    [G,wells] = wellfractureNNCs3D(G,fracplanes,wells,the_input.frac.tol);
end


%% fluid
if strcmp(the_input.flash.all, 'on')
    [fluid, info] = getShaleCompFluidCase(the_input.flash.casename);
    G1cell = cartGrid([1 1],[1 1]);
    G1cell = computeGeometry(G1cell);
    EOSModel = EquationOfStateModel(G1cell, fluid, the_input.flash.eosname);
    [~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(the_input.flash.p_sc, the_input.flash.T_sc, info.initial, EOSModel);
    the_input.fluid.rho = [1000, rhoO_S, rhoG_S];
end

flowfluid = initSimpleADIFluid('phases', the_input.fluid.phases, ...
                               'mu' ,the_input.fluid.mu, ...
                               'rho', the_input.fluid.rho, ...
                               'n'  , the_input.fluid.n);
 

                          
if strcmp(the_input.mechanism.compressibility, 'on')
    fluid.bW = @(p) exp((p - the_input.flash.pRef)*the_input.fluid.c_w);
    fluid.bO = @(p) exp((p - the_input.flash.pRef)*the_input.fluid.c_o);
    fluid.bG = @(p) exp((p - the_input.flash.pRef)*the_input.fluid.c_g);
end


%% mechanism
if strcmp(the_input.mechanism.sorption, 'on')
    G.rock.shaleMechanisms.sorption = 1;
end

if strcmp(the_input.mechanism.diffusion, 'on')
    G.rock.shaleMechanisms.diffusion = 1;
    G.rock.Di=[2.8,2.5,1.9]*10^-7;
    G.rock.tau = 2;
end


%% model
gravity reset off

model = NatVarsShaleModel(G, [], flowfluid, fluid, 'water', true);

TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, the_input.frac.tol);
model.operators = TPFAoperators;

%% wells
%{
W   = addWell([], G.Matrix, G.Matrix.rock, cellinj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', sum(model.operators.pv)/the_input.schedule.tot_time*0.7, 'Radius', the_input.well.radius, 'Comp_i', [1, 0, 0], 'Name', 'Injector', ...
    'Sign', 1);
%}

W   = addWell([], G.Matrix, G.Matrix.rock, cellprod, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 50*barsa, 'Radius', the_input.well.radius, 'Comp_i', [1, 1, 0] ./ 2, 'Name', 'Producer', ...
    'Sign', -1);

%W(1).components = info.injection;
W(1).components = info.initial;

if checkPlot(the_input, the_input.plot.wells)
    plotEDFMgrid(G);
    hold on;
    plotWell(G,W);view(60, 20);
end


%% state
if strcmp(the_input.flash.pres_mode, 'file')
    pres = load(the_input.flash.pres_file).pres(1:prod(the_input.G.cartdim));
    for i=1:length(the_input.frac.index)
        data = getfield(G.FracGrid, ['Frac', num2str(i)]);
        index = data.matrix_connection.cells;
        pres_frac = pres(index);
        pres = vertcat(pres, pres_frac);
    end
end
if strcmp(the_input.state.sat_mode, 'file')
    sat = load(the_input.state.sat_file).sat(1:prod(the_input.G.cartdim), :);
    for i=1:length(the_input.frac.index)
        data = getfield(G.FracGrid, ['Frac', num2str(i)]);
        index = data.matrix_connection.cells;
        sat_frac = sat(index, :);
        sat = vertcat(sat, sat_frac);
    end
end

state = initCompositionalState(G, pres, info.temp, sat, info.initial, model.EOSModel);

%% schedule
schedule = simpleSchedule(the_input.schedule.tot_time, 'W', W);


%% simulate
[ws, states, report] = simulateScheduleAD(state, model, schedule);


%% save
pres = states{1}.pressure;
sat = states{1}.s;
qOs = ws{1}.qOs;
save(the_input.flash.pres_file, "pres");
save(the_input.state.sat_file, "sat");
save(the_input.well.prod_file, 'qOs')
save(the_input.G.G_file, 'G')