function [fluid, info] = getBenchmarkMixtureOG(name, varargin)
descr = name;
switch(lower(name))
    case 'eec'

        T_c = [33.145, 190.564, 304.1282, 126.192, 639,658];
        P_c = [12.9583, 45.3886, 72.8645, 33.5554, 19.245,17.962]*atm;
        mw = [2.01588, 16.0428, 44.0098, 28.0135, 156.312,170.338]/1000;
        acc = [-0.215993, 0.0115478, 0.223621, 0.0377215, 0.530316,0.576385];
        Z_c = [0.305, 0.286, 0.27588, 0.289, 0.242,0.238];

        V_c = Z_c.*8.314.*T_c./P_c;

        names = {'Hydrogen', 'Methane', 'Carbondioxide','Nitrogen','n-Undecane','n-Dodecane'};
        fluid = CompositionalMixture(names, T_c, P_c, V_c, acc, mw);

        bic = [0,       0,       0,       0,      0,               0;...
               0.202,   0,       0,       0,      0,               0; ...
               0.1202,  0.1,     0,       0,      0,               0; ...
               -0.036,  0.036,   -0.02,   0,      0,               0; ...
               0.2921,  0.04799, 0.101,   0.1,    0,               0; ...
               0.3560,  0.05180, 0.101,   0.1,    0.00008,         0];


        bic = bic + tril(bic, -1)';

        fluid = fluid.setBinaryInteraction(bic);
        info = makeInfo('injection', [1, 0, 0, 0, 0,  0], ...
            'initial',   [0.00, 0.31, 0.31, 0.32, 0.03, 0.03], ...
            'pressure', 50*barsa, ...
            'T', 323.15*Kelvin);
end
fluid.name = descr;
end

function info = makeInfo(varargin)
info = struct('injection', [], 'initial', [],...
    'pressure',  [], 'temp',    [], 'T', []);
info = merge_options(info, varargin{:});
info.temp = info.T;
end