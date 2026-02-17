%% Setup paths
paths = {'Models';
         'Models/aero';
         'Models/act';
         'Models/mass';
         'Models/atmo';
         'Models/earth';
         'Models/eqmo';
         'Models/prop';
         'Models/vehicle';
         'Trim'
        };
for idx = 1:length(paths)
    addpath(paths{idx});
end
clear paths idx

%% Run vars scripts 
aero_vars;
act_vars;
mass_vars;
atmo_vars;
earth_vars;
prop_vars;
vehicle_vars;

%% Initialize Autocode Folder
% Create Autocode folder if it does not exist
if ~exist('./Autocode','dir')
    mkdir('Autocode')
end
CODEpath = [cd '/Autocode'];
% Define single code generation location for all auto-generated files
Simulink.fileGenControl('set','CacheFolder',CODEpath,'CodeGenFolder',CODEpath);