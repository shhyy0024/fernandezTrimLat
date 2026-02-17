%% Close all open models
disp('Closing Simulink models...')
while 1
    if isempty(gcs),break,end
    try
        close_system
    catch
        choice = questdlg(['Save "',gcs,'" before closing?'],...
            'Message','Yes','No','Cancel','Cancel');
        switch choice
            case 'Yes'
                save_system(gcs)
                close_system(gcs)
            case 'No'
                bdclose(gcs)
            case 'Cancel'
                error('Script terminated by user')
        end
    end
end
disp('Done...')
%% Delete any autosave file
delete *.autosave

%% Remove paths
warning off
rmpath(genpath(cd))
warning on

%% Reset location of auto generated files
Simulink.fileGenControl('reset');

%% Clear variables 
clear all
clear mex
warning('on','all')
disp('Shutdown Complete')

