if usejava('desktop')
   clear; 
else
    ver;
end
close all;

% name of file with experiment definitions
experdef='20200904';
experdef_file=['experdef_',experdef];

% name of main experiment (needs to have been computed to convergence)
res_path='res_20200904_lowratefinal';

% For the transition of the baseline model
experseq = { 'trans01_eta';
              'trans02_eta';
              'trans03_eta';
              'trans04_eta';
              'trans05_eta';
              'trans06_eta';
              'trans07_eta';
              'trans08_eta'};
          
% experseq = { 'trans01_secur_DTI';
%               'trans02_secur_DTI';
%               'trans03_secur_DTI';
%               'trans04_secur_DTI';
%               'trans05_secur_DTI';
%               'trans06_secur_DTI';
%               'trans07_secur_DTI';
%               'trans08_secur_DTI'};          
         
         
% write names of res files into cell array
% For the one where we take all periods until convergence
out_array=['ST_',res_path,'_names.mat'];

         
guess_path=res_path;
guess_mode='guess';
nexp=length(experseq);
maxit=30;
checkTrans=1;
tolavg=1e-5;
npp=16;


res_array=cell(nexp+1,1);
res_array{nexp+1}=res_path; % final step

% initial value and forecasting functions
load(res_path,'mobj');
nextvfct=mobj.Vfct;
nextffct=mobj.Ffct;
clear mobj;


% loop over all experiments
for ie=nexp:-1:1
    % create experiment definition
    expername=experseq{ie};
    main_create_env(experdef_file,expername,guess_mode,guess_path);
    
    
    % Name of the file:
    exper_path=['env_',expername,'.mat'];
    outname=['res_',experdef,'_',expername,'.mat'];
    
   
    mobj=main_run_exper(exper_path,maxit,outname,checkTrans,tolavg,npp);
    % set guess for next period
    guess_path=outname;
    res_array{ie}=outname;
end


save(out_array,'res_array');
