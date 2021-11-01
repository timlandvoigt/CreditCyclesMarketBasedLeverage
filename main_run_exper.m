function mobj=main_run_exper(varargin)

if nargin>0
    exper_path=varargin{1};
    maxit=varargin{2};
    outname=varargin{3};
    if nargin>3
        checkTrans=varargin{4};
        if nargin>4
            tol_avg=varargin{5};
             if nargin>5
                    no_par_processes=varargin{6};
                    if nargin>6
                        printmode=varargin{7};
                        damp=varargin{8};
                    end
             end
        end
    end    
end
% ===========================================
% program control
% ===========================================


% path to file with experiment definition
if ~exist('exper_path','var')
%    exper_path='env_bench.mat';
    exper_path='env_lowratefinal.mat';

end

if ~exist('no_par_processes','var')
%    exper_path='env_bench.mat';
    no_par_processes=2;

end

if ~exist('maxit','var')
    maxit =100;
end
if ~exist('tol_avg','var')
    % mean convergence
    tol_avg=1e-5;
end

if ~exist('printmode','var')
    % print mode
    printmode=1;
end

if ~exist('checkTrans','var')
    % print mode
    checkTrans=0;
end


if ~exist('damp','var')
    % dampening
    damp=0.5;
end

open_parpool; 


% ===========================================
% load environment structure
% ===========================================

load(exper_path);

% ===========================================
% outer convergence loop
% ===========================================

% policy iteration
% polIter is in 
[mobj,failedPoints,dist,distT]=mobj.polIter(maxit,1,printmode,damp,tol_avg,checkTrans);

% make date string
curr_date=clock;
date_str='';
for i=1:5
    date_str=[date_str,'_',num2str(curr_date(i))];
end

if ~exist('outname','var')
    outname=['res',date_str,'.mat'];
end


save(outname,'mobj','stv','failedPoints','dist','distT');

% Close parallel pool
end