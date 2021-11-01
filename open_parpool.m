
disp(['PPN: ',num2str(no_par_processes)]);


if usejava('desktop')
    disp('Local Job')
    mycluster = 'local';
    mypool = parcluster(mycluster);
    cp=gcp('nocreate');
    if ~isempty(cp)
        if  cp.NumWorkers~=no_par_processes
            delete(cp);
            if no_par_processes>0
                parpool(mypool, no_par_processes);
            end
        end
    else
        if no_par_processes>0
            parpool(mypool, no_par_processes);
        end
    end
    
else
    
   cp=gcp('nocreate');
   if ~isempty(cp)
       if contains(getenv('QUEUE'), 'aws-')
           disp('AWS Job')
           mycluster = 'local';
           myfolderbase = strcat(getenv('TMP'), '/');
       else
           disp('HPCC Job')
           mycluster = 'WhartonHPC';
           myfolderbase = '~/matlabtmp/';
       end
       mypool = parcluster(mycluster);
       tmpFolder = strcat(myfolderbase, getenv('JOB_ID'), '-', getenv('SGE_TASK_ID'));
       mkdir(tmpFolder);
       mypool.JobStorageLocation = tmpFolder;
       parpool(mypool, no_par_processes);
   else
       disp('Cluster job with existing pool');
   end
end
    






