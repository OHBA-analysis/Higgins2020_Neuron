
%%%%%%%%%%%%%%%%%%
%% SET UP FILE NAMES

rawdatadir = 
datadir     = [rawdatadir,'/'];%[tilde 'homedir/vols_data/notts_ukmp/raw_data/'];
africadir   = [spmfilesdir '/africa/'];

subjectdirs = dir([datadir '/3*']);
subjectdirs = sort_nat({subjectdirs.name});
africadir   = [spmfilesdir 'africa/'];

subjectdirs = dir([datadir '/3*']);
subjectdirs = sort_nat({subjectdirs.name});

ctf_files_eo            = cell(size(subjectdirs));
ctf_files_vmg           = cell(size(subjectdirs));
ctf_files_vml           = cell(size(subjectdirs));
ctf_files_vms           = cell(size(subjectdirs));
spm_files_eo            = cell(size(subjectdirs));
spm_files_vmg           = cell(size(subjectdirs));
spm_files_vml           = cell(size(subjectdirs));
spm_files_vms           = cell(size(subjectdirs));
structural_files        = cell(size(subjectdirs));
pos_files               = cell(size(subjectdirs));

do_mri_convert = false;

for s = 1:length(subjectdirs)
    fprintf(['Processing subject ',int2str(s),'\n']); 
    dsfile = dir([datadir subjectdirs{s} '/*Eyes_Open*.ds']); 
    if ~isempty(dsfile)
        ctf_files_eo{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_eo{s}    = [spmfilesdir subjectdirs{s} '_eo.mat'];
    end
    
    dsfile = dir([datadir subjectdirs{s} '/*VisMotor_Gamma*.ds']); 
    if ~isempty(dsfile)
        ctf_files_vmg{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_vmg{s}    = [spmfilesdir subjectdirs{s} '_vmg.mat'];
    end
    
    dsfile = dir([datadir subjectdirs{s} '/*VisMotor_Long*.ds']); 
    if ~isempty(dsfile)
        ctf_files_vml{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_vml{s}    = [spmfilesdir subjectdirs{s} '_vml.mat'];
    end
 
    
    dsfile = dir([datadir subjectdirs{s} '/*VisMotor_Short*.ds']); 
    if ~isempty(dsfile)
        ctf_files_vms{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_vms{s}    = [spmfilesdir subjectdirs{s} '_vms.mat'];
    end
        
    % structural files
    try
        
        niifile = [datadir subjectdirs{s} '/' subjectdirs{s} '_CRG']; 
        
        if do_mri_convert || isempty(dir([niifile '.nii']))
            runcmd(['rm -f ' niifile '.nii']);
            runcmd(['rm -f ' niifile '.nii.gz']);
            runcmd(['rm -f ' datadir subjectdirs{s} '/y_' subjectdirs{s} '_CRG.nii']);
            runcmd(['rm -f ' niifile '_*.gii']);
            runcmd(['rm -f ' niifile '*_*.*']);
            runcmd(['rm -f ' datadir subjectdirs{s} '/' 'rhino*']);
            runcmd(['rm -f ' datadir subjectdirs{s} '/' 'y_rhino*']);
            
            mrifile = dir([datadir subjectdirs{s} '/' num2str(subjectdirs{s}) '_CRG.mri']); 
            
            if isempty(mrifile)
                me=MException('Preprocess:noStuctural',['No mri or nii structural file for subject ' subjectdirs{s}]);
                throw(me);
            end           
            
            %[ niifile ] = mri2analyze( [datadir subjectdirs{s} '/' mrifile.name] );
            
            
            niifilename=[niifile '.nii'];
            convert_mri([datadir subjectdirs{s} '/' mrifile.name], niifilename);
           
        end
        
        niifilename=[niifile '.nii'];
        structural_files{s} = niifilename;
        
    catch me
        warning(me.message)
        structural_files{s} = [];
    end
    
    % list of head position files
    
    pfname=[datadir 'pos_files/' subjectdirs{s} '.pos'];
    pf = dir(pfname); 
    if ~isempty(pf)
        pos_files{s}=pfname;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select spm files to work with

spm_files=spm_files_eo;
ctf_files=ctf_files_eo;
event_type=[];
    

subjects_to_do = find(~cellfun(@isempty,spm_files) & ~cellfun(@isempty,structural_files));

%% consolidate into is object:

is = [];
is.spm_files = spm_files_eo;
is.ctf_files = ctf_files; %rawdata files
is.pos_files = pos_files;
is.structural_files = structural_files;
is.subjects_to_do = subjects_to_do;

clear spm_files ctf_files pos_files structural_files spm_files_eo ctf_files_eo spm_files_vml ctf_files_vml ...
    spm_files_vms ctf_files_vms spm_files_vmg ctf_files_vmg event_type subjectdirs subjects_to_do