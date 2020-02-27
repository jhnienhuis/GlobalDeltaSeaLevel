function create_netcdf(ncname,out,funits,fmeta)

fnames = fieldnames(out);

n = length(out.(fnames{1}));


ncid = netcdf.create(ncname,'WRITE');
% Define a dimension in the file.
dimid = netcdf.defDim(ncid,'n',n);

for ii=1:length(fnames),
    
    if length(out.(fnames{ii}))~=n, continue, end
        
    %nccreate(ncname,fnames{ii},'Dimensions',{'n',n})
    

    
    % Define a new variable in the file.
    varid = netcdf.defVar(ncid,fnames{ii},'double',dimid);
    
    %ncwrite(ncname,fnames{ii},out.(fnames{ii}))
    netcdf.endDef(ncid);
    netcdf.putVar(ncid,varid,out.(fnames{ii}));
    netcdf.reDef(ncid);
    netcdf.putAtt(ncid,varid,'units',funits{ii});
    netcdf.putAtt(ncid,varid,'about',fmeta{ii});
    
end
    

netcdf.close(ncid);
