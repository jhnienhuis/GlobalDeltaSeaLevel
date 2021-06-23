function create_netcdf(ncname,out,funits,fmeta)

fnames = fieldnames(out);

ncid = netcdf.create(ncname,'WRITE');

dim_a = netcdf.defDim(ncid,'n0',size(out.(fnames{1}),1));
dim_b = netcdf.defDim(ncid,'n1',size(out.(fnames{1}),2));


for ii=1:length(fnames),
    
    if isa(out.(fnames{ii}),'integer'),
        out.(fnames{ii}) = int32(out.(fnames{ii}));
        cla = 'NC_INT';
    else,
        cla = class(out.(fnames{ii}));
    end
        
    
    if size(out.(fnames{ii}),1)==10848, 
        dim1 = dim_a;
    else,
        dim1 = netcdf.defDim(ncid,['n' num2str(ii) '_1'],size(out.(fnames{ii}),1));
    end
    
    if size(out.(fnames{ii}),2)==1, 
        dim2 = dim_b;
    else,
        dim2 = netcdf.defDim(ncid,['n' num2str(ii) '_2'],size(out.(fnames{ii}),2));
    end
        
    varid = netcdf.defVar(ncid,fnames{ii},cla,[dim1 dim2]);
        

    
    %ncwrite(ncname,fnames{ii},out.(fnames{ii}))
    netcdf.endDef(ncid);
    netcdf.putVar(ncid,varid,out.(fnames{ii}));
    netcdf.reDef(ncid);
    netcdf.putAtt(ncid,varid,'units',funits{ii});
    netcdf.putAtt(ncid,varid,'about',fmeta{ii});
    
end
    

netcdf.close(ncid);
