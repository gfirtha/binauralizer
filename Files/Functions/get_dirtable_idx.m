function idx = get_dirtable_idx( dir_tables, source )
if isempty(dir_tables)
    idx = [];
else
    a = cellfun( @(x) strcmp( source.source_type.Shape,x),...
        cellfun( @(x) x.type.Shape, dir_tables, 'UniformOutput', false) );
    b = (source.source_type.R(1) == cell2mat(cellfun( @(x) x.type.R(1), dir_tables, 'UniformOutput', false)));
    idx = find(a&b);
end
end

