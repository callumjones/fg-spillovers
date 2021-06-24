function gen_set_params(SET,dyn_in)

M_ = dyn_in.M_ ;

%tochg = [SET.tv_str_chg.ant_chg_params SET.tv_str_chg.unant_chg_params] ;

tmp = cellstr(dyn_in.M_.param_names) ;

index = 1:M_.param_nbr; 

tmpstr = '' ;
for i=1:M_.param_nbr
    tmpstr = [tmpstr, tmp{i},' = dyn_in.M_.params(', int2str(i), '); \n'] ;
end

old = [' \n', ... % put it here
' '];

% for i=1:length(tochg)
%     idx = index(find(strcmp(tmp, tochg{i}))) ;
%     new = strrep(old, tochg{i}, ['dyn_in.M_.params(', int2str(idx), ')']) ;
%     old = new ;
% end

tmpstr2 = '' ;
for i=1:M_.param_nbr
    tmpstr2 = [tmpstr2, ' out(', int2str(i), ')=',tmp{i},';'] ;
end

old = [tmpstr, old, tmpstr2] ;

fid = fopen('set_params.m','wt');

fprintf(fid, old) ;

fclose(fid);