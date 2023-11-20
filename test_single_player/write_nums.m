function write_nums(var,filename)
    fid = fopen(filename, 'w');
    for i = 1:size(var,1)
         b = fprintf(fid, "%d,",var(i));
    end
    fclose(fid);
end

