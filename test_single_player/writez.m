function writez(var,filename)
    fid = fopen(filename, 'w');
    for i = 1:size(var,1)
         b = fprintf(fid, "%.9f,",var(i));
    end
    fclose(fid);
end

