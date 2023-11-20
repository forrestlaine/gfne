function M = readz(filename)
    fid = fopen(filename);
    M = textscan(fid, '%f', 'Delimiter', ',');
    M = M{1};
    fclose(fid);
end

