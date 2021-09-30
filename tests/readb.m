function a0 = readb(filename)
fileid = fopen(filename);
a0 = fscanf(fileid, '%f');
fclose(fileid);
