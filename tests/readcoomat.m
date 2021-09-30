function A0 = readcoomat(filename)
fileid = fopen(filename);
A = fscanf(fileid,'%d %d %f');
fclose(fileid);
A = reshape(A, [3, numel(A)/3])';
A(:,1) = A(:,1) + 1;
A(:,2) = A(:,2) + 1;
A0 = sparse(A(:,1),A(:,2),A(:,3));
