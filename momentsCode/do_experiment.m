fprintf('Now Testing Gehler Moments 9 ... \n');
do_graham_moments_gehler;
fid = fopen('result.txt','w');
fprintf(fid,strcat('Moments 9 : \n','Mean ',mat2str(mea,5),'\nStd ',mat2str(std1,4),'\n\n'));
fclose(fid);

fprintf('Now Testing Gehler Moments 19 ...  \n');
do_graham_moments_gehler_19;
fid = fopen('result.txt','a');
fprintf(fid,strcat('Moments 19 : \n','Mean ',mat2str(mea,5),'\nStd ',mat2str(std1,4),'\n\n'));
fclose(fid);

fprintf('Now Testing Gehler Edge Moments 9 ... \n');
do_graham_moments_gehler_edges;
fid = fopen('result.txt','a');
fprintf(fid,strcat('Edge Moments 9 : \n','Mean ',mat2str(mea,5),'\nStd ',mat2str(std1,4),'\n\n'));
fclose(fid);

fprintf('Now Testing Gehler Edge Moments 19 ... \n');
do_graham_moments_gehler_edges_19;
fid = fopen('result.txt','a');
fprintf(fid,strcat('Edge Moments 19 : \n','Mean ',mat2str(mea,5),'\nStd ',mat2str(std1,4),'\n\n'));
fclose(fid);

fprintf('Now Testing Gehler Chrom & Bright ... \n');
do_graham_moments_gehler_chrom_combo;
fid = fopen('result.txt','a');
fprintf(fid,strcat('Chrome & Bright 8+8 : \n','Mean ',mat2str(mea,5),'\nStd ',mat2str(std1,4),'\n\n'));
fclose(fid);

fprintf('Now Testing Gehler Chrom & Bright & Zeta ...  \n');
do_graham_moments_gehler_chromall_combo;
fid = fopen('result.txt','a');
fprintf(fid,strcat('Chrome & Bright & Zeta 8+8 : \n','Mean ',mat2str(mea,5),'\nStd ',mat2str(std1,4),'\n\n'));
fclose(fid);
