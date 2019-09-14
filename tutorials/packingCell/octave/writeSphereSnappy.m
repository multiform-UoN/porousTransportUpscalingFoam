function writeSphereSnappy (fid,centre,rad,id)


fprintf(fid,'sphere%i\n',id);
fprintf(fid,'{\n');
fprintf(fid,'    surface sphere;\n');
fprintf(fid,'    scale (%f %f %f);\n', rad,rad,rad);
fprintf(fid,'    transform\n');
fprintf(fid,'    {\n');
fprintf(fid,'         coordinateSystem\n');
fprintf(fid,'         {\n');
fprintf(fid,'             type cartesian;\n');
fprintf(fid,'             origin (%f %f %f);\n',centre(1),centre(2),centre(3));
fprintf(fid,'             coordinateRotation\n');
fprintf(fid,'             {\n');
fprintf(fid,'                 type axesRotation;\n');
fprintf(fid,'                 e1 (1 0 0);\n');
fprintf(fid,'                 e2 (0 0 1);\n');
fprintf(fid,'             }\n');
fprintf(fid,'         }\n');
fprintf(fid,'    }\n');

fprintf(fid,'}\n\n');

  
endfunction