% Aircraft design tool
%
% Mario Bras (mbras@uvic.ca) and Ricardo Marques (ricardoemarques@uvic.ca) 2019
%
% This file is subject to the license terms in the LICENSE file included in this distribution

function save_project(data, filename)

str = prettyjson(jsonencode(data));

fid = fopen(filename,'wt');
fprintf(fid, str);
fclose(fid);
