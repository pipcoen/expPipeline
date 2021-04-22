addpath(genpath('../allenCCF'));
addpath(genpath('../npy-matlab'));

st = loadStructureTree('../Allen/structure_tree_safe_2017.csv'); % a table of what all the labels mean

tv = readNPY('../Allen/template_volume_10um.npy'); % grey-scale "background signal intensity"

av = readNPY('../Allen/annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below

allen_ccf_npx(tv,av,st);

 file_save_location = 'C:\Histology\Mouse1'; % where will the probe locations be saved
probe_name = 'test'; % name probe to avoid overwriting

f = allenAtlasBrowser(tv, av, st, file_save_location, probe_name);

plotBrainGrid();
