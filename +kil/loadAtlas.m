allenPath = prc.pathFinder('allenAtlas');
tv = readNPY([allenPath 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allenPath 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allenPath 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
