function updateEphysRecord(penetrationIdx, calcLine, calcDepth)
%%
ephysListPath = prc.pathFinder('ephysrecord');
ephysRec = table2struct(readtable(ephysListPath));
columns = fields(ephysRec);
selectedRow = find([ephysRec.penetrationIdx]==penetrationIdx)+1;
lineFitColumn = char('A' + find(contains(columns, 'calcLine'))-1);
depthFitColumn = char('A' + find(contains(columns, 'calcDepth'))-1);
xlswrite(ephysListPath, calcLine, 'sheet1', [lineFitColumn num2str(selectedRow)]);
xlswrite(ephysListPath, calcDepth, 'sheet1', [depthFitColumn num2str(selectedRow)]);
