function sArr = chkThenRemoveFields(sArr, flds)
aFld = fields(sArr);
sArr = rmfield(sArr, aFld(ismember(aFld, flds)));
end