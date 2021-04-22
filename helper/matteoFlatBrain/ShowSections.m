function ax = ShowSections( TheVolume, iX, iZ, DecFac )

if nargin < 4
    DecFac = 1;
end

nn = size(TheVolume);

figure; clf; ax = [];
ax(1) = subplot(1,2,1); imagesc( squeeze( TheVolume( iX/DecFac,:,:) ) ); axis image; hold on
ax(2) = subplot(1,2,2); imagesc( squeeze( TheVolume( :,:,iZ/DecFac))' ); axis image; hold on
xlabel(ax(1),'Dimension 3 (ML)');
ylabel(ax(1),'Dimension 2 (DV)');
xlabel(ax(2),'Dimension 1 (AP)');
axes(ax(1)); hold on; plot(iZ/DecFac*[1 1], [1 nn(2)],'k-')
axes(ax(2)); hold on; plot(iX/DecFac*[1 1], [1 nn(2)],'k-')