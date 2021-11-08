function pVal = chiTest(nRes, nTotal)

if (length(nRes)~= 2)||(length(nRes)~=length(nTotal))
    disp('Error: bad vector length')
elseif (nRes(1)>nTotal(1))|| (nRes(2)>nTotal(2))
    disp('Error: bad counts (nRes>nTotal)')
else  
   
    % Observed data
    n1 = nRes(1);
    n2 = nRes(2);
    nTotal1 = nTotal(1);
    nTotal2 = nTotal(2);
    
    % Pooled estimate of proportion
    p0 = (n1+n2) / (nTotal1+nTotal2);
    
    % Expected counts under H0 (null hypothesis)
    n10 = nTotal1 * p0;
    n20 = nTotal2 * p0;
    observed = [n1 nTotal1-n1 n2 nTotal2-n2];
    expected = [n10 nTotal1-n10 n20 nTotal2-n20];
    
    % Standard Chi-square test
    chi2stat = sum((observed-expected).^2 ./ expected);
    pVal = 1 - chi2cdf(chi2stat,1);
    
end