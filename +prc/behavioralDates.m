function dateRange = behavioralDates(subjectID)
switch subjectID
    case 'PC010'
        dateRange = {'rng', '2017-09-05', '2017-10-13'}; %Regular trials within inactivation
    case 'PC011'
        dateRange = {'rng', '2017-06-14', '2017-08-16'};
    case 'PC012'
        dateRange = {'rng', '2017-09-06', '2017-10-17'}; %Regular trials within inactivation
    case 'PC013'
        dateRange = {'rng', '2017-09-05', '2017-10-17'}; %Regular trials within inactivation
    case 'PC015'
        dateRange = {'rng', '2017-09-25', '2017-10-17'};
    case 'PC022'
        dateRange = {'rng', '2017-11-02', '2017-11-21'};
    case 'PC025'
        dateRange = {'rng', '2018-03-08', '2018-03-30'};
    case 'PC027'
        dateRange = {'rng', '2018-02-05', '2018-03-30'};
    case 'PC029'
        dateRange = {'rng', '2018-02-26', '2018-03-30'};
    case 'DJ008'
        dateRange = {'rng', '2018-05-11', '2018-05-30'};
    case 'DJ007'
        dateRange = {'rng', '2018-06-07', '2018-06-18'};
    case 'DJ010'
        dateRange = {'rng', '2018-06-09', '2018-06-18'};
end