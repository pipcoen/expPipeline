function dateRange = keyDates(subjectID, dataTag)
switch lower(dataTag{1})
    case 'behavior'
        switch subjectID
            case 'all'; dateRange = {'PC011','PC012','PC013','PC015','PC022','PC027','PC029','DJ007','DJ008','DJ010'}';
            case 'PC011'; dateRange = {'rng', '2017-06-14', '2017-08-16'};
            case 'PC012'; dateRange = {'rng', '2017-09-06', '2017-10-17'}; %Regular trials within inactivation
            case 'PC013'; dateRange = {'rng', '2017-09-05', '2017-10-17'}; %Regular trials within inactivation
            case 'PC015'; dateRange = {'rng', '2017-09-25', '2017-10-17'};
            case 'PC022'; dateRange = {'rng', '2018-02-20', '2018-05-11'}; %Regular trials within inactivation
            case 'PC027'; dateRange = {'rng', '2018-02-05', '2018-05-09'}; %Regular trials within inactivation
            case 'PC029'; dateRange = {'rng', '2018-02-27', '2018-06-21'};
            case 'DJ008'; dateRange = {'rng', '2018-05-11', '2018-05-30'};
            case 'DJ007'; dateRange = {'rng', '2018-05-11', '2018-06-19'}; %Intermediate days without conflicts are excluded
            case 'DJ010'; dateRange = {'rng', '2018-06-09', '2018-06-18'};
            otherwise, dateRange = dataTag;
                %PC022 dateRange = {'rng', '2017-11-02', '2017-11-21'}; %Data without any inactivation, but fewer trials
        end
    case 'uniscan'
        switch subjectID{1}
            case 'all'; dateRange = {'PC022','PC027','PC029','DJ008'}';
            case 'PC022'; dateRange = {'rng', '2018-02-06', '2018-03-20'}; %Regular trials within inactivation
            case 'PC027'; dateRange = {'rng', '2018-02-05', '2018-03-21'}; %Regular trials within inactivation
            case 'PC029'; dateRange = {'rng', '2018-06-12', datestr(now, 'yyyy-mm-dd')};
            case 'DJ008'; dateRange = {'rng', '2018-06-12', datestr(now, 'yyyy-mm-dd')};
            otherwise, dateRange = dataTag;
                %PC022 dateRange = {'rng', '2017-11-02', '2017-11-21'}; %Data without any inactivation, but fewer trials
        end
    otherwise, dateRange = dataTag;
end
end