function dateRange = keyDates(subjectID, dataTag)
if iscell(dataTag{1}); dateRange = dataTag; return; end
switch lower(dataTag{1})
    case 'behavior'
        switch subjectID{1}
            case 'all'; dateRange = {'PC011','PC012','PC013','PC015','PC022','PC027','PC029','DJ007','DJ008','DJ010'}';
            case 'PC011'; dateRange = {'rng', '2017-06-14', '2017-08-16'};
            case 'PC012'; dateRange = {'rng', '2017-09-06', '2017-10-17'}; %Regular trials within inactivation
            case 'PC013'; dateRange = {'rng', '2017-09-05', '2017-10-17'}; %Regular trials within inactivation
            case 'PC015'; dateRange = {'rng', '2017-09-25', '2017-10-17'};
            case 'PC022'; dateRange = {'rng', '2018-02-20', '2018-05-11'}; %Regular trials within inactivation
            case 'PC027'; dateRange = {'rng', '2018-02-05', '2018-05-09'}; %Regular trials within inactivation
            case 'PC029'; dateRange = {'rng', '2018-02-27', '2018-06-21'};
            case 'PC031'; dateRange = {'rng', '2019-02-19', '2019-03-01'};
            case 'DJ007'; dateRange = {'rng', '2018-05-11', '2018-06-19'}; %Intermediate days without conflicts are excluded
            case 'DJ008'; dateRange = {'rng', '2018-05-11', '2018-05-30'};
            case 'DJ010'; dateRange = {'rng', '2018-06-09', '2018-06-18'};
            otherwise, dateRange = dataTag;
                %PC022 dateRange = {'rng', '2017-11-02', '2017-11-21'}; %Data without any inactivation, but fewer trials
        end
    case 'aud5'
        switch subjectID{1}
            case 'all'; dateRange = {'PC013','PC015','PC022'}';
            case 'PC013'; dateRange = {'rng', '2017-10-18', '2017-11-04'}; %Regular trials within inactivation
            otherwise, dateRange = dataTag;
                %PC022 dateRange = {'rng', '2017-11-02', '2017-11-21'}; %Data without any inactivation, but fewer trials
        end
    case 'uniscan'
        switch subjectID{1}
            case 'all'; dateRange = {'PC027','PC029','DJ008','DJ006','DJ007'}';
                %             case 'PC022'; dateRange = {'rng', '2018-02-06', '2018-03-20'}; %Missing some of the locations and power was only 1.5mW
            case 'PC027'; dateRange = {'rng', '2018-02-05', '2018-03-21'}; %Power was only 1.5mW
            case 'PC029'; dateRange = {'rng', '2018-06-12', '2018-07-17'};
            case 'DJ008'; dateRange = {'rng', '2018-06-12', '2018-07-17'};
            case 'DJ006'; dateRange = {'rng', '2018-08-06', '2018-09-16'};
            case 'DJ007'; dateRange = {'rng', '2018-08-06', '2018-10-18'};
            otherwise, dateRange = dataTag;
                %PC022 dateRange = {'rng', '2017-11-02', '2017-11-21'}; %Data without any inactivation, but fewer trials
        end
        
    case 'biscan'
        switch subjectID{1}
            case 'all'; dateRange = {'PC010','PC012','PC013'}';
            case 'PC010'; dateRange = {'rng', '2017-06-20', '2017-07-05'}; %Power was only 1.5mW
            case 'PC012'; dateRange = {'rng', '2017-06-21', '2017-07-05'};
            case 'PC013'; dateRange = {'rng', '2017-06-20', '2017-07-06'};
            otherwise, dateRange = dataTag;
                %PC022 dateRange = {'rng', '2017-11-02', '2017-11-21'}; %Data without any inactivation, but fewer trials
        end
    otherwise, dateRange = dataTag;
end
end