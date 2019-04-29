function dateRange = keyDates(subjectID, dataTag)
if iscell(dataTag{1}); dateRange = dataTag; return; end
behaviorMice = {'PC011','PC012','PC013','PC015','PC022','PC027','PC029','PC030','PC031','PC032','PC033','PC034','PC043',...
                    'DJ007','DJ008','DJ010'};
                
switch lower(dataTag{1})
    case {'behavior'; 'behaviour'}
        switch subjectID{1}
            case 'all'; dateRange = behaviorMice';
            case 'PC011'; dateRange = {'rng', '2017-06-14', '2017-08-16'};
            case 'PC012'; dateRange = {'rng', '2017-09-06', '2017-10-17'}; %Regular trials within inactivation
            case 'PC013'; dateRange = {'rng', '2017-09-05', '2017-10-17'}; %Regular trials within inactivation
            case 'PC015'; dateRange = {'rng', '2017-09-25', '2017-10-17'};
            case 'PC022'; dateRange = {'rng', '2018-02-20', '2018-05-11'}; %Regular trials within inactivation
            case 'PC027'; dateRange = {'rng', '2018-02-05', '2018-05-09'}; %Regular trials within inactivation
            case 'PC029'; dateRange = {'rng', '2018-02-27', '2018-06-21'};
            case 'PC030'; dateRange = {'rng', '2019-02-12', '2019-03-29'};
            case 'PC031'; dateRange = {'rng', '2019-02-19', '2019-03-01'};
            case 'PC032'; dateRange = {'rng', '2019-02-12', '2019-03-21'};
            case 'PC033'; dateRange = {'rng', '2019-02-27', '2019-03-15'};
            case 'PC034'; dateRange = {'rng', '2019-02-19', '2019-03-15'}; %Issues with eye after this point
            case 'PC043'; dateRange = {'rng', '2019-03-06', '2019-03-29'}; %Issues with eye after this point
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
        
    case 'm2ephys'
        switch subjectID{1}
            case 'all'; dateRange = {'DJ007','PC029','PC033','PC032'}';
            case 'DJ007'; dateRange = {'rng', '2018-11-28', '2018-12-02'}; %Power was only 1.5mW
            case 'PC029'; dateRange = {'rng', '2018-10-18', '2018-10-19'};
            case 'PC032'; dateRange = {'rng', '2019-04-03', '2019-04-10'};
            case 'PC033'; dateRange = {'rng', '2019-03-27', '2019-03-31'};
            otherwise, dateRange = dataTag;
                %PC022 dateRange = {'rng', '2017-11-02', '2017-11-21'}; %Data without any inactivation, but fewer trials
        end
    otherwise, dateRange = dataTag;
end
end