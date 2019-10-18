function dateRange = keyDates(subjectID, dataTag)
                
switch lower(dataTag{1})
    case {'behavior'; 'behaviour'}
        switch subjectID{1}
            case 'PC011'; dateRange = {'2017-06-14:2017-08-16'};
            case 'PC012'; dateRange = {'2017-09-06:2017-10-17'}; %Regular trials within inactivation
            case 'PC013'; dateRange = {'2017-09-05:2017-10-17'}; %Regular trials within inactivation
            case 'PC015'; dateRange = {'2017-09-25:2017-10-17'};
            case 'PC022'; dateRange = {'2018-02-20:2018-05-11'}; %Regular trials within inactivation
            case 'PC027'; dateRange = {'2018-02-05:2018-05-09'}; %Regular trials within inactivation
            case 'PC029'; dateRange = {'2018-02-27:2018-06-21'};
            case 'PC030'; dateRange = {'2019-02-12:2019-03-29'};
            case 'PC031'; dateRange = {'2019-02-19:2019-03-01'};
            case 'PC032'; dateRange = {'2019-02-12:2019-03-21'};
            case 'PC033'; dateRange = {'2019-02-27:2019-03-15'};
            case 'PC034'; dateRange = {'2019-02-19:2019-03-15'}; %Issues with eye after this point
            case 'PC043'; dateRange = {'2019-03-06:2019-03-29'}; 
            case 'PC045'; dateRange = {'2019-07-30:2019-08-26'}; 
            case 'PC046'; dateRange = {'2019-07-24:2019-08-26'}; 
            case 'PC050'; dateRange = {'2019-07-24:2019-09-11'}; 
            case 'PC051'; dateRange = {'2019-07-24:2019-09-13'}; 
            case 'DJ007'; dateRange = {'2018-05-11:2018-05-19'; '2018-05-11:2018-06-19'}; %Intermediate days without conflicts are excluded
            case 'DJ008'; dateRange = {'2018-05-11:2018-05-30'};
            case 'DJ010'; dateRange = {'2018-06-09:2018-06-18'};
            otherwise, dateRange = [];
        end
        
    case 'aud5'
        switch subjectID{1}
            case 'PC013'; dateRange = {'2017-10-18:2017-11-04'}; %Regular trials within inactivation
            otherwise, dateRange = [];
        end
    case 'uniscan'
        switch subjectID{1}
            case 'PC027'; dateRange = {'2018-02-05:2018-03-21'}; %Power was only 1.5mW
            case 'PC029'; dateRange = {'2018-06-12:2018-07-17'};
            case 'DJ008'; dateRange = {'2018-06-12:2018-07-17'};
            case 'DJ006'; dateRange = {'2018-08-06:2018-09-16'};
            case 'DJ007'; dateRange = {'2018-08-06:2018-10-18'};
            otherwise, dateRange = []; 
        end
        
    case 'biscan'
        switch subjectID{1}
            case 'PC010'; dateRange = {'2017-06-20:2017-07-05'}; %Power was only 1.5mW
            case 'PC012'; dateRange = {'2017-06-21:2017-07-05'};
            case 'PC013'; dateRange = {'2017-06-20:2017-07-06'};
            otherwise, dateRange = []; 
        end
        
    case 'm2ephysgood'
        switch subjectID{1}
            case 'PC045'; dateRange = {'2019-08-27:2019-09-06'};
            case 'PC046'; dateRange = {'2019-08-27:2019-09-06'};
            case 'PC048'; dateRange = {'2019-09-17:2019-09-26'};
            case 'PC050'; dateRange = {'2019-09-17:2019-09-26'};
%             case 'PC032'; dateRange = {'2019-04-03:2019-04-14'};
            otherwise, dateRange = []; 
        end
        
    case 'm2ephys'
        switch subjectID{1}
            case 'DJ007'; dateRange = {'2018-11-28:2018-12-02'}; %Power was only 1.5mW
            case 'PC029'; dateRange = {'2018-10-18:2018-10-19'};
            case 'PC032'; dateRange = {'2019-04-03:2019-04-14'};
            case 'PC033'; dateRange = {'2019-03-27:2019-03-31'};
            case 'PC030'; dateRange = {'2019-05-07:2019-05-13'};
            case 'PC043'; dateRange = {'2019-05-07:2019-05-16'};
            otherwise, dateRange = []; 
        end
        
    case 'presurg'
        switch subjectID{1}
            case 'DJ007'; dateRange = {'2018-11-13:2018-11-27'}; %Power was only 1.5mW
            case 'PC029'; dateRange = {'2018-10-03:2018-10-17'};
            case 'PC032'; dateRange = {'2019-03-13:2019-04-02'};
            case 'PC033'; dateRange = {'2019-03-12:2019-03-26'};
            case 'PC030'; dateRange = {'2019-04-21:2019-05-06'};
            case 'PC043'; dateRange = {'2019-04-21:2019-05-06'};
            otherwise, dateRange = []; 
        end
        
    otherwise, dateRange = dataTag;
end
end