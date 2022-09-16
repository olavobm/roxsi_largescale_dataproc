clear all
close all
clc
dd=dir('*spot*');
%%
for m=1:length(dd)
    fname_out=dd(m).name;
    dir_csv = string(dd(m).name)+'/parsed' + filesep + "*.csv";
    d=dir(dir_csv);
    for i=1:length(d)
        fname=d(i).name(1:end-4);
        if fname~="system"
            dir_fn = string(dd(m).name)+'/parsed'  + filesep + string(d(i).name);
            fprintf('Processing %s\n',dir_fn)
            a=readtable(dir_fn,'ReadVariableNames',true,'VariableNamingRule','preserve');
            a = renamevars(a,["# year"],["year"]);
            if fname=="location"|fname=="displacement"|fname=="sst"
                time=datetime(a.year,a.month,a.day,a.hour-7,a.min,a.sec,a.msec);
            else
                time=datetime(a.year,a.month,a.day,a.hour-7,a.min,a.sec,a.milisec);
            end
            a = addvars(a,time,'Before','year');
            s.(fname) = a;    
        else
        end
    end
    save('-v7.3',fname_out,'s')
end
