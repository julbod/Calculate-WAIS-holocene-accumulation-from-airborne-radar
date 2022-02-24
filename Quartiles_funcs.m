%% Calculates quartiles for standard boxplot or to obtain summary stats
% Written by Julien Bodart (UoE) -  04.05.2020
%
%%
function output = Quartiles_funcs(x)
    disp(['#################   Quartile Calculations   #################']);

    output = [];

    % compute the standard deviation
    sigma = std(x,'omitnan');
    disp(['Stand Dev: ',num2str(sigma)]);

    % compute the median
    med = median(x,'omitnan');
    disp(['Median: ' , num2str(med)]);
    
    % compute 25th percentile (first quartile)
    Q1 = median(x(find(x<median(x,'omitnan'))),'omitnan');
    disp(['Quartile 1 - 25th: ',num2str(Q1)]);

    % compute 50th percentile (second quartile)
    Q2 = median(x, 'omitnan');
    disp(['Quartile 2 - median: ',num2str(Q2)]);

    % compute 75th percentile (third quartile)
    Q3 = median(x(find(x>median(x,'omitnan'))), 'omitnan');
    disp(['Quartile 3 - 75th: ',num2str(Q3)]);
    
    % compute min and max
    mn = min (x);
    disp(['Minimum: ',num2str(mn)]);
    
    mx = max (x);
    disp(['Maximum: ',num2str(mx)]);

    % compute Interquartile Range (IQR)
    IQR = Q3-Q1;
    disp(['Interquartile Range: ',num2str(IQR)]);    
    
    % compute Semi Interquartile Deviation (SID)
    % The importance and implication of the SID is that if you 
    % start with the median and go 1 SID unit above it 
    % and 1 SID unit below it, you should (normally) 
    % account for 50% of the data in the original data set
    SID = IQR/2;
    disp(['Semi-interquartile Deviation - add/subtract this to median = Q1/Q3: ',num2str(SID)]);   

    % determine extreme Q1 outliers (e.g., x < Q1 - 3*IQR)
    %iy = find(x<Q1-3*IQR);
    %if length(iy)>0,
    %    outliersQ1 = x(iy);
    %else
    %    outliersQ1 = [];
    %end
    
    %disp(['Index of Outlier 1: ',num2str(iy(1))]);   
    %disp(['Value of Outlier 1: ',num2str(outliersQ1(1))]);
    %disp(['Attention: There might be more - find all indices/values by entering: num2str(iy())']);   

    % determine extreme Q3 outliers (e.g., x > Q1 + 3*IQR)
    %iy = find(x>Q1+3*IQR);
    %if length(iy)>0,
    %    outliersQ3 = x(iy);
    %else
    %    outliersQ3 = [];
    %end
    
    % compute total number of outliers
    %Noutliers = length(outliersQ1)+length(outliersQ3);
    %disp(['Total number of outliers: ',num2str(Noutliers)]);   

    %output = [output; sigma; med; Q1; Q2; Q3; mn; mx; IQR; SID; Noutliers; outliersQ1; outliersQ3];
    output = [output; sigma; med; Q1; Q2; Q3; mn; mx; IQR; SID];

end
