function set1 = set1_cols

% Colourscheme based on the ColorBrewer Set 1

% http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=8
set1 = { [228,26,28], [55,126,184], [77,175,74], [152,78,163],...
              [255,127,0], [247,129,191], [255,255,51],[166,86,40]};
          
          %K=12 option:
set1 = {[141,211,199],[188,128,189],[190,186,218],[251,128,114], ...
            [128,177,211], [253,180,98],[179,222,105],[252,205,229], ...
            [217,217,217],[255,255,179],[204,235,197],[255,237,111]};
        

% normalise and convert back to cell array
set1 = cellfun(@(x) (0.9*x)./255, set1,'UniformOutput', false );

% actually better just to use inbuilt matlab colors:
colors_alt=jet(12);
rng(10);
colors_alt = colors_alt(randperm(12),:);
for k=1:12
    set1{k}=colors_alt(k,:);
end
set1{8} = set1{8} + 0.5*[0,1,1];
set1{7} = [141,211,199]./255;
set1{11} = [190,186,218]./255;
