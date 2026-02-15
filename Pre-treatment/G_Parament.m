function [normsty,sumD] = G_Parament(data,V,f)

sumD=0;
for v=1:V
    sumD=sumD+size(data{v},2);
end

if sumD<1e3
    normsty="var";
elseif sumD>1e4
    normsty="norm";
else
    if f==0
        normsty="range";
    else
        normsty="norm";
    end
end

end