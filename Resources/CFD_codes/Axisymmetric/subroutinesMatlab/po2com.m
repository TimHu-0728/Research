function b = po2com(a)
str = sprintf('%g',a);
b = strrep(str,'.',',');
end