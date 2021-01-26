


line1="1 99789U 09999B   21022.64108674  .00000009  00000-0  59233-6 0 00008"
line2="2 99789 097.4984 077.3576 0009198 294.5224 091.6031 15.08170883000346"


function chksum(line1)
  result = false; c = 0;

  str = [string(line1[k]) for i = 1:length(line1)]

  for k = 1:68
    if str[k] > '0' && str[k] <= '9'
      c = c + str[k] - 48;
    elseif str[k] == '-'
      c = c + 1;
    end
  end
​  if mod(c,10) == (str[69] - 48)
    result = true
  end

  return result
end
