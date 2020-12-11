function signal = upsample_flag(input, flag)

signal = zeros([1, length(flag)]);
a = 1;
sub_i = 1;

for i = 1:length(flag)
    if flag(i) > 0
        if mod(i,2^flag(i)) ~= 0 && a < length(input)-2^flag(i)
            diff = (input(a+1) - input(a))/2^flag(i);
            signal(i) = input(a) + diff*sub_i;
            sub_i = sub_i+1;
        elseif mod(i, 2^flag(i)) == 0
                signal(i) = input(a);
                a = a+1;
                sub_i = 1;
        end
    else
        signal(i) = input(a);
        a = a+1;
        sub_i = 1;
    end
end

end
