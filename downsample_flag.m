function signal = downsample_flag(input, flag)

% signal = zeros([1, length(input)]);

a = 1;

for i = 1:length(input)
    if flag(i) > 0
        if mod(i,2^flag(i)) == 0
            signal(a) = input(i);
            a = a+1;
        end
    else
        signal(a) = input(i);
        a = a+1;
    end
end
