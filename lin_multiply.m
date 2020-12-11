function mult_out = lin_multiply(input)

mult_out = 1;

for i = 1:length(input)
    mult_out = mult_out.*input(i);
end
