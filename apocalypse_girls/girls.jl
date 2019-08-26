using PyPlot;

number_of_pairs = 100:100:100000;

boys = zeros(length(number_of_pairs));

for n = 1:length(number_of_pairs)
	for i = 1 : number_of_pairs[n]
		while rand() < 0.5
			boys[n] += 1;
		end;
	end;
end;
