using PyPlot;

samples = 1000;

PI = 3.14;

funcs = [s -> sin(2 * PI * s / samples)];#, s -> sin(2 * PI * 10 * s / samples), s -> sin(2 * PI * s / samples), s -> sin(2 * PI * 30 * s / samples + 1 / 8 * PI), s -> cos(2*PI*50*s/samples)];
signal = zeros(length(funcs), samples);

for s = 1:samples
	for c = 1:length(funcs)
		signal[c, s] = funcs[c](s - 1);
	end;
end;

t = 0 : 2*PI/samples : 2*PI - 0.001;

figure();
for c = 1:length(funcs)
	plot(t, signal[c, :]);
end;

figure();
sig = zeros(samples);
for c = 1:length(funcs)
	sig .+= signal[c, :];
end; 

#sig = signal[1, :];# + signal[2, :] + signal[3, :] + signal[4, :] + signal[5, :];
plot(t, sig);

X = complex(zeros(samples));


for r = 1:samples
	re = 0;
	imag = 0;
	for c = 1:samples
		re += sig[c]*cos(2 * PI * (r - 1) * (c - 1) / samples);
		imag -= sig[c]*sin(2 * PI * (r - 1) * (c - 1) / samples);
	end;
	X[r] = re/samples + imag*im/samples;
end;


XX = zeros(length(X), samples);

for i = 1 : length(X)
	for j = 1 : samples
		XX[i, j] = sqrt(real(X[i])^2 + imag(X[i])^2) * sin((i - 1) * t[j] + atan(imag(X[i]), real(X[i])));
	end;
end;

dft_signal = zeros(samples);
for i = 1 : length(X)
	dft_signal .+= XX[i, :];
end;

#for i = 1 : samples
#	akku = 0;
#	for j = 1 : length(X)
#		akku += sqrt(real(X[j])*real(X[j]) + imag(X[j])*imag(X[j])) * sin((j - 1) * t[i] + atan(imag(X[j]), real(X[j])));
#	end;
#	XX[i] = akku;
#end;

P = real(X).^2 + imag(X).^2;

f = x -> (x * (1/(t[2]-t[1]))/samples);

figure();
plot(f.(0:999), P);
figure();
plot(t, dft_signal);