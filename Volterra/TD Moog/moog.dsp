import("stdfaust.lib");
// Constants
// Sampling rate
fs = ma.SR;
ts = 1/fs;
// Cutoff frequency
fc = hslider("Cutoff", 600, 100, 5000, 0.1);
// Normalized pulsation
nu = fc*2*ma.PI*ts;

// One-layer filter
// Linear kernel
F1 = *(nu/(1+nu)): +~*(1/(1+nu));
// Third order kernel
t3 = -1/3;
F3 = _<: F1 ,( _ <: (F1 : ^(3) : *(-1)), ^(3) : + : *(t3) : F1) : +;

// module used each time 
M(u1, u3) = F1(u1), F1(u1^3-F1(u1)^3+u3);
M2 = ((_<:(F1, (_^(3), F1^3:>-))), _):_,((_,_):>+:F1);
M3 = (((_<:(((F1<:(_,^(3)))),^(3))):_,((_,_):>-)),_):_,((_,_):>+:F1);
process4 = _,((_,_):>+);

S(u1, u3) = u1 + t3*u3;
S2 = (_,*(t3)):>+;

// four layers filter
H = _,0:M:M:M:M:S;

phaser(f) = ((f/fs):(+:fmod(_,1))~_);

sawtooth(f) = phaser(f) - 0.5;

square(f) = (sawtooth(f) > 0) - 0.5;

drive = hslider("Drive", 1, 1, 100, 0.1);
freq = hslider("Frequency", 150, 50, 5000, 0.1);
delta = hslider("Delta", 1, 0, 10, 0.01);
gain = hslider("Gain", 0, 0, 10, 0.01);

signal = (square(freq)+square(freq+delta)) / 2;

moog = gain * drive * signal: H / drive;

process = moog;
