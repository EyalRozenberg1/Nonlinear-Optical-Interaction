function [figNum] = PrintExpResults(Print, figNum, date, f, Pulse_wdt)
% Printing experimental results

set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultAxesFontName','times')
set(0,'DefaultTextFontSize',18)
set(0,'DefaultTextFontName','times')
set(gcf,'color','white')

T       = 0;
P2wVsT  = 0;
Pw_1    = 0;
Pw_vec  = 0;

Attenuation = 1; % energy attenuation due to pusle gaussian shape
                 % not in use in case it is covered in the simulation (not in use -> set to 1)

if(strcmp(date,'26_10_17')) % Set gradient and acanning over Tm
% Notes: - Using Lab gradient oven (the one I have created using the new const temp oven) for Tm measurments
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
% for constant temperature of Ts=Te=146.2, Tm=~150.5 --> P2w=14.38
Pw_1      = 50e-3; %[W]
Ts        = (151.2-10:0.1:151.2+10) + offset; % Celsius Degrees - start of crystal
Te        = (141.2-10:0.1:141.2+10) + offset; % Celsius Degrees - end of crystal
Tm_diff   = (Ts+Te)/2;

Tm        = [140   140   140.5 140.5 140.5 141   141   141   141   141   ...
             141   141   141   141.5 142   142   142   142   142   142   ...
             142   142   142.5 142.5 143   143   143   143   143   143   ...
             143   143   143.5 143.5 144   144   144   144   144   144   ...
             144   144   144   144.5 144.5 144.5 145   145   145   145   ...
             145   145   145   145.5 146   146   146   146   146   146   ...
             146   146   146   146   146.5 147   147   147   147   147   ...
             147   147   147   147.5 147.5 148   148   148   148   148   ...
             148   148   148   148   148   148.5 149   149   149   149   ...
             149   149   149   149   149   149.5 149.5 150   150   150   ...
             150   150   150   150   150   150   150.5 151   151   151   ...
             151   151   151   151   151   151.5 151.5 152   152   152   ...
             152   152   152   152   152   152   152.5 153   153   153   ...
             153   153   153   153   153   153.5 154   154   154   154   ...
             154   154   154   154   154   154   154.5 154.5 155   155   ...
             155   155   155   155   155   155   155   156   156   156   ...
             157   157   157   157   157   157   157   158   158   158   ...
             158   158   158   158   158   158   158.5 158.5 159   159   ...
             159   159   159   159   159   159.5 159.5 159.5 159.5 159.5 ...
             159.5 160   160   160   160   160   160   160   160   160   160]; % Celsius Degrees - middle of crystal (measured with Escort)

P2wVsT    = [0.47  0.45  0.45  0.47  0.45  0.41  0.39  0.36  0.35  0.38  ...
             0.48  0.51  0.53  0.50  0.48  0.44  0.40  0.42  0.44  0.52  ...
             0.59  0.61  0.62  0.58  0.51  0.48  0.46  0.49  0.58  0.60  ...
             0.65  0.70  0.71  0.68  0.58  0.56  0.56  0.59  0.64  0.73  ...
             0.75  0.79  0.85  0.78  0.74  0.72  0.67  0.68  0.72  0.84  ...
             0.92  0.98  1.03  1.00  0.94  0.86  0.83  0.82  0.89  0.99  ...
             1.05  1.12  1.20  1.25  1.24  1.17  1.13  1.10  1.08  1.11  ...
             1.14  1.12  1.18  1.62  1.67  1.64  1.65  1.65  1.66  1.57  ...
             1.53  1.55  1.56  1.74  2.03  2.40  2.50  2.66  2.67  2.75  ...
             2.77  2.80  3.06  3.32  3.93  4.24  5.30  6.70  7.00  7.27  ...
             8.87  9.00  10.92 12.34 13.01 13.93 14.00 13.97 13.40 13.00 ...
             12.82 13.16 12.86 12.22 12.01 11.64 10.96 10.73 10.77 11.08 ...
             10.84 11.03 10.77 10.79 11.17 11.11 11.08 10.88 10.71 10.30 ...
             10.20 9.99  9.70  9.67  9.57  9.03  8.71  8.50  8.44  8.41  ...
             8.28  8.02  7.80  7.70  7.17  7.00  6.60  6.61  6.85  6.82  ...
             6.74  6.44  6.52  6.11  5.64  5.56  5.54  5.44  5.42  5.66  ...
             5.45  4.45  4.42  4.48  4.50  4.53  4.59  4.55  4.13  3.77  ...
             3.75  3.61  3.50  3.65  3.77  3.80  3.83  3.77  3.67  3.44  ...
             2.74  2.48  2.30  2.34  2.52  2.60  2.68  2.68  2.48  2.38  ...
             2.24  2.14  2.01  1.94  2.06  2.05  2.13  2.19  2.25  2.20  2.12]*1e-3/(transmission*Attenuation); %[W] ?? samples

end                                                   
                 
                 
if(strcmp(date,'24_10_17')) % Set gradient and acanning over Tm
% Notes: - Using Lab gradient oven (the one I have created using the new const temp oven) for Tm measurments
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
% for constant temperature of Ts=Te=147.2, Tm=~150.5 --> P2w=14.26
Pw_1      = 50e-3; %[W]
Ts        = (152.2-10:0.1:152.2+10) + offset; % Celsius Degrees - start of crystal
Te        = (142.2-10:0.1:142.2+10) + offset; % Celsius Degrees - end of crystal
Tm_diff   = (Ts+Te)/2;

Tm        = [140  140.5  141   141   141   141  141  141  141  141  ...
             141  141    141.5 142   142   142  142  142  142  142  ...
             142  142    142   142.5 142.5 143  143  143  143  143  ...
             143  143    143.5 143.5 144   144  144  144  144  144  ...
             144  144    144   144.5 145   145  145  145  145  145  ...
             145  145    145   146   146]; % Celsius Degrees - middle of crystal (measured with Escort)

P2wVsT    = [0.46  0.46  0.42  0.41  0.38  0.37  0.39  0.41  0.45  0.50  ...
             0.54  0.53  0.52  0.48  0.42  0.41  0.44  0.47  0.53  0.54  ...
             0.59  0.63  0.60  0.56  0.53  0.48  0.49  0.50  0.55  0.68  ...
             0.72  0.72  0.71  0.66  0.62  0.57  0.57  0.60  0.75  0.82  ...
             0.85  0.85  0.83  0.73  0.69  0.69  0.81  0.95  0.96  1.03  ...
             1.04  1.01  0.97  0.87  0.86  1.09  1.15  1.22  ]*1e-3/(transmission*Attenuation); %[W] ?? samples

end                                  
                 
if(strcmp(date,'22_10_17')) % Set gradient and acanning over Tm
% Notes: - Using Lab gradient oven (the one I have created using the new const temp oven) for Tm measurments
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
% for constant temperature of Ts=Te=148.6, Tm=~150 --> P2w=14.24
Pw_1      = 50e-3; %[W]
Ts        = (151.1-5:0.1:151.1+5) + offset; % Celsius Degrees - start of crystal
Te        = (146.1-5:0.1:146.1+5) + offset; % Celsius Degrees - end of crystal
Tm_diff   = (Ts+Te)/2;
Tm        = [147   147   147   147   147   147   147   147   148   148   ...
             148   148   148   148   148   148   148   148   148   148.5 ...
             149   149   149   149   149   149   149   149   149   149.5 ...
             149.5 149.5 150   150   150   150   150   150   150   150   ...
             150   150   150   150.5 150.5 151   151   151   151   151   ...
             151   151   151   151   151   151   151.5 152   152   152   ...
             152   152   152   152   152   152   152   152.5 152.5 152.5 ...
             153   153   153   153   153   153   153   153   153   153.5 ...
             153.5 153.5 154   154   154   154   154   154   154   154   ...
             154   154.5 154.5 155   155   155   155   155   155   155   155]; % Celsius Degrees - middle of crystal (measured with Escort)

P2wVsT    = [1.52  1.58  1.66  1.80  1.85  1.87  1.89  1.85  1.75  1.77  ...
             2.01  1.78  1.88  2.21  2.37  2.73  2.93  3.00  3.12  3.17  ...
             3.23  3.20  3.31  3.60  3.71  4.22  4.97  5.15  5.68  7.09  ...
             8.03  8.74  8.88  10.15 11.44 12.47 12.79 13.12 13.60 13.65 ...
             13.74 12.92 12.85 12.82 12.79 12.62 11.95 11.75 11.18 10.95 ...
             10.68 10.40 10.29 10.90 10.69 10.80 10.87 10.89 10.75 10.22 ...
             10.20 10.13 10.00 9.97  9.79  9.61  9.45  9.27  8.75  8.81  ...
             7.97  8.05  7.81  7.80  7.85  7.86  7.87  7.83  7.72  7.52  ...
             7.08  6.63  6.49  6.27  6.23  6.27  6.30  6.31  6.31  6.32  ...
             6.20  5.97  5.67  5.30  5.16  5.10  5.09  5.04  5.10  5.17  5.20]*1e-3/(transmission*Attenuation); %[W] ?? samples

end                 
                 
if(strcmp(date,'17_10_17')) % Gradient diff
% Notes: - Using Lab gradient oven (the one I have created using the new const temp oven) for gradient measurments
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 50e-3; %[W]
Ts        = (153.3:-0.1:144.1) + offset; % Celsius Degrees - start of crystal
Tm        = 150 + offset; % Celsius Degrees - middle of crystal (measured with Escort)
Te        = (139.3:0.1:148.5)  + offset; % Celsius Degrees - end of crystal

P2wVsT    = [13.11  13.32  13.53  13.54  13.54  13.46  13.44  13.53  13.55  13.56  ...
             13.57  13.57  13.61  13.60  13.44  13.57  13.57  13.58  13.59  13.64  ...
             13.58  13.58  13.36  13.19  13.17  13.21  13.17  13.15  13.19  13.19  ...
             13.18  13.09  13.07  13.14  13.10  13.07  13.08  13.10  13.12  13.17  ...
             13.15  13.13  13.15  13.39  13.32  13.45  13.50  13.46  13.43  13.35  ...
             13.15  13.34  13.41  13.50  13.56  13.60  13.66  13.80  13.76  13.73  ...
             13.80  13.65  13.80  13.73  13.81  13.81  13.70  13.75  13.78  13.82  ...
             13.84  13.85  13.55  13.52  13.41  13.60  13.63  13.58  13.58  13.68  ...
             13.71  13.53  13.41  13.35  13.34  13.36  13.38  13.03  13.16  13.23  ...
             13.45  15.53  13.47]*1e-3/(transmission*Attenuation); %[W] ?? samples

end

if(strcmp(date,'10_10_17')) % Scanning temperature from low to high (as usual)
% Notes: - Using Lab gradient oven (the one I have created using the new const temp oven) for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = -1.5+3.8; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 10e-3; %[W]
T        = (142:0.1:152) + offset; % Celsius Degrees

P2wVsT   = [0.11  0.10  0.11  0.10  0.09  0.10  0.12  0.14  0.14  0.14  ...
            0.13  0.13  0.14  0.14  0.14  0.13  0.13  0.18  0.20  0.21  ...
            0.22  0.23  0.23  0.23  0.25  0.29  0.40  0.46  0.49  0.49  ...
            0.54  0.64  0.67  0.74  0.85  0.95  1.36  1.55  2.52  2.95  ...
            3.35  4.21  4.56  4.32  3.80  2.96  2.48  2.25  2.51  2.62  ...
            2.47  2.39  2.27  2.03  1.77  1.66  1.61  1.62  1.60  1.57  ...
            1.52  1.47  1.37  1.20  1.16  1.13  1.13  1.13  1.13  1.12  ...
            1.07  1.04  0.99  0.85  0.75  0.74  0.76  0.78  0.81  0.80  ...
            0.75  0.70  0.58  0.49  0.44  0.44  0.50  0.54  0.57  0.59  ...
            0.59  0.50  0.43  0.30  0.28  0.27  0.35  0.41  0.43  0.44  0.44]*1e-3/(transmission*Attenuation); %[W] ?? samples

end


if(strcmp(date,'01_10_17')) % Scanning temperature from low to high (as usual)
% Notes: - Using Lab gradient oven (the one I have created using the gradient with the teflon chunks) for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
%        - ** After temperature ~147 the laser was very unstable !! **
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = -1.5; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 10e-3; %[W]
T        = (140:0.1:149) + offset; % Celsius Degrees


P2wVsT   = [0.13  0.14  0.17  0.23  0.32  0.43  0.52  0.56  0.60  0.56  ...
            0.50  0.46  0.39  0.36  0.56  0.82  1.04  1.17  1.22  1.20  ...
            1.17  1.11  1.03  0.79  0.65  0.59  0.73  1.03  1.31  1.61  ...
            1.88  2.20  2.17  1.90  1.04  0.86  0.62  0.71  1.00  1.20  ...
            1.36  1.73  2.10  2.48  2.65  2.45  1.91  1.51  1.06  0.85  ...
            0.60  0.53  0.51  0.75  1.25  1.69  1.94  2.00  2.14  2.00  ...
            1.67  1.32  1.04  0.87  0.82  0.85  0.94  1.02  1.16  1.31  ...
            1.37  1.51  1.59  1.62  1.74  1.70  1.68  1.70  1.51  1.29  ...
            0.92  0.52  0.40  0.31  0.28  0.34  0.47  0.53  0.67  0.77  0.91]*1e-3/(transmission*Attenuation); %[W] ?? samples

end


if(strcmp(date,'26_09_17')) % Scanning temperature from low to high (as usual)
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = -1.5; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 50e-3; %[W]
T        = (145:0.1:156) + offset; % Celsius Degrees


P2wVsT   = [1.30  1.41  1.53  1.63  1.68  1.67  1.62  1.57  1.535 1.56  ...
            1.64  1.79  1.96  2.12  2.26  2.30  2.31  2.285 2.25  2.24  ...
            2.27  2.40  2.615 2.86  3.13  3.37  3.57  3.68  3.74  3.72  ...
            3.69  3.75  3.98  4.48  5.18  5.95  6.41  6.67  6.33  6.90  ...
            8.13  8.76  9.85  10.67 10.16 11.19 11.60 11.81 11.90 12.31 ...
            12.86 11.95 12.19 12.51 11.89 11.56 11.44 11.34 11.20 10.87 ...
            10.41 9.95  9.52  9.18  8.92  8.71  8.51  8.31  8.11  7.88  ...
            7.75  7.53  7.30  7.10  6.85  6.68  6.52  6.37  6.20  6.04  ...
            5.87  5.71  5.54  5.31  5.18  5.10  4.98  4.88  4.79  4.67  ...
            4.55  4.45  4.33  4.24  4.13  4.02  3.92  3.82  3.73  3.67  ...
            3.63  3.58  3.53  3.46  3.36  3.26  3.16  3.08  3.04  3.02  3.02]*1e-3/(transmission*Attenuation); %[W] 111 samples

end


if(strcmp(date,'24_09_17')) % Scanning temperature from low to high (as usual)
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = -1.5; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 30e-3; %[W]
T        = (145:0.1:155) + offset; % Celsius Degrees


P2wVsT   = [0.68  0.68  0.65  0.65  0.65  0.69  0.74  0.79  0.84  0.90  ...
            0.92  0.91  0.89  0.88  0.88  0.90  0.94  1.02  1.11  1.21  ...
            1.285 1.32  1.33  1.31  1.29  1.265 1.34  1.44  1.615 1.82  ...
            2.01  2.185 2.30  2.36  2.355 2.33  2.355 2.56  3.08  3.76  ...
            4.45  4.98  5.35  4.80  5.81  6.805 7.05  7.305 7.90  9.21  ...
            10.32 9.06  8.68  7.98  7.14  7.00  6.99  6.94  6.82  6.59  ...
            6.32  5.94  5.57  5.25  4.99  4.83  4.70  4.63  4.52  4.37  ...
            4.13  3.93  3.74  3.61  3.53  3.50  3.42  3.34  3.26  3.15  ...
            2.98  2.89  2.82  2.75  2.72  2.71  2.64  2.54  2.43  2.30  ...
            2.26  2.18  2.14  2.17  2.11  2.09  1.99  1.89  1.82  1.73  1.63]*1e-3/(transmission*Attenuation); %[W] 101 samples

end


if(strcmp(date,'19_09_17')) % Scanning temperature from low to high (as usual)
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = -1.5; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 10e-3; %[W]
T        = (145:0.1:155) + offset; % Celsius Degrees


P2wVsT   = [0.08  0.11  0.13  0.15  0.16  0.15  0.13  0.10  0.09  0.09  ...
            0.11  0.15  0.19  0.21  0.215 0.20  0.18  0.15  0.14  0.15  ...
            0.18  0.22  0.27  0.30  0.31  0.30  0.28  0.27  0.27  0.29  ...
            0.33  0.39  0.45  0.51  0.56  0.60  0.63  0.65  0.67  0.71  ...
            0.79  0.935 1.165 1.45  1.75  2.03  2.48  2.85  3.48  4.10  ...
            4.92  3.875 2.95  2.755 2.605 2.44  2.30  2.17  2.04  1.90  ...
            1.71  1.55  1.43  1.36  1.30  1.25  1.22  1.18  1.12  1.02  ...
            0.89  0.79  0.73  0.69  0.68  0.70  0.71  0.70  0.64  0.56  ...
            0.46  0.38  0.35  0.35  0.39  0.44  0.46  0.44  0.39  0.31  ...
            0.25  0.20  0.18  0.21  0.26  0.30  0.31  0.29  0.24  0.18  0.12]*1e-3/(transmission*Attenuation); %[W] 101 samples
%         0.79  0.935 1.165 1.45  1.75  2.03  2.64  ...]

end

if(strcmp(date,'15_09_17')) % Scanning temperature from low to high (as usual)
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-PD-3W 2w-PD-UV
%        - Waist  of intensity = 49e-6 [m]
transmission = 0.7; % 
%        - laser type: Luce
offset   = -1.1; % Celsius Degrees Ghosh Deltak: -0.8
% P2w Vs Temperature
Pw_1     = 25e-3; %[W]
T        = (148:0.1:153) + offset; % Celsius Degrees


P2wVsT   = [0.048  0.051  0.060  0.071  0.082  0.088  0.088  0.083  0.075  0.081  ...
            0.117  0.206  0.347  0.545  0.813  1.147  1.474  1.797  2.012  2.114  ...
            2.136  2.044  1.842  1.625  1.296  1.043  0.835  0.672  0.550  0.428  ...
            0.344  0.285  0.239  0.195  0.162  0.154  0.145  0.143  0.147  0.143  ...
            0.137  0.123  0.107  0.091  0.075  0.067  0.064  0.066  0.069  0.073  0.073]*1e-3/(transmission*Attenuation); %[W] 51 samples

%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Filter PD on
Temp    =  150;
Pw_vec  =  [2      6.64   9.35   11.59  13.21  14.53  16.34  17.53  18.88  20.42]*1e-3; %[W] 10 samples
P2wVsPw =  [0.048  0.198  0.356  0.518  0.661  0.783  0.968  1.109  1.268  1.426]*1e-3/(transmission*Attenuation); %[W]
end


if(strcmp(date,'13_09_17')) % Scanning temperature from high to low
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has 0.7 transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity = 49e-6 [m]
transmission = 0.7; % 
%        - laser type: Luce
offset   = -1.2; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 24.9e-3; %[W]
T        = (147:0.1:155) + offset; % Celsius Degrees


P2wVsT   = [0.047  0.045  0.045  0.046  0.049  0.053  0.057  0.058  0.058  0.055  ...
            0.052  0.050  0.053  0.061  0.072  0.085  0.094  0.092  0.088  0.081  ...
            0.088  0.124  0.218  0.361  0.573  0.845  1.173  1.501  1.704  1.902  ...
            1.987  1.919  1.847  1.680  1.452  1.243  1.067  0.930  0.857  0.766  ...
            0.694  0.650  0.560  0.485  0.410  0.340  0.265  0.211  0.168  0.138  ...
            0.121  0.107  0.095  0.084  0.075  0.068  0.065  0.065  0.069  0.074  ...
            0.078  0.079  0.077  0.072  0.067  0.063  0.060  0.059  0.058  0.057  ...
            0.055  0.052  0.049  0.046  0.045  0.045  0.046  0.049  0.050  0.050  0.049]*1e-3/(transmission*Attenuation); %[W] 81 samples

end

if(strcmp(date,'06_09_17_2'))% quick measurment to avoid Pw degradation (as happened in 06_09_17)
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity  = 49e-6 [m]
transmission = 0.7; % 
%        - laser type: Luce
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 30e-3; % to 29.9e-6 %[W]
T        = (147:0.1:151) + offset; % Celsius Degrees


P2wVsT   = [0.061  0.055  0.051  0.059  0.075  0.097  0.123  0.132  0.130  0.110  ...
            0.094  0.091  0.148  0.282  0.536  0.930  1.392  1.870  2.416  2.97   ...
            3.10   3.30   3.14   3.04   2.581  2.145  1.675  1.359  0.970  0.796  ...
            0.575  0.463  0.366  0.320  0.265  0.226  0.191  0.175  0.166  0.176  0.179]*1e-3/(transmission*Attenuation); %[W] 41 samples

%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Filter PD on
Temp    =  149.1;
Pw_vec  =  [5.94  9.18  11.3  13.79 15.44 16.83 18.50 20.4  21.60 22.45 24.0]*1e-3; %[W] 11 samples
P2wVsPw =  [0.192 0.405 0.575 0.823 1.001 1.160 1.388 1.669 1.838 2.040 2.20]*1e-3/(transmission*Attenuation); %[W]
end


if(strcmp(date,'06_09_17'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity  = 49e-6 [m]
transmission = 0.7; % 
%        - laser type: Luce
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 28.78e-3;  % or less due to degradation from 30 e-3 to ~26e-3 duering measurment%[W]
T        = (146:0.1:152) + offset; % Celsius Degrees


P2wVsT   = [0.051  0.0482 0.0451 0.0443 0.0471 0.0525 0.0591 0.0645 0.0654 0.0625 ...
            0.0559 0.0505 0.0496 0.0575 0.0721 0.092  0.109  0.114  0.111  0.095  ...
            0.081  0.088  0.147  0.275  0.510  0.845  1.203  1.625  2.035  2.430  ...
            2.559  2.600  2.563  2.323  2.065  1.655  1.261  0.984  0.744  0.549]*1e-3/(transmission*Attenuation); %[W] 40 samples

end


if(strcmp(date,'05_09_17'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity  = 49e-6 [m]
transmission = 0.7; % 
%        - laser type: Luce
offset   = -0.3; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 29.20e-3; % or 30e-3;  %[W]
T        = (147:0.1:151) + offset; % Celsius Degrees


P2wVsT   = [0.056  0.050  0.049  0.056  0.073  0.094  0.113  0.124  0.118  0.103  ...
            0.085  0.092  0.147  0.282  0.535  0.910  1.370  1.888  2.350  2.800  ...
            3.10   3.25   3.20   3.00   2.650  2.130  1.700  1.300  0.960  0.690  ...
            0.505  0.405  0.340  0.285  0.245  0.213  0.183  0.177  0.173  0.165  0.174]*1e-3/(transmission*Attenuation); %[W] ?? samples

%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Filter PD on
Temp    =  149.1;


Pw_vec  =  [3.845  5.85   8.76   11.02  13.34 14.62  16.21  17.88  19.40 20.9  ...
            22.15  23.09  24.24  25.50  26.2  27.44  28.30  29.50  30.1  31.7]*1e-3; %[W]

P2wVsPw =  [0.107 0.197  0.394  0.589  0.825  0.970  1.160  1.400  1.610  1.86  ...
            2.060 2.2    2.4    2.670  2.81   3      3.2    3.4    3.61   3.83]*1e-3/(transmission*Attenuation); %[W]
end


if(strcmp(date,'29_08_17_3'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity  = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las

%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Filter PD on
Temp    =  149.5; % Was set according to 23_08_17_1 as it is undepleted measurment
Pw_vec  =  [0.73  0.77  0.85  0.87  0.91  0.94  1.00  1.04  1.09  1.15  1.23  1.34]*1e-3; %[W] 12 Samples     1.46  1.56  1.96  2.46  
P2wVsPw =  [0.18  0.20  0.23  0.24  0.25  0.27  0.29  0.31  0.34  0.36  0.40  0.45]*1e-3/(transmission*Attenuation); %[W]               0.51  0.56  0.76  1.02  
end


if(strcmp(date,'29_08_17_2'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity  = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 10.03e-3;  %[W]
T        = (147.5:0.1:151.5) + offset; % Celsius Degrees


P2wVsT   = [0.37  0.44  0.50  0.55  0.575 ...
            0.58  0.585 0.59  0.65  0.76  0.95  1.18  1.42  1.65  1.82  ...
            1.99  2.26  2.74  3.50  4.45  4.70  3.42  2.36  2.45  2.525 ...
            2.53  2.45  2.30  2.10  1.88  1.655 1.47  1.30  1.24  1.23  ...
            1.23  1.24  1.22  1.14  1.01  0.87]*1e-3/(transmission*Attenuation); %[W] 41 samples

%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Filter PD on
Temp    =  149.5; % Was set according to 23_08_17_1 as it is undepleted measurment
Pw_vec  =  [4.385 4.55  4.68  4.77  4.93  5.225 5.41  5.65  5.95  6.22  ...
            6.39  6.63  6.91  7.33]*1e-3; %[W] 14 Samples
P2wVsPw =  [2.02  2.11  2.17  2.21  2.29  2.43  2.515 2.62  2.76  2.87  ...
            2.92  3.03  3.14  3.31]*1e-3/(transmission*Attenuation); %[W]
end


if(strcmp(date,'29_08_17'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity  = 120e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0.7; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 50e-3;  %[W]
T        = (147:0.1:152) + offset; % Celsius Degrees

P2wVsT   = [3.49  3.50  3.44  3.45  3.70  4.04  4.53  5.20  5.62  6.05  ...
            6.35  6.20  5.99  7.40  8.50  9.80  10.29 10.30 10.44 11.30 ...
            12.20 12.77 13.05 13.38 13.24 12.91 13.38 13.63 13.23 12.91 ...
            12.55 12.18 11.93 11.73 11.53 11.32 11.10 10.87 10.60 10.23 ...
            9.80  9.40  9.04  8.72  8.40  8.15  7.93  7.75  7.55  7.35  7.12]*1e-3/(transmission*Attenuation); %[W] 51 samples

%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Filter PD on
Temp    =  149.7; % Was set according to 23_08_17_1 as it is undepleted measurment
Pw_vec  =  [1.32  2.01  2.81  4.45  5.93  8.00  8.51  9.44  10.02]*1e-3; %[W] 9 Samples
P2wVsPw =  [0.34  0.56  0.77  1.01  1.17  1.49  1.64  2.03  2.36 ]*1e-3/(transmission*Attenuation); %[W]

% FIRST MEASURMENT WAS UNDER pd'limit which is 200uW therefor 0.16 is deleted
% Pw_vec  =  [0.79  1.32  2.01  2.81  4.45  5.93  8.00  8.51  9.44  10.02]*1e-3; %[W] 10 Samples
% P2wVsPw =  [0.16  0.34  0.56  0.77  1.01  1.17  1.49  1.64  2.03  2.36 ]*1e-3/(transmission*Attenuation); %[W]
end



if(strcmp(date,'23_08_17_3'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-Thermal 2w-PD-3W
%        - Waist  of intensity  = 95e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 50e-3;  %[W]
T        = (147:0.1:152) + offset; % Celsius Degrees

P2wVsT   = [2.63  3.02  3.44  3.85  4.11  4.20  4.14  3.97  3.84  3.94  ...
            4.40  5.29  6.42  7.20  7.65  7.95  7.86  7.84  8.23  9.57  ...
            10.60 10.97 11.40 11.70 11.63 11.40 11.05 11.26 11.20 10.87 ...
            11.17 10.63 10.10 9.35  9.05  8.73  8.43  8.19  7.99  7.71  ...
            7.50  6.05  5.90  5.80  5.66  5.50  5.33  5.18  5.06  4.97  4.85]*1e-3/(transmission*Attenuation); %[W] 51 samples

%        - Power Measurments: w-PD-3W 2w-PD-3W
Temp    =  149.7; % Was set according to 23_08_17_1 as it is undepleted measurment
Pw_vec  =  [3.78  5.02  7.14  9.09  10.69]*1e-3; %[W] 5 Samples
P2wVsPw =  [1.46  1.87  2.21  2.99  3.90]*1e-3/(transmission*Attenuation); %[W]
end


if(strcmp(date,'23_08_17_2'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity  = 95e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 10e-3;  %[W]
T        = (147:0.1:152) + offset; % Celsius Degrees

P2wVsT   = [0.28  0.30  0.32  0.34  0.35  0.36  0.38  0.40  0.43  0.46  ...
            0.50  0.56  0.64  0.73  0.82  0.91  0.99  1.07  1.17  1.32  ...
            1.54  1.78  2.00  2.31  2.79  3.37  3.94  4.03  3.20  2.71  ...
            2.65  2.45  2.27  2.12  1.97  1.83  1.68  1.54  1.44  1.37  ...
            1.32  1.29  1.26  1.20  1.14  1.06  0.98  0.91  0.84  0.80  0.78]*1e-3/(transmission*Attenuation); %[W] 51 samples

end


if(strcmp(date,'23_08_17_1'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-PD-3W 2w-PD-3W
%        - Waist  of intensity  = 95e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 5.6e-3;  %[W]
T        = (146:0.1:152) + offset; % Celsius Degrees

P2wVsT   = [0.07  0.06  0.05  0.04  0.04  0.05  0.08  0.09  0.11  0.11  ...
            0.10  0.09  0.08  0.08  0.09  0.11  0.14  0.17  0.19  0.21  ...
            0.20  0.20  0.20  0.20  0.23  0.29  0.36  0.44  0.52  0.59  ...
            0.65  0.74  0.85  1.03  1.26  1.54  1.82  2.14  1.83  1.01  ...
            0.94  0.97  0.99  0.97  0.91  0.81  0.69  0.57  0.46  0.39  ...
            0.35  0.35  0.36  0.39  0.40  0.39  0.34  0.28  0.21  0.16  0.14]*1e-3/(transmission*Attenuation); %[W] 61 samples

end


if(strcmp(date,'22_08_17'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm eksma 2A641
%        - FL532-10 laser line filter has ?? transmittance
%        - Power Measurments: w-thermal 2w-thermal
%        - Waist  of intensity  = 95e-6 [m]
transmission = 0.7; % 
%        - laser type: Alpha Las
offset   = 0; % Celsius Degrees
% P2w Vs Temperature
Pw_1     = 50e-3;  %[W] ~51e-3 [W]
T        = (148.5:0.1:151.5) + offset; % Celsius Degrees

P2wVsT   = [3  4  5  7  8  9  9  8  9  10 ...
            12 14 15 16 16 15 15 15 14 13 ...
            12 12 12 11 11 10 10 10 9  9  9]*1e-3/(transmission*Attenuation); %[W] 31 samples

% P2w Vs Pw
Temp    =  149.8 + offset;
Pw_vec  =  [17 22 26 28 31 33 37 40 43 47 51]*1e-3; %[W] 11 Samples
P2wVsPw =  [2  4  5  6  8  9  10 11 12 13 15]*1e-3/(transmission*Attenuation); %[W]
end


if(strcmp(date,'03_07_17'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 4x4x50 mm Gamdan
%        - 532 bandpass filter has 44.5% transmittance
%        - Waist  of intensity  = <45.5e-6 [m] (maybe sqrt(2)*45.5e-6)
transmission = 0.445;
%        - laser type: Luce

% P2w Vs Temperature
Pw_1     = 30e-3;  %[W]
T        = 148.5:0.1:150.5;

P2wVsT   = [0.0426  0.0625  0.139  0.299  0.533  0.933  1.425  1.945  2.494  2.99 ...
            3.34    3.49    3.51   3.28   2.97   2.520  2.065  1.598  1.281  1.015  0.810]*1e-3/(transmission*Attenuation); %[W] 21 samples

% P2w Vs Pw
Temp    =  149.6;
Pw_vec  =  [5     10    15.02 20   25.17 30   34.95 40   45   50 ]*1e-3; %[W] 13 Samples
P2wVsPw =  [0.152 0.505 1.065 1.77 2.6   3.53 4.655 5.77 7.05 8.4]*1e-3/(transmission*Attenuation); %[W]
end

if(strcmp(date,'02_07_17'))
% Notes: - Using Eksma oven for constant temperature
%        - Not the same crystal as 22_05_17
%        - Crystal 5x5x50 mm
%        - 532 bandpass filter has 44.5% transmittance
%        - Error: uncalicrated lambda/2 plate - need to be done again
%        - Waist  of intensity  = <45.5e-6 [m] (maybe sqrt(2)*45.5e-6)
transmission = 0.445;
%        - laser type: Luce

% P2w Vs Temperature
Pw_1     = 30e-3;  %[W]
T        = 148.5:0.1:150.5;

P2wVsT =[0.0484 0.0775 0.1675 0.353 0.630 1.015 1.440 1.865 2.400 2.700 ...
         2.900  2.83   2.700  2.5    2.10 1.750 1.370 1.080 0.900 0.730 0.640]*1e-3/(transmission*Attenuation); %[W]

% P2w Vs Pw
Temp    =  149.6;
Pw_vec  =  [4.92  8.32  12.99 15.65 18.22 22.53 26.23 30.1 35.6 40.2 43.0 46.6 49.8]*1e-3; %[W]
P2wVsPw =  [0.127 0.293 0.636 0.875 1.158 1.650 2.224 2.75 3.85 4.63 5.18 5.97 6.69]*1e-3/(transmission*Attenuation); %[W]     
end

if(strcmp(date,'22_05_17'))
% Notes: - Using Eksma oven for constant temperature
%        - Crystal 5x5x50 mm
%        - 532 bandpass filter has 44.5% transmittance
%        - Error: uncalicrated lambda/2 plate - need to be done again
% transmission = 0.445;
%        - Waist  of intensity  = <45.5e-6 [m] (maybe sqrt(2)*45.5e-6)
transmission = 1;%0.32
%        - laser type: Luce

% P2w Vs Temperature
Pw_1     = 30e-3;  %[W]
T        = 146.5:0.1:152.5;

P2wVsT =[0.0536 0.0518 0.0505 0.0507 0.0530 0.0568 0.0604 0.0630 0.0628 0.0597 ...
         0.0557 0.0534 0.0553 0.0623 0.0740 0.0859 0.0945 0.0940 0.0850 0.0723 ...
         0.0700 0.0963 0.177  0.325  0.537  0.838  1.190  1.560  1.900  2.208  ...
         2.330  2.400  2.320  2.130  1.870  1.610  1.365  1.110  0.930  0.825  ...
         0.740  0.695  0.670  0.635  0.595  0.549  0.500  0.446  0.400  0.350  ...
         0.314  0.289  0.265  0.243  0.230  0.227  0.225  0.222  0.221  0.215  0.212]*1e-3/(transmission*Attenuation); %[W]

% P2w Vs Pw
Temp    =  149.6;
Pw_vec  =  [3.530 4.060 10.04 20.10 30.70 40.50]*1e-3; %[W]
% Pw_vec  =  [3.530 4.060 10.04 20.10 30.70 40.50 50.10 61.20 71.00 80.00 ...
%             90.60 100.8 110.8 120.4 130.0 141.1 150.0 161.0 170.7 180.3 ...
%             190.6 200.4 210.1 219.9 230.3 240.0 249.4 260.6 270.8 280.5 ...
%             289.8 300.0 312.0 320.0 330.0 340.0 350.0 360.0 371.0 382.0 ...
%             390.0 400.0 410.0 420.0 430.0 440.0 449.0 460.0 470.0 481.0 ...
%             491.0 500.0]*1e-3; %[W]

P2wVsPw = [0.088 0.100 0.349 1.160 2.404 4.010]*1e-3/(transmission*Attenuation); %[W]
       
% P2wVsPw = [0.088 0.100 0.349 1.160 2.404 4.010 5.800 7.960 10.13 12.33 ...
%            15.00 17.53 20.15 23.00 25.50 28.75 31.10 34.40 37.50 40.40 ...
%            44.10 47.40 48.50 51.50 54.80 58.00 61.50 65.20 68.90 71.70 ...
%            75.50 78.70 82.50 85.20 87.80 92.30 94.80 98.00 101.3 105.2 ...
%            107.3 111.7 114.1 117.1 121.1 123.0 126.3 129.5 133.1 136.7 ...
%            138.3 140.7]*1e-3/(transmission*Attenuation); %[W]       
end


if(strcmp(date,'09_10_16'))
% P2w Vs Temperature
Pw_1  = 500e-3;  %[W]

T   =   [128.0 128.5 129.0 129.5 130.0 130.2 130.4 130.6 130.8 131 ...
         131.2 131.4 131.6 131.8 132.0 132.2 132.4 132.6 132.8 133 ...
         133.2 133.4 133.6 133.8 134.0 134.2 134.4 134.6 134.8 135 ...
         135.2 135.4 135.6 135.8 136.0 136.2 136.4 136.6 136.8 137 ...
         137.1 137.2 137.3 137.4 137.5 137.6 137.7 137.8 137.9 138 ...
         138.1 138.2 138.3 138.4 138.5 138.6 138.7 138.8 138.9 139 ...
         139.2 139.4 139.6 139.8 140.0 140.2 140.4 140.6 140.8 141 ...
         141.2 141.4 141.6 141.8 142.0 142.2 142.4 142.6 142.8 143 ...
         143.2 143.4 143.6 143.8 144.0 144.2 144.4 144.6 144.8 145 ...
         145.2 145.4 145.6 145.8 146.0 146.2 146.4 146.6 146.8 147];

P2wVsT =[0  1  1  2  0  2  5  6  3  6  ...
         6  8  10 10 6  11 22 24 24 17 ...
         12 9  7  17 25 20 16 18 28 38 ...
         45 51 46 49 55 64 76 85 84 76 ...
         73 74 76 80 86 89 94 94 93 88 ...
         80 79 68 65 65 64 61 57 59 59 ...
         53 47 37 27 20 12 7  6  6  5  ...
         7  10 14 15 15 16 16 15 17 18 ...
         22 25 26 22 18 11 8  5  5  5  ...
         6  7  9  12 13 13 10 8  8  8]*10^-3; %[W]

% P2w Vs Pw
Temp =  137.7;
Pw_vec =  [903   875  819  782 732 710 680 630 581 540 ...
         507   472  437  405 383 315 286 248 195 146 ...
         107.5 80   53 ]*10^-3; % [W] 

P2wVsPw=[226 217 198 183 171 160 150 131 119 106 ...
         96  85  76  68  61  45  38  30  19  11  ...
         6   3   1 ]*10^-3; % [W]
end

if(strcmp(date,'06_11_16'))
% P2w Vs Temperature
Pw_1  = 200e-3;  %[W]

T   =   132:0.1:141;

P2wVsT =[0   1  1  2    1   2  2  2   3  4  ...
		 5   5  6  6    6   6  6  5.5 6  7  ...
		 9   11 13 13   14  14 13 13  13 13 ...
		 13  14 17 20   20  21 23 22  22 20 ...
		 18  17 16 17   18  20 22 22  23 22 ...
		 19  19 15 14   11  10 10 10  10 11 ...
		 11  11 11 10.5 9   8  6  4.5 4  3  ...
		 1   1  1  1    0.5 1   1 1   1  1  ... 
		 1.5 2  2  3    3   3   3 3   3  3  ...
		 3]*10^-3; %[W]

% P2w Vs Pw
Temp =  136.8;
Pw_vec =  [200 191  178 167 162 153 141 131 120 107 ...
		 97  89.5 79  71  60  47  0]*10^-3; % [W] 

P2wVsPw=[21  20   19  17  15  14  12  11  10  8 ...
		 6   5    4   3   2   1   0]*10^-3; % [W]
end

if(strcmp(date,'29_11_16'))
% P2w Vs Temperature
Pw_1  = 183e-3;  %[W]

T   =   129:0.1:138.4;

P2wVsT =[1.4    1.63    2.07    2.3     2.53    3.18    3.78    4.4     4.9     5.54    ...
         4.4    4.83    5.9     6.7     7.48    8       9.7     9.3     10.05   10.79   ...
         9.3    10.33   11      11.6    13      14.5    15.07   16      16.26   15.44   ...
         15.5   15      14.6    14.3    15      15.6    16.35   18.5    17.2    16.5    ...
         15.7   13.9    12.9    12.7    13.15   13.6    14      14.5    14.2    13.33   ...
         12.85  10.3    8       7.44    6.6     4.7     5.77    5.8     6       6.3     ...
         6.41   6.63    6.55    5.81    4.85    3.75    2.75    2.05    1.37    0.96    ...
         0.918  0.905   1.055   1.33    1.96    2.34    3.04    3.58    3.95    3.97    ...
         4.16   3.86    3.56    3.14    2.85    2.5     1.95    1.88    1.99    2.16    ...
         2.6    2.86    3.08    3.58    3.55]*10^-3; %[W]

% P2w Vs Pw
Temp    =  132.2;
Pw_vec  =  [200     188     181     171     161     150     139     128     119     109     ...
            99      89      80.5    70      61      50      39.8    30      21.72   10.54   ...
            0]*10^-3; % [W] 

P2wVsPw =  [21.9    17.8    15.8    14.8    13.3    11.4    10      9       7.9     6.7     ...
            5.6     4.5     3.7     3.1     2.325   1.6     1.08    0.626   0.36    0.106   ...
            0]*10^-3; % [W]
end


if(Print.Exp.P2wVsT)
    Pw_1        = Pw_1/(f*Pulse_wdt);
    P2wVsT      = P2wVsT/(f*Pulse_wdt);
    Efficiency  = P2wVsT/Pw_1;
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    set(FigHandle, 'Position', [70, 50, 980, 600]);
%     plot(T,P2wVsT/max(P2wVsT),'*-','Linewidth',2);
    plot(T,Efficiency,'*-','Linewidth',2);
    str = sprintf('Conversion Efficiency Vs T, Pw = %s [W]',num2str(Pw_1));
    title(str);
    ylabel('Conversion Efficiency');
    xlabel('T [^\circC]');

    [M,I] = max(Efficiency);
    text(T(I), Efficiency(I), ['\leftarrow = ', num2str(M),' ,T = ',num2str(T(I))]);
end


if(Print.Exp.P2wVsPw)
    Pw_vec = Pw_vec/(f*Pulse_wdt);
    P2wVsPw= P2wVsPw/(f*Pulse_wdt);
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    set(FigHandle, 'Position', [70, 50, 980, 600]);
    plot(Pw_vec,P2wVsPw,'*-','Linewidth',2);
    str = sprintf('P2w Vs Pw, T= %s [C]', num2str(Temp));
    title(str);
    xlabel('P\omega [W]');
    ylabel('P2\omega [W]');
    hold on;
    coefs = polyfit(Pw_vec, P2wVsPw, 2);
    P2wVsPwPoly = polyval(coefs,0:1:max(Pw_vec));
    plot(0:1:max(Pw_vec),P2wVsPwPoly,'r');
    xlim([min(Pw_vec) max(Pw_vec)])
    a1 = coefs(2);
    a2 = coefs(1);
    polyfit_str = ['P2\omega = ' num2str(a1) ' P\omega + ' num2str(a2) ' P\omega^2'];
    %text(7,7,polyfit_str)
    hold off;
    legend('Experiment',polyfit_str);
    
end

if(Print.Exp.effcncyVsPw)
    Pw_vec = Pw_vec/(f*Pulse_wdt);
    P2wVsPw= P2wVsPw/(f*Pulse_wdt);
    % conversion percentage efficiency vs Pw
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    set(FigHandle, 'Position', [70, 50, 980, 600]);
    Prcntg = 100*P2wVsPw./Pw_vec;
    Prcntg(find(isnan(Prcntg) == 1))=0;
    plot(Pw_vec,Prcntg,'*-','Linewidth',2);
    str = sprintf('Conversion Efficiency Vs Pw, T= %s [C]', num2str(Temp));
    title(str);
    xlabel('P\omega [W]');
    ylabel('% P2\omega/P\omega Steady state');
    hold on;
    coefs2 = polyfit(Pw_vec, Prcntg, 1);
    PrcntgPoly = polyval(coefs2,0:1:max(Pw_vec));
    plot(0:1:max(Pw_vec),PrcntgPoly,'r');
    xlim([min(Pw_vec) max(Pw_vec)])
    b1 = coefs2(1);
    polyfit_str2 = ['P2\omega/P\omega = ' num2str(b1) ' P\omega'];
    %text(7,7,polyfit_str2)
    hold off;
    legend('Experiment',polyfit_str2);
    [M,I] = max(Prcntg);
    text(Pw_vec(I), Prcntg(I), [num2str(M),' %']);
end

if (Print.Exp.P2w_vs_GradDiff)
    Pw_1        = Pw_1/(f*Pulse_wdt);
    P2wVsT      = P2wVsT/(f*Pulse_wdt);
    Efficiency  = P2wVsT/Pw_1;
    
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    set(FigHandle, 'Position', [70, 50, 980, 600]);
    DeltaT = Te - Ts;
    
    
    plot(DeltaT, Efficiency,'b*-','Linewidth',2);
    xlim([min(DeltaT) max(DeltaT)]);
    xlabel('\DeltaT [^\circC]');
    ylabel('Conversion Efficiency');

    str = sprintf('Conversion Efficiency Vs T, Pw = %s [W]',num2str(Pw_1));
    title(str);
    ylabel('Conversion Efficiency');
    xlabel('T [^\circC]');

    [M,I] = max(Efficiency);
    text(DeltaT(I), Efficiency(I), ['\leftarrow = ', num2str(M),' ,T = ',num2str(DeltaT(I))]);
     
end
