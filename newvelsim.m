function [B1 Gz] = newvelsim (seg,ref,n)
%
%This function is to create a B1 and Gz pulse with velocity encoded
%gradient and a composite RF pulse with double refocussers interleaved
%between them - Refer Qin et al MRM 2016
%
% Inputs: seg - choose from sinc,sinc7(7lobe sinc),bhard,sech followed by FA of
%         10,20,30 or 90
%         example - bhard20 for hard 20 degree pulse 
%         ref - choose from sinc,sinc7,bhard,sech followed by FA of
%         180
%         example - bhard180 for hard reforcusser pulse
%         n - choose according to FA chosen for seg pulse 
%         8 for 9 segs(20 deg) , 5 for 6 segs (30deg), 17 for 18 segs (10 deg)
% Outputs: B1 and Gz of same length
% Author : Sai Abitha Srinivas 

dt = 1e-3;
B1 = [];
Gz = [];
Mzfinal = [];
mlev_phases = repmat([1 1 -1 -1   -1 1 1 -1   -1 -1 1 1    1 -1 -1 1], [1,4]);
hard90 = (0.0195e-4 )* ones(1/dt,1);
B1ref_dur = 0.8;
B1ref_len = B1ref_dur/dt;
B1ref_max = (1/B1ref_dur)* 0.117e-4; % Tesla ... amplitude of 1 ms hard 180
Gvs = 4*1e-4;  % T/cm
ramp_len = 0.2/dt;
trap = [
    [0:ramp_len-1]'/(ramp_len-1);
    ones(0,1);
    [ramp_len-1:-1:0]'/(ramp_len-1);
    ]*Gvs;
trap_len = length(trap);


%Sinc1lobe
Sinc1 = sinc(linspace(-4,4, 1/dt));
hw = hamming(1/dt)';
sinc1 = Sinc1 .* hw;
sinc1_area = sum((sinc1));

%Sinc1lobe10 degrees
%variable = sinc10

sinc10 = sinc1 * (0.0117/18) / sinc1_area ; 
sinc10_max = max((sinc10));
sinc10_area = sum((sinc10));
sinc10 = sinc10';

sinc20 = sinc1 * (0.0117/9) / sinc1_area ; 
sinc20_max = max((sinc20));
sinc20_area = sum((sinc20));
sinc20 = sinc20';

%Sinc1lobe 30 degrees
%variable = sinc30
sinc30 = sinc1 * (0.0117/6) / sinc1_area ; 
sinc30_max = max((sinc30));
sinc30_area = sum((sinc30));
sinc30 = sinc30';

%Sinc1lobe 90 degrees
%variable = sinc90
sinc90 = sinc1 * (0.0117/2) / sinc1_area ; 
sinc90_max = max((sinc90));
sinc90_area = sum((sinc90));
sinc90 = sinc90';

%Sinc1lobe 180 degrees
%variable = sinc90
sinc180 = sinc1 * (0.0117) / sinc1_area ; 
sinc180_max = max((sinc180));
sinc180_area = sum((sinc180));
sinc180 = sinc180';

%Sinc7lobe 
Sinc7 = sinc(linspace(-8,8, 4/dt));
hw2 = hamming(4/dt)';
sinc7 = Sinc7 .* hw2;
sinc7_area = sum((sinc7));

%Sinc7lobe 10 degrees
%variable = sinc710
sinc710 = sinc7 * (0.0117/18) / sinc7_area ; 
sinc710_max = max((sinc710));
sinc710_area = sum((sinc710));
sinc710 = sinc710';

%Sinc7lobe 20 degrees
%variable = sinc720
sinc720 = sinc7 * (0.0117/9) / sinc7_area ; 
sinc720_max = max((sinc720));
sinc720_area = sum((sinc720));
sinc720 = sinc720';

%Sinc7lobe 30 degrees
%variable = sinc730
sinc730 = sinc7 * (0.0117/6) / sinc7_area ; 
sinc730_max = max((sinc730));
sinc730_area = sum((sinc730));
sinc730 = sinc730';

%Sinc7lobe 90 degrees 
%variable = sinc790
sinc790 = sinc7 * (0.0117/2) / sinc7_area ; 
sinc790_max = max((sinc790));
sinc790_area = sum((sinc790));

%Sinc7lobe 180 degrees
%variable = sinc7180
sinc7180 = sinc7 * (0.0117) / sinc7_area ; 
sinc7180_max = max((sinc7180));
sinc7180_area = sum((sinc7180));
sinc7180 = sinc7180';

%hard 10 degrees
bhard10 = hard90( 1 : 333);
%hard 20 degrees
bhard20 = hard90( 1 : 667);
%hard 30 degrees
bhard30 = hard90( 1 : 1000);
%hard 90 degrees
bhard90 = (B1ref_max/2) * ones(B1ref_len,1);
%hard 180 degrees
bhard180 = B1ref_max * ones(B1ref_len,1);

%sech
mysech = genSech180(0.5, 1)';
%Sech 10 degrees
Sech10 = 1000*(0.117e-4/18) * mysech/sum(mysech);
%Sech 20 degrees
Sech20 = 1000*(0.117e-4/9) * mysech/sum(mysech);
%Sech 30 degrees
Sech30 = 1000*(0.117e-4/6) * mysech/sum(mysech);
%Sech 90 degrees
Sech90 = 1000*(0.117e-4/2) * mysech/sum(mysech);
%Sech 180 degrees
Sech180 = 1000*(0.117e-4) * mysech/sum(mysech);

x = seg;
y = ref;

for i = 1: n 
    
    Gz = [Gz;
        zeros(size(x));
        trap;
        zeros(size(y));                                                       
        -trap;
        trap;
        zeros(size(y));
        -trap;
        ];
    
    B1 = [B1;
        x;
        zeros(size(trap));
        y * mlev_phases(2*n-1);  
        zeros(size(trap));
        zeros(size(trap));
        y * mlev_phases(2*n);  
        zeros(size(trap));       
        ];

end

B1 = [B1; x];
Gz = [ Gz ;zeros(size(x))]; 


