
function [timeout,obsAuto]=createAutoCor(a,bb,s,timeL)

% define measerment error
bb_un = bb-(bb*0.01);
bb_up = bb+(bb*0.01);


timecourAuto=sim_autocor(a,bb,1,timeL);

obsAuto=corr(timecourAuto(1:end-1),timecourAuto(2:end));


if   obsAuto>bb_un&obsAuto<bb_up


%nothing here

else


while  or(obsAuto<bb_un,obsAuto>bb_up)

timecourAuto=sim_autocor(0,bb,1,timeL);

obsAuto=corr(timecourAuto(1:end-1),timecourAuto(2:end));      

end


end

timeout=timecourAuto;
