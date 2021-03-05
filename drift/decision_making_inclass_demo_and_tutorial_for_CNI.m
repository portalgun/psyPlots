% Decision-making: Drift Diffusion Model (DDM)
% Created by Yunshu Fan 03-22-2017

% This tutorial contains three parts:
% Part 1: in-class demo, includes solving a read-life decision problem and
% a perceptual decision task used in decision-making researches both using the
% drif-diffusion model

% Part 2: After class tutorial: For the same perceptual decision task,
% instead of simulate one trials, now we simulate multiple trials, and look
% at two important behavioral measurements: decision accuracy and decision time distribution

% Part3: Assignment

% My suggestions
% Whenever you see a new DDM with some modification to the basic parameters, always simulate the model yourself, and plot
% psychometric function and RT distribution yourself to see if that new
% model produce noticable difference from the basic model.

%% In-class demo:
%% Solving a real-life decision problem using DDM
% Choose a decision problem to solve: apple or banana
% The quesion is: choosing between A(apple   ) and B( banana  )

% Define momentary evidence: How strong does each piece of information
% supports A over B, on a scale of -1 to 1, -1 being strongly support B,
% and 1 being strongly support A.

%% Set decision bounds: 
% The total amount of evidence needed for choosing A:
bound_A = 4;
% The total amount of evidence needed for choosing B:
bound_B = - bound_A;
%% plot them out:
figure(1); clf;
subplot(2,1,1) % momentary evidence
hold on
plot([0,20],[1,1]*0,'-k','linewidth',1.5) % reference when evidence supports A and B equally
ylim([-1,1])
ylabel(sprintf('momentary evidence\nsupport A over B'))
xlabel('time')
set(gca,'xgrid','on')
set(gca,'ygrid','on')

subplot(2,1,2) % decision variable
hold on
plot([0,20],[1,1]*bound_A,'-k','linewidth',1.5) % bound A
plot([0,20],[1,1]*bound_B,'-k','linewidth',1.5) % bound B
ylim([bound_B-1.5,bound_A+1.5])
ylabel(sprintf('accumulated evidence\nsupport A over B'))
xlabel('time')
set(gca,'xgrid','on')
set(gca,'ygrid','on')

%% Start gethering information-converting into momentary evidence and
% accumulating into decision variable.
stop = false;
t = 0;
momentary_evidence = 0; % before you start gethering information, there is no evidence
decision_variable = 0;

while stop == false
    % time increases
    t = t+1;
    % convert information into momentary evidence
    momentary_evidence(t+1) = input('momentary evidence = ');
    % accumulate momentary evidence into decision variable
    decision_variable(t+1) = decision_variable(t) + momentary_evidence(t+1);
    
    % =======plotting==========
    subplot(2,1,1)
    plot([t-1,t],momentary_evidence(t:t+1),'-vr','linewidth',1.5)
    
    subplot(2,1,2)
    plot([t-1,t],decision_variable(t:t+1),'-vr','linewidth',1.5)
%     plot([t-1,t],[1,1]*bound_A,'-k','linewidth',1.5) % bound A
%     plot([t-1,t],[1,1]*bound_B,'-k','linewidth',1.5) % bound B
    % ======================
    
    % If decision variable hits bound A or bound B, decision stops,
    % otherwise, continue
    if (decision_variable(t+1) >= bound_A) || (decision_variable(t+1) <= bound_B)
        stop = true;
        decision_time = t;
        if decision_variable(t+1) >= bound_A
            choice = 1;
        else
            choice = 2;
        end
    else
        stop = false;
    end
    
end


%% questions to think about
% If we set higher bounds, i.e. needs more total evidence to commit to a
% decision, how will that change decision time? When do we want to do that?

% How about setting lower bounds? (hint: I want to make a decision fast.)

%% model a decision process in the random-dot motion discrimination task

% The stimulus in the task is a field of randomly positioned dots. A
% proportion of the dots moving coherently to left or right. We call the
% proportion "motion coherence".

% The model assumes that the motion information gets into the brain is
% converted into momentary evidence, which measures how strong the information supports
% leftward motion than rightward motion.

% You can think of the "momentary evidence" as taking a glimpse at the
% moving dots, and get a "feeling" of how strong do you think the dots move
% to the left based on that glimpse, and use a number to rate how strong your feeling is. And
% then take another glimpse.
% A big positive number indicates that you feel what you saw strongly
% supports leftward motion; a big negative number indicates that you feel
% what you saw strongly supports rightward motion; zero indicates that you
% feel ambivalent.

% Usually because of the noise in the stimulus as well as in our brain, the
% momentary evidence at one time would be different from that at another time,
% so we model the momentary evidence as a random variable drawing from a
% distribution. 
% The distribution is constructed in the following way:
% Based on our experience, when the coherence is high, we tend to feel that the evidence 
% is strong; when the coherence is low, the evidence is weak. 
% To model this, we set the mean of the momentary evidence to be proportional to motion coherence.
% For simplicity, we set the variance of the momentary evidence to 1.


% Now, let's simulate a single trial:

% set motion coherence in the trial. (+): dots move to left, (-): dots move to right
coh = 0.2; % 20% of the dots move coherently to left
m = coh *0.3 ; % mean of momentary evidence = coherence * scaling factor
v = 1; % variance of momentary evidence

% set bounds
bound_left = 25;
bound_right = -25;

decision_variable = 0; % starting from neutral
t = 0;
stop = false;
while stop == false
    t = t+1;
    % simulate momentary evidence
    momentary_evidence(t+1) = normrnd(m,v,1,1);

    % accumulate momentary evidence into decision variable
    decision_variable(t+1) = decision_variable(t) + momentary_evidence(t+1);
    
    % compare decision variable with bounds to determine if the decision is
    % stop at this step or not
    
    if (decision_variable(t+1) >= bound_left) || (decision_variable(t+1) <= bound_right)
        stop = true;
        decision_time = t;
        if decision_variable(t+1) >= bound_left
            choice = 1; % 1:choose left, 2:choose right.
        else
            choice = 2;
        end
    else
        stop = false;
    end
    
end

if choice == 1
    fprintf('choose A, decision time = %.d\n',decision_time)
else
    fprintf('choose B, decision time = %.d\n',decision_time)
end

% plot momentary evidence and decision variable
figure(3)
subplot(2,1,1)
hold all
plot([0,50],[1,1]*0,'k','linewidth',1.5) % only plot the first 50 steps
plot(0:50, momentary_evidence(1:51))
title('momentary evidence')
ylabel('support left over right')
ylim([-1,1]*3.5)

subplot(2,1,2)
hold all
plot([0,1500],[1,1]*bound_left,'k','linewidth',1.5)
plot([0,1500],[1,1]*bound_right,'k','linewidth',1.5)
plot(0:decision_time, decision_variable(1:decision_time+1))
title('decision variable')
ylim([-1,1]*(bound_left+10))

text(700,bound_left+5,'choose left')
text(700,bound_right-5,'choose right')

%% tutorial
% Let's do 100 trials by repeat the procedure above 100 times, then look at choice accuracy and decision time distribution
clear;
coh = 0.2; % 20% of the dots move coherently to left
m = coh * 0.1; % mean of momentary evidence
v = 1; % variance of momentary evidence

% set bounds
bound_left = 25;
bound_right = -25;

ntrials = 100; % do 100 trials
tmax = 10000; % set the maximum time you are willing to take before you quit the decision
momentary_evidence_ALL = nan(ntrials, tmax);
decision_variable_ALL = nan(ntrials, tmax);
choice_ALL = nan(ntrials,1);
decision_time_ALL = nan(ntrials,1);

for ii = 1:ntrials
    t = 0;
    stop = false;
    momentary_evidence = 0;
    decision_variable = 0;
    while stop == false
        t = t+1;
        % if t > tmax, start from beginning
        if t > tmax
            t = 1;
            momentary_evidence = 0;
            decision_variable = 0;
        end
        
        % simulate momentary evidence
        momentary_evidence(t+1) = normrnd(m,v,1,1);

        % accumulate momentary evidence into decision variable
        decision_variable(t+1) = decision_variable(t) + momentary_evidence(t+1);

        % compare decision variable with bounds to determine if the decision is
        % stop at this step or not

        if (decision_variable(t+1) >= bound_left) || (decision_variable(t+1) <= bound_right)
            stop = true;
            decision_time = t;
            if decision_variable(t+1) >= bound_left
                choice = 1; % 1:choose left, 2:choose right.
            else
                choice = 2;
            end
        else
            stop = false;
        end

    end
    
    momentary_evidence_ALL(ii,1:t+1) = momentary_evidence;
    decision_variable_ALL(ii,1:t+1) = decision_variable;
    choice_ALL(ii,1) = choice;
    decision_time_ALL(ii,1) = decision_time;
end

% compute proprotion of leftward and rightward choices
n_left = sum(choice_ALL == 1);
p_left = n_left/ntrials;
p_right = 1-p_left;
% compute accuracy based on dots motion direction
% When coh > 0, dots moving to left, leftward choices are correct; 
% when coh < 0, dots moving to right, rightward choices are correct;
if coh>0
    p_correct = p_left;
else
    p_correct = p_right;
end

% look at decision time distribution
mean_decision_time = mean(decision_time_ALL);
figure(4)
clf
subplot(2,2,1)
bar(3, p_correct)
ylim([0,1])
xlim([0,6])
ylabel('Prob. correct');
title('accuracy')
set(gca,'xticklabel',{})

subplot(2,2,2)
hist(decision_time_ALL);
title('decision time distribution')
xlabel('time')
ylabel('# of trials')

% look at momentary evidence and decision variable averaged across trials
% as a function of elapsed time
momentary_evidence_avg = nanmean(momentary_evidence_ALL,1);
momentary_evidence_avg(isnan(momentary_evidence_avg)) = [];

decision_variable_avg = nanmean(decision_variable_ALL,1);
decision_variable_avg(isnan(decision_variable_avg)) = [];


subplot(2,2,3)
hold on;
t_plot = 50; % plot the first 50ms
plot(0:t_plot,momentary_evidence_avg(1:t_plot+1),'r','linewidth',1.5)
plot([0,t_plot+1],[0,0],'k','linewidth',1.5)
ylim([-1,1])
title('average momentary evidence')
xlabel('elapsed time in a trial')

subplot(2,2,4)
hold on;
if length(decision_variable_avg)>=1000
    t_plot = 1000; % plot the first 10s00ms
else
    t_plot = length(decision_variable_avg)-1;
end
plot(0:t_plot,decision_variable_avg(1:t_plot+1),'r','linewidth',1.5)
plot([0,t_plot+1],[0,0],'k','linewidth',1.5)
title('average decision variable')
xlabel('elapsed time in a trial')

% the average momentary evidence is relatively steady over time, while the
% average decision variable ramps up (or down, if the coherence is
% negative)

% If you try multiple levels of coherence, e.g. 0.1, 0.3, 0.5, and plot
% them together. You will see the patterns look like the neural activities
% in brain areas MT (momentary evidence) and LIP (decision variable)

%% Assignment
% Part 1: How does bound height influence accuracy and decision time?
% Use coherence = 0.2. Set bound heights to different values, such as 10, 20, 40. 
% Plot the accuracy as a function of bound height.
% Plot the decision time distribution and mean decision time for different bound heights.
% You should get something like: the higher the bounds, the more accurate and slow a decision is, on average.
% This phenomenon is called "speed-accuracy trade-off".
% Intuitively, the higher the bounds are, the more evidence you will need
% to collect, which increases accuracy, but collecting more evidence will increase decision
% time.

% Think about how this would inform you about how you would do for the
% daily decisions you need to make. For me, I'm usually a perfectionist. In
% order to acheive high accuracy, I tend to set my bound very high, which also means that it
% takes a long time for me to make a decision and take action.
% (That's also my excuse or reason for procrastination.)
% When I need to make a quick decision, or just want push myself out of the
% procrastination zone, I'll tell myself that I'm in a "speed mode", which
% means I'm happy to sacrifice a little accuracy, but in doing that, I can
% get my work done on time. And that also inform me to be more careful
% while I'm doing the work, so that I can use my attention to compansate for the accuracy loss due
% to prioritizing speed.

% Part 2: How does stimulus strength (i.e. motion coherence) influence
% accuracy and decision time?
% Set the bound heights to a constant value, such as 25. Set coherence =
% 0.1, 0.3, 0.5, 0.7
% Plot(1): Plot the accuracy as a function of coherence.
% Plot(2): Plot the mean decision time OF CORRECT TRIALS as a function of coherence.
% You should get something like this: accuracy increases as coherence
% increases; mean decision time of correct trials decreases as coherence
% increases.

% Plot(1) is called "psychometric function". Plot(2) is called "chronometric function"
% When we study decision-making in the lab, we ask subjects to do this decision task for 
% different coherence levels, under different conditions. Then we measure
% the corresponding psychometric and chronometric functions from their
% behaviors.
% By finding the values of bound height that give us the psychometric and chronometric
% functions that best match to the ones we measured from subjects, we can
% find out what bound heights they used and what factors influence how they
% set their bounds.
