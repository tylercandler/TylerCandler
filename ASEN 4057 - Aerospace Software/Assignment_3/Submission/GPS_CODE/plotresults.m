function [] = plotresults(e_i,e_q, p_i, p_q,l_i, l_q, carrierfq, codefq, hand_Tracked_Code_Freq, hand_Tracked_Int_Freq, hand_Correlation_results, hand_Prompt_I, hand_Prompt_Q)

%Simple script to plot the results of the bistatic processing of the
%direct channel;  create two figures and plot important aspects

% figure(901)
% subplot(2,1,1),
plot(hand_Correlation_results,p_i .^2 + p_q .^ 2, 'g.-')
hold(hand_Correlation_results,'on')
grid
% subplot(2,1,1),
plot(hand_Correlation_results,e_i .^2 + e_q .^ 2, 'bx-')
% subplot(2,1,1),
plot(hand_Correlation_results,l_i .^2 + l_q .^ 2, 'r+-')
xlabel(hand_Correlation_results,'milliseconds')
ylabel(hand_Correlation_results,'amplitude')
title(hand_Correlation_results,'Correlation Results')
legend(hand_Correlation_results,'prompt','early','late')
% subplot(2,2,3),
plot(hand_Prompt_I,p_i)

grid
xlabel(hand_Prompt_I,'milliseconds')
ylabel(hand_Prompt_I,'amplitude')
title(hand_Prompt_I,'Prompt I Channel')
% subplot(2,2,4),
plot(hand_Prompt_Q,p_q)
hold(hand_Prompt_Q,'on');
grid
xlabel(hand_Prompt_Q,'milliseconds')
ylabel(hand_Prompt_Q,'amplitude')
title(hand_Prompt_Q,'Prompt Q Channel')

% figure(902)
% subplot(2,1,1),
plot(hand_Tracked_Code_Freq,1.023e6 - codefq)
hold(hand_Tracked_Code_Freq,'on');
grid
xlabel(hand_Tracked_Code_Freq,'milliseconds')
ylabel(hand_Tracked_Code_Freq,'Hz')
title(hand_Tracked_Code_Freq,'Tracked Code Frequency (Deviation from 1.023MHz)')
% subplot(2,1,2),
plot(hand_Tracked_Int_Freq,carrierfq)
hold(hand_Tracked_Int_Freq,'on');
grid
xlabel(hand_Tracked_Int_Freq,'milliseconds')
ylabel(hand_Tracked_Int_Freq,'Hz')
title(hand_Tracked_Int_Freq,'Tracked Intermediate Frequency')
end
