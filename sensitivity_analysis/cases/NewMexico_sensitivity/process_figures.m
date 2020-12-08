for i=1:5
    
    fig_id = 100+(i-1);
    h = figure(fig_id)
    grid on
    xlim([1.95 2.15])
    xlabel('t [s]')
    title('10 LHS realizations of PH2 yaw model parameter, run 941')
    ylabel(strcat('F_N at section ',num2str(i)));

    savefig(strcat('Results/NewMexico/example_FN_PH2_run941_section',num2str(i)));
    saveas(h,strcat('Results/NewMexico/example_FN_PH2_run941_section',num2str(i)),'png');
end