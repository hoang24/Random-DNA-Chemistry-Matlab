function plotChem(timeVec, conS)
    figure;
    hold on;
    grid;
    title('Random DNA Strand Circuit');
    xlabel('time (s)');
    ylabel('concentration (M)');
    plot(timeVec, cell2mat(conS(2:end, :)));
    for i = 1 : size(conS, 2)
        legendInfo{i} = num2str(conS{1, i});
        legend(legendInfo);
    end
    hold off;
end