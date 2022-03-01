function save_figures( FolderName, number_of_figures_excluded )
% save figures into given folder


% % Ripped from Jan's answer and others on MATLAB help online. SAM 3/5/21
% FolderName = tempdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)-number_of_figures_excluded

    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    set(0, 'CurrentFigure', FigHandle);
    savefig(fullfile(FolderName, [FigName '.fig']));

end


end