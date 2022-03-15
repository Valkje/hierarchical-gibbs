function savePDF(filepath)
%SAVEPDF saves a figure as a PDF file, properly cropped.
%   MATLAB's savefig always saves a PDF figure as if it was placed on an A4
%   paper, which is not a very useful default. This function crops the PDF
%   to the size of the figure.

h = gcf;

set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
    'PaperSize', [pos(3), pos(4)])

print(h, filepath, '-dpdf', '-r0');

end

