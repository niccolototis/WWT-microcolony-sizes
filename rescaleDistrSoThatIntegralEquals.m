function [pdf_values]=rescaleDistrSoThatIntegralEquals(pdf_values,x2_widths,mustintegrateTo)
if any(pdf_values) % If I have all zeros, do nothing
    M_NisTilde=false;
    integralOfTheMarg=computeMarg(pdf_values(:),M_NisTilde,[],x2_widths(:),[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
    proportionalityK=mustintegrateTo/integralOfTheMarg;
    pdf_values=pdf_values*proportionalityK;
end
end
