function [pdf_values]=rescaleDistrSoThatSumEquals(pdf_values,mustSumTo)
if any(pdf_values) % If I have all zeros, do nothing
    SumOfTheMarg=sum(pdf_values);
    proportionalityK=mustSumTo/SumOfTheMarg;
    pdf_values=pdf_values*proportionalityK;
end
end
