#include <RooRealVar.h>
#include "OriginMethod.h"
#include "RooAddPdf.h"
#include "RooVoigtian.h"
#include "RooMyPdf.h"

using namespace RooFit;

int main()
{
    RooRealVar x("x", "x", 15, 10, 20);
    RooRealVar width("width", "width", 10, 8, 12);
    RooRealVar mean("width", "width", 10, 8, 12);
    RooRealVar sigma("width", "width", 0.5, 0.5, 0.5);
    RooMyPdf my_pdf("my_pdf", "MyPDF", x, mean, width);
    RooVoigtian voigtian_pdf("voigtian_pdf", "VoigtianPdf", x, mean, width, sigma);
    RooPlot *xfram = x.frame();
    my_pdf.plotOn(xfram);
    voigtian_pdf.plotOn(xfram);
    
}