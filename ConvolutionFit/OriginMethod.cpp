#include <RooRealVar.h>
#include <TApplication.h>
#include <TCanvas.h>
#include "OriginMethod.h"
#include "RooAddPdf.h"
#include "RooVoigtian.h"
#include "RooMyPdf.h"
#include "RooPlot.h"

using namespace RooFit;

int main(int argc, char **argv)
{
    TApplication *myapp = new TApplication("App", &argc, argv);
    RooRealVar x("x", "x", 15, 0, 30);
    RooRealVar width("width", "width", 10, 10, 10);
    RooRealVar mean("mean", "mean", 10, 10, 10);
    RooRealVar sigma("sigma", "sigma", 0.5, 0.5, 0.5);
    RooMyPdf my_pdf_1("my_pdf_1", "MyPDF", x, mean, width, 0);
    RooMyPdf my_pdf_2("my_pdf_2", "MyPDF", x, mean, width, 1);
    RooVoigtian voigtian_pdf("voigtian_pdf", "VoigtianPdf", x, mean, width, sigma);
    RooPlot *xfram = x.frame();
    voigtian_pdf.plotOn(xfram, LineColor(1));
    my_pdf_1.plotOn(xfram, LineColor(2));
    my_pdf_2.plotOn(xfram, LineColor(3));
    xfram->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->cd();
    my_pdf_1.h->Draw();
    c2->SetLogy(true);
    c2->Update();
    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->cd();
    my_pdf_2.h->Draw();
    c3->SetLogy(true);
    c3->Update();
    myapp->Run();
}