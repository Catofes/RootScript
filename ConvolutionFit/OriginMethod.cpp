#include <RooRealVar.h>
#include <TApplication.h>
#include <TCanvas.h>
#include "OriginMethod.h"
#include "RooAddPdf.h"
#include "RooVoigtian.h"
#include "RooMyPdf.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooMinimizer.h"


using namespace RooFit;

int main(int argc, char **argv)
{
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    RooFit::PrintLevel(5);
    TApplication *myapp = new TApplication("App", &argc, argv);
    RooRealVar x("x", "x", 15, 0, 30);
    RooRealVar width("width", "width", 10, 8, 12);
    RooRealVar mean("mean", "mean", 10, 8, 12);
    RooRealVar sigma("sigma", "sigma", 0.5, 0.5, 0.5);
    RooMyPdf my_pdf_1("my_pdf_1", "MyPDF", x, mean, width, 0);
    RooMyPdf my_pdf_2("my_pdf_2", "MyPDF", x, mean, width, 1);
    RooMyPdf my_pdf_3("my_pdf_3", "MyPDF", x, mean, width, 2);
    RooMyPdf my_pdf_4("my_pdf_4", "MyPDF", x, mean, width, 3);
    my_pdf_4.cuda_gaus_prepare();
    RooVoigtian voigtian_pdf("voigtian_pdf", "VoigtianPdf", x, mean, width, sigma);
    RooPlot *xfram = x.frame();
    voigtian_pdf.plotOn(xfram, LineColor(1));
//    my_pdf_1.plotOn(xfram, LineColor(2));
    my_pdf_2.plotOn(xfram, LineColor(3));
    my_pdf_3.plotOn(xfram, LineColor(4));
//    my_pdf_4.plotOn(xfram, LineColor(5));
//    auto h_voigtian_pdf = voigtian_pdf.createHistogram("x", 1000000);
//    auto h_my_pdf_1 = my_pdf_1.createHistogram("x", 100000);
//    auto h_my_pdf_2 = my_pdf_2.createHistogram("x", 100000);
//    auto h_my_pdf_3 = my_pdf_3.createHistogram("x", 100000);
//
//    std::cout << "Diff h1" << h_voigtian_pdf->Chi2Test(h_my_pdf_1) << std::endl;
//    std::cout << "Diff h2" << h_voigtian_pdf->Chi2Test(h_my_pdf_2) << std::endl;
//    std::cout << "Diff h3" << h_voigtian_pdf->Chi2Test(h_my_pdf_3) << std::endl;

//    auto data = voigtian_pdf.generate(x, 10000);
//    auto h = data->createHistogram("x", 100, 0, 30);
//    h->Draw();
//    my_pdf_3.fitTo(*data);
//    my_pdf_3.plotOn(xfram);

    xfram->Draw();

//    TCanvas *c2 = new TCanvas("gaus_evaluate", "gaus_evaluate");
//    c2->cd();
//    my_pdf_1.h->Draw();
//    c2->SetLogy(true);
//    c2->Update();
//
//    TCanvas *c3 = new TCanvas("normal_evaluate", "normal_evaluate");
//    c3->cd();
//    my_pdf_2.h->Draw();
//    c3->SetLogy(true);
//    c3->Update();
//
//    TCanvas *c4 = new TCanvas("cuda_normal_evaluate", "cuda_normal_evaluate");
//    c4->cd();
//    my_pdf_3.h->Draw();
//    c4->SetLogy(true);
//    c4->Update();
//
//    TCanvas *c5 = new TCanvas("cuda_gaus_evaluate", "cuda_gaus_evaluate");
//    c5->cd();
//    my_pdf_4.h->Draw();
//    c5->SetLogy(true);
//    c5->Update();

    myapp->Run();
}