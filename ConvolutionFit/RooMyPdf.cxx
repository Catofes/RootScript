/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

// Your description goes here... 

#include <complex>
#include <RooMath.h>
#include "Riostream.h"

#include "RooMyPdf.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TMath.h"
#include "iostream"

ClassImp(RooMyPdf)

RooMyPdf::RooMyPdf(const char *name, const char *title,
                   RooAbsReal &_x, RooAbsReal &_mean, RooAbsReal &_width, int _method) :
        RooAbsPdf(name, title),
        x("x", "x", this, _x),
        mean("mean", "mean", this, _mean),
        width("width", "Breit-Wigner Width", this, _width),
        method(_method)
{
}


RooMyPdf::RooMyPdf(const RooMyPdf &other, const char *name) :
        RooAbsPdf(other, name),
        x("x", this, other.x),
        mean("mean", this, other.mean),
        width("width", this, other.width),
        method(other.method)
{
}


Double_t RooMyPdf::sub_f(Double_t t) const
{
    Double_t w = (width > 0) ? width : -width;
    Double_t arg = t - mean;
    return (1. / (arg * arg + 0.25 * w * w));
}

Double_t RooMyPdf::sub_sigma(Double_t t) const
{
    return 0.5;
}

// F(t) = f(t) * \frac{1}{\sigma (t)} \exp{-\frac{(t-x)^2}{2*\sigma^2(t)}}
Double_t RooMyPdf::sub_evaluate(Double_t t) const
{
    Double_t s = (sub_sigma(t) > 0) ? sub_sigma(t) : -sub_sigma(t);
    Double_t arg = t - x;
    Double_t coef = -0.5 / (s * s);
    return sub_f(t) * exp(coef * arg * arg) * 1 / s;
}

Double_t RooMyPdf::gaus_evaluate() const
{
    Double_t upper = x.max();
    Double_t lower = x.min();
    Double_t result = 0;

    double x0[5] = {0,
                    sqrt(245 - 14 * sqrt(70)) / 21.,
                    -sqrt(245 - 14 * sqrt(70)) / 21.,
                    sqrt(245 + 14 * sqrt(70)) / 21.,
                    -sqrt(245 + 14 * sqrt(70)) / 21.,
    };

    double a0[5] = {
            128 / 225.,
            (322 + 13 * sqrt(70)) / 900.,
            (322 + 13 * sqrt(70)) / 900.,
            (322 - 13 * sqrt(70)) / 900.,
            (322 - 13 * sqrt(70)) / 900.,
    };

    for (int i = 0; i < 5; i++) {
        result += a0[i] * sub_evaluate((upper - lower) / 2. * x0[i] + (upper + lower) / 2.);
    }
    return (upper - lower) / 2. * result;
}

Double_t RooMyPdf::normal_evaluate() const
{
    int cut = 10000;
    double upper = x.max();
    double lower = x.min();
    double step = (upper - lower) / cut;
    double result = 0;
    for (int i = 0; i < cut; i++) {
        double point = lower + (i + 0.5) * step;
        result += sub_evaluate(point) * step;
    }
    return result;
}

Double_t RooMyPdf::evaluate() const
{
    switch (method) {
        case 0:
            return gaus_evaluate();
        case 1:
            return normal_evaluate();
    }
}

