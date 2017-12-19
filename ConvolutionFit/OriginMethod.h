#include "TF1.h"

class OriginConvolutionFit
{
public:
    OriginConvolutionFit(TF1 *sigma, TF1 *f);

    ~OriginConvolutionFit()
    {};

private:
    TF1 *_sigma;
    TF1 *_f;
};

