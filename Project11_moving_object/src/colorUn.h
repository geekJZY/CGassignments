#ifndef COLOR_UN
#define COLOR_UN

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "definition.h"

using namespace std;
using namespace cv;

///Define colorUnite
class colorUnite{
private:

public:
    Vec4f ColorH0, ColorH1, ColorH2;
    Vec4f  ColorH3;
    Vec3d Ph;
    Vec3d norm;

    colorUnite(void);
    colorUnite(Vec3d Ph_, Vec3d norm_, Vec4f colorH0_, Vec4f colorH1_, Vec4f colorH2_, Vec4f colorH3_);
};

#endif
