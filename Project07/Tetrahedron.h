#ifndef TETRAHEDRON
#define TETRAHEDRON

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "Quartics.h"
#include "colorUn.h"

using namespace std;
using namespace cv;

///Define Tetrahedron
class Tetrahedron{
private:

public:
    Vec4f  ColorH3;
    vector<vector<Vec3d>> vecPoints;
    vector<Vec3d> As;
    vector<Vec3d> norms;
    Vec4f ColorH0, ColorH1, ColorH2;
    Quadrics quad;

    Tetrahedron(vector<Vec3d> points, Vec4f colorH0_, Vec4f colorH1_, Vec4f colorH2_, Vec4f colorH3_, Quadrics quad_);
    colorUnite rayTrace(Vec3d pe, Vec3d nep);
};

#endif // TETRAHEDRON
