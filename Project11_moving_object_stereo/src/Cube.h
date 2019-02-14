#ifndef CUBE
#define CUBE

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "definition.h"
#include "Quartics.h"
#include "colorUn.h"

using namespace std;
using namespace cv;

///Define Cube
class Cube{
private:
    int traceTimeComp(vector<double> traceVecTime);
    Mat textureMap;
    Mat normalMap;
public:
    Vec4f  ColorH3;
    vector<vector<Vec3d>> vecPoints;
    vector<Vec3d> As;
    vector<Vec3d> norms;
    Vec4f ColorH0, ColorH1, ColorH2;
    Quadrics quad;
    int texture;
    float nRR;

    Vec3d pointSave;
    vector<Vec4f> colorSave;
    Vec3d normalSave;
    float nRR_Save;

    Cube(vector<Vec3d> points, Vec4f colorH0_, Vec4f colorH1_, Vec4f colorH2_, Vec4f colorH3_,int texture_, float nRR_, Quadrics quad_, Mat textureMap_, Mat normalMap_);
    double rayTrace(Vec3d pe, Vec3d nep);
};

#endif // CUBE
