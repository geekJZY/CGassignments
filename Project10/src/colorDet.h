#ifndef COLORDET
#define COLORDET

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "Quartics.h"
#include "plane.h"
#include "lightSource.h"
#include "definition.h"
#include "Cube.h"

using namespace std;
using namespace cv;

class colorDet{
private:
    int traceTimeComp(vector<double> traceVecTime);
    Vec3d jitter(Vec3d n);

    Mat projectImg;
    Vec3d planeA;
    Vec3d planeX;
    Vec3d planeC;
    float x_solid;
    float y_solid;
    Vec3d projP;
    Vec3d ph;

public:
    vector<Planes> planes;
    vector<Quadrics> spheres;
    vector<Cube> cubes;
    vector<lightSource> lighters;
    Vec3d Peye;
    Vec4f ColorH3;
    double kb, k0;

    colorDet(vector<Planes> _planes, vector<Quadrics> _spheres, vector<Cube> _cubes,vector<lightSource> _lighters, Vec3d _Peye, Vec4f ColorH3_, double kb_, double k0_);
    pair<int, Vec3d> spaceTracer(Vec3d nep, Vec3d pe);
    Vec4f colorReturn(int objFrom, int objInx, Vec3d ph, Vec3d nep, int tierCnt);
    void projectColorInit(Mat projectImg_, Vec3d planeA_, Vec3d planeX_, Vec3d planeC_, float x_solid_, float y_solid_, Vec3d projP_);
    Vec4f projectColorDet(Vec3d ph);
};

#endif // COLORDET
