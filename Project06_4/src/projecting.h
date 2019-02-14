#ifndef PROJECTING
#define PROJECTING

#include <cmath>
#include <opencv2/opencv.hpp>
#include "colorDet.h"

#define PI 3.1415927

using namespace std;
using namespace cv;


class projecting{
public:
    Vec3d planeA;
    Vec3d planeC;
    Vec3d projP;
    Vec3d n0;
    Vec3d n1;
    float Sx;
    float Sy;
    Mat projImg;
    projecting(Vec3d planeA, Vec3d planeC, Vec3d n0, Vec3d n1, Vec3d projP, float Sx, float Sy, Mat projImg);
    Vec4f colorCal(Vec3d pos, colorDet &colorDeter);
};

#endif
