#ifndef PLANES
#define PLANES

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "definition.h"
#define PI 3.1415927

using namespace std;
using namespace cv;

///Define Objects

class Planes{
private:
    Vec4f ColorH0, ColorH1, ColorH2;
    enum PlaneType{
        TEXTURE,
        HOMO
    };
    PlaneType planeType;
    Mat textureMap;
    Mat normalMap;
    float s0, s1; //params for texture map
public:
    int texture;
    float nRR;
    Vec3d pc;
    Vec3d n1, n2, n0;
    Vec4f ColorH3;
    Planes(Vec3d n2_, Vec3d pc_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_);
    Planes(Vec3d n2_, Vec3d pc_, Vec3d n1_, int texture_, float nRR_, Mat textureMap_, Mat normalMap_, float s0_, float s1_, Vec4f ColorH3_);
    vector<Vec4f> TexColorDeter(Vec3d ph);
    bool checkIn(Vec3d pos);
    Vec3d nH(void);
    Vec3d nH_Tex(Vec3d Ph);
    double rayTracer(Vec3d pe, Vec3d npe);
};

#endif
