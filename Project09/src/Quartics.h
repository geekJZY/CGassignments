#ifndef QUARTICS
#define QUARTICS

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "definition.h"
#define PI 3.1415927

using namespace std;
using namespace cv;

///Define Objects
class Quadrics{
private:
    Vec4f ColorH0, ColorH1, ColorH2;
    Mat textureMap;
public:
    enum SphereType{
        TEXTURE,
        SOLID_TEXTURE,
        SOLID_TEXTURE_PER,
        HOMO,
        INFINITE_SPHERE
    };
    SphereType sphereType;
    float nRR;
    int texture;
    int a02,a12,a22,a21,a00;
    double s0, s1, s2;
    float x_solid, y_solid;
    Vec3d n0, n1, n2;
    Vec4f  ColorH3;
    Vec3d pc;
    Vec3d planeA, planeY, planeX, planeC, projPoint;
    bool checkIn(Vec3d pos);
    double rayTracer(Vec3d pe, Vec3d npe);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, int texture_, float nRR_, Vec4f ColorH0_, Vec4f ColorH1_, Vec4f ColorH2_, Vec4f ColorH3_);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, int texture_, float nRR_, Mat textureMap_, Vec4f ColorH3_);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3d pc_, int texture_, float nRR_, Mat textureMap_, Vec3d planeA_, Vec3d planeX_, Vec3d planeC_, Vec3d projPoint_, float x_solid, float y_solid, Vec4f ColorH3_, SphereType _sphereType = SOLID_TEXTURE);

    vector<Vec4f> TexColorDeter(Vec3d Ph);
    Vec3d nH(Vec3d Ph);
};

#endif
