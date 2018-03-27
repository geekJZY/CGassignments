#ifndef LIGHT_SOURCE
#define LIGHT_SOURCE

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#define PI 3.1415927

using namespace std;
using namespace cv;


///Define Light Source
class lightSource{
public:
    enum Mode
    {
      SPOT,
      DIRECTIONAL
    };

    Vec3f Pl, nLH_D; //light point OR light Direction
    Vec4f ColorL;   //color of light
    //params for color calculation
    double ks, kd, k0, kb;
    Mode lightMode;
    lightSource(Vec3f Pl, Vec3f nLH_, Vec4f ColorL, double ks, double kd, double k0, double kb, Mode lightMode);
    //H0 ambient, H1 diffuse, H2 specular
    Vec3b colorDeter(Vec3f Pe, Vec3f Ph, Vec3f nH, Vec4f ColorH0, Vec4f ColorH1, Vec4f ColorH2, Vec4f ColorH3);
};

#endif
