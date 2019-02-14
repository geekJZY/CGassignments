#ifndef LIGHT_SOURCE
#define LIGHT_SOURCE

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "definition.h"
#define PI 3.1415927

using namespace std;
using namespace cv;


///Define Light Source
class lightSource{
public:
    enum Mode
    {
      SPOT,
      DIRECTIONAL,
      POINTLIGHT
    };

    Vec3d Pl, nLH_D; //light point OR light Direction
    double spotRange;
    double ks, kd;
    Vec4f ColorL;   //color of light
    //params for color calculation
    Mode lightMode;
    lightSource(Vec3d Pl, Vec3d nLH, double spotRange, double ks, double kd, Vec4f ColorL, Mode lightMode);
    //H0 ambient, H1 diffuse, H2 specular
    Vec4f colorSDcal(Vec3d Pe, Vec3d Ph, Vec3d nH, Vec4f ColorH0, Vec4f ColorH1, Vec4f ColorH2, Vec4f ColorH3);
};
#endif
