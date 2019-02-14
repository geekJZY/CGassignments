#ifndef OBJECTMOVING
#define OBJECTMOVING

#include <cmath>
#include <opencv2/opencv.hpp>
#include "colorDet.h"
#include "Quartics.h"

#define PI 3.1415927

using namespace std;
using namespace cv;


class ObjectMoving{
public:
    vector<Quadrics> spheres;
    ObjectMoving(vector<Quadrics> spheres);
    void rotateTowards(Vec3d dir, int idxS);
};


Vec3d rotateFun(Vec3d n, Vec3d v, double theta);

#endif
