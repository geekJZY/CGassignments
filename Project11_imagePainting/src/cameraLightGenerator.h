#ifndef CAMERALIGHTGENE
#define CAMERALIGHTGENE

#include <cmath>
#include <opencv2/opencv.hpp>

#define PI 3.1415927

using namespace std;
using namespace cv;


class lightShooter{
public:
    Vec3d Peye, n2;
    float f;
    Vec3d n0, n1;
    vector<Vec3d> pointsAroundC;
    lightShooter(Vec3d pe,Vec3d v2, float f_);
    pair<vector<Vec3d>, vector<Vec3d>> lightGene(Vec3d pos);
};

#endif
