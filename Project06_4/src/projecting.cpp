#include "projecting.h"

projecting::projecting(Vec3d planeA, Vec3d planeC, Vec3d n0, Vec3d n1, Vec3d projP, float Sx, float Sy, Mat projImg){
    this -> planeA = planeA;
    this -> projP = projP;
    this -> n0 = n0;
    this -> n1 = n1;
    this -> Sx = Sx;
    this -> Sy = Sy;
    this -> projImg = projImg;
    this -> planeC = planeC;
    if((fabs(n0.dot(n1)) > 1e-4)||(fabs(planeA.dot(n1)) > 1e-4)||(fabs(n0.dot(planeA)) > 1e-4)) cout << "orientation is not orthogonal" << endl;
}

Vec4f projecting::colorCal(Vec3d pos, colorDet &colorDeter){
    Vec3d projLight, interPoint;
    projLight = (pos - projP)/norm((pos - projP));
    if(fabs(planeA.dot(projLight)) < 1e-4) return Vec4f(0,0,0,0);
    double t = planeA.dot(planeC - projP)/planeA.dot(projLight);
    Vec3d interP = projP + t * projLight;
    //calculate the intersecting point
    if(((interP-planeC).dot(n0) < Sx)&&((interP-planeC).dot(n0) >= 0)&&((interP-planeC).dot(n1) < Sy)&&((interP-planeC).dot(n1) >= 0)){
        int x = max(min(int(((interP-planeC).dot(n0)/Sx)*projImg.cols),projImg.cols-1),0);
        int y = max(min(int(((interP-planeC).dot(n1)/Sy)*projImg.rows),projImg.rows-1),0);
        Vec3b colorRGB = projImg.at<Vec3b>(y,x);
        //cout << "in range" << endl;
        Vec4f colorReturn = Vec4f(colorRGB[0], colorRGB[1], colorRGB[2], 255);

        pair<int, Vec3d> rayTracing = colorDeter.spaceTracer(projLight, projP);
        if(norm(rayTracing.second - pos) < 1e-4){
            //cout << colorReturn << endl;
            return colorReturn;
        }
        else return Vec4f(0,0,0,0);
    }
    else return Vec4f(0,0,0,0);
}
