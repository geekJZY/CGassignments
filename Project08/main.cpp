
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "src/colorUn.h"
#include "src/lightSource.h"
#include "src/Quartics.h"
#include "src/plane.h"
#include "src/Tetrahedron.h"
#include "src/Cube.h"
#include "src/colorDet.h"

#define PI 3.1415927

using namespace cv;
using namespace std;


int main(int argc, char *argv[]){

    vector<Scalar> color{Scalar(255,224,147),Scalar(93,66,255),Scalar(255,255,255),Scalar(230,180,80)};
    Mat img(Size(640,480),CV_8UC3,color[0]);
    Mat texturePlane = imread("Texture.jpg");
    Mat wallPaper = imread("wall_paper.jpg");
    Mat sphereTexture = imread("ball_image_small.jpg");
    Mat sphereTexture2 = imread("ball_image_10_small.jpg");
    Mat universeTexture = imread("universe3.jpg");
    Mat normalTexture = imread("planenorm.png");
    Mat projImg = imread("projImg.jpg");
    double Sx = 16, Sy = 12;
    Vec3d pos;
    Vec3d Ph;
    int objectIdx;
    int xSample = 2, ySample = 2;
    double randX,randY;
    int M = img.cols, N = img.rows;

    //the efficient of background and outline
    double kb = 0.5, k0 = 0.8;

    //Graphics
    vector<Quadrics> spheres;
    vector<Planes> planes;
    vector<Cube> cubes;
    vector<lightSource> lighters;


    ///define Tetrahedron
//    Quadrics quad(1,1,1,0,-1,2*sqrt(3),2*sqrt(3),2*sqrt(3),Vec3d(0,12,2),TX_OPAQUE, 1,Vec4f(0,0,0,0), Vec4f(0,0,0,0), Vec4f(0,0,0,0),Vec4f(0,0,0,255));
//    vector<Vec3d> points;
//    points.push_back(Vec3d(2,10,4));
//    points.push_back(Vec3d(-2,10,4));
//    points.push_back(Vec3d(-2,14,4));
//    points.push_back(Vec3d(2,14,4));
//    points.push_back(Vec3d(2,14,0));
//    points.push_back(Vec3d(-2,14,0));
//    points.push_back(Vec3d(-2,10,0));
//    points.push_back(Vec3d(2,10,0));
//    cubes.push_back(Cube(points, Vec4f(53,36,155,255),Vec4f(113,86,255,255),Vec4f(1930,1660,2550,2550),Vec4f(0,0,0,255), TX_OPAQUE, 1.33, quad, sphereTexture2, sphereTexture));

    spheres.push_back(Quadrics(1,1,1,0,-1,5,5,5,Vec3d(0,0,5),TX_OPAQUE,sphereTexture,Vec3d(0,-1,0),Vec3d(-1,0,0),Vec3d(6,-6,11),Vec3d(0,-12,5),12,12,Vec4f(0,0,0,255),Quadrics::SOLID_TEXTURE));
    spheres.push_back(Quadrics(1,1,1,0,-1,3,3,3,Vec3d(0,12,3),TX_OPAQUE,sphereTexture2,Vec3d(0,-1,0),Vec3d(-1,0,0),Vec3d(4,6,7),Vec3d(0,4,4),8,8,Vec4f(0,0,0,255),Quadrics::SOLID_TEXTURE));

    //spheres.push_back(Quadrics(1,1,1,0,-1,100,100,100,Vec3d(0,0,0),TX_INFINITE_SPHERE,universeTexture,Vec3d(0,-1,0),Vec3d(-1,0,0),Vec3d(100,100,100),Vec3d(0,4,4),120,60,Vec4f(0,0,0,255),Quadrics::INFINITE_SPHERE));
    planes.push_back(Planes(Vec3d(0,0,1),Vec3d(-10,-10,0),Vec3d(0,1,0),TX_OPAQUE, 1,texturePlane,normalTexture,10,10,Vec4f(10,10,10,255)));
//    planes.push_back(Planes(Vec3d(1,0,0),Vec3d(-10,-10,0),Vec3d(0,1,0),TX_OPAQUE, 1,wallPaper,normalTexture,10,10,Vec4f(10,10,10,255)));
//    planes.push_back(Planes(Vec3d(0,1,0),Vec3d(-10,-10,0),Vec3d(1,0,0),TX_OPAQUE, 1,wallPaper,normalTexture,10,10,Vec4f(10,10,10,255)));
//    planes.push_back(Planes(Vec3d(-1,0,0),Vec3d(40,40,50),Vec3d(0,1,0),TX_OPAQUE, 1,wallPaper,normalTexture,10,10,Vec4f(10,10,10,255)));
//    planes.push_back(Planes(Vec3d(0,-1,0),Vec3d(40,40,50),Vec3d(1,0,0),TX_OPAQUE, 1,wallPaper,normalTexture,10,10,Vec4f(10,10,10,255)));
//    planes.push_back(Planes(Vec3d(0,0,-1),Vec3d(40,40,50),Vec3d(1,0,0),TX_OPAQUE, 1,wallPaper,normalTexture,10,10,Vec4f(10,10,10,255)));

    //area light
    lighters.push_back(lightSource(Vec3d(0,30,30), Vec3d(0.5773,-0.5773,-0.5773), 0.9,0.5,0.5,Vec4f(60,60,60,60), lightSource::POINTLIGHT));

    //camera
    Vec3d pe(40,25,26);
    Vec3d vUp(-1,0,1);
    Vec3d v2(20,10,10);
    Vec3d n0Vec(0,-1,0);
    Vec3d n2 = v2/norm(v2);

    Vec3d n0 = n0Vec/norm(n0Vec);
    Vec3d n1 = n0.cross(n2);

    Vec3d pcMat = pe + v2;
    Vec3d pc(pcMat[0], pcMat[1], pcMat[2]);
    Vec3d th_P;

    Vec3i colorSum;
    double cameraX,cameraY; //the coordinate in the camera
    Vec3d npe;

    //colorDeter
    colorDet colorDeter(planes, spheres, cubes, lighters, pe, Vec4f(10,10,10,100), kb, k0);

    pair<int, Vec3d> rayTracing;

    cout << n0 << endl << n1 << endl << n2 << endl;
    colorDeter.projectColorInit(projImg, Vec3d(0,-1.44,-1.92)/norm(Vec3d(0,-1.44,-1.92)), Vec(-1,0,0), Vec(13,1,1), 2, 2, Vec(15,0,0));
    //checkIn
    for(int j = 0; j < planes.size(); j ++){
        if(planes[j].checkIn(pe)){
            cout << "fail in checking in " << endl;
            return 0;
        }
    }

    for(int m = 0; m < M; m ++)
        for(int n = 0; n < N; n ++)
        {
            //with anti-aliasing
            randX = rand()%1000/1000.0;
            randY = rand()%1000/1000.0;
            colorSum = Vec3i(0,0,0);
            for(int xCnt = 0; xCnt < xSample; xCnt ++)
                for(int yCnt = 0; yCnt < ySample; yCnt ++)
                {
                    cameraX = (Sx*(m+(xCnt+randX)/xSample)/M-Sx/2.0);
                    cameraY = (Sy*(n+(yCnt+randY)/ySample)/N-Sy/2.0);
                    pos[0] = pc[0] + n1[0]*cameraY + n0[0]*cameraX;
                    pos[1] = pc[1] + n1[1]*cameraY + n0[1]*cameraX;
                    pos[2] = pc[2] + n1[2]*cameraY + n0[2]*cameraX;
                    npe = (pe - pos)/norm(pe - pos);

                    rayTracing = colorDeter.spaceTracer(npe, pe);
                    Vec4f temp = colorDeter.colorReturn(-1, rayTracing.first, rayTracing.second, npe, 1) + colorDeter.projectColorDet(rayTracing.second);
                    colorSum += Vec3b(min(255, int(round(temp[0]*255/temp[3]))),min(255, int(round(temp[1]*255/temp[3]))),min(255, int(round(temp[2]*255/temp[3]))));
                }
            img.at<Vec3b>(n,m) = Vec3b((colorSum)/(xSample*ySample));

        }

//    flip(img, img, 1);
    imshow("lena", img);
    imwrite("img.png",img);
    resize(img,img,Size(160,120));
    imwrite("00.jpg",img);
    waitKey(0);
    return 0;
}
