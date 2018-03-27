
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

#include "colorUn.h"
#include "lightSource.h"
#include "Quartics.h"
#include "Tetrahedron.h"
#include "Cube.h"

#define PI 3.1415927

using namespace cv;
using namespace std;

int main(int argc, char *argv[]){

    vector<Scalar> color{Scalar(255,224,147),Scalar(93,66,255),Scalar(255,255,255),Scalar(230,180,80)};
    Mat img(Size(640,480),CV_8UC3,color[0]);
    double Sx = 16, Sy = 12;
    Vec3d pos;
    Vec3d Ph;
    int objectIdx;
    int xSample = 2, ySample = 2;
    double randX,randY;
    int M = img.cols, N = img.rows;

    //the efficient of background and outline
    double kb = 0.5, k0 = 0.5;
    vector<float> shadowFlag;
    Quadrics quad(1,1,1,0,-1,sqrt(3),sqrt(3),sqrt(3),Vec3d(0,0,0),Vec4f(0,0,0,0), Vec4f(0,0,0,0), Vec4f(0,0,0,0),Vec4f(0,0,0,255));
    colorUnite colorIns;

    ///define Tetrahedron
    vector<Vec3d> points;
    points.push_back(Vec3d(1,-1,1));
    points.push_back(Vec3d(-1,-1,1));
    points.push_back(Vec3d(-1,1,1));
    points.push_back(Vec3d(1,1,1));
    points.push_back(Vec3d(1,1,-1));
    points.push_back(Vec3d(-1,1,-1));
    points.push_back(Vec3d(-1,-1,-1));
    points.push_back(Vec3d(1,-1,-1));
    Cube cube(points, Vec4f(53,36,155,255),Vec4f(113,86,255,255),Vec4f(1930,1660,2550,2550),Vec4f(0,0,0,255), quad);

    lightSource light(Vec3d(0,40,50), Vec3d(0.5773,-0.5773,-0.5773),Vec4f(60,60,60,60), 0.9,0.5,0.5,0.5, lightSource::SPOT);

    //camera
    Vec3d pe(20,10,10);
    Vec3d vUp(-1,0,1);
    Vec3d v2(60,30,30);
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

    cout << n0 << endl << n1 << endl << n2 << endl;

    //checkIn
    if(quad.checkIn(pe)){
        cout << "fail in checking in " << endl;
        return 0;
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

                    if(quad.rayTracer(pe, npe) > 0){
                        colorIns = cube.rayTrace(pe, npe);
                        if(norm(colorIns.norm) > 1e-5){
                            colorSum += light.colorDeter(pe, colorIns.Ph, colorIns.norm, colorIns.ColorH0, colorIns.ColorH1, colorIns.ColorH2, colorIns.ColorH3);
//                            cout << colorIns.norm << endl;
                        }
                    }

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


