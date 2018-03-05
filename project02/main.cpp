#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;



class Quadrics{
public:
    int a02,a12,a22,a21,a00;
    double s0, s1, s2;
    Vec3f n0, n1, n2;
    Vec3f pc;
    bool checkIn(Vec3f pos);
    double rayTracer(Vec3f, Vec3f);
    Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3f pc_);
};

Quadrics::Quadrics(int a02_, int a12_, int a22_, int a21_, int a00_, double s0_, double s1_, double s2_, Vec3f pc_){
    a02 = a02_;
    a12 = a12_;
    a22 = a22_;
    a21 = a21_;
    a00 = a00_;
    s0 = s0_;
    s1 = s1_;
    s2 = s2_;
    pc = pc_;
    n0 = Vec3f(1,0,0);
    n1 = Vec3f(0,1,0);
    n2 = Vec3f(0,0,1);
}

bool Quadrics::checkIn(Vec3f pos){
    double F = pow((pos-pc).dot(n0)/s0,2)*a02+pow((pos-pc).dot(n1)/s1,2)*a12+pow((pos-pc).dot(n2)/s2,2)*a22+a00;
    if(F<0) return true;
    else return false;
}

double Quadrics::rayTracer(Vec3f pe, Vec3f npe){
    double A = a02*pow(n0.dot(npe)/s0,2)+a12*pow(n1.dot(npe)/s1,2)+a22*pow(n2.dot(npe)/s2,2);
    double B = a02*(2*n0.dot(npe))*(n0.dot(pe-pc))/(s0*s0)+
               a12*(2*n1.dot(npe))*(n1.dot(pe-pc))/(s1*s1)+
               a22*(2*n2.dot(npe))*(n2.dot(pe-pc))/(s2*s2)+
               a21*n2.dot(npe)/s2;
    double C = a02*pow(n0.dot(pe-pc)/s0,2) +
               a12*pow(n1.dot(pe-pc)/s1,2) +
               a22*pow(n2.dot(pe-pc)/s2,2) +
               a21*n2.dot(pe-pc)/s2 +
               a00;
//    cout << "A is " << A <<" B is " << B <<" C is " << C << endl;
    double delta = pow(B,2)-4*A*C;
    if(delta < 0) return -1;
    else return ((-B-sqrt(pow(B,2)-4*A*C))/(2*A));
}

class Planes{
public:
    Vec3f pc;
    Vec3f n2;
    Planes(Vec3f n2_, Vec3f pc_);
    bool checkIn(Vec3f pos);
    double rayTracer(Vec3f pe, Vec3f npe);
};

Planes::Planes(Vec3f n2_, Vec3f pc_){
    n2 = n2_;
    pc = pc_;
}

bool Planes::checkIn(Vec3f pos){
    double F = n2.dot(pos - pc);
    if(F<0) return true;
    else return false;
}

double Planes::rayTracer(Vec3f pe, Vec3f npe){
    return n2.dot(pc-pe)/n2.dot(npe);
}

int main(int argc, char *argv[])
{
    vector<Scalar> color{Scalar(255,224,147),Scalar(93,66,255),Scalar(255,255,255),Scalar(230,180,80)};
    Mat img(Size(640,480),CV_8UC3,color[0]);
    double Sx = 16, Sy = 12;
    Vec3f pos;
    int xSample = 2, ySample = 2;
    vector<int> sample(5);
    double randX,randY;
    int M = img.cols, N = img.rows;

    //Graphics
    Quadrics ellipsoid(1,1,1,0,-1,8,8,4,Vec3f(0,0,0));
    Quadrics sphere(1,1,1,0,-1,5,5,5,Vec3f(0,0,5));
    Planes plane(Vec3f(-1,-1,0),Vec3f(10,10,0));
    cout << "n0 is " << ellipsoid.n0 << " n1 is " << ellipsoid.n1 << " n2 is " << ellipsoid.n2 << endl;

    //camera
    Vec3f pe(-36,-36,0);
    Vec3f vUp(1,1,3);
    Vec3f v2(-18,-18,0);
    Vec3f n0Vec(0,0,1);
    Vec3f n2 = v2/norm(v2);
//    Mat n1 = normlizeMat(vUp);
//    Mat n0 = mulMat(n2,n1);
    Vec3f n0 = n0Vec/norm(n0Vec);
    Vec3f n1 = n0.cross(n2);
    v2 = n2 * sqrt(550);
    Vec3f pcMat = pe + v2;
    Vec3f pc(pcMat[0], pcMat[1], pcMat[2]);
    Vec3f th_P;

    double cameraX,cameraY; //the coordinate in the camera
    Vec3f npe;

    cout << n0 << endl << n1 << endl << n2 << endl;

    //checkIn
    if(ellipsoid.checkIn(pe)||sphere.checkIn(pe)){
        cout << "fail in checking in " << endl;
        return 0;
    }
//    Mat npeTry = normlizeMat(Mat(Point3d(1,1,0)));
//    cout << "npeTry is "<< npeTry << endl;
//    cout << "TracePlane is " << plane.rayTracer(pe,npeTry) << endl;
//    cout << "TracePlane point is " << Mat(pe) + plane.rayTracer(pe,npeTry) * npeTry << endl;
    for(int m = 0; m < M; m ++)
        for(int n = 0; n < N; n ++)
        {
            //with anti-aliasing
            randX = rand()%1000/1000.0;
            randY = rand()%1000/1000.0;
            for(int i = 0; i < sample.size(); i ++)
            {
                sample[i] = 0;
            }
            for(int xCnt = 0; xCnt < xSample; xCnt ++)
                for(int yCnt = 0; yCnt < ySample; yCnt ++)
                {
                    cameraX = (Sx*(m+(xCnt+randX)/xSample)/M-Sx/2.0);
                    cameraY = (Sy*(n+(yCnt+randY)/ySample)/N-Sy/2.0);
                    pos[0] = pc[0] + n1[0]*cameraY + n0[0]*cameraX;
                    pos[1] = pc[1] + n1[1]*cameraY + n0[1]*cameraX;
                    pos[2] = pc[2] + n1[2]*cameraY + n0[2]*cameraX;
                    npe = (pe - pos)/norm(pe - pos);

//                    if(fabs(npe.at<double>(0)+1)<1e-6){
//                        cout << "pos is " << pos << endl;
//                        cout << "npe is " << npe << endl;
//                        cout << "rayTracer is " << ellipsoid.rayTracer(pe,npe,n0,n1,n2) << endl;
//                        while(!getchar());
//                    }
                    double tEllipsoid, tSphere;
                    tEllipsoid = ellipsoid.rayTracer(pe,npe);
                    tSphere = sphere.rayTracer(pe,npe);

                    if(tEllipsoid > 0 && tSphere < 0) sample[1]+=1;
                    else if(tEllipsoid < 0 && tSphere > 0) sample[2]+=1;
                    else if(tEllipsoid > 0 && tSphere > 0){
                        if(tEllipsoid > tSphere) sample[2]+=1;
                        else sample[1]+=1;
                    }
                    else{
                        double tplane = plane.rayTracer(pe,npe);
                        if(tplane > 0){
                            th_P = pe + tplane * npe;
                            if(th_P[2] > 0) sample[3] ++;
                            else sample[0] ++;
                        }
                    }
                }
            img.at<Vec3b>(n,m)[0] = (color[0].val[0] * sample[0] + color[1].val[0] * sample[1] + color[2].val[0] * sample[2] + color[3].val[0] * sample[3])/(xSample*ySample);
            img.at<Vec3b>(n,m)[1] = (color[0].val[1] * sample[0] + color[1].val[1] * sample[1] + color[2].val[1] * sample[2] + color[3].val[1] * sample[3])/(xSample*ySample);
            img.at<Vec3b>(n,m)[2] = (color[0].val[2] * sample[0] + color[1].val[2] * sample[1] + color[2].val[2] * sample[2] + color[3].val[2] * sample[3])/(xSample*ySample);


            //without anti-aliasing
//            pos.x = Sx*(m+0.5)/M + V0.x;
//            pos.y = Sy*(n+0.5)/N + V0.y;
//            pos.z = 10;
//
//            if(sphere1(pos))
//            {
//                img.at<Vec3b>(n,m)[0] = color[2].val[0];
//                img.at<Vec3b>(n,m)[1] = color[2].val[1];
//                img.at<Vec3b>(n,m)[2] = color[2].val[2];
//            }
//            else if(sphere2(pos))
//            {
//                img.at<Vec3b>(n,m)[0] = color[3].val[0];
//                img.at<Vec3b>(n,m)[1] = color[3].val[1];
//                img.at<Vec3b>(n,m)[2] = color[3].val[2];
//            }
//            else if(Ellipsoid(pos))
//            {
//                img.at<Vec3b>(n,m)[0] = color[1].val[0];
//                img.at<Vec3b>(n,m)[1] = color[1].val[1];
//                img.at<Vec3b>(n,m)[2] = color[1].val[2];
//            }
//            else if(plane(pos))
//            {
//                img.at<Vec3b>(n,m)[0] = color[4].val[0];
//                img.at<Vec3b>(n,m)[1] = color[4].val[1];
//                img.at<Vec3b>(n,m)[2] = color[4].val[2];
//            }
        }

    imshow("lena", img);
    imwrite("img.png",img);
    resize(img,img,Size(160,120));
    imwrite("00.jpg",img);
    waitKey(0);
    return 0;
}



