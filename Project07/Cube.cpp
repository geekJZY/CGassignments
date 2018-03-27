#include "Cube.h"
#include "Tetrahedron.h"

Cube::Cube(vector<Vec3d> points, Vec4f colorH0_, Vec4f colorH1_, Vec4f colorH2_, Vec4f colorH3_, Quadrics quad_)
    : quad(quad_)
{
    vector<Vec3d> temp(4);

    ColorH0 = colorH0_;
    ColorH1 = colorH1_;
    ColorH2 = colorH2_;
    ColorH3 = colorH3_;

    temp[0] = points[0];
    temp[1] = points[3];
    temp[2] = points[2];
    temp[3] = points[1];
    vecPoints.push_back(temp);

    temp[0] = points[0];
    temp[1] = points[7];
    temp[2] = points[4];
    temp[3] = points[3];
    vecPoints.push_back(temp);

    temp[0] = points[0];
    temp[1] = points[1];
    temp[2] = points[6];
    temp[3] = points[7];
    vecPoints.push_back(temp);

    temp[0] = points[5];
    temp[1] = points[4];
    temp[2] = points[7];
    temp[3] = points[6];
    vecPoints.push_back(temp);

    temp[0] = points[5];
    temp[1] = points[2];
    temp[2] = points[3];
    temp[3] = points[4];
    vecPoints.push_back(temp);

    temp[0] = points[5];
    temp[1] = points[6];
    temp[2] = points[1];
    temp[3] = points[2];
    vecPoints.push_back(temp);

    for(int i = 0; i < 6; i ++){
        As.push_back((vecPoints[i][1]-vecPoints[i][0]).cross(vecPoints[i][3]-vecPoints[i][0]));
        norms.push_back(As[i]/norm(As[i]));
    }
}

colorUnite Cube::rayTrace(Vec3d pe, Vec3d nep){
    double u,v;
    Vec3d ph;
    vector<double> ts(6);
    double temp;
    int timeIdx;
    ///First: get the intersection point with all planes
    for(int i = 0; i < 6; i ++){
        temp = As[i].dot(nep);
        if(fabs(temp)<1e-4){
            ts[i] = -1;
            continue;
        }
        ts[i] = As[i].dot(vecPoints[i][0]-pe)/temp;
    }

    ///Judge which intersection
    for(int i = 0; i < 6; i ++){
        timeIdx = traceTimeComp(ts);
        if(timeIdx < 0) break;
        ph = pe + ts[timeIdx] * nep;
        if(quad.checkIn(ph)){
            u = (ph - vecPoints[timeIdx][0]).dot(vecPoints[timeIdx][1]-vecPoints[timeIdx][0])/norm(vecPoints[timeIdx][1]-vecPoints[timeIdx][0])/norm(vecPoints[timeIdx][1]-vecPoints[timeIdx][0]);
            v = (ph - vecPoints[timeIdx][0]).dot(vecPoints[timeIdx][3]-vecPoints[timeIdx][0])/norm(vecPoints[timeIdx][3]-vecPoints[timeIdx][0])/norm(vecPoints[timeIdx][3]-vecPoints[timeIdx][0]);
            if(u >= 0 && u <= 1 && v >= 0 && v <= 1){
//                cout << " u is " << u << "v is " << v << endl;
                return colorUnite(ph, norms[timeIdx],ColorH0,ColorH1,ColorH2,ColorH3);
            }
        }
        ts[timeIdx] = -1;
    }
    return colorUnite(Vec3f(0,0,0), Vec3f(0,0,0),ColorH0,ColorH1,ColorH2,ColorH3);
}
