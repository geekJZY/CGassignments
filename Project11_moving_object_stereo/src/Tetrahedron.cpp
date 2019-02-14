#include "Tetrahedron.h"

int Tetrahedron::traceTimeComp(vector<double> traceVecTime){
    double timeMin = -1;
    int minIdx = -1;
    int idx = 0;
    while(traceVecTime[idx] <= 0){
        idx ++;
        if(idx == traceVecTime.size()) return -1;
    }
    timeMin = traceVecTime[idx];
    minIdx = idx;
    for(;idx < traceVecTime.size(); idx ++){
        if((traceVecTime[idx] < timeMin)&&(traceVecTime[idx] > 0)){
            timeMin = traceVecTime[idx];
            minIdx = idx;
        }
    }
    return minIdx;
}

Tetrahedron::Tetrahedron(vector<Vec3d> points, Vec4f colorH0_, Vec4f colorH1_, Vec4f colorH2_, Vec4f colorH3_, Quadrics quad_)
    :quad(quad_)
{
    vector<Vec3d> temp(3);

    ColorH0 = colorH0_;
    ColorH1 = colorH1_;
    ColorH2 = colorH2_;
    ColorH3 = colorH3_;

    temp[0] = points[0];
    temp[1] = points[1];
    temp[2] = points[2];
    vecPoints.push_back(temp);

    temp[0] = points[0];
    temp[1] = points[1];
    temp[2] = points[3];
    vecPoints.push_back(temp);

    temp[0] = points[0];
    temp[1] = points[2];
    temp[2] = points[3];
    vecPoints.push_back(temp);

    temp[0] = points[1];
    temp[1] = points[2];
    temp[2] = points[3];
    vecPoints.push_back(temp);

    for(int i = 0; i < 4; i ++){
        As.push_back((vecPoints[i][1]-vecPoints[i][0]).cross(vecPoints[i][2]-vecPoints[i][0]));
        norms.push_back(As[i]/norm(As[i]));
    }
}

colorUnite Tetrahedron::rayTrace(Vec3d pe, Vec3d nep){
    Vec3d A0, A1, A2;
    Vec3d ph;
    vector<double> ts(4);
    double temp;
    int timeIdx;
    ///First: get the intersection point with all planes
    for(int i = 0; i < 4; i ++){
        temp = As[i].dot(nep);
        if(fabs(temp)<1e-4){
            ts[i] = -1;
            continue;
        }
        ts[i] = As[i].dot(vecPoints[i][0]-pe)/temp;
    }

    ///Judge which intersection
    for(int i = 0; i < 4; i ++){
        timeIdx = traceTimeComp(ts);
        if(timeIdx < 0) break;
        ph = pe + ts[timeIdx] * nep;
        if(quad.checkIn(ph)){
            A0 = (ph - vecPoints[timeIdx][1]).cross(ph - vecPoints[timeIdx][2]);
            A1 = (ph - vecPoints[timeIdx][2]).cross(ph - vecPoints[timeIdx][0]);
            A2 = (ph - vecPoints[timeIdx][0]).cross(ph - vecPoints[timeIdx][1]);
            if(As[timeIdx].dot(A1) > 0 && As[timeIdx].dot(A2) > 0 && As[timeIdx].dot(A0) > 0)
                return colorUnite(ph, norms[timeIdx],ColorH0,ColorH1,ColorH2,ColorH3);
        }
        ts[timeIdx] = -1;
    }



    return colorUnite(Vec3f(0,0,0), Vec3f(0,0,0),ColorH0,ColorH1,ColorH2,ColorH3);
}
