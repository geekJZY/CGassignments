#include "Cube.h"
#include "Tetrahedron.h"

Cube::Cube(vector<Vec3d> points, Vec4f colorH0_, Vec4f colorH1_, Vec4f colorH2_, Vec4f colorH3_,int texture_, float nRR_, Quadrics quad_, Mat textureMap_, Mat normalMap_)
    : quad(quad_),
    textureMap(textureMap_),
    texture(texture_),
    nRR(nRR_),
    normalMap(normalMap_)
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

double Cube::rayTrace(Vec3d pe, Vec3d nep){
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
                int xI, yI;
                vector<Vec4f> texColors(3);
                pointSave = ph;
                xI = min(textureMap.cols-1, int(round(u * textureMap.cols)));
                yI = min(textureMap.rows-1, int(round(v * textureMap.rows)));

                texColors[1][0] = textureMap.at<Vec3b>(yI, xI)[0];
                texColors[1][1] = textureMap.at<Vec3b>(yI, xI)[1];
                texColors[1][2] = textureMap.at<Vec3b>(yI, xI)[2];
                texColors[1][3] = 255;

                for(int i = 0; i < 3; i ++){
                    texColors[0][i] = texColors[1][i]/5;
                    texColors[2][i] = min(texColors[1][i], float(255.0));
                    texColors[2][i] *= 10;
                }
                texColors[0][3] = 255;
                texColors[2][3] = 2550;
                colorSave = texColors;

                xI = min(normalMap.cols-1, int(round(u * normalMap.cols)));
                yI = min(normalMap.rows-1, int(round(v * normalMap.rows)));
                Vec3d n0_R, n1_G, n2_B;
                float a,b,c;
//                n2_B = norms[timeIdx];
//                n0_R = (vecPoints[timeIdx][1]-vecPoints[timeIdx][0])/norm(vecPoints[timeIdx][1]-vecPoints[timeIdx][0]);
//                n1_G = (vecPoints[timeIdx][3]-vecPoints[timeIdx][0])/norm(vecPoints[timeIdx][3]-vecPoints[timeIdx][0]);
//                a = 2*normalMap.at<Vec3b>(yI, xI)[2]/255.0 - 1;
//                b = 2*normalMap.at<Vec3b>(yI, xI)[1]/255.0 - 1;
//                c = 2*normalMap.at<Vec3b>(yI, xI)[0]/255.0 - 1;
//                normalSave = a*n0_R + b*n1_G + c*n2_B;
                normalSave = norms[timeIdx];

                return ts[timeIdx];
            }
        }
        ts[timeIdx] = -1;
    }
    return -1;
}

int Cube::traceTimeComp(vector<double> traceVecTime){
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
