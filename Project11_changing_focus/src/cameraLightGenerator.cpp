#include "cameraLightGenerator.h"

lightShooter::lightShooter(Vec3d pe,Vec3d v2, float f_){
    Peye = pe;
    n2 = v2/norm(v2);
    f = f_;
    float r = 0.38;

    //generate n0 n1
    Vec3d n0 = Vec3d(1,0,0);
    if(n0.dot(n2)>0.99) n0 = Vec3d(0,1,0);
    n0 = n0 - n0.dot(n2)*n2;
    n0 = n0 / norm(n0);
    Vec3d n1 = n0.cross(n2);

    //generate points around pe
    pointsAroundC.push_back(pe+n0*r*sin(PI/4)+n1*r*cos(PI/4));
    pointsAroundC.push_back(pe+n0*r*sin(PI/2)+n1*r*cos(PI/2));
    pointsAroundC.push_back(pe+n0*r*sin(PI/4*3)+n1*r*cos(PI/4*3));
    pointsAroundC.push_back(pe+n0*r*sin(PI/4*4)+n1*r*cos(PI/4*4));
    pointsAroundC.push_back(pe+n0*r*sin(PI/4*5)+n1*r*cos(PI/4*5));
    pointsAroundC.push_back(pe+n0*r*sin(PI/4*6)+n1*r*cos(PI/4*6));
    pointsAroundC.push_back(pe+n0*r*sin(PI/4*7)+n1*r*cos(PI/4*7));
    pointsAroundC.push_back(pe+n0*r*sin(PI/4*8)+n1*r*cos(PI/4*8));
}

pair<vector<Vec3d>, vector<Vec3d>> lightShooter::lightGene(Vec3d pos){
    vector<Vec3d> rayToShoot(pointsAroundC.size());

    //calculate the point which it would cross
    Vec3d ray_N = (Peye - pos)/norm(Peye - pos);
    Vec3d A = this -> n2;
    double S1 = -(Peye - pos).dot(A);
    double S2 = 1/(1/this->f - 1/S1);
//    cout << "S1 is " << S1 << " S2 is " << S2 << endl;
//    cout << "S2 should be " << (Peye - Vec3d(20,18,3)).dot(n2) << "f should be " << 1/(1/((Peye - Vec3d(20,18,3)).dot(n2))+1/S1) << endl;
//    getchar();
    Vec3d b = this->Peye - S2*n2;
    double t = A.dot(b - Peye)/A.dot(ray_N);
    Vec3d pCrs = t*ray_N + this->Peye;
    for(int i = 0; i < pointsAroundC.size(); i ++){
        rayToShoot[i] = (pCrs - pointsAroundC[i])/norm(pCrs - pointsAroundC[i]);
    }


    return make_pair(pointsAroundC, rayToShoot);
}
