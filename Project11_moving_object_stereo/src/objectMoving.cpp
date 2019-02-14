#include "objectMoving.h"

Vec3d rotateFun(Vec3d n, Vec3d v, double theta){
    return cos(theta)*v + (1 - cos(theta))*(n.dot(v))*n + sin(theta)*v.cross(n);
}

ObjectMoving::ObjectMoving(vector<Quadrics> spheres){
    this -> spheres = spheres;
}

void ObjectMoving::rotateTowards(Vec3d dir, int idxS){
    Vec3d originalC;
    Vec3d nowC;
    Vec3d normVec;
    originalC = spheres[idxS].pc;
    nowC = spheres[idxS].pc + dir;
    spheres[idxS].pc = nowC;
    normVec[0] = dir[1];
    normVec[1] = -dir[0];
    normVec = normVec/norm(normVec);

    if(norm(dir) < 1e-6) return;
    double S = norm(dir);
    double theta = S/spheres[idxS].s0;
    spheres[idxS].planeA = rotateFun(normVec, spheres[idxS].planeA, theta);
    spheres[idxS].planeY = rotateFun(normVec, spheres[idxS].planeY, theta);
    spheres[idxS].planeX = rotateFun(normVec, spheres[idxS].planeX, theta);

    spheres[idxS].planeC = rotateFun(normVec, spheres[idxS].planeC - originalC, theta) + nowC;
    spheres[idxS].projPoint = rotateFun(normVec, spheres[idxS].projPoint - originalC, theta) + nowC;
}
