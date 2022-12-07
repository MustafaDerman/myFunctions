#include "raisim/RaisimServer.hpp"
#include "D:/Raisim_v115/raisimLib-master/build/examples/TrajGen.h"

using namespace Eigen;

Matrix <double, 6, 1> CalcZMP(Vector3d LF1, Vector3d LF2, Vector3d LF3, Vector3d LF4, Vector3d RF1, Vector3d RF2, Vector3d RF3, Vector3d RF4) {

    Matrix <double, 6, 1> ZMP(0, 0, 0, 0, 0, 0);
    Vector3d LFoot_Force, RFoot_Force;
    LFoot_Force.setZero(), RFoot_Force.setZero();
    double LxZMP, LyZMP, RxZMP, RyZMP;
    LxZMP = 0;                          // Left foot x-ZMP
    LyZMP = 0;                          // Left foot y-ZMP
    RxZMP = 0;                          // Right foot x-ZMP
    RyZMP = 0;                          // Right foot y-ZMP
    Vector3d L3;
    L3 << 0.031, -0.0615, -0.1245;      // Foot center coordinates, used to calculate total ZMP
    Vector3d L_ZMP, R_ZMP;
    L_ZMP.setZero(), R_ZMP.setZero();

    LFoot_Force = LF1 + LF2 + LF3 + LF4;
    RFoot_Force = RF1 + RF2 + RF3 + RF4;
    if ((LFoot_Force(2) + RFoot_Force(2)) < 5) {
        L_ZMP << 0, 0, 0;
        R_ZMP << 0, 0, 0;
    }
    else {
        LxZMP = (0.131 * LF1(2) + 0.131 * LF2(2) - 0.131 * LF3(2) - 0.131 * LF4(2)) / LFoot_Force(2);
        LyZMP = (-0.055 * LF1(2) + 0.055 * LF2(2) - 0.055 * LF3(2) + 0.055 * LF4(2)) / LFoot_Force(2);
        L_ZMP << LxZMP, LyZMP, 0;

        RxZMP = (0.131 * RF1(2) + 0.131 * RF2(2) - 0.131 * RF3(2) - 0.131 * RF4(2)) / RFoot_Force(2);
        RyZMP = (0.055 * RF1(2) - 0.055 * RF2(2) + 0.055 * RF3(2) - 0.055 * RF4(2)) / RFoot_Force(2);
        R_ZMP << RxZMP, RyZMP, 0;
    }
    ZMP(0) = L_ZMP(0);
    ZMP(1) = L_ZMP(1);
    ZMP(2) = L_ZMP(2);
    ZMP(3) = R_ZMP(0);
    ZMP(4) = R_ZMP(1);
    ZMP(5) = R_ZMP(2);

    return ZMP;
}

Matrix <double, 6, 1> PosPID(double Kp, double Ki, double Kd,
    Matrix <double, 6, 1> err,
    Matrix <double, 6, 1> errPrev, Matrix <double, 6, 1> integral, double dt) {

    VectorXd Force(6);
    Force = Kp * err + Ki * integral + Kd * (err - errPrev) / dt;

    return Force;
}

Matrix <double, 6, 1> Squat(Vector3d Des, double n, int i, int stallTime, double RealTime, double dt) {

    VectorXd refPos(6);
    Vector3d ankleTraj, kneeTraj, hipTraj;
    ankleTraj = FuncPoly6th(RealTime, stallTime / 1000, stallTime / 1000 + 4.5,
        0, 0, 0,
        0, 0, 0, M_PI / 9);
    kneeTraj = FuncPoly5th(RealTime, stallTime / 1000, stallTime / 1000 + 4.0,
        0, 0, 0,
        M_PI / 4, 0, 0);
    hipTraj = FuncPoly6th(RealTime, stallTime / 1000, stallTime / 1000 + 5.0,
        0, 0, 0,
        M_PI / 6, 0, 0, M_PI / 3);

    if (i < stallTime) {
        refPos(0) = -Des(0); //Left Ankle
        refPos(1) = -Des(1); //Left Knee
        refPos(2) = -Des(2); //Left Hip
        refPos(3) = -Des(2); //Right Hip
        refPos(4) = -Des(1); //Right Knee
        refPos(5) = -Des(0); //Right Ankle
    }
    else {
        refPos(0) = -Des(0) - 0 * ankleTraj(0) - 0 * (10 * M_PI / 180) * sin(2 * M_PI * 0.2 * n);
        refPos(1) = -Des(1) + 0 * kneeTraj(0) + (15 * M_PI / 180) * sin(2 * M_PI * 0.2 * n);
        refPos(2) = -Des(2) - 0 * hipTraj(0) - (15 * M_PI / 180) * sin(2 * M_PI * 0.2 * n);
        refPos(3) = -Des(2) - 0 * hipTraj(0) - 0 * (10 * M_PI / 180) * sin(2 * M_PI * 0.2 * n);
        refPos(4) = -Des(1) + 0 * kneeTraj(0) + 0 * (15 * M_PI / 180) * sin(2 * M_PI * 0.2 * n);
        refPos(5) = -Des(0) - 0 * ankleTraj(0) - 0 * (10 * M_PI / 180) * sin(2 * M_PI * 0.2 * n);
    }

    return refPos;
}

Matrix <double, 6, 1> Walking(Vector3d DesLeft, Vector3d DesRight, int i, int stallTime) {

    VectorXd refPos(6);

    if (i < stallTime) {
        refPos(0) = DesLeft(0); //Left Ankle
        refPos(1) = DesLeft(1); //Left Knee
        refPos(2) = DesLeft(2); //Left Hip
        refPos(3) = DesRight(2); //Right Hip
        refPos(4) = DesRight(1); //Right Knee
        refPos(5) = DesRight(0); //Right Ankle
    }
    else {
        refPos(0) = DesLeft(0);
        refPos(1) = DesLeft(1);
        refPos(2) = DesLeft(2);
        refPos(3) = DesRight(2);
        refPos(4) = DesRight(1);
        refPos(5) = DesRight(0);
    }

    return refPos;

}