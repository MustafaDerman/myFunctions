#include "raisim/RaisimServer.hpp"
#include "math.h"

Eigen::Matrix3d RotatePitch(double a) {

    double RotY[3][3];
    Eigen::Matrix3d pitchMatrix;
    RotY[0][0] = cos(a);
    RotY[0][1] = 0;
    RotY[0][2] = sin(a);

    RotY[1][0] = 0;
    RotY[1][1] = 1;
    RotY[1][2] = 0;

    RotY[2][0] = -sin(a);
    RotY[2][1] = 0;
    RotY[2][2] = cos(a);

    pitchMatrix << RotY[0][0], RotY[0][1], RotY[0][2], RotY[1][0], RotY[1][1], RotY[1][2], RotY[2][0], RotY[2][1], RotY[2][2];

    return pitchMatrix;
}

Eigen::Vector2d FoKinThreeDof(double q1, double q2, double q3) {
    Eigen::Vector2d Out;
    double L0 = 0.117;
    double L1 = 0.395;
    double L2 = 0.455;
    double L3 = 0.1815;
    double L4 = 0.007877;

    double Xe = L3 - L1 * sin(q2 + q3) - L2 * sin(q3) - L0 * sin(q1 + q2 + q3);
    double Ze = L4 - L1 * cos(q2 + q3) - L2 * cos(q3) - L0 * cos(q1 + q2 + q3);

    Out << Xe, Ze;
    return Out;
}

Eigen::Vector3d InKinThreeDofFromHip(double Xe, double Ze) {
    Eigen::VectorXd OutInv(3);
    double L0 = 0.117; // Foot base to ankle
    double L1 = 0.395; // Ankle to Knee
    double L2 = 0.455; // Knee to Hip
    double L3 = 0*0.1815; // Hip to Root on X
    double L4 = 0.007877; // Hip to Root on Z
    double beta = -0.174533/3;

    double c2 = ((Xe - L3) * (Xe - L3) + (Ze - L4 + L0) * (Ze - L4 + L0) - L1 * L1 - L2 * L2) / (2 * L1 * L2);
    double s2 = sqrt(1 - c2 * c2);

    double q2 = atan2(s2, c2);

    double A = -L1 * c2 - L2;
    double B = -L1 * s2;

    double c3 = (A * (Ze - L4 + L0) + B * (Xe - L3)) / (A * A + B * B);
    double s3 = -sqrt(1 - c3 * c3);
    double q3 = atan2(s3, c3);

    double q1 = beta - q2 - q3;

    OutInv << q1, q2, q3;
    return OutInv;
}