
#include "raisim/RaisimServer.hpp"
#include "Kinematics.h"
#include "robotDynamics.h"

#define Grav 9.81

Eigen::Vector3d computeNewtEul1Leg(Eigen::Vector3d W0, Eigen::Vector3d dW0, Eigen::Vector3d dV0,
    Eigen::Matrix3d tempOrient, Eigen::Vector3d grForce, Eigen::Vector3d ZMP,
    Eigen::Vector3d legjPos, Eigen::Vector3d legjVel, Eigen::Vector3d legjAcc, int legIndex) {

    Eigen::Vector3d legjTorques(0, 0, 0);
    Eigen::Matrix3d R01, R12, R23, RG3, I1, I2, I3;
    Eigen::Vector3d L0, L1, L2, L3, LC1, LC2, LC3;
    R01 = RotatePitch(legjPos(0)); // Coordinate transformation from 2(KneeJoint) to 1(HipJoint)
    R12 = RotatePitch(legjPos(1)); // Coordinate transformation from 3(AnkleJoint) to 2(KneeJoint)
    R23 = RotatePitch(legjPos(2)); // Coordinate transformation from End-Eff(Foot) to 3(AnkleJoint)

    switch (legIndex) {
    case 1: // Left_Leg
        L0 << 0.1815, 0.2395, 0.007877;
        L1 << 0, 0, -0.455;
        L2 << 0, 0, -0.395;
        // L3 + Offset to the Center of Foot
        L3 << 0.031, -0.0615, -0.117;
        // L3 + ZMP
        L3 = L3 + ZMP;
        LC1 << 0, 0, -0.2275;
        LC2 << 0, 0, -0.1975;
        LC3 << 0.027568, -0.054691, -0.0498;
        // <inertia ixx = "0.003" ixy = "3.297e-4" ixz = "3.002e-4" iyy = "0.01" iyz = "-5.956e-4" izz = "0.011" / >
        I3 << 0.003, 3.297e-4, 3.002e-4,
            3.297e-4, 0.01, -5.956e-4,
            3.002e-4, -5.956e-4, 0.011;
        break;
    case 2: // Right_Leg
        L0 << 0.1815, -0.2395, 0.007877;
        L1 << 0, 0, -0.455;
        L2 << 0, 0, -0.395;
        // L3 + Offset to the Center of Foot
        L3 << 0.031, 0.0615, -0.117;
        // L3 + ZMP
        L3 = L3 + ZMP;
        LC1 << 0, 0, -0.2275;
        LC2 << 0, 0, -0.1975;
        LC3 << 0.027568, 0.054691, -0.0498;
        //<inertia ixx = "0.003" ixy = "-3.297e-4" ixz = "3.002e-4" iyy = "0.01" iyz = "5.956e-4" izz = "0.011" / >
        I3 << 0.003, -3.297e-4, 3.002e-4,
            -3.297e-4, 0.01, 5.956e-4,
            3.002e-4, 5.956e-4, 0.011;
        break;
    }
#define m1 6
#define m2 5
#define m3 2
    //<inertia ixx = "0.052" ixy = "0.004" ixz = "0.012" iyy = "0.057" iyz = "-0.005" izz = "0.016" / >
    I1 << 0.052, 0.004, 0.012,
        0.004, 0.057, -0.005,
        0.012, -0.005, 0.016;
    //<inertia ixx = "0.005" ixy = "0.001" ixz = "-0.0008785" iyy = "0.01" iyz = "0.002" izz = "0.011" / >
    I2 << 0.005, 0.001, -0.0008785,
        0.001, 0.01, 0.002,
        -0.0008785, 0.002, 0.011;
    //<inertia ixx = "0.0009458" ixy = "0" ixz = "0" iyy = "0.001" iyz = "0.0001202 " izz = "0.0002289" / >

    Eigen::Vector3d JVel1, JVel2, JVel3;
    JVel1 << 0, legjVel(0), 0;
    JVel2 << 0, legjVel(1), 0;
    JVel3 << 0, legjVel(2), 0;
    Eigen::Vector3d JAcc1, JAcc2, JAcc3;
    JAcc1 << 0, legjAcc(0), 0;
    JAcc2 << 0, legjAcc(1), 0;
    JAcc3 << 0, legjAcc(2), 0;
    Eigen::Vector3d W1, dW1, dV1, dVc1, F1, N1, f1, n1;
    Eigen::Vector3d W2, dW2, dV2, dVc2, F2, N2, f2, n2;
    Eigen::Vector3d W3, dW3, dV3, dVc3, F3, N3, f3, n3;
    Eigen::Vector3d n4(0, 0, 0), f4(0, 0, 0);
    Eigen::Vector3d Temp(0, 0, 0);
    Eigen::Matrix3d Eye3;
    Eye3 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    RG3 = tempOrient * R01 * R12 * R23;
    f4 = RG3.transpose() * -grForce;

    //from i=0 to i=(N-1); N=3(# of Bodies) -> Outward Iteration
    // Outward (0 -> 2)
// i = 0 (Floating-Base to Hip-Joint)
// Angular Velocity and Angular Acceleration
    Temp = R01.transpose() * W0;
    W1 = Temp + JVel1;
    dW1 = R01.transpose() * dW0 + Temp.cross(JVel1) + JAcc1;
    // Translational Acceleration
    Temp = W0.cross(L0);
    Temp = W0.cross(Temp) + dW0.cross(L0);
    dV1 = R01.transpose() * (Temp + dV0);
    // Translational Acceleration of CoM
    Temp = W1.cross(LC1);
    Temp = W1.cross(Temp);
    dVc1 = dW1.cross(LC1) + dV1 + Temp;
    // Force and Torque acting on Link due to Inertial Forces
    F1 = m1 * dVc1;
    Temp = I1 * W1;
    Temp = W1.cross(Temp);
    N1 = I1 * dW1 + Temp;

    // i = 1 (Hip-Joint to Knee-Joint)
    // Angular Velocity and Angular Acceleration
    Temp = R12.transpose() * W1;
    W2 = Temp + JVel2;
    dW2 = R12.transpose() * dW1 + Temp.cross(JVel2) + JAcc2;
    // Translational Acceleration
    Temp = W1.cross(L1);
    Temp = W1.cross(Temp) + dW1.cross(L1);
    dV2 = R12.transpose() * (Temp + dV1);
    // Translational Acceleration of CoM
    Temp = W2.cross(LC2);
    Temp = W2.cross(Temp);
    dVc2 = dW2.cross(LC2) + dV2 + Temp;
    // Force and Torque acting on Link due to Inertial Forces
    F2 = m2 * dVc2;
    Temp = I2 * W2;
    Temp = W2.cross(Temp);
    N2 = I2 * dW2 + Temp;

    // i = 2 (Knee-Joint to Ankle-Joint)
    // Angular Velocity and Angular Acceleration
    Temp = R23.transpose() * W2;
    W3 = Temp + JVel3;
    dW3 = R23.transpose() * dW2 + Temp.cross(JVel3) + JAcc3;
    // Translational Acceleration
    Temp = W2.cross(L2);
    Temp = W2.cross(Temp) + dW2.cross(L2);
    dV3 = R23.transpose() * (Temp + dV2);
    // Translational Acceleration of CoM
    Temp = W3.cross(LC3);
    Temp = W3.cross(Temp);
    dVc3 = dW3.cross(LC3) + dV3 + Temp;
    // Force and Torque acting on Link due to Inertial Forces
    F3 = m3 * dVc3;
    Temp = I3 * W3;
    Temp = W3.cross(Temp);
    N3 = I3 * dW3 + Temp;

    // from i=N=3 to (i=1); N = 3 -> Invard Iteration
    // Inward (3 -> 1)
// i = 3 (Ankle-Joint)
    f3 = Eye3 * f4 + F3;
    Temp = Eye3 * f4;
    Temp = L3.cross(Temp);
    Temp = LC3.cross(F3) + Temp;
    n3 = N3 + Eye3 * n4 + Temp;

    // i = 2 (Knee-Joint)
    f2 = R23 * f3 + F2;
    Temp = R23 * f3;
    Temp = L2.cross(Temp);
    Temp = LC2.cross(F2) + Temp;
    n2 = N2 + R23 * n3 + Temp;

    // i = 1 (Hip-Joint)
    f1 = R12 * f2 + F1;
    Temp = R12 * f2;
    Temp = L1.cross(Temp);
    Temp = LC1.cross(F1) + Temp;
    n1 = N1 + R12 * n2 + Temp;

    legjTorques(2) = n3(1);
    legjTorques(1) = n2(1);
    legjTorques(0) = n1(1);

    return legjTorques;
}

Eigen::Matrix <double, 6, 1> computeNewtEul2Legs(Eigen::Vector3d rootAbsacc, Eigen::Matrix3d rootOrient,
    Eigen::Vector3d rootAngvel, Eigen::Vector3d rootAngacc,
    Eigen::Vector3d LFoot_Force, Eigen::Vector3d RFoot_Force, Eigen::Vector3d L_ZMP, Eigen::Vector3d R_ZMP,
    Eigen::Matrix <double, 6, 1> jointPos, Eigen::Matrix <double, 6, 1> jointVel, Eigen::Matrix <double, 6, 1> jointAcc) {

    Eigen::Matrix <double, 6, 1> ffTorques(0, 0, 0, 0, 0, 0);
    Eigen::Matrix3d trarootOrient, Eye3;
    Eigen::Vector3d W0, dW0, dV0;
    Eigen::Vector3d tempLegjTorques, tempLegjPos, tempLegjVel, tempLegjAcc;
    Eye3 << 1, 0, 0, 0, 1, 0, 0, 0, 1;

    rootAbsacc(2) = rootAbsacc(2) + Grav;

    // Leg 1: LF
    tempLegjPos << jointPos(0), jointPos(1), jointPos(2);
    tempLegjVel << jointVel(0), jointVel(1), jointVel(2);
    tempLegjAcc << jointAcc(0), jointAcc(1), jointAcc(2);
    trarootOrient = rootOrient.transpose();
    dV0 = trarootOrient * rootAbsacc;
    W0 = trarootOrient * rootAngvel;
    dW0 = trarootOrient * rootAngacc;
    tempLegjTorques = computeNewtEul1Leg(W0, dW0, dV0, rootOrient, LFoot_Force, L_ZMP,
        tempLegjPos, tempLegjVel, tempLegjAcc, 1);
    ffTorques(0) = tempLegjTorques(0);
    ffTorques(1) = tempLegjTorques(1);
    ffTorques(2) = tempLegjTorques(2);

    // Leg 2: RF
    tempLegjPos << jointPos(3), jointPos(4), jointPos(5);
    tempLegjVel << jointVel(3), jointVel(4), jointVel(5);
    tempLegjAcc << jointAcc(3), jointAcc(4), jointAcc(5);
    trarootOrient = rootOrient.transpose();
    dV0 = trarootOrient * rootAbsacc;
    W0 = trarootOrient * rootAngvel;
    dW0 = trarootOrient * rootAngacc;
    tempLegjTorques = computeNewtEul1Leg(W0, dW0, dV0, rootOrient, RFoot_Force, R_ZMP,
        tempLegjPos, tempLegjVel, tempLegjAcc, 2);
    ffTorques(3) = tempLegjTorques(0);
    ffTorques(4) = tempLegjTorques(1);
    ffTorques(5) = tempLegjTorques(2);

    return ffTorques;
}