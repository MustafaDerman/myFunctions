#include "raisim/RaisimServer.hpp"

#define Grav 9.81

bool FuncInterval(double RealTime, double ta, double ts, double dt)
{
    double Lim1 = ta - dt * 0.25;
    double Lim2 = ts + dt * 0.25;
    if ((RealTime > Lim1) && (RealTime < Lim2)) // t = [Lim1 Lim2]
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool FuncGreater(double RealTime, double ta, double dt)
{
    double Lim1 = ta + dt * 0.25;
    if (RealTime > Lim1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Eigen::Vector3d FuncPoly5th(double RealTime, double t_start, double t_end,
    double Z0, double dZ0, double ddZ0,
    double Ze, double dZe, double ddZe)
{
    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel, Acc;
    Eigen::Vector3d trajOut;
    double n0 = Z0;
    double n1 = dZ0;
    double n2 = ddZ0 / 2;
    double n3 = -(20 * Z0 - 20 * Ze + 12 * dZ0 * tw + 8 * dZe * tw + 3 * ddZ0 * pow(tw, 2) - ddZe * pow(tw, 2)) / (2 * pow(tw, 3));
    double n4 = (30 * Z0 - 30 * Ze + 16 * dZ0 * tw + 14 * dZe * tw + 3 * ddZ0 * pow(tw, 2) - 2 * ddZe * pow(tw, 2)) / (2 * pow(tw, 4));
    double n5 = -(12 * Z0 - 12 * Ze + 6 * dZ0 * tw + 6 * dZe * tw + ddZ0 * pow(tw, 2) - ddZe * pow(tw, 2)) / (2 * pow(tw, 5));
    if (RealTime >= t_start && RealTime <= t_end)
    {
        Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3) + n4 * pow(Rt, 4) + n5 * pow(Rt, 5);
        Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2) + 4 * n4 * pow(Rt, 3) + 5 * n5 * pow(Rt, 4);
        Acc = 2 * n2 + 6 * n3 * Rt + 12 * n4 * pow(Rt, 2) + 20 * n5 * pow(Rt, 3);
    }
    else if (RealTime < t_start)
    {
        Pos = Z0;
        Vel = dZ0;
        Acc = ddZ0;
    }
    else if (RealTime > t_end)
    {
        Pos = Ze;
        Vel = dZe;
        Acc = ddZe;
    }
    trajOut << Pos, Vel, Acc;
    return trajOut;
}

Eigen::Vector3d FuncPoly6th(double RealTime, double t_start, double t_end,
    double Z0, double dZ0, double ddZ0,
    double Ze, double dZe, double ddZe, double Fh)
{
    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel, Acc;
    Eigen::Vector3d trajOut;
    double n0 = Z0;
    double n1 = dZ0;
    double n2 = ddZ0 / 2;
    double n3 = -(84 * Z0 - 128 * Fh + 44 * Ze + 32 * dZ0 * tw - 12 * dZe * tw + 5 * ddZ0 * pow(tw, 2) + ddZe * pow(tw, 2)) / (2 * pow(tw, 3));
    double n4 = (222 * Z0 - 384 * Fh + 162 * Ze + 76 * dZ0 * tw - 46 * dZe * tw + 9 * ddZ0 * pow(tw, 2) + 4 * ddZe * pow(tw, 2)) / (2 * pow(tw, 4));
    double n5 = -(204 * Z0 - 384 * Fh + 180 * Ze + 66 * dZ0 * tw - 54 * dZe * tw + 7 * ddZ0 * pow(tw, 2) + 5 * ddZe * pow(tw, 2)) / (2 * pow(tw, 5));
    double n6 = (32 * Z0 - 64 * Fh + 32 * Ze + 10 * dZ0 * tw - 10 * dZe * tw + ddZ0 * pow(tw, 2) + ddZe * pow(tw, 2)) / pow(tw, 6);
    if (RealTime >= t_start && RealTime <= t_end)
    {
        Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3) + n4 * pow(Rt, 4) + n5 * pow(Rt, 5) + n6 * pow(Rt, 6);
        Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2) + 4 * n4 * pow(Rt, 3) + 5 * n5 * pow(Rt, 4) + 6 * n6 * pow(Rt, 5);
        Acc = 2 * n2 + 6 * n3 * Rt + 12 * n4 * pow(Rt, 2) + 20 * n5 * pow(Rt, 3) + 30 * n6 * pow(Rt, 4);
    }
    else if (RealTime < t_start)
    {
        Pos = Z0;
        Vel = dZ0;
        Acc = ddZ0;
    }
    else if (RealTime > t_end)
    {
        Pos = Ze;
        Vel = dZe;
        Acc = ddZe;
    }
    trajOut << Pos, Vel, Acc;
    return trajOut;
}

Eigen::VectorXd zmpBasedTraj(double RealTime, double T_transfer, double Ts, double Td,
    int Nphase, double px, double Vx_mean, double Zc, double XOffset, double dt) {

    // Note: Nphse tek sayı olunca son yarım adımda problem oluyor.
    double t0 = 0.0; // Beginning of the Universe(or Simulation)
    double t1 = t0 + 3; // Halt
    double t2 = t1 + T_transfer; // First half step
    double tstart_walk = t2 + Td; // Double phase after half step
    double tend_walk = tstart_walk + (Ts + Td) * Nphase; // Walking end
    double t3 = tend_walk + T_transfer; // Last half step

    double w = sqrt(Grav / Zc); // Natural freq of equivalent pendulum

    Eigen::Vector3d X_Com_Traj, LFoot_xTraj, RFoot_xTraj, LFoot_yTraj, RFoot_yTraj, LFoot_zTraj, RFoot_zTraj;
    Eigen::VectorXd Out(9); // Returns the output values

    double COMx, dCOMx, ddCOMx, COPx, LFoot_xPos, LFoot_zPos, RFoot_xPos, RFoot_zPos;

    /* Single and double support phase calculation */
    int k, kx;
    double StepIndex = (RealTime - tstart_walk) / (Ts + Td);
    if (FuncInterval(RealTime, t0, tstart_walk, dt) == true)
    {
        k = 0;
        kx = 0;
    }
    else if (FuncGreater(RealTime, tend_walk, dt) == true)
    {
        k = Nphase - 1;
        kx = floor(k / 2);
    }
    else
    {
        k = floor(StepIndex);
        kx = floor(StepIndex / 2);
    }
    /* Single and double support phase calculation */

    // COM_x parameters
    double COM_xd = (px + Vx_mean * Ts / 2);
    double Cx0 = (2 * px - COM_xd);
    double dCx0 = w * (px - Cx0) * 1 / tanh(w * Ts / 2);
    double ddCx0 = w * w * (Cx0 - px);
    double d_Cx0 = COM_xd;
    double d_dCx0 = dCx0;

    double Kx = d_dCx0 + w * (COM_xd - px) * 1 / tanh(w * Td / 2);
    double Strx = 2 * Td * Kx;
    double Cxd0 = Cx0 + Strx / 2;

    /* CoM Trajectory Start */
    if (FuncInterval(RealTime, t0, t1, dt) == true) // Put robot on the ground
    {
        COMx = 0; dCOMx = 0; ddCOMx = 0;
    }
    else if (FuncInterval(RealTime, t1, t2, dt) == true) // Initialize CoM position (First half step)
    {
        X_Com_Traj = FuncPoly5th(RealTime, t1, t2, XOffset, 0, 0, COM_xd, dCx0, -ddCx0);
        COMx = X_Com_Traj(0); dCOMx = X_Com_Traj(1); ddCOMx = X_Com_Traj(2);
    }
    else if (FuncInterval(RealTime, t2, tstart_walk, dt) == true) // Initialite CoM trajectory (Just before walking starts)
    {
        X_Com_Traj = FuncPoly5th(RealTime, t2, tstart_walk, COM_xd, dCx0, -ddCx0, Cxd0, dCx0, ddCx0);
        COMx = X_Com_Traj(0); dCOMx = X_Com_Traj(1); ddCOMx = X_Com_Traj(2);
    }
    else if (FuncInterval(RealTime, tstart_walk, tend_walk, dt) == true) // Start walking
    {
        double Rt = RealTime - tstart_walk - (Ts + Td) * (k);
        if (FuncInterval(Rt, 0, Ts, dt) == true) // Single Support Phase
        {
            double s_T = Rt;
            COMx = (Cx0 - px) * cosh(w * s_T) + (dCx0 / w) * sinh(w * s_T) + px;
            dCOMx = w * (Cx0 - px) * sinh(w * s_T) + dCx0 * cosh(w * s_T);
            ddCOMx = w * w * (Cx0 - px) * cosh(w * s_T) + w * dCx0 * sinh(w * s_T);
            COMx = COMx + (k + 1) * Strx / 2;
        }
        else if (FuncInterval(Rt, Ts, Ts + Td, dt) == true) // Double Support Phase
        {
            double d_T = Rt - Ts;
            COMx = (d_Cx0 - px) * cosh(w * d_T) + ((d_dCx0 - Kx) / w) * sinh(w * d_T) + Kx * d_T + px;
            dCOMx = w * (d_Cx0 - px) * sinh(w * d_T) + (d_dCx0 - Kx) * cosh(w * d_T) + Kx;
            ddCOMx = w * w * (d_Cx0 - px) * cosh(w * d_T) + w * (d_dCx0 - Kx) * sinh(w * d_T);
            COMx = COMx + (k + 1) * Strx / 2;
        }
    }
    else if (FuncInterval(RealTime, tend_walk, t3, dt) == true) // Stop CoM trajectory(After walking ends)
    {
        double PosX_start = Cxd0 + (k + 1) * Strx / 2;
        double PosX_end = Strx / 2 + (k + 1) * Strx / 2;

        X_Com_Traj = FuncPoly5th(RealTime, tend_walk, t3, PosX_start, dCx0, ddCx0, PosX_end, 0, 0);
        COMx = X_Com_Traj(0); dCOMx = X_Com_Traj(1); ddCOMx = X_Com_Traj(2);
    }
    else
    {
        COMx = Strx / 2 + (k + 1) * Strx / 2; dCOMx = 0; ddCOMx = 0;
    }

    COPx = COMx - Zc * ddCOMx / Grav;
    /* CoM Trajectory End */

    /* Z-axis Feet Trajectory START */
    double Fh = 0.06; // 0.06
    double RFoot_zVel, RFoot_zAcc, LFoot_zVel, LFoot_zAcc;
    if (FuncInterval(RealTime, t0, t1, dt) == true) // Put robot on the ground
    {
        RFoot_zPos = 0; RFoot_zVel = 0; RFoot_zAcc = 0;
        LFoot_zPos = 0; LFoot_zVel = 0; LFoot_zAcc = 0;

    }
    else if (FuncInterval(RealTime, t1, t2, dt) == true) // Z - axis Foot trajectory initialization(First half step)
    {
        LFoot_zTraj = FuncPoly6th(RealTime, t1, t2, 0, 0, 0, 0, 0, 0, Fh);
        LFoot_zPos = LFoot_zTraj(0); LFoot_zVel = LFoot_zTraj(1); LFoot_zAcc = LFoot_zTraj(2);
        RFoot_zPos = 0; RFoot_zVel = 0; RFoot_zAcc = 0;

    }
    else if (FuncInterval(RealTime, t2, tstart_walk, dt) == true)
    {
        RFoot_zPos = 0; RFoot_zVel = 0; RFoot_zAcc = 0;
        LFoot_zPos = 0; LFoot_zVel = 0; LFoot_zAcc = 0;

    }
    else if (FuncInterval(RealTime, tstart_walk, tend_walk, dt) == true) // Z - axis Foot trajectory afater initialization
    {
        if (((k + 1) % 2) == 0) // Left foot swing, right foot stand
        {
            double ts_r = tstart_walk + (Ts + Td) * k;
            double ts_l = tstart_walk + (Ts + Td) * k;
            LFoot_zTraj = FuncPoly6th(RealTime, ts_l, ts_l + Ts, 0, 0, 0, 0, 0, 0, Fh);
            RFoot_zTraj = FuncPoly6th(RealTime, ts_r, ts_r + Ts, 0, 0, 0, 0, 0, 0, 0);
            LFoot_zPos = LFoot_zTraj(0); LFoot_zVel = LFoot_zTraj(1); LFoot_zAcc = LFoot_zTraj(2);
            RFoot_zPos = RFoot_zTraj(0); RFoot_zVel = RFoot_zTraj(1); RFoot_zAcc = RFoot_zTraj(2);
        }
        else // Right foot swing, left foot stand
        {
            double ts_r = tstart_walk + (Ts + Td) * k;
            double ts_l = tstart_walk + (Ts + Td) * k;
            LFoot_zTraj = FuncPoly6th(RealTime, ts_l, ts_l + Ts, 0, 0, 0, 0, 0, 0, 0);
            RFoot_zTraj = FuncPoly6th(RealTime, ts_r, ts_r + Ts, 0, 0, 0, 0, 0, 0, Fh);
            LFoot_zPos = LFoot_zTraj(0); LFoot_zVel = LFoot_zTraj(1); LFoot_zAcc = LFoot_zTraj(2);
            RFoot_zPos = RFoot_zTraj(0); RFoot_zVel = RFoot_zTraj(1); RFoot_zAcc = RFoot_zTraj(2);
        }
    }
    else if (FuncInterval(RealTime, tend_walk, t3, dt) == true)
    {
        LFoot_zPos = 0; LFoot_zVel = 0; LFoot_zAcc = 0;
        RFoot_zTraj = FuncPoly6th(RealTime, tend_walk, t3, 0, 0, 0, 0, 0, 0, Fh);
        RFoot_zPos = RFoot_zTraj(0); RFoot_zVel = RFoot_zTraj(1); RFoot_zAcc = RFoot_zTraj(2);
    }
    else
    {
        LFoot_zPos = 0; LFoot_zVel = 0; LFoot_zAcc = 0;
        RFoot_zPos = 0; RFoot_zVel = 0; RFoot_zAcc = 0;
    }
    /* Z-axis Feet Trajectory END */

    /* X-axis and Y-axis Feet Trajectory START */
    double RFoot_xVel, RFoot_xAcc, LFoot_xVel, LFoot_xAcc;
    if (FuncInterval(RealTime, t0, t1, dt) == true)
    {
        RFoot_xPos = 0; RFoot_xVel = 0; RFoot_xAcc = 0;
        LFoot_xPos = 0; LFoot_xVel = 0; LFoot_xAcc = 0;
    }
    else if (FuncInterval(RealTime, t1, t2, dt) == true) // X - axis Foot trajectory initialization(First half step)
    {
        LFoot_xTraj = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, Strx / 2, 0, 0);
        LFoot_xPos = LFoot_xTraj(0); LFoot_xVel = LFoot_xTraj(1); LFoot_xAcc = LFoot_xTraj(2);
        RFoot_xPos = 0; RFoot_xVel = 0; RFoot_xAcc = 0;
    }
    else if (FuncInterval(RealTime, t2, tstart_walk, dt) == true)
    {
        RFoot_xPos = 0; RFoot_xVel = 0; RFoot_xAcc = 0;
        LFoot_xPos = Strx / 2; LFoot_xVel = 0; LFoot_xAcc = 0;
    }
    else if (FuncInterval(RealTime, tstart_walk, tend_walk, dt) == true) // X - axis Foot trajectory after initialization
    {
        if (((k + 1) % 2) == 0) // Left foot swing, right foot stand
        {
            double ts_r = tstart_walk + (Ts + Td) * k;
            double ts_l = tstart_walk + (Ts + Td) * k;
            LFoot_xTraj = FuncPoly5th(RealTime, ts_l, ts_l + Ts, Strx / 2 + Strx * kx, 0, 0, Strx / 2 + Strx * (kx + 1), 0, 0);
            RFoot_xTraj = FuncPoly5th(RealTime, ts_r, ts_r + Ts, Strx * (kx + 1), 0, 0, Strx * (kx + 1), 0, 0);
            LFoot_xPos = LFoot_xTraj(0); LFoot_xVel = LFoot_xTraj(1); LFoot_xAcc = LFoot_xTraj(2);
            RFoot_xPos = RFoot_xTraj(0); RFoot_xVel = RFoot_xTraj(1); RFoot_xAcc = RFoot_xTraj(2);
        }
        else // Right foot swing, left foot stand
        {
            double ts_r = tstart_walk + (Ts + Td) * k;
            double ts_l = tstart_walk + (Ts + Td) * k;
            LFoot_xTraj = FuncPoly5th(RealTime, ts_l, ts_l + Ts, Strx / 2 + Strx * kx, 0, 0, Strx / 2 + Strx * kx, 0, 0);
            RFoot_xTraj = FuncPoly5th(RealTime, ts_r, ts_r + Ts, Strx * kx, 0, 0, Strx * (kx + 1), 0, 0);
            LFoot_xPos = LFoot_xTraj(0); LFoot_xVel = LFoot_xTraj(1); LFoot_xAcc = LFoot_xTraj(2);
            RFoot_xPos = RFoot_xTraj(0); RFoot_xVel = RFoot_xTraj(1); RFoot_xAcc = RFoot_xTraj(2);
        }
    }
    else if (FuncInterval(RealTime, tend_walk, t3, dt) == true)
    {
        LFoot_xPos = Strx / 2 + Strx * (kx + 1); LFoot_xVel = 0; LFoot_xAcc = 0;
        RFoot_xTraj = FuncPoly5th(RealTime, tend_walk, t3, Strx * (kx + 1), 0, 0, Strx / 2 + Strx * (kx + 1), 0, 0);
        RFoot_xPos = RFoot_xTraj(0); RFoot_xVel = RFoot_xTraj(1); RFoot_xAcc = RFoot_xTraj(2);
    }
    else
    {
        LFoot_xPos = Strx / 2 + Strx * (kx + 1); LFoot_xVel = 0; LFoot_xAcc = 0;
        RFoot_xPos = Strx / 2 + Strx * (kx + 1); RFoot_xVel = 0; RFoot_xAcc = 0;
    }
    /* X-axis and Y-axis Feet Trajectory END */

    Out << COMx, dCOMx, ddCOMx, RFoot_xPos, RFoot_zPos, LFoot_xPos, LFoot_zPos, COPx, k;
    return Out;
}