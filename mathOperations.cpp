// This file contains the mathematical operations needed in a robotic controller project.
//

/* Differentiate Function */
double diff(double In, double prevIn, double dt)
{
    double Out = (In - prevIn)/dt; // derivative of the input signal
    return Out;
}
/* Differentiate Function */

/* Integral Function (Numerical Integral Calculation) */
double numIntegral(double In, double prevIn, double prevOut, double dt)
{    
    double Out = prevOut + (In + prevIn)*dt/2; // integral of the input signal
    return Out;
}
/* Integral Function (Numerical Integral Calculation) */

/* Low Pass Filter Function */
double LPF(double In, double PrevOut, double Freq, double dt)
{
    return (In*Freq+ (PrevOut/dt))/(Freq + (1/dt));
}
/* Low Pass Filter Function */

/* High Pass Filter Function */
double HPF(double In,double PrevIn, double PrevOut, double Freq, double dt)
{
    return (1/(1+(dt/Freq)))*(PrevOut+ (In- PrevIn));
}
/* High Pass Filter Function */

/* Unwrap function for a circular motion (i.e. rotary encoder) */
double angleConv(double angle){
    return fmod(angle, 2*M_PI);
}
double angleDiff(double a, double b){
    double dif = fmod(b - a + M_PI, 2*M_PI);
    if (dif < 0)
        dif += 2*M_PI;
    return dif - M_PI;
}
double unwrap(double previous_angle, double new_angle){
    return previous_angle - angleDiff(new_angle, angleConv(previous_angle));
}
/* Unwrap function for a circular motion (i.e. rotary encoder) */
