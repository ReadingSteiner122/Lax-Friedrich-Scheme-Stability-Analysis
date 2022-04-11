#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <exception>
#include <functional>
#include <algorithm>

using namespace std;

const double a = 1.0;
double CFL = 1;
const int NumOfPnts = 501;
const double xL = -1.0, xR = 1.0;
const double dx = (xR - xL)/(NumOfPnts - 1);
vector<double> x(NumOfPnts, 0.0f);
vector<double> u_prev(NumOfPnts+2, 0.0f);
vector<double> u_cur(NumOfPnts+2, 0.f);
vector<double> reset(NumOfPnts+2, 0.f);

const double dt = CFL * dx/a;

const double NumOfSteps1 = 4/dt;
const int NumOfSteps = int(NumOfSteps1);
const double c = a * dt/dx;
const double L = xR - xL;

double u0(double x){
    //return 1.0 * exp(-8.0 * pow(x, 2));
    if(3*x>=1 || 3*x<=-1)
        return 0;
    else{
        if(x<0)
            return 1 + 3*x;
        else
            return 1 - 3*x;
    }
}

void write_u(ofstream &f, const vector<double> &x, const vector<double> &y)
{
    const int n = x.size()-1;

    f << x[1];
    for (int i = 0; i < n; i++)
        f  << y[i] << '\t'<< x[i]<<endl;
    f << endl;
}

void write_x(ofstream &f, const vector<double> &x)
{
    f << x[0];
    for(int k = 1; k < NumOfPnts; k++)
        f << '\t' << x[k];
    f << endl;
}
void exact()
{
    string str = "Exacta2var";
    str = str + to_string(NumOfPnts);
    str = str + ".dat";
    ofstream fout(str);
    u_prev.assign(reset.begin(), reset.end());

    fout << NumOfSteps << "\t" << NumOfPnts << endl;
    //write_x(fout, x);

    //IC
    for(int i = 1; i <= NumOfPnts; i++)
        u_prev[i] = u0(x[i-1]);
    
    //Periodical BC
    u_prev[0] = u_prev[NumOfPnts];
    u_prev[NumOfPnts+1] = u_prev[1];

    write_u(fout, u_prev, x);

    //Iterate over time
    double t = 0;
    for(int k = 1; k < NumOfSteps; k++)
    {
        t += dt;
        for(int i = 1; i <= NumOfPnts; i++)
        {
            double x_origin = x[i-1] - a * t;
            while (x_origin < xL)
                x_origin += L;
            while (x_origin > xR)
                x_origin -= L;
            u_cur[i] = u0(x_origin);
        }
        
        //Periodical BC
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];
        //write_u(fout, u_cur, x);
       
    }
    write_u(fout, u_cur, x);
    fout.close();
    cout << "Exact Soln Done!" << endl;
}

void expLaxW(){
    ofstream fout("Explicit-Lax-Wendroff-a2.dat");
    u_prev.assign(reset.begin(), reset.end());
    fout << NumOfSteps << "\t" <<NumOfPnts << endl;
    //write_x(fout, x);

    for(int i=1; i<=NumOfPnts; i++)
        u_prev[i] = u0(x[i-1]);
    
    u_prev[0] = u_prev[NumOfPnts];
    u_prev[NumOfPnts+1] = u_prev[1];

    //write_u(fout, u_prev);

    double  b_0 = 0.5 * CFL * (1 + CFL);
    double b_1 = -0.5 * CFL * (1 - CFL);
    double b_2 = 1 - pow(CFL, 2);

    for(int i=1; i<=NumOfSteps; i++){
        for(int j=1; j<=NumOfPnts; j++)
            u_cur[j] = b_0 * u_prev[j - 1] + b_1 * u_prev[j + 1] + b_2 * u_prev[j];
        
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];

        write_u(fout, u_cur, x);
        u_prev.swap(u_cur);
    }
    fout.close();
    cout<<"Solution according to Explicit-Lax-Wendroff Scheme Calculated!"<<endl;
}
void expLaxF(){
    u_prev.assign(reset.begin(), reset.end());
    ofstream fout("Explicit-Lax-Friedrichs-a2.dat");
    fout << NumOfSteps << "\t" <<NumOfPnts << endl;
    //write_x(fout, x);

    for(int i=1; i<=NumOfPnts; i++)
        u_prev[i] = u0(x[i-1]);
    
    u_prev[0] = u_prev[NumOfPnts];
    u_prev[NumOfPnts+1] = u_prev[1];

    //write_u(fout, u_prev);

    double  b_0 = 0.5 * (1 + CFL);
    double b_1 = 0.5 * (1 - CFL);

    for(int i=1; i<=NumOfSteps; i++){
        for(int j=1; j<=NumOfPnts; j++)
            u_cur[j] = b_0 * u_prev[j - 1] + b_1 * u_prev[j + 1];
        
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];

        //write_u(fout, u_cur, x);
        u_prev.swap(u_cur);
    }
    write_u(fout, u_cur, x);
    fout.close();
    cout<<"Solution according to Explicit-Lax-Friedrich Scheme Calculated!"<<endl;
}

void FTBS(){
    u_prev.assign(reset.begin(), reset.end());
    string str = "FTBS-Schemevar4cfl-";
    str = str + to_string(NumOfPnts) + to_string(CFL);
    str = str + ".dat";
    ofstream fout(str);
    fout << NumOfSteps << "\t" <<NumOfPnts << endl;
    //write_x(fout, x);

    for(int i=1; i<=NumOfPnts; i++)
        u_prev[i] = u0(x[i-1]);
    
    u_prev[0] = u_prev[NumOfPnts];
    u_prev[NumOfPnts+1] = u_prev[1];

    //write_u(fout, u_prev, x);

    //double  b_0 = 0.5 * (1 + CFL);
    //double b_1 = 0.5 * (1 - CFL);

    for(int i=1; i<=NumOfSteps; i++){
        for(int j=1; j<=NumOfPnts; j++)
            u_cur[j] = (1 - CFL) * u_prev[j] + CFL * u_prev[j-1];
        
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];

        //write_u(fout, u_cur, x);
        u_prev.assign(u_cur.begin(), u_cur.end());
    }
    write_u(fout, u_cur, x);
    fout.close();
    cout<<"Solution according to FTBS Scheme Calculated!"<<endl;
}

void FTCS(){
    u_prev.assign(reset.begin(), reset.end());
    string str = "FTCS-Schemevar4cfl-";
    str = str + to_string(NumOfPnts) + to_string(CFL);
    str = str + ".dat";
    ofstream fout(str);
    fout << NumOfSteps << "\t" <<NumOfPnts << endl;
    //write_x(fout, x);

    for(int i=1; i<=NumOfPnts; i++)
        u_prev[i] = u0(x[i-1]);
    
    u_prev[0] = u_prev[NumOfPnts];
    u_prev[NumOfPnts+1] = u_prev[1];

    //write_u(fout, u_prev, x);

    //double  b_0 = 0.5 * (1 + CFL);
    //double b_1 = 0.5 * (1 - CFL);

    for(int i=1; i<=NumOfSteps; i++){
        for(int j=1; j<=NumOfPnts; j++)
            u_cur[j] = u_prev[j] - 0.5 * CFL * (u_prev[j+1] - u_prev[j-1]);
        
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];

        //write_u(fout, u_cur, x);
        u_prev.swap(u_cur);
    }
    write_u(fout, u_cur, x);
    fout.close();
    cout<<"Solution according to FTCS Scheme Calculated!"<<endl;
}

void LeapFrog(){
    u_prev.assign(reset.begin(), reset.end());
    string str = "Leapfrogvar4cfl-";
    str = str + to_string(NumOfPnts) + to_string(CFL);
    str = str + ".dat";
    ofstream fout(str);
    vector<double> u_prev_v(NumOfPnts+2, 0.f);
    fout << NumOfSteps << "\t" << NumOfPnts << endl;
    //write_x(fout, x);
    for(int i=1; i<=NumOfPnts; i++)
        u_prev_v[i] = u0(x[i-1]);

    u_prev_v[0] = u_prev_v[NumOfPnts];
    u_prev_v[NumOfPnts+1] = u_prev_v[1];
    //write_u(fout, u_prev_v, x);
    //Calculate u at t=del(t) using the FTBS Scheme
    for(int j=1; j<=NumOfPnts; j++)
            u_prev[j] = (1 - CFL) * u_prev_v[j] + CFL * u_prev_v[j-1];
    

    //write_u(fout, u_prev_v, x);
    //write_u(fout, u_prev, x);

    for(int i=1; i<=NumOfSteps-1; i++){
        for(int j=1; j<=NumOfPnts; j++)
            u_cur[j] = u_prev_v[j] - CFL * (u_prev[j+1] - u_prev[j-1]);
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];

        //write_u(fout, u_cur, x);

        u_prev_v.assign(u_prev.begin(), u_prev.end());
        u_prev.assign(u_cur.begin(), u_cur.end());
    }
    write_u(fout, u_cur, x);
    fout.close();
    cout<<"Solution according to Leapfrog Scheme Calculated!"<<endl;

}

int main(){
    cout<<"Please Enter a CFL Number: "<<endl;
    cin>>CFL;
    cout<<"Initialized Grid"<<endl;
    x[0] = xL;
    for(int i = 1; i<NumOfPnts; i++)
        x[i] = x[i-1] + dx;
    //cout<<"Done"<<endl;
    //expLaxF();
    //cout<<"Explicit Lax-Friedrichs Calculated"<<endl;
    exact();
    double arr_CFL[] = {0.5, 1, 3};
    //expLaxW();
    for(int i=0; i<3; i++){
        CFL = arr_CFL[i];
        FTBS();
        FTCS();
        LeapFrog();
    }
    return 0;
}