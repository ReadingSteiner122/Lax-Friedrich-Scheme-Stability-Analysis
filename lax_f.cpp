#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <exception>
#include <functional>
#include <algorithm>

using namespace std;

const double a = 1.0;
double CFL = 1;
const int NumOfPnts = 101;
const double xL = -1.0, xR = 1.0;
const double dx = (xR - xL)/(NumOfPnts - 1);
vector<double> x(NumOfPnts, 0.0f);
vector<double> u_prev(NumOfPnts+2, 0.0f);
vector<double> u_cur(NumOfPnts+2, 0.f);

const int NumOfSteps = 10000;
const double dt = CFL * dx/a;
const double c = a * dt/dx;


double u0(double x){
    return exp(-0.8 * pow(x, 2));
}

void write_u(ofstream &f, const vector<double> &x)
{
    const int n = x.size()-1;

    f << x[1];
    for (int i = 2; i < n; i++)
        f << '\t' << x[i];
    f << endl;
}

void write_x(ofstream &f, const vector<double> &x)
{
    f << x[0];
    for(int k = 1; k < NumOfPnts; k++)
        f << '\t' << x[k];
    f << endl;
}

void expLaxF(){
    ofstream fout("Explicit-Lax-Friedrichs.txt");
    fout << NumOfSteps << "\t" <<NumOfPnts << endl;
    write_x(fout, x);

    for(int i=1; i<=NumOfPnts; i++)
        u_prev[i] = u0(x[i-1]);
    
    u_prev[0] = u_prev[NumOfPnts];
    u_prev[NumOfPnts+1] = u_prev[1];

    write_u(fout, u_prev);

    double  b_0 = 0.5 * (1 + c);
    double b_1 = 0.5 * (1 - c);

    for(int i=1; i<=NumOfSteps; i++){
        for(int j=1; j<=NumOfPnts; j++)
            u_cur[j] = b_0 * u_prev[i - 1] + b_1 * u_prev[i + 1];
        
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];

        write_u(fout, u_cur);
        u_prev.swap(u_cur);
    }
    fout.close();
    cout<<"Solution according to Explicit-Lax-Friedrich Scheme Calculated!"<<endl;
}

void impLaxF(){
     ofstream fout("Implicit-Lax-Friedrichs.txt");
    fout << NumOfSteps << "\t" <<NumOfPnts << endl;
    write_x(fout, x);

    for(int i=1; i<=NumOfPnts; i++)
        u_prev[i] = u0(x[i-1]);
    
    u_prev[0] = u_prev[NumOfPnts];
    u_prev[NumOfPnts+1] = u_prev[1];

    write_u(fout, u_prev);

    double  b_0 = 0.5 * (1 + c);
    double b_1 = 0.5 * (1 - c);

    for(int i=1; i<=NumOfSteps; i++){
        for(int j=1; j<=NumOfPnts; j++)
            u_cur[j] = b_0 * u_prev[i - 1] + b_1 * u_prev[i + 1];
        
        u_cur[0] = u_cur[NumOfPnts];
        u_cur[NumOfPnts+1] = u_cur[1];

        write_u(fout, u_cur);
        u_prev.swap(u_cur);
    }
    fout.close();
    cout<<"Solution according to Implicit-Lax-Friedrich Scheme Calculated!"<<endl;
}

int main(){
    cout<<"Please Enter a CFL Number: "<<endl;
    cin>>CFL;
    cout<<"Initialize Grid"<<endl;
    x[0] = xL;
    for(int i = 1; i<NumOfPnts; i++)
        x[i] = x[i-1] + dx;
    expLaxF();
    return 0;
}