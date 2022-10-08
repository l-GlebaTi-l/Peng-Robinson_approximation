#include <iostream>     /// for console "cout <<" and "cin >>"
#include <iomanip>      /// only for "cout << fixed << setprecision(5);" in void build()
#include <math.h>       /// for isnan();
#include <fstream>      /// for r-w files

#define IS_DEBUGG 0

#define MAX_COMP_NAME 10
#define MAX_COMPONENTS 20
#define PARAM_GRAPH_SCALE 5    //  меняет количество значений x и y в графике - прямая (положительная) зависимость
#define PARAM_GRAPH_LINES_INVISIBLE 200
#define PARAM_GRAPH_LINES_DASHES 100
#define PARAM_ADD_PIXELS_TO_GRAPH_LEFT 120
#define PARAM_ADD_PIXELS_TO_GRAPH_DOWN 60
#define PARAM_GRAPH_VALUE_SCALE 15
int PARAM_MOLAR_OR_SPECIFIC_VOLUME = 0;

using namespace std;

#include "data/graphbuilder.h"

#define R 0.08314   /// бар*л / моль*К
#define M_PI (3.141592653589793)
#define M_2PI (2.*M_PI)

#define OmegaA0 0.457235529
#define OmegaB0 0.077796074

#define FOR(I,N) for(int I = 0; I < N; I++)
#define FOR1(I,N) for(int I = 1; I < N+1; I++)

struct fluid_components;
class Fluid;

int n_in_char(const char* str);
int Solution_PR_auto(long double *w, long double T, Fluid fluid, int jp);
int Solution_PR_auto(long double *w, long double T, long double p, fluid_components fluid);
int Solution_Cubic(long double *x, long double a, long double b, long double c);
void Build_PR_Liquid_Approximation(unsigned char *graph, long double a, long double b, long double pz, long double T, long double P0, long double P1, long double W0, long double W1, unsigned char color_R, unsigned char color_G, unsigned char color_B, int Ix, int Iy);
void Build_PR_GAS_Approximation(unsigned char *graph, long double beta, long double bg, long double T, long double P0, long double P1, long double W0, long double W1, unsigned char color_R, unsigned char color_G, unsigned char color_B, int Ix, int Iy);

struct fluid_components{
    char* comp = new char[MAX_COMP_NAME + 1];
    long double C_in_deposit = 0;    // concentration
    long double Tc = 0;      // T crit
    long double Pc = 0;      // P crit
    long double Wc = 0;      // V crit
    long double acf = 0;     // acentric factor
    long double MW = 0;      // molar weight
};

struct component_concentrations{
    int n_p = 0;
    char **name;
    long double *p;
    long double **C;    //  C[i][j], первый массив i выделен под количество компонент, второй массив j - значение концентрации для компоненты i при давлении p[j]
};

struct root_G_L{    //  в большей степени введён для оптимизации расчётов (полноценно ещё не внедрён)
    long double p = -1;  //  давление, при котором получены корни
    long double G = -1;  //  газовый корень
    long double L = -1;  //  жидкий корень
};

class Fluid{
protected:
    void compute_gamma(){           //  вызов происходит в Load_GAS_and_Liquid_Concentrations() после считывания концентраций
        gamma_G = new long double[CT_G.n_p];
        for(int j = 0; j < CT_G.n_p; j++){
            gamma_G[j] = 0;
            for(int i = 0; i < N_comp; i++){
                gamma_G[j] += CT_G.C[i][j] * (1000 / component[i].MW);
            }
            if(IS_DEBUGG == 1) cout << j << " gamma_G = " << gamma_G[j] << endl;
        }
        if(IS_DEBUGG == 1) cout << endl;

        gamma_L = new long double[CT_L.n_p];
        for(int j = 0; j < CT_L.n_p; j++){
            gamma_L[j] = 0;
            for(int i = 0; i < N_comp; i++){
                gamma_L[j] += CT_L.C[i][j] * (1000 / component[i].MW);
            }
            if(IS_DEBUGG == 1) cout << j << " gamma_L = " << gamma_L[j] << endl;
        }
        if(IS_DEBUGG == 1) cout << endl;
    }
    void convert_concentration_to_molar(){
        if(N_comp > 1){
            long double *MC = new long double[N_comp];  /// MC - Molar Concentration
            long double f_MW = 0;

            for(int i = 0; i < N_comp; i++){
//                MC[i] = (component[i].C * 100.) / component[i].MW;
                f_MW += MC[i];
            }
            for(int i = 0; i < N_comp; i++){
//                component[i].C = MC[i] / f_MW;
            }

            delete(MC);
        }
    }
    long double calculate_Wcrit_V_D_V_method(){    /// рассчёт критического молярного объёма через формулу-вывод из уравнения состояния Ван-Дер-Ваальса
        if((isnan(Tc) == 1 || isnan(Pc) == 1) && N_comp > 1){cout << "Input fluid Tc and Pc: "; cin >> Tc >> Pc;}
        if(N_comp == 1){
            Tc = component[0].Tc;
            Pc = component[0].Pc;
            Wc = 3.*R*Tc / (8.*Pc);
            return Wc;
        }
        Wc = 3.*R*Tc / (8.*Pc);
        return Wc;
    }
    long double calculate_Wcrit(){
        if((isnan(Tc) == 1 || isnan(Pc) == 1) && N_comp > 1){cout << "Input fluid Tc and Pc: "; cin >> Tc >> Pc;}
        long double W[3];
        int n_roots;

        if(N_comp == 1) n_roots = Solution_PR_auto(W, component[0].Tc, component[0].Pc, component[0]);
        else            return calculate_Wcrit_V_D_V_method();

        if(n_roots == 3) Wc = min3(W);
        else             Wc = W[0];
        return Wc;
    }
    int Load_Concentrations(char* filename, component_concentrations &CT){
        CT.C = new long double*[N_comp];
        CT.name = new char*[N_comp];
        if(N_comp > 1){
            ifstream file (filename);
            if(!file.is_open()){
                cout << "NO FLUID CONCENTRATIONS DATA FOR THIS TEMPERATURE! Save into \"" << filename << "\"" << endl;
                return -1;
            }
            long double temp_ld;
            char temp_ch[100];

            file >> temp_ch;
            for(int i = 0; i < N_comp; i++){
                CT.name[i] = new char[MAX_COMP_NAME];
                file >> CT.name[i];
            }
            while(!file.eof()){     //  расчёт, сколько значений p в файле данных
                for(int i = 0; i < 1 + N_comp; i++){
                    file >> temp_ld;
                }
                CT.n_p++;
            }
            {   //  выделение памяти под массивы концентраций и соответствующих давлений
            CT.p = new long double[CT.n_p];
            for(int i = 0; i < N_comp; i++){
                CT.C[i] = new long double[CT.n_p];
            }
            }

            file.close();
            file.open(filename);
            for(int i = 0; i < N_comp+1; i++){
                file >> temp_ch;
            }
            for(int pj = 0; pj < CT.n_p; pj++){
                file >> CT.p[pj];
                for(int i = 0; i < N_comp; i++){
                    file >> CT.C[i][pj];
                }
            }
            file.close();
            return 2;
        }
        else if(N_comp == 1){
            CT.C[0] = new long double[1];
            CT.C[0][0] = 1.;
            return 1;
        }
        else return -1;
    }
public:
    fluid_components *component;
    component_concentrations CT_G;            //  двумерный массив концентраций с массивом соответствующих давлений
    component_concentrations CT_L;            //  двумерный массив концентраций с массивом соответствующих давлений
    int N_comp = 0;
    long double Tc;
    long double Pc;
    long double Wc;             //  пока что выступает только в качестве определения фазы (жидкая или газ). Является плохо приближенным методом, практически неприменим в закритическом состоянии.
    long double *gamma_G;         //  коэффициент перевода молярного объёма в удельный. Считается в момент считывания концентраций.
    long double *gamma_L;
    long double p_min = 0;
    long double p_max = 250;
    Fluid(char* filename){
        component = new fluid_components[MAX_COMPONENTS];

        ifstream file (filename);

        if (!file.is_open()){
            cout << "FILE NOT OPEN!\n" << endl;
            return;
        }

        while(!file.eof()){
            file >> component[N_comp].comp >> component[N_comp].C_in_deposit >> component[N_comp].Tc >> component[N_comp].Pc >> component[N_comp].acf >> component[N_comp].MW;
            N_comp++;
        }
        file.close();

//        int ans;
//        cout << "Convert concentration to molar?" << endl;
//        if(IS_DEBUGG == 0)  cin >> ans;
//        else                ans = 1;
//
//        if(ans == 1) convert_concentration_to_molar();
//        cout << endl;

        if(filename == "substance/fluid.txt") {Tc = 667.65; Pc = 41.6;}   //  значения критической точки смеси из моей работы

        if(N_comp > 1)       calculate_Wcrit_V_D_V_method();       /// нахождение критического молярного объёма выведенной формулой из уравнения Ван-дер-Ваальса
        else if(N_comp == 1) calculate_Wcrit();       /// нахождение критического молярного объёма прямым рассчётом (через УРС ПР)
        else cout << "NO FLUID DATA!" << endl;
    }
    Fluid(char* filename, int N_Components){
        component = new fluid_components[N_Components];

        ifstream file (filename);

        if (!file.is_open()){
            cout << "FILE NOT OPEN!\n" << endl;
            return;
        }

        while(!file.eof()){
            file >> component[N_comp].comp >> component[N_comp].C_in_deposit >> component[N_comp].Tc >> component[N_comp].Pc >> component[N_comp].acf >> component[N_comp].MW;
            N_comp++;
        }
        file.close();

//        int ans;
//        cout << "Convert concentration to molar?" << endl;
//        if(IS_DEBUGG == 0)  cin >> ans;
//        else                ans = 1;
//
//        if(ans == 1) convert_concentration_to_molar();
//        cout << endl;

        if(filename == "substance/fluid.txt") {Tc = 667.65; Pc = 41.6;}   //  значения критической точки смеси из моей работы

        calculate_Wcrit_V_D_V_method();       /// нахождение критического молярного объёма выведенной формулой из уравнения Ван-дер-Ваальса
//        calculate_Wcrit();       /// нахождение критического молярного объёма прямым рассчётом (через УРС ПР)

    }
    int Load_GAS_and_Liquid_Concentrations(int temperature_in_C){
        char filename_G[100] = "substance/", filename_L[100] = "substance/";

        {
            char clear_filename[100];
            sprintf(clear_filename, "%d", temperature_in_C);
            int n = n_in_char(clear_filename);
            char tchar[10] = "_gas.TXT";
            for(int i = n; i < n+9; i++)
                clear_filename[i] = tchar[i-n];
            n = n_in_char(clear_filename);

            for(int i = 0; i < n+1; i++)
                filename_G[10 + i] = clear_filename[i];
        }

        if(Load_Concentrations(filename_G, CT_G) == -1) return -1;

        {
            char clear_filename[100];
            sprintf(clear_filename, "%d", temperature_in_C);
            int n = n_in_char(clear_filename);
            char tchar[10] = "_oil.TXT";
            for(int i = n; i < n+9; i++)
                clear_filename[i] = tchar[i-n];
            n = n_in_char(clear_filename);

            for(int i = 0; i < n+1; i++)
                filename_L[10 + i] = clear_filename[i];
        }

        if(Load_Concentrations(filename_L, CT_L) == -1) return -1;

        compute_gamma();
        return 0;
    }
    void print_Concentrations(){
        cout << "GAS CONCENTRATIONS:" << endl;
        cout << CT_G.n_p << endl;
        for(int j = 0; j < N_comp; j++){
            cout << CT_G.name[j] << " =";
            for(int i = 0; i < CT_G.n_p; i++){
                cout << " " << CT_G.C[j][i];
            }
            cout << endl;
        }
        cout << endl;
        for(int i = 0; i < CT_G.n_p; i++)
            cout << CT_G.p[i] << " ";
        cout << endl;

        cout << "LIQUID CONCENTRATIONS:" << endl;
        cout << CT_L.n_p << endl;
        for(int j = 0; j < N_comp; j++){
            cout << CT_L.name[j] << " =";
            for(int i = 0; i < CT_L.n_p; i++){
                cout << " " << CT_L.C[j][i];
            }
            cout << endl;
        }
        cout << endl;
        for(int i = 0; i < CT_L.n_p; i++)
            cout << CT_L.p[i] << " ";
        cout << endl;
    }
    void print_pressure_j_from_data(){
        cout << "Available pressures:" << endl;
        for(int i = 0; i < CT_G.n_p; i++)
            cout << i << ") " << CT_G.p[i] << endl;
    }
    void add_comp(fluid_components new_component){
        component[N_comp].comp = new_component.comp;
        component[N_comp].C_in_deposit = new_component.C_in_deposit;
        component[N_comp].Tc = new_component.Tc;
        component[N_comp].Pc = new_component.Pc;
        component[N_comp].Wc = new_component.Wc;
        component[N_comp].acf = new_component.acf;
        component[N_comp].MW = new_component.MW;
        N_comp++;
    }
    void set_crit_param(long double Tcrit, long double Pcrit){
        Tc = Tcrit;
        Pc = Pcrit;
    }
    void set_crit_param(long double Tcrit, long double Pcrit, long double Wcrit){
        Tc = Tcrit;
        Pc = Pcrit;
        Wc = Wcrit;
    }
};

int Solution_Cubic(long double *x, long double a, long double b, long double c)
{
    long double q, r, r2, q3;

    q=(a*a-3.*b)/9.;
    r=(a*(2.*a*a-9.*b)+27.*c)/54.;
    r2=r*r; q3=q*q*q;

    if(r2<q3) {
        long double t=acos(r/sqrt(q3));
        a/=3.; q=-2.*sqrt(q);
        x[0]=q*cos(t/3.)-a;
        x[1]=q*cos((t+M_2PI)/3.)-a;
        x[2]=q*cos((t-M_2PI)/3.)-a;
        return(3);
    }
    else {
        long double aa,bb;
        aa = -1 * sign(r) * pow(module(r)+sqrt(r2-q3),1./3.);
        if(aa!=0.) bb=q/aa;
        else bb=0.;
        a/=3.; q=aa+bb; r=aa-bb;
        x[0]=q-a;
        x[1]=(-0.5)*q-a;
        x[2]=(sqrt(3.)*0.5)*fabs(r);
        if(x[2] > -0.0000000000001 && x[2] < 0.0000000000001) return(2);
        return(1);
    }
}

int Solution_PR_auto(long double *w, long double T, Fluid fluid, int jp)   /// returns num of P-R roots,  jp is currect pressure number in component_concentrations.p
{
    int nWG, nWL;
    long double W_gas, W_liquid;
    {           /// FOR GAS
        long double At = 0, Bt = 0, A[fluid.N_comp], B[fluid.N_comp];

        FOR(i,fluid.N_comp){
            long double Pr, Tr;
            Pr = fluid.CT_G.p[jp]/fluid.component[i].Pc;
            Tr = T/fluid.component[i].Tc;

            if(fluid.component[i].acf <= 0.49)    A[i] = OmegaA0 * stp(1 + (0.37464 + 1.54226*fluid.component[i].acf - 0.26992*stp(fluid.component[i].acf))*(1 - sqrt(Tr))) * Pr / (Tr * Tr);
            else                        A[i] = OmegaA0 * stp(1 + (0.37964 + 1.48503*fluid.component[i].acf - 0.16442*stp(fluid.component[i].acf) + 0.01666*stp(fluid.component[i].acf, 3))*(1 - sqrt(Tr))) * Pr / (Tr * Tr);
            B[i] = OmegaB0 * Pr/Tr;
        }

        long double sigma = 0;
        FOR(i,fluid.N_comp){
            FOR(j,fluid.N_comp){
                if((i < 2)&&(j > 1) || ((i > 1)&&(j < 2))) sigma = 0.1;
                else sigma = 0;
                At += fluid.CT_G.C[i][jp]*fluid.CT_G.C[j][jp]*(1. - sigma)*sqrt(A[i]*A[j]);  // sigma - look for diplom work Kodzaev Nikita, p.19
            }
        }
        FOR(i,fluid.N_comp){
            Bt += fluid.CT_G.C[i][jp]*B[i];
        }

        long double C =  R*T / fluid.CT_G.p[jp];    //  it's need to get coeff in front of w^3 be = 1
    //    long double cZ = p / R*T;

    //    z^3*C*C*C - (1 - B)*C*C*C*z^2 + (A - 3*B*B - 2*B)*C*C*C*z - (A*B - B*B - B*B*B)*C*C*C = 0;
    //    w^3       - (1 - B)*C * w^2   + (A - 3*B*B - 2*B)*C*C * w - (A*B - B*B - B*B*B)*C*C*C = 0;

        nWG = Solution_Cubic(w, ((Bt - 1)*C), ((At - 3*Bt*Bt - 2*Bt)*C*C), ((Bt*Bt*Bt + Bt*Bt - At*Bt)*C*C*C));
        if(nWG == 1)        W_gas = w[0];
        else if(nWG == 3)   W_gas = max3(w);
        else                cout << "ERROR IN P-R !!!!!!!!!!!!!!!!!!! 2 roots in gas" << endl;
    }
    {           /// FOR LIQUID
        long double At = 0, Bt = 0, A[fluid.N_comp], B[fluid.N_comp];

        FOR(i,fluid.N_comp){
            long double Pr, Tr;
            Pr = fluid.CT_L.p[jp]/fluid.component[i].Pc;
            Tr = T/fluid.component[i].Tc;

            if(fluid.component[i].acf <= 0.49)    A[i] = OmegaA0 * stp(1 + (0.37464 + 1.54226*fluid.component[i].acf - 0.26992*stp(fluid.component[i].acf))*(1 - sqrt(Tr))) * Pr / (Tr * Tr);
            else                        A[i] = OmegaA0 * stp(1 + (0.37964 + 1.48503*fluid.component[i].acf - 0.16442*stp(fluid.component[i].acf) + 0.01666*stp(fluid.component[i].acf, 3))*(1 - sqrt(Tr))) * Pr / (Tr * Tr);
            B[i] = OmegaB0 * Pr/Tr;
        }

        long double sigma = 0;
        FOR(i,fluid.N_comp){
            FOR(j,fluid.N_comp){
                if((i < 2)&&(j > 1) || ((i > 1)&&(j < 2))) sigma = 0.1;
                else sigma = 0;
                At += fluid.CT_L.C[i][jp]*fluid.CT_L.C[j][jp]*(1. - sigma)*sqrt(A[i]*A[j]);  // sigma - look for diplom work Kodzaev Nikita, p.19
            }
        }
        FOR(i,fluid.N_comp){
            Bt += fluid.CT_L.C[i][jp]*B[i];
        }

        long double C =  R*T / fluid.CT_L.p[jp];    //  it's need to get coeff in front of w^3 be = 1
    //    long double cZ = p / R*T;

    //    z^3*C*C*C - (1 - B)*C*C*C*z^2 + (A - 3*B*B - 2*B)*C*C*C*z - (A*B - B*B - B*B*B)*C*C*C = 0;
    //    w^3       - (1 - B)*C * w^2   + (A - 3*B*B - 2*B)*C*C * w - (A*B - B*B - B*B*B)*C*C*C = 0;

        nWL = Solution_Cubic(w, ((Bt - 1)*C), ((At - 3*Bt*Bt - 2*Bt)*C*C), ((Bt*Bt*Bt + Bt*Bt - At*Bt)*C*C*C));
        if(nWL == 1)        W_liquid = w[0];
        else if(nWL == 3)   W_liquid = min3(w);
        else                cout << "ERROR IN P-R !!!!!!!!!!!!!!!!!!! 2 roots in liquid" << endl;
    }

//    if(nWG != nWL) cout << "WARNING ---------> NOT TWO-PHASE FLUID AREA!!!!!!!!!!!!!!!!!" << endl;
    if(nWG != 1 && nWG != 2 && nWG != 3 && nWL != 1 && nWL != 2 && nWL != 3) cout << "WARNING ---------> UNKNOWN ERROR IN P-R!!!!!!!!!!!!!!!!!" << endl;

    if(PARAM_MOLAR_OR_SPECIFIC_VOLUME == 1){
        w[0] = W_gas * fluid.gamma_G[jp];
        w[1] = W_liquid * fluid.gamma_L[jp];
        w[2] = (w[0]+w[1])/2;
    }
    else{
        w[0] = W_gas;
        w[1] = W_liquid;
        w[2] = (w[0]+w[1])/2;
    }
//    cout << "P-R by p = " << fluid.CT_G.p[jp] << ", G roots = " << nWG << ", L roots = " << nWL << ", W_gas = " << W_gas << ", W_liquid = " << W_liquid << endl;

    return 3;
}

int Solution_PR_auto(long double *w, long double T, long double p, fluid_components fluid)   /// returns num of P-R roots
{
    long double A, B;

    long double Pr, Tr;
    Pr = p/fluid.Pc;
    Tr = T/fluid.Tc;

    if(fluid.acf <= 0.49)   A = OmegaA0 * stp(1 + (0.37464 + 1.54226*fluid.acf - 0.26992*stp(fluid.acf))*(1 - sqrt(Tr))) * Pr / (Tr * Tr);
    else                    A = OmegaA0 * stp(1 + (0.37964 + 1.48503*fluid.acf - 0.16442*stp(fluid.acf) + 0.01666*stp(fluid.acf, 3))*(1 - sqrt(Tr))) * Pr / (Tr * Tr);
    B = OmegaB0 * Pr/Tr;

    long double C =  R*T / p;    //  it's need to get coeff in front of w^3 be = 1
//    long double cZ = p / R*T;

//    z^3*C*C*C - (1 - B)*C*C*C*z^2 + (A - 3*B*B - 2*B)*C*C*C*z - (A*B - B*B - B*B*B)*C*C*C = 0;
//    w^3       - (1 - B)*C * w^2   + (A - 3*B*B - 2*B)*C*C * w - (A*B - B*B - B*B*B)*C*C*C = 0;

    int nW = 0; //  num of roots (w)
    nW = Solution_Cubic(w, ((B - 1)*C), ((A - 3*B*B - 2*B)*C*C), ((B*B*B + B*B - A*B)*C*C*C));

    if(nW != 1 && nW != 2 && nW != 3) cout << "WARNING ---------> UNKNOWN ERROR IN P-R!!!!!!!!!!!!!!!!!" << endl;

    return nW;
}

void Build_PR_Graph_PW_with_approximations(Fluid fluid, long double gT, long double a, long double b, long double pz, long double beta, long double bg, long double max_error_G, long double max_error_L, int Ix, int Iy)
{
    int Ix_G = Ix, Iy_G = Iy,
        Ix_L = Ix, Iy_L = Iy;
    int ans;
    if(fluid.N_comp > 1)
        ans = 1;
    else if(fluid.N_comp == 1)
        ans = 2;
    else{cout << "Error while building graph - zero components in substance!!!" << endl; return;}

    long double GP0, GP1, GW0, GW1, LP0, LP1, LW0, LW1;
    if(ans == 1){
        if(0) switch((int)(gT - 273.15)){
        case 0:
            if(IS_DEBUGG == 1) cout << "T = 0" << endl;
            GP0 = 0.5;
            GP1 = 5.;
            GW0 = 0.;
            GW1 = 25.;
            LP0 = 0.5;
            LP1 = 5.;
            LW0 = 0.174;
            LW1 = 0.19;
            break;
        case 100:
            if(IS_DEBUGG == 1) cout << "T = 100" << endl;
            GP0 = 0.;
            GP1 = 15.;
            GW0 = 0.;
            GW1 = 40.;
            LP0 = 0.;
            LP1 = 15.;
            LW0 = 0.18;
            LW1 = 0.22;
            break;
        case 200:
            if(IS_DEBUGG == 1) cout << "T = 200" << endl;
            GP0 = 0.;
            GP1 = 30.;
            GW0 = 0.;
            GW1 = 40.;
            LP0 = 0.;
            LP1 = 30.;
            LW0 = 0.2;
            LW1 = 0.3;
            break;
        case 300:
            if(IS_DEBUGG == 1) cout << "T = 300" << endl;
            GP0 = 0.;
            GP1 = 45.;
            GW0 = 0.;
            GW1 = 15.;
            LP0 = 0.;
            LP1 = 45.;
            LW0 = 0.26;
            LW1 = 0.36;
            break;
        default:
            cout << "DEFAULT TEMPERATURE, NO DATA FOR IT" << endl;
    //        cout << "Enter range of graph as (Pmin, Pmax, Wmin, Wmax) for gas and then for liquid:" << endl;
        }
        else{
            GP0 = (int)fluid.CT_G.p[0];
            if(fluid.CT_G.p[fluid.CT_G.n_p-1] > (int)fluid.CT_G.p[fluid.CT_G.n_p-1])
                GP1 = ((int)fluid.CT_G.p[fluid.CT_G.n_p-1]) + 1;
            else
                GP1 = ((int)fluid.CT_G.p[fluid.CT_G.n_p-1]);
            LP0 = GP0;
            LP1 = GP1;

            long double W0[3], W1[3];
            Solution_PR_auto(W0, gT, fluid, fluid.CT_G.n_p-1);
            Solution_PR_auto(W1, gT, fluid, 0);
            if(W0[0] > W1[0])
                swap(W0[0], W1[0]);
            if(W0[1] > W1[1])
                swap(W0[1], W1[1]);
            GW0 = W0[0];
            LW0 = W0[1];
            GW1 = W1[0];
            LW1 = W1[1];

            GW1 = GW1 + (GW1 - GW0)*0.1;
            if(GW0 - (GW1/3) < 0)
                GW0 = 0;
            else{
                GW0 = GW0 - (W1[0] - GW0)*0.1;
                if(GW0 < 0) GW0 = 0;
            }

            LW1 = LW1 + (LW1 - LW0)*0.1;
            if(LW0 - (LW1/3) < 0)
                LW0 = 0;
            else{
                LW0 = LW0 - (W1[1] - LW0)*0.1;
                if(LW0 < 0) LW0 = 0;
            }
        }
    }
    else{    ///  поиск верхних и нижних границ жидкого и газового графов
//        P0 = fluid.p_min;
//        P1 = fluid.p_max;
        long double dp, Pt, W[3];
        int Nw;
        dp = (fluid.p_max - fluid.p_min) / Ix;

        GW0 = -1; GW1 = -1; LW0 = -1; LW1 = -1;
        GP0 = -1; GP1 = -1; LP0 = -1; LP1 = -1;
        for(int i = 0; i < Ix; i++){
            Pt = fluid.p_min + i*dp;
            Nw = Solution_PR_auto(W, gT, Pt, fluid.component[0]);

            if(Nw == 3){
                if(GP0 == -1)
                    GP0 = Pt;
                GP1 = Pt;
                if(LP0 == -1)
                    LP0 = Pt;
                LP1 = Pt;
                if(max3(W) > GW1)
                    GW1 = max3(W);
                if(max3(W) < GW0)
                    GW0 = max3(W);
                if(min3(W) > LW1)
                    LW1 = min3(W);
                if(min3(W) < LW0)
                    LW0 = min3(W);
            }
            else{
                if(W[0] > fluid.Wc){
                    if(GP0 == -1)
                        GP0 = Pt;
                    GP1 = Pt;
                    if(W[0] > GW1)
                        GW1 = W[0];
                    if(W[0] < GW0)
                        GW0 = W[0];
                }
                else{
                    if(LP0 == -1)
                        LP0 = Pt;
                    LP1 = Pt;
                    if(W[0] > LW1)
                        LW1 = W[0];
                    if(W[0] < LW0)
                        LW0 = W[0];
                }
            }
        }

        if(GP0 < GP1){
            GP0 = int(GP0);
            if(GP1 > ((int)(GP1)))
                GP1 = ((int)GP1) + 1;

            GW1 = GW1 + (GW1 - GW0)*0.1;
            if(GW0 - (GW1/3) < 0)
                GW0 = 0;
            else{
                GW0 = GW0 - (GW1 - GW0)*0.1;
                if(GW0 < 0) GW0 = 0;
            }
        }

        if(LP0 < LP1){
            LP0 = int(LP0);
            if(LP1 > ((int)(LP1)))
                LP1 = ((int)LP1) + 1;

            LW1 = LW1 + (LW1 - LW0)*0.1;
            if(LW0 - (LW1/3) < 0)
                LW0 = 0;
            else{
                LW0 = LW0 - (LW1 - LW0)*0.1;
                if(LW0 < 0) LW0 = 0;
            }
        }
    }
    resolution_shirina = Ix;
    resolution_vysota = Iy;

    unsigned char *graph_G = new unsigned char[Ix*4*Iy];
    unsigned char *graph_L = new unsigned char[Ix*4*Iy];


    Graph_get_ready(graph_G, Ix, Iy);
    Graph_get_ready(graph_L, Ix, Iy);

    if(ans == 1){
        {   /// FOR LIQUID
            int Nw;     //  num of w-roots in P-R
            long double W[3];   //  Pt - P temp, Wt - volume temp
            long double lineW, lineP, koorWdoub_L;
            int koorW, koorW0, koorP, koorP0;

            lineW = LW1 - LW0;
            lineP = LP1 - LP0;
            Nw = Solution_PR_auto(W, gT, fluid, 0);
            koorP = ((fluid.CT_G.p[0] - LP0) / lineP) * Ix;
            koorWdoub_L = ((min3(W) - LW0) / lineW) * Iy;

            for(int jp = 1; jp < fluid.CT_G.n_p; jp++){
                Nw = Solution_PR_auto(W, gT, fluid, jp);

                koorP0 = koorP;
                koorP = ((fluid.CT_G.p[jp] - LP0) / lineP) * Ix;

                if(Nw == 3){
                    koorW0 = koorWdoub_L;
                    koorWdoub_L = ((min3(W) - LW0) / lineW) * Iy;     /// значение "y" координат изображения
                    koorW = koorWdoub_L;
                    Build_Line(graph_L, Ix, Iy, koorP0, koorW0, koorP, koorW, 0, 0, 255);
                }
                else{
                    cout << "UNKNOWN ERROR IN BUILDING PR-GRAPH" << endl;
                }
            }
        }
        {   /// FOR GAS
            int Nw;     //  num of w-roots in P-R
            long double W[3];   //  Pt - P temp, Wt - volume temp
            long double lineW, lineP, koorWdoub_G;
            int koorW, koorW0, koorP, koorP0;

            lineW = GW1 - GW0;
            lineP = GP1 - GP0;
            Nw = Solution_PR_auto(W, gT, fluid, 0);
            koorP = ((fluid.CT_G.p[0] - GP0) / lineP) * Ix;
            koorWdoub_G = ((max3(W) - GW0) / lineW) * Iy;

            for(int jp = 1; jp < fluid.CT_G.n_p; jp++){
                Nw = Solution_PR_auto(W, gT, fluid, jp);

                koorP0 = koorP;
                koorP = ((fluid.CT_G.p[jp] - GP0) / lineP) * Ix;

                if(Nw == 3){
                    koorW0 = koorWdoub_G;
                    koorWdoub_G = ((max3(W) - GW0) / lineW) * Iy;     /// значение "y" координат изображения
                    koorW = koorWdoub_G;
                    Build_Line(graph_G, Ix, Iy, koorP0, koorW0, koorP, koorW, 0, 255, 0);
                }
                else{
                    cout << "UNKNOWN ERROR IN BUILDING PR-GRAPH" << endl;
                }
            }
        }
    }
    else{
        {   /// FOR LIQUID
            if(LP0 >= LP1){
                /// do nothing. No liquid graph.
            }
            else{
                long double dp, Pt, W[3], L_root;
                int Nw;
                int last_i = -10;
                long double lineW, koorWdoub;
                lineW = LW1 - LW0;
                dp = (LP1 - LP0) / Ix;

                int koorW, koorW0;

                for(int i = 0; i < Ix; i++){
                    Pt = LP0 + i*dp;

                    Nw = Solution_PR_auto(W, gT, Pt, fluid.component[0]);
                    if(Nw == 3)                         L_root = min3(W);
                    else if(Nw == 1 && W[0] < fluid.Wc) L_root = W[0];
                    else continue;

                    if(last_i == i-1){
                        koorW0 = koorWdoub;
                    }
                    else{
                        koorWdoub = ((L_root - LW0) / lineW) * Iy;
                        last_i = i;
                        continue;
                    }
                    koorWdoub = ((L_root - LW0) / lineW) * Iy;
                    koorW = koorWdoub;

                    Build_Line(graph_L, Ix, Iy, i, koorW0, koorW, 0, 0, 255);

                    last_i = i;
                }
            }
        }
        {   /// FOR GAS
            if(GP0 >= GP1){
                /// do nothing. No gas graph.
            }
            else{
                long double dp, Pt, W[3], G_root;
                int Nw;
                int last_i = -10;
                long double lineW, koorWdoub;
                lineW = GW1 - GW0;
                dp = (GP1 - GP0) / Ix;

                int koorW, koorW0;

                for(int i = 0; i < Ix; i++){
                    Pt = GP0 + i*dp;

                    Nw = Solution_PR_auto(W, gT, Pt, fluid.component[0]);
                    if(Nw == 3)                         G_root = max3(W);
                    else if(Nw == 1 && W[0] > fluid.Wc) G_root = W[0];
                    else continue;

                    if(last_i == i-1){
                        koorW0 = koorWdoub;
                    }
                    else{
                        koorWdoub = ((G_root - GW0) / lineW) * Iy;
                        last_i = i;
                        continue;
                    }
                    koorWdoub = ((G_root - GW0) / lineW) * Iy;
                    koorW = koorWdoub;

                    Build_Line(graph_G, Ix, Iy, i, koorW0, koorW, 0, 255, 0);

                    last_i = i;
                }
            }
        }
    }

    Build_PR_GAS_Approximation(graph_G, beta, bg, gT, GP0, GP1, GW0, GW1, 0, 120, 0, Ix, Iy);
    Build_PR_Liquid_Approximation(graph_L, a, b, pz, gT, LP0, LP1, LW0, LW1, 0, 0, 120, Ix, Iy);

    Graph_numeric(graph_G, Ix_G, Iy_G, GP0, GP1, GW0, GW1);
    Graph_numeric(graph_L, Ix_L, Iy_L, LP0, LP1, LW0, LW1);

//    eng.say(0, 0, 32, "error=", graph_G, Ix, Iy);
    Graph_print_error(graph_G, max_error_G, Ix_G, Iy_G);
    Graph_print_error(graph_L, max_error_L, Ix_L, Iy_L);

    if(ans == 1){
        char filename[256];
        sprintf(filename, "results/Fluid %.2f°C GAS.bmp", (float)(gT - 273.15));
        saveBMP(filename, graph_G, Ix_G, Iy_G);

        sprintf(filename, "results/Fluid %.2f°C Liquid.bmp", (float)(gT - 273.15));
        saveBMP(filename, graph_L, Ix_L, Iy_L);
    }
    else{
        char filename[256];
        sprintf(filename, "results/%s %.2f°K %d_%d GAS.bmp", fluid.component[0].comp, (float)gT, (int)fluid.p_min, (int)fluid.p_max);
        saveBMP(filename, graph_G, Ix_G, Iy_G);

        sprintf(filename, "results/%s %.2f°K %d_%d Liquid.bmp", fluid.component[0].comp, (float)gT, (int)fluid.p_min, (int)fluid.p_max);
        saveBMP(filename, graph_L, Ix_L, Iy_L);
    }
//    saveBMP("Graph P-R and approximation GAS.bmp", graph_G, Ix_G, Iy_G);
//    saveBMP("Graph P-R and approximation LIQUID.bmp", graph_L, Ix_L, Iy_L);

    delete(graph_G);
    delete(graph_L);
}

void Build_PR_Liquid_Approximation(unsigned char *graph, long double a, long double b, long double pz, long double T, long double P0, long double P1, long double W0, long double W1, unsigned char color_R, unsigned char color_G, unsigned char color_B, int Ix, int Iy)
{
    resolution_shirina = Ix;
    resolution_vysota = Iy;
    int dx, dy;
//    cout << w1 << " ?= " << ((a*R*T)/(p[0] + pz) + b) << endl;

    long double dp, lineW;

    dp = (P1 - P0) / Ix;
    lineW = W1 - W0;

    long double W, p;   //  wt - w-temp;
    int koorW, koorWt;

    W = ((a*R*T)/(P0 + pz) + b);
    koorW = ((W - W0) / lineW) * Iy;

    int LINE_type_param = sqrt(Ix*Ix + Iy*Iy) / PARAM_GRAPH_LINES_DASHES;
    int LINE_type_counter = 0;
    for(int i = 1; i < Ix; i++){
        p = P0 + i*dp;

        W = ((a*R*T)/(p + pz) + b);

        koorWt = koorW;
        koorW = ((W - W0) / lineW) * Iy;
        Build_Line(graph, Ix, Iy, i, koorWt, koorW, LINE_type_param, LINE_type_counter, color_R, color_G, color_B);
    }
}

void Build_PR_GAS_Approximation(unsigned char *graph, long double beta, long double bg, long double T, long double P0, long double P1, long double W0, long double W1, unsigned char color_R, unsigned char color_G, unsigned char color_B, int Ix, int Iy)
{
    resolution_shirina = Ix;
    resolution_vysota = Iy;
    int dx, dy;
//    cout << w1 << " ?= " << ((a*R*T)/(p[0] + pz) + b) << endl;

    long double dp, lineW;

    dp = (P1 - P0) / Ix;
    lineW = W1 - W0;

    long double W, p;   //  wt - w-temp;
    int koorW, koorWt;

    W = ((beta*R*T)/(P0) + bg);
    koorW = ((W - W0) / lineW) * Iy;

    int LINE_type_param = sqrt(Ix*Ix + Iy*Iy) / PARAM_GRAPH_LINES_DASHES;
    int LINE_type_counter = 0;
    for(int i = 1; i < Ix; i++){
        p = P0 + i*dp;

        W = ((beta*R*T)/(p) + bg);

        koorWt = koorW;
        koorW = ((W - W0) / lineW) * Iy;
        Build_Line(graph, Ix, Iy, i, koorWt, koorW, LINE_type_param, LINE_type_counter, color_R, color_G, color_B);
    }
}

void Build_PR_Graph_PW(Fluid fluid, long double gT, long double P0, long double P1, long double W0, long double W1, int Ix, int Iy)     /// 1 if you want to build PR result
{
    resolution_shirina = Ix;
    resolution_vysota = Iy;
    int dx, dy;

    unsigned char *graph = new unsigned char[Ix*4*Iy];

    int Nw;     //  num of w-roots in P-R
    long double W[3], dw, dp, Wt, Pt;   //  Pt - P temp, Wt - volume temp

    for(int j = 0; j < Iy; j++){    /// clearning the image
        for(int i = 0; i < Ix; i++){
            graph[j*Ix*4 + i*4    ] = 255;    /// B
            graph[j*Ix*4 + i*4 + 1] = 255;    /// G
            graph[j*Ix*4 + i*4 + 2] = 255;    /// R
            graph[j*Ix*4 + i*4 + 3] = 0;
        }
    }

    for(int i = 0; i < Ix; i++){
        graph[i*4    ] = 0;    /// B
        graph[i*4 + 1] = 0;    /// G
        graph[i*4 + 2] = 0;    /// R
        graph[i*4 + 3] = 0;
        graph[(Iy - 1)*Ix*4 + i*4    ] = 0;    /// B
        graph[(Iy - 1)*Ix*4 + i*4 + 1] = 0;    /// G
        graph[(Iy - 1)*Ix*4 + i*4 + 2] = 0;    /// R
        graph[(Iy - 1)*Ix*4 + i*4 + 3] = 0;
    }
    for(int j = 0; j < Iy; j++){    ///
        graph[j*Ix*4    ] = 0;    /// B
        graph[j*Ix*4 + 1] = 0;    /// G
        graph[j*Ix*4 + 2] = 0;    /// R
        graph[j*Ix*4 + 3] = 0;
        graph[j*Ix*4 + (Ix - 1)*4    ] = 0;    /// B
        graph[j*Ix*4 + (Ix - 1)*4 + 1] = 0;    /// G
        graph[j*Ix*4 + (Ix - 1)*4 + 2] = 0;    /// R
        graph[j*Ix*4 + (Ix - 1)*4 + 3] = 0;
    }

    {                /// расчерчивание графика полосами
        if(Ix >= Iy){
            dx = (int)((double)Ix / ((double)PARAM_GRAPH_SCALE * (double)Ix / (double)Iy));
            dy = (int)((double)Iy / ((double)PARAM_GRAPH_SCALE*1.5));
        }
        else{
            dx = (int)(Ix / (int)((double)PARAM_GRAPH_SCALE));
            dy = (int)((double)Iy*1.5 / ((double)PARAM_GRAPH_SCALE * (double)Iy / (double)Ix));
        }

        for(int i = 0; i < Ix; i += dx){
            for(int j = 0; j < Iy; j++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
        for(int j = 0; j < Iy; j += dy){
            for(int i = 0; i < Ix; i++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
    }

    long double lineW, koorWdoub[3];
    dp = (P1 - P0) / Ix;

    if(0){
        cout << "Find automatic Wmin and Wmax?" << endl;
        for(int i = 0; i < Ix; i++){
            Pt = P0 + i*dp;

            Nw = Solution_PR_auto(W, gT, Pt, fluid.component[0]);
            if(Nw == 3){
                if(max3(W) > W1)
                    W1 = max3(W);
                if(min3(W) < W0)
                    W0 = min3(W);
            }
            else{
                if(W[0] > W1)
                    W1 = max3(W);
                if(W[0] < W0)
                    W0 = min3(W);
            }
        }
        cout << "Wmin and Wmax = " << W0 << " " << W1 << endl;
    }

    lineW = W1 - W0;

    int koorW, koorW0, kFrom, kTo;



    Nw = Solution_PR_auto(W, gT, P0, fluid.component[0]);
    koorWdoub[0] = ((W[0] - W0) / lineW) * Iy;
    if(Nw == 3){
        koorWdoub[1] = ((max3(W) - W0) / lineW) * Iy;
        koorWdoub[2] = ((min3(W) - W0) / lineW) * Iy;
    }

    for(int i = 1; i < Ix; i++){
        Pt = P0 + i*dp;

        Nw = Solution_PR_auto(W, gT, Pt, fluid.component[0]);

        if(Nw == 1){
            koorW0 = koorWdoub[0];
            koorWdoub[0] = ((W[0] - W0) / lineW) * Iy;
            koorW = koorWdoub[0];
            if(fluid.N_comp == 1){
                if(W[0] > fluid.Wc)  Build_Line(graph, Ix, Iy, i, koorW0, koorW, 0, 255, 0);
                else                    Build_Line(graph, Ix, Iy, i, koorW0, koorW, 0, 0, 255);
            }
            else{
                if(W[0] > fluid.Wc)  Build_Line(graph, Ix, Iy, i, koorW0, koorW, 0, 255, 0);
                else                    Build_Line(graph, Ix, Iy, i, koorW0, koorW, 0, 0, 255);
            }
        }
        else if(Nw == 3){
            koorW0 = koorWdoub[1];
            koorWdoub[1] = ((max3(W) - W0) / lineW) * Iy;     /// значение "y" координат изображения
            koorW = koorWdoub[1];
            Build_Line(graph, Ix, Iy, i, koorW0, koorW, 0, 255, 0);

            koorW0 = koorWdoub[2];
            koorWdoub[2] = ((min3(W) - W0) / lineW) * Iy;     /// значение "y" координат изображения
            koorW = koorWdoub[2];
            Build_Line(graph, Ix, Iy, i, koorW0, koorW, 0, 0, 255);
        }
        else if(Nw == 2){
            cout << "ERROR IN BUILDING PR-GRAPH, 2 ROOTS" << endl;
        }
        else{
            cout << "UNKNOWN ERROR IN BUILDING PR-GRAPH" << endl;
        }
    }

    {                /// Увеличиваем изображение для накидывания размерностей. К примеру, по 120 пикселей с каждой стороны (чтобы нормально влезли числа, пока без увеличения на параметры изображения)
        int addx = PARAM_ADD_PIXELS_TO_GRAPH_LEFT, addy = PARAM_ADD_PIXELS_TO_GRAPH_DOWN;     /// на сколько увеличить изображение слева и снизу
        int dIx, dIy;
        dIx = Ix + addx;
        dIy = Iy + addy;
        unsigned char *diagram = new unsigned char[dIx*4*dIy];

        for(int j = 0; j < dIy; j++){    /// clearning the image
            for(int i = 0; i < dIx; i++){
                diagram[j*dIx*4 + i*4    ] = 255;    /// B
                diagram[j*dIx*4 + i*4 + 1] = 255;    /// G
                diagram[j*dIx*4 + i*4 + 2] = 255;    /// R
                diagram[j*dIx*4 + i*4 + 3] = 0;
            }
        }

        for(int j = 0; j < Iy; j++){
            for(int i = 0; i < Ix*4; i++){
                diagram[(j+addy)*dIx*4 + addx*4 + i] = graph[j*Ix*4 + i];
            }
        }

        char str[20];
        int n;
        double PT, WT, dP, dW; //  P temp, W temp
        PT = P0;
        WT = W0;
        dP = (P1 - P0) / ((double)Ix / (double)dx);
        dW = (W1 - W0) / ((double)Iy / (double)dy);

        for(int cX = addx; cX < dIx; cX += dx){
            sprintf(str, "%.f", PT);
            eng.say(cX, addy - (PARAM_GRAPH_VALUE_SCALE*1.25), PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
            PT += dP;
        }
        for(int cY = addy; cY < dIy; cY += dy){
            sprintf(str, "%.3f", WT);
            n = n_in_char(str);
            eng.say(addx - n*PARAM_GRAPH_VALUE_SCALE, cY, PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
            WT += dW;
        }
        eng.say(dIx - PARAM_GRAPH_VALUE_SCALE*2.5*5, 0, PARAM_GRAPH_VALUE_SCALE*2.5, "P,bar", diagram, dIx, dIy);
        eng.say(0, dIy - PARAM_GRAPH_VALUE_SCALE*2.5, PARAM_GRAPH_VALUE_SCALE*2.5, "T,K", diagram, dIx, dIy);

        saveBMP("Graph_Peng_Robinson.bmp", diagram, dIx, dIy);
        delete(diagram);
    }

//    saveBMP("Graph_Peng_Robinson.bmp", graph, Ix, Iy);

    delete(graph);
}

void Build_PR_Liquid_Approximation(long double a, long double b, long double pz, long double T, long double P0, long double P1, long double W0, long double W1, int Ix, int Iy)
{
    resolution_shirina = Ix;
    resolution_vysota = Iy;
    int dx, dy;
//    cout << w1 << " ?= " << ((a*R*T)/(p[0] + pz) + b) << endl;
//    Ix = 2000, Iy = 500;
    unsigned char *graph = new unsigned char[Ix*4*Iy];

    for(int j = 0; j < Iy; j++){    /// clearning the image
        for(int i = 0; i < Ix; i++){
            graph[j*Ix*4 + i*4    ] = 255;    /// B
            graph[j*Ix*4 + i*4 + 1] = 255;    /// G
            graph[j*Ix*4 + i*4 + 2] = 255;    /// R
            graph[j*Ix*4 + i*4 + 3] = 0;
        }
    }

    for(int i = 0; i < Ix; i++){    /// рамка графа
        graph[i*4    ] = 0;    /// B
        graph[i*4 + 1] = 0;    /// G
        graph[i*4 + 2] = 0;    /// R
        graph[i*4 + 3] = 0;
        graph[(Iy - 1)*Ix*4 + i*4    ] = 0;    /// B
        graph[(Iy - 1)*Ix*4 + i*4 + 1] = 0;    /// G
        graph[(Iy - 1)*Ix*4 + i*4 + 2] = 0;    /// R
        graph[(Iy - 1)*Ix*4 + i*4 + 3] = 0;
    }
    for(int j = 0; j < Iy; j++){    ///
        graph[j*Ix*4    ] = 0;    /// B
        graph[j*Ix*4 + 1] = 0;    /// G
        graph[j*Ix*4 + 2] = 0;    /// R
        graph[j*Ix*4 + 3] = 0;
        graph[j*Ix*4 + (Ix - 1)*4    ] = 0;    /// B
        graph[j*Ix*4 + (Ix - 1)*4 + 1] = 0;    /// G
        graph[j*Ix*4 + (Ix - 1)*4 + 2] = 0;    /// R
        graph[j*Ix*4 + (Ix - 1)*4 + 3] = 0;
    }

    {                /// расчерчивание графика полосами
        if(Ix >= Iy){
            dx = (int)((double)Ix / ((double)PARAM_GRAPH_SCALE * (double)Ix / (double)Iy));
            dy = (int)((double)Iy / ((double)PARAM_GRAPH_SCALE*1.5));
        }
        else{
            dx = (int)(Ix / (int)((double)PARAM_GRAPH_SCALE));
            dy = (int)((double)Iy*1.5 / ((double)PARAM_GRAPH_SCALE * (double)Iy / (double)Ix));
        }

        for(int i = 0; i < Ix; i += dx){
            for(int j = 0; j < Iy; j++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
        for(int j = 0; j < Iy; j += dy){
            for(int i = 0; i < Ix; i++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
    }


    long double dp, lineW;

    dp = (P1 - P0) / Ix;
    lineW = W1 - W0;

    long double W, p;   //  wt - w-temp;
    int koorW, koorWt;

    W = ((a*R*T)/(P0 + pz) + b);
    koorW = ((W - W0) / lineW) * Iy;

    for(int i = 1; i < Ix; i++){
        p = P0 + i*dp;

        W = ((a*R*T)/(p + pz) + b);

        koorWt = koorW;
        koorW = ((W - W0) / lineW) * Iy;
        Build_Line(graph, Ix, Iy, i, koorWt, koorW, 0, 0, 255);
    }


    {                /// Увеличиваем изображение для накидывания размерностей. К примеру, по 120 пикселей с каждой стороны (чтобы нормально влезли числа, пока без увеличения на параметры изображения)
        int addx = PARAM_ADD_PIXELS_TO_GRAPH_LEFT, addy = PARAM_ADD_PIXELS_TO_GRAPH_DOWN;     /// на сколько увеличить изображение слева и снизу
        int dIx, dIy;
        dIx = Ix + addx;
        dIy = Iy + addy;
        unsigned char *diagram = new unsigned char[dIx*4*dIy];

        for(int j = 0; j < dIy; j++){    /// clearning the image
            for(int i = 0; i < dIx; i++){
                diagram[j*dIx*4 + i*4    ] = 255;    /// B
                diagram[j*dIx*4 + i*4 + 1] = 255;    /// G
                diagram[j*dIx*4 + i*4 + 2] = 255;    /// R
                diagram[j*dIx*4 + i*4 + 3] = 0;
            }
        }

        for(int j = 0; j < Iy; j++){
            for(int i = 0; i < Ix*4; i++){
                diagram[(j+addy)*dIx*4 + addx*4 + i] = graph[j*Ix*4 + i];
            }
        }

        char str[20];
        int n;
        double PT, WT, dP, dW; //  P temp, W temp
        PT = P0;
        WT = W0;
        dP = (P1 - P0) / ((double)Ix / (double)dx);
        dW = (W1 - W0) / ((double)Iy / (double)dy);

        for(int cX = addx; cX < dIx; cX += dx){
            sprintf(str, "%.f", PT);
            eng.say(cX, addy - (PARAM_GRAPH_VALUE_SCALE*1.25), PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
            PT += dP;
        }
        for(int cY = addy; cY < dIy; cY += dy){
            sprintf(str, "%.3f", WT);
            n = n_in_char(str);
            eng.say(addx - n*PARAM_GRAPH_VALUE_SCALE, cY, PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
            WT += dW;
        }
        eng.say(dIx - PARAM_GRAPH_VALUE_SCALE*2.5*5, 0, PARAM_GRAPH_VALUE_SCALE*2.5, "P,bar", diagram, dIx, dIy);
        eng.say(0, dIy - PARAM_GRAPH_VALUE_SCALE*2.5, PARAM_GRAPH_VALUE_SCALE*2.5, "T,K", diagram, dIx, dIy);

        saveBMP("Liquid_Approximation.bmp", diagram, dIx, dIy);
        delete(diagram);
    }

//    saveBMP("Graph_Liquid_Approximation.bmp", graph, Ix, Iy);

    delete(graph);
}

void Build_PR_GAS_Approximation(long double beta, long double bg, long double T, long double P0, long double P1, long double W0, long double W1, int Ix, int Iy)
{
    resolution_shirina = Ix;
    resolution_vysota = Iy;
    int dx, dy;
//    cout << w1 << " ?= " << ((a*R*T)/(p[0] + pz) + b) << endl;
//    Ix = 2000, Iy = 500;
    unsigned char *graph = new unsigned char[Ix*4*Iy];

    for(int j = 0; j < Iy; j++){    /// clearning the image
        for(int i = 0; i < Ix; i++){
            graph[j*Ix*4 + i*4    ] = 255;    /// B
            graph[j*Ix*4 + i*4 + 1] = 255;    /// G
            graph[j*Ix*4 + i*4 + 2] = 255;    /// R
            graph[j*Ix*4 + i*4 + 3] = 0;
        }
    }

    for(int i = 0; i < Ix; i++){
        graph[i*4    ] = 0;    /// B
        graph[i*4 + 1] = 0;    /// G
        graph[i*4 + 2] = 0;    /// R
        graph[i*4 + 3] = 0;
        graph[(Iy - 1)*Ix*4 + i*4    ] = 0;    /// B
        graph[(Iy - 1)*Ix*4 + i*4 + 1] = 0;    /// G
        graph[(Iy - 1)*Ix*4 + i*4 + 2] = 0;    /// R
        graph[(Iy - 1)*Ix*4 + i*4 + 3] = 0;
    }
    for(int j = 0; j < Iy; j++){    ///
        graph[j*Ix*4    ] = 0;    /// B
        graph[j*Ix*4 + 1] = 0;    /// G
        graph[j*Ix*4 + 2] = 0;    /// R
        graph[j*Ix*4 + 3] = 0;
        graph[j*Ix*4 + (Ix - 1)*4    ] = 0;    /// B
        graph[j*Ix*4 + (Ix - 1)*4 + 1] = 0;    /// G
        graph[j*Ix*4 + (Ix - 1)*4 + 2] = 0;    /// R
        graph[j*Ix*4 + (Ix - 1)*4 + 3] = 0;
    }

    {                /// расчерчивание графика полосами
        if(Ix >= Iy){
            dx = (int)((double)Ix / ((double)PARAM_GRAPH_SCALE * (double)Ix / (double)Iy));
            dy = (int)((double)Iy / ((double)PARAM_GRAPH_SCALE*1.5));
        }
        else{
            dx = (int)(Ix / (int)((double)PARAM_GRAPH_SCALE));
            dy = (int)((double)Iy*1.5 / ((double)PARAM_GRAPH_SCALE * (double)Iy / (double)Ix));
        }

        for(int i = 0; i < Ix; i += dx){
            for(int j = 0; j < Iy; j++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
        for(int j = 0; j < Iy; j += dy){
            for(int i = 0; i < Ix; i++){
                graph[j*Ix*4 + i*4] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 1] = PARAM_GRAPH_LINES_INVISIBLE;
                graph[j*Ix*4 + i*4 + 2] = PARAM_GRAPH_LINES_INVISIBLE;
            }
        }
    }

    long double dp, lineW;

    dp = (P1 - P0) / Ix;
    lineW = W1 - W0;

    long double W, p;   //  wt - w-temp;
    int koorW, koorWt;

    W = ((beta*R*T)/(P0) + bg);
    koorW = ((W - W0) / lineW) * Iy;

    for(int i = 1; i < Ix; i++){
        p = P0 + i*dp;

        W = ((beta*R*T)/(p) + bg);

        koorWt = koorW;
        koorW = ((W - W0) / lineW) * Iy;
        Build_Line(graph, Ix, Iy, i, koorWt, koorW, 0, 255, 0);
    }

    {                /// Увеличиваем изображение для накидывания размерностей. К примеру, по 120 пикселей с каждой стороны (чтобы нормально влезли числа, пока без увеличения на параметры изображения)
        int addx = PARAM_ADD_PIXELS_TO_GRAPH_LEFT, addy = PARAM_ADD_PIXELS_TO_GRAPH_DOWN;     /// на сколько увеличить изображение слева и снизу
        int dIx, dIy;
        dIx = Ix + addx;
        dIy = Iy + addy;
        unsigned char *diagram = new unsigned char[dIx*4*dIy];

        for(int j = 0; j < dIy; j++){    /// clearning the image
            for(int i = 0; i < dIx; i++){
                diagram[j*dIx*4 + i*4    ] = 255;    /// B
                diagram[j*dIx*4 + i*4 + 1] = 255;    /// G
                diagram[j*dIx*4 + i*4 + 2] = 255;    /// R
                diagram[j*dIx*4 + i*4 + 3] = 0;
            }
        }

        for(int j = 0; j < Iy; j++){
            for(int i = 0; i < Ix*4; i++){
                diagram[(j+addy)*dIx*4 + addx*4 + i] = graph[j*Ix*4 + i];
            }
        }

        char str[20];
        int n;
        double PT, WT, dP, dW; //  P temp, W temp
        PT = P0;
        WT = W0;
        dP = (P1 - P0) / ((double)Ix / (double)dx);
        dW = (W1 - W0) / ((double)Iy / (double)dy);

        for(int cX = addx; cX < dIx; cX += dx){
            sprintf(str, "%.f", PT);
            eng.say(cX, addy - (PARAM_GRAPH_VALUE_SCALE*1.25), PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
            PT += dP;
        }
        for(int cY = addy; cY < dIy; cY += dy){
            sprintf(str, "%.3f", WT);
            n = n_in_char(str);
            eng.say(addx - n*PARAM_GRAPH_VALUE_SCALE, cY, PARAM_GRAPH_VALUE_SCALE, str, diagram, dIx, dIy);
            WT += dW;
        }
        eng.say(dIx - PARAM_GRAPH_VALUE_SCALE*2.5*5, 0, PARAM_GRAPH_VALUE_SCALE*2.5, "P,bar", diagram, dIx, dIy);
        eng.say(0, dIy - PARAM_GRAPH_VALUE_SCALE*2.5, PARAM_GRAPH_VALUE_SCALE*2.5, "T,K", diagram, dIx, dIy);

        saveBMP("GAS_Approximation.bmp", diagram, dIx, dIy);
        delete(diagram);
    }

//    saveBMP("Graph_GAS_Approximation.bmp", graph, Ix, Iy);

    delete(graph);
}

int build(int enter_param)
{
    cout << fixed << setprecision(5);   /// num of numbers after point

    char* name = new char[128];
    int ans;
    if(enter_param == 0){
        cout << "Which file to open?\n1)fluid.txt\n2)fluid H2O.txt\n3)fluid methane.txt\n4)fluid ethane.txt\n5)fluid dodecane.txt\n6)else" << endl;
        cin >> ans;
    }
    else ans = enter_param;

    if(ans == 1) name = "substance/fluid.txt";
    else if(ans == 2) name = "substance/fluid H2O.txt";
    else if(ans == 3) name = "substance/fluid methane.txt";
    else if(ans == 4) name = "substance/fluid ethane.txt";
    else if(ans == 5) name = "substance/fluid dodecane.txt";
    else if(ans == 6) cin >> name;
    cout << endl;

    Fluid fluiddd(name);

    cout << "elem" << "\t" << "Concentration" << "\t" << "Tcrit" << "\t\t" << "Pcrit" << "\t\t" << "acentric" << "\t" << "molar weight" << endl;
    FOR(i, fluiddd.N_comp){
        cout << fluiddd.component[i].comp << "\t" << fluiddd.component[i].C_in_deposit << "\t\t" << fluiddd.component[i].Tc << "\t" << fluiddd.component[i].Pc << "\t" << fluiddd.component[i].acf << "\t\t" << fluiddd.component[i].MW << endl;
    }
    cout << endl;

    if(fluiddd.N_comp > 1)
        ans = 1;

    long double T, p[3];
    int jp[3];
    {
        if(ans != 1){
            cout << "Enter temperature in K: "; cin >> T;
            cout << "Enter p1, p2, p3 in bar, which will use to build approximation: "; cin >> p[0] >> p[1] >> p[2];}
        else{
            int intT;
            cout << "Enter temperature in C: "; cin >> intT;
            if(fluiddd.Load_GAS_and_Liquid_Concentrations(intT) != 0) return 0;
//            fluiddd.print_Concentrations();
            fluiddd.print_pressure_j_from_data();
            T = (long double)intT + 273.15;
            cout << "Enter numbers of pj to take from that list, which will used to build approximation: "; cin >> jp[0] >> jp[1] >> jp[2];
            p[0] = fluiddd.CT_G.p[jp[0]];
            p[1] = fluiddd.CT_G.p[jp[1]];
            p[2] = fluiddd.CT_G.p[jp[2]];
            cout << "p1,p2,p3 = " << p[0] << " " << p[1] << " " << p[2] << endl;
        }
    }

    long double w[3][3];
    int nW[3];

//    if(0){                  /// Решение УРС ПР с классическими коэффициентами a и b (не приведёнными)
//        cout << "\t\tMetric method" << endl << endl;
//        Solution_PR(T, p[0], fluid, fN);
//        Solution_PR(T, p[1], fluid, fN);
//        Solution_PR(T, p[2], fluid, fN);
//    }

    cout << endl << "\t\tNon-metric method" << endl;
    if(ans == 1){
        nW[0] = Solution_PR_auto(w[0], T, fluiddd, jp[0]);
        nW[1] = Solution_PR_auto(w[1], T, fluiddd, jp[1]);
        nW[2] = Solution_PR_auto(w[2], T, fluiddd, jp[2]);
    }
    else{
        nW[0] = Solution_PR_auto(w[0], T, p[0], fluiddd.component[0]);
        nW[1] = Solution_PR_auto(w[1], T, p[1], fluiddd.component[0]);
        nW[2] = Solution_PR_auto(w[2], T, p[2], fluiddd.component[0]);
    }



    cout << endl << "Roots of P-R are:" << endl;
    if(fluiddd.N_comp > 0)
        for(int i = 0; i < 3; i++){
            cout << "for p = " << p[i] << " : ";
            if(nW[i] == 1){
                if(w[i][0] > fluiddd.Wc)   cout << "gas root = " << w[i][0] << endl;
                else                        cout << "liquid root = " << w[i][0] << endl;
            }
            else if(nW[i] == 3) cout << "gas root = " << max3(w[i]) << ", liquid root = " << min3(w[i]) << endl;
            else if(nW[i] == 2){
                    if(w[i][0] > w[i][1]) cout << "gas root = " << w[i][0] << ", liquid root = " << w[i][1] << endl;
                    else                  cout << "gas root = " << w[i][1] << ", liquid root = " << w[i][0] << endl;
            }
        }
    else    {cout << "CRITICAL ERROR!!! Number of components is 0 for no reason!!!" << endl; return -1;}
    cout << endl;


    long double Wg[3], Wl[3];       //  gas and liquid roots
    long double a, b, pz, beta, bg; //  coefficients of approximations

    for(int i = 0; i < 3; i++){
        if(nW[i] == 1){ if(w[i][0] > fluiddd.Wc) Wg[i] = w[i][0];
                        else                     Wl[i] = w[i][0];}
        else if(nW[i] == 3) {Wg[i] = max3(w[i]); Wl[i] = min3(w[i]);}
        else if(nW[i] == 2){
                if(w[i][0] > w[i][1]){ Wg[i] = w[i][0]; Wl[i] = w[i][1];}
                else                 { Wg[i] = w[i][1]; Wl[i] = w[i][0];}
        }
    }

    a = -((Wl[1] - Wl[2])*(p[1] - p[2])*(Wl[0]*Wl[0] - Wl[0]*Wl[1] - Wl[0]*Wl[2] + Wl[1]*Wl[2])*(p[0]*p[0] - p[0]*p[1] - p[0]*p[2] + p[1]*p[2])) / (R*T*(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1])*(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1]));
    b = (Wl[0]*Wl[1]*p[0] - Wl[0]*Wl[1]*p[1] - Wl[0]*Wl[2]*p[0] + Wl[0]*Wl[2]*p[2] + Wl[1]*Wl[2]*p[1] - Wl[1]*Wl[2]*p[2])/(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] -Wl[2]*p[0] + Wl[2]*p[1]);
    pz = (Wl[0]*p[0]*p[1] - Wl[0]*p[0]*p[2] - Wl[1]*p[0]*p[1] + Wl[1]*p[1]*p[2] + Wl[2]*p[0]*p[2] - Wl[2]*p[1]*p[2])/(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1]);

    if(isnan(a) != 1 && isnan(b) != 1 && isnan(pz) != 1){
        cout << endl << "Liquid approximation coefficients (on p1, p2, p3):" << endl;
        cout << "a = " << a << endl << "b = " << b << endl << "p* = " << pz << endl;
    }
    else cout << "-----Couldn't build the liquid approximation. Try with another p1, p2, p3!" << endl;


    beta = (p[0]*p[1]*(Wg[0] - Wg[1])) / (R*T*(p[1] - p[0]));
    bg = (Wg[0]*p[0] - Wg[1]*p[1]) / (p[0] - p[1]);

    if(isnan(beta) != 1 && isnan(bg) != 1){
        cout << endl << "Gas approximation coefficients (on p1, p2):" << endl;
        cout << "beta = " << beta << endl << "bg = " << bg << endl << endl;
    }
    else cout << "-----Couldn't build the gas approximation. Try with another p1, p2!" << endl;


//    long double pg[6] = {1.1, 5., 12., 19., 26., 30.};

//    cout << " liquid approximation roots:" << endl;
//    if(1)   for(int i = 0; i < 6; i++)
//                cout << "with p = " << pg[i] << " Wl = " << ((a*R*T)/(pg[i] + pz) + b) << endl;
//    cout << endl;
////    long double pg[6] = {1., 5., 25., 50., 100., 250.};
//
//    cout << " gas approximation roots:" << endl;
//    if(1)   for(int i = 0; i < 6; i++)
//                cout << "with p = " << pg[i] << " Wg = " << ((beta*R*T)/(pg[i]) + bg) << endl;
//    cout << endl;



                    ///     Подсчёт погрешности аппроксимации
    long double max_G_error = -1, max_L_error = -1;
    if(ans == 1){
        long double roots_W[3], root_L = -1, root_G = -1, root_app_L, root_app_G;
        int roots_n;

        for(int i = 0; i < fluiddd.CT_G.n_p; i++){
            root_app_L = ((a*R*T)/(fluiddd.CT_G.p[i] + pz) + b);
            root_app_G = ((beta*R*T)/(fluiddd.CT_G.p[i]) + bg);
            roots_n = Solution_PR_auto(roots_W, T, fluiddd, i);

            if(roots_n == 1){ if(roots_W[0] > fluiddd.Wc) root_G = roots_W[0];
                            else                     root_L = roots_W[0];}
            else if(roots_n == 3) {root_G = max3(roots_W); root_L = min3(roots_W);}
            else if(roots_n == 2){
                    if(roots_W[0] > roots_W[1]){ root_G = roots_W[0]; root_L = roots_W[1];}
                    else                 { root_G = roots_W[1]; root_L = roots_W[0];}
            }

            if(root_L != -1)
                if(max_L_error < module((root_L - root_app_L) / root_L))
                    max_L_error = module((root_L - root_app_L) / root_L);
            if(root_G != -1)
                if(max_G_error < module((root_G - root_app_G) / root_G))
                    max_G_error = module((root_G - root_app_G) / root_G);
            root_L = -1;
            root_G = -1;
        }
        if(max_L_error != -1)
                cout << "Liquid approximation inaccuracy (in 100%) = " << max_L_error * 100 << "%" << endl;
        if(max_G_error != -1)
                cout << "GAS approximation inaccuracy (in 100%) = " << max_G_error * 100 << "%" << endl;
        cout << endl;
    }
    else{
//        long double pg[6] = {1., 5., 25., 50., 100., 250.};
        long double p_min, p_max, d_pressure, pg;
        int pressure_N;
        cout << "Enter pressure range to calculate the error (min, max): ";
        cin >> p_min >> p_max;
        pressure_N = p_max - p_min;
        if(pressure_N < 1) cout << "The range is too small!" << endl;
        else{
            d_pressure = (p_max - p_min) / pressure_N;
            pg = p_min;
            root_G_L *roots_PR = new root_G_L[pressure_N];

            long double roots_W[3], root_L = -1, root_G = -1, root_app_L, root_app_G;
            int roots_n;

            for(int i = 0; i < pressure_N; i++){
                root_app_L = ((a*R*T)/(pg + pz) + b);
                root_app_G = ((beta*R*T)/(pg) + bg);
                roots_n = Solution_PR_auto(roots_W, T, pg, fluiddd.component[0]);

                if(roots_n == 1){ if(roots_W[0] > fluiddd.Wc) root_G = roots_W[0];
                                else                     root_L = roots_W[0];}
                else if(roots_n == 3) {root_G = max3(roots_W); root_L = min3(roots_W);}
                else if(roots_n == 2){
                        if(roots_W[0] > roots_W[1]){ root_G = roots_W[0]; root_L = roots_W[1];}
                        else                 { root_G = roots_W[1]; root_L = roots_W[0];}
                }

                cout << " p = " << pg << ":" << endl;
                if(root_L != -1){
                    cout << " liquid P-R root = " << root_L << ", approximation = " << root_app_L << endl;
                    if(max_L_error < module((root_L - root_app_L) / root_L))
                        max_L_error = module((root_L - root_app_L) / root_L);
                }
                if(root_G != -1){
                    cout << " gas P-R root = " << root_G << ", approximation = " << root_app_G << endl;
                    if(max_G_error < module((root_G - root_app_G) / root_G))
                        max_G_error = module((root_G - root_app_G) / root_G);
                }
                root_L = -1;
                root_G = -1;
                pg += d_pressure;
            }
            /////////////////////////////////////////////
//            for(int i = 0; i < 6; i++){
//                root_app_L = ((a*R*T)/(pg[i] + pz) + b);
//                root_app_G = ((beta*R*T)/(pg[i]) + bg);
//                roots_n = Solution_PR_auto(roots_W, T, pg[i], fluiddd.component[0]);
//
//                if(roots_n == 1){ if(roots_W[0] > fluiddd.Wc) root_G = roots_W[0];
//                                else                     root_L = roots_W[0];}
//                else if(roots_n == 3) {root_G = max3(roots_W); root_L = min3(roots_W);}
//                else if(roots_n == 2){
//                        if(roots_W[0] > roots_W[1]){ root_G = roots_W[0]; root_L = roots_W[1];}
//                        else                 { root_G = roots_W[1]; root_L = roots_W[0];}
//                }
//
//                cout << " p = " << pg[i] << ":" << endl;
//                if(root_L != -1){
//                    cout << " liquid P-R root = " << root_L << ", approximation = " << root_app_L << endl;
//                    if(max_L_error < module((root_L - root_app_L) / root_L))
//                        max_L_error = module((root_L - root_app_L) / root_L);
//                }
//                if(root_G != -1){
//                    cout << " gas P-R root = " << root_G << ", approximation = " << root_app_G << endl;
//                    if(max_G_error < module((root_G - root_app_G) / root_G))
//                        max_G_error = module((root_G - root_app_G) / root_G);
//                }
//                root_L = -1;
//                root_G = -1;
//            }
            if(max_L_error != -1)
                    cout << "Liquid approximation error (in 100%) = " << max_L_error * 100 << "%" << endl;
            if(max_G_error != -1)
                    cout << "GAS approximation error (in 100%) = " << max_G_error * 100 << "%" << endl;
            cout << endl;
        }
    }




                ///         GRAPH BUILDER       ///
    if(ans == 1){
        int ImageX, ImageY;                                 //  diogramm resolution
        ImageX = 1000;
        ImageY = 1000;

        cout << "Building Peng-Robinson with approximation graph..." << endl;
        Build_PR_Graph_PW_with_approximations(fluiddd, T, a, b, pz, beta, bg, max_G_error, max_L_error, ImageX, ImageY);     /// Peng-Robins PV Graph.

    //    if(isnan(a) != 1 && isnan(b) != 1 && isnan(pz) != 1){           //  is "a, b, pz" is numbers?
    //        cout << "Building Liquid approximation graph..." << endl;
    //        Build_PR_Liquid_Approximation(a, b, pz, T, Gp0, Gp1, GW0, GW1, ImageX, ImageY);
    //    }
    //
    //    if(isnan(beta) != 1 && isnan(bg) != 1){
    //        cout << "Building GAS approximation graph..." << endl;
    //        Build_PR_GAS_Approximation(beta, bg, T, Gp0, Gp1, GW0, GW1, ImageX, ImageY);
    //    }
    }
    else{
        long double Gp0, Gp1, GW0, GW1;                     //  Gp0 Gp1 GW0 GW1 - span of graph (chart range).
        int ImageX, ImageY;                                 //  diogramm resolution
        Gp0 = 1.;
        Gp1 = 250.;
        GW0 = 0.;
        GW1 = 40.;
        ImageX = 1000;
        ImageY = 350;

        cout << "Building PW graph..." << endl;
        Build_PR_Graph_PW(fluiddd, T, Gp0, Gp1, GW0, GW1, ImageX, ImageY);     /// Peng-Robins PV Graph.

        if(isnan(a) != 1 && isnan(b) != 1 && isnan(pz) != 1){           //  is "a, b, pz" is numbers?
            cout << "Building Liquid approximation graph..." << endl;
            Build_PR_Liquid_Approximation(a, b, pz, T, Gp0, Gp1, GW0, GW1, ImageX, ImageY);
        }

        if(isnan(beta) != 1 && isnan(bg) != 1){
            cout << "Building GAS approximation graph..." << endl;
            Build_PR_GAS_Approximation(beta, bg, T, Gp0, Gp1, GW0, GW1, ImageX, ImageY);
        }
    }


    return 0;
}

int build_with_best_approximation_finder()
{
    cout << fixed << setprecision(6);   /// num of numbers after point

    char* name = new char[128];
    int ans;
    cout << "Which file to open?\n1)fluid.txt\n2)fluid H2O.txt\n3)fluid methane.txt\n4)fluid ethane.txt\n5)fluid dodecane.txt\n6)else\n0)exit" << endl;
    cin >> ans;

    if(ans == 1) name = "substance/fluid.txt";
    else if(ans == 2) name = "substance/fluid H2O.txt";
    else if(ans == 3) name = "substance/fluid methane.txt";
    else if(ans == 4) name = "substance/fluid ethane.txt";
    else if(ans == 5) name = "substance/fluid dodecane.txt";
    else if(ans == 6) cin >> name;
    else if(ans == 0) return -1;
    cout << endl;

    Fluid fluiddd(name);

    cout << "elem" << "\t" << "Concentration" << "\t" << "Tcrit" << "\t\t" << "Pcrit" << "\t\t" << "acentric" << "\t" << "molar weight" << endl;
    FOR(i, fluiddd.N_comp){
        cout << fluiddd.component[i].comp << "\t" << fluiddd.component[i].C_in_deposit << "\t\t" << fluiddd.component[i].Tc << "\t" << fluiddd.component[i].Pc << "\t" << fluiddd.component[i].acf << "\t\t" << fluiddd.component[i].MW << endl;
    }
    cout << endl;

    if(fluiddd.N_comp > 1)
        ans = 1;

    long double T;
    {
        if(ans != 1){
            cout << "Enter temperature in K: "; cin >> T;
            cout << "Enter pressure range [min, max] (will be used to calculate the error and build the graph): ";
            cin >> fluiddd.p_min >> fluiddd.p_max;
        }
        else{
            int intT;
            cout << "Enter temperature in C: "; cin >> intT;
            if(fluiddd.Load_GAS_and_Liquid_Concentrations(intT) != 0) return 0;
//            fluiddd.print_Concentrations();
//            fluiddd.print_pressure_j_from_data();
            T = (long double)intT + 273.15;
        }
    }





    cout << endl;
    long double best_a, best_b, best_pz, best_beta, best_bg;
    long double MIN_max_G_error = 373., MIN_max_L_error = 373.;

    cout << "Best approximation finding and max error calculation process, wait..." << endl;
    if(ans == 1){                       /// Best approximation finder for fluid       !!!!!Notice: will only work with Solutions_PR_auto() in return = 3 version
        long double p[3];
        int jp[3];
        long double w[3][3];
        int nW[3];
        long double L_best_p[3], G_best_p[3];       //  for GAS is need only p[0], p[1], but p[3] is still used after for Liquid calculation
        int L_best_jp[3], G_best_jp[3];
        for(jp[0] = 0; jp[0] < fluiddd.CT_G.n_p - 2; jp[0]++){
            for(jp[1] = jp[0]+1; jp[1] < fluiddd.CT_G.n_p - 1; jp[1]++){
                for(jp[2] = jp[1]+1; jp[2] < fluiddd.CT_G.n_p; jp[2]++){    /// скорее всего я не установил jp[2] = jp[1]; после усовершенствования самого расчёта. Строчка выше, по идее, тоже должна быть без "-1".
                    p[0] = fluiddd.CT_G.p[jp[0]];
                    p[1] = fluiddd.CT_G.p[jp[1]];
                    p[2] = fluiddd.CT_G.p[jp[2]];
                    nW[0] = Solution_PR_auto(w[0], T, fluiddd, jp[0]);
                    nW[1] = Solution_PR_auto(w[1], T, fluiddd, jp[1]);
                    nW[2] = Solution_PR_auto(w[2], T, fluiddd, jp[2]);

                    long double Wg[3], Wl[3];       //  gas and liquid roots
                    long double a, b, pz, beta, bg; //  coefficients of approximations

                    for(int i = 0; i < 3; i++){
                        Wg[i] = w[i][0];
                        Wl[i] = w[i][1];
                    }

                    if(p[1] != p[2]){
                        a = -((Wl[1] - Wl[2])*(p[1] - p[2])*(Wl[0]*Wl[0] - Wl[0]*Wl[1] - Wl[0]*Wl[2] + Wl[1]*Wl[2])*(p[0]*p[0] - p[0]*p[1] - p[0]*p[2] + p[1]*p[2])) / (R*T*(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1])*(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1]));
                        b = (Wl[0]*Wl[1]*p[0] - Wl[0]*Wl[1]*p[1] - Wl[0]*Wl[2]*p[0] + Wl[0]*Wl[2]*p[2] + Wl[1]*Wl[2]*p[1] - Wl[1]*Wl[2]*p[2])/(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] -Wl[2]*p[0] + Wl[2]*p[1]);
                        pz = (Wl[0]*p[0]*p[1] - Wl[0]*p[0]*p[2] - Wl[1]*p[0]*p[1] + Wl[1]*p[1]*p[2] + Wl[2]*p[0]*p[2] - Wl[2]*p[1]*p[2])/(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1]);
                    }

                    beta = (p[0]*p[1]*(Wg[0] - Wg[1])) / (R*T*(p[1] - p[0]));
                    bg = (Wg[0]*p[0] - Wg[1]*p[1]) / (p[0] - p[1]);


                    long double roots_W[3], root_L = -1, root_G = -1, root_app_L, root_app_G, max_G_error = -1, max_L_error = -1;
                    int roots_n;

                    for(int i = 0; i < fluiddd.CT_G.n_p; i++){
                        root_app_L = ((a*R*T)/(fluiddd.CT_G.p[i] + pz) + b);
                        root_app_G = ((beta*R*T)/(fluiddd.CT_G.p[i]) + bg);
                        roots_n = Solution_PR_auto(roots_W, T, fluiddd, i);

                        if(roots_n == 3){root_G = roots_W[0]; root_L = roots_W[1];}
                        else cout << "ERROR!!!!!!!!!!!!! USING NOT RIGHT Solutions_PR_Auto() VERSION" << endl;

                        if(p[1] != p[2])
                            if(max_L_error < module((root_L - root_app_L) / root_L))
                                max_L_error = module((root_L - root_app_L) / root_L);
                        if(max_G_error < module((root_G - root_app_G) / root_G))
                            max_G_error = module((root_G - root_app_G) / root_G);
                    }
                    if(max_L_error != -1){
                        if(max_L_error < MIN_max_L_error){
                            MIN_max_L_error = max_L_error;
                            best_a = a;
                            best_b = b;
                            best_pz = pz;
                            L_best_p[0] = p[0];
                            L_best_p[1] = p[1];
                            L_best_p[2] = p[2];
                            L_best_jp[0] = jp[0];
                            L_best_jp[1] = jp[1];
                            L_best_jp[2] = jp[2];
                            if(IS_DEBUGG == 1)
                                cout << "Liquid approximation inaccuracy (in 100%) = " << max_L_error * 100 << "%" << endl;
                        }
                    }
                    if(max_G_error != -1){
                        if(max_G_error < MIN_max_G_error){
                            MIN_max_G_error = max_G_error;
                            best_beta = beta;
                            best_bg = bg;
                            G_best_p[0] = p[0];
                            G_best_p[1] = p[1];
                            G_best_p[2] = p[2];
                            G_best_jp[0] = jp[0];
                            G_best_jp[1] = jp[1];
                            G_best_jp[2] = jp[2];
                            if(IS_DEBUGG == 1)
                                cout << "GAS approximation inaccuracy (in 100%) = " << max_G_error * 100 << "%" << endl;
                        }
                    }

    //                if(jp[1] == fluiddd.CT_G.n_p - 2){
    //                    do one more for(jp[1] = fluiddd.CT_G.n_p - 1 only for GAS){
    //
    //                    }
    //                }
                }
            }
        }

        cout << endl;
        cout << "Best Luqiud approximation on pi = " << L_best_p[0] << " " << L_best_p[1] << " " << L_best_p[2] << "\nApprox. coefficients are:\na = " << best_a << "\nb = " << best_b << "\np* = " << best_pz << endl << "inaccuracy (in 100%) = " << MIN_max_L_error * 100 << "%" << endl;
        cout << "L jp: " << L_best_jp[0] << " " << L_best_jp[1] << " " << L_best_jp[2] << endl << endl;
        cout << "Best GAS approximation on pi = " << G_best_p[0] << " " << G_best_p[1] << "\nApprox. coefficients are:\nbeta = " << best_beta << "\nbg = " << best_bg << endl << "inaccuracy (in 100%) = " << MIN_max_G_error * 100 << "%" << endl;
        cout << "G jp: " << G_best_jp[0] << " " << G_best_jp[1] << endl << endl;
    }
    else{                               /// Best approximation finder for pure substance
        long double p[3];
        long double L_best_p[3], G_best_p[3];       //  for GAS is need only p[0], p[1], but p[3] is still used after for Liquid calculation
        long double d_pressure, pg, roots_W[3];
        int pressure_N, roots_n;
//        cout << "Enter pressure range to calculate the error (min, max): ";
//        cin >> p_min >> p_max;
        pressure_N = fluiddd.p_max - fluiddd.p_min;
        if(pressure_N < 1){cout << "The range is too small!" << endl;}
        else{
            d_pressure = (fluiddd.p_max - fluiddd.p_min) / pressure_N;
            pressure_N++;
            pg = fluiddd.p_min;
            root_G_L *roots_PR = new root_G_L[pressure_N];

            for(int i = 0; i < pressure_N; i++){
                roots_PR[i].p = pg;
                roots_n = Solution_PR_auto(roots_W, T, pg, fluiddd.component[0]);
                if(roots_n == 1){ if(roots_W[0] > fluiddd.Wc) roots_PR[i].G = roots_W[0];
                                else                     roots_PR[i].L = roots_W[0];}
                else if(roots_n == 3) {roots_PR[i].G = max3(roots_W); roots_PR[i].L = min3(roots_W);}
                else if(roots_n == 2){
                        if(roots_W[0] > roots_W[1]){ roots_PR[i].G = roots_W[0]; roots_PR[i].L = roots_W[1];}
                        else                 { roots_PR[i].G = roots_W[1]; roots_PR[i].L = roots_W[0];}
                }
                pg += d_pressure;
            }

            long double Wg[3], Wl[3];       //  gas and liquid roots
            long double a, b, pz, beta, bg; //  coefficients of approximations
            long double roots_error, root_G = -1, max_G_error = -1, max_L_error = -1;
//            auto Tm_Start = steady_clock::now();
            for(int i = 0; i < pressure_N-1; i++){
                for(int j = i+1; j < pressure_N; j++){
                    for(int k = j; k < pressure_N; k++){
                        p[0] = roots_PR[i].p;
                        p[1] = roots_PR[j].p;
                        p[2] = roots_PR[k].p;


                            Wg[0] = roots_PR[i].G;
                            Wl[0] = roots_PR[i].L;
                            Wg[1] = roots_PR[j].G;
                            Wl[1] = roots_PR[j].L;
                            Wg[2] = roots_PR[k].G;
                            Wl[2] = roots_PR[k].L;

                        if(p[1] != p[2] && Wl[0] != -1 && Wl[1] != -1 && Wl[2] != -1){
                            a = -((Wl[1] - Wl[2])*(p[1] - p[2])*(Wl[0]*Wl[0] - Wl[0]*Wl[1] - Wl[0]*Wl[2] + Wl[1]*Wl[2])*(p[0]*p[0] - p[0]*p[1] - p[0]*p[2] + p[1]*p[2])) / (R*T*(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1])*(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1]));
                            b = (Wl[0]*Wl[1]*p[0] - Wl[0]*Wl[1]*p[1] - Wl[0]*Wl[2]*p[0] + Wl[0]*Wl[2]*p[2] + Wl[1]*Wl[2]*p[1] - Wl[1]*Wl[2]*p[2])/(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] -Wl[2]*p[0] + Wl[2]*p[1]);
                            pz = (Wl[0]*p[0]*p[1] - Wl[0]*p[0]*p[2] - Wl[1]*p[0]*p[1] + Wl[1]*p[1]*p[2] + Wl[2]*p[0]*p[2] - Wl[2]*p[1]*p[2])/(-Wl[0]*p[1] + Wl[0]*p[2] + Wl[1]*p[0] - Wl[1]*p[2] - Wl[2]*p[0] + Wl[2]*p[1]);
                        }

                        if(Wg[0] != -1 && Wg[1] != -1){
                            beta = (p[0]*p[1]*(Wg[0] - Wg[1])) / (R*T*(p[1] - p[0]));
                            bg = (Wg[0]*p[0] - Wg[1]*p[1]) / (p[0] - p[1]);
                        }

                        max_G_error = -1;
                        max_L_error = -1;

                        for(int ijk = 0; ijk < pressure_N; ijk++){
                            if(p[1] != p[2] && Wl[0] != -1 && Wl[1] != -1 && Wl[2] != -1 && roots_PR[ijk].L != -1){
                                roots_error = module((roots_PR[ijk].L - ((a*R*T)/(roots_PR[ijk].p + pz) + b)) / roots_PR[ijk].L);
                                if(max_L_error < roots_error)
                                    max_L_error = roots_error;
                            }
                            if(Wg[0] != -1 && Wg[1] != -1 && roots_PR[ijk].G != -1){
                                roots_error = module((roots_PR[ijk].G - ((beta*R*T)/(roots_PR[ijk].p) + bg)) / roots_PR[ijk].G);
                                if(max_G_error < roots_error)
                                    max_G_error = roots_error;
                            }
                        }
                        if(max_L_error != -1){
                            if(max_L_error < MIN_max_L_error){
                                MIN_max_L_error = max_L_error;
                                best_a = a;
                                best_b = b;
                                best_pz = pz;
                                L_best_p[0] = p[0];
                                L_best_p[1] = p[1];
                                L_best_p[2] = p[2];
                                if(IS_DEBUGG == 1)
                                    cout << "Liquid approximation inaccuracy (in 100%) = " << max_L_error * 100 << "%" << endl;
                            }
                        }
                        if(max_G_error != -1){
                            if(max_G_error < MIN_max_G_error){
                                MIN_max_G_error = max_G_error;
                                best_beta = beta;
                                best_bg = bg;
                                G_best_p[0] = p[0];
                                G_best_p[1] = p[1];
                                G_best_p[2] = p[2];
                                if(IS_DEBUGG == 1)
                                    cout << "GAS approximation inaccuracy (in 100%) = " << max_G_error * 100 << "%" << endl;
                            }
                        }
                    }
                }
            }

//    auto Tm_End = steady_clock::now();
//    cout << "Elapsed time: " << duration_cast<microseconds>(Tm_End - Tm_Start).count()/1000000. << endl;
            cout << endl;
            if(MIN_max_L_error != 373.)
                cout << "Best Luqiud approximation on pi = " << L_best_p[0] << " " << L_best_p[1] << " " << L_best_p[2] << "\nApprox. coefficients are:\na = " << best_a << "\nb = " << best_b << "\np* = " << best_pz << endl << "inaccuracy (in 100%) = " << MIN_max_L_error * 100 << "%" << endl << endl;
            if(MIN_max_G_error != 373.)
                cout << "Best GAS approximation on pi = " << G_best_p[0] << " " << G_best_p[1] << "\nApprox. coefficients are:\nbeta = " << best_beta << "\nbg = " << best_bg << endl << "inaccuracy (in 100%) = " << MIN_max_G_error * 100 << "%" << endl << endl;
        }
    }
    cout << "______________________________________________________" << endl;

if(0) if(ans != 1){   /// исключительно для курсовой - вывод значений в 6 точках
        long double pg[6]/* = {16., 50., 100., 150., 200., 250.}*/;
        cout << "Enter p[i]. i from 0 to 6: ";
        for(int i = 0; i < 6; i++)
            cin >> pg[i];
    {   /// LIQUID
        int nW[6];

        long double Wg[6], Wl[6], w[6][3];       //  gas and liquid roots
        for(int i = 0; i < 6; i++){
//            cin >> pg[i];
            nW[i] = Solution_PR_auto(w[i], T, pg[i], fluiddd.component[0]);
        }

        for(int i = 0; i < 6; i++){
            if(nW[i] == 1){ if(w[i][0] > fluiddd.Wc) Wg[i] = w[i][0];
                            else                     Wl[i] = w[i][0];}
            else if(nW[i] == 3) {Wg[i] = max3(w[i]); Wl[i] = min3(w[i]);}
            else if(nW[i] == 2){
                    if(w[i][0] > w[i][1]){ Wg[i] = w[i][0]; Wl[i] = w[i][1];}
                    else                 { Wg[i] = w[i][1]; Wl[i] = w[i][0];}
            }
        }

        cout << " liquid roots:" << endl;
        for(int i = 0; i < 6; i++)
            cout << "with p = " << pg[i] << " apprx = " << ((best_a*R*T)/(pg[i] + best_pz) + best_b) << ", P-R = " << Wl[i] << endl;
        cout << endl;
    }
    {   /// GAS
//        long double pg[6] = {1., 5., 25., 50., 150., 250.};
        int nW[6];

        long double Wg[6], Wl[6], w[6][3];       //  gas and liquid roots
        for(int i = 0; i < 6; i++)
            nW[i] = Solution_PR_auto(w[i], T, pg[i], fluiddd.component[0]);

        for(int i = 0; i < 6; i++){
            if(nW[i] == 1){ if(w[i][0] > fluiddd.Wc) Wg[i] = w[i][0];
                            else                     Wl[i] = w[i][0];}
            else if(nW[i] == 3) {Wg[i] = max3(w[i]); Wl[i] = min3(w[i]);}
            else if(nW[i] == 2){
                    if(w[i][0] > w[i][1]){ Wg[i] = w[i][0]; Wl[i] = w[i][1];}
                    else                 { Wg[i] = w[i][1]; Wl[i] = w[i][0];}
            }
        }

        cout << " gas approximation roots:" << endl;
        for(int i = 0; i < 6; i++)
            cout << "with p = " << pg[i] << " Wg = " << ((best_beta*R*T)/(pg[i]) + best_bg) << ", P-R = " << Wg[i] << endl;
        cout << endl;
    }

///////////////////
}

                ///         GRAPH BUILDER       ///
    {
        int ImageX, ImageY;                                 //  diogramm resolution. 1600 500 as NIST(150%) //  930 620 fluid // 1000 500 pure substance
        ImageX = 930;
        ImageY = 620;

        cout << "Building Peng-Robinson with approximation graph..." << endl;
        Build_PR_Graph_PW_with_approximations(fluiddd, T, best_a, best_b, best_pz, best_beta, best_bg, MIN_max_G_error, MIN_max_L_error, ImageX, ImageY);     /// Peng-Robins PV Graph.
    }

    return 0;
}

int main()
{
    cout << "Which volume to use: 0 = Molar, 1 = Specific?" << endl;
    cin >> PARAM_MOLAR_OR_SPECIFIC_VOLUME;
    cout << endl;
//    if(IS_DEBUGG == 0)
//        build_with_best_approximation_finder();
//    else
        while(1){
            if(build_with_best_approximation_finder() == -1) break;
            cout << endl << endl << endl;
        }
    cout << "End of program. Press Enter to exit." << endl;
    getchar();
    getchar();
    return 0;
}










