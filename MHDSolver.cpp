//
// Created by Иван on 4/29/2024.
//
#include "MHDSolver.h"

/*          0      1      2      3    4   5   6   7
 * state:  rho,  rho*u, rho*v, rho*w, e, Bx, Bz, By
 * */

// Давление
double pressure(const double &gam_hcr, const double &e, const double &rho, const double &u, const double &v, const double &w, const double &Bx, const double &By, const double &Bz){
    return (gam_hcr - 1) * (e - 0.5 * rho * (u * u + v * v + w * w) - 0.5*(Bx * Bx + By * By + Bz * Bz));
}

// Энергия
double energy(const double &gam_hcr, const double &p, const double &rho, const double &u, const double &v, const double &w, const double &Bx, const double &By, const double &Bz){
    return p/(gam_hcr-1) + 0.5*rho*(u * u + v * v + w * w) + 0.5 * (Bx * Bx + By * By + Bz * Bz);
}

// Суммарное давление
double ptotal(const double &p, const double &Bx, const double &By, const double &Bz) {
    return p + 0.5 * (Bx*Bx + By*By + Bz*Bz);
}

// Вектор состояния из параметров
std::vector<double>
state_from_primitive_vars(const double &rho, const double &u, const double &v, const double &w, const double &p,
                          const double &Bx, const double &By, const double &Bz, const double &gam_hcr) {
    std::vector<double> U(8,0);

    double mx = rho * u;
    double my = rho * v;
    double mz = rho * w;
    double e = energy(gam_hcr, p, rho, u, v, w, Bx, By, Bz);

    U[0] = rho;
    U[1] = mx;
    U[2] = my;
    U[3] = mz;
    U[4] = e;
    U[5] = Bx;
    U[6] = By;
    U[7] = Bz;

    return U;
}

// Определяем МГД-поток
std::vector<double> MHD_flux(const std::vector<double>& U, const double &gam_hcr) {
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];
    double Bx = U[5];
    double By = U[6];
    double Bz = U[7];

    double p = pressure(gam_hcr, e, rho, u, v, w, Bx, By, Bz);
    double pT = ptotal(p, Bx, By, Bz);

    std::vector<double> F(8, 0);
    F[0] = rho * u;
    F[1] = rho * u * u + pT - Bx * Bx /*/(4*PI)*/;
    F[2] = rho * v * u - Bx * By /*/ (4 * PI)*/;
    F[3] = rho * w * u - Bx * Bz /*/ (4 * PI)*/;
    F[4] = (e + pT) * u - Bx * (u * Bx + v * By + w * Bz) /*/ (4 * PI)*/;
    F[5] = 0;
    F[6] = By * u - Bx * v;
    F[7] = Bz * u - Bx * w;

    return F;
}

// Скорость быстрой волны в точке
double cfast(const std::vector<double>& U, const double& gam_hcr) {
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];
    double Bx = U[5];
    double By = U[6];
    double Bz = U[7];

    //|B|^2
    double BB = Bx * Bx + By * By + Bz * Bz;
    //p
    double p = pressure(gam_hcr, e, rho, u, v, w, Bx, By, Bz);

    return std::sqrt((gam_hcr * p + BB + std::sqrt((gam_hcr * p + BB) * (gam_hcr * p + BB) - 4 * gam_hcr * p * Bx * Bx)) / (2 * rho));
}

// Определяем HLL поток F
std::vector<double> HLL_flux(const std::vector<double>& U_L, const std::vector<double>& U_R, const double &gam_hcr) {
    /*          0      1      2      3    4   5   6   7
     * state:  rho,  rho*u, rho*v, rho*w, e, Bx, Bz, By
     */
    double rho_L = U_L[0];
    double u_L = U_L[1]/rho_L;
    double v_L = U_L[2]/rho_L;
    double w_L = U_L[3]/rho_L;
    double e_L = U_L[4];
    double Bx_L = U_L[5];
    double By_L = U_L[6];
    double Bz_L = U_L[7];

    double rho_R = U_R[0];
    double u_R = U_R[1]/rho_R;
    double v_R = U_R[2]/rho_R;
    double w_R = U_R[3]/rho_R;
    double e_R = U_R[4];
    double Bx_R = U_R[5];
    double By_R = U_R[6];
    double Bz_R = U_R[7];

    //быстрые магнитозвуковые скорости на левом и правом концах
    double cf_L = cfast(U_L, gam_hcr);
    double cf_R = cfast(U_R, gam_hcr);

    //скорость левого сигнала, рассчитываемая как минимальное значение скорости левого состояния (uL) и быстрой магнитозвуковой скорости (cfL).
    double SL = std::min(u_L - cf_L, u_R - cf_R);//std::min(uL, uR) - std::max(cfL, cfR);//
    // Скорость правого сигнала, рассчитываемая как максимальное значение скорости правого состояния (uR) и быстрой магнитозвуковой скорости (cfR).ы
    double SR = std::max(u_L + cf_L, u_R + cf_R);

    if (SL > 0) {
        return MHD_flux(U_L, gam_hcr);
    }
    else if (SL <= 0 && SR >= 0) {
        return 1 / (SR - SL) * (SR * MHD_flux(U_L, gam_hcr) - SL * MHD_flux(U_R, gam_hcr) + SL * SR * (U_R - U_L));
    }
    //if (SR < 0)
    else  {
        return MHD_flux(U_R, gam_hcr);
    }
}

/* Инициализация начального состояния системы (каждой точке пространства ставим в соотвествие величину начального отклонения)
 * args: problem
 * return: initial state
 * */
std::vector<std::vector<double>> initializeState(const MHDProblem &problem){
    std::vector<std::vector<double>> new_state(problem.num_space_steps+1, std::vector<double>(8,0));
    double x_i = problem.x0;

    for(int i = 0; i <= problem.num_space_steps; ++i){
        new_state[i] = problem.initStateFunc(x_i);
        x_i += problem.h;
    }

    return new_state;
}

bool HLLScheme(const MHDProblem &problem, const std::string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting HLL Scheme..." << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        // Физические параметры
        double gam_hcr = problem.gam_hcr;
        double x0 = problem.x0;    // начало отсчёта по пространству
        double L = problem.L;      // характерный пространственный размер (длина струны)
        double X = x0 + L;         // координата правого края
        double t0 = problem.t0;    // время отсчёта
        double T = problem.T;      // время окончания
        double tau = problem.tau;  // шаг по времени
        double h = problem.h;      // шаг по пространству
        double gam_courant = problem.gam_courant;  // число Куранта

        // Шаги по времени и пространству
        int num_time_steps = problem.num_time_steps;
        int num_space_steps = problem.num_space_steps;

        double t_i = t0;
        double x_i = x0;

        std::vector<std::vector<double>> state_j = initializeState(problem);
        std::vector<std::vector<double>> state_jp(state_j);
        std::vector<std::vector<double>> fluxes(num_space_steps, std::vector<double>(8,0));

        // Запись первого слоя в файл
        fpoints << t_i << std::endl;
        for (int i = 0; i <= num_space_steps; ++i)
        {
            writeVectorToFile(fpoints, state_j[i]);
        }
        fpoints << "--------" << std::endl;

        // Эволюция системы во времени
        for(int j = 1; j <= num_time_steps; ++j) {
            t_i += tau;

            // Граничные условия слева
            state_jp[0] = problem.leftBoundaryFunction(t_i);

            // Граничные условия справа
            state_jp[num_space_steps] = problem.rightBoundaryFunction(t_i);

            // HLL-потоки
            for(int i = 0; i < num_space_steps; ++i){
                fluxes[i] = HLL_flux(state_j[i], state_j[i+1], gam_hcr);
            }

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                state_jp[i] = state_j[i] - tau/h * (fluxes[i]-fluxes[i-1]);
            }
            // Запись в файл
            fpoints << t_i << std::endl;
            for (int i = 0; i <= num_space_steps; ++i)
            {
                writeVectorToFile(fpoints, state_jp[i]);
            }
            fpoints << "--------" << std::endl;

            // j  jp   -> jp j
            state_j.swap(state_jp);
        }
        fpoints.close();

        std::cout << "LOG[INFO]: HLL-calculation is completed." << std::endl;

        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file." << std::endl;
        return false;
    }
};

bool HLLCScheme(const MHDProblem &problem, const std::string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting HLLC Scheme..." << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        // Физические параметры
        double gam_hcr = problem.gam_hcr;
        double x0 = problem.x0;    // начало отсчёта по пространству
        double L = problem.L;      // характерный пространственный размер (длина струны)
        double X = x0 + L;         // координата правого края
        double t0 = problem.t0;    // время отсчёта
        double T = problem.T;      // время окончания
        double tau = problem.tau;  // шаг по времени
        double h = problem.h;      // шаг по пространству
        double gam_courant = problem.gam_courant;  // число Куранта

        // Шаги по времени и пространству
        int num_time_steps = problem.num_time_steps;
        int num_space_steps = problem.num_space_steps;

        double t_i = t0;
        double x_i = x0;

        std::vector<std::vector<double>> state_j = initializeState(problem);
        std::vector<std::vector<double>> state_jp(state_j);
        std::vector<std::vector<double>> fluxes(num_space_steps, std::vector<double>(8,0));

        // Запись первого слоя в файл
        fpoints << t_i << std::endl;
        for (int i = 0; i <= num_space_steps; ++i)
        {
            writeVectorToFile(fpoints, state_j[i]);
        }
        fpoints << "--------" << std::endl;

        // Эволюция системы во времени
        for(int j = 1; j <= num_time_steps; ++j) {
            t_i += tau;

            // Граничные условия слева
            state_jp[0] = problem.leftBoundaryFunction(t_i);

            // Граничные условия справа
            state_jp[num_space_steps] = problem.rightBoundaryFunction(t_i);

            // HLL-потоки
            for(int i = 0; i < num_space_steps; ++i){
                fluxes[i] = HLL_flux(state_j[i], state_j[i+1], gam_hcr);
            }

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                state_jp[i] = state_j[i] - tau/h * (fluxes[i]-fluxes[i-1]);
            }
            // Запись в файл
            fpoints << t_i << std::endl;
            for (int i = 0; i <= num_space_steps; ++i)
            {
                writeVectorToFile(fpoints, state_jp[i]);
            }
            fpoints << "--------" << std::endl;

            // j  jp   -> jp j
            state_j.swap(state_jp);
        }
        fpoints.close();

        std::cout << "LOG[INFO]: HLLC-calculation is completed." << std::endl;

        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file." << std::endl;
        return false;
    }
};