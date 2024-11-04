//
// Created by Иван on 10/21/2024.
//

#include "GasSolver2D.h"

/*            0      1      2      3    4
 * state U:  rho,  rho*u, rho*v, rho*w, e
 * */

/**
 * pressure
 * @arg gam_hcr (c_p / c_v  or heat capacity ratio)
 * @arg e (energy component)
 * @arg rho (mass density)
 * @arg u (x-velocity)
 * @arg v (y-velocity)
 * @arg w (z-velocity)
 */
double gasPressure(const double &gam_hcr, const double &e, const double &rho, const double &u, const double &v, const double &w){
    return (gam_hcr - 1) * (e - 0.5 * rho * (u * u + v * v + w * w));
}

double gasPressure(const std::vector<double> &U, const double &gam_hcr){
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];
    return (gam_hcr - 1) * (e - 0.5 * rho * (u * u + v * v + w * w) );
}

double cSound(const std::vector<double> & U, const double & gam_hcr){
    return std::sqrt(gam_hcr * gasPressure(U, gam_hcr) / U[0]);
}

/**
 * F-flux is a flux along the x-axis
 * @arg U (state)
 * @arg gam_hcr (c_p / c_v  or heat capacity ratio)
 */
std::vector<double> F_flux(const std::vector<double>& U, const double &gam_hcr) {
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];

    double p = gasPressure(gam_hcr, e, rho, u, v, w);

    std::vector<double> F(5, 0);
    F[0] = rho * u;
    F[1] = rho * u * u + p;
    F[2] = rho * v * u;
    F[3] = rho * w * u;
    F[4] = (e + p) * u;

    return F;
}

/**
 * G-flux is a flux along the y-axis
 * @arg U (state)
 * @arg gam_hcr (c_p / c_v  or heat capacity ratio)
 */
std::vector<double> G_flux(const std::vector<double>& U, const double &gam_hcr) {
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];

    double p = gasPressure(gam_hcr, e, rho, u, v, w);

    std::vector<double> G(5, 0);
    G[0] = rho * v;
    G[1] = rho * u * v;
    G[2] = rho * v * v + p;
    G[3] = rho * v * w;
    G[4] = (e + p) * v;

    return G;
}

/**
 * HLL-X flux is a HLL flux along the x-axis
 * @arg U_L (left state)
 * @arg U_R (right state)
 * @arg gam_hcr (c_p / c_v  or heat capacity ratio)
 */
std::vector<double> HLL_fluxFx(const std::vector<double>& U_L, const std::vector<double>& U_R, const double &gam_hcr) {
    /*          0      1      2      3    4
     * state:  rho,  rho*u, rho*v, rho*w, e
     */
    double rho_L = U_L[0];
    double u_L = U_L[1]/rho_L;

    double rho_R = U_R[0];
    double u_R = U_R[1]/rho_R;

    //sound speeds
    double cf_L = cSound(U_L, gam_hcr);
    double cf_R = cSound(U_R, gam_hcr);

    //
    double SL = std::min(u_L - cf_L, u_R - cf_R);//std::min(uL, uR) - std::max(cfL, cfR);//
    //
    double SR = std::max(u_L + cf_L, u_R + cf_R);

    if (SL > 0) {
        return F_flux(U_L, gam_hcr);
    }
    else if (SL <= 0 && SR >= 0) {
        return 1 / (SR - SL) * (SR * F_flux(U_L, gam_hcr) - SL * F_flux(U_R, gam_hcr) + SL * SR * (U_R - U_L));
    }
        //if (SR < 0)
    else  {
        return F_flux(U_R, gam_hcr);
    }
}

/**
 * HLL-Y flux is a HLL flux along the y-axis
 * @arg U_L (left state)
 * @arg U_R (right state)
 * @arg gam_hcr (c_p / c_v  or heat capacity ratio)
 */
std::vector<double> HLL_fluxGy(const std::vector<double>& U_L, const std::vector<double>& U_R, const double &gam_hcr) {
    /*          0      1      2      3    4
     * state:  rho,  rho*u, rho*v, rho*w, e
     */
    double rho_L = U_L[0];
    double v_L = U_L[2]/rho_L;

    double rho_R = U_R[0];
    double v_R = U_R[2]/rho_R;

    //sound speeds
    double cf_L = cSound(U_L, gam_hcr);
    double cf_R = cSound(U_R, gam_hcr);

    //
    double SL = std::min(v_L - cf_L, v_R - cf_R);
    //
    double SR = std::max(v_L + cf_L, v_R + cf_R);

    if (SL > 0) {
        return G_flux(U_L, gam_hcr);
    }
    else if (SL <= 0 && SR >= 0) {
        return 1 / (SR - SL) * (SR * G_flux(U_L, gam_hcr) - SL * G_flux(U_R, gam_hcr) + SL * SR * (U_R - U_L));
    }
        //if (SR < 0)
    else  {
        return G_flux(U_R, gam_hcr);
    }
}

// Энергия
double energy(const double &gam_hcr, const double &p, const double &rho, const double &u, const double &v, const double &w){
    return p/(gam_hcr-1) + 0.5*rho*(u * u + v * v + w * w);
}
// Вектор состояния из параметров
std::vector<double>
state_from_primitive_vars(const double &rho, const double &u, const double &v, const double &w, const double &p, const double &gam_hcr) {
    std::vector<double> U(5,0.);

    double mx = rho * u;
    double my = rho * v;
    double mz = rho * w;
    double e = energy(gam_hcr, p, rho, u, v, w);

    U[0] = rho;
    U[1] = mx;
    U[2] = my;
    U[3] = mz;
    U[4] = e;
    return U;
}

void solverHLL2D(const World& world){
    NodePool np = world.getNodePool();
    ElementPool ep = world.getElementPool();
    std::vector<std::vector<double>> states(np.nodeCount, std::vector<double>(5, 0.0));
    /*
     *     states[i] <---> np.nodes[i]
     *
     * */

    //test1:   [0, 1] x [0, 1]
    /*   quadrant            rho        u         v         p
     *   x>0.5, y>0.5:          0.5313     0.0       0.0       0.4
     *   x<0.5, y>0.5:          1.0        0.7276    0.0       1.0
     *   x<0.5, y<0.5:          0.8        0.0       0.0       1.0
     *   x>0.5, y<0.5:          1.0        0.0       0.7276    1.0
     *   FREE-FLOW boundary conditions
     * */
    //printf("elCount = %d", ep.elCount);

    //initializer
    for(int i = 0; i < np.nodeCount; ++i){
        Node tmp_node = np.nodes[i];
        if(tmp_node.x < 0.5 && tmp_node.y > 0.5){
            states[i] = state_from_primitive_vars();
        }
        else if(tmp_node.x > 0.5 && tmp_node.y > 0.5){

        }
        else if(tmp_node.x < 0.5 && tmp_node.y < 0.5){

        }
        else{

        }
    }


    //solver
    for(int i = 0; i < ep.elCount; ++i){
        //пока делаем для n=4, для n=3 попозже...
        std::vector<int> nodeInds = ep.elements[i].nodeIndexes;

    }

}