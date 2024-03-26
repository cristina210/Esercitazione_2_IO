#include "iostream"
#include "Eigen/Eigen"
// "err_rel_qr" è una funzione che calcola la soluzione di un sistema lineare utilizzando la decomposizione QR e successivamente l'errore relativo
double err_rel_qr(Eigen::MatrixXd A,Eigen::Vector2d b){
    Eigen::Matrix<double, 2, 1> sol_ex = {-1.0e+0, -1.0e+00};
    Eigen::Vector2d xQR = A.householderQr().solve(b); // in https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR.html
    double errRel=(xQR-sol_ex).norm()/sol_ex.norm();
    return errRel;
}
// "err_rel_palu" è una funzione che calcola la soluzione di un sistema lineare utilizzando la decomposizione PALU e successivamente l'errore relativo
double err_rel_palu(Eigen::MatrixXd A,Eigen::Vector2d b){
    Eigen::Matrix<double, 2, 1> sol_ex = {-1.0e+0, -1.0e+00};
    Eigen::Vector2d xQR = A.partialPivLu().solve(b);  // in https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html
    double errRel=(xQR-sol_ex).norm();
    return errRel;
}
int main()
{   std::cout<<"Errori relativi nella risoluzione del sistema lineare ->\n";
    Eigen::MatrixXd A1  {
        {5.547001962252291e-01, -3.770900990025203e-02},
        {8.320502943378437e-01, -9.992887623566787e-01}
    };
    Eigen::Matrix<double, 2, 1> b1 = {-5.169911863249772e-01, 1.672384680188350e-01};
    double err_qr_1 = err_rel_qr(A1,b1);
    std::cout<<"della prima matrice con la decomposizione QR:   "<<err_qr_1<<std::endl;
    double err_palu_1 = err_rel_palu(A1,b1);
    std::cout<<"della prima matrice con la decomposizione PALU:  "<<err_palu_1<<std::endl;
    Eigen::MatrixXd A2  {
        {5.547001962252291e-01, -5.540607316466765e-01},
        {8.320502943378437e-01, -8.324762492991313e-01}
    };
    Eigen::Matrix<double, 2, 1> b2 = {-6.394645785530173e-04, 4.259549612877223e-04};
    double err_qr_2 = err_rel_qr(A2,b2);
    std::cout<<"della seconda matrice con la decomposizione QR:  "<<err_qr_2<<std::endl;
    double err_palu_2 = err_rel_palu(A2,b2);
    std::cout<<"della seconda matrice con la decomposizione PALU:  "<<err_palu_2<<std::endl;
    Eigen::MatrixXd A3  {
        {5.547001962252291e-01, -5.547001955851905e-01},
        { 8.320502943378437e-01,-8.320502947645361e-01}
    };
    Eigen::Matrix<double, 2, 1> b3 = {-6.400391328043042e-10, 4.266924591433963e-10};
    double err_qr_3 = err_rel_qr(A3,b3);
    std::cout<<"della terza matrice con la decomposizione QR:  "<<err_qr_3<<std::endl;
    double err_palu_3 = err_rel_palu(A3,b3);
    std::cout<<"della terza matrice con la decomposizione PALU:  "<<err_palu_3<<std::endl;
    return 0;
}
