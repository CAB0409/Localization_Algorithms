#include <iostream>
#include <armadillo>
#include "typeExg_arma_eig.h"
#include "engine.h"

using namespace arma;
using namespace std;

void print_vec_dimensions(vec input) {
    std::cout << "size: " << input.size() << std::endl;
    std::cout << "columns: " << input.n_cols << std::endl;
    std::cout << "rows: " << input.n_rows << std::endl;
}

void print_mat_dimensions(mat input) {
    std::cout << "size: " << input.size() << std::endl;
    std::cout << "columns: " << input.n_cols << std::endl;
    std::cout << "rows: " << input.n_rows << std::endl;
}
/*
    Convert armadillo matrix to eigen to take a more
    accuate inverse.  Converts matrix back to armadillo
*/
template <typename T>
void eigen_inv(arma::Mat<T> &matOut) {
    //Init Eigen Matrix
    Eigen::MatrixXd eigen_mat(matOut.size(),matOut.size());
    //Convert armadillo matrix to eigen
    arma2eigen(matOut, eigen_mat);
    //Take inverse and replace old matrix
    eigen_mat = eigen_mat.inverse();
    //Convert eigen matrix back to armadillo matrix
    eigen2arma(eigen_mat,matOut);
}

mat TDOALocCRLB(mat S, vec u, mat Q) {

    int M = Q.size() + 1;

    //error check here for M

    M = S.n_rows;
    //ro = sqrt(sum((u*ones(1,M)-S).^2));
    rowvec ro = sqrt(sum(pow(u*ones(1,M)-S.t(),2),0));

    S = S.t();
    mat d_u = (S.cols(1,S.n_cols-1)-u*ones(1,M-1)) / (ro(span(1,ro.size()-1)).t()*ones(1,S.n_rows)).t();
    d_u = d_u.t() - (ones(M-1,1)*((S.col(0)-u).t()/ro(0)));

    mat J = d_u.t() * inv(Q) * d_u;

    return inv(J);

}

//SensorPositions, r, Q
//TODO input Q is not correct
vec TDOALoc(mat S, vec r, mat Q) {

    int RptCnt = 3; // number of repetitions in stage 1 to recompute W1
    int M = r.size() + 1;

    //R=sqrt(sum(S.^2))';
    mat cum_r = cumsum(pow(S,2).t());
    mat subrow_r = cum_r.row(1);
    //correct
    vec R = vectorise(sqrt(subrow_r));//element wise sqrt and convert matrix to vector

    //h1 = r.^2 - R(2:end).^2 + R(1)^2;
    //close enough
    vec h1 = pow(r,2) - pow(R.rows(1,R.n_rows-1),2) + pow(R(0),2);

    //=========== construct related vector and matrix ============
    //G1 = -2*[S(:,2:end)'-ones(M-1,1)*S(:,1)' ,  r];
    //correct
    mat G1 = join_rows(-2*(S.submat(1,0,S.n_rows-1, S.n_cols-1)-ones(M-1,1)*S.row(0)), -2*r);

    //============= first stage ===================================  
    //TODO These values were random and need to be verified from now to the end
    mat B(M-1,M-1,fill::eye);
    //correct
    mat W1 = B*Q*B.t();
    //Take inverse
    eigen_inv(W1);
    //correct
    mat u1_1 = G1.t()*W1*G1;
    eigen_inv(u1_1);
    //correct
    vec u1 = vectorise(u1_1*G1.t()*W1*h1);//Initial position guess. compute position deviations

    mat u_mid = zeros<mat>(G1.size(),G1.size());

    vec ri_hat = zeros<vec>(10);

    for (int j=0; j<RptCnt; j++) {
        // ri_hat = sqrt(sum((S-u1(1:end-1)*ones(1,M)).^2));
        //This is correct  
        ri_hat = vectorise(sqrt(sum(pow(S.t() - u1.rows(0,u1.n_rows-2)*ones(1,M),2))));

        B = 2*diagmat(ri_hat(span(1,M-1))); 
        
        W1 = B*Q*B.t(); // up date weights
        //take inverse
        eigen_inv(W1);

/*        cout << "W1_2 " << endl;
        W1.print();*/

        u_mid = G1.t()*W1*G1;
        eigen_inv(u_mid);

        //TODO I will hack a solution for now since the inv is not solving correctly
        u1 = (u_mid*G1.t()*W1*h1);
/*        cout << "u1: " << endl;
        u1.print();*/
     
    }
    //vector containing one zero
    dvec zero_vec = zeros<vec>(1);
    //correct
    dvec u1p = u1 - vectorise(join_rows(S.row(0), zero_vec));

    //========== second stage =====================================
    //correct
    dvec h2 = pow(u1p,2);
  
    dmat G_eye(u1p.size()-1,u1p.size()-1,fill::eye);
    //correct
    dmat G2 = join_cols(G_eye, ones(1,u1p.size()-1));
    //correct
    dmat B2 = 2 * diagmat(u1p);
 
    //W2 = inv(B2')*(G1'*W1*G1)*inv(B2);
    //correct equation but W1 is not right
    dmat W2 = inv(B2.t())*(G1.t()*W1*G1)*inv(B2);

    // u2 = inv(G2'*W2*G2)*G2'*W2*h2;
    //off a bit but correct
    dmat u2 = inv(G2.t()*W2*G2)*G2.t()*W2*h2;

    //=========== mapping ========================================
    //correct but u2 may be off
    mat SourceLocation = sign(diagmat(u1p(span(0,u2.n_rows-1))))*sqrt(abs(u2))+S.row(0).t();
    
    //correct
    if( as_scalar(u1(u1.size()-1)) < 0 || (as_scalar(min(u2)) < 0)) {
        SourceLocation = u1(span(0,u2.n_rows-1));
    }

    return SourceLocation;
}
int main()
{
    int L = 10000; // Number of ensemble runs

    arma_rng::set_seed_random();  // set the seed to a random value

    vec uo = {-50, 250}; //Source position

    vec x = {0, -5, 4, -2, 7, -7, 2, -4, 3, 1}; // Sensor position matrix
    vec y = {0, 8, 6, 4, 3, 5, 5, 2, 3, 8};
    mat S = join_rows(x,y);

    int M = 10;
    int N = 2;

    mat ones_mat = ones<mat>(1,M);
    //correct
    vec ro = vectorise(sqrt(sum(pow((uo*ones_mat - S.t()),2))));
  
    //correct
    vec rdo = ro(span(1,ro.size()-1)) - ro(0); 

    mat R(M-1,M-1,fill::eye); //Covariance 
    mat A(M-1,M-1,fill::ones);
    //correct
    R = (R + A)/2;
    //correct
    vec NsePwrVecdB = regspace<vec>(-60, 4, -24);

    //Initialize variables
    double nsePwr = 0;
    double simulationMSE = 0;
    mat Q = zeros<mat>(M-1,M-1);
    vec crlb = zeros<vec>(M);
    vec rdNse = zeros<vec>(M);
    vec rd = zeros<vec>(M-1);
    vec u = zeros<vec>(uo.size());
    vec mse = zeros<vec>(M);

    //-----
    //START LOOP HERE
    // ------
    cout << "Starting simulation" << endl;
    for(int i=0; i<NsePwrVecdB.size(); i++) {
        //correct
        nsePwr = pow(10,NsePwrVecdB(0)/10);
        //correct
        Q = nsePwr * R;//add tdoa correlated gaussian random noises with covariance matrix 

        //Sum of the elements on the main diagonal of matrix
        //lowest crlb contains most likely location
        crlb(i) = trace(TDOALocCRLB(S,uo,Q));

        //Monte Carlo Simulation
        simulationMSE = 0;
        // for(int j=0; j<L; j++) {

            //random noise
            rdNse = sqrt(nsePwr/2) * randn(M,1);
            // cout << " rdnse: " << endl;
            // rdNse.print();

            //add noise to 
            rd = rdo + rdNse(span(1,rdNse.size()-1)) - rdNse(1);
            // cout << " rd: " << endl;
            // rd.print();e
            //contains source location lat/long
            u = TDOALoc(S,rd,Q);
            simulationMSE =  + pow(norm(u-uo,2),2);
        // }   
        //TODO change this equation

        mse(i) = simulationMSE/L;
    }

/*figure(1); plot(NsePwrVecdB/2,10*log10(mse),'xk','MarkerSize',8); hold on;
plot(NsePwrVecdB/2,10*log10(crlb),'k'); grid on; hold off;

xlabel('10 log(c\sigma)'); 
ylabel('10 log(MSE)');
legend('New Method','CRLB');
ylim([0 60]);*/




    return 0;
}
