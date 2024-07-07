
#include "MMA.hpp"

// Eigen lib
#include <Eigen/Core>
#include <Eigen/Dense>


void MMA::mmasub(
    int const m,
    int const n,
    int const iter,
    Eigen::VectorXd const& xval,
    Eigen::VectorXd const& xmin,
    Eigen::VectorXd const& xmax,
    Eigen::VectorXd const& df0dx,
    Eigen::VectorXd const& fval,
    Eigen::MatrixXd const& dfdx,
    Eigen::VectorXd& xold1,
    Eigen::VectorXd& xold2,
    Eigen::VectorXd& low,
    Eigen::VectorXd& upp,
    Eigen::VectorXd& ox) {
  // Initialization of MMA parameters
  double const epsimin= 1.e-7;
  double const raa0= 1.e-5;
  double const move= 0.5;
  double const albefa= 0.1;

  double const asyinit= 0.5;
  double const asyincr= 1.2;
  double const asydecr= 0.7;

  // Computation of asymptotes low and upp
  if (iter < 2) {
    low= xval - asyinit * (xmax - xmin);
    upp= xval + asyinit * (xmax - xmin);
  }
  else {
    Eigen::VectorXd zzz= (xval - xold1).cwiseProduct(xold1 - xold2);
    Eigen::VectorXd factor= Eigen::VectorXd::Ones(n);
    for (int i= 0; i < n; i++) {
      if (zzz(i) > 0.0) factor(i)= asyincr;
      if (zzz(i) < 0.0) factor(i)= asydecr;
    }
    low= xval - factor.cwiseProduct(xold1 - low);
    low= low.cwiseMax(xval - 10.0 * (xmax - xmin));
    low= low.cwiseMin(xval - 0.01 * (xmax - xmin));
    upp= xval + factor.cwiseProduct(upp - xold1);
    upp= upp.cwiseMin(xval + 10.0 * (xmax - xmin));
    upp= upp.cwiseMax(xval + 0.01 * (xmax - xmin));
  }

  // Computation of bounds alpha and beta:
  Eigen::VectorXd alpha= ((low + albefa * (xval - low)).cwiseMax(xval - move * (xmax - xmin))).cwiseMax(xmin);
  Eigen::VectorXd beta= ((upp - albefa * (upp - xval)).cwiseMin(xval + move * (xmax - xmin))).cwiseMin(xmax);

  // Precomputations
  Eigen::VectorXd xmamiinv= ((xmax - xmin).cwiseMax(1.e-5 * Eigen::VectorXd::Ones(n))).cwiseInverse();
  Eigen::VectorXd ux2= (upp - xval).cwiseAbs2();
  Eigen::VectorXd xl2= (xval - low).cwiseAbs2();
  Eigen::VectorXd uxinv= (upp - xval).cwiseInverse();
  Eigen::VectorXd xlinv= (xval - low).cwiseInverse();

  // Computation of p0 and q0
  Eigen::VectorXd p0= (df0dx).cwiseMax(Eigen::VectorXd::Zero(n));
  Eigen::VectorXd q0= (-df0dx).cwiseMax(Eigen::VectorXd::Zero(n));
  Eigen::VectorXd pq0= 0.001 * (p0 + q0) + raa0 * xmamiinv;
  p0= (p0 + pq0).cwiseProduct(ux2);
  q0= (q0 + pq0).cwiseProduct(xl2);

  // Computation of P and Q
  Eigen::MatrixXd P= (dfdx).cwiseMax(Eigen::MatrixXd::Zero(m, n));
  Eigen::MatrixXd Q= (-dfdx).cwiseMax(Eigen::MatrixXd::Zero(m, n));
  Eigen::MatrixXd PQ= 0.001 * (P + Q) + raa0 * Eigen::VectorXd::Ones(m) * xmamiinv.transpose();
  P= P + PQ;
  Q= Q + PQ;
  P= P * ux2.asDiagonal();
  Q= Q * xl2.asDiagonal();

  // Computation of a0, a, b, c, d
  double a0= 1.0;
  Eigen::VectorXd a= Eigen::MatrixXd::Zero(int(fval.size()), 1);
  Eigen::VectorXd b= P * uxinv + Q * xlinv - fval;
  Eigen::VectorXd c= 1000.0 * Eigen::MatrixXd::Ones(int(fval.size()), 1);
  Eigen::VectorXd d= Eigen::MatrixXd::Zero(int(fval.size()), 1);

  // Solve the subproblem by a primal dual Newton method
  MMA::subsolve_barrier_method(m, n, epsimin, low, upp, alpha, beta, p0, q0, P, Q, a0, a, b, c, d, ox);

  // Update the variable history
  xold2= xold1;
  xold1= xval;
}


void MMA::subsolve_barrier_method(
    int const m,
    int const n,
    double const epsimin,
    Eigen::VectorXd const& low,
    Eigen::VectorXd const& upp,
    Eigen::VectorXd const& alpha,
    Eigen::VectorXd const& beta,
    Eigen::VectorXd const& p0,
    Eigen::VectorXd const& q0,
    Eigen::MatrixXd const& P,
    Eigen::MatrixXd const& Q,
    double const a0,
    Eigen::VectorXd const& a,
    Eigen::VectorXd const& b,
    Eigen::VectorXd const& c,
    Eigen::VectorXd const& d,
    Eigen::VectorXd& ox) {
  // Initialize W
  Eigen::VectorXd x= 0.5 * (alpha + beta);
  Eigen::VectorXd y= Eigen::VectorXd::Ones(m);
  double z= 1.0;
  Eigen::VectorXd lam= Eigen::VectorXd::Ones(m);
  Eigen::VectorXd xsi= ((x - alpha).cwiseInverse()).cwiseMax(Eigen::VectorXd::Ones(n));
  Eigen::VectorXd eta= ((beta - x).cwiseInverse()).cwiseMax(Eigen::VectorXd::Ones(n));
  Eigen::VectorXd mu= (0.5 * c).cwiseMax(Eigen::VectorXd::Ones(m));
  double zet= 1.0;
  Eigen::VectorXd s= Eigen::VectorXd::Ones(m);

  // Update W with the barrier method
  for (double epsilon= 1.0; epsilon > epsimin; epsilon*= 0.1)
    MMA::barrier_method_iteration(m, n, low, upp, alpha, beta, p0, q0, P, Q, a0, a, b, c, d, epsilon, x, y, z, lam, xsi, eta, mu, zet, s);

  // Output the first block-component of W : X
  ox= x;
}


void MMA::barrier_method_iteration(
    int const m,
    int const n,
    Eigen::VectorXd const& low,
    Eigen::VectorXd const& upp,
    Eigen::VectorXd const& alpha,
    Eigen::VectorXd const& beta,
    Eigen::VectorXd const& p0,
    Eigen::VectorXd const& q0,
    Eigen::MatrixXd const& P,
    Eigen::MatrixXd const& Q,
    double const a0,
    Eigen::VectorXd const& a,
    Eigen::VectorXd const& b,
    Eigen::VectorXd const& c,
    Eigen::VectorXd const& d,
    double const epsilon,
    Eigen::VectorXd& x,
    Eigen::VectorXd& y,
    double& z,
    Eigen::VectorXd& lam,
    Eigen::VectorXd& xsi,
    Eigen::VectorXd& eta,
    Eigen::VectorXd& mu,
    double& zet,
    Eigen::VectorXd& s) {
  // Assemble F
  Eigen::VectorXd F(3 * n + 4 * m + 2);
  MMA::assemble_F_vector(m, n, low, upp, alpha, beta, p0, q0, P, Q, a0, a, b, c, d, epsilon, x, y, z, lam, xsi, eta, mu, zet, s, F);
  double Fnorm= F.norm();
  double Fmax= F.cwiseAbs().maxCoeff();

  // Solve for W in F*W= 0 with Newton method
  for (int iter= 0; iter < 200; iter++) {
    if (Fmax <= 0.9 * epsilon) break;

    // Compute delta_W
    Eigen::VectorXd Dx, Dy, Dlam, Dxsi, Deta, Dmu, Ds;
    double Dz, Dzet;
    MMA::newton_method_find_direction(m, n, x, y, z, lam, xsi, eta, mu, zet, s, low, upp, alpha, beta, p0, q0, P, Q, a0, a, b, c, d, epsilon,
                                                     Dx, Dy, Dz, Dlam, Dxsi, Deta, Dmu, Dzet, Ds);

    // Find the best step and update W and therefore Fmax and Fnorm
    MMA::newton_method_find_coefficient(m, n, low, upp, alpha, beta, p0, q0, P, Q, a0, a, b, c, d, epsilon, Dx, Dy, Dz, Dlam, Dxsi, Deta, Dmu, Dzet, Ds,
                                                       x, y, z, lam, xsi, eta, mu, zet, s, Fmax, Fnorm);
  }
}


void MMA::newton_method_find_direction(
    int const m,
    int const n,
    Eigen::VectorXd const& x,
    Eigen::VectorXd const& y,
    double const z,
    Eigen::VectorXd const& lam,
    Eigen::VectorXd const& xsi,
    Eigen::VectorXd const& eta,
    Eigen::VectorXd const& mu,
    double const zet,
    Eigen::VectorXd const& s,
    Eigen::VectorXd const& low,
    Eigen::VectorXd const& upp,
    Eigen::VectorXd const& alpha,
    Eigen::VectorXd const& beta,
    Eigen::VectorXd const& p0,
    Eigen::VectorXd const& q0,
    Eigen::MatrixXd const& P,
    Eigen::MatrixXd const& Q,
    double const a0,
    Eigen::VectorXd const& a,
    Eigen::VectorXd const& b,
    Eigen::VectorXd const& c,
    Eigen::VectorXd const& d,
    double const epsilon,
    Eigen::VectorXd& Dx,
    Eigen::VectorXd& Dy,
    double& Dz,
    Eigen::VectorXd& Dlam,
    Eigen::VectorXd& Dxsi,
    Eigen::VectorXd& Deta,
    Eigen::VectorXd& Dmu,
    double& Dzet,
    Eigen::VectorXd& Ds) {
  // To get delta_w, we need to solve J*delta_w = -F

  // Build necessary components to simplify the expression J*delta_w = -F
  Eigen::VectorXd ux2, xl2, plam, qlam, gvec, dpsidx;
  MMA::assemble_precomp(low, upp, p0, q0, P, Q, x, lam, ux2, xl2, plam, qlam, gvec, dpsidx);
  Eigen::MatrixXd GG= P * (ux2.cwiseInverse().asDiagonal()) - Q * (xl2.cwiseInverse().asDiagonal());

  Eigen::VectorXd epsvecn= epsilon * Eigen::VectorXd::Ones(n);
  Eigen::VectorXd epsvecm= epsilon * Eigen::VectorXd::Ones(m);

  Eigen::VectorXd delx= dpsidx - epsvecn.cwiseQuotient(x - alpha) + epsvecn.cwiseQuotient(beta - x);
  Eigen::VectorXd dely= c + d.cwiseProduct(y) - lam - epsvecm.cwiseQuotient(y);
  double delz= a0 - a.transpose() * lam - epsilon / z;
  Eigen::VectorXd dellam= gvec - a * z - y - b + epsvecm.cwiseQuotient(lam);
  Eigen::VectorXd diagx= plam.cwiseQuotient((upp - x).cwiseProduct(ux2));
  diagx= diagx + qlam.cwiseQuotient((x - low).cwiseProduct(xl2));
  diagx= 2.0 * diagx + xsi.cwiseQuotient(x - alpha) + eta.cwiseQuotient(beta - x);
  Eigen::VectorXd diagy= d + mu.cwiseQuotient(y);
  Eigen::VectorXd diaglamyi= s.cwiseQuotient(lam) + diagy.cwiseInverse();

  // In the m < n case, the simplified expression of J*delta_w = -F gives
  Eigen::VectorXd bb(m + 1);
  bb << (dellam + dely.cwiseQuotient(diagy) - GG * (delx.cwiseQuotient(diagx))), delz;
  Eigen::MatrixXd AA_1(m, m + 1);
  AA_1 << Eigen::MatrixXd(diaglamyi.asDiagonal()) + GG * (diagx.cwiseInverse().asDiagonal()) * (GG.transpose()), a;
  Eigen::MatrixXd AA_2(1, m + 1);
  AA_2 << a.transpose(), (-zet / z);
  Eigen::MatrixXd A(m + 1, m + 1);
  A << AA_1, AA_2;

  // Equation to solve
  Eigen::MatrixXd solution= A.fullPivLu().solve(bb);

  // Retrieve delta_w from the simplified expression
  Dlam= solution.topRows(m);
  Dz= solution(m, 0);
  Dx= -delx.cwiseQuotient(diagx) - (GG.transpose() * Dlam).cwiseQuotient(diagx);
  Dy= -dely.cwiseQuotient(diagy) + Dlam.cwiseQuotient(diagy);
  Dxsi= -xsi + epsvecn.cwiseQuotient(x - alpha) - (xsi.cwiseProduct(Dx)).cwiseQuotient(x - alpha);
  Deta= -eta + epsvecn.cwiseQuotient(beta - x) + (eta.cwiseProduct(Dx)).cwiseQuotient(beta - x);
  Dmu= -mu + epsvecm.cwiseQuotient(y) - (mu.cwiseProduct(Dy)).cwiseQuotient(y);
  Dzet= -zet + epsilon / z - zet * Dz / z;
  Ds= -s + epsvecm.cwiseQuotient(lam) - (s.cwiseProduct(Dlam)).cwiseQuotient(lam);
}


void MMA::newton_method_find_coefficient(
    int const m,
    int const n,
    Eigen::VectorXd const& low,
    Eigen::VectorXd const& upp,
    Eigen::VectorXd const& alpha,
    Eigen::VectorXd const& beta,
    Eigen::VectorXd const& p0,
    Eigen::VectorXd const& q0,
    Eigen::MatrixXd const& P,
    Eigen::MatrixXd const& Q,
    double const a0,
    Eigen::VectorXd const& a,
    Eigen::VectorXd const& b,
    Eigen::VectorXd const& c,
    Eigen::VectorXd const& d,
    double const epsilon,
    Eigen::VectorXd const& Dx,
    Eigen::VectorXd const& Dy,
    double const Dz,
    Eigen::VectorXd const& Dlam,
    Eigen::VectorXd const& Dxsi,
    Eigen::VectorXd const& Deta,
    Eigen::VectorXd const& Dmu,
    double const Dzet,
    Eigen::VectorXd const& Ds,
    Eigen::VectorXd& x,
    Eigen::VectorXd& y,
    double& z,
    Eigen::VectorXd& lam,
    Eigen::VectorXd& xsi,
    Eigen::VectorXd& eta,
    Eigen::VectorXd& mu,
    double& zet,
    Eigen::VectorXd& s,
    double& Fmax,
    double& Fnorm) {
  // Find the largets coefficient step t such that conditions 5.9j - 5.9.n are satisfied
  Eigen::VectorXd W_prim(2 * n + 4 * m + 2);
  W_prim << y, z, lam, xsi, eta, mu, zet, s;
  Eigen::VectorXd dW_prim(2 * n + 4 * m + 2);
  dW_prim << Dy, Dz, Dlam, Dxsi, Deta, Dmu, Dzet, Ds;

  double t_temp_w= (-1.01 * dW_prim.cwiseQuotient(W_prim)).maxCoeff();
  double t_temp_alpha= (-1.01 * Dx.cwiseQuotient(x - alpha)).maxCoeff();
  double t_temp_beta= (1.01 * Dx.cwiseQuotient(beta - x)).maxCoeff();
  double t= std::max(std::max(std::max(t_temp_w, t_temp_alpha), t_temp_beta), 1.0);
  t= 1.0 / t;

  // Save the initial value of W for later
  Eigen::VectorXd x_old= x;
  Eigen::VectorXd y_old= y;
  double z_old= z;
  Eigen::VectorXd lam_old= lam;
  Eigen::VectorXd xsi_old= xsi;
  Eigen::VectorXd eta_old= eta;
  Eigen::VectorXd mu_old= mu;
  double zet_old= zet;
  Eigen::VectorXd s_old= s;

  // Needed variables to update Fmax
  Eigen::VectorXd F(3 * n + 4 * m + 2);
  double F_old= Fnorm;
  Fnorm= 2.0 * F_old;

  // We search t such that F(w_new) < F(w_old) using an iterative bisection algorithm
  for (int iter= 0; iter < 50; iter++) {
    if (Fnorm <= F_old) break;

    // Update W
    x= x_old + t * Dx;
    y= y_old + t * Dy;
    z= z_old + t * Dz;
    lam= lam_old + t * Dlam;
    xsi= xsi_old + t * Dxsi;
    eta= eta_old + t * Deta;
    mu= mu_old + t * Dmu;
    zet= zet_old + t * Dzet;
    s= s_old + t * Ds;

    MMA::assemble_F_vector(m, n, low, upp, alpha, beta, p0, q0, P, Q, a0, a, b, c, d, epsilon, x, y, z, lam, xsi, eta, mu, zet, s, F);

    Fnorm= F.norm();

    // Change t
    t= t / 2.0;
  }
  Fmax= F.cwiseAbs().maxCoeff();
}


void MMA::assemble_F_vector(
    int const m,
    int const n,
    Eigen::VectorXd const& low,
    Eigen::VectorXd const& upp,
    Eigen::VectorXd const& alpha,
    Eigen::VectorXd const& beta,
    Eigen::VectorXd const& p0,
    Eigen::VectorXd const& q0,
    Eigen::MatrixXd const& P,
    Eigen::MatrixXd const& Q,
    double const a0,
    Eigen::VectorXd const& a,
    Eigen::VectorXd const& b,
    Eigen::VectorXd const& c,
    Eigen::VectorXd const& d,
    double const epsilon,
    Eigen::VectorXd const& x,
    Eigen::VectorXd const& y,
    double const z,
    Eigen::VectorXd const& lam,
    Eigen::VectorXd const& xsi,
    Eigen::VectorXd const& eta,
    Eigen::VectorXd const& mu,
    double const zet,
    Eigen::VectorXd const& s,
    Eigen::VectorXd& F) {
  // Precomputation
  Eigen::VectorXd ux2, xl2, plam, qlam, gvec, dpsidx;
  MMA::assemble_precomp(low, upp, p0, q0, P, Q, x, lam, ux2, xl2, plam, qlam, gvec, dpsidx);

  Eigen::VectorXd epsvecn= epsilon * Eigen::VectorXd::Ones(n);
  Eigen::VectorXd epsvecm= epsilon * Eigen::VectorXd::Ones(m);

  // Compute coefficients
  Eigen::VectorXd residx= dpsidx - xsi + eta;
  Eigen::VectorXd residy= c + d.cwiseProduct(y) - mu - lam;
  double residz= a0 - zet - a.transpose().dot(lam);
  Eigen::VectorXd residlam= gvec - a * z - y + s - b;
  Eigen::VectorXd residxsi= xsi.cwiseProduct(x - alpha) - epsvecn;
  Eigen::VectorXd resideta= eta.cwiseProduct(beta - x) - epsvecn;
  Eigen::VectorXd residmu= mu.cwiseProduct(y) - epsvecm;
  double residzet= zet * z - epsilon;
  Eigen::VectorXd resids= lam.cwiseProduct(s) - epsvecm;

  // Assemble vector
  F << residx, residy, residz, residlam, residxsi, resideta, residmu, residzet, resids;
}


void MMA::assemble_precomp(
    Eigen::VectorXd const& low,
    Eigen::VectorXd const& upp,
    Eigen::VectorXd const& p0,
    Eigen::VectorXd const& q0,
    Eigen::MatrixXd const& P,
    Eigen::MatrixXd const& Q,
    Eigen::VectorXd const& x,
    Eigen::VectorXd const& lam,
    Eigen::VectorXd& ux2,
    Eigen::VectorXd& xl2,
    Eigen::VectorXd& plam,
    Eigen::VectorXd& qlam,
    Eigen::VectorXd& gvec,
    Eigen::VectorXd& dpsidx) {
  ux2= (upp - x).cwiseAbs2();
  xl2= (x - low).cwiseAbs2();
  plam= p0 + P.transpose() * lam;
  qlam= q0 + Q.transpose() * lam;
  gvec= P * ((upp - x).cwiseInverse()) + Q * ((x - low).cwiseInverse());
  dpsidx= plam.cwiseQuotient(ux2) - qlam.cwiseQuotient(xl2);
}
