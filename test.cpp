// trial 1 for MPC
#include<acado_toolkit.hpp>
#include<acado_gnuplot.hpp>

int main(){

  USING_NAMESPACE_ACADO;

  // variables

  DifferentialState       p_x;
  DifferentialState       p_y;
  DifferentialState       p_z;
  DifferentialState       q_w;
  DifferentialState       q_x;
  DifferentialState       q_y;
  DifferentialState       q_z;
  DifferentialState       v_x;
  DifferentialState       v_y;
  DifferentialState       v_z;

  Control                 T_w;
  Control                 w_x;
  Control                 w_y;
  Control                 w_z;

  DifferentialEquation      f;
  Function                h, hN;

  const double t_start = 0.0;
  const double t_end   = 1.0;
  const double dt      = 0.1;
  const int N          = round(t_end/dt);
  const double g_z     = 9.8066;
  const double w_max_yaw = 1;
  const double w_max_xy = 3;
  const double T_min = 2;
  const double T_max = 20;

  f << dot(p_x) == v_x;
  f << dot(p_y) == v_y;
  f << dot(p_z) == v_z;
  f << dot(q_w) == 0.5*(-w_x*q_x - w_y*q*y - w_z*q_z);
  f << dot(q_x) == 0.5*(w_x*q_w + w_z*q_y - w_y*q_z);
  f << dot(q_y) == 0.5*(w_y*q_w - w_z*q_x - w_x*q_z);
  f << dot(q_z) == 0.5*(w_z*q_w + w_y*q_x - w_z*q_y);
  f << dot(v_x) == 2*(q_w*q_y + q_x*q_z)*T;
  f << dot(v_y) == 2*(q_y*q_z - q_w*q_x)*T;
  f << dot(v_z) == (1 - 2*(q_x*q_x + q_y*q_y))*T - g_z;

  // Optimization Problem setup
  OCP ocp(t_start, t_end, N);

  h << p_x << p_y << p_z
    << q_w << q_x << q_y << q_z
    << v_x << v_y << v_z
    << T << w_x << w_y << w_z;

  hN << p_x << p_y << p_z
      << q_w << q_x << q_y << q_z
      << v_x << v_y << v_z

  DMatrix Q(h.getDim(), h.getDim());
  Q.setIdentity();
  Q(0,0) = 100;
  Q(1,1) = 100;
  Q(2,2) = 100;
  Q(3,3) = 100;
  Q(4,4) = 100;
  Q(5,5) = 100;
  Q(6,6) = 100;
  Q(7,7) = 10;
  Q(8,8) = 10;
  Q(9,9) = 10;
  Q(10,10) = 1;
  Q(11,11) = 1;
  Q(12,12) = 1;
  Q(13,13) = 1;

  DMatrix QN(hN.getDim(),hN.getDim());
  QN.setIdentity();
  QN(0,0) = Q(0,0);
  QN(1,1) = Q(1,1);
  QN(2,2) = Q(2,2);
  QN(3,3) = Q(3,3);
  QN(4,4) = Q(4,4);
  QN(5,5) = Q(5,5);
  QN(6,6) = Q(6,6);
  QN(7,7) = Q(7,7);
  QN(8,8) = Q(8,8);
  QN(9,9) = Q(9,9);

  //Reference for hover=2m qw=1

  DVector r(h.getDim());
  r.setZero();
  r(0) = 2.0;
  r(3) = 1.0;

  DVector rN(hN.getDim());
  rN.setZero();
  rN(0) = 2.0;
  rN(3) = 1.0;

  ocp.minimizeLSQ(Q,h,r);
  ocp.minimizeLSQEndTerm(QN, hN, rN);

  //System Dynamics

  ocp.subjectTo(f);

  //constraints
  ocp.subjectTo(-w_max_xy <= w_x <= w_max_xy);
  ocp.subjectTo(-w_max_xy <= w_y <= w_max_xy);
  ocp.subjectTo(-w_max_yaw <= w_z <= w_max_yaw);
  ocp.subjectTo(T_min <= T <= T_max);

  ocp.setNOD(10);

  // Set initial state
  ocp.subjectTo( AT_START, p_x ==  0.0 );
  ocp.subjectTo( AT_START, p_y ==  0.0 );
  ocp.subjectTo( AT_START, p_z ==  0.0 );
  ocp.subjectTo( AT_START, q_w ==  1.0 );
  ocp.subjectTo( AT_START, q_x ==  0.0 );
  ocp.subjectTo( AT_START, q_y ==  0.0 );
  ocp.subjectTo( AT_START, q_z ==  0.0 );
  ocp.subjectTo( AT_START, v_x ==  0.0 );
  ocp.subjectTo( AT_START, v_y ==  0.0 );
  ocp.subjectTo( AT_START, v_z ==  0.0 );
  ocp.subjectTo( AT_START, w_x ==  0.0 );
  ocp.subjectTo( AT_START, w_y ==  0.0 );
  ocp.subjectTo( AT_START, w_z ==  0.0 );

  /*  // Setup some visualization
  GnuplotWindow window1( PLOT_AT_EACH_ITERATION );
  window1.addSubplot( p_x,"position x" );
  window1.addSubplot( p_y,"position y" );
  window1.addSubplot( p_z,"position z" );
  window1.addSubplot( v_x,"verlocity x" );
  window1.addSubplot( v_y,"verlocity y" );
  window1.addSubplot( v_z,"verlocity z" );

  GnuplotWindow window3( PLOT_AT_EACH_ITERATION );
  window3.addSubplot( w_x,"rotation-acc x" );
  window3.addSubplot( w_y,"rotation-acc y" );
  window3.addSubplot( w_z,"rotation-acc z" );
  window3.addSubplot( T,"Thrust" );
*/

    // Define an algorithm to solve it.
  OptimizationAlgorithm algorithm(ocp);
  algorithm.set( INTEGRATOR_TOLERANCE, 1e-6 );
  algorithm.set( KKT_TOLERANCE, 1e-3 );
  //algorithm << window1;
  //algorithm << window3;
  algorithm.solve();




}
