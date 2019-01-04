/*
*    This file is part of ACADO Toolkit.
*
*    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
*    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
*    Developed within the Optimization in Engineering Center (OPTEC) under
*    supervision of Moritz Diehl. All rights reserved.
*
*    ACADO Toolkit is free software; you can redistribute it and/or
*    modify it under the terms of the GNU Lesser General Public
*    License as published by the Free Software Foundation; either
*    version 3 of the License, or (at your option) any later version.
*
*    ACADO Toolkit is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public
*    License along with ACADO Toolkit; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*
*/


/**
*    Author David Ariens, Rien Quirynen
*    Date 2009-2013
*    http://www.acadotoolkit.org/matlab 
*/

#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado/utils/matlab_acado_utils.hpp>

USING_NAMESPACE_ACADO

#include <mex.h>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
 { 
 
    MatlabConsoleStreamBuf mybuf;
    RedirectStream redirect(std::cout, mybuf);
    clearAllStaticCounters( ); 
 
    mexPrintf("\nACADO Toolkit for Matlab - Developed by David Ariens and Rien Quirynen, 2009-2013 \n"); 
    mexPrintf("Support available at http://www.acadotoolkit.org/matlab \n \n"); 

    if (nrhs != 4){ 
      mexErrMsgTxt("This problem expects 4 right hand side argument(s) since you have defined 4 MexInput(s)");
    } 
 
    TIME autotime;
    DifferentialState p_x;
    DifferentialState p_y;
    DifferentialState p_z;
    DifferentialState q_w;
    DifferentialState q_x;
    DifferentialState q_y;
    DifferentialState q_z;
    DifferentialState v_x;
    DifferentialState v_y;
    DifferentialState v_z;
    Control T;
    Control w_x;
    Control w_y;
    Control w_z;
    double *mexinput0_temp = NULL; 
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !(mxGetM(prhs[0])==1 && mxGetN(prhs[0])==1) ) { 
      mexErrMsgTxt("Input 0 must be a noncomplex scalar double.");
    } 
    mexinput0_temp = mxGetPr(prhs[0]); 
    double mexinput0 = *mexinput0_temp; 

    int mexinput1_count = 0;
    if (mxGetM(prhs[1]) == 1 && mxGetN(prhs[1]) >= 1) 
       mexinput1_count = mxGetN(prhs[1]);
    else if (mxGetM(prhs[1]) >= 1 && mxGetN(prhs[1]) == 1) 
       mexinput1_count = mxGetM(prhs[1]);
    else 
       mexErrMsgTxt("Input 1 must be a noncomplex double vector of dimension 1xY.");

    double *mexinput1_temp = NULL; 
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) { 
      mexErrMsgTxt("Input 1 must be a noncomplex double vector of dimension 1xY.");
    } 
    mexinput1_temp = mxGetPr(prhs[1]); 
    DVector mexinput1(mexinput1_count);
    for( int i=0; i<mexinput1_count; ++i ){ 
        mexinput1(i) = mexinput1_temp[i];
    } 

    int mexinput2_count = 0;
    if (mxGetM(prhs[2]) == 1 && mxGetN(prhs[2]) >= 1) 
       mexinput2_count = mxGetN(prhs[2]);
    else if (mxGetM(prhs[2]) >= 1 && mxGetN(prhs[2]) == 1) 
       mexinput2_count = mxGetM(prhs[2]);
    else 
       mexErrMsgTxt("Input 2 must be a noncomplex double vector of dimension 1xY.");

    double *mexinput2_temp = NULL; 
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) { 
      mexErrMsgTxt("Input 2 must be a noncomplex double vector of dimension 1xY.");
    } 
    mexinput2_temp = mxGetPr(prhs[2]); 
    DVector mexinput2(mexinput2_count);
    for( int i=0; i<mexinput2_count; ++i ){ 
        mexinput2(i) = mexinput2_temp[i];
    } 

    double *mexinput3_temp = NULL; 
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ) { 
      mexErrMsgTxt("Input 3 must be a noncomplex double vector of dimension XxY.");
    } 
    mexinput3_temp = mxGetPr(prhs[3]); 
    DMatrix mexinput3(mxGetM(prhs[3]), mxGetN(prhs[3]));
    for( int i=0; i<mxGetN(prhs[3]); ++i ){ 
        for( int j=0; j<mxGetM(prhs[3]); ++j ){ 
           mexinput3(j,i) = mexinput3_temp[i*mxGetM(prhs[3]) + j];
        } 
    } 

    Function acadodata_f2;
    acadodata_f2 << p_x;
    acadodata_f2 << p_y;
    acadodata_f2 << p_z;
    acadodata_f2 << q_w;
    acadodata_f2 << q_x;
    acadodata_f2 << q_y;
    acadodata_f2 << q_z;
    acadodata_f2 << v_x;
    acadodata_f2 << v_y;
    acadodata_f2 << v_z;
    acadodata_f2 << T;
    acadodata_f2 << w_x;
    acadodata_f2 << w_y;
    acadodata_f2 << w_z;
    DMatrix acadodata_M1;
    acadodata_M1.read( "quad_mpc_4_data_acadodata_M1.txt" );
    Function acadodata_f3;
    acadodata_f3 << p_x;
    acadodata_f3 << p_y;
    acadodata_f3 << p_z;
    acadodata_f3 << q_w;
    acadodata_f3 << q_x;
    acadodata_f3 << q_y;
    acadodata_f3 << q_z;
    acadodata_f3 << v_x;
    acadodata_f3 << v_y;
    acadodata_f3 << v_z;
    DMatrix acadodata_M2;
    acadodata_M2.read( "quad_mpc_4_data_acadodata_M2.txt" );
    DifferentialEquation acadodata_f1;
    acadodata_f1 << dot(p_x) == v_x;
    acadodata_f1 << dot(p_y) == v_y;
    acadodata_f1 << dot(p_z) == v_z;
    acadodata_f1 << dot(q_w) == ((-w_x)*q_x-q_y*w_y-q_z*w_z)*5.00000000000000000000e-01;
    acadodata_f1 << dot(q_x) == (q_w*w_x+q_y*w_z-q_z*w_y)*5.00000000000000000000e-01;
    acadodata_f1 << dot(q_y) == (q_w*w_y-q_x*w_z+q_z*w_x)*5.00000000000000000000e-01;
    acadodata_f1 << dot(q_z) == (q_w*w_z+q_x*w_y+q_y*w_z)*5.00000000000000000000e-01;
    acadodata_f1 << dot(v_x) == (q_w*q_y+q_x*q_z)*2.00000000000000000000e+00*T;
    acadodata_f1 << dot(v_y) == (-q_w*q_x+q_y*q_z)*2.00000000000000000000e+00*T;
    acadodata_f1 << dot(v_z) == ((1.00000000000000000000e+00-2.00000000000000000000e+00*q_x*q_x-2.00000000000000000000e+00*q_y*q_y)*T-9.80659999999999953957e+00);

    OCP ocp1(0, 0.2, 2);
    ocp1.minimizeLSQ(acadodata_M1, acadodata_f2, mexinput1);
    ocp1.minimizeLSQEndTerm(acadodata_M2, acadodata_f3, mexinput1);
    ocp1.subjectTo(acadodata_f1);
    ocp1.subjectTo((-3.00000000000000000000e+00) <= w_x <= 3.00000000000000000000e+00);
    ocp1.subjectTo((-3.00000000000000000000e+00) <= w_y <= 3.00000000000000000000e+00);
    ocp1.subjectTo((-1.00000000000000000000e+00) <= w_z <= 1.00000000000000000000e+00);
    ocp1.subjectTo(2.00000000000000000000e+00 <= T <= 2.00000000000000000000e+01);
    ocp1.setNOD( 10 );


    OutputFcn acadodata_f4;

    DynamicSystem dynamicsystem1( acadodata_f1,acadodata_f4 );
    Process process2( dynamicsystem1,INT_RK45 );

    RealTimeAlgorithm algo1(ocp1, 0.1);
    algo1.set( MAX_NUM_ITERATIONS, 3 );
    algo1.set( INTEGRATOR_TOLERANCE, 1.000000E-06 );
    algo1.set( KKT_TOLERANCE, 1.000000E-03 );
    algo1.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON );
    algo1.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING );

    StaticReferenceTrajectory referencetrajectory(mexinput3);
    Controller controller3( algo1,referencetrajectory );

    SimulationEnvironment algo2(0, 5, process2, controller3);
     algo2.init(mexinput2);
    returnValue returnvalue = algo2.run();


    VariablesGrid out_processout; 
    VariablesGrid out_feedbackcontrol; 
    VariablesGrid out_feedbackparameter; 
    VariablesGrid out_states; 
    VariablesGrid out_algstates; 
    algo2.getSampledProcessOutput(out_processout);
    algo2.getProcessDifferentialStates(out_states);
    algo2.getFeedbackControl(out_feedbackcontrol);
    const char* outputFieldNames[] = {"STATES_SAMPLED", "CONTROLS", "PARAMETERS", "STATES", "ALGEBRAICSTATES", "CONVERGENCE_ACHIEVED"}; 
    plhs[0] = mxCreateStructMatrix( 1,1,6,outputFieldNames ); 
    mxArray *OutSS = NULL;
    double  *outSS = NULL;
    OutSS = mxCreateDoubleMatrix( out_processout.getNumPoints(),1+out_processout.getNumValues(),mxREAL ); 
    outSS = mxGetPr( OutSS );
    for( int i=0; i<out_processout.getNumPoints(); ++i ){ 
      outSS[0*out_processout.getNumPoints() + i] = out_processout.getTime(i); 
      for( int j=0; j<out_processout.getNumValues(); ++j ){ 
        outSS[(1+j)*out_processout.getNumPoints() + i] = out_processout(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"STATES_SAMPLED",OutSS );
    mxArray *OutS = NULL;
    double  *outS = NULL;
    OutS = mxCreateDoubleMatrix( out_states.getNumPoints(),1+out_states.getNumValues(),mxREAL ); 
    outS = mxGetPr( OutS );
    for( int i=0; i<out_states.getNumPoints(); ++i ){ 
      outS[0*out_states.getNumPoints() + i] = out_states.getTime(i); 
      for( int j=0; j<out_states.getNumValues(); ++j ){ 
        outS[(1+j)*out_states.getNumPoints() + i] = out_states(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"STATES",OutS );
    mxArray *OutC = NULL;
    double  *outC = NULL;
    OutC = mxCreateDoubleMatrix( out_feedbackcontrol.getNumPoints(),1+out_feedbackcontrol.getNumValues(),mxREAL ); 
    outC = mxGetPr( OutC );
    for( int i=0; i<out_feedbackcontrol.getNumPoints(); ++i ){ 
      outC[0*out_feedbackcontrol.getNumPoints() + i] = out_feedbackcontrol.getTime(i); 
      for( int j=0; j<out_feedbackcontrol.getNumValues(); ++j ){ 
        outC[(1+j)*out_feedbackcontrol.getNumPoints() + i] = out_feedbackcontrol(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"CONTROLS",OutC );
    mxArray *OutP = NULL;
    double  *outP = NULL;
    OutP = mxCreateDoubleMatrix( out_feedbackparameter.getNumPoints(),1+out_feedbackparameter.getNumValues(),mxREAL ); 
    outP = mxGetPr( OutP );
    for( int i=0; i<out_feedbackparameter.getNumPoints(); ++i ){ 
      outP[0*out_feedbackparameter.getNumPoints() + i] = out_feedbackparameter.getTime(i); 
      for( int j=0; j<out_feedbackparameter.getNumValues(); ++j ){ 
        outP[(1+j)*out_feedbackparameter.getNumPoints() + i] = out_feedbackparameter(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"PARAMETERS",OutP );
    mxArray *OutZ = NULL;
    double  *outZ = NULL;
    OutZ = mxCreateDoubleMatrix( out_algstates.getNumPoints(),1+out_algstates.getNumValues(),mxREAL ); 
    outZ = mxGetPr( OutZ );
    for( int i=0; i<out_algstates.getNumPoints(); ++i ){ 
      outZ[0*out_algstates.getNumPoints() + i] = out_algstates.getTime(i); 
      for( int j=0; j<out_algstates.getNumValues(); ++j ){ 
        outZ[(1+j)*out_algstates.getNumPoints() + i] = out_algstates(i, j); 
       } 
    } 

    mxSetField( plhs[0],0,"ALGEBRAICSTATES",OutZ );
    mxArray *OutConv = NULL;
    if ( returnvalue == SUCCESSFUL_RETURN ) { OutConv = mxCreateDoubleScalar( 1 ); }else{ OutConv = mxCreateDoubleScalar( 0 ); } 
    mxSetField( plhs[0],0,"CONVERGENCE_ACHIEVED",OutConv );


    clearAllStaticCounters( ); 
 
} 

