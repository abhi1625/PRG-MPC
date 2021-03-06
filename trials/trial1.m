BEGIN_ACADO;

    acadoSet('problemname','simplerocket');
    
    DifferentialState v;
    DifferentialState s;
    DifferentialState m;
    DifferentialState L;
    
    Control u;
    
    f = acado.DifferentialEquation();
    
    f.add(dot(s) == v);
    f.add(dot(v) == (u-0.02*v*v)/(m));
    f.add(dot(m) == -0.01*u*u);
    f.add(dot(L) == u*u);
    
    ocp = acado.OCP(0.0, 10.0, 20);
    
    ocp.minimizeMayerTerm(L);
    
    ocp.subjectTo(f);
    ocp.subjectTo('AT_START',s == 0.0);
    ocp.subjectTo('AT_START',v == 0.0);
    ocp.subjectTo('AT_START',L == 0.0);
    ocp.subjectTo('AT_START',m == 1.0);
    ocp.subjectTo('AT_END',s == 10.0);
    
    ocp.subjectTo('AT_END',v == 0.0);
    ocp.subjectTo(-0.01 <= v<= 1.3);
    
    algo = acado.OptimizationAlgorithm(ocp);
    algo.set('KKT.TOLERANCE', 1e-10);
    
END_ACADO;
    