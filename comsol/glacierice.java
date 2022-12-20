/*
 * glacierice.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Sep 25 2012, 11:29 by COMSOL 4.2.1.110. */
public class glacierice {

  public static void main(String[] args) {
    run();
  }

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/data/docs2/rheadley/src/svn/ESD_Simulator/comsol");

    model.modelNode().create("mod1");

    model.geom().create("geom1", 3);

    model.mesh().create("mesh1", "geom1");

    model.physics().create("spf", "LaminarFlow", "geom1");

    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");

    /*model.geom("geom1").run();*/

    model.view().create("view2", "geom1");
    model.view("view2").camera().set("preserveaspect", "off");

    model.param().set("rho_ice", "917[kg/m^3]");
    model.param().descr("rho_ice", "ice density");
    model.param().set("g", "9.8[m/s^2]");
    model.param().descr("g", "gravitational acceleration");
    model.param().set("spy", "365.25*24*60*60");
    model.param().descr("spy", "seconds per year");
    model.param().set("At", "6.8e-24");
    model.param()
         .descr("At", "from Paterson for temperate ice [s^(-1)*Pa^(-3)] (used by Brian)");
    model.param().set("n", "3");
    model.param().descr("n", "Glen's flow constant");
    model.param().set("epsilon", "1e-30");
    model.param()
         .descr("epsilon", "to make sure that eta doesn't go to zero");
    model.param().set("etalin", "3.46e12");
    model.param().descr("etalin", "viscosity for viscous approximation");
    model.param().set("nsl", "3");
    model.param().descr("nsl", "sliding exponent (used by brian)");
    model.param().set("Asl", "1e-15/(365*24*60*60)[m*Pa^-3/s]");
    model.param().descr("Asl", "sliding constant (from ice cascade for now)");
    model.param().set("mb", "1");
    model.param()
         .descr("mb", "test mass balance; eventually make this a function that gets input into the program");
    model.param().set("delt", ".01[yr]");
    model.param().descr("delt", "time step (years)");
    model.param().set("Tzb", "0[degC]");
    model.param().descr("Tzb", "test Tzb");
    model.param().set("xmin", "1.0e5[m]");
    model.param().descr("xmin", "");
    model.param().set("xmax", "1.5e5[m]");
    model.param().descr("xmax", "");
    model.param().set("ymin", "1.1e5[m]");
    model.param().descr("ymin", "");
    model.param().set("ymax", "1.6e5[m]");
    model.param().descr("ymax", "");

    model.variable().create("var1");
    model.variable("var1").model("mod1");
    model.variable("var1")
         .set("eta", "1/2*(At^-(1/n))*(epseff)^((1-n)/(n))*1[Pa*s]");
    model.variable("var1").descr("eta", "effective viscosity");
    model.variable("var1").set("epsxx", "ux");
    model.variable("var1").descr("epsxx", "");
    model.variable("var1").set("epsxy", "1/2*(uy+vx)");
    model.variable("var1").descr("epsxy", "");
    model.variable("var1").set("epsxz", "1/2*uz");
    model.variable("var1").descr("epsxz", "");
    model.variable("var1").set("epsyy", "vy");
    model.variable("var1").descr("epsyy", "");
    model.variable("var1").set("epsyz", "1/2*vz");
    model.variable("var1").descr("epsyz", "");
    model.variable("var1").set("epszz", "wz");
    model.variable("var1").descr("epszz", "");
    model.variable("var1")
         .set("epseff", "(epsxx^2+epsyy^2+epsxx*epsyy+epsxy^2+epsxz^2+epsyz^2+epsilon^2)^(1/2)*1[s]");
    model.variable("var1").descr("epseff", "effective strain rate");
    model.variable("var1").set("etax", "dtang(eta,x)");
    model.variable("var1").descr("etax", "");
    model.variable("var1").set("etay", "dtang(eta,y)");
    model.variable("var1").descr("etay", "");
    model.variable("var1").set("epsxzf", "1/2*(uz+wx)");
    model.variable("var1").descr("epsxzf", "");
    model.variable("var1").set("epsyzf", "1/2*(vz+wy)");
    model.variable("var1").descr("epsyzf", "");
    model.variable("var1")
         .set("fullepseff", "(epsxx^2+epsyy^2+epsxx*epsyy+epsxy^2+epsxzf^2+epsyzf^2+epsilon^2)^(1/2)*1[s]");
    model.variable("var1").descr("fullepseff", "effective strain rate");
    model.variable("var1")
         .set("fulleta", "1/2*(At^-(1/n))*(fullepseff)^((1-n)/(n))*1[Pa*s]");
    model.variable("var1").descr("fulleta", "");
    model.variable("var1").set("streff", "(fullepseff/At)^(1/n)");
    model.variable("var1")
         .descr("streff", "effective stress/second invarient (make sure this is right?!)");
    model.variable().create("var2");
    model.variable("var2").model("mod1");
    model.variable("var2").set("S", "surface(z)");
    model.variable("var2").descr("S", "surface elevation");
    model.variable("var2").set("B", "bed(z)");
    model.variable("var2").descr("B", "bed elevation");
    model.variable("var2").set("H", "S-B");
    model.variable("var2").descr("H", "ice thickness");
    model.variable("var2").set("dzdx", "dtang(z,x)");
    model.variable("var2").descr("dzdx", "");
    model.variable("var2").set("dzdy", "dtang(z,y)");
    model.variable("var2").descr("dzdy", "");
    model.variable("var2").name("Geometry");
    model.variable().create("var3");
    model.variable("var3").model("mod1");
    model.variable("var3").set("tauxz", "2*fulleta*.5*(uz+wx)");
    model.variable("var3").descr("tauxz", "");
    model.variable("var3").set("tauyz", "2*fulleta*.5*(vz+wy)");
    model.variable("var3").descr("tauyz", "");
    model.variable("var3").set("uslx", "bed(((Asl/0.8)*tauxz^nsl))");
    model.variable("var3").descr("uslx", "");
    model.variable("var3").set("usly", "bed(((Asl/0.8)*tauyz^nsl))");
    model.variable("var3").descr("usly", "");
    model.variable("var3").set("taubxSIA", "rho_ice*g*H*dSdx");
    model.variable("var3").descr("taubxSIA", "");
    model.variable("var3").set("taubySIA", "rho_ice*g*H*dSdy");
    model.variable("var3").descr("taubySIA", "");
    model.variable("var3").set("usl", "ifsl*sqrt(uslx^2+usly^2)");
    model.variable("var3").descr("usl", "");
    model.variable("var3").name("Sliding");
    model.variable("var3").set("Tb", "basalT(x,y)");
    model.variable("var3").set("ifsl", "(Tb>=Tpmp)*1+0");
    model.variable("var3").descr("ifsl", "test if greater than freezing");
    model.variable().create("var4");
    model.variable("var4").model("mod1");
    model.variable("var4").name("Outputs");
    model.variable("var4").set("MB", "mb(x,y)/spy");
    model.variable("var4").descr("MB", "");
    model.variable("var4").set("Hp", "p/(rho_ice*g)");
    model.variable("var4")
         .descr("Hp", "ice thickness (full thickness at bed)");
    model.variable("var4").set("fluxx", "genproj1(utot)");
    model.variable("var4").descr("fluxx", "Flux using general extrusion");
    model.variable("var4").set("fluxy", "genproj1(vtot)");
    model.variable("var4").descr("fluxy", "");
    model.variable("var4").set("dqdx", "dtang(fluxx,x)");
    model.variable("var4").descr("dqdx", "x-derivative of flux2");
    model.variable("var4").set("dqdy", "dtang(fluxy,y)");
    model.variable("var4").descr("dqdy", "y-derivative of flux2");
    model.variable("var4").set("uave", "genproj1(u)/Hp");
    model.variable("var4").descr("uave", "depth-averaged velocity");
    model.variable("var4").set("p1x", "dtang(uave,x)*Hp");
    model.variable("var4").descr("p1x", "test flux part 1");
    model.variable("var4").set("p2x", "dtang(Hp,x)*uave");
    model.variable("var4").descr("p2x", "test flux part 2");
    model.variable("var4").set("dqdx2", "p1x+p2x");
    model.variable("var4").descr("dqdx2", "compare flux components");
    model.variable("var4").set("Tpmpold", "T0-beta*p");
    model.variable("var4").set("Tpmp", "T0-8.7e-4*H*1[K/m]");
    model.variable("var4").descr("Tpmp", "same as in ICE.f90");    
    model.variable("var4")
         .descr("Tpmp", "calculate pressure-melting point at base of glacier");
    model.variable("var4").set("T0", "0[degC]");
    model.variable("var4").descr("T0", "");
    model.variable("var4").set("beta", "9.8e-8[K/Pa]");
    model.variable("var4")
         .descr("beta", "dependence of melting on pressure, Pattyn, 2003, table 1");
    model.variable("var4").set("delh", "delt*(MB-(dqdx+dqdy))");
    model.variable("var4")
         .descr("delh", "change in surface elevation; based on Bryan's. It's different because my flux terms include sliding implicitly");
    model.variable("var4").set("utot", "u+uslx*ifsl");
    model.variable("var4").set("vtot", "u+usly*ifsl");

    model.func().create("int1", "Interpolation");
    model.func("int1").set("source", "file");
    model.func("int1")
         .set("filename", "/data/docs2/rheadley/src/svn/ESD_Simulator/output/comsol/export_to_comsol_tempb.txt");
    model.func("int1").setIndex("funcs", "basalT", 0, 0);
    model.func("int1").set("argunit", "m");
    model.func("int1").set("fununit", "degC");
    model.func().create("int2", "Interpolation");
    model.func().create("int3", "Interpolation");
    model.func().create("int4", "Interpolation");
    model.func("int2").set("source", "file");
    model.func("int2")
         .set("filename", "/data/docs2/rheadley/src/svn/ESD_Simulator/output/comsol/export_to_comsol_ice_mass_balance.txt");
    model.func("int2").setIndex("funcs", "mb", 0, 0);
    model.func("int3").set("source", "file");
    model.func("int3")
         .set("filename", "/data/docs2/rheadley/src/svn/ESD_Simulator/output/comsol/export_to_comsol_ithick.txt");
    model.func("int3").setIndex("funcs", "H", 0, 0);
    model.func("int4").set("source", "file");
    model.func("int4")
         .set("filename", "/data/docs2/rheadley/src/svn/ESD_Simulator/output/comsol/export_to_comsol_bed.txt");
    model.func("int4").setIndex("funcs", "zb", 0, 0);
    model.func().create("int5", "Interpolation");
    model.func("int5").set("source", "file");
    model.func("int5")
         .set("filename", "/data/docs2/rheadley/src/svn/ESD_Simulator/output/comsol/export_to_comsol_isurf.txt");
    model.func("int5").setIndex("funcs", "isurf", 0, 0);

    model.geom("geom1").feature().create("ps1", "ParametricSurface");
    model.geom("geom1").feature("ps1").set("parmin1", "xmin");
    model.geom("geom1").feature("ps1").set("parmax1", "xmax");
    model.geom("geom1").feature("ps1").set("parmin2", "ymin");
    model.geom("geom1").feature("ps1").set("parmax2", "ymax");
    model.geom("geom1").feature("ps1").setIndex("coord", "s1", 0);
    model.geom("geom1").feature("ps1").setIndex("coord", "s2", 1);
    model.geom("geom1").feature("ps1")
         .setIndex("coord", "isurf(s1,s2)-H(s1,s2)", 2);
    model.geom("geom1").feature("ps1").set("rtol", "1.0E-4");
    model.geom("geom1").feature("ps1").set("maxknots", "10");
    model.geom("geom1").run("ps1");
    model.geom("geom1").feature().duplicate("ps2", "ps1");
    model.geom("geom1").feature("ps2").setIndex("coord", "isurf(s1,s2)", 2);
    model.geom("geom1").run("ps2");
    model.geom("geom1").run("ps2");
    model.geom("geom1").feature().create("blk1", "Block");
    model.geom("geom1").feature("blk1").setIndex("size", "(xmax-xmin)", 0);
    model.geom("geom1").feature("blk1").setIndex("size", "(ymax-ymin)", 1);
    model.geom("geom1").feature("blk1").setIndex("size", "4000", 2);
    model.geom("geom1").feature("blk1").setIndex("pos", "xmin", 0);
    model.geom("geom1").feature("blk1").setIndex("pos", "xmax", 1);
    model.geom("geom1").feature("blk1").setIndex("pos", "-1500", 2);
    model.geom("geom1").run("blk1");
    model.geom("geom1").feature().create("uni1", "Union");
    model.geom("geom1").run("blk1");
    model.geom("geom1").runPre("uni1");
    model.geom("geom1").feature("blk1").setIndex("pos", "ymin", 1);
    model.geom("geom1").run("blk1");
    model.geom("geom1").feature("uni1").selection("input")
         .set(new String[]{"blk1", "ps1", "ps2"});
    model.geom("geom1").feature("uni1").set("repairtol", "1.0E-4");
    model.geom("geom1").run("blk1");
    model.geom("geom1").feature().create("del1", "Delete");
    model.geom("geom1").feature().remove("del1");
    model.geom("geom1").run("uni1");
    model.geom("geom1").run("uni1");
    model.geom("geom1").feature().create("del1", "Delete");
    model.geom("geom1").feature("del1").selection("input").init(3);
    model.geom("geom1").feature("del1").selection("input")
         .set("uni1", new int[]{1, 3});
    model.geom("geom1").run("del1");
    model.geom("geom1").run();

    model.cpl().create("genext1", "GeneralExtrusion", "geom1");
    model.cpl("genext1").selection().geom("geom1", 2);
    model.cpl("genext1").selection().set(new int[]{4});
    model.cpl("genext1").set("opname", "surface");
    model.cpl("genext1").setIndex("dstmap", "", 2);
    model.cpl("genext1").set("usesrcmap", "on");
    model.cpl("genext1").setIndex("srcmap", "", 2);
    model.cpl().duplicate("genext2", "genext1");
    model.cpl("genext2").selection().set(new int[]{3});
    model.cpl("genext2").set("opname", "bed");
    model.cpl().create("genproj1", "GeneralProjection", "geom1");
    model.cpl("genproj1").selection().set(new int[]{1});

    model.selection().create("sel1", "Explicit");
    model.selection("sel1").geom(2);
    model.selection("sel1").set(new int[]{3});
    model.selection().duplicate("sel2", "sel1");
    model.selection().duplicate("sel3", "sel1");
    model.selection("sel1").name("Bed");
    model.selection("sel2").set(new int[]{4});
    model.selection("sel2").name("Surface");
    model.selection("sel3").set(new int[]{1, 2, 5, 6});
    model.selection("sel3").name("Edges");

    model.mesh("mesh1").feature().create("bl1", "BndLayer");
    model.mesh("mesh1").feature("bl1").feature()
         .create("blp", "BndLayerProp");
    model.mesh("mesh1").feature().create("ftet1", "FreeTet");
    model.mesh("mesh1").feature().create("conv1", "Convert");
    model.mesh("mesh1").feature("size").set("hauto", "5");	
    /* In this, size 6 is coarse, 5 is normal, 4 is fine, etc. */
    model.mesh("mesh1").feature("bl1").feature("blp").selection()
         .named("sel1");
    model.mesh("mesh1").run();

    model.physics("spf").feature("fp1").set("rho_mat", 1, "userdef");
    model.physics("spf").feature("fp1").set("mu_mat", 1, "userdef");
    model.physics("spf").prop("CompressibilityProperty")
         .set("Compressibility", 1, "Incompressible");
    model.physics("spf").feature("fp1").set("rho", 1, "rho_ice");
    model.physics("spf").feature("fp1").set("mu", 1, "fulleta");
    model.physics("spf").feature("wall1").name("Bottom");
    model.physics("spf").feature().create("vf1", "VolumeForce", 3);
    model.physics("spf").feature().create("open1", "OpenBoundary", 2);
    model.physics("spf").feature().create("open2", "OpenBoundary", 2);
    model.physics("spf").feature("init1")
         .set("u", new String[]{"10/spy", "10/spy", "0"});
    model.physics("spf").feature("init1").set("p", 1, "1e5");
    model.physics("spf").feature("vf1").selection().set(new int[]{1});
    model.physics("spf").feature("vf1")
         .set("F", new String[]{"0", "0", "-rho_ice*g"});
    model.physics("spf").feature("open1").selection().named("sel2");
    model.physics("spf").feature("open1").name("Top Free Surface");
    model.physics("spf").feature("open2")
         .set("BoundaryCondition", 1, "NoViscousStress");
    model.physics("spf").feature("open2").selection().named("sel3");
    model.physics("spf").feature("open2").name("Sidewalls Open");
    model.physics("spf").feature().create("wall2", "Wall", 2);
    model.physics("spf").feature("wall2")
         .set("BoundaryCondition", 1, "MovingWall");
    model.physics("spf").feature("wall2").active(false);

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").feature().create("st1", "StudyStep");
    model.sol("sol1").feature("st1").set("study", "std1");
    model.sol("sol1").feature("st1").set("studystep", "stat");
    model.sol("sol1").feature().create("v1", "Variables");
    model.sol("sol1").feature().create("s1", "Stationary");
    model.sol("sol1").feature("s1").feature().create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("dtech", "auto");
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 25);
    model.sol("sol1").feature("s1").feature().create("i1", "Iterative");
    model.sol("sol1").feature("s1").feature("i1").set("linsolver", "gmres");
    model.sol("sol1").feature("s1").feature("i1").set("prefuntype", "left");
    model.sol("sol1").feature("s1").feature("i1").set("itrestart", 50);
    model.sol("sol1").feature("s1").feature("i1").set("rhob", 20);
    model.sol("sol1").feature("s1").feature("i1").set("maxlinit", 200);
    model.sol("sol1").feature("s1").feature("i1").set("nlinnormuse", "on");
    model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "i1");
    model.sol("sol1").feature("s1").feature("i1").feature()
         .create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1")
         .set("prefun", "gmg");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1")
         .set("mcasegen", "any");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature().create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("iter", 2);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("linerelax", 0.3);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("linealgorithm", "matrix");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("linemethodmatrix", "coupled");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("linevar", "mod1_u");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("maxline", 15);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("seconditer", 1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr")
         .feature("sl1").set("relax", 0.3);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature().create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("iter", 2);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("linerelax", 0.4);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("linealgorithm", "matrix");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("linemethodmatrix", "coupled");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("linevar", "mod1_u");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("maxline", 15);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("seconditer", 2);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po")
         .feature("sl1").set("relax", 0.5);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs")
         .feature().create("d1", "Direct");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs")
         .feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("dtech", "auto");
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 25);
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    model.sol("sol1").attach("std1");

    model.result().create("pg1", 3);
    model.result("pg1").set("data", "dset1");
    model.result("pg1").feature().create("slc1", "Slice");
    model.result("pg1").feature("slc1").set("expr", new String[]{"spf.U"});
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").name("Velocity (spf)");
    model.result().dataset().create("surf1", "Surface");
    model.result().dataset("surf1").set("data", "dset1");
    model.result().create("pg2", 3);
    model.result("pg2").set("data", "surf1");
    model.result("pg2").set("frametype", "spatial");
    model.result("pg2").feature().create("surf1", "Surface");
    model.result("pg2").feature("surf1").set("expr", new String[]{"1"});
    model.result("pg2").feature("surf1").set("coloring", "uniform");
    model.result("pg2").feature("surf1").set("color", "gray");
    model.result("pg2").feature().create("con1", "Contour");
    model.result("pg2").feature("con1").set("expr", new String[]{"p"});
    model.result("pg2").feature("con1").set("number", 40);
    model.result("pg2").name("Pressure (spf)");

    model.sol("sol1").feature("s1").feature("dDef").active(true);
    model.sol("sol1").runAll();

    model.result("pg1").run();
    model.result("pg1").run();
    model.result("pg1").feature().remove("slc1");
    model.result("pg1").run();
    model.result("pg1").feature().create("surf1", "Surface");
    model.result("pg1").feature("surf1").set("expr", "spf.U*spy");
    model.result("pg1").run();
    model.result().export().create("data1", "Data");
    model.result().export("data1").setIndex("expr", "delh", 0);
    model.result().export("data1").setIndex("unit", "m", 0);
    model.result().export("data1").setIndex("expr", "usl*spy", 1);
    model.result().export("data1")
         .set("filename", "/data/docs2/rheadley/src/svn/ESD_Simulator/output/comsol/import_from_comsol.txt");
    model.result().dataset().create("surf2", "Surface");
    model.result().dataset("surf2").selection().named("sel1");
    model.result().export("data1").set("data", "surf2");
    model.result().export("data1").set("header", false); 
    model.result().export("data1").set("resolution", "fine");
    model.result().export("data1").run();

    return model;
  }

}