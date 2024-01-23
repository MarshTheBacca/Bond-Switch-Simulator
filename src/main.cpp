#include <iostream>
#include <sstream>
#include "outputfile.h"
#include "linked_network.h"
#include "lammps_c_interface.h"
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

int main(){

    //Set up logfile
    Logfile logfile("./netmc.log");
    logfile.datetime("Simulation begun at: ");
    logfile.write("NETwork Monte Carlo");
    logfile.write("Written By: David OM, Wilson Group, 2019");
    logfile.separator();

    //Read input file
    logfile.write("Reading input file");
    ++logfile.currIndent;
    ifstream inputFile("./netmc.inpt", ios::in);
    if(!inputFile.good()) logfile.criticalError("Cannot find input file 'netmc.inpt' in current directory");
    //I/O
    string skip,line,prefixFolderOut, prefixFileOut; //prefix for output files
    string prefixFolderIn, prefixFileIn;
    bool boolIn, boolRestart;
    getline(inputFile,skip);

    getline(inputFile,line);
    istringstream(line)>>prefixFolderOut;

    getline(inputFile,line);
    istringstream(line)>>prefixFileOut;

    getline(inputFile,line);
    istringstream(line)>>boolIn;
    getline(inputFile,line);
    istringstream(line)>>boolRestart;

    getline(inputFile,line);
    istringstream(line)>>prefixFolderIn;
    getline(inputFile,line);
    istringstream(line)>>prefixFileIn;



    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("I/O read");
    //Network properties
    int nRings,minRingSize,maxRingSize;
    string lattice;
    string crystal;
    string moveType;
    int minCnd,maxCnd;
    bool fixed_rings;
    string fixed_rings_file;
    getline(inputFile,line);
    istringstream(line)>>nRings;
    getline(inputFile,line);
    istringstream(line)>>minRingSize;
    getline(inputFile,line);
    istringstream(line)>>maxRingSize;
    getline(inputFile,line);
    istringstream(line)>>minCnd;
    getline(inputFile,line);
    istringstream(line)>>maxCnd;
    getline(inputFile,line);
    istringstream(line)>>lattice;
    getline(inputFile,line);
    istringstream(line)>>crystal;
    getline(inputFile,line);
    istringstream(line)>>fixed_rings;
    getline(inputFile,line);
    istringstream(line)>>fixed_rings_file;

    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Network properties read");

    bool isOpenMP, isSimpleGraphene, isTriangleRaft, isBilayer, isTersoffGraphene, isBN;
    int MC_Routine;
    getline(inputFile, line);
    istringstream(line)>>isOpenMP;
    getline(inputFile, line);
    istringstream(line)>>isSimpleGraphene;
    getline(inputFile, line);
    istringstream(line)>>isTriangleRaft;
    getline(inputFile, line);
    istringstream(line)>>isBilayer;
    getline(inputFile, line);
    istringstream(line)>>isTersoffGraphene;
    getline(inputFile, line);
    istringstream(line)>>isBN;
    getline(inputFile, line);
    istringstream(line)>>MC_Routine;
    getline(inputFile,skip);
    getline(inputFile,skip);



    if (MC_Routine==1 && !isSimpleGraphene) exit(2);
    if (MC_Routine==2 && !isTriangleRaft) exit(2);
    if (MC_Routine==3 && !isBilayer) exit(2);
    if (MC_Routine==4 && !isTersoffGraphene) exit(2);
    if (MC_Routine==5 && !isBN) exit(2);

//    cout << isSimpleGraphene << isTersoffGraphene << isTriangleRaft << isBilayer << endl;


    //Monte carlo
    string runType;
    int randomSeed;
    getline(inputFile,line);
    cout << line << endl;
    istringstream(line)>>runType;
    cout << "Run Type : " << runType << endl;
    getline(inputFile,line);
    istringstream(line)>>moveType;
    getline(inputFile,line);
    istringstream(line)>>randomSeed;
    bool spiral;
    int spiralRadius;
    string MCWeighting;
    getline(inputFile,line);
    istringstream(line)>>spiral;
    getline(inputFile,line);
    istringstream(line)>>spiralRadius;
    getline(inputFile,line);
    istringstream(line)>>MCWeighting;

    cout << "Spiral : " << spiral << endl;
    cout << "Spiral Radius : " << spiralRadius << endl;
    cout << "MC Weighting : " << MCWeighting << endl;

    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Monte Carlo parameters read");
    //Energy search
    int mcSteps,equilSteps;
    double mcStartT,mcEndT,mcIncT,mcThermT;
    getline(inputFile,line);
    istringstream(line)>>mcStartT;
    getline(inputFile,line);
    istringstream(line)>>mcEndT;
    getline(inputFile,line);
    istringstream(line)>>mcIncT;
    getline(inputFile,line);
    istringstream(line)>>mcThermT;
    getline(inputFile,line);
    istringstream(line)>>mcSteps;
    getline(inputFile,line);
    istringstream(line)>>equilSteps;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Energy search parameters read");
    //Cost search
    int costEquilSteps,costInitSteps,costSteps;
    double costLowLimP,costUpLimP,costIncP,costLowLimR,costUpLimR,costIncR,costT;
    getline(inputFile,line);
    istringstream(line)>>costEquilSteps;
    getline(inputFile,line);
    istringstream(line)>>costLowLimP;
    getline(inputFile,line);
    istringstream(line)>>costUpLimP;
    getline(inputFile,line);
    istringstream(line)>>costIncP;
    getline(inputFile,line);
    istringstream(line)>>costLowLimR;
    getline(inputFile,line);
    istringstream(line)>>costUpLimR;
    getline(inputFile,line);
    istringstream(line)>>costIncR;
    getline(inputFile,line);
    istringstream(line)>>costInitSteps;
    getline(inputFile,line);
    istringstream(line)>>costSteps;
    getline(inputFile,line);
    istringstream(line)>>costT;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Cost search parameters read");
    //Cost function
    double costPK,costRK;
    getline(inputFile,line);
    istringstream(line)>>costPK;
    getline(inputFile,line);
    istringstream(line)>>costRK;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Cost function parameters read");
    //Potential model
    int potConvex;
    double potAK,potBK,potCK;
    getline(inputFile,line);
    istringstream(line)>>potBK;
    getline(inputFile,line);
    istringstream(line)>>potAK;
    getline(inputFile,line);
    istringstream(line)>>potCK;
    getline(inputFile,line);
    istringstream(line)>>potConvex;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Potential parameters read");
    //Geometry optimisation
    int goptLocalIt, goptGlobalIt, goptLocalSize;
    double goptTau, goptTol;
    getline(inputFile,line);
    istringstream(line)>>goptLocalIt;
    getline(inputFile,line);
    istringstream(line)>>goptGlobalIt;
    getline(inputFile,line);
    istringstream(line)>>goptTau;
    getline(inputFile,line);
    istringstream(line)>>goptTol;
    getline(inputFile,line);
    istringstream(line)>>goptLocalSize;
    getline(inputFile,skip);
    getline(inputFile,skip);
    logfile.write("Geometry optimisation parameters read");
    //Analysis
    int analysisFreq,writeStructures,structureFreq;
    getline(inputFile,line);
    istringstream(line)>>analysisFreq;
    getline(inputFile,line);
    istringstream(line)>>writeStructures;
    getline(inputFile,line);
    istringstream(line)>>structureFreq;
    logfile.write("Write parameters read");

    float lj_cutoff;
    getline(inputFile,skip);
    getline(inputFile,skip);
    getline(inputFile,line);
    istringstream(line)>>lj_cutoff;
    getline(inputFile,line);


    inputFile.close();
    --logfile.currIndent;
    logfile.write("Input file closed");
    logfile.separator();

    cout << "LAMMPS BOND ANGLE :: " << potAK << " " << potBK << endl;

    OutputFile lammpsBondAngle(  prefixFolderIn+"/PARM_Si.lammps");
    lammpsBondAngle.write("bond_style harmonic");
    lammpsBondAngle.write("bond_coeff 1 "+ to_string(potBK)+" 1.000");
    lammpsBondAngle.write("angle_style cosine/squared");
    lammpsBondAngle.write("angle_coeff 1 "+ to_string(potAK)+" 120");


    //Initialise network
    logfile.write("Initialising network");
    ++logfile.currIndent;
    logfile.write("Lattice type:", lattice);
    logfile.write("Number of rings:",nRings);
    logfile.write("Min ring size:",minRingSize);
    logfile.write("Max ring size:",maxRingSize);
    logfile.write("Min node coordination:",minCnd);
    logfile.write("Max node coordination:",maxCnd);
    logfile.write("Monte Carlo run type:",runType);
    logfile.write("Monte Carlo move type:",moveType);
    if(runType=="energy") {
        logfile.write("Monte carlo initial log temperature:", mcStartT);
        logfile.write("Monte carlo final log temperature:", mcEndT);
        logfile.write("Monte carlo log temperature increment:", mcIncT);
        logfile.write("Monte carlo steps per increment:", mcSteps);
        logfile.write("Equilibrium steps:", equilSteps);
    }
    else if(runType=="cost"){
        logfile.write("Monte carlo p lower limit:", costLowLimP);
        logfile.write("Monte carlo p upper limit:", costUpLimP);
        logfile.write("Monte carlo r lower limit:", costLowLimR);
        logfile.write("Monte carlo r upper limit:", costLowLimR);
        logfile.write("Monte carlo effective temperature:", costT);
        logfile.write("Monte carlo steps per increment:", costSteps);
    }
    bool mixedLattice=false; //whether mixed coordination lattice
    if(lattice.substr(0,3)=="mix" || lattice.substr(0,5)=="cairo" || lattice=="alt_square" || moveType=="mix") mixedLattice=true;

//    LinkedNetwork network(nRings, lattice, minCnd, maxCnd, minRingSize, maxRingSize);

//    if (!boolIn) {
    LinkedNetwork network(prefixFolderIn, prefixFileIn, prefixFolderOut,minCnd, maxCnd, minRingSize, maxRingSize, isSimpleGraphene, isTriangleRaft, isBilayer, isTersoffGraphene, isBN, boolRestart);
//    }


    minCnd = network.minACnxs;
    maxCnd = network.maxACnxs;
    minRingSize = network.minBCnxs;
    maxRingSize = network.maxBCnxs;

    //    if(boolIn){


//    network.write("./output_files/pre");
//    network.writeXYZ("./output_files/pre");

    // Find Ignored Rings
    network.pushPrefix(prefixFolderIn, prefixFolderOut);
    network.findIgnoredRings(fixed_rings, prefixFolderIn+"/"+fixed_rings_file);

    network.initialisePotentialModel(network.networkA, potAK,potBK,potCK,potConvex);
    cout << "Initialise Geometry Opt" << endl;
    network.initialiseGeometryOpt(goptLocalIt,goptTau,goptTol,goptLocalSize);
    cout << "Initialsie Monte Carlo" << endl;
    network.initialiseMonteCarlo(network.networkA, pow(10,mcStartT),randomSeed,mixedLattice);
    cout << "Initialise Cost Function" << endl;
    network.initialiseCostFunction(costT,randomSeed,costPK,costRK);
    network.isOpenMP=isOpenMP;
    network.MCWeighting=MCWeighting;
    cout << "305 " << MCWeighting << endl;
    network.isSimpleGraphene=isSimpleGraphene;
    network.isTersoffGraphene=isTersoffGraphene;
    network.isTriangleRaft=isTriangleRaft;
    network.isBilayer=isBilayer;
    network.isBN=isBN;
    cout << "309" << endl;
    network.spiralRadius = spiralRadius;
    if (spiral) network.makerFixed();
    cout << "310" << endl;
    network.MC_Routine=MC_Routine;
    if(crystal!="default") network.makeCrystal(crystal,lattice);
    if(lattice=="goldberg" || lattice=="inv_cubic") network.optimalProjection("sphere");



    if(mixedLattice){//get statistics on number of 3/4 coordinate nodes for mixed lattice
        cout << "Mixed Lattice " << endl;
        double mixA=network.networkT.nodes.n,mixA3=0,mixA4=0;
        for(int i=0; i<mixA; ++i){
            if(network.networkT.nodes[i].netCnxs.n==3) ++mixA3;
            else if(network.networkT.nodes[i].netCnxs.n==4) ++mixA4;
        }
        VecF<double> ringStats = network.getNodeDistribution("B");
        double mixAv = 0.0;
        for(int i=0; i<=maxRingSize; ++i) mixAv += i*ringStats[i];
        logfile.write("Mixed lattice total rings:",network.networkB.nodes.n);
        logfile.write("Mixed average ring size:",mixAv);
        logfile.write("Mixed lattice total nodes:",mixA);
        logfile.write("Mixed lattice 3 coordinate nodes:",mixA3/mixA);
        logfile.write("Mixed lattice 4 coordinate nodes:",mixA4/mixA);
    }
    --logfile.currIndent;
    logfile.write("Network initialised");
    logfile.separator();
    cout << "Network initialised " << endl;

    //Initialise output files
    logfile.write("Initialising analysis output files");
    ++logfile.currIndent;
    int status;
    string Folder;
    cout << "Initialising analysis output files" << endl;

//    Folder = prefixFolderOut+"_ti_"+to_string(int(100*mcStartT))+"_tf_"+ to_string(int(100*mcEndT));
//    status = mkdir(Folder.c_str(), 0777);
    //    if (status==0) cout << "Folder Generation Issue" << endl;

    Folder = prefixFolderOut;
    status = mkdir(Folder.c_str(), 0777);
    //status = mkdir(Folder.c_str());


    string prefixOut = Folder+"/"+prefixFileOut;
    OutputFile outEnergyStats(  prefixOut+"_e_compare.out");
    OutputFile outRingStats(    prefixOut+"_ringstats.out");
    OutputFile outCorr(         prefixOut+"_correlations.out");
    OutputFile outEnergy(       prefixOut+"_energy.out");
    OutputFile outEntropy(      prefixOut+"_entropy.out");
    OutputFile outTemperature(  prefixOut+"_temperature.out");
    OutputFile outGeometry(     prefixOut+"_geometry.out");
    OutputFile outEmatrix(      prefixOut+"_ematrix.out");
    OutputFile outGeomHist(     prefixOut+"_geomhist.out");
    OutputFile outAreas(        prefixOut+"_areas.out");
    OutputFile outClusterA(     prefixOut+"_cluster_a.out");
    OutputFile outClusterB(     prefixOut+"_cluster_b.out");
    OutputFile outCndStats(     prefixOut+"_cndstats.out");
    outGeometry.initVariables(6,4,60,20);
    outAreas.initVariables(6,4,60,30);
    outEmatrix.initVariables(1,4,60,int(log10(nRings*12))+2);
    outGeomHist.initVariables(6,4,60,20);
    outClusterB.initVariables(1,4,60,10);
    logfile.write("Ring statistics file created");
    logfile.write("Correlations file created");
    logfile.write("Energy file created");
    logfile.write("Entropy file created");
    logfile.write("Temperature file created");
    logfile.write("Geometry file created");
    logfile.write("Geometry histogram file created");
    logfile.write("Edge distribution file created");
    logfile.write("Coordination statistics file created");
    --logfile.currIndent;
    logfile.write("Files initialised");
    logfile.separator();

    bool lammps = true;

    //Initialise total analysis variables - only update in main simulation not equilibration
    VecF<double> lenHist(10000),angHist(10000);
    lenHist=0.0;
    angHist=0.0;

    //Run monte carlo
    cout<<"                             " << network.mc.getEnergy()<<endl;
    int accepted=0,optIterations=0;
    VecF<int> optCodes(5);
    optCodes=0;
    int trackFreq=100;
    VecF<int> moveStatus;
    if(runType=="energy") {//energy run
        //Run monte carlo thermalisation
        logfile.write("Running Monte Carlo thermalisation");
        cout << "Running Monte Carlo Thermalisation" << endl;
        ++logfile.currIndent;
        double energy=network.mc.getEnergy();
        cout << "Initial Energy : " << energy << endl;

        double SimpleGrapheneEnergy=0.0;
        double TersoffGrapheneEnergy=0.0;
        double TriangleRaftEnergy=0.0;
        double BilayerEnergy = 0.0;
        double BNEnergy = 0.0;
//        cout << isSimpleGraphene << isTersoffGraphene << isTriangleRaft << isBilayer << endl;

//        logfile.write("Simple Graphene : "+isSimpleGraphene);
//        logfile.write("Tersoff Graphene : "+isTersoffGraphene);
//        logfile.write("Triangle Raft : "+isTriangleRaft);
//        logfile.write("Bilayer : "+isBilayer);



        if (isSimpleGraphene)   {
//            cout << "Simple Graphene : " << isSimpleGraphene << " " << network.SimpleGraphene.GlobalPotentialEnergy() << endl;
            SimpleGrapheneEnergy=network.SimpleGraphene.GlobalPotentialEnergy();
        }
        if (isTersoffGraphene) {
//            cout << "Tersoff Graphene : " << isTersoffGraphene << endl;
            TersoffGrapheneEnergy = network.TersoffGraphene.GlobalPotentialEnergy();
        }
        if (isTriangleRaft)     {
//            cout << "Triangle Raft    : " << isTriangleRaft << " " << network.Triangle_Raft.GlobalPotentialEnergy() << endl;
            TriangleRaftEnergy=network.Triangle_Raft.GlobalPotentialEnergy();
        }
        if (isBilayer)          {
//            cout << "Bilayer           : " << isBilayer << endl;
            BilayerEnergy = network.Bilayer.GlobalPotentialEnergy();
        }
        if (isBN)               {
            BNEnergy = network.BN.GlobalPotentialEnergy();
        }

        logfile.write("Calculated energies for running systems");
        double mcT = pow(10, mcThermT);
        network.mc.setTemperature(mcT);
        if (spiral){
            cout << "Pre-run : attempting MC steps on all atoms within radius" << endl;
            bool disallowed_node;
            for (int i=0;i<network.networkA.nodes.n;++i){
                if (network.rFixed[i]>1000/spiralRadius){
                    disallowed_node=false;
                    for (int j=0;j<network.networkB.nodes[network.Fixed_Ring[0]].dualCnxs.n;++j){
                        if (i==network.networkB.nodes[network.Fixed_Ring[0]].dualCnxs[j]) disallowed_node=true;
                    }
                    if (!disallowed_node){
                        cout << "Node " << i << " withing radius at distance " << network.rFixed[i] << endl;
                        moveStatus = network.SpiralmonteCarloSwitchMoveLAMMPS(i, SimpleGrapheneEnergy, TersoffGrapheneEnergy, TriangleRaftEnergy, BilayerEnergy, BNEnergy,0);

                        double dt = logfile.timeElapsed();
                        string track =
                                to_string(accepted) + "/" + to_string(i) + " moves accepted/completed in " + to_string(dt) +
                                " seconds";
                        logfile.write(track);
                        cout << "t" << " " << i << " " << accepted << "  Energy = " << energy << endl;

                        VecF<double> ringStats = network.getNodeDistribution("B");
                        double r = network.getAssortativity("B");
                        double aEst = network.getAboavWeaireEstimate("B");
                        VecF<double> aw = network.getAboavWeaire("B");
                        VecF<double> s = network.getEntropy("B");
                        VecF<double> corr(6);
                        VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                        double rr = network.getAssortativity("A");
                        corr[0] = r;
                        corr[1] = aEst;
                        corr[2] = aw[0];
                        corr[3] = aw[1];
                        corr[4] = aw[2];
                        corr[5] = rr;


                        VecF<double> eCompare(6);

                        eCompare[0] = ringStats[6];
                        eCompare[1] = SimpleGrapheneEnergy;
                        eCompare[2] = TersoffGrapheneEnergy;
                        eCompare[3] = TriangleRaftEnergy;
                        eCompare[4] = BilayerEnergy;
                        eCompare[5] = BNEnergy;


                        VecF<double> emptyL,emptyA; //dummy histograms
                        VecF<double> geomStats = network.getOptimisationGeometry(network.networkA, emptyL,emptyA);
                        VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
//                VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
                        VecF<double> cndStats = network.getNodeDistribution("A");
//                network.getRingAreas(a,aSq);
                        outRingStats.writeRowVector(ringStats);
                        outCorr.writeRowVector(corr);
                        outEnergy.write(energy);
                        outEntropy.writeRowVector(s);
                        outTemperature.write(mcT);
                        outGeometry.writeRowVector(geomStats);
                        outAreas.writeRowVector(a);
                        outAreas.writeRowVector(aSq);

                        outEnergyStats.writeRowVector(eCompare);


                        for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
//                outClusterB.writeRowVector(clusters);
                        outCndStats.writeRowVector(cndStats);
                        //sleep(100);
                        if(mixedLattice) {
                            VecF<double> cluster(3);
                            cluster[0] = network.getMaxCluster("A", 3);
                            cluster[1] = network.getMaxCluster("A", 4);
                            cluster[2] = network.getAssortativity("A");
                            outClusterA.writeRowVector(cluster);
                        }



                    }
                }
            }
        }


        for (int i = 1; i <= equilSteps; ++i) {
            cout << endl << endl << endl << endl;
            cout << "                                               Equilibriating : " << i << " -- energy : " << SimpleGrapheneEnergy << " " << TersoffGrapheneEnergy << " " << TriangleRaftEnergy << " " << BilayerEnergy << " " << BNEnergy << endl;
            if(lammps){
                moveStatus = network.monteCarloSwitchMoveLAMMPS(SimpleGrapheneEnergy, TersoffGrapheneEnergy, TriangleRaftEnergy, BilayerEnergy, BNEnergy, 0);
            }
            else if(!mixedLattice) {
                moveStatus = network.monteCarloSwitchMove(network.networkT,energy);
                cout << moveStatus[0] << endl;
            }

            else moveStatus = network.monteCarloMixMove(energy);


            accepted += moveStatus[0];
            optCodes[moveStatus[1]] += 1;
            optIterations += moveStatus[2];
//        cout << i << endl;
            if (i % trackFreq == 0) {
                double dt = logfile.timeElapsed();
                string track =
                        to_string(accepted) + "/" + to_string(i) + " moves accepted/completed in " + to_string(dt) +
                        " seconds";
                logfile.write(track);
                cout << "t" << " " << i << " " << accepted << "  Energy = " << energy << endl;
            }
            if (i % analysisFreq == 0) {
                VecF<double> ringStats = network.getNodeDistribution("B");
                double r = network.getAssortativity("B");
                double aEst = network.getAboavWeaireEstimate("B");
                VecF<double> aw = network.getAboavWeaire("B");
                VecF<double> s = network.getEntropy("B");
                VecF<double> corr(6);
                VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                double rr = network.getAssortativity("A");
                corr[0] = r;
                corr[1] = aEst;
                corr[2] = aw[0];
                corr[3] = aw[1];
                corr[4] = aw[2];
                corr[5] = rr;


                VecF<double> eCompare(6);

                eCompare[0] = ringStats[6];
                eCompare[1] = SimpleGrapheneEnergy;
                eCompare[2] = TersoffGrapheneEnergy;
                eCompare[3] = TriangleRaftEnergy;
                eCompare[4] = BilayerEnergy;
		eCompare[5] = BNEnergy;


                VecF<double> emptyL,emptyA; //dummy histograms
                VecF<double> geomStats = network.getOptimisationGeometry(network.networkA, emptyL,emptyA);
                VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
//                VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
                VecF<double> cndStats = network.getNodeDistribution("A");
//                network.getRingAreas(a,aSq);
                outRingStats.writeRowVector(ringStats);
                outCorr.writeRowVector(corr);
                outEnergy.write(energy);
                outEntropy.writeRowVector(s);
                outTemperature.write(mcT);
                outGeometry.writeRowVector(geomStats);
                outAreas.writeRowVector(a);
                outAreas.writeRowVector(aSq);

                outEnergyStats.writeRowVector(eCompare);


                for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
//                outClusterB.writeRowVector(clusters);
                outCndStats.writeRowVector(cndStats);
                //sleep(100);
                if(mixedLattice) {
                    VecF<double> cluster(3);
                    cluster[0] = network.getMaxCluster("A", 3);
                    cluster[1] = network.getMaxCluster("A", 4);
                    cluster[2] = network.getAssortativity("A");
                    outClusterA.writeRowVector(cluster);
                }
            }
            if (i%structureFreq==0 && writeStructures==1) {
                network.syncCoordinates();
                //network.write(prefixOut+"_therm_"+to_string(i));
                cout << prefixOut+"_therm_"+to_string(i) << endl;
                network.writeXYZ(prefixOut+"_therm_"+to_string(i));
            }

        }
        --logfile.currIndent;
        logfile.write("Monte Carlo equilibration complete");
        logfile.separator();

        //Perform monte carlo simulation
        logfile.write("Running Monte Carlo simulation");
        ++logfile.currIndent;
        int nT = (mcEndT - mcStartT) / mcIncT;

        for (int t = 0; t <= nT; ++t) {

            mcT = pow(10, mcStartT + t * mcIncT);
            network.mc.setTemperature(mcT);
            cout << "Temperature : " << mcT << endl;
            logfile.write("Temperature:", mcT);
            ++logfile.currIndent;
            for (int i = 1; i <= mcSteps; ++i) {
                cout << "                                                       Running : " << i << " energy : " << SimpleGrapheneEnergy << " " << TersoffGrapheneEnergy << " " << TriangleRaftEnergy << " " << BilayerEnergy << " " << BNEnergy << endl;
                if(lammps){
                    moveStatus = network.monteCarloSwitchMoveLAMMPS(SimpleGrapheneEnergy, TersoffGrapheneEnergy, TriangleRaftEnergy, BilayerEnergy, BNEnergy, 0);
                }
                else if(!mixedLattice)   moveStatus = network.monteCarloSwitchMove(network.networkT, energy);
                else                moveStatus = network.monteCarloMixMove(energy);
//                    cout << "Mixed Lattice" << endl;

                accepted += moveStatus[0];
                optCodes[moveStatus[1]] += 1;
                optIterations += moveStatus[2];
//            cout << i << endl;
                if (i % trackFreq == 0) {
                    double dt = logfile.timeElapsed();
                    string track =
                            to_string(accepted) + "/" + to_string(i) + " moves accepted/completed in " + to_string(dt) +
                            " seconds";
                    logfile.write(track);
                    cout << t << " " << i << " " << accepted << "  Energy = " << energy << endl;
                }
                if (i % analysisFreq == 0) {

                    VecF<double> ringStats = network.getNodeDistribution("B");
                    double r = network.getAssortativity("B");
                    double aEst = network.getAboavWeaireEstimate("B");
                    VecF<double> aw = network.getAboavWeaire("B");
                    VecF<double> s = network.getEntropy("B");
                    VecF<double> corr(6);
                    VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                    double rr = network.getAssortativity("A");
                    corr[0] = r;
                    corr[1] = aEst;
                    corr[2] = aw[0];
                    corr[3] = aw[1];
                    corr[4] = aw[2];
                    corr[5] = rr;
                    VecF<double> geomStats = network.getOptimisationGeometry(network.networkA, lenHist,angHist);
                    VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
//                    VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
                    VecF<double> cndStats = network.getNodeDistribution("A");
//                    network.getRingAreas(a,aSq);
                    VecF<double> eCompare(6);
                    eCompare[0] = ringStats[6];
                    eCompare[1] = SimpleGrapheneEnergy;
                    eCompare[2] = TersoffGrapheneEnergy;
                    eCompare[3] = TriangleRaftEnergy;
                    eCompare[4] = BilayerEnergy;
		    eCompare[5] = BNEnergy;

                    outRingStats.writeRowVector(ringStats);
                    outCorr.writeRowVector(corr);
                    outEnergy.write(energy);
                    outEntropy.writeRowVector(s);
                    outTemperature.write(mcT);
                    outGeometry.writeRowVector(geomStats);
                    outAreas.writeRowVector(a);
                    outAreas.writeRowVector(aSq);
                    outEnergyStats.writeRowVector(eCompare);
//                    outClusterB.writeRowVector(clusters);
                    for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
                    outCndStats.writeRowVector(cndStats);
                    if(mixedLattice) {
                        VecF<double> cluster(3);
                        cluster[0] = network.getMaxCluster("A", 3);
                        cluster[1] = network.getMaxCluster("A", 4);
                        cluster[2] = network.getAssortativity("A");
                        outClusterA.writeRowVector(cluster);
                    }
                }
//                cout<<i<<" "<<network.checkConsistency()<<endl;
                if (i%structureFreq==0 && writeStructures==1) {
                    network.syncCoordinates();
                    //network.write(prefixOut+"_t"+to_string(t)+"_"+to_string(i));
                    network.writeXYZ(prefixOut+"_t"+to_string(t)+"_"+to_string(i));
                }
            }
            --logfile.currIndent;
        }
        --logfile.currIndent;
        logfile.write("Monte Carlo simulation complete");
        logfile.separator();

    }
    else if(runType=="cost"){

        //randomise lattice
        double energy;
        if(costEquilSteps>0){
            network.mc.setTemperature(1e10);
            logfile.write("Running Monte Carlo randomisation");
            ++logfile.currIndent;
            for (int i = 1; i <=costEquilSteps; ++i) {
                moveStatus = network.monteCarloSwitchMove(network.networkT, energy);
                if (i % trackFreq == 0) {
                    double dt = logfile.timeElapsed();
                    string track =
                            to_string(i) + " randomisation moves completed in " + to_string(dt) +" seconds";
                    logfile.write(track);
                    cout << "r" << " " << i << endl;
                }
            }
            --logfile.currIndent;
            logfile.write("Monte Carlo randomisation complete");
            logfile.separator();
        }

        //equilibrate
        logfile.write("Running Monte Carlo equilibration");
        ++logfile.currIndent;
        double cost;
        double costP=costLowLimP;
        double costR=costLowLimR;
        for(int i=1; i<=costInitSteps; ++i){
            moveStatus = network.monteCarloCostSwitchMove(cost,energy,costP,costR);
            accepted += moveStatus[0];
            optCodes[moveStatus[1]] += 1;
            optIterations += moveStatus[2];
//            cout << i << endl;
            if (i % trackFreq == 0) {
                double dt = logfile.timeElapsed();
                string track =
                        to_string(accepted) + "/" + to_string(i) + " moves accepted/completed in " + to_string(dt) +
                        " seconds";
                logfile.write(track);
                cout << "e" << " " << i << " " << accepted << endl;
            }
            if (i % analysisFreq == 0) {
                VecF<double> ringStats = network.getNodeDistribution("B");
                double r = network.getAssortativity("B");
                double aEst = network.getAboavWeaireEstimate("B");
                VecF<double> aw = network.getAboavWeaire("B");
                VecF<double> s = network.getEntropy("B");
                VecF<double> corr(5);
                VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                corr[0] = r;
                corr[1] = aEst;
                corr[2] = aw[0];
                corr[3] = aw[1];
                corr[4] = aw[2];
                VecF<double> emptyL,emptyA; //dummy histograms
                VecF<double> geomStats = network.getOptimisationGeometry(network.networkT, emptyL,emptyA);
                VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
//                network.getRingAreas(a,aSq);
                outRingStats.writeRowVector(ringStats);
                outCorr.writeRowVector(corr);
                outEnergy.write(energy);
                outEntropy.writeRowVector(s);
                outTemperature.write(costT);
                outGeometry.writeRowVector(geomStats);
                outAreas.writeRowVector(a);
                outAreas.writeRowVector(aSq);
                for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
                outClusterB.writeRowVector(clusters);
                if(mixedLattice) {
                    VecF<double> cluster(3);
                    cluster[0] = network.getMaxCluster("A", 3);
                    cluster[1] = network.getMaxCluster("A", 4);
                    cluster[2] = network.getAssortativity("A");
                    outClusterA.writeRowVector(cluster);
                }
            }
        }
        --logfile.currIndent;
        logfile.write("Monte Carlo equilibration complete");
        logfile.separator();

        //search
        logfile.write("Running Monte Carlo simulation");
        ++logfile.currIndent;
        int nR = ceil((costUpLimR - costLowLimR) / costIncR);
        int nP = ceil((costUpLimP - costLowLimP) / costIncP);
        for(int r=0; r<=nR; ++r){
            costR=costLowLimR+costIncR*r;
            logfile.write("Searching r =",costR);
            ++logfile.currIndent;
            for (int p = 0; p <= nP; ++p) {
                if(r%2==0) costP=costLowLimP+costIncP*p;
                else costP=costUpLimP-costIncP*p;
                network.mcCost.setEnergy(numeric_limits<double>::infinity());
                logfile.write("Searching p =",costP);
                ++logfile.currIndent;
                for(int i=1; i<=costSteps; ++i){
                    moveStatus = network.monteCarloCostSwitchMove(cost,energy,costP,costR);
                    accepted += moveStatus[0];
                    optCodes[moveStatus[1]] += 1;
                    optIterations += moveStatus[2];
//                    cout << i << endl;
                    if (i % trackFreq == 0) {
                        double dt = logfile.timeElapsed();
                        string track =
                                to_string(accepted) + "/" + to_string(i) + " moves accepted/completed in " + to_string(dt) +
                                " seconds";
                        logfile.write(track);
                        cout << costR << " " << costP << " " << i << " " << accepted << endl;
                    }
                    if (i % analysisFreq == 0) {
                        VecF<double> ringStats = network.getNodeDistribution("B");
                        double r = network.getAssortativity("B");
                        double aEst = network.getAboavWeaireEstimate("B");
                        VecF<double> aw = network.getAboavWeaire("B");
                        VecF<double> s = network.getEntropy("B");
                        VecF<double> corr(5);
                        VecF<double> a(maxRingSize+1),aSq(maxRingSize+1);
                        corr[0] = r;
                        corr[1] = aEst;
                        corr[2] = aw[0];
                        corr[3] = aw[1];
                        corr[4] = aw[2];
                        VecF<double> geomStats = network.getOptimisationGeometry(network.networkT,lenHist,angHist);
                        VecF< VecF<int> > edgeDist = network.getEdgeDistribution("B");
                        VecF<int> clusters = network.getMaxClusters("B",minRingSize,maxRingSize);
//                        network.getRingAreas(a,aSq);
                        outRingStats.writeRowVector(ringStats);
                        outCorr.writeRowVector(corr);
                        outEnergy.write(energy);
                        outEntropy.writeRowVector(s);
                        outTemperature.write(costT);
                        outGeometry.writeRowVector(geomStats);
                        outAreas.writeRowVector(a);
                        outAreas.writeRowVector(aSq);
                        for(int j=0; j<edgeDist.n; ++j) outEmatrix.writeRowVector(edgeDist[j]);
                        outClusterB.writeRowVector(clusters);
                        if(mixedLattice) {
                            VecF<double> cluster(3);
                            cluster[0] = network.getMaxCluster("A", 3);
                            cluster[1] = network.getMaxCluster("A", 4);
                            cluster[2] = network.getAssortativity("A");
                            outClusterA.writeRowVector(cluster);
                        }
                    }
                    if (i%structureFreq==0 && writeStructures==0) {
                        network.syncCoordinates();
                        network.write(prefixOut+"_t"+to_string(r)+"_"+to_string(i));
                        network.writeXYZ(prefixOut+"_t"+to_string(r)+"_"+to_string(i));
                   }
                }
                --logfile.currIndent;
            }
            --logfile.currIndent;
        }

        --logfile.currIndent;
        logfile.write("Monte Carlo simulation complete");
        logfile.separator();
    }

    //Write total analysis
    for(int i=0; i<10000; ++i){
        VecF<double> hist(4);
        hist[0]=i*4.0/10000.0;
        hist[1]=lenHist[i];
        hist[2]=i*2*M_PI/10000.0;
        hist[3]=angHist[i];
        outGeomHist.writeRowVector(hist);
    }

    if (isSimpleGraphene) {
        cout << "Writing Simple Graphene Results" << endl;
        network.SimpleGraphene.write_data(0);
        network.SimpleGraphene.write_restart(0);
    }
    if (isTriangleRaft) {
        cout << "Writing Triangle Raft Results" << endl;
        network.Triangle_Raft.write_data(1);
        network.Triangle_Raft.write_restart(1);
    }
    if (isBN) {
        cout << "Writing BN Results" << endl;
        network.BN.write_data(4);
        network.BN.write_restart(4);
    }

    //Check network
    logfile.write("Simulation diagnostics");
    ++logfile.currIndent;
    cout << "Check Consistency" << endl;
    bool consistent=network.checkConsistency();
    cout << "Check Convexity" << endl;
    bool convex=network.checkConvexity();
    logfile.write("Network consistent:",consistent);
    logfile.write("Rings convex:",convex);
    logfile.write("Monte Carlo acceptance:",(double)accepted/mcSteps);
    logfile.write("Geometry optimisation codes:");
    ++logfile.currIndent;
    logfile.write("Converged: ",optCodes[0]);
    logfile.write("Converged (zero force): ",optCodes[1]);
    logfile.write("Unconverged: ",optCodes[2]);
    logfile.write("Failed (overlapping): ",optCodes[3]);
    logfile.write("Failed (initially non-convex): ",optCodes[4]);
    --logfile.currIndent;
    logfile.write("Geometry optimisation average iterations:",optIterations/vSum(optCodes));
    --logfile.currIndent;
    logfile.write("Diagnostics complete");
    logfile.separator();

    //Write files
    logfile.write("Writing files");
    ++logfile.currIndent;
    network.wrapCoordinates();
    network.syncCoordinates();
    network.write(prefixOut);
    network.writeXYZ(prefixOut);
    cout << prefixOut<< endl;
//    network.writeBilayer(prefixOut);
//    network.writeBilayer(Folder, lj_cutoff);
    //network.writeBN(prefixOut);
    logfile.write("Writing network files");
    --logfile.currIndent;
    logfile.write("Files written");
    logfile.separator();

    //Close files
    logfile.datetime("Simulation complete at: ");
    return 0;
}
