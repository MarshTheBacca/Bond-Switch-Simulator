#include "linked_network.h"
#include <iostream>
#include <chrono>
#include <ctime>
#include <unistd.h>
#include <omp.h>

//Default constructor
LinkedNetwork::LinkedNetwork() {
    //
}

// Import as bilayer?



//Construct with defined A network, B generated as the dual of A
LinkedNetwork::LinkedNetwork(int nodesA, string latticeA, int minA, int maxA, int minB, int maxB) {

    //Initialise lattices
    if(latticeA=="square" || latticeA=="triangular"){
        networkA=Network(nodesA,latticeA,maxA);
        networkB=networkA.constructDual(maxB);
    }
    else if(latticeA=="hexagonal"){
        networkB=Network(nodesA,"triangular",maxB);
        networkA=networkB.constructDual(maxA);
        rescale(sqrt(3.0));
    }
    else if(latticeA=="cairo"){
        networkB=Network(nodesA,"snubsquare",maxB);
        networkA=networkB.constructDual(maxA);
        rescale(6.0/(3.0+sqrt(3.0)));
        rescale((3.0+sqrt(3.0))/(8.0-2*sqrt(3.0)));
    }
    else if(latticeA=="alt_square"){
        networkA=Network(nodesA,"altsquare",maxA);
        networkB=networkA.constructDual(maxB);
    }
    else if(latticeA=="cubic"){
        networkA=Network(nodesA,"cubic",maxA);
        networkB=networkA.constructDual(maxB);
    }
    else if(latticeA=="inv_cubic"){
        networkB=Network(nodesA,"cubic",maxB);
        networkA=networkB.constructDual(maxA);
    }
    else if(latticeA=="geodesic"){
        networkA=Network(nodesA,"geodesic",maxA);
        networkB=networkA.constructDual(maxB);
    }
    else if(latticeA=="goldberg"){
        networkB=Network(nodesA,"geodesic",maxB);
        networkA=networkB.constructDual(maxA);
    }
    else if(latticeA.substr(0,3)=="mix"){
        double mix=stod(latticeA.substr(4,latticeA.length()));
        networkB=Network(nodesA,"mixTS",maxB,mix);
        networkA=networkB.constructDual(maxA);
    }
    networkB.generateAuxConnections(networkA,0);

    if (minA<=minACnxs){ minACnxs=minA;}
    else {cout << "Initial network does not fit within allowed coordination numbers" << endl << minA << " vs " << minACnxs << endl;}
    if (maxA>=maxACnxs){ maxACnxs=maxA;}
    else {cout << "Initial network does not fit within allowed coordination numbers" << endl << maxA << " vs " << maxACnxs << endl;}
    if (minB<=minBCnxs){ minBCnxs=minB;}
    else {cout << "Initial network does not fit within allowed coordination numbers" << endl << minB << " vs " << minBCnxs << endl;}
    if (maxB>=maxBCnxs){ maxBCnxs=maxB;}
    else {cout << "Initial network does not fit within allowed coordination numbers" << endl << maxB << " vs " << maxBCnxs << endl;}
//    minACnxs=minA;
//    maxACnxs=maxA;
//    minBCnxs=minB;
//    maxBCnxs=maxB;
}

//Construct by loading networks from files
LinkedNetwork::LinkedNetwork(string prefixFolderIn, string prefixFileIn, string prefixFolderOut, int minA, int maxA, int minB, int maxB, bool isSimpleGraphene, bool isTriangleRaft, bool isBilayer, bool isTersoffGraphene, bool isBN, bool restartLammps) {
    string prefix = prefixFolderIn+'/'+prefixFileIn;
    networkA=Network(prefix+"_A", maxB, maxB);
    //cout << "Network A Nodes 63 : " << networkA.nodes.n << endl;
    int ACnxs=0, BCnxs=0;
    minACnxs=10;
    maxACnxs=0;
    for(int i=0; i<networkA.nodes.n; ++i) {
        ACnxs = networkA.nodes[i].netCnxs.n;
        if (ACnxs < minACnxs) minACnxs = ACnxs;
        if (ACnxs > maxACnxs) maxACnxs = ACnxs;
    }

    //cout << "Network A Max A Cnxs : " << maxACnxs << endl;
    networkB=Network(prefix+"_B", maxB, maxB);
    minBCnxs=10;
    maxBCnxs=0;
    for(int i=0; i<networkB.nodes.n; ++i) {
        BCnxs = networkB.nodes[i].netCnxs.n;
        if (BCnxs < minBCnxs) minBCnxs = BCnxs;
        if (BCnxs > maxBCnxs) maxBCnxs = BCnxs;
    }


    if (minA<=minACnxs){ minACnxs=minA;}
    else {cout << "Initial network does not fit within allowed min node coordination numbers" << endl << minA << " vs " << minACnxs << endl;}
    if (maxA>=maxACnxs){ maxACnxs=maxA;}
    else {cout << "Initial network does not fit within allowed max node coordination numbers" << endl << maxA << " vs " << maxACnxs << endl;}
    if (minB<=minBCnxs){ minBCnxs=minB;}
    else {
        cout << "Initial network does not fit within allowed min ring coordination numbers" << endl << minB << " vs " << minBCnxs << endl;
        minBCnxs = minB;
    }
    if (maxB>=maxBCnxs){ maxBCnxs=maxB;}
    else {
        cout << "Initial network does not fit within allowed max ring coordination numbers" << endl << maxB << " vs " << maxBCnxs << endl;
        maxBCnxs=maxB;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "*****************************************************" << endl << "            CREATING NETWORK T " << endl;
    cout << "*****************************************************" << endl;
    if (restartLammps && isTriangleRaft){
        networkT = Network(prefix+"_Si2O3", maxB, maxB);
    }
    else{
        networkT = Network(networkA.nodes,networkA.pb, networkA.rpb, "t", maxACnxs, maxACnxs);
    }


    //  Linked Network can contain lammps objects for minimisation
    //  These will not necessarily have the same coordinate systems, bond orders, etc etc.


    if (isSimpleGraphene) SimpleGraphene          =LammpsObject(0, prefixFolderIn, prefixFolderOut);
    if (isTersoffGraphene) TersoffGraphene         =LammpsObject(3, prefixFolderIn, prefixFolderOut);
    if (isTriangleRaft) Triangle_Raft           =LammpsObject(1, prefixFolderIn, prefixFolderOut);
    if (isBilayer) Bilayer                 =LammpsObject(2, prefixFolderIn, prefixFolderOut);
    if (isBN) BN                         =LammpsObject(4, prefixFolderIn, prefixFolderOut);


}

void LinkedNetwork::pushPrefix(string prefixin, string prefixout){
    prefixIn=prefixin;
    prefixOut=prefixout;
}


void LinkedNetwork::findIgnoredRings(bool fixed_rings, string filename){
    if (fixed_rings){
        string line;
        istringstream ss("");
        cout << filename << endl;
        ifstream auxFile(filename+".dat", ios::in);
        int NoIgnoredRings;
        getline(auxFile,line);
        cout << line << endl;
        istringstream(line)>>NoIgnoredRings;
        Fixed_Ring.setSize(NoIgnoredRings);
        cout << "Fixing " << NoIgnoredRings <<" Rings ..." << endl;

//    VecR<int> Ignored(NoIgnoredRings);
        for (int i=0;i<NoIgnoredRings;++i){
            getline(auxFile,line);
            istringstream(line)>>Fixed_Ring[i];
            cout << Fixed_Ring[i] << " ";
        }
        cout << endl;

    }
    else{
        Fixed_Ring.setSize(0);
    }
}

//Set up potential model with single angle and bond parameter set
void LinkedNetwork::initialisePotentialModel(Network network, double ak, double bk, double ck, int convex) {


    /*
    //Make copy of lattice A coordinates
    if(network.geometryCode=="2DE"){
        crds=VecF<double>(2*network.nodes.n);
        for(int i=0; i<network.nodes.n; ++i){
            crds[2*i]=network.nodes[i].crd[0];
            crds[2*i+1]=network.nodes[i].crd[1];
        }
    }
    else if (network.geometryCode=="2DEtr"){
        crds=VecF<double>(2*network.nodes.n);
        for(int i=0; i<network.nodes.n; ++i){
            crds[2*i]=network.nodes[i].crd[0];
            crds[2*i+1]=network.nodes[i].crd[1];
        }
    }
    else{
        crds=VecF<double>(3*network.nodes.n);
        for(int i=0; i<network.nodes.n; ++i){
            crds[3*i]=network.nodes[i].crd[0];
            crds[3*i+1]=network.nodes[i].crd[1];
            crds[3*i+2]=network.nodes[i].crd[2];
        }
    }
    */

    //Make copy of lattice A coordinates
    if(networkA.geometryCode=="2DE"){
        crds=VecF<double>(2*networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i){
            crds[2*i]=networkA.nodes[i].crd[0];
            crds[2*i+1]=networkA.nodes[i].crd[1];
        }
    }
    else{
        crds=VecF<double>(3*networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i){
            crds[3*i]=networkA.nodes[i].crd[0];
            crds[3*i+1]=networkA.nodes[i].crd[1];
            crds[3*i+2]=networkA.nodes[i].crd[2];
        }
    }


    //Initialise potential model parameters
    //Angle parameters
    potParamsA=VecF<double>(6); //for 3 and 4 coordinate
    potParamsA[0]=ak;
    potParamsA[1]=cos(2*M_PI/3.0);
    potParamsA[2]=ak;
    potParamsA[3]=cos(2*M_PI/4.0);

    //Bond parameters
    potParamsB=VecF<double>(3);
    potParamsB[0]=bk;
    potParamsB[1]=1.0;

    //Geometry constraint parameters
    potParamsC=VecF<double>(2);
    potParamsC[0]=ck; //k, r0 updated through optimal projection

    //Line intersection parameters
    potParamsD=VecF<int>(2);
    potParamsD[0]=1;
    potParamsD[1]=convex;
}

//Set up geometry optimisation parameters
void LinkedNetwork::initialiseGeometryOpt(int iterations, double tau, double tolerance, int localExtent) {

    goptParamsA=VecF<int>(2);
    goptParamsA[0]=iterations;
    goptParamsA[1]=localExtent;
    goptParamsB=VecF<double>(2);
    goptParamsB[0]=tau;
    goptParamsB[1]=tolerance;
}

//Set up monte carlo and random number generators
void LinkedNetwork::initialiseMonteCarlo(Network network, double temperature, int seed, bool globalOpt) {
    cout << "Running Global Optimisation " << endl;
    if(globalOpt) globalGeometryOptimisation(true,true, networkA);
    cout << "Collecting global energy " << endl;
    double energy=globalPotentialEnergy(potParamsD[0],potParamsD[1], networkA);
    cout << "Global Optimisation Energy : " << energy << endl;
    mc=Metropolis(seed,temperature,energy);
    mtGen.seed(seed);
}

//Set up cost function parameters and monte carlo
void LinkedNetwork::initialiseCostFunction(double temperature, int seed, double pk, double rk) {

    costParams=VecF<double>(2);
    costParams[0]=pk;
    costParams[1]=rk;
    double cost=numeric_limits<double>::infinity();
    mcCost=Metropolis(seed+1,temperature,cost);
}

//Rescale lattice dimensions
void LinkedNetwork::rescale(double scaleFactor) {
    networkA.rescale(scaleFactor);
    networkB.rescale(scaleFactor);
}

//Project lattice onto different geometry
void LinkedNetwork::project(string projType, double param) {
    if(projType=="sphere"){
        networkA.project(projType,param);
        networkB.project(projType,param);
    }
}

//Perform defined monte carlo moves to make crystal
void LinkedNetwork::makeCrystal(string crystalCode, string lattice) {

    if(lattice=="hexagonal") {
        if (crystalCode == "var1" || crystalCode == "var2") {
            int n = sqrt(networkB.nodes.n);
            for (int i = 0; i < n; i += 2) {
                for (int j = 0; j < n; j += 2) {
                    int u = n * i + j + (i % 4) / 2;
                    int v = n * i + (j + (i % 4) / 2 + 1) % n;
                    VecR<int> common = vCommonValues(networkB.nodes[u].dualCnxs, networkB.nodes[v].dualCnxs);
                    int a = common[0];
                    int b = common[1];
                    VecF<int> switchIdsA, switchIdsB, switchIdsT;
                    generateSwitchIds34(33, switchIdsA, switchIdsB, switchIdsT, a, b, u, v);
                    switchCnx33(switchIdsA, switchIdsB, switchIdsT);
                    localGeometryOptimisation(a, b, 5, false, false);
                }
            }
            if (crystalCode == "var2") {
                for (int i = 0; i < n; i += 2) {
                    for (int j = 0; j < n; j += 2) {
                        int u = n * i + j + ((i + 2) % 4) / 2;
                        int v = n * i + (j + ((i + 2) % 4) / 2 + 1) % n;
                        VecR<int> common = vCommonValues(networkB.nodes[u].dualCnxs, networkB.nodes[v].dualCnxs);
                        int a = common[0];
                        int b = common[1];
                        VecF<int> switchIdsA, switchIdsB, switchIdsT;
                        generateSwitchIds34(33, switchIdsA, switchIdsB, switchIdsT, a, b, u, v);
                        switchCnx33(switchIdsA, switchIdsB, switchIdsT);
                        localGeometryOptimisation(a, b, 5, false, false);
                    }
                }
            }
            globalGeometryOptimisation(false, false, networkA);
            double energy=globalPotentialEnergy(false,false, networkA);
            mc.setEnergy(energy);
        }
    }

}

//Project lattice onto different geometry with optimal parameters
void LinkedNetwork::optimalProjection(string projType) {

    if(projType=="sphere"){
        //Geometry optimise changing sphere radius until hit minimum in energy

        /* Find initial radius to nearest 1
         * 1) may be multiple minima so try all values in range
         * 2) find lowest value and take limits as radii either side */
        int searchLim=20;
        VecF<double> saveCrdsA=crds;
        VecF<double> energies(20);
        energies[0]=numeric_limits<double>::infinity();
        double radius=1.0;
        for(int i=1; i<searchLim; ++i){
            networkA.project(projType,radius);
            potParamsC[1]=radius;
            for(int j=0; j<networkA.nodes.n; ++j){
                crds[3*j]=networkA.nodes[j].crd[0];
                crds[3*j+1]=networkA.nodes[j].crd[1];
                crds[3*j+2]=networkA.nodes[j].crd[2];
            }
            globalGeometryOptimisation(false,false,networkA);
            energies[i]=globalPotentialEnergy(false,false,networkA);
            cout<<radius<<" "<<energies[i]<<endl;
            radius+=1.0;
            crds=saveCrdsA;
        }
        int id0;
        double e0=energies[0],e1;
        double lowerLim, upperLim, minRadius, minEnergy;
        for(int i=1; i<searchLim; ++i){
            if(energies[i]<e0){
                e0=energies[i];
                id0=i;
            }
        }
        if(id0==searchLim-1) throw string("Initial spherical minimisation reached search limit");
        else{
            lowerLim=id0-1.0;
            minRadius=id0;
            upperLim=id0+1.0;
            networkA.project(projType,minRadius);
            potParamsC[1]=minRadius;
            for(int i=0; i<networkA.nodes.n; ++i){
                crds[3*i]=networkA.nodes[i].crd[0];
                crds[3*i+1]=networkA.nodes[i].crd[1];
                crds[3*i+2]=networkA.nodes[i].crd[2];
            }
            saveCrdsA=crds;
            globalGeometryOptimisation(false,false, networkA);
            minEnergy=globalPotentialEnergy(false,false, networkA);
        }

        /* Refine minimimum
         * 1) could have multiple minima if small - hence reset coordinates as can cause instability
         * 2) search between lower and upper limits until pass through minimum
         * 3) decrease search increment and search again */
        double radiusInc=0.1;
        for(int i=0; i<3; ++i){
            radius=lowerLim;
            e0=numeric_limits<double>::infinity();
            for(;;){
                networkA.project(projType,radius);
                for(int i=0; i<networkA.nodes.n; ++i){
                    crds[3*i]=networkA.nodes[i].crd[0];
                    crds[3*i+1]=networkA.nodes[i].crd[1];
                    crds[3*i+2]=networkA.nodes[i].crd[2];
                }
                potParamsC[1]=radius;
                globalGeometryOptimisation(false,false,networkA);
                e1=globalPotentialEnergy(false,false,networkA);
                cout<<radius<<" "<<e1<<" "<<minEnergy<<endl;
                if(e1>e0 && e0<=minEnergy){
                    lowerLim=radius-2*radiusInc;
                    minRadius=radius-radiusInc;
                    upperLim=radius;
                    minEnergy=e0;
                    break;
                }
                else e0=e1;
                radius+=radiusInc;
                if(radius>upperLim) break;
            }
            radiusInc/=10.0;
        }

        //Optimise with minimum radius
//        potParamsC[0]=0.0;
        potParamsC[1]=minRadius;
        globalGeometryOptimisation(false,false, networkA);
    }
}

//Select nodes forming a random edge in lattice A, and linked nodes in lattice B, only for 3/4 coordinate nodes
int LinkedNetwork::pickSpiralCnx34(int &a, int &b, int &u, int &v, mt19937 &gen) {
    bool includesFixed=false;

    int n0 = a;
    int cnd0 = networkA.nodes[n0].netCnxs.n;
    int n1, n1Ref;
    double rn1 = -1;
    VecF<double> aCnxs(3), rFixedLocal(3);

    //aCnxs = networkA.nodes[n0].netCnxs;
    for (int i=0;i<cnd0;++i){
        rFixedLocal[i] = rFixed[networkA.nodes[n0].netCnxs[i]];
    }
    /*
    //aCnxs = networkA.nodes[n0].netCnxs;
    for (int i=0;i<cnd0;++i){
        for (int j=0; j<networkB.nodes[Fixed_Ring[0]].dualCnxs.n; ++j) {
            if (networkA.nodes[n0].netCnxs[i] == networkB.nodes[Fixed_Ring[0]].dualCnxs[j]) {
                includesFixed = true;
            }
        }
        if (!includesFixed) rFixedLocal[i] = rFixed[networkA.nodes[n0].netCnxs[i]];
        else                rFixedLocal[i]=-1;
        includesFixed = false;
    }
    */

    if (rFixedLocal[0]<0 && rFixedLocal[1]<0 && rFixedLocal[2]<0) return 33;
    while (rn1<0){
        uniform_int_distribution<int> randomCnx(0, cnd0 - 1);

        n1Ref = randomCnx(gen);
        rn1 = rFixedLocal[n1Ref];
        n1 = networkA.nodes[n0].netCnxs[n1Ref];
    }
    a = n0, b = n1;

    VecR<int> common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[b].dualCnxs);
    uniform_int_distribution<int> randomDirection(0, 1);
    if(common.n==2){
        int randIndex=randomDirection(gen);
        u=common[randIndex];
        v=common[1-randIndex];
    }
    else {
        wrapCoordinates();
        syncCoordinates();
        write("debug");
        cout << a << " " << b << " " << common << endl;
        throw string("Error in random connection - incorrect dual ids");
    }

    int cnxType=33;

    return cnxType;
}

//Select nodes forming a random edge in lattice A, and linked nodes in lattice B, only for 3/4 coordinate nodes
int LinkedNetwork::pickDiscreteCnx34(int &a, int &b, int &u, int &v, mt19937 &gen) {

    //pick random node and one of its random connections

    cout << "Running Discrete Distribution" << endl;
    discrete_distribution<> randomNode(weights.begin(), weights.end());
    cout << "Discrete distribution created" << endl;
//    uniform_int_distribution<int> randomNode(0, networkA.nodes.n - 1);
    bool picking_acceptable_ring = true;
    int cnxType;
    bool includesFixed;
    while (picking_acceptable_ring) {
        includesFixed=false;


        // MAYBE EXCLUDE SOME ATOMS ? //


        int n0 = randomNode(gen);
//        n0 = 0;
//        cout << endl;
//        cout << "n0 : " << n0 << endl;

        int cnd0 = networkA.nodes[n0].netCnxs.n;

//        cout << "cnd0 : " << cnd0 << endl;

        uniform_int_distribution<int> randomCnx(0, cnd0 - 1);
        int n1 = networkA.nodes[n0].netCnxs[randomCnx(gen)];
        int cnd1 = networkA.nodes[n1].netCnxs.n;

/*        cout << "Connecting Node " << n0 << " and " << n1 << endl;
        cout << n0 << " connected to ";
        for (int i=0;i<3;++i)
            cout <<  networkA.nodes[n0].netCnxs[i] << " ";
        cout << endl;

        cout << n1 << " connected to ";
        for (int i=0;i<3;++i)
            cout <<  networkA.nodes[n1].netCnxs[i] << " ";
        cout << endl;
*/
        //find connection type and assign a,b

        if (cnd0 == 3 && cnd1 == 3) {
            cnxType = 33;
            a = n0;
            b = n1;
        } else if (cnd0 == 4 && cnd1 == 4) {
            cnxType = 44;
            a = n0;
            b = n1;
        } else if (cnd0 == 3 && cnd1 == 4) {
            cnxType = 43;
            a = n1;
            b = n0;
        } else if (cnd0 == 4 && cnd1 == 3) {
            cnxType = 43;
            a = n0;
            b = n1;
        } else throw string("Error in random connection - incorrect coordinations");

        /*    cout << "a/b " << a << " " << b << endl;
            cout << endl;
            cout << endl;
            cout << n0 << " connected to ";
            for (int i=0;i<networkA.nodes[n0].dualCnxs.n;++i)
                cout <<  networkA.nodes[n0].dualCnxs[i] << " ";
            cout << endl;

            cout << n1 << " connected to ";
            for (int i=0;i<networkA.nodes[n1].dualCnxs.n;++i)
                cout <<  networkA.nodes[n1].dualCnxs[i] << " ";
            cout << endl;
        */

        //get nodes in dual in random orientation
        uniform_int_distribution<int> randomDirection(0, 1);
/*
        cout << n0 << " dual : " << endl;
        for (int i=0;i<networkA.nodes[a].dualCnxs.n;++i){
            cout << networkA.nodes[a].dualCnxs[i] << " ";
        }
        cout << endl;

        cout << n1 << " dual : " << endl;
        for (int i=0;i<networkA.nodes[b].dualCnxs.n;++i){
            cout << networkA.nodes[b].dualCnxs[i] << " ";
        }
        cout << endl;
*/
        VecR<int> common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[b].dualCnxs);

/*        cout << "Common : " << endl;
        for (int i=0;i<common.n;++i){
            cout << common[i] << " ";
        }
        cout << endl;
*/
        if(common.n==2){
            int randIndex=randomDirection(gen);
            u=common[randIndex];
            v=common[1-randIndex];
        }
        else {
            wrapCoordinates();
            syncCoordinates();
            write("debug");
            cout << a << " " << b << " " << common << endl;
            throw string("Error in random connection - incorrect dual ids");
        }
        /*
        cout << endl;
        cout << endl;
        cout << "Fixed Ring " << Fixed_Ring[0] << endl;
        for (int i=0; i<networkB.nodes[Fixed_Ring[0]].dualCnxs.n;++i){
            cout << networkB.nodes[Fixed_Ring[0]].dualCnxs[i] << " ";
        }
        cout << endl;
        */
        int fixed_ring;
        for (int i=0; i<Fixed_Ring.n;++i){
            fixed_ring = Fixed_Ring[i];
            for (int j=0; j<networkB.nodes[fixed_ring].dualCnxs.n; ++j){

                if (a==networkB.nodes[fixed_ring].dualCnxs[j] || b==networkB.nodes[fixed_ring].dualCnxs[j]){
//                    cout << "illegal ring " << endl;
                    includesFixed=true;
//                    cout << a << " " << b << " " <<  networkB.nodes[i].dualCnxs[j] << endl;

                }
                else{
//                    cout << "Dual " << networkB.nodes[fixed_ring].dualCnxs[j] << " isn't "<< a << " or " << b << endl;
                }
            }

//            if (u==Fixed_Ring[i] || v==Fixed_Ring[i]) {includesFixed=true;}
        }
        /*
        cout << endl;
        cout << endl;
         */
        if (!includesFixed) {
            picking_acceptable_ring=false;
//            cout << "Breaking the loop " << endl;
        }
    }

    return cnxType;
}

//Select nodes forming a random edge in lattice A, and linked nodes in lattice B, only for 3/4 coordinate nodes
int LinkedNetwork::pickRandomCnx34(int &a, int &b, int &u, int &v, mt19937 &gen) {

    //pick random node and one of its random connections
    uniform_int_distribution<int> randomNode(0, networkA.nodes.n - 1);
    bool picking_acceptable_ring = true;
    int cnxType;
    bool includesFixed;
    while (picking_acceptable_ring) {
        includesFixed=false;


        // MAYBE EXCLUDE SOME ATOMS ? //


        int n0 = randomNode(gen);
//        n0 = 0;
//        cout << endl;
//        cout << "n0 : " << n0 << endl;

        int cnd0 = networkA.nodes[n0].netCnxs.n;

//        cout << "cnd0 : " << cnd0 << endl;

        uniform_int_distribution<int> randomCnx(0, cnd0 - 1);
        int n1 = networkA.nodes[n0].netCnxs[randomCnx(gen)];
        int cnd1 = networkA.nodes[n1].netCnxs.n;

/*        cout << "Connecting Node " << n0 << " and " << n1 << endl;
        cout << n0 << " connected to ";
        for (int i=0;i<3;++i)
            cout <<  networkA.nodes[n0].netCnxs[i] << " ";
        cout << endl;

        cout << n1 << " connected to ";
        for (int i=0;i<3;++i)
            cout <<  networkA.nodes[n1].netCnxs[i] << " ";
        cout << endl;
*/
        //find connection type and assign a,b

        if (cnd0 == 3 && cnd1 == 3) {
            cnxType = 33;
            a = n0;
            b = n1;
        } else if (cnd0 == 4 && cnd1 == 4) {
            cnxType = 44;
            a = n0;
            b = n1;
        } else if (cnd0 == 3 && cnd1 == 4) {
            cnxType = 43;
            a = n1;
            b = n0;
        } else if (cnd0 == 4 && cnd1 == 3) {
            cnxType = 43;
            a = n0;
            b = n1;
        } else throw string("Error in random connection - incorrect coordinations");

        /*    cout << "a/b " << a << " " << b << endl;
            cout << endl;
            cout << endl;
            cout << n0 << " connected to ";
            for (int i=0;i<networkA.nodes[n0].dualCnxs.n;++i)
                cout <<  networkA.nodes[n0].dualCnxs[i] << " ";
            cout << endl;

            cout << n1 << " connected to ";
            for (int i=0;i<networkA.nodes[n1].dualCnxs.n;++i)
                cout <<  networkA.nodes[n1].dualCnxs[i] << " ";
            cout << endl;
        */

        //get nodes in dual in random orientation
        uniform_int_distribution<int> randomDirection(0, 1);
/*
        cout << n0 << " dual : " << endl;
        for (int i=0;i<networkA.nodes[a].dualCnxs.n;++i){
            cout << networkA.nodes[a].dualCnxs[i] << " ";
        }
        cout << endl;

        cout << n1 << " dual : " << endl;
        for (int i=0;i<networkA.nodes[b].dualCnxs.n;++i){
            cout << networkA.nodes[b].dualCnxs[i] << " ";
        }
        cout << endl;
*/
        VecR<int> common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[b].dualCnxs);

/*        cout << "Common : " << endl;
        for (int i=0;i<common.n;++i){
            cout << common[i] << " ";
        }
        cout << endl;
*/
        if(common.n==2){
            int randIndex=randomDirection(gen);
            u=common[randIndex];
            v=common[1-randIndex];
        }
        else {
            wrapCoordinates();
            syncCoordinates();
            write("debug");
            cout << a << " " << b << " " << common << endl;
            throw string("Error in random connection - incorrect dual ids");
        }
        /*
        cout << endl;
        cout << endl;
        cout << "Fixed Ring " << Fixed_Ring[0] << endl;
        for (int i=0; i<networkB.nodes[Fixed_Ring[0]].dualCnxs.n;++i){
            cout << networkB.nodes[Fixed_Ring[0]].dualCnxs[i] << " ";
        }
        cout << endl;
        */
        int fixed_ring;
        for (int i=0; i<Fixed_Ring.n;++i){
            fixed_ring = Fixed_Ring[i];
            for (int j=0; j<networkB.nodes[fixed_ring].dualCnxs.n; ++j){

                if (a==networkB.nodes[fixed_ring].dualCnxs[j] || b==networkB.nodes[fixed_ring].dualCnxs[j]){
//                    cout << "illegal ring " << endl;
                    includesFixed=true;
//                    cout << a << " " << b << " " <<  networkB.nodes[i].dualCnxs[j] << endl;

                }
                else{
//                    cout << "Dual " << networkB.nodes[fixed_ring].dualCnxs[j] << " isn't "<< a << " or " << b << endl;
                }
            }

//            if (u==Fixed_Ring[i] || v==Fixed_Ring[i]) {includesFixed=true;}
        }
        /*
        cout << endl;
        cout << endl;
         */
        if (!includesFixed) {
            picking_acceptable_ring=false;
//            cout << "Breaking the loop " << endl;
        }
    }

    return cnxType;
}

// pick weighted Cnx
// pick sequential Cnxs

//Select nodes forming a random edge in lattice A, and linked nodes in lattice B, for coordinations >=2
int LinkedNetwork::pickRandomCnx(int& a, int& b, int& u, int& v, mt19937& gen) {

    //pick random node and one of its random connections
    uniform_int_distribution<int> randomNode(0,networkA.nodes.n-1);
    int n0=randomNode(gen);
    int cnd0=networkA.nodes[n0].netCnxs.n;
    uniform_int_distribution<int> randomCnx(0,cnd0-1);
    int n1=networkA.nodes[n0].netCnxs[randomCnx(gen)];
    int cnd1=networkA.nodes[n1].netCnxs.n;

    //find connection type and assign a,b
    int cnxType=10*cnd0+cnd1;
    a=n0;
    b=n1;
    //get nodes in dual in random orientation
    uniform_int_distribution<int> randomDirection(0,1);
    VecR<int> common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[b].dualCnxs);
    if(common.n==2){
        int randIndex=randomDirection(gen);
        u=common[randIndex];
        v=common[1-randIndex];
    }
    else if(common.n==3){
        //check a-b in ring
        VecR<int> uv(0,common.n);
        for(int i=0; i<common.n; ++i){
            VecR<int> ringCnxs=networkB.nodes[common[i]].dualCnxs;
            for(int j=0; j<ringCnxs.n; ++j){
                int k=(j+ringCnxs.n-1)%ringCnxs.n;
                int l=(j+1)%ringCnxs.n;
                if(ringCnxs[j]==a && ringCnxs[k]==b){
                    uv.addValue(common[i]);
                    break;
                }
                else if(ringCnxs[j]==a && ringCnxs[l]==b){
                    uv.addValue(common[i]);
                    break;
                }
            }
        }
        if(uv.n==2){
            int randIndex=randomDirection(gen);
            u=uv[randIndex];
            v=uv[1-randIndex];
        }
        else{
            wrapCoordinates();
            syncCoordinates();
            write("debug");
            cout<<a<<" "<<b<<" "<<common<<endl;
            throw string("Error in random connection - incorrect dual ids");
        }
    }

    return cnxType;
}

//Generate all ids of nodes in lattices A and B needed for switch move, only for 3/4 coordinate nodes
int LinkedNetwork::generateSwitchIds34(int cnxType, VecF<int> &switchIdsA, VecF<int> &switchIdsB, VecF<int> &switchIdsT, int a, int b, int u, int v) {

    //lots of error checking to remove any potential pathological cases
    if(a==b || u==v){
//        cout<<"Note: skip in switch generation"<<endl;
        return 1;
    }

    if(cnxType==33){
        /* Switch connectivities in lattice and dual
         * 3-3 coordination connection
         * a,b,c,d,e,f are nodes in lattice A
         * u,v,w,x are nodes in lattice B
         *  E      F            V
         *   \    /           / | \
         *   A---B           W  |  X
         *  /     \           \ | /
         * C       D            U
         */

        int errorFlag=0;
        int c,d,e,f;
        int w,x;
        VecR<int> common,common1;
//        common=vCommonValues(networkA.nodes[a].netCnxs,networkB.nodes[u].dualCnxs);
//        common.delValue(b);
//        if(common.n!=1) errorFlag=1;
//        c=common[0];
//        common=vCommonValues(networkA.nodes[b].netCnxs,networkB.nodes[u].dualCnxs);
//        common.delValue(a);
//        if(common.n!=1) errorFlag=2;
//        d=common[0];
//        common=vCommonValues(networkA.nodes[a].netCnxs,networkB.nodes[v].dualCnxs);
//        common.delValue(b);
//        if(common.n!=1) errorFlag=3;
//        e=common[0];
//        common=vCommonValues(networkA.nodes[b].netCnxs,networkB.nodes[v].dualCnxs);
//        common.delValue(a);
//        if(common.n!=1) errorFlag=4;
//        f=common[0];

        c=findAssociatedNodeAB(a,u,b);
        d=findAssociatedNodeAB(b,u,a);
        e=findAssociatedNodeAB(a,v,b);
        f=findAssociatedNodeAB(b,v,a);
        w=findAssociatedNodeAA(a,c,u);
        x=findAssociatedNodeAA(b,d,u);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        int alpha=0, beta=0, gamma=0, delta=0, eta=0;
        for (int i=0; i<networkT.nodes[a].netCnxs.n;++i){
            for (int j=0; j<networkT.nodes[b].netCnxs.n;++j){
                if (networkT.nodes[a].netCnxs[i]==networkT.nodes[b].netCnxs[j]) {alpha=networkT.nodes[a].netCnxs[i];}
            }
            for (int j=0; j<networkT.nodes[c].netCnxs.n;++j){
                if (networkT.nodes[a].netCnxs[i]==networkT.nodes[c].netCnxs[j]) {beta=networkT.nodes[a].netCnxs[i];}
            }
            for (int j=0; j<networkT.nodes[e].netCnxs.n;++j){
                if (networkT.nodes[a].netCnxs[i]==networkT.nodes[e].netCnxs[j]) {gamma=networkT.nodes[a].netCnxs[i];}
            }
        }
        for (int i=0; i<networkT.nodes[b].netCnxs.n;++i){
            for (int j=0; j<networkT.nodes[d].netCnxs.n;++j){
                if (networkT.nodes[b].netCnxs[i]==networkT.nodes[d].netCnxs[j]) {delta=networkT.nodes[b].netCnxs[i];}
            }
            for (int j=0; j<networkT.nodes[f].netCnxs.n;++j){
                if (networkT.nodes[b].netCnxs[i]==networkT.nodes[f].netCnxs[j]) {eta=networkT.nodes[b].netCnxs[i];}
            }
        }
        if (alpha==0) {cout << "alpha broken" << endl;}
        if (beta==0) {cout << "beta broken" << endl;}
        if (gamma==0) {cout << "gamma broken" << endl;}
        if (delta==0) {cout << "delta broken" << endl;}
        if (eta==0) {cout << "eta broken" << endl;}
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//        common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[c].dualCnxs);
//        common1=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[e].dualCnxs);
//        common.delValue(u);
//        common1.delValue(v);
//        if(common.n!=1 || common1.n!=1 || common[0]!=common1[0]) errorFlag=5;
//        w=common[0];
//        common=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[d].dualCnxs);
//        common1=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[f].dualCnxs);
//        common.delValue(u);
//        common1.delValue(v);
//        if(common.n!=1 || common1.n!=1 || common[0]!=common1[0]) errorFlag=5;
//        x=common[0];

        //Additional error checking including preventing two nodes connecting multiple times
//        if(c==d || e==f) errorFlag=6;
//        if(vContains(networkB.nodes[w].netCnxs,x)) errorFlag=7;
//
//        if(errorFlag!=0){
////            cout<<"Note: skip in switch generation 33 with code "<<errorFlag<<endl;
//            return 1;
//        }

        //Additional error checking
        if (c == d || e == f) errorFlag = 6; //can simply be triangle edge sharing pair (not an error)
        //Prevent rings having only two or fewer neighbours
        VecR<int> vCnxs=vUnique(networkB.nodes[v].netCnxs);
        if(vCnxs.n<=3) errorFlag = 10;
        VecR<int> uCnxs=vUnique(networkB.nodes[u].netCnxs);
        uCnxs.delValue(v);
        if(uCnxs.n<=2){
            for(int i=0; i<uCnxs.n; ++i){
                if(uCnxs[i]==x){
                    errorFlag = 10;
                    break;
                }
            }
        }

//        //Prevent two large rings from completely surrounding group
//        VecR<int> uv(2),uvCount(2);
//        uv[0]=u;
//        uv[1]=v;
//        uvCount=vValCount(networkB.nodes[w].netCnxs,uv);
//        if(uvCount[0]>=2 && uvCount[1]>=2){
//            int adjCount=0;
//            VecR<int> cnxs=networkB.nodes[w].netCnxs;
//            for(int i=0; i<cnxs.n; ++i){
//                if(cnxs[i]==u){
//                    int j=(i+cnxs.n-1)%cnxs.n;
//                    int k=(i+1)%cnxs.n;
//                    if(cnxs[j]==v) ++adjCount;
//                    if(cnxs[k]==v) ++adjCount;
//                }
//            }
//            if(adjCount>1) errorFlag = 11;
//        }
//        uvCount=vValCount(networkB.nodes[x].netCnxs,uv);
//        if(uvCount[0]>=2 && uvCount[1]>=2){
//            int adjCount=0;
//            VecR<int> cnxs=networkB.nodes[x].netCnxs;
//            for(int i=0; i<cnxs.n; ++i){
//                if(cnxs[i]==u){
//                    int j=(i+cnxs.n-1)%cnxs.n;
//                    int k=(i+1)%cnxs.n;
//                    if(cnxs[j]==v) ++adjCount;
//                    if(cnxs[k]==v) ++adjCount;
//                }
//            }
//            if(adjCount>1) errorFlag = 11;
//        }
//        if(a==22 && b==18){
//            cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<errorFlag<<endl;
//        }
        if (errorFlag != 0)  return 1;

        //check move will not violate dual connectivity limits
        if(networkB.nodes[u].netCnxs.n==minBCnxs || networkB.nodes[v].netCnxs.n==minBCnxs
            || networkB.nodes[w].netCnxs.n==maxBCnxs || networkB.nodes[x].netCnxs.n==maxBCnxs) return 1;
        else {
            switchIdsA = VecF<int>(6);
            switchIdsB = VecF<int>(4);
            switchIdsT = VecF<int>(11);
            switchIdsA[0] = a;
            switchIdsA[1] = b;
            switchIdsA[2] = c;
            switchIdsA[3] = d;
            switchIdsA[4] = e;
            switchIdsA[5] = f;
            switchIdsB[0] = u;
            switchIdsB[1] = v;
            switchIdsB[2] = w;
            switchIdsB[3] = x;
            switchIdsT[0] = a;
            switchIdsT[1] = b;
            switchIdsT[2] = c;
            switchIdsT[3] = d;
            switchIdsT[4] = e;
            switchIdsT[5] = f;
            switchIdsT[6] = alpha;
            switchIdsT[7] = beta;
            switchIdsT[8] = gamma;
            switchIdsT[9] = delta;
            switchIdsT[10] = eta;


            return 0;
        }
    }
    else if(cnxType==44){
        /* Switch connectivities in lattice and dual
         * 4-4 coordination connection
         * a,b,c,d,e,f,g,h are nodes in lattice A
         * u,v,w,x,y,z are nodes in lattice B
         * a-b, a-c, a-e, a-g
         * b-a, b-d, b-f, b-h
         * a-b share u-v
         * c-a-b-d share u
         * e-a-b-f share v
         * g-a-e share w
         * d-b-h share x
         * g-a-c share y
         * f-b-h share z
         * u-v, u-y, u-x
         * v-u, v-w, v-z
         * w-y, x-z*/
        int errorFlag=0;
        int c,d,e,f,g,h;
        int w,x,y,z;

        VecR<int> common,common1;
        common=vCommonValues(networkA.nodes[a].netCnxs,networkB.nodes[u].dualCnxs);
        common.delValue(b);
        if(common.n!=1) errorFlag=1;
        c=common[0];
        common=vCommonValues(networkA.nodes[b].netCnxs,networkB.nodes[u].dualCnxs);
        common.delValue(a);
        if(common.n!=1) errorFlag=2;
        d=common[0];
        common=vCommonValues(networkA.nodes[a].netCnxs,networkB.nodes[v].dualCnxs);
        common.delValue(b);
        if(common.n!=1) errorFlag=3;
        e=common[0];
        common=vCommonValues(networkA.nodes[b].netCnxs,networkB.nodes[v].dualCnxs);
        common.delValue(a);
        if(common.n!=1) errorFlag=4;
        f=common[0];

        common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[e].dualCnxs);
        common.delValue(v);
        if(common.n!=1) errorFlag=5;
        w=common[0];
        common=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[d].dualCnxs);
        common.delValue(u);
        if(common.n!=1) errorFlag=6;
        x=common[0];
        common=vCommonValues(networkA.nodes[a].dualCnxs,networkA.nodes[c].dualCnxs);
        common.delValue(u);
        if(common.n!=1) errorFlag=7;
        y=common[0];
        common=vCommonValues(networkA.nodes[b].dualCnxs,networkA.nodes[f].dualCnxs);
        common.delValue(v);
        if(common.n!=1) errorFlag=8;
        z=common[0];

        common=networkA.nodes[a].netCnxs;
        common.delValue(b);
        common.delValue(c);
        common.delValue(e);
        if(common.n!=1) errorFlag=9;
        g=common[0];
        common=networkA.nodes[b].netCnxs;
        common.delValue(a);
        common.delValue(d);
        common.delValue(f);
        if(common.n!=1) errorFlag=10;
        h=common[0];

        //Additional error checking including preventing two nodes connecting multiple times
        if(c==d || e==f) errorFlag=6;
        if(vContains(networkB.nodes[w].netCnxs,x)) errorFlag=7;
        if(vContains(networkB.nodes[y].netCnxs,x)) errorFlag=7;
        if(vContains(networkB.nodes[y].auxCnxs,x)) errorFlag=7;
        if(vContains(networkB.nodes[z].netCnxs,w)) errorFlag=7;
        if(vContains(networkB.nodes[z].auxCnxs,w)) errorFlag=7;

        if(errorFlag!=0){
//            cout<<"Note: skip in switch generation 44"<<endl;
            return 1;
        }

        //check move will not violate dual connectivity limits
        if(networkB.nodes[u].netCnxs.n==minBCnxs || networkB.nodes[v].netCnxs.n==minBCnxs
           || networkB.nodes[w].netCnxs.n==maxBCnxs || networkB.nodes[x].netCnxs.n==maxBCnxs) return 1;
        else {
            switchIdsA = VecF<int>(8);
            switchIdsB = VecF<int>(6);
            switchIdsA[0] = a;
            switchIdsA[1] = b;
            switchIdsA[2] = c;
            switchIdsA[3] = d;
            switchIdsA[4] = e;
            switchIdsA[5] = f;
            switchIdsA[6] = g;
            switchIdsA[7] = h;
            switchIdsB[0] = u;
            switchIdsB[1] = v;
            switchIdsB[2] = w;
            switchIdsB[3] = x;
            switchIdsB[4] = y;
            switchIdsB[5] = z;
            return 0;
        }
    }
    else if(cnxType==43) {
        /* Switch connectivities in lattice and dual
         * 4-3 coordination connection
         * a,b,c,d,e,f,g are nodes in lattice A
         * u,v,w,x,y are nodes in lattice B
         * a-b, a-c, a-e, a-g
         * b-a, b-d, b-f
         * a-b share u-v
         * c-a-b-d share u
         * e-a-b-f share v
         * g-a-e share w
         * d-b-f share x
         * g-a-c share y
         * u-v, u-y, u-x
         * v-u, v-w, v-x
         * w-y */
        int errorFlag = 0;
        int c, d, e, f, g;
        int w, x, y;

        VecR<int> common, common1;
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(b);
        if (common.n != 1) errorFlag = 1;
        c = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(a);
        if (common.n != 1) errorFlag = 2;
        d = common[0];
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(b);
        if (common.n != 1) errorFlag = 3;
        e = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(a);
        if (common.n != 1) errorFlag = 4;
        f = common[0];

        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[e].dualCnxs);
        common.delValue(v);
        if (common.n != 1) errorFlag = 5;
        w = common[0];
        common = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[d].dualCnxs);
        common1 = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[f].dualCnxs);
        common.delValue(u);
        common1.delValue(v);
        if (common.n != 1 || common1.n != 1 || common[0] != common1[0]) errorFlag = 5;
        x = common[0];
        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[c].dualCnxs);
        common.delValue(u);
        if (common.n != 1) errorFlag = 7;
        y = common[0];

        common = networkA.nodes[a].netCnxs;
        common.delValue(b);
        common.delValue(c);
        common.delValue(e);
        if (common.n != 1) errorFlag = 9;
        g = common[0];

        //Additional error checking including preventing two nodes connecting multiple times
        if (c == d || e == f) errorFlag = 6; //can simply be triangle edge sharing pair (not an error)
        if(vContains(networkB.nodes[w].netCnxs,x)) errorFlag=7;
        if(vContains(networkB.nodes[y].netCnxs,x)) errorFlag=7;
        if(vContains(networkB.nodes[y].auxCnxs,x)) errorFlag=7;

        if (errorFlag != 0) {
//            cout << "Note: skip in switch generation 43 with error flag" << " " << errorFlag << endl;
            return 1;
        }

        //check move will not violate dual connectivity limits
        if (networkB.nodes[u].netCnxs.n == minBCnxs || networkB.nodes[v].netCnxs.n == minBCnxs
            || networkB.nodes[w].netCnxs.n == maxBCnxs || networkB.nodes[x].netCnxs.n == maxBCnxs)
            return 1;
        else {
            switchIdsA = VecF<int>(7);
            switchIdsB = VecF<int>(5);
            switchIdsA[0] = a;
            switchIdsA[1] = b;
            switchIdsA[2] = c;
            switchIdsA[3] = d;
            switchIdsA[4] = e;
            switchIdsA[5] = f;
            switchIdsA[6] = g;
            switchIdsB[0] = u;
            switchIdsB[1] = v;
            switchIdsB[2] = w;
            switchIdsB[3] = x;
            switchIdsB[4] = y;
            return 0;
        }
    }
    return 0;
}

//Generate all ids of nodes in lattices A and B needed for mix move, only for 3/4 coordinate nodes
int LinkedNetwork::generateMixIds34(int cnxType, VecF<int> &mixIdsA, VecF<int> &mixIdsB, int a, int b, int u, int v) {

    if(cnxType!=43){//cannot decrement either 3 cnd nodes
        return 1;
    }
    else {
        /* Mix connectivities in lattice and dual
         * a,b,c,d,e,f,g are nodes in lattice A
         * u,v,w,x,y are nodes in lattice B
         * a-b, a-c, a-e, a-g
         * b-a, b-d, b-f
         * a-b share u-v
         * c-a-b-d share u
         * e-a-b-f share v
         * g-a-e share w
         * d-b-f share x
         * g-a-c share y
         * u-v, u-y, u-x
         * v-u, v-w, v-x
         * w-y */

        int errorFlag = 0;
        int c, d, e, f, g;
        int w, x, y;

        VecR<int> common, common1;
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(b);
        if (common.n != 1) errorFlag = 1;
        c = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(a);
        if (common.n != 1) errorFlag = 2;
        d = common[0];
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(b);
        if (common.n != 1) errorFlag = 3;
        e = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(a);
        if (common.n != 1) errorFlag = 4;
        f = common[0];

        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[e].dualCnxs);
        common.delValue(v);
        if (common.n != 1) errorFlag = 5;
        w = common[0];
        common = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[d].dualCnxs);
        common1 = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[f].dualCnxs);
        common.delValue(u);
        common1.delValue(v);
        if (common.n != 1 || common1.n != 1 || common[0] != common1[0]) errorFlag = 5;
        x = common[0];
        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[c].dualCnxs);
        common.delValue(u);
        if (common.n != 1) errorFlag = 7;
        y = common[0];

        common = networkA.nodes[a].netCnxs;
        common.delValue(b);
        common.delValue(c);
        common.delValue(e);
        if (common.n != 1) errorFlag = 9;
        g = common[0];

        //Additional error checking including preventing two nodes connecting multiple times
        if (c == d || e == f) errorFlag = 6; //can simply be triangle edge sharing pair (not an error)
        if(vContains(networkB.nodes[u].netCnxs,w)) errorFlag=7;
        if(vContains(networkB.nodes[u].auxCnxs,v)) errorFlag=7;
        if(vContains(networkB.nodes[w].netCnxs,x)) errorFlag=7;
        if(vContains(networkB.nodes[w].auxCnxs,x)) errorFlag=7;

        if (errorFlag != 0) {
//            cout << "Note: skip in switch generation 43 with error flag" << " " << errorFlag << endl;
            return 1;
        }

        //check move will not violate dual connectivity limits
        if (networkB.nodes[v].netCnxs.n == minBCnxs || networkB.nodes[w].netCnxs.n == maxBCnxs) return 1;
        else {
            mixIdsA = VecF<int>(7);
            mixIdsB = VecF<int>(5);
            mixIdsA[0] = a;
            mixIdsA[1] = b;
            mixIdsA[2] = c;
            mixIdsA[3] = d;
            mixIdsA[4] = e;
            mixIdsA[5] = f;
            mixIdsA[6] = g;
            mixIdsB[0] = u;
            mixIdsB[1] = v;
            mixIdsB[2] = w;
            mixIdsB[3] = x;
            mixIdsB[4] = y;
            return 0;
        }
    }
}

//Generate all ids of nodes in lattices A and B needed for mix move, for all coordinations >=2
int LinkedNetwork::generateMixIds(int cnxType, VecF<int> &mixIdsA, VecF<int> &mixIdsB, int a, int b, int u, int v) {

    if(cnxType<30){//cannot mix if first node is 2 coordinate
        return 1;
    }
    else {
        /* Get required node ids in lattice and dual
         * a,b,c,d,e,f are nodes in lattice A
         * u,v,w,x,y,z,xx are nodes in lattice B
         *  E      F            V
         *   \    /          /  |  \
         * --A---B--    XX--X   |   Z
         *  /     \         W   |   Y
         * C       D         \  |  /
         *                      U
        */

        int errorFlag = 0;
        int c, d, e, f;
        int w, x, y, z, xx;

        VecR<int> common, common1;
        //c is defined as sharing a,u but not being b
        c=findAssociatedNodeAB(a, u, b);
        //d is defined as sharing b,u but not being a
        d=findAssociatedNodeAB(b, u, a);
        //e is defined as sharing a,v but not being b
        e=findAssociatedNodeAB(a, v, b);
        //f is defined as sharing b,v but not being a
        f=findAssociatedNodeAB(b, v, a);

        //w is defines as sharing a,c but not being u
        w=findAssociatedNodeAA(a,c,u);
        //x is defines as sharing a,e but not being v
        x=findAssociatedNodeAA(a,e,v);
        //y is defines as sharing b,d but not being u
        y=findAssociatedNodeAA(b,d,u);
        //z is defines as sharing b,f but not being v
        z=findAssociatedNodeAA(b,f,v);

        //Find next node connected to x by inspecting around a
        int nCnxs=networkA.nodes[a].dualCnxs.n;
        for(int i=0; i<nCnxs; ++i){
            if(networkA.nodes[a].dualCnxs[i]==u){
                int j=(i+1)%nCnxs;
                int k=(i-1+nCnxs)%nCnxs;
                if(networkA.nodes[a].dualCnxs[j]==v) xx=networkA.nodes[a].dualCnxs[(j+2)%nCnxs];
                else if(networkA.nodes[a].dualCnxs[k]==v) xx=networkA.nodes[a].dualCnxs[(k-2+nCnxs)%nCnxs];
                else cout<<"Error in xx detection"<<endl;
            }
        }

        //Additional error checking
        if (c == d || e == f) errorFlag = 6; //can simply be triangle edge sharing pair (not an error)
        //Prevent rings having only two or fewer neighbours
        VecR<int> vCnxs=vUnique(networkB.nodes[v].netCnxs);
        if(vCnxs.n<=3) errorFlag = 10;
        VecR<int> uCnxs=vUnique(networkB.nodes[u].netCnxs);
        uCnxs.delValue(v);
        if(uCnxs.n<=2){
            for(int i=0; i<uCnxs.n; ++i){
                if(uCnxs[i]==x){
                    errorFlag = 10;
                    break;
                }
            }
        }
        if (errorFlag != 0)  return 1;

        //check move will not violate connectivity limits
        if (networkB.nodes[v].netCnxs.n == minBCnxs || networkB.nodes[x].netCnxs.n == maxBCnxs
            || networkA.nodes[a].netCnxs.n == minACnxs || networkA.nodes[b].netCnxs.n == maxACnxs) return 1;
        else {
            mixIdsA = VecF<int>(6);
            mixIdsB = VecF<int>(7);
            mixIdsA[0] = a;
            mixIdsA[1] = b;
            mixIdsA[2] = c;
            mixIdsA[3] = d;
            mixIdsA[4] = e;
            mixIdsA[5] = f;
            mixIdsB[0] = u;
            mixIdsB[1] = v;
            mixIdsB[2] = w;
            mixIdsB[3] = x;
            mixIdsB[4] = y;
            mixIdsB[5] = z;
            mixIdsB[6] = xx;
            return 0;
        }
    }
}

int LinkedNetwork::findAssociatedNodeAB(int idA, int idB, int idDel) {

    //Find node that shares idA and idB but is not idDel
    int associated=-1;
    VecR<int> common;
    common = vCommonValues(networkA.nodes[idA].netCnxs, networkB.nodes[idB].dualCnxs);
    common.delValue(idDel);
    if(common.n==1) associated=common[0];
    else{//rare high temperature occurrence as a result of 2-cnd nodes giving ring inside ring
        VecR<int> common1(0,common.n);
        int nCnxs=networkA.nodes[idA].netCnxs.n;
        for(int i=0; i<nCnxs; ++i){
            if(networkA.nodes[idA].netCnxs[i]==idDel){
                int l=networkA.nodes[idA].netCnxs[(i+1)%nCnxs];
                int r=networkA.nodes[idA].netCnxs[(i-1+nCnxs)%nCnxs];
                for(int j=0; j<common.n; ++j){
                    if(common[j]==l) common1.addValue(common[j]);
                    else if(common[j]==r) common1.addValue(common[j]);
                }
                break;
            }
        }
        if(common1.n==1) associated=common1[0];
        else{//even rarer case
            int nCnxs=networkB.nodes[idB].dualCnxs.n;
            for(int i=0; i<nCnxs; ++i){
                if(networkB.nodes[idB].dualCnxs[i]==idDel){
                    int j;
                    if(networkB.nodes[idB].dualCnxs[(i+1)%nCnxs]==idA) j=(i+2)%nCnxs;
                    else if(networkB.nodes[idB].dualCnxs[(i+nCnxs-1)%nCnxs]==idA) j=(i+nCnxs-2)%nCnxs;
                    if(vContains(common1,networkB.nodes[idB].dualCnxs[j])){
                        associated=networkB.nodes[idB].dualCnxs[j];
                        break;
                    }
                }
            }
        }
    }

    if(associated==-1){
        wrapCoordinates();
        syncCoordinates();
        write("debug");
        cout<<idA<<" "<<idB<<" "<<" "<<idDel<<" "<<common<<endl;
        cout<<"ERROR IN ASSOCIATED NODE"<<endl;
        exit(9);
    }

    return associated;
}

int LinkedNetwork::findAssociatedNodeAA(int idA, int idB, int idDel) {

    //Find node that shares idA and idB but is not idDel
    int associated=-1;
    VecR<int> common;
    common = vCommonValues(networkA.nodes[idA].dualCnxs, networkA.nodes[idB].dualCnxs);
    /*
    cout << " A " << idA << endl;
    for (int i=0; i<networkA.nodes[idA].dualCnxs.n;++i){
        cout << networkA.nodes[idA].dualCnxs[i] << " ";
    }
    cout << endl;

    cout << " B " << idB << endl;
    for (int i=0; i<networkA.nodes[idB].dualCnxs.n;++i){
        cout << networkA.nodes[idB].dualCnxs[i] << " ";
    }
    cout << endl;
    */
    common.delValue(idDel);
    if(common.n==1) associated=common[0];
    else {//rare case with periodic interactions
        VecR<int> common1(0,common.n);
        for(int i=0; i<common.n; ++i){
            if(vContains(networkB.nodes[common[i]].netCnxs,idDel)) common1.addValue(common[i]);
        }
        if(common1.n==1) associated=common1[0];
        else{//rare case of large ring surrounding group
            //check a,b adjacent on ring
            VecR<int> common2(0,common1.n);
            for(int i=0; i<common1.n; ++i){
                VecR<int> ring=networkB.nodes[common[i]].dualCnxs;
                for(int j=0; j<ring.n; ++j){
                    int k=(j+1)%ring.n;
                    if(ring[j]==idA && ring[k]==idB){
                        common2.addValue(common1[i]);
                        break;
                    }
                    else if(ring[j]==idB && ring[k]==idA){
                        common2.addValue(common1[i]);
                        break;
                    }
                }
            }
            if(common2.n==1) associated=common2[0];
        }
    }

    if(associated==-1){
        wrapCoordinates();
        syncCoordinates();
        write("debug");
        cout<<idA<<" "<<idB<<" "<<idDel<<" "<<common<<endl;
        cout<<"ERROR IN ASSOCIATED NODE"<<endl;
    }

    return associated;
}

//Switch connectivities in lattice between 2x3 coordinate nodes
void LinkedNetwork::switchCnx33(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT) {
    bool globalVerbose=false;
    //unpck parameters
    int a,b,c,d,e,f;
    int u,v,w,x;
    int alpha, beta, gamma, delta, eta;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];

    alpha =     switchIdsT[6];
    beta =      switchIdsT[7];
    gamma =     switchIdsT[8];
    delta =     switchIdsT[9];
    eta =       switchIdsT[10];
    if (globalVerbose) cout << "1556" << endl;
    //Apply changes to descriptors due to breaking connections
    //For network A node distribution and edge distribution will remain unchanged

    //For network B node and edge distribution will change
    int nu, nv, nw, nx;
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) --networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) --networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }

    if (globalVerbose) cout << "1595" << endl;

    //A-A connectivities
    //swap a->e/d, b->d/e, d->b/a, e->a/b
    networkA.nodes[a].netCnxs.swapValue(e,d);
    networkA.nodes[b].netCnxs.swapValue(d,e);
    networkA.nodes[d].netCnxs.swapValue(b,a);
    networkA.nodes[e].netCnxs.swapValue(a,b);

    if (globalVerbose) cout << "1604" << endl;

    //A-B connectvities
    //swap a->v/x, b->u/w
    networkA.nodes[a].dualCnxs.swapValue(v,x);
    networkA.nodes[b].dualCnxs.swapValue(u,w);

    if (globalVerbose) cout << "1611" << endl;
    if (globalVerbose) cout << "u : " << u << " v : " << v << " x : " << x << " w : " << w << endl;
    //B-B connectivities
    //have to account for the fact that two nodes may connect multiple times
    //break insert w:u-(x)-v, x:u-(w)-v
    if (globalVerbose) {
        cout << "u cnxs" << endl;
        for (int i = 0; i < networkB.nodes[u].netCnxs.n; ++i) cout << networkB.nodes[u].netCnxs[i] << " ";
        cout << endl;
        cout << "v cnxs" << endl;
        for (int i=0;i<networkB.nodes[v].netCnxs.n;++i) cout << networkB.nodes[v].netCnxs[i] << " ";
        cout << endl;
        cout << "w cnxs" << endl;
        for (int i=0;i<networkB.nodes[w].netCnxs.n;++i) cout << networkB.nodes[w].netCnxs[i] << " ";
        cout << endl;
        cout << "x cnxs" << endl;
        for (int i=0;i<networkB.nodes[x].netCnxs.n;++i) cout << networkB.nodes[x].netCnxs[i] << " ";
        cout << endl;

    }

    networkB.nodes[u].netCnxs.swapValue(v,-1,w,x);
    networkB.nodes[v].netCnxs.swapValue(u,-1,w,x);
    networkB.nodes[u].netCnxs.delValue(-1);
    networkB.nodes[v].netCnxs.delValue(-1);
    if (globalVerbose) cout << "1620" << endl;
    if(vAdjCount(networkB.nodes[w].netCnxs,u,v)>1){
        VecR<int> eCnxs=networkA.nodes[e].dualCnxs;
        eCnxs.delValue(v);
        eCnxs.delValue(w);
        int ww=eCnxs[0];
        networkB.nodes[w].netCnxs.swapValue(v,-1,ww,u);
        networkB.nodes[w].netCnxs.insertValue(x,-1,u);
        networkB.nodes[w].netCnxs.swapValue(-1,v);
    }
    else{
/*
        cout << x << " between " << u << " and " << v << endl;
        for (int i=0;i<networkB.nodes[w].netCnxs.n;++i){
            cout << networkB.nodes[w].netCnxs[i] << " ";
        }
        cout << endl;
*/
        if (globalVerbose) cout << "1644" << endl;
        networkB.nodes[w].netCnxs.insertValue(x,u,v);
    }
    if(vAdjCount(networkB.nodes[x].netCnxs,u,v)>1){
        VecR<int> fCnxs=networkA.nodes[f].dualCnxs;
        fCnxs.delValue(v);
        fCnxs.delValue(x);
        int xx=fCnxs[0];
        networkB.nodes[x].netCnxs.swapValue(v,-1,xx,u);
        networkB.nodes[x].netCnxs.insertValue(w,-1,u);
        networkB.nodes[x].netCnxs.swapValue(-1,v);
    }
    else {
        if (globalVerbose) cout << "1657" << endl;
        if (globalVerbose) cout << "Insert " << w << " between " << u << " " << v << endl;


        networkB.nodes[x].netCnxs.insertValue(w, u, v);
    }
    if (globalVerbose) cout << "1650" << endl;

    //B-A connectivities
    //break u->b, v->a, insert w:a-(b)-e, z:b-(a)-d
    networkB.nodes[u].dualCnxs.delValue(b);
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b,a,e);
    networkB.nodes[x].dualCnxs.insertValue(a,b,d);

    if (globalVerbose) cout << "1659" << endl;

    //Apply changes to descriptors due to making connections
    //Network B
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) ++networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) ++networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    cout << "                                                                       Triangle raft switch" << endl;

    /*
    alpha=0, beta=0, gamma=0, delta=0, eta=0;
    for (int i=0; i<networkT.nodes[a].netCnxs.n;++i){
        for (int j=0; j<networkT.nodes[b].netCnxs.n;++j){
            if (networkT.nodes[a].netCnxs[i]==networkT.nodes[b].netCnxs[j]) {alpha=networkT.nodes[a].netCnxs[i];}
        }
        for (int j=0; j<networkT.nodes[c].netCnxs.n;++j){
            if (networkT.nodes[a].netCnxs[i]==networkT.nodes[c].netCnxs[j]) {beta=networkT.nodes[a].netCnxs[i];}
        }
        for (int j=0; j<networkT.nodes[e].netCnxs.n;++j){
            if (networkT.nodes[a].netCnxs[i]==networkT.nodes[e].netCnxs[j]) {gamma=networkT.nodes[a].netCnxs[i];}
        }
    }
    for (int i=0; i<networkT.nodes[b].netCnxs.n;++i){
        for (int j=0; j<networkT.nodes[d].netCnxs.n;++j){
            if (networkT.nodes[b].netCnxs[i]==networkT.nodes[d].netCnxs[j]) {delta=networkT.nodes[b].netCnxs[i];}
        }
        for (int j=0; j<networkT.nodes[f].netCnxs.n;++j){
            if (networkT.nodes[b].netCnxs[i]==networkT.nodes[f].netCnxs[j]) {eta=networkT.nodes[b].netCnxs[i];}
        }
    }
    */
    /*
    if (alpha==0) {cout << "alpha broken" << endl;}
    if (beta==0) {cout << "beta broken" << endl;}
    if (gamma==0) {cout << "gamma broken" << endl;}
    if (delta==0) {cout << "delta broken" << endl;}
    if (eta==0) {cout << "eta broken" << endl;}
    cout << alpha << " " << beta << " " << gamma << " " << delta << " " << eta << endl;
     */
    //T-T connectivities
    //swap a->e/d, b->d/e, d->b/a, e->a/b
    /*
    cout << "a  ";
    for (int i=0; i<networkT.nodes[a].netCnxs.n;++i){cout << networkT.nodes[a].netCnxs[i] << " ";}
    cout << endl;
    cout << "b  ";
    for (int i=0; i<networkT.nodes[b].netCnxs.n;++i){cout << networkT.nodes[b].netCnxs[i] << " ";}
    cout << endl;
    */
/*
    for (int i=0; i<networkT.nodes[a].netCnxs.n;++i){           if(networkT.nodes[a].netCnxs[i]==gamma){cout << "gamma in a" << endl;}    }
    for (int i=0; i<networkT.nodes[gamma].netCnxs.n;++i){       if(networkT.nodes[gamma].netCnxs[i]==a){cout << "a in gamma" << endl;}    }

    for (int i=0; i<networkT.nodes[b].netCnxs.n;++i){           if(networkT.nodes[b].netCnxs[i]==delta){cout << "delta in b" << endl;}    }
    for (int i=0; i<networkT.nodes[delta].netCnxs.n;++i){       if(networkT.nodes[delta].netCnxs[i]==b){cout << "b in delta" << endl;}    }
*/

    networkT.nodes[a].netCnxs.swapValue(gamma,delta);
    networkT.nodes[b].netCnxs.swapValue(delta,gamma);
    networkT.nodes[delta].netCnxs.swapValue(b,a);
    networkT.nodes[gamma].netCnxs.swapValue(a,b);

//    cout << "Si - O connections switch" << endl;
/*
    for (int i=0; i<networkT.nodes[gamma].netCnxs.n;++i){       if(networkT.nodes[gamma].netCnxs[i]==beta){cout << "beta in gamma" << endl;}    }
    for (int i=0; i<networkT.nodes[beta].netCnxs.n;++i){        if(networkT.nodes[beta].netCnxs[i]==gamma){cout << "gamma in beta" << endl;}    }

    for (int i=0; i<networkT.nodes[delta].netCnxs.n;++i){       if(networkT.nodes[delta].netCnxs[i]==eta){cout << "eta in delta" << endl;}    }
    for (int i=0; i<networkT.nodes[eta].netCnxs.n;++i){         if(networkT.nodes[eta].netCnxs[i]==delta){cout << "delta in eta" << endl;}    }
*/
    networkT.nodes[gamma].netCnxs.swapValue(beta, eta);
    networkT.nodes[beta].netCnxs.swapValue(gamma,delta);
    networkT.nodes[delta].netCnxs.swapValue(eta, beta);
    networkT.nodes[eta].netCnxs.swapValue(delta, gamma);
//    cout << "O - O connections switch" << endl;
    //T-B connectvities
    //swap a->v/x, b->u/w
/*
    cout << u << " " << v << " " << w << " " << x << endl;
    cout << "a  ";
    for (int i=0; i<networkT.nodes[a].dualCnxs.n;++i){cout << networkT.nodes[a].dualCnxs[i] << " ";}
    cout << endl;
    cout << "b  ";
    for (int i=0; i<networkT.nodes[b].dualCnxs.n;++i){cout << networkT.nodes[b].dualCnxs[i] << " ";}
    cout << endl;
*/
    networkT.nodes[a].dualCnxs.swapValue(v,x);
    networkT.nodes[b].dualCnxs.swapValue(u,w);

//    cout << "                                                                       Finished" << endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}

//Switch connectivities in lattice between 2x4 coordinate nodes
void LinkedNetwork::switchCnx44(VecF<int> switchIdsA, VecF<int> switchIdsB) {

    //unpck parameters
    int a,b,c,d,e,f,g,h;
    int u,v,w,x,y,z;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    g=switchIdsA[6];
    h=switchIdsA[7];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];
    y=switchIdsB[4];
    z=switchIdsB[5];

    //Apply changes to descriptors due to breaking connections
    //For network A node distribution and edge distribution will remain unchanged

    //For network B node and edge distribution will change
    int nu, nv, nw, nx;
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) --networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) --networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }

    //A-A connectivities
    //break a->e, b->d, insert a:b-(d)-c b:a-(e)-f, swap d->b/a, e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.delValue(d);
    networkA.nodes[a].netCnxs.insertValue(d,b,c);
    networkA.nodes[b].netCnxs.insertValue(e,a,f);
    networkA.nodes[d].netCnxs.swapValue(b,a);
    networkA.nodes[e].netCnxs.swapValue(a,b);

    //A-B connectvities
    //swap a->v/x, b->u/w
    networkA.nodes[a].dualCnxs.swapValue(v,x);
    networkA.nodes[b].dualCnxs.swapValue(u,w);

    //B-B connectivities
    //have to account for the fact that two nodes may connect multiple times
    //break u:x-(v)-y, v:w-(u)-z, insert w:y-(x)-v, x:u-(w)-z
    networkB.nodes[u].netCnxs.swapValue(v,-1,x,y);
    networkB.nodes[v].netCnxs.swapValue(u,-1,w,z);
    networkB.nodes[u].netCnxs.delValue(-1);
    networkB.nodes[v].netCnxs.delValue(-1);
    networkB.nodes[w].netCnxs.insertValue(x,y,v);
    networkB.nodes[x].netCnxs.insertValue(w,u,z);

    //B-A connectivities
    //break u->b, v->a, insert w:a-(b)-e, x:b-(a)-d
    networkB.nodes[u].dualCnxs.delValue(b);
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b,a,e);
    networkB.nodes[x].dualCnxs.insertValue(a,b,d);

    //B-B secondary connectivities
    //break v-y, u-z, swap y->v/x, z->u/w, add x-y, w-z
    networkB.nodes[v].auxCnxs.delValue(y);
    networkB.nodes[u].auxCnxs.delValue(z);
    networkB.nodes[y].auxCnxs.swapValue(v,x);
    networkB.nodes[z].auxCnxs.swapValue(u,w);
    networkB.nodes[x].auxCnxs.addValue(y);
    networkB.nodes[w].auxCnxs.addValue(z);

    //Apply changes to descriptors due to making connections
    //Network B
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) ++networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) ++networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }
}

//Switch connectivities in lattice between 4 and 3 coordinate nodes
void LinkedNetwork::switchCnx43(VecF<int> switchIdsA, VecF<int> switchIdsB) {

    //unpck parameters
    int a,b,c,d,e,f,g;
    int u,v,w,x,y;
    a=switchIdsA[0];
    b=switchIdsA[1];
    c=switchIdsA[2];
    d=switchIdsA[3];
    e=switchIdsA[4];
    f=switchIdsA[5];
    g=switchIdsA[6];
    u=switchIdsB[0];
    v=switchIdsB[1];
    w=switchIdsB[2];
    x=switchIdsB[3];
    y=switchIdsB[4];

    //Apply changes to descriptors due to breaking connections
    //For network A node distribution will remain unchanged but edge distribution will change
    int na, nb, nd, ne;
    na=networkA.nodes[a].netCnxs.n;
    nb=networkA.nodes[b].netCnxs.n;
    nd=networkA.nodes[d].netCnxs.n;
    ne=networkA.nodes[e].netCnxs.n;
    --networkA.edgeDistribution[na][ne];
    --networkA.edgeDistribution[nb][nd];
    --networkA.edgeDistribution[ne][na];
    --networkA.edgeDistribution[nd][nb];

    //For network B node and edge distribution will change
    int nu, nv, nw, nx;
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) --networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) --networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) --networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }

    //A-A connectivities
    //break a->e, insert a:b-(d)-c, swap b->d/e, d->b/a, e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.swapValue(d,e);
    networkA.nodes[a].netCnxs.insertValue(d,b,c);
    networkA.nodes[d].netCnxs.swapValue(b,a);
    networkA.nodes[e].netCnxs.swapValue(a,b);

    //A-B connectvities
    //swap a->v/x, b->u/w
    networkA.nodes[a].dualCnxs.swapValue(v,x);
    networkA.nodes[b].dualCnxs.swapValue(u,w);

    //B-B connectivities
    //have to account for the fact that two nodes may connect multiple times
    //break u:x-(v)-y, v:w-(u)-x, insert w:y-(x)-v, x:u-(w)-v
    networkB.nodes[u].netCnxs.swapValue(v,-1,x,y);
    networkB.nodes[v].netCnxs.swapValue(u,-1,w,x);
    networkB.nodes[u].netCnxs.delValue(-1);
    networkB.nodes[v].netCnxs.delValue(-1);
    networkB.nodes[w].netCnxs.insertValue(x,y,v);
    networkB.nodes[x].netCnxs.insertValue(w,u,v);

    //B-A connectivities
    //break u->b, v->a, insert w:a-(b)-e, x:b-(a)-d
    networkB.nodes[u].dualCnxs.delValue(b);
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b,a,e);
    networkB.nodes[x].dualCnxs.insertValue(a,b,d);

    //B-B secondary connectivities
    //break v-y, swap y->v/x, add x-y
    networkB.nodes[v].auxCnxs.delValue(y);
    networkB.nodes[y].auxCnxs.swapValue(v,x);
    networkB.nodes[x].auxCnxs.addValue(y);

    //Apply changes to descriptors due to making connections
    //Network A
    ++networkA.edgeDistribution[na][nd];
    ++networkA.edgeDistribution[nb][ne];
    ++networkA.edgeDistribution[nd][na];
    ++networkA.edgeDistribution[ne][nb];
    //Network B
    nu=networkB.nodes[u].netCnxs.n;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for(int i=0; i<nu; ++i){
        int id=networkB.nodes[u].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if(id!=v && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nu]; //prevent double counting
    }
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if(id!=u && id!=w && id!=x) ++networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if(id!=u && id!=v && id!=x) ++networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if(id!=u && id!=v && id!=w) ++networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }
}

//Mix connectivities to exchange 3/4 coordination nodes
void LinkedNetwork::mixCnx34(VecF<int> mixIdsA, VecF<int> mixIdsB) {

    //unpck parameters
    int a,b,c,d,e,f,g;
    int u,v,w,x,y;
    a=mixIdsA[0];
    b=mixIdsA[1];
    c=mixIdsA[2];
    d=mixIdsA[3];
    e=mixIdsA[4];
    f=mixIdsA[5];
    g=mixIdsA[6];
    u=mixIdsB[0];
    v=mixIdsB[1];
    w=mixIdsB[2];
    x=mixIdsB[3];
    y=mixIdsB[4];

    //Apply changes to descriptors due to breaking connections
    //For network A node distribution will remain unchanged but edge distribution will change
    int na, nb;
    na=networkA.nodes[a].netCnxs.n;
    nb=networkA.nodes[b].netCnxs.n;
    for(int i=0; i<na; ++i){
        int id=networkA.nodes[a].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[na][nCnx];
        if(id!=b) --networkA.edgeDistribution[nCnx][na]; //prevent double counting
    }
    for(int i=0; i<nb; ++i){
        int id=networkA.nodes[b].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[nb][nCnx];
        if(id!=a) --networkA.edgeDistribution[nCnx][nb]; //prevent double counting
    }

    //For network B node and edge distribution will change
    int nv,nw;
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if(id!=w) --networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if(id!=v) --networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }

    //A-A connectivities
    //break a->e, insert b:a-(e)-f, swap e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.insertValue(e,a,f);
    networkA.nodes[e].netCnxs.swapValue(a,b);

    //A-B connectvities
    //break a->v, insert b:u-(w)-v
    networkA.nodes[a].dualCnxs.delValue(v);
    networkA.nodes[b].dualCnxs.insertValue(w,u,v);

    //B-B connectivities
    //break v->u, insert w:y-(u)-v, swap u->v/w
    networkB.nodes[v].netCnxs.delValue(u);
    networkB.nodes[w].netCnxs.insertValue(u,y,v);
    networkB.nodes[u].netCnxs.swapValue(v,w);

    //B-A connectivities
    //break v->a, insert w:a-(b)-e
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b,a,e);

    //B-B secondary connectivities
    //break y-v, u-z, swap u->w/v, v->y/u, w->u/x, add x-w
    networkB.nodes[y].auxCnxs.delValue(v);
    networkB.nodes[u].auxCnxs.swapValue(w,v);
    networkB.nodes[w].auxCnxs.swapValue(u,x);
    networkB.nodes[v].auxCnxs.swapValue(y,u);
    networkB.nodes[x].auxCnxs.addValue(w);

    //Apply changes to descriptors due to making connections
    //Network A
    na=networkA.nodes[a].netCnxs.n;
    nb=networkA.nodes[b].netCnxs.n;
    for(int i=0; i<na; ++i){
        int id=networkA.nodes[a].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[na][nCnx];
        if(id!=b) ++networkA.edgeDistribution[nCnx][na]; //prevent double counting
    }
    for(int i=0; i<nb; ++i){
        int id=networkA.nodes[b].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[nb][nCnx];
        if(id!=a) ++networkA.edgeDistribution[nCnx][nb]; //prevent double counting
    }
    //Network B
    nv=networkB.nodes[v].netCnxs.n;
    nw=networkB.nodes[w].netCnxs.n;
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if(id!=w) ++networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nw; ++i){
        int id=networkB.nodes[w].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if(id!=v) ++networkB.edgeDistribution[nCnx][nw]; //prevent double counting
    }
}

//Mix connectivities of adjacent nodes to give -1/+1 in coordinations
void LinkedNetwork::mixCnx(VecF<int> mixIdsA, VecF<int> mixIdsB) {

    //unpack parameters
    /* a,b,c,d,e,f are nodes in lattice A
     * u,v,w,x,y,z are nodes in lattice B
     *  E      F            V
     *   \    /          /  |  \
     * --A---B--    XX--X   |   Z
     *  /     \         W   |   Y
     * C       D         \  |  /
     *                      U
     *
     *       E F               V
     *       |/              /  \
     * --A---B--        XX--X   Z
     *  /     \         W   |   Y
     * C       D         \  |  /
     *                      U
    */
    int a,b,c,d,e,f;
    int u,v,w,x,y,z,xx;
    a=mixIdsA[0];
    b=mixIdsA[1];
    c=mixIdsA[2];
    d=mixIdsA[3];
    e=mixIdsA[4];
    f=mixIdsA[5];
    u=mixIdsB[0];
    v=mixIdsB[1];
    w=mixIdsB[2];
    x=mixIdsB[3];
    y=mixIdsB[4];
    z=mixIdsB[5];
    xx=mixIdsB[6];


    //Apply changes to descriptors due to breaking connections
    //For network A node distribution and edge distribution will change as a result of a,b mixing
    int na, nb;
    na=networkA.nodes[a].netCnxs.n;
    nb=networkA.nodes[b].netCnxs.n;
    --networkA.nodeDistribution[na];
    --networkA.nodeDistribution[nb];
    for(int i=0; i<na; ++i){
        int id=networkA.nodes[a].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[na][nCnx];
        if(id!=b) --networkA.edgeDistribution[nCnx][na]; //prevent double counting
    }
    for(int i=0; i<nb; ++i){
        int id=networkA.nodes[b].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[nb][nCnx];
        if(id!=a) --networkA.edgeDistribution[nCnx][nb]; //prevent double counting
    }

    //For network B node and edge distribution will change as a result of v,x mixing
    int nv,nx;
    nv=networkB.nodes[v].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nx];
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if(id!=x) --networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if(id!=v) --networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }

    //A-A connectivities
    //break a->e, insert b:a-(e)-f, swap e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.insertValue(e,a,f);
    networkA.nodes[e].netCnxs.swapValue(a,b);

    //A-B connectvities
    //break a->v, insert b:u-(x)-v
    networkA.nodes[a].dualCnxs.delValue(v);
    networkA.nodes[b].dualCnxs.insertValue(x,u,v);

    //B-B connectivities
    //break v->u, insert x:xx-(u)-v, swap u->v/x
    networkB.nodes[v].netCnxs.swapValue(u,-1,x,z); //swap and delete as can occur multiple times if 2-cnd node
    networkB.nodes[v].netCnxs.delValue(-1);
    networkB.nodes[x].netCnxs.insertValue(u,xx,v);
    networkB.nodes[u].netCnxs.swapValue(v,x,w,y);

    //B-A connectivities
    //break v->a, insert x:a-(b)-e
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[x].dualCnxs.insertValue(b,a,e);

    //B-B secondary connectivities
    //DONT NEED FOR THIS GENERAL MIX MOVE

    //Apply changes to descriptors due to making connections
    //Network A
    na=networkA.nodes[a].netCnxs.n;
    nb=networkA.nodes[b].netCnxs.n;
    ++networkA.nodeDistribution[na];
    ++networkA.nodeDistribution[nb];
    for(int i=0; i<na; ++i){
        int id=networkA.nodes[a].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[na][nCnx];
        if(id!=b) ++networkA.edgeDistribution[nCnx][na]; //prevent double counting
    }
    for(int i=0; i<nb; ++i){
        int id=networkA.nodes[b].netCnxs[i];
        int nCnx=networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[nb][nCnx];
        if(id!=a) ++networkA.edgeDistribution[nCnx][nb]; //prevent double counting
    }
    //Network B
    nv=networkB.nodes[v].netCnxs.n;
    nx=networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nx];
    for(int i=0; i<nv; ++i){
        int id=networkB.nodes[v].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if(id!=x) ++networkB.edgeDistribution[nCnx][nv]; //prevent double counting
    }
    for(int i=0; i<nx; ++i){
        int id=networkB.nodes[x].netCnxs[i];
        int nCnx=networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if(id!=v) ++networkB.edgeDistribution[nCnx][nx]; //prevent double counting
    }
}

//Check mix move did not introduce any edges which form part of three rings
bool LinkedNetwork::checkThreeRingEdges(int id) {

    bool edgeCheck=true;
    int nCnxs=networkB.nodes[id].dualCnxs.n;
    for(int i=0; i<nCnxs; ++i){
        int j=networkB.nodes[id].dualCnxs[i];
        int k=networkB.nodes[id].dualCnxs[(i+1)%nCnxs];
        VecR<int> common=vCommonValues(networkA.nodes[j].dualCnxs,networkA.nodes[k].dualCnxs);
        if(common.n>2){
            edgeCheck=false;
            break;
        }
    }
    return edgeCheck;
}

//Rearrange nodes after connection switch to maintain convexity
bool LinkedNetwork::convexRearrangement(int cnxType, VecF<int> switchIdsA, VecF<int> switchIdsB) {

    bool convex;

    if(cnxType==33 || cnxType==34 || cnxType==44){
        //Unpack nodes
        int a,b,c,d,e,f;
        a=switchIdsA[0];
        b=switchIdsA[1];
        c=switchIdsA[2];
        d=switchIdsA[3];
        e=switchIdsA[4];
        f=switchIdsA[5];

        //Maintains which nodes are convex
        VecF<bool> convexNodes;
        if(cnxType==33) convexNodes=VecF<bool>(6);
        else if(cnxType==34) convexNodes=VecF<bool>(7);
        else if(cnxType==44) convexNodes=VecF<bool>(8);

        //Initial guess places a at the centre of cd, b at the centre of ef
        VecF<double> va(2),vb(2),vc(2),vd(2),ve(2),vf(2);
        va[0]=crds[2*a];
        va[1]=crds[2*a+1];
        vb[0]=crds[2*b];
        vb[1]=crds[2*b+1];
        vc[0]=crds[2*c];
        vc[1]=crds[2*c+1];
        vd[0]=crds[2*d];
        vd[1]=crds[2*d+1];
        ve[0]=crds[2*e];
        ve[1]=crds[2*e+1];
        vf[0]=crds[2*f];
        vf[1]=crds[2*f+1];
        VecF<double> vce(2),vdf(2),vcd(2),vef(2);
        vce=ve-vc;
        vdf=vf-vd;
        vcd=vd-vc;
        vef=vf-ve;
        vce[0]-=networkA.pb[0]*nearbyint(vce[0]*networkA.rpb[0]);
        vce[1]-=networkA.pb[1]*nearbyint(vce[1]*networkA.rpb[1]);
        vdf[0]-=networkA.pb[0]*nearbyint(vdf[0]*networkA.rpb[0]);
        vdf[1]-=networkA.pb[1]*nearbyint(vdf[1]*networkA.rpb[1]);
        vcd[0]-=networkA.pb[0]*nearbyint(vcd[0]*networkA.rpb[0]);
        vcd[1]-=networkA.pb[1]*nearbyint(vcd[1]*networkA.rpb[1]);
        vef[0]-=networkA.pb[0]*nearbyint(vef[0]*networkA.rpb[0]);
        vef[1]-=networkA.pb[1]*nearbyint(vef[1]*networkA.rpb[1]);
        va=vd-vcd/2.0+vdf/10.0;
        vb=ve+vef/2.0-vce/10.0;
        va[0]-=networkA.pb[0]*nearbyint(va[0]*networkA.rpb[0]);
        va[1]-=networkA.pb[1]*nearbyint(va[1]*networkA.rpb[1]);
        vb[0]-=networkA.pb[0]*nearbyint(vb[0]*networkA.rpb[0]);
        vb[1]-=networkA.pb[1]*nearbyint(vb[1]*networkA.rpb[1]);
        crds[2*a]=va[0];
        crds[2*a+1]=va[1];
        crds[2*b]=vb[0];
        crds[2*b+1]=vb[1];
        for(int i=0; i<switchIdsA.n; ++i) convexNodes[i]=checkConvexity(switchIdsA[i]);
        convex=(convexNodes==true);

        //Guess move a,b towards each other
        VecF<double> vab(2);
        if(!convex) {
            vab = (vb - va) * 0.5 * 0.1;
            vab[0] -= networkA.pb[0] * nearbyint(vab[0] * networkA.rpb[0]);
            vab[1] -= networkA.pb[1] * nearbyint(vab[1] * networkA.rpb[1]);
            for (int i = 0; i < 9; ++i) {
                va += vab;
                vb -= vab;
                va[0] -= networkA.pb[0] * nearbyint(va[0] * networkA.rpb[0]);
                va[1] -= networkA.pb[1] * nearbyint(va[1] * networkA.rpb[1]);
                vb[0] -= networkA.pb[0] * nearbyint(vb[0] * networkA.rpb[0]);
                vb[1] -= networkA.pb[1] * nearbyint(vb[1] * networkA.rpb[1]);
                crds[2*a]=va[0];
                crds[2*a+1]=va[1];
                crds[2*b]=vb[0];
                crds[2*b+1]=vb[1];
                for (int i = 0; i < switchIdsA.n; ++i) convexNodes[i] = checkConvexity(switchIdsA[i]);
                convex = (convexNodes == true);
                if (convex) break;
            }
        }
    }
    else throw(string("Not yet implemented!"));

    return convex;
}


void LinkedNetwork::makerFixed(){
    int count=0;
    rFixed.resetMaxSize(networkA.nodes.n);
    rFixed.setSize(networkA.nodes.n);
    VecF<double> crdFixed = networkB.nodes[Fixed_Ring[0]].crd;
    VecF<double> v(2);
    cout << "***********************" << endl;
    cout << networkB.nodes[Fixed_Ring[0]].crd[0] << " " << networkB.nodes[Fixed_Ring[0]].crd[1] << endl;

    for (int i=0; i<networkA.nodes.n; ++i){
        for (int j=0; j<2; j++){
            v[j] = networkA.nodes[i].crd[j]-crdFixed[j];
//            if (v[j]>networkA.pb[j]/2) v[j] -= networkA.pb[j]/2;
//            else if (v[j]<-networkA.pb[j]/2) v[j] += networkA.pb[j]/2;
        }
        double r = sqrt(pow(v[0],2)+pow(v[1],2));

        rFixed[i] = 1000/r;
        if (r<spiralRadius){
            count +=1;
//            cout << i <<" " << r << endl;
//            cout << networkA.nodes[i].crd[0] << " " << networkA.nodes[i].crd[1] << endl;
        }
    }


    for (int i=0;i<rFixed.n;++i){
        weights.push_back(rFixed[i]);
    }
    cout << "---------- " << count << " Spiral atoms" << endl;
    return;
}

//Single monte carlo switching move
VecF<int> LinkedNetwork::SpiralmonteCarloSwitchMoveLAMMPS(int a, double& SimpleGrapheneEnergy, double& TersoffGrapheneEnergy, double& TriangleRaftEnergy, double& BilayerEnergy, double& BNEnergy, int Selected) {
    bool globalVerbose=false;
    if (globalVerbose) cout << "Start Spiral MonteCarlo" << endl;
    Selected=0;
    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */

//    cout << " CONFIRM THAT SIMPLE GRAPHENE ENERGIES FROM LAMMPS AND MOLSIM ALIGN" << endl;
//    if (SimpleGraphene.GlobalPotentialEnergy() - globalPotentialEnergy(false, false, networkA)!=0) {
//        cout << "LAMMPS : " << SimpleGraphene.GlobalPotentialEnergy() << " vs opt.tpp " << globalPotentialEnergy(false, false, networkA) << endl;
//    }



    //Select valid random connection - that will not violate connection limits
    int b,u,v;
    VecF<int> switchIdsA, switchIdsB, switchIdsT;
    int validMove;
    int cnxType;

    for (int i=0;i<10;++i){
        //cout << "Failed pickSpiral" << endl;
        cnxType = pickSpiralCnx34(a,b,u,v,mtGen);

        if (globalVerbose) cout << "2459" << endl;

        validMove=generateSwitchIds34(cnxType,switchIdsA,switchIdsB,switchIdsT,a,b,u,v);
        if (validMove==0) break;
    }


    if (globalVerbose) cout << "2463 " << validMove << endl;
    if(validMove==1) {
        cout << "Cannot find any valid switch moves" << endl;
        //throw string("Cannot find any valid switch moves");
        VecF<int> status(3);
        status[0]=0;
        status[1]=0;
        status[2]=0;
        cout << "returning ..." <<endl;
        return status;
    }
    if (globalVerbose) cout << "2469" << endl;
    //Save current state
    //double saveEnergy=energy;
    double saveEnergySimpleGraphene=SimpleGrapheneEnergy;
    double saveEnergyTersoffGraphene=TersoffGrapheneEnergy;
    double saveEnergyTriangleRaft=TriangleRaftEnergy;
    double saveEnergyBilayer = BilayerEnergy;
    double saveEnergyBN = BNEnergy;
    if (globalVerbose) cout << "2476" << endl;
    if (isSimpleGraphene) {
        if (abs(SimpleGrapheneEnergy - SimpleGraphene.GlobalPotentialEnergy()) > 0.001) {
            cout << "SimpleGraphene" << endl;
            cout << "Saved : " << SimpleGrapheneEnergy << " vs Calculated " << SimpleGraphene.GlobalPotentialEnergy()
                 << endl;
        }
        else{
            VecF<double> gpe= globalPotentialEnergy(potParamsD[0], potParamsD[1], networkA);
            cout << gpe << endl;
            cout << "Simple Graphene : " << SimpleGraphene.GlobalPotentialEnergy() << "(" << gpe[2] <<")" << endl;

        }
    }
    if (isTriangleRaft) {
        if (abs(TriangleRaftEnergy - Triangle_Raft.GlobalPotentialEnergy()) > 0.001) {
            cout << "Triangle Raft" << endl;
            cout << "Saved : " << TriangleRaftEnergy << " vs Calculated " << Triangle_Raft.GlobalPotentialEnergy()
                 << endl;
        }
//        else{
//            cout << "Triangle Raft : " << Triangle_Raft.GlobalPotentialEnergy() << endl;
//        }
    }
    if (globalVerbose) cout << "2496" << endl;
    cout << crds[0] << endl;
    VecF<double> saveCrds=crds;
    cout << "Saved crds" << endl;
    double* saveCrdsSimpleGraphene ;
    double* saveCrdsTersoffGraphene;
    double* saveCrdsTriangleRaft   ;
    double* saveCrdsBilayer        ;
    double* saveCrdsBN             ;
    if (globalVerbose) cout << "2503" << endl;

    if (isSimpleGraphene)   saveCrdsSimpleGraphene      =SimpleGraphene.fetchCrds(2);
    if (isTersoffGraphene)  saveCrdsTersoffGraphene    =TersoffGraphene.fetchCrds(2);
    if (isTriangleRaft)     saveCrdsTriangleRaft        =Triangle_Raft.fetchCrds(2);
    if (isBilayer)          saveCrdsBilayer            =Bilayer.fetchCrds(3);
    if (isBN)               saveCrdsBN                  = BN.fetchCrds(2);


    VecF<int> saveNodeDistA=networkA.nodeDistribution;
    VecF<int> saveNodeDistB=networkB.nodeDistribution;

    VecF< VecF<int> > saveEdgeDistA=networkA.edgeDistribution;
    VecF< VecF<int> > saveEdgeDistB=networkB.edgeDistribution;

    if (globalVerbose) cout << "2518" << endl;

    VecF<Node> saveNodesA(switchIdsA.n), saveNodesB(switchIdsB.n), saveNodesT(switchIdsT.n);
    for(int i=0; i<saveNodesA.n; ++i) saveNodesA[i]=networkA.nodes[switchIdsA[i]];
    for(int i=0; i<saveNodesB.n; ++i) saveNodesB[i]=networkB.nodes[switchIdsB[i]];
    for(int i=0; i<saveNodesT.n; ++i) saveNodesT[i]=networkT.nodes[switchIdsT[i]];


    //Switch and geometry optimise
    cout << "Switch" << endl;
    VecF<int> optStatus_networkA(2);
    VecF<int> optStatus_SimpleGraphene(2);
    VecF<int> optStatus_TersoffGraphene(2);
    VecF<int> optStatus_TriangleRaft(2);
    VecF<int> optStatus_Bilayer(2);
    VecF<int> optStatus_BN(2);
    if (globalVerbose) cout << "2542" << endl;

    if(cnxType==33) {
        if (globalVerbose) cout << "2545" << endl;
        // works for network version of lammps objects
        switchCnx33(switchIdsA,switchIdsB, switchIdsT);
        if (globalVerbose) cout << "2548" << endl;

        // works for lammps objects
        if (isSimpleGraphene)   SimpleGraphene.switchGraphene(switchIdsA, networkA);
        if (isTriangleRaft)     Triangle_Raft.switchTriangleRaft(switchIdsA,switchIdsB, switchIdsT, networkT);
        if (isBilayer)          Bilayer.switchBilayer(switchIdsA,switchIdsB, switchIdsT);
        if (isBN)               BN.switchGraphene(switchIdsA, networkA);

    }
    else if(cnxType==44) switchCnx44(switchIdsA,switchIdsB);
    else if(cnxType==43) switchCnx43(switchIdsA,switchIdsB);
    else {
        cout << "Not yet implemented (?)" << endl;
        //throw string("Not yet implemented!");
    }

    if (globalVerbose) cout << "2544" << endl;

    //Rearrange nodes after switch
    bool geometryOK=true;
    geometryOK=checkThreeRingEdges(u);
    if(geometryOK) geometryOK=checkThreeRingEdges(v);
    if(geometryOK) {
//        if (potParamsD[1] == 0) localGeometryOptimisation(a, b, 1, false, false);
        if (potParamsD[1] == 0){

//            optStatus= globalGeometryOptimisation(false, false, network);
//            energy=globalPotentialEnergy(potParamsD[0],potParamsD[1], network);
//            cout << "                               optimise 1 " << energy  << endl;
        }
        else {
            geometryOK = convexRearrangement(cnxType, switchIdsA, switchIdsB);
            for (int i = 0; i < switchIdsA.n; ++i) {
                geometryOK = checkConvexity(switchIdsA[i]);
                if (!geometryOK) break;
            }
        }
//        if(!geometryOK) cout<<"yyy"<<endl;
    }
    else{
//        cout<<a<<" "<<b<<" xxx"<<endl;
        optStatus_SimpleGraphene = VecF<int>(3);
    }
    if(!geometryOK) optStatus_SimpleGraphene[0]=4;

    if (globalVerbose) cout << "2573" << endl;

    //Geometry optimisation of local region
    //cout << "Geometry Optimise" << endl;
    if(geometryOK){

//        syncCoordinates();
        //cout << "..............optimising" << endl;
//        cout << "---- Simple Graphene Minimisation" <<endl;
        optStatus_SimpleGraphene = SimpleGraphene.GlobalPotentialMinimisation();

        //cout << "feed simple graphene coords to Triangle raft" << endl;

        double* localCrdsSimpleGraphene       ;
        double* localCrdsTersoff              ;
        double* localCrdsTriangleRaft         ;
        double* localCrdsBilayer              ;
        double* localCrdsBN                   ;


        if (isSimpleGraphene) localCrdsSimpleGraphene         =SimpleGraphene.fetchCrds(2);
        if (isTersoffGraphene) localCrdsTersoff                =TersoffGraphene.fetchCrds(2);
        if (isTriangleRaft) localCrdsTriangleRaft           =Triangle_Raft.fetchCrds(2);
        if (isBilayer) localCrdsBilayer                =Bilayer.fetchCrds(3);
        if (isBN)   localCrdsBN                         = BN.fetchCrds(2);

        int natoms, nSi, nO;
        if (isBilayer) {
            int natoms = Bilayer.natoms;
            int nSi = (int) (round(natoms / 3) + 0.5);
            int nO = natoms - nSi;
        }
        double x,y;
        if (isTriangleRaft){
            x = localCrdsTriangleRaft[2*switchIdsT[6]] + localCrdsTriangleRaft[2*switchIdsT[7]]+localCrdsTriangleRaft[2*switchIdsT[9]];
            y = localCrdsTriangleRaft[2*switchIdsT[6]+1] + localCrdsTriangleRaft[2*switchIdsT[7]+1]+localCrdsTriangleRaft[2*switchIdsT[9]+1];
            localCrdsTriangleRaft[2*switchIdsA[0]] = x/3;
            localCrdsTriangleRaft[2*switchIdsA[0]+1] = y/3;

            x = localCrdsTriangleRaft[2*switchIdsT[6]] + localCrdsTriangleRaft[2*switchIdsT[8]]+localCrdsTriangleRaft[2*switchIdsT[10]];
            y = localCrdsTriangleRaft[2*switchIdsT[6]+1] + localCrdsTriangleRaft[2*switchIdsT[8]+1]+localCrdsTriangleRaft[2*switchIdsT[10]+1];
            localCrdsTriangleRaft[2*switchIdsA[1]] = x/3;
            localCrdsTriangleRaft[2*switchIdsA[1]+1] = y/3;
        }




//            if (isTriangleRaft and isSimpleGraphene) {
//                localCrdsTriangleRaft[2 * switchIdsA[i]] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
//                localCrdsTriangleRaft[2 * switchIdsA[i] + 1] =
//                        localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
//            }


        //    if (isSimpleGraphene and isTersoffGraphene) {
        //        localCrdsTersoff[2 * switchIdsA[i]] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * CScaling;
        //        localCrdsTersoff[2 * switchIdsA[i] + 1] = localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * CScaling;
        //    }
        //    if (isSimpleGraphene and isBilayer) {
        //        localCrdsBilayer[3 * (switchIdsA[i] * 2)] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
        //        localCrdsBilayer[3 * (switchIdsA[i] * 2) + 1] =
        //                localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
        //        localCrdsBilayer[3 * (switchIdsA[i] * 2 + 1)] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
        //        localCrdsBilayer[3 * (switchIdsA[i] * 2 + 1) + 1] =
        //                localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
        //        localCrdsBilayer[3 * (switchIdsA[i] + nSi)] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
        //        localCrdsBilayer[3 * (switchIdsA[i] + nSi) + 1] =
        //                localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
        //    }


        if (isTriangleRaft) Triangle_Raft.pushCrds(2, localCrdsTriangleRaft);
        if (isTersoffGraphene) TersoffGraphene.pushCrds(2,localCrdsTersoff);
        if (isBilayer) Bilayer.pushCrds(3,localCrdsBilayer);
        if (isBN) BN.pushCrds(2, localCrdsBN);

        int nThreads=0;
        //omp_set_dynamic(0);
        //omp_set_num_threads(4);


        /*
        #pragma omp parallel for num_threads(4)
        {
            cout << "threads=" << omp_get_num_threads() << endl;
            for (int i = 0; i < 4; ++i) {

                if (i == 0) {
                    cout << endl << "---- Triangle Raft Minimisation (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
                } else if (i == 1) {
                    cout << endl << "---- Tersoff Graphene Minimisation (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
                } else if (i == 2) {
                    cout << endl << "---- Bilayer Minimisation  (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
                } else if (i == 3) {
                    cout << endl << "---- networkA Minimisation  (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_networkA = globalGeometryOptimisation(false, false, networkA);
                }
            }
        }
        */

//        cout << endl << "---- Triangle Raft Minimisation"<< endl << endl;
//        optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
//        cout << endl << "---- Tersoff Graphene Minimisation " << endl << endl;
//        optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
//        cout << endl << "---- Bilayer Minimisation " << endl << endl;
//        optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
//        cout << endl << "---- networkA Minimisation" << endl << endl;
//        optStatus_networkA = globalGeometryOptimisation(false, false, networkA);
//        optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
        if (isOpenMP) {
//            cout << "---- networkA Minimisation" << endl;
            optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
//            cout << "                                                                       start parallel" << endl;
#pragma omp parallel num_threads(3)
            {
//            cout << "Number of Threads : " << omp_get_num_threads() << endl;



                if (omp_get_thread_num() == 0 and isTriangleRaft) {
                    //cout << omp_get_thread_num();
//                    cout << "---- Triangle Raft Minimisation        (thread " << omp_get_thread_num() << ")" << endl;
                    optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
                } else if (omp_get_thread_num() == 1 and isTersoffGraphene) {
                    //cout << omp_get_thread_num();
//                    cout << "---- Tersoff Graphene Minimisation     (thread " << omp_get_thread_num() << ")" << endl;
                    optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
                } else if (omp_get_thread_num() == 2 and isBilayer) {
                    //cout << omp_get_thread_num();
                    //                    cout << "---- Bilayer Minimisation              (thread " << omp_get_thread_num() << ")" << endl;
                    optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
                }
//                else if (omp_get_thread_num()==3) {
//                    cout << "---- networkA Minimisation             (thread " << omp_get_thread_num() << ")" << endl;
//                    optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
//                }


            }

//            if (isTriangleRaft) cout << "---- Triangle Raft Minimisation        "<< endl;
//            if (isTersoffGraphene) cout << "---- Tersoff Graphene Minimisation     "<< endl;
//            if (isBilayer) cout << "---- Bilayer Minimisation              "<< endl;

        }
        else{
            if (isTriangleRaft) {
//                cout << "---- Triangle Raft Minimisation" << endl;
                optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
            }
            if (isTersoffGraphene) {
//                cout << "---- Tersoff Graphene Minimisation " << endl;
                optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
            }
            if (isBilayer) {
//                cout << "---- Bilayer Minimisation " << endl;
                optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
            }
            if (isBN) {
                optStatus_BN = BN.GlobalPotentialEnergy();
            }
//            cout << "---- networkA Minimisation" <<  endl;
//            optStatus_networkA = globalGeometryOptimisation(false, false, networkA);
            optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
        }

//        optStatus=localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);

//        SimpleGrapheneEnergy=globalPotentialEnergy(potParamsD[0],potParamsD[1], network);
//        networkAEnergy      = globalPotentialEnergy(potParamsD[0], potParamsD[1], networkA);
        if (isSimpleGraphene) SimpleGrapheneEnergy = SimpleGraphene.GlobalPotentialEnergy();
        if (isTriangleRaft) TriangleRaftEnergy = Triangle_Raft.GlobalPotentialEnergy();
        if (isTersoffGraphene) TersoffGrapheneEnergy = TersoffGraphene.GlobalPotentialEnergy();
        if (isBilayer) BilayerEnergy = Bilayer.GlobalPotentialEnergy();
        if (isBN) BNEnergy = BN.GlobalPotentialEnergy();
//        sleep(100);
        syncCoordinates();
//        write("./output_files/debug");
    }
    else {
        SimpleGrapheneEnergy=numeric_limits<double>::infinity();
        TriangleRaftEnergy=numeric_limits<double>::infinity();
        TersoffGrapheneEnergy=numeric_limits<double>::infinity();
        BilayerEnergy=numeric_limits<double>::infinity();
        BNEnergy = numeric_limits<double>::infinity();
    }



    //Accept or reject
    //int accept=mc.acceptanceCriterion(energy);
    int accept;
    if (MC_Routine==1)        accept=mc.acceptanceCriterion(SimpleGrapheneEnergy, saveEnergySimpleGraphene,1.0);
//    else if (Selected==1)   accept=mc.acceptanceCriterion(TersoffGrapheneEnergy);
    else if (MC_Routine==2)   accept=mc.acceptanceCriterion(TriangleRaftEnergy, saveEnergyTriangleRaft,7.3448);
//    else if (Selected==3)   accept=mc.acceptanceCriterion(BilayerEnergy);
    else if (MC_Routine==5)   accept=mc.acceptanceCriterion(BNEnergy, saveEnergyBN, 7.0);

    //int accept=mc.acceptanceCriterion(SimpleGrapheneEnergy);

    if(accept==0){
        cout << "               Rejected MC Move        Ei = ";
        if (MC_Routine==1)      cout << saveEnergySimpleGraphene << " Ef = " << SimpleGrapheneEnergy << endl;
        else if (MC_Routine==2) cout << saveEnergyTriangleRaft   << " Ef = " << TriangleRaftEnergy   << endl;
        else if (MC_Routine==5) cout << saveEnergyBN << " Ef = " << BNEnergy << endl;



        //Triangle_Raft.revertTriangleRaft(switchIdsA, switchIdsB, switchIdsT);
        crds=saveCrds;
        networkA.nodeDistribution=saveNodeDistA;
        networkA.edgeDistribution=saveEdgeDistA;
        networkB.nodeDistribution=saveNodeDistB;
        networkB.edgeDistribution=saveEdgeDistB;

        for(int i=0; i<saveNodesA.n; ++i) networkA.nodes[switchIdsA[i]]=saveNodesA[i];
        for(int i=0; i<saveNodesB.n; ++i) networkB.nodes[switchIdsB[i]]=saveNodesB[i];
        for(int i=0; i<saveNodesT.n; ++i) networkT.nodes[switchIdsT[i]]=saveNodesT[i];

//        if (isSimpleGraphene) SimpleGrapheneEnergy=saveEnergySimpleGraphene;
//        if (isTriangleRaft) TriangleRaftEnergy  =saveEnergyTriangleRaft;
//        if (isTersoffGraphene) TersoffGrapheneEnergy = saveEnergyTersoffGraphene;
//        if (isBilayer) BilayerEnergy = saveEnergyBilayer;

//        crds=saveCrdsSimpleGraphene;


        if (isSimpleGraphene)       SimpleGraphene.pushCrds(2, saveCrdsSimpleGraphene);
        if (isTriangleRaft)         Triangle_Raft.pushCrds(2, saveCrdsTriangleRaft);
        if (isBilayer)              Bilayer.pushCrds(3, saveCrdsBilayer);
        if (isTersoffGraphene)      TersoffGraphene.pushCrds(2, saveCrdsTersoffGraphene);
        if (isBN)                   BN.pushCrds(2,saveCrdsBN);
        if (isSimpleGraphene)       {
            SimpleGraphene.revertGraphene(switchIdsA, networkA);
            SimpleGraphene.GlobalPotentialMinimisation();
        }
        if (isTriangleRaft)         {
            Triangle_Raft.revertTriangleRaft(switchIdsA, switchIdsB, switchIdsT);
            Triangle_Raft.GlobalPotentialMinimisation();
        }
        if (isBilayer)              {
            Bilayer.revertBilayer(switchIdsA, switchIdsB, switchIdsT);
            Bilayer.GlobalPotentialMinimisation();
        }
        if (isBN)                   {
            BN.revertGraphene(switchIdsA, networkA);
            BN.GlobalPotentialMinimisation();
        }



        if (isSimpleGraphene) SimpleGrapheneEnergy=SimpleGraphene.GlobalPotentialEnergy();
        if (isTriangleRaft) TriangleRaftEnergy  = Triangle_Raft.GlobalPotentialEnergy();
        if (isTersoffGraphene) TersoffGrapheneEnergy = TersoffGraphene.GlobalPotentialEnergy();
        if (isBilayer) BilayerEnergy = Bilayer.GlobalPotentialEnergy();
        if (isBN) BNEnergy = BN.GlobalPotentialEnergy();

    }
    else{
/*
        cout << "               Accepted MC Move        Ei = " << saveEnergy << " Ef = " << energy << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(crds[i]-saveCrds[i],2);
        }
        cout << "Change in coordinates on acceptance (new vs saved) : " << sum << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(crds[i]-networkT.nodes[(i-i%2)/2].crd[i%2],2);
        }
        cout << "Change in coordinates on acceptance (new vs networkT) : " << sum << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(saveCrds[i]-networkT.nodes[(i-i%2)/2].crd[i%2],2);
            //cout << "Saved : " << saveCrds[i] << " vs " << (i-i%2)/2 << " " << i%2 << " -> " << network.nodes[(i-i%2)/2].crd[i%2] << endl;
        }
        cout << "Change in coordinates on acceptance (saved vs networkT) : " << sum << endl;

*/
        cout << "               Accepted MC Move        Ei = ";
        if (MC_Routine==1)      cout << saveEnergySimpleGraphene << " Ef = " << SimpleGrapheneEnergy << endl;
        else if (MC_Routine==2) cout << saveEnergyTriangleRaft   << " Ef = " << TriangleRaftEnergy   << endl;
        else if (MC_Routine==5) cout << saveEnergyBN << " Ef = " << BNEnergy << endl;
        else cout << "MC ROUTINE : " << MC_Routine << endl;

        syncCoordinates();
        //writeXYZ("intermediate");
    }

    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force) 2=unsuccessful(it limit) 3=unsuccessful(intersection) 4=unsuccessful(non-convex)
     * [2] optimisation iterations */
    VecF<int> status(3);
    status[0]=accept;
    status[1]=optStatus_SimpleGraphene[0];
    status[2]=optStatus_SimpleGraphene[1];

    return status;
}

//Single monte carlo switching move
VecF<int> LinkedNetwork::monteCarloSwitchMoveLAMMPS(double& SimpleGrapheneEnergy, double& TersoffGrapheneEnergy, double& TriangleRaftEnergy, double& BilayerEnergy, double&BNEnergy, int Selected) {

    Selected=0;
    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */

//    cout << " CONFIRM THAT SIMPLE GRAPHENE ENERGIES FROM LAMMPS AND MOLSIM ALIGN" << endl;
//    if (SimpleGraphene.GlobalPotentialEnergy() - globalPotentialEnergy(false, false, networkA)!=0) {
//        cout << "LAMMPS : " << SimpleGraphene.GlobalPotentialEnergy() << " vs opt.tpp " << globalPotentialEnergy(false, false, networkA) << endl;
//    }



    //Select valid random connection - that will not violate connection limits
    int a,b,u,v;
    VecF<int> switchIdsA, switchIdsB, switchIdsT;
    int validMove;
    int cnxType;
    cout << "Find Move" << endl;
    for(int i=0; i<networkA.nodes.n*networkA.nodes.n; ++i){//catch in case cannot find any valid moves
        if (MCWeighting=="Weighted" || MCWeighting=="weighted") cnxType= pickDiscreteCnx34(a,b,u,v,mtGen);
        else                                                    cnxType= pickRandomCnx34(a, b, u, v, mtGen);
        validMove=generateSwitchIds34(cnxType,switchIdsA,switchIdsB,switchIdsT,a,b,u,v);
        if(validMove==0) break;
    }
    if(validMove==1) {
        cout << "Cannot find any valid switch moves" << endl;
        throw string("Cannot find any valid switch moves");
    }
    else{
        //    cout << "Found a valid switch ! " << endl;
    }

    //Save current state
    cout << "Save energies" << endl;
    //double saveEnergy=energy;
    double saveEnergySimpleGraphene=SimpleGrapheneEnergy;
    double saveEnergyTersoffGraphene=TersoffGrapheneEnergy;
    double saveEnergyTriangleRaft=TriangleRaftEnergy;
    double saveEnergyBilayer = BilayerEnergy;
    double saveEnergyBN      = BNEnergy;

    if (isSimpleGraphene) {
        if (abs(SimpleGrapheneEnergy - SimpleGraphene.GlobalPotentialEnergy()) > 0.001) {
            cout << "SimpleGraphene" << endl;
            cout << "Saved : " << SimpleGrapheneEnergy << " vs Calculated " << SimpleGraphene.GlobalPotentialEnergy()
                 << endl;
        }
        else{
//            VecF<double> gpe= globalPotentialEnergy(potParamsD[0], potParamsD[1], networkA);
//            cout << gpe << endl;
//            cout << "Simple Graphene : " << SimpleGraphene.GlobalPotentialEnergy() << "(" << gpe[2] <<")" << endl;
        }
    }
    if (isTriangleRaft) {
        if (abs(TriangleRaftEnergy - Triangle_Raft.GlobalPotentialEnergy()) > 0.001) {
            cout << "Triangle Raft" << endl;
            cout << "Saved : " << TriangleRaftEnergy << " vs Calculated " << Triangle_Raft.GlobalPotentialEnergy()
                 << endl;
        }
        else{
            cout << "Triangle Raft : " << Triangle_Raft.GlobalPotentialEnergy() << endl;
        }
    }
    cout << "Save Crds" << endl;
    VecF<double> saveCrds=crds;
    double* saveCrdsSimpleGraphene ;
    double* saveCrdsTersoffGraphene;
    double* saveCrdsTriangleRaft   ;
    double* saveCrdsBilayer        ;
    double* saveCrdsBN             ;

    if (isSimpleGraphene)   saveCrdsSimpleGraphene      =SimpleGraphene.fetchCrds(2);
    if (isTersoffGraphene)  saveCrdsTersoffGraphene    =TersoffGraphene.fetchCrds(2);
    if (isTriangleRaft)     saveCrdsTriangleRaft        =Triangle_Raft.fetchCrds(2);
    if (isBilayer)          saveCrdsBilayer            =Bilayer.fetchCrds(3);
    if (isBN)               saveCrdsBN                  = BN.fetchCrds(2);



    VecF<int> saveNodeDistA=networkA.nodeDistribution;
    VecF<int> saveNodeDistB=networkB.nodeDistribution;

    VecF< VecF<int> > saveEdgeDistA=networkA.edgeDistribution;
    VecF< VecF<int> > saveEdgeDistB=networkB.edgeDistribution;

    VecF<Node> saveNodesA(switchIdsA.n), saveNodesB(switchIdsB.n), saveNodesT(switchIdsT.n);
    for(int i=0; i<saveNodesA.n; ++i) saveNodesA[i]=networkA.nodes[switchIdsA[i]];
    for(int i=0; i<saveNodesB.n; ++i) saveNodesB[i]=networkB.nodes[switchIdsB[i]];
    for(int i=0; i<saveNodesT.n; ++i) saveNodesT[i]=networkT.nodes[switchIdsT[i]];


    //Switch and geometry optimise
    cout << "Switch" << endl;
    VecF<int> optStatus_networkA(2);
    VecF<int> optStatus_SimpleGraphene(2);
    VecF<int> optStatus_TersoffGraphene(2);
    VecF<int> optStatus_TriangleRaft(2);
    VecF<int> optStatus_Bilayer(2);
    VecF<int> optStauts_BN(2);


    if(cnxType==33) {
        // works for network version of lammps objects
        switchCnx33(switchIdsA,switchIdsB, switchIdsT);
        // works for lammps objects

        if (isSimpleGraphene)   {
            cout << "switch Simple Graphene" << endl;
            SimpleGraphene.switchGraphene(switchIdsA, networkA);
        }
        if (isTriangleRaft)     Triangle_Raft.switchTriangleRaft(switchIdsA,switchIdsB, switchIdsT, networkT);
        if (isBilayer)          Bilayer.switchBilayer(switchIdsA,switchIdsB, switchIdsT);
        if (isBN)               {
            cout << "switch BN" << endl;
            BN.switchGraphene(switchIdsA, networkA);
        }
    }
    else if(cnxType==44) switchCnx44(switchIdsA,switchIdsB);
    else if(cnxType==43) switchCnx43(switchIdsA,switchIdsB);
    else throw string("Not yet implemented!");



    //Rearrange nodes after switch
    bool geometryOK=true;
    geometryOK=checkThreeRingEdges(u);
    if(geometryOK) geometryOK=checkThreeRingEdges(v);
    if(geometryOK) {
//        if (potParamsD[1] == 0) localGeometryOptimisation(a, b, 1, false, false);
        if (potParamsD[1] == 0){

//            optStatus= globalGeometryOptimisation(false, false, network);
//            energy=globalPotentialEnergy(potParamsD[0],potParamsD[1], network);
//            cout << "                               optimise 1 " << energy  << endl;
        }
        else {
            geometryOK = convexRearrangement(cnxType, switchIdsA, switchIdsB);
            for (int i = 0; i < switchIdsA.n; ++i) {
                geometryOK = checkConvexity(switchIdsA[i]);
                if (!geometryOK) break;
            }
        }
//        if(!geometryOK) cout<<"yyy"<<endl;
    }
    else{
//        cout<<a<<" "<<b<<" xxx"<<endl;
        optStatus_SimpleGraphene = VecF<int>(3);
    }
    if(!geometryOK) optStatus_SimpleGraphene[0]=4;

    //Geometry optimisation of local region
    cout << "Geometry Optimise" << endl;
    if(geometryOK){

//        syncCoordinates();
        //cout << "..............optimising" << endl;
//        cout << "---- Simple Graphene Minimisation" <<endl;
        optStatus_SimpleGraphene = SimpleGraphene.GlobalPotentialMinimisation();

        //cout << "feed simple graphene coords to Triangle raft" << endl;

        double* localCrdsSimpleGraphene       ;
        double* localCrdsTersoff              ;
        double* localCrdsTriangleRaft         ;
        double* localCrdsBilayer              ;
        double* localCrdsBN                   ;

        if (isSimpleGraphene) localCrdsSimpleGraphene         =SimpleGraphene.fetchCrds(2);
        if (isTersoffGraphene) localCrdsTersoff                =TersoffGraphene.fetchCrds(2);
        if (isTriangleRaft) localCrdsTriangleRaft           =Triangle_Raft.fetchCrds(2);
        if (isBilayer) localCrdsBilayer                =Bilayer.fetchCrds(3);
        if (isBN) localCrdsBN = BN.fetchCrds(2);

        int natoms, nSi, nO;
        if (isBilayer) {
            int natoms = Bilayer.natoms;
            int nSi = (int) (round(natoms / 3) + 0.5);
            int nO = natoms - nSi;
        }

        for (int i=0; i<6;++i){
//            if (isTriangleRaft and isSimpleGraphene) {
//                localCrdsTriangleRaft[2 * switchIdsA[i]] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
//                localCrdsTriangleRaft[2 * switchIdsA[i] + 1] =
//                        localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
//            }


            if (isSimpleGraphene and isTersoffGraphene) {
                localCrdsTersoff[2 * switchIdsA[i]] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * CScaling;
                localCrdsTersoff[2 * switchIdsA[i] + 1] = localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * CScaling;
            }
            if (isSimpleGraphene and isBilayer) {
                localCrdsBilayer[3 * (switchIdsA[i] * 2)] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] * 2) + 1] =
                        localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] * 2 + 1)] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] * 2 + 1) + 1] =
                        localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] + nSi)] = localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] + nSi) + 1] =
                        localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
            }
        }

//        if (isTriangleRaft) Triangle_Raft.pushCrds(2, localCrdsTriangleRaft);
        if (isTersoffGraphene) TersoffGraphene.pushCrds(2,localCrdsTersoff);
        if (isBilayer) Bilayer.pushCrds(3,localCrdsBilayer);

        int nThreads=0;
        //omp_set_dynamic(0);
        //omp_set_num_threads(4);


        /*
        #pragma omp parallel for num_threads(4)
        {
            cout << "threads=" << omp_get_num_threads() << endl;
            for (int i = 0; i < 4; ++i) {

                if (i == 0) {
                    cout << endl << "---- Triangle Raft Minimisation (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
                } else if (i == 1) {
                    cout << endl << "---- Tersoff Graphene Minimisation (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
                } else if (i == 2) {
                    cout << endl << "---- Bilayer Minimisation  (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
                } else if (i == 3) {
                    cout << endl << "---- networkA Minimisation  (thread " << omp_get_thread_num() << ")" << endl
                         << endl;
                    optStatus_networkA = globalGeometryOptimisation(false, false, networkA);
                }
            }
        }
        */

//        cout << endl << "---- Triangle Raft Minimisation"<< endl << endl;
//        optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
//        cout << endl << "---- Tersoff Graphene Minimisation " << endl << endl;
//        optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
//        cout << endl << "---- Bilayer Minimisation " << endl << endl;
//        optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
//        cout << endl << "---- networkA Minimisation" << endl << endl;
//        optStatus_networkA = globalGeometryOptimisation(false, false, networkA);
//        optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
        if (isOpenMP) {
//            cout << "---- networkA Minimisation" << endl;
            optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
//            cout << "                                                                       start parallel" << endl;
            #pragma omp parallel num_threads(3)
            {
//            cout << "Number of Threads : " << omp_get_num_threads() << endl;



                if (omp_get_thread_num() == 0 and isTriangleRaft) {
                    //cout << omp_get_thread_num();
//                    cout << "---- Triangle Raft Minimisation        (thread " << omp_get_thread_num() << ")" << endl;
                    optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
                } else if (omp_get_thread_num() == 1 and isTersoffGraphene) {
                    //cout << omp_get_thread_num();
//                    cout << "---- Tersoff Graphene Minimisation     (thread " << omp_get_thread_num() << ")" << endl;
                    optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
                } else if (omp_get_thread_num() == 2 and isBilayer) {
                    //cout << omp_get_thread_num();
                    //                    cout << "---- Bilayer Minimisation              (thread " << omp_get_thread_num() << ")" << endl;
                    optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
                }
//                else if (omp_get_thread_num()==3) {
//                    cout << "---- networkA Minimisation             (thread " << omp_get_thread_num() << ")" << endl;
//                    optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
//                }


            }

//            if (isTriangleRaft) cout << "---- Triangle Raft Minimisation        "<< endl;
//            if (isTersoffGraphene) cout << "---- Tersoff Graphene Minimisation     "<< endl;
//            if (isBilayer) cout << "---- Bilayer Minimisation              "<< endl;

        }
        else{

            if (isTriangleRaft) {
//                cout << "---- Triangle Raft Minimisation" << endl;
                optStatus_TriangleRaft = Triangle_Raft.GlobalPotentialMinimisation();
            }
            if (isTersoffGraphene) {
//                cout << "---- Tersoff Graphene Minimisation " << endl;
                optStatus_TersoffGraphene = TersoffGraphene.GlobalPotentialMinimisation();
            }
            if (isBilayer) {
//                cout << "---- Bilayer Minimisation " << endl;
                optStatus_Bilayer = Bilayer.GlobalPotentialMinimisation();
            }
            if (isBN){
                optStauts_BN = BN.GlobalPotentialMinimisation();
            }
//            cout << "---- networkA Minimisation" <<  endl;
//            optStatus_networkA = globalGeometryOptimisation(false, false, networkA);
            optStatus_networkA = localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
        }

//        optStatus=localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);

//        SimpleGrapheneEnergy=globalPotentialEnergy(potParamsD[0],potParamsD[1], network);
//        networkAEnergy      = globalPotentialEnergy(potParamsD[0], potParamsD[1], networkA);
        if (isSimpleGraphene) SimpleGrapheneEnergy = SimpleGraphene.GlobalPotentialEnergy();
        if (isTriangleRaft) TriangleRaftEnergy = Triangle_Raft.GlobalPotentialEnergy();
        if (isTersoffGraphene) TersoffGrapheneEnergy = TersoffGraphene.GlobalPotentialEnergy();
        if (isBilayer) BilayerEnergy = Bilayer.GlobalPotentialEnergy();
        if (isBN) BNEnergy = BN.GlobalPotentialEnergy();
//        sleep(100);
        syncCoordinates();
//        write("./output_files/debug");
    }
    else {
        SimpleGrapheneEnergy=numeric_limits<double>::infinity();
        TriangleRaftEnergy=numeric_limits<double>::infinity();
        TersoffGrapheneEnergy=numeric_limits<double>::infinity();
        BilayerEnergy=numeric_limits<double>::infinity();
        BNEnergy=numeric_limits<double>::infinity();
    }



    cout << "Accept or reject" << endl;
    //int accept=mc.acceptanceCriterion(energy);
    int accept;
    if (MC_Routine==1)        accept=mc.acceptanceCriterion(SimpleGrapheneEnergy, saveEnergySimpleGraphene,1.0);
//    else if (Selected==1)   accept=mc.acceptanceCriterion(TersoffGrapheneEnergy);
    else if (MC_Routine==2)   accept=mc.acceptanceCriterion(TriangleRaftEnergy, saveEnergyTriangleRaft,7.3448);
//    else if (Selected==3)   accept=mc.acceptanceCriterion(BilayerEnergy);
    else if (MC_Routine==5)   accept=mc.acceptanceCriterion(BNEnergy, saveEnergyBN, 7.0);
    //int accept=mc.acceptanceCriterion(SimpleGrapheneEnergy);


    if(accept==0){
        cout << "               Rejected MC Move        Ei = ";
        if (MC_Routine==1)      cout << saveEnergySimpleGraphene << " Ef = " << SimpleGrapheneEnergy << endl;
        else if (MC_Routine==2) cout << saveEnergyTriangleRaft   << " Ef = " << TriangleRaftEnergy   << endl;
        else if (MC_Routine==5) cout << saveEnergyBN << " Ef = " << BNEnergy << endl;


        //Triangle_Raft.revertTriangleRaft(switchIdsA, switchIdsB, switchIdsT);
        crds=saveCrds;
        networkA.nodeDistribution=saveNodeDistA;
        networkA.edgeDistribution=saveEdgeDistA;
        networkB.nodeDistribution=saveNodeDistB;
        networkB.edgeDistribution=saveEdgeDistB;

        for(int i=0; i<saveNodesA.n; ++i) networkA.nodes[switchIdsA[i]]=saveNodesA[i];
        for(int i=0; i<saveNodesB.n; ++i) networkB.nodes[switchIdsB[i]]=saveNodesB[i];
        for(int i=0; i<saveNodesT.n; ++i) networkT.nodes[switchIdsT[i]]=saveNodesT[i];

        if (isSimpleGraphene) SimpleGrapheneEnergy=saveEnergySimpleGraphene;
        if (isTriangleRaft) TriangleRaftEnergy  =saveEnergyTriangleRaft;
        if (isTersoffGraphene) TersoffGrapheneEnergy = saveEnergyTersoffGraphene;
        if (isBilayer) BilayerEnergy = saveEnergyBilayer;
        if (isBN) BNEnergy = saveEnergyBN;

//        crds=saveCrdsSimpleGraphene;


        if (isSimpleGraphene)       SimpleGraphene.pushCrds(2, saveCrdsSimpleGraphene);
        if (isTriangleRaft)         Triangle_Raft.pushCrds(2, saveCrdsTriangleRaft);
        if (isBilayer)              Bilayer.pushCrds(3, saveCrdsBilayer);
        if (isTersoffGraphene)      TersoffGraphene.pushCrds(2, saveCrdsTersoffGraphene);
        if (isBN)                   BN.pushCrds(2, saveCrdsBN);


        if (isSimpleGraphene)       {
                SimpleGraphene.revertGraphene(switchIdsA, networkA);
                SimpleGraphene.GlobalPotentialMinimisation();
            }
        if (isTriangleRaft)         {
                Triangle_Raft.revertTriangleRaft(switchIdsA, switchIdsB, switchIdsT);
                Triangle_Raft.GlobalPotentialMinimisation();
            }
        if (isBilayer)              {
                Bilayer.revertBilayer(switchIdsA, switchIdsB, switchIdsT);
                Bilayer.GlobalPotentialMinimisation();
            }
        if (isBN)                   {
            BN.revertGraphene(switchIdsA, networkA);
            BN.GlobalPotentialMinimisation();
        }

    }
    else{
/*
        cout << "               Accepted MC Move        Ei = " << saveEnergy << " Ef = " << energy << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(crds[i]-saveCrds[i],2);
        }
        cout << "Change in coordinates on acceptance (new vs saved) : " << sum << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(crds[i]-networkT.nodes[(i-i%2)/2].crd[i%2],2);
        }
        cout << "Change in coordinates on acceptance (new vs networkT) : " << sum << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(saveCrds[i]-networkT.nodes[(i-i%2)/2].crd[i%2],2);
            //cout << "Saved : " << saveCrds[i] << " vs " << (i-i%2)/2 << " " << i%2 << " -> " << network.nodes[(i-i%2)/2].crd[i%2] << endl;
        }
        cout << "Change in coordinates on acceptance (saved vs networkT) : " << sum << endl;

*/
        cout << "               Accepted MC Move        Ei = ";
        if (MC_Routine==1)      cout << saveEnergySimpleGraphene << " Ef = " << SimpleGrapheneEnergy << endl;
        else if (MC_Routine==2) cout << saveEnergyTriangleRaft   << " Ef = " << TriangleRaftEnergy   << endl;
        else if (MC_Routine==5) cout << saveEnergyBN << " Ef = " << BNEnergy << endl;
        else cout << "MC ROUTINE : " << MC_Routine << endl;
        syncCoordinates();
        cout << "sync" << endl;
        //writeXYZ("intermediate");
    }

    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force) 2=unsuccessful(it limit) 3=unsuccessful(intersection) 4=unsuccessful(non-convex)
     * [2] optimisation iterations */
    VecF<int> status(3);
    status[0]=accept;
    status[1]=optStatus_SimpleGraphene[0];
    status[2]=optStatus_SimpleGraphene[1];

    cout << "status" << endl;

    return status;
}

//Single monte carlo switching move
VecF<int> LinkedNetwork::monteCarloSwitchMove(Network network, double& energy) {

    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */

    double sum=0;

    //Select valid random connection - that will not violate connection limits
    int a,b,u,v;
    VecF<int> switchIdsA, switchIdsB, switchIdsT;
    int validMove;
    int cnxType;
    for(int i=0; i<networkA.nodes.n*networkA.nodes.n; ++i){//catch in case cannot find any valid moves
        cnxType= pickRandomCnx34(a, b, u, v, mtGen);
        validMove=generateSwitchIds34(cnxType,switchIdsA,switchIdsB,switchIdsT,a,b,u,v);
        if(validMove==0) break;
    }
    if(validMove==1) {
        cout << "Cannot find any valid switch moves" << endl;
        throw string("Cannot find any valid switch moves");
    }
    else{
    //    cout << "Found a valid switch ! " << endl;
    }

    //Save current state
    //cout << "Save currrent state" << endl;
    double saveEnergy=energy;
    //networkT = Network(networkA.nodes, networkA.pb, networkA.rpb, "t", networkA.maxNetCnxs, networkA.maxDualCnxs);
//    globalGeometryOptimisation(false,false,networkT);
//    double saveEnergyTR=energy;
    //cout << "Triangle Raft Energy : " << saveEnergyTR<<endl;
    VecF<double> saveCrds=crds;
    VecF<int> saveNodeDistA=networkA.nodeDistribution;
    VecF<int> saveNodeDistB=networkB.nodeDistribution;
    //VecF<int> saveNodeDistT=networkT.nodeDistribution;

    VecF< VecF<int> > saveEdgeDistA=networkA.edgeDistribution;
    VecF< VecF<int> > saveEdgeDistB=networkB.edgeDistribution;
    //VecF< VecF<int> > saveEdgeDistT=networkT.edgeDistribution;

    VecF<Node> saveNodesA(switchIdsA.n), saveNodesB(switchIdsB.n), saveNodesT(switchIdsT.n);
    for(int i=0; i<saveNodesA.n; ++i) saveNodesA[i]=networkA.nodes[switchIdsA[i]];
    for(int i=0; i<saveNodesB.n; ++i) saveNodesB[i]=networkB.nodes[switchIdsB[i]];
    for(int i=0; i<saveNodesT.n; ++i) saveNodesT[i]=networkT.nodes[switchIdsT[i]];


    //Switch and geometry optimise
    cout << "Switch" << endl;
    VecF<int> optStatus(3);
    if(cnxType==33) switchCnx33(switchIdsA,switchIdsB, switchIdsT);
    else if(cnxType==44) switchCnx44(switchIdsA,switchIdsB);
    else if(cnxType==43) switchCnx43(switchIdsA,switchIdsB);
    else throw string("Not yet implemented!");



    //Rearrange nodes after switch
    bool geometryOK=true;
    geometryOK=checkThreeRingEdges(u);
    if(geometryOK) geometryOK=checkThreeRingEdges(v);
    if(geometryOK) {
//        if (potParamsD[1] == 0) localGeometryOptimisation(a, b, 1, false, false);
        if (potParamsD[1] == 0){

//            optStatus= globalGeometryOptimisation(false, false, network);
//            energy=globalPotentialEnergy(potParamsD[0],potParamsD[1], network);
//            cout << "                               optimise 1 " << energy  << endl;
        }
        else {
            geometryOK = convexRearrangement(cnxType, switchIdsA, switchIdsB);
            for (int i = 0; i < switchIdsA.n; ++i) {
                geometryOK = checkConvexity(switchIdsA[i]);
                if (!geometryOK) break;
            }
        }
//        if(!geometryOK) cout<<"yyy"<<endl;
    }
    else{
//        cout<<a<<" "<<b<<" xxx"<<endl;
        optStatus = VecF<int>(3);
    }
    if(!geometryOK) optStatus[0]=4;

    //Geometry optimisation of local region
    cout << "Geometry Optimise" << endl;
    if(geometryOK){

        optStatus= globalGeometryOptimisation(false, false, network);
//        optStatus=localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],potParamsD[1]);
        energy=globalPotentialEnergy(potParamsD[0],potParamsD[1], network);
//        cout << "                               optimise 2 " << energy  << endl;

//        syncCoordinates(network);
//        write("./output_files/debug");
    }
    else energy=numeric_limits<double>::infinity();

    //Accept or reject
    int accept=mc.acceptanceCriterion(energy,saveEnergy,1.00);
    if(accept==0){
        cout << "               Rejected MC Move        Ei = " << saveEnergy << " Ef = " << energy << endl;
        energy=saveEnergy;
        crds=saveCrds;
        networkA.nodeDistribution=saveNodeDistA;
        networkA.edgeDistribution=saveEdgeDistA;
        networkB.nodeDistribution=saveNodeDistB;
        networkB.edgeDistribution=saveEdgeDistB;
        for(int i=0; i<saveNodesA.n; ++i) networkA.nodes[switchIdsA[i]]=saveNodesA[i];
        for(int i=0; i<saveNodesB.n; ++i) networkB.nodes[switchIdsB[i]]=saveNodesB[i];
        for(int i=0; i<saveNodesT.n; ++i) networkT.nodes[switchIdsT[i]]=saveNodesT[i];
    }
    else{
/*
        cout << "               Accepted MC Move        Ei = " << saveEnergy << " Ef = " << energy << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(crds[i]-saveCrds[i],2);
        }
        cout << "Change in coordinates on acceptance (new vs saved) : " << sum << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(crds[i]-networkT.nodes[(i-i%2)/2].crd[i%2],2);
        }
        cout << "Change in coordinates on acceptance (new vs networkT) : " << sum << endl;
        sum = 0;
        for (int i=0;i<crds.n;++i){
            sum += pow(saveCrds[i]-networkT.nodes[(i-i%2)/2].crd[i%2],2);
            //cout << "Saved : " << saveCrds[i] << " vs " << (i-i%2)/2 << " " << i%2 << " -> " << network.nodes[(i-i%2)/2].crd[i%2] << endl;
        }
        cout << "Change in coordinates on acceptance (saved vs networkT) : " << sum << endl;

*/



        syncCoordinates();
        //writeXYZ("intermediate");
    }

    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force) 2=unsuccessful(it limit) 3=unsuccessful(intersection) 4=unsuccessful(non-convex)
     * [2] optimisation iterations */
    VecF<int> status(3);
    status[0]=accept;
    status[1]=optStatus[0];
    status[2]=optStatus[1];

    return status;
}

//Single monte carlo mixing move
VecF<int> LinkedNetwork::monteCarloMixMove(double& energy) {

    /* Single MC mix move (exchange 3<->4 coordination)
     * 1) select random connection
     * 2) mix connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */

    //Select valid random connection - that will not violate connection limits
    int a,b,u,v;
    VecF<int> mixIdsA, mixIdsB;
    int validMove;
    int cnxType;
    for(int i=0; i<networkA.nodes.n*networkA.nodes.n; ++i){//catch in case cannot find any valid moves
//        cnxType=pickRandomCnx34(a, b, u, v, mtGen);
        cnxType=pickRandomCnx(a,b,u,v,mtGen);
//        validMove=generateMixIds34(cnxType,mixIdsA,mixIdsB,a,b,u,v);
        validMove=generateMixIds(cnxType,mixIdsA,mixIdsB,a,b,u,v);
        if(validMove==0) break;
    }
    if(validMove==1) throw string("Cannot find any valid switch moves");

    //Save current state
    double saveEnergy=energy;
    VecF<double> saveCrdsA=crds;
    VecF<int> saveNodeDistA=networkA.nodeDistribution;
    VecF<int> saveNodeDistB=networkB.nodeDistribution;
    VecF< VecF<int> > saveEdgeDistA=networkA.edgeDistribution;
    VecF< VecF<int> > saveEdgeDistB=networkB.edgeDistribution;
    VecF<Node> saveNodesA(mixIdsA.n), saveNodesB(mixIdsB.n);
    for(int i=0; i<saveNodesA.n; ++i) saveNodesA[i]=networkA.nodes[mixIdsA[i]];
    for(int i=0; i<saveNodesB.n; ++i) saveNodesB[i]=networkB.nodes[mixIdsB[i]];

    //Switch and geometry optimise
    bool geometryOK=true;
    VecF<int> optStatus;
//    mixCnx34(mixIdsA,mixIdsB);
    mixCnx(mixIdsA,mixIdsB);
    geometryOK= checkThreeRingEdges(u);
    if(geometryOK) geometryOK= checkThreeRingEdges(v);
    //Unrestricted local optimisation of switched atoms
    if(geometryOK) {
        optStatus = localGeometryOptimisation(a, b, 1, false, false); //bond switch atoms only
        if (potParamsD[1] == 1) {
            for (int i = 0; i < mixIdsA.n; ++i) {
                geometryOK = checkConvexity(mixIdsA[i]);
                if (!geometryOK) break;
            }
        }
    }
    else optStatus = VecF<int>(3);
    if(!geometryOK) optStatus[0]=4;

    //Restricted optimisation of local region
    if(geometryOK) {
        optStatus = localGeometryOptimisation(a, b, goptParamsA[1], potParamsD[0], potParamsD[1]); //wider area
        energy = globalPotentialEnergy(potParamsD[0],potParamsD[1], networkA);
    }
    else energy=numeric_limits<double>::infinity();

    //Accept or reject
    int accept=mc.acceptanceCriterion(energy, saveEnergy, 1.00);
    if(accept==0){
        energy=saveEnergy;
        crds=saveCrdsA;
        networkA.nodeDistribution=saveNodeDistA;
        networkA.edgeDistribution=saveEdgeDistA;
        networkB.nodeDistribution=saveNodeDistB;
        networkB.edgeDistribution=saveEdgeDistB;
        for(int i=0; i<saveNodesA.n; ++i) networkA.nodes[mixIdsA[i]]=saveNodesA[i];
        for(int i=0; i<saveNodesB.n; ++i) networkB.nodes[mixIdsB[i]]=saveNodesB[i];
    }

    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force) 2=unsuccessful(it limit) 3=unsuccessful(intersection)
     * [2] optimisation iterations */
    VecF<int> status(3);
    status[0]=accept;
    status[1]=optStatus[0];
    status[2]=optStatus[1];

    return status;
}

//Single monte carlo switching move based on cost function
VecF<int> LinkedNetwork::monteCarloCostSwitchMove(double &cost, double &energy, double pTarget, double rTarget) {
    //***** CHECK IF NEED TO INCLUDE CONVEX REARRANGEMENT****//
    cout<<"WARNING DEPRECATED - MAY NOT FUNCTION AS INTENDED"<<endl;

    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate cost
     * 4) accept or reject */

    //Select valid random connection - that will not violate connection limits
    int a,b,u,v;
    VecF<int> switchIdsA, switchIdsB, switchIdsT;
    int validMove;
    int cnxType;
    for(int i=0; i<networkA.nodes.n*networkA.nodes.n; ++i){//catch in case cannot find any valid moves
        cnxType= pickRandomCnx34(a, b, u, v, mtGen);
        validMove=generateSwitchIds34(cnxType,switchIdsA,switchIdsB,switchIdsT,a,b,u,v);
        if(validMove==0) break;
    }
    if(validMove==1) throw string("Cannot find any valid switch moves");

    //Save current state
    double saveCost=cost;
    double saveEnergy=energy;
    VecF<double> saveCrdsA=crds;
    VecF<int> saveNodeDistA=networkA.nodeDistribution;
    VecF<int> saveNodeDistB=networkB.nodeDistribution;
    VecF< VecF<int> > saveEdgeDistA=networkA.edgeDistribution;
    VecF< VecF<int> > saveEdgeDistB=networkB.edgeDistribution;
    VecF<Node> saveNodesA(switchIdsA.n), saveNodesB(switchIdsB.n);
    for(int i=0; i<saveNodesA.n; ++i) saveNodesA[i]=networkA.nodes[switchIdsA[i]];
    for(int i=0; i<saveNodesB.n; ++i) saveNodesB[i]=networkB.nodes[switchIdsB[i]];

    //Switch, evaluate cost and geometry optimise
    VecF<int> optStatus;
    if(cnxType==33) switchCnx33(switchIdsA,switchIdsB,switchIdsT);
    else if(cnxType==44) switchCnx44(switchIdsA,switchIdsB);
    else throw string("Not yet implemented!");
    cost=costFunction(pTarget,rTarget);
    //Unrestricted optimisation of switched atoms
    localGeometryOptimisation(a,b,1,false,false); //bond switch atoms only
    //Restricted optimistation of local region
    optStatus=localGeometryOptimisation(a,b,goptParamsA[1],potParamsD[0],true); //wider area

    //Make sure not infinite energy the accept or reject
    int accept;
    energy=globalPotentialEnergy(potParamsD[0],potParamsD[1], networkA);
    if(energy==numeric_limits<double>::infinity()) accept=0;
    if(accept!=0) accept=mcCost.acceptanceCriterion(cost,saveEnergy, 1.00);
    if(accept==0){
        energy=saveEnergy;
        crds=saveCrdsA;
        networkA.nodeDistribution=saveNodeDistA;
        networkA.edgeDistribution=saveEdgeDistA;
        networkB.nodeDistribution=saveNodeDistB;
        networkB.edgeDistribution=saveEdgeDistB;
        for(int i=0; i<saveNodesA.n; ++i) networkA.nodes[switchIdsA[i]]=saveNodesA[i];
        for(int i=0; i<saveNodesB.n; ++i) networkB.nodes[switchIdsB[i]]=saveNodesB[i];
    }

    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force) 2=unsuccessful(it limit) 3=unsuccessful(intersection)
     * [2] optimisation iterations */
    VecF<int> status(3);
    status[0]=accept;
    status[1]=optStatus[0];
    status[2]=optStatus[1];

    return status;
}

//Cost fuctional based on ring statistics and assortative mixing
double LinkedNetwork::costFunction(double &pTarget, double &rTarget) {

    double p6=getNodeDistribution("B")[6];
    double r=getAssortativity("B");
    double cost=costParams[0]*abs(p6-pTarget)+costParams[1]*abs(r-rTarget);
    return cost;
}


//Calculate potential energy of entire system
double LinkedNetwork::globalPotentialEnergy(bool useIntx, bool restrict, Network network) {

    /* Potential model
     * Bonds as harmonic
     * Angles as harmonic
     * Local line intersections */
    bool globalVerbose=false;
    if (globalVerbose) cout << "Number of nodes to minimise : " << network.nodes.n << endl;
    if (globalVerbose) cout << "Max NetworkA net cnxs : " << networkA.maxNetCnxs << endl;
    int maxCnxs = network.maxNetCnxs;

    if (globalVerbose) cout << "Max Coordination Number     : " << network.maxNetCnxs << endl;
    //Set up potential model
    //Bond and angle harmonics

//    VecR<int> bonds(0,network.nodes.n*maxCnxs*2),angles(0,network.nodes.n*maxCnxs*3);
    if (globalVerbose) cout << "create blank identities" << endl;
    if (globalVerbose) cout << network.nodes.n*6 << " " << network.nodes.n*maxCnxs*3<< endl;

    VecR<int> bonds(0,network.nodes.n*6),angles(0,network.nodes.n*maxCnxs*3);
//    VecR<double> bondParams(0,network.nodes.n*maxCnxs*3),angleParams(0,network.nodes.n*maxCnxs*3);
    if (globalVerbose) cout << "create blank params" << endl;
    VecR<double> bondParams(0,network.nodes.n*6),angleParams(0,network.nodes.n*maxCnxs*3) ;



    if (globalVerbose) cout << "created blanks" << endl;
    if (network.nodes.n>networkA.nodes.n) {
//        cout << "harmonics only" << endl;
        for (int i = 0; i < network.nodes.n; ++i) {
            generateHarmonicsOnly(i, bonds, bondParams, network);
        }
    }
    else{


        for(int i=0; i<network.nodes.n; ++i){
            generateHarmonics(i,bonds,bondParams,angles,angleParams, network);
        }
    }
    if (globalVerbose) cout << "generated harmonics" << endl;
    //Intersections
    VecR<int> intersections(0,network.nodes.n*1000);
    if(useIntx){
        for(int i=0; i<networkB.nodes.n; ++i){
            generateRingIntersections(i,intersections);
        }
    }

    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),angs(angles.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),angP(angleParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<angs.n; ++i) angs[i]=angles[i];
    for(int i=0; i<angP.n; ++i) angP[i]=angleParams[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];

    //Potential model based on geometry code
    double potEnergy;
    if(network.geometryCode=="2DE"){
        if(!restrict) {
            HI2DP potModel(network.pb[0], network.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx) potModel.setIntersections(intx, intxP);
            potEnergy = potModel.function(crds);
        }
        else {
            HRI2DP potModel(network.pb[0], network.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx) potModel.setIntersections(intx, intxP);
            potEnergy = potModel.function(crds);
        }
    }
    else if(network.geometryCode=="2DS"){
        VecF<int> constrained(network.nodes.n);
        for(int i=0; i<network.nodes.n; ++i) constrained[i]=i;
        HI3DS potModel;
        potModel.setBonds(bnds,bndP);
        potModel.setAngles(angs,angP);
        potModel.setGeomConstraints(constrained,potParamsC);
        if(useIntx) potModel.setIntersections(intx,intxP);
//        for(int i=0; i<intx.n/4; ++i){
//            cout<<intx[4*i]<<" "<<intx[4*i+1]<<" "<<intx[4*i+2]<<" "<<intx[4*i+3]<<endl;
//        }
        potEnergy=potModel.function(crds);
    }

    else if(network.geometryCode=="2DEtr"){
        HLJ2DP potModel(network.pb[0], network.pb[1]);
        potModel.setBonds(bnds, bndP);
        potModel.setAngles(angs, angP);
        if (useIntx) potModel.setIntersections(intx, intxP);
        potEnergy = potModel.function(crds);
    }
    //Convexity
    if(restrict) {
        bool convex = checkConvexity();
        if (!convex) potEnergy = numeric_limits<double>::infinity();
    }
    if (globalVerbose) cout << endl << endl << endl;
    if (globalVerbose) cout << bonds.n << "   " << angles.n << endl << endl << endl;

//    OutputFile bondfile("bonds.dat");
//    bondfile.initVariables(6,4,60,20);
//    for (int i=0; i<bonds.n/2;++i){
//        bondfile.write(to_string(bonds[2*i]+1)+"    "+to_string(bonds[2*i+1]+1));
//    }
//
//
//    OutputFile anglefile("angles.dat");
//    anglefile.initVariables(6,4,60,20);
//    for (int i=0; i<angles.n/3;++i){
//        anglefile.write(to_string(angles[3*i]+1)+"    "+to_string(angles[3*i+1]+1) + "    "+to_string(angles[3*i+2]+1));
//    }
//

    return potEnergy;
}
/*
//Calculate potential energy of entire system
double LinkedNetwork::globalPotentialEnergyTR(bool useIntx, bool restrict, Network network) {

    //* Potential model
    //* Bonds as harmonic
    //* Angles as harmonic
    //* Local line intersections

    //Set up potential model
    //Bond and angle harmonics
    VecR<int> bonds(0,network.nodes.n*maxACnxs*2),angles(0,network.nodes.n*maxACnxs*3);
    VecR<double> bondParams(0,network.nodes.n*maxACnxs*3),angleParams(0,network.nodes.n*maxACnxs*3);
    for(int i=0; i<network.nodes.n; ++i){
        generateHarmonics(i,bonds,bondParams,angles,angleParams, network);
    }
    //Intersections
    VecR<int> intersections(0,network.nodes.n*1000);
    if(useIntx){
        for(int i=0; i<networkB.nodes.n; ++i){
            generateRingIntersections(i,intersections);
        }
    }

    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),angs(angles.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),angP(angleParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<angs.n; ++i) angs[i]=angles[i];
    for(int i=0; i<angP.n; ++i) angP[i]=angleParams[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];

    //Potential model based on geometry code
    double potEnergy;
    if(network.geometryCode=="2DE"){
        if(!restrict) {
            HI2DP potModel(network.pb[0], network.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx) potModel.setIntersections(intx, intxP);
            potEnergy = potModel.function(crds);
        }
        else {
            HRI2DP potModel(network.pb[0], network.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx) potModel.setIntersections(intx, intxP);
            potEnergy = potModel.function(crds);
        }
    }
    else if(network.geometryCode=="2DS"){
        VecF<int> constrained(network.nodes.n);
        for(int i=0; i<network.nodes.n; ++i) constrained[i]=i;
        HI3DS potModel;
        potModel.setBonds(bnds,bndP);
        potModel.setAngles(angs,angP);
        potModel.setGeomConstraints(constrained,potParamsC);
        if(useIntx) potModel.setIntersections(intx,intxP);
//        for(int i=0; i<intx.n/4; ++i){
//            cout<<intx[4*i]<<" "<<intx[4*i+1]<<" "<<intx[4*i+2]<<" "<<intx[4*i+3]<<endl;
//        }
        potEnergy=potModel.function(crds);
    }

    else if(network.geometryCode=="2DEtr"){
        HLJ2DP potModel(network.pb[0], network.pb[1]);
        potModel.setBonds(bnds, bndP);
        potModel.setAngles(angs, angP);
        if (useIntx) potModel.setIntersections(intx, intxP);
        potEnergy = potModel.function(crds);
    }
    //Convexity
    if(restrict) {
        bool convex = checkConvexity();
        if (!convex) potEnergy = numeric_limits<double>::infinity();
    }

    return potEnergy;
}

*/

//Geometry optimise entire system
VecF<int> LinkedNetwork::globalGeometryOptimisation(bool useIntx, bool restrict, Network network) {
//    cout << "start_GEO" << endl;
    auto start_GEO = chrono::system_clock::now();

    /* Potential model
     * Bonds as harmonic
     * Angles as approximated harmonic
     * Local line intersections */

//    cout << "Number of nodes to minimise : " << network.nodes.n << endl;
    int maxCnxs = network.maxNetCnxs;
//    cout << "Max Coordination Number     : " << network.maxNetCnxs << endl;
    //Set up potential model
    //Bond and angle harmonics
    VecR<int> bonds(0,network.nodes.n*maxCnxs*2),angles(0,network.nodes.n*maxCnxs*3);
    VecR<double> bondParams(0,network.nodes.n*maxCnxs*3),angleParams(0,network.nodes.n*maxCnxs*3);
//    cout << "Made empty vectors" << endl;
    if (network.nodes.n>networkA.nodes.n){
        for(int i=0; i<network.nodes.n; ++i){
            generateHarmonicsOnly(i,bonds,bondParams,network);
        }
    }
    else{
        for(int i=0; i<network.nodes.n; ++i){
            generateHarmonics(i,bonds,bondParams,angles,angleParams, network);
        }
    }
//    cout << "Length of Bonds : " << bonds.n << endl;
//    cout << "Generated Harmonics" << endl;
    //Intersections

    VecR<int> intersections(0,network.nodes.n*1000);
    /*
    useIntx=false;
    if(useIntx){
        for(int i=0; i<networkB.nodes.n; ++i){
            generateRingIntersections(i,intersections);
        }
    }
    */

    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),angs(angles.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),angP(angleParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<angs.n; ++i) angs[i]=angles[i];
    for(int i=0; i<angP.n; ++i) angP[i]=angleParams[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];



    VecF<int> optStatus(2);

    //Geometry optimise
    if(network.geometryCode=="2DE"){
        if(!restrict) {
            HI2DP potModel(network.pb[0], network.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//only optimise if no line intersections
                SteepestDescentArmijoMultiDim<HI2DP> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;

            }
        }
        else{
            HRI2DP potModel(network.pb[0], network.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//only optimise if no line intersections
                SteepestDescentArmijoMultiDim<HRI2DP> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;

            }
        }
    }
    else if(network.geometryCode=="2DS"){
        VecF<int> constrained(network.nodes.n);
        for(int i=0; i<network.nodes.n; ++i) constrained[i]=i;
        if(!restrict) {
            HI3DS potModel;
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            potModel.setGeomConstraints(constrained, potParamsC);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//only optimise if no line intersections
                SteepestDescentArmijoMultiDim<HI3DS> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;

            }
        }
        else{
            HRI3DS potModel;
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            potModel.setGeomConstraints(constrained, potParamsC);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//only optimise if no line intersections
                SteepestDescentArmijoMultiDim<HRI3DS> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;

            }
        }
    }
    else if (network.geometryCode=="2DEtr"){
        cout << "Minimise Triangle Raft" << endl;
        HLJ2DP potModel(network.pb[0], network.pb[1]);
        potModel.setBonds(bnds, bndP);
//        cout << "        set bonds" << endl;
//        potModel.setRepulsions();
        if (potModel.function(crds) < numeric_limits<double>::infinity()) {//only optimise if no line intersections
//            cout << "        e < infinity" << endl;
            SteepestDescentArmijoMultiDim<HLJ2DP> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
            optStatus = optimiser(potModel, crds);



//            SteepestDescentMultiDim<HLJ2DP> optimiser2(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
//            optStatus = optimiser2(potModel, crds);
        } else {
            optStatus[0] = 3;
            optStatus[1] = 0;

        }
    }
    auto end_GEO = chrono::system_clock::now();
    chrono::duration<double> elapsed_GEO = end_GEO-start_GEO;
    cout << "GEO took " << elapsed_GEO.count() << endl;
    return optStatus;
}


/*
//Geometry optimise entire system
VecF<int> LinkedNetwork::globalGeometryOptimisationTR(bool useIntx, bool restrict) {

    // * Potential model
    // * Bonds as harmonic
    // * Angles as approximated harmonic
    // * Local line intersections

    //Set up potential model
    //Bond and angle harmonics
    VecR<int> bonds(0,networkT.nodes.n*maxACnxs*2),angles(0,networkT.nodes.n*maxACnxs*3);
    VecR<double> bondParams(0,networkT.nodes.n*maxACnxs*3),angleParams(0,networkT.nodes.n*maxACnxs*3);
    for(int i=0; i<networkT.nodes.n; ++i){
        generateHarmonicsTR(i,bonds,bondParams,angles,angleParams);
    }

    //Intersections
    VecR<int> intersections(0,networkT.nodes.n*1000);

    if(useIntx){
        for(int i=0; i<networkB.nodes.n; ++i){
            generateRingIntersections(i,intersections);
        }
    }


    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),angs(angles.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),angP(angleParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<angs.n; ++i) angs[i]=angles[i];
    for(int i=0; i<angP.n; ++i) angP[i]=angleParams[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];

    HLJ2DP potModel(networkT.pb[0], networkT.pb[1]);
    potModel.setBonds(bnds, bndP);
    if (potModel.function(crdsT) < numeric_limits<double>::infinity()) {//only optimise if no line intersections
        SteepestDescentArmijoMultiDim<HLJ2DP> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
        optimiser(potModel, crds);
    }

}

*/



//Geometry optimise subsection of system by only including interactions in a specified range
VecF<int> LinkedNetwork::localGeometryOptimisation(int centreA, int centreB, int extent, bool useIntx, bool restrict) {

    /* Find three regions
     * 1) local (full interactions)
     * 2) fixed inner (interactions with local and fixed but immobile)
     * 3) fixed outer (interactions with fixed inner but immobile) */
    VecR<int> local, fixedInner, fixedOuter;
    networkA.findLocalRegion(centreA,centreB,extent,local,fixedInner,fixedOuter);

//    /* Find unique dual nodes associated with nodes in local and fixed regions */
//    VecR<int> localDual(0,local.n*100);
//    for(int i=0; i<local.n; ++i){
//        int id=local[i];
//        for(int j=0; j<networkA.nodes[id].dualCnxs.n; ++j){
//            int rId=networkA.nodes[id].dualCnxs[j];
//            localDual.addValue(rId);
//        }
//    }
//    for(int i=0; i<fixedInner.n; ++i){
//        int id=fixedInner[i];
//        for(int j=0; j<networkA.nodes[id].dualCnxs.n; ++j){
//            int rId=networkA.nodes[id].dualCnxs[j];
//            localDual.addValue(rId);
//        }
//    }
//    localDual=vUnique(localDual);

    //Harmonics
    int reserveSize=(local.n+fixedInner.n)*maxACnxs*3;
    VecR<int> bonds(0,reserveSize),angles(0,reserveSize);
    VecR<double> bondParams(0,reserveSize),angleParams(0,reserveSize);
    for(int i=0; i<local.n; ++i){
        generateHarmonics(local[i],bonds,bondParams,angles,angleParams, networkA);
    }
    for(int i=0; i<fixedInner.n; ++i){
        generateHarmonics(fixedInner[i],bonds,bondParams,angles,angleParams, networkA);
    }

    /*
    OutputFile bondfile("bonds.dat");
    bondfile.initVariables(6,4,60,20);
    for (int i=0; i<bonds.n/2;++i){
        bondfile.write(to_string(bonds[2*i])+"    "+to_string(bonds[2*i+1]));
    }


    OutputFile anglefile("angles.dat");
    anglefile.initVariables(6,4,60,20);
    for (int i=0; i<angles.n/3;++i){
        anglefile.write(to_string(angles[3*i])+"    "+to_string(angles[3*i+1]) + "    "+to_string(angles[3*i+2]));
    }
    */

    //Intersections - Expensive, turn off and use in global potential energy
    VecR<int> intersections(0,local.n*1000);
//    if(useIntx) {
//        for(int i=0; i<localDual.n; ++i){
//            generateRingIntersections(localDual[i],intersections);
//        }
//    }

    //Assign to fixed size arrays
    VecF<int> bnds(bonds.n),fixd(fixedInner.n+fixedOuter.n),angs(angles.n),intx(intersections.n);
    VecF<double> bndP(bondParams.n),angP(angleParams.n),intxP; //placeholder for intersections
    for(int i=0; i<bnds.n; ++i) bnds[i]=bonds[i];
    for(int i=0; i<bndP.n; ++i) bndP[i]=bondParams[i];
    for(int i=0; i<angs.n; ++i) angs[i]=angles[i];
    for(int i=0; i<angP.n; ++i) angP[i]=angleParams[i];
    for(int i=0; i<fixedInner.n; ++i) fixd[i]=fixedInner[i];
    for(int i=0; i<fixedOuter.n; ++i) fixd[i+fixedInner.n]=fixedOuter[i];
    for(int i=0; i<intersections.n; ++i) intx[i]=intersections[i];

    ///cout << endl << endl << endl;
    //cout << bonds.n << "   " << angles.n << endl << endl << endl;

    //Geometry optimise
    VecF<int> optStatus(2);
    if(networkA.geometryCode=="2DE"){
        if(!restrict) {
            HI2DP potModel(networkA.pb[0], networkA.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            potModel.setFixedAtoms(fixd);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//optimise if no line intersections
                SteepestDescentArmijoMultiDim<HI2DP> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;
            }
        }
        else{
            HRI2DP potModel(networkA.pb[0], networkA.pb[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            potModel.setFixedAtoms(fixd);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//optimise if no line intersections
                SteepestDescentArmijoMultiDim<HRI2DP> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;
            }
        }
    }
    else if(networkA.geometryCode=="2DS"){
        VecF<int> constrained(networkA.nodes.n);
        for(int i=0; i<networkA.nodes.n; ++i) constrained[i]=i;
        if(!restrict) {
            HI3DS potModel;
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            potModel.setFixedAtoms(fixd);
            potModel.setGeomConstraints(constrained, potParamsC);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//optimise if no line intersections
                SteepestDescentArmijoMultiDim<HI3DS> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;
            }
        }
        else {
            HRI3DS potModel;
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            potModel.setFixedAtoms(fixd);
            potModel.setGeomConstraints(constrained, potParamsC);
            if (useIntx) potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) < numeric_limits<double>::infinity()) {//optimise if no line intersections
                SteepestDescentArmijoMultiDim<HRI3DS> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;
            }
        }
    }
    return optStatus;
}

//Generate harmonic potentials corresponding to bonds and angles associated with a given node in lattice A
void LinkedNetwork::generateHarmonics(int id, VecR<int> &bonds, VecR<double> &bondParams, VecR<int> &angles, VecR<double> &angleParams, Network network) {

    //Potential parameters
    double bk=potParamsB[0], br=potParamsB[1]; //bond k and r0
    double bk0=bk, br0=sqrt(3)*br;
    //cout << br << " " << br0 << endl;
    double ak, act; //angle k, cos theta0

    //Harmonics
    int cnd; //coordination
    int id0 = id, id1, id2;
    cnd = network.nodes[id0].netCnxs.n;
//    if (cnd == 3) {
//        ak = potParamsA[0];
//        act = potParamsA[1];
//    } else if (cnd == 4) {
//        ak = potParamsA[2];
//        act = potParamsA[3];
//    }
    ak=potParamsA[0];
    act=cos(2.0*M_PI/cnd);
    for (int i = 0; i < cnd; ++i) {
        id1 = network.nodes[id0].netCnxs[i];
        id2 = network.nodes[id0].netCnxs[(i + 1) % cnd];
        //bonds
        if (id0 < id1) {
            //prevent double counting
            //if (id0>networkA.nodes.n){cout<<"large"<<endl;}
            //if (id1>networkA.nodes.n){cout<<"large"<<endl;}
            if (id0>=networkA.nodes.n && id1>=networkA.nodes.n){
                bonds.addValue(id0);
                bonds.addValue(id1);
                bondParams.addValue(bk);
                bondParams.addValue(br0);
            }
            else{
                bonds.addValue(id0);
                bonds.addValue(id1);
                bondParams.addValue(bk);
                bondParams.addValue(br);
            }

        }
        //angles
        angles.addValue(id0);
        angles.addValue(id1);
        angles.addValue(id2);
        angleParams.addValue(ak);
        angleParams.addValue(act);
    }
}

//Generate harmonic potentials corresponding to bonds and angles associated with a given node in lattice A
void LinkedNetwork::generateHarmonicsOnly(int id, VecR<int> &bonds, VecR<double> &bondParams, Network network) {

    //Potential parameters
    double bk=potParamsB[0], br=potParamsB[1]; //bond k and r0
    double bk0=bk, br0=sqrt(3)*br;
    //cout << br << " " << br0 << endl;
    //double ak, act; //angle k, cos theta0

    //Harmonics
    int cnd; //coordination
    int id0 = id, id1, id2;
    cnd = network.nodes[id0].netCnxs.n;
//    if (cnd == 3) {
//        ak = potParamsA[0];
//        act = potParamsA[1];
//    } else if (cnd == 4) {
//        ak = potParamsA[2];
//        act = potParamsA[3];
//    }
    //ak=potParamsA[0];
    //act=cos(2.0*M_PI/cnd);

    int sio_bond_count=0;
    int o_o_bond_count=0;


    for (int i = 0; i < cnd; ++i) {
        id1 = network.nodes[id0].netCnxs[i];
        id2 = network.nodes[id0].netCnxs[(i + 1) % cnd];
        //bonds
        if (id0 < id1) {
            if (network.nodes.n>networkA.nodes.n){
                //prevent double counting
                //if (id0>networkA.nodes.n){cout<<"large"<<endl;}
                //if (id1>networkA.nodes.n){cout<<"large"<<endl;}
                if (id0>=networkA.nodes.n && id1>=networkA.nodes.n){
                    bonds.addValue(id0);
                    bonds.addValue(id1);
                    bondParams.addValue(bk);
                    bondParams.addValue(br0/2);
                    sio_bond_count +=1;
//                    cout << id0 << " " << id1 << " " << bk << " " << br0/2 << endl;
                }
                else{
                    bonds.addValue(id0);
                    bonds.addValue(id1);
                    bondParams.addValue(bk);
                    bondParams.addValue(br/2);
                    o_o_bond_count +=1;
//                    cout << id0 << " " << id1 << " " << bk << " " << br/2 << endl;

                }
            }
            else{
                bonds.addValue(id0);
                bonds.addValue(id1);
                bondParams.addValue(bk);
                bondParams.addValue(br);
//                cout << id0 << " " << id1 << " " << bk << " " << br << endl;

            }


        }

//        //angles
//        angles.addValue(id0);
//        angles.addValue(id1);
//        angles.addValue(id2);
//        angleParams.addValue(ak);
//        angleParams.addValue(act);
    }

}

//Generate ring edge intersections for a specific ring
void LinkedNetwork::generateRingIntersections(int rId, VecR<int> &intersections) {

    int rCnd=networkB.nodes[rId].netCnxs.n;
    int nCnd=networkB.nodes[rId].dualCnxs.n;
    int e0,e1,e2,e3;
    for(int i=0; i<rCnd; ++i){//loop over neighbouring rings
        int rId0=networkB.nodes[rId].netCnxs[i];
        if(rId<rId0){//prevent double counting
            for(int j=0; j<nCnd; ++j){//loop over nodes
                e0=networkB.nodes[rId].dualCnxs[j];
                e1=networkB.nodes[rId].dualCnxs[(j+1)%nCnd];
                int nCnd0=networkB.nodes[rId0].dualCnxs.n;
                for(int k=0; k<nCnd0; ++k){
                    e2=networkB.nodes[rId0].dualCnxs[k];
                    e3=networkB.nodes[rId0].dualCnxs[(k+1)%nCnd0];
                    if(e0!=e2 && e0!=e3 && e1!=e2 && e1!=e3) {
                        intersections.addValue(e0);
                        intersections.addValue(e1);
                        intersections.addValue(e2);
                        intersections.addValue(e3);
                    }
                }
            }
        }
    }
}

//Generate intersections required to maintain convexity for a given node
void LinkedNetwork::generateConvexIntersections(int nId, VecR<int> &intersections) {

    int cnd=networkA.nodes[nId].netCnxs.n;
    int id0,id1,id2;
    for (int i = 0; i < cnd; ++i) {
        id0 = networkA.nodes[nId].netCnxs[i];
        id1 = networkA.nodes[nId].netCnxs[(i + 1) % cnd];
        id2 = networkA.nodes[nId].netCnxs[(i + 2) % cnd];
        intersections.addValue(id0);
        intersections.addValue(id1);
        intersections.addValue(nId);
        intersections.addValue(id2);
    }

}


//Update networks with geometry optimised coordinates
void LinkedNetwork::syncCoordinates() {

    //Sync T coordinates
    for (int i=0; i<networkA.nodes.n;++i){
        networkT.nodes[i].crd[0]=crds[2*i];
        networkT.nodes[i].crd[1]=crds[2*i+1];
    }
/*
    float sum = 0.0;
    for (int i=0;i<crds.n;++i){
        sum += pow(crds[i]-networkT.nodes[(i-i%2)/2].crd[i%2],2);
    }
    cout << "Sync Coordinates sum : " << sum << endl;
//    for (int i=0;i<network.nodes.n;++i) {
//        cout << "Atom "<<i<<" for ref : " << network.nodes[i].crd[0] << " " << network.nodes[i].crd[1] << endl;
//    }
*/
    //Sync A coordinates
    if(networkA.geometryCode=="2DE"){
        for(int i=0; i<networkA.nodes.n; ++i){
            networkA.nodes[i].crd[0]=crds[2*i];
            networkA.nodes[i].crd[1]=crds[2*i+1];
        }
    }
    else{
        for(int i=0; i<networkA.nodes.n; ++i){
            networkA.nodes[i].crd[0]=crds[3*i];
            networkA.nodes[i].crd[1]=crds[3*i+1];
            networkA.nodes[i].crd[2]=crds[3*i+2];
        }
    }

    //Sync B coordinates
    if(networkB.geometryCode=="2DE"){
        for(int i=0; i<networkB.nodes.n; ++i){
            VecF<double> x(networkB.nodes[i].dualCnxs.n);
            VecF<double> y(networkB.nodes[i].dualCnxs.n);
            for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
                x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
                y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
            }
            VecF<double> origin(2);
            origin[0]=x[0];
            origin[1]=y[0];
            x-=origin[0];
            y-=origin[1];
            for(int j=0; j<x.n; ++j) x[j]-=networkB.pb[0]*nearbyint(x[j]*networkB.rpb[0]);
            for(int j=0; j<y.n; ++j) y[j]-=networkB.pb[1]*nearbyint(y[j]*networkB.rpb[1]);
            VecF<double> c(2);
            c[0]=origin[0]+vMean(x);
            c[1]=origin[1]+vMean(y);
            networkB.nodes[i].crd=c;
        }
    }
    else{
        for (int i = 0; i < networkB.nodes.n; ++i) {
            VecF<double> x(networkB.nodes[i].dualCnxs.n);
            VecF<double> y(networkB.nodes[i].dualCnxs.n);
            VecF<double> z(networkB.nodes[i].dualCnxs.n);
            for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
                x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
                y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
                z[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[2];
            }
            VecF<double> origin(3);
            origin[0] = x[0];
            origin[1] = y[0];
            origin[2] = z[0];
            x -= origin[0];
            y -= origin[1];
            z -= origin[2];
            for (int j = 0; j < x.n; ++j) x[j] -= networkB.pb[0] * nearbyint(x[j] * networkB.rpb[0]);
            for (int j = 0; j < y.n; ++j) y[j] -= networkB.pb[1] * nearbyint(y[j] * networkB.rpb[1]);
            for (int j = 0; j < z.n; ++j) z[j] -= networkB.pb[2] * nearbyint(z[j] * networkB.rpb[2]);
            VecF<double> c(3);
            c[0] = origin[0] + vMean(x);
            c[1] = origin[1] + vMean(y);
            c[2] = origin[2] + vMean(z);
            networkB.nodes[i].crd = c;
        }
    }

 }

//Update networks with geometry optimised coordinates
void LinkedNetwork::syncCoordinatesTD() {

    //Sync A coordinates
    if(networkA.geometryCode=="2DE"){
        for(int i=0; i<networkA.nodes.n; ++i){
            networkA.nodes[i].crd[0]=crds[2*i];
            networkA.nodes[i].crd[1]=crds[2*i+1];
        }
        for (int i=0; i<networkT.nodes.n;++i){
            networkT.nodes[i].crd[0]=crds[2*i];
            networkT.nodes[i].crd[1]=crds[2*i+1];
        }
    }
    else{
        for(int i=0; i<networkA.nodes.n; ++i){
            networkA.nodes[i].crd[0]=crds[3*i];
            networkA.nodes[i].crd[1]=crds[3*i+1];
            networkA.nodes[i].crd[2]=crds[3*i+2];
        }
    }

    //Sync B coordinates
    if(networkB.geometryCode=="2DE"){
        for(int i=0; i<networkB.nodes.n; ++i){
            VecF<double> x(networkB.nodes[i].dualCnxs.n);
            VecF<double> y(networkB.nodes[i].dualCnxs.n);
            for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
                x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
                y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
            }
            VecF<double> origin(2);
            origin[0]=x[0];
            origin[1]=y[0];
            x-=origin[0];
            y-=origin[1];
            for(int j=0; j<x.n; ++j) x[j]-=networkB.pb[0]*nearbyint(x[j]*networkB.rpb[0]);
            for(int j=0; j<y.n; ++j) y[j]-=networkB.pb[1]*nearbyint(y[j]*networkB.rpb[1]);
            VecF<double> c(2);
            c[0]=origin[0]+vMean(x);
            c[1]=origin[1]+vMean(y);
            networkB.nodes[i].crd=c;
        }
    }
    else{
        for (int i = 0; i < networkB.nodes.n; ++i) {
            VecF<double> x(networkB.nodes[i].dualCnxs.n);
            VecF<double> y(networkB.nodes[i].dualCnxs.n);
            VecF<double> z(networkB.nodes[i].dualCnxs.n);
            for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
                x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
                y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
                z[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[2];
            }
            VecF<double> origin(3);
            origin[0] = x[0];
            origin[1] = y[0];
            origin[2] = z[0];
            x -= origin[0];
            y -= origin[1];
            z -= origin[2];
            for (int j = 0; j < x.n; ++j) x[j] -= networkB.pb[0] * nearbyint(x[j] * networkB.rpb[0]);
            for (int j = 0; j < y.n; ++j) y[j] -= networkB.pb[1] * nearbyint(y[j] * networkB.rpb[1]);
            for (int j = 0; j < z.n; ++j) z[j] -= networkB.pb[2] * nearbyint(z[j] * networkB.rpb[2]);
            VecF<double> c(3);
            c[0] = origin[0] + vMean(x);
            c[1] = origin[1] + vMean(y);
            c[2] = origin[2] + vMean(z);
            networkB.nodes[i].crd = c;
        }
    }
}

//Get normalised probability distribution of nodes of each size in given lattice
VecF<double> LinkedNetwork::getNodeDistribution(string lattice) {

    if(lattice=="A") return networkA.getNodeDistribution();
    else return networkB.getNodeDistribution();
}

//Get unnormalised probability distribution of node connections in given lattice
VecF< VecF<int> > LinkedNetwork::getEdgeDistribution(string lattice) {

    if(lattice=="A") return networkA.edgeDistribution;
    else return networkB.edgeDistribution;
}

//Get Aboav-Weaire fitting parameters
VecF<double> LinkedNetwork::getAboavWeaire(string lattice) {

    if(lattice=="A") return networkA.aboavWeaireParams();
    else return networkB.aboavWeaireParams();
}

//Get assortativity
double LinkedNetwork::getAssortativity(string lattice) {

    if(lattice=="A") return networkA.assortativity();
    else return networkB.assortativity();

}

//Get alpha estimate
double LinkedNetwork::getAboavWeaireEstimate(string lattice) {

    if(lattice=="A") return networkA.aboavWeaireEstimate();
    else return networkB.aboavWeaireEstimate();
}

//Get entropy
VecF<double> LinkedNetwork::getEntropy(string lattice) {

    if(lattice=="A") return networkA.entropy();
    else return networkB.entropy();
}

//Get cluster information
double LinkedNetwork::getMaxCluster(string lattice, int nodeCnd){

    if(lattice=="A") return networkA.maxCluster(nodeCnd);
    else return networkB.maxCluster(nodeCnd);
}

//Get cluster information
VecF<int> LinkedNetwork::getMaxClusters(string lattice, int minCnd, int maxCnd){

    if(lattice=="A") return networkA.maxClusters(minCnd,maxCnd,3,2);
    else return networkB.maxClusters(minCnd,maxCnd,3,2);
}

//Check linked networks for consistency
bool LinkedNetwork::checkConsistency() {

    bool checkCnx=checkCnxConsistency();
    bool checkDesc=checkDescriptorConsistency();
    bool check=checkCnx*checkDesc;

    return check;
}

//Check linked networks have mutual network and dual connections
bool LinkedNetwork::checkCnxConsistency() {

    //check number of network connections is equal to number of dual connections
    bool netDualEquality=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        if(networkA.nodes[i].netCnxs.n!=networkA.nodes[i].dualCnxs.n) netDualEquality=false;
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        if(networkB.nodes[i].netCnxs.n!=networkB.nodes[i].dualCnxs.n) netDualEquality=false;
    }

    //check mutual network connections
    bool mutualNetCnx=true;
    int id0, id1;
    for(int i=0; i<networkA.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkA.nodes[i].netCnxs.n; ++j){
            id1=networkA.nodes[i].netCnxs[j];
            mutualNetCnx=vContains(networkA.nodes[id1].netCnxs,id0);
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkB.nodes[i].netCnxs.n; ++j){
            id1=networkB.nodes[i].netCnxs[j];
            mutualNetCnx=vContains(networkB.nodes[id1].netCnxs,id0);
        }
    }

    //check mutual dual connections
    bool mutualDualCnx=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkA.nodes[i].dualCnxs.n; ++j){
            id1=networkA.nodes[i].dualCnxs[j];
            mutualDualCnx=vContains(networkB.nodes[id1].dualCnxs,id0);
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkB.nodes[i].dualCnxs.n; ++j){
            id1=networkB.nodes[i].dualCnxs[j];
            mutualDualCnx=vContains(networkA.nodes[id1].dualCnxs,id0);
        }
    }

    //check network connections are neighbours by lying on same ring (some highly strained cases could give a false positive)
    bool nbNetCnx=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        int nCnxs=networkA.nodes[i].netCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkA.nodes[i].netCnxs[j];
            id1=networkA.nodes[i].netCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkA.nodes[id0].dualCnxs,networkA.nodes[id1].dualCnxs);
            if(common.n==0) nbNetCnx=false;
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int nCnxs=networkB.nodes[i].netCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkB.nodes[i].netCnxs[j];
            id1=networkB.nodes[i].netCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkB.nodes[id0].dualCnxs,networkB.nodes[id1].dualCnxs);
            if(common.n==0) nbNetCnx=false;
        }
    }

    //check dual connections are neighbours by lying on same ring (some highly strained cases could give a false positive)
    bool nbDualCnx=true;
    for(int i=0; i<networkA.nodes.n; ++i){
        int nCnxs=networkA.nodes[i].dualCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkA.nodes[i].dualCnxs[j];
            id1=networkA.nodes[i].dualCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkB.nodes[id0].dualCnxs,networkB.nodes[id1].dualCnxs);
            common.delValue(i);
            if(common.n==0) nbDualCnx=false;
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int nCnxs=networkB.nodes[i].dualCnxs.n;
        for(int j=0; j<nCnxs; ++j){
            id0=networkB.nodes[i].dualCnxs[j];
            id1=networkB.nodes[i].dualCnxs[(j+1)%nCnxs];
            VecR<int> common=vCommonValues(networkA.nodes[id0].dualCnxs,networkA.nodes[id1].dualCnxs);
            common.delValue(i);
            if(common.n==0) nbDualCnx=false;
        }
    }

    //check expected number of auxiliary connections
    bool numAux=true;
    for(int i=0; i<networkB.nodes.n; ++i){
        int expAux = 0;
        for(int j=0; j<networkB.nodes[i].dualCnxs.n; ++j) expAux+=networkA.nodes[networkB.nodes[i].dualCnxs[j]].netCnxs.n-3;
        if(networkB.nodes[i].auxCnxs.n!=expAux) numAux=false;
    }
    numAux=true;

    //check mutual auxiliary connections
    bool mutualAuxCnx=true;
    for(int i=0; i<networkB.nodes.n; ++i){
        id0=i;
        for(int j=0; j<networkB.nodes[i].auxCnxs.n; ++j){
            id1=networkB.nodes[i].auxCnxs[j];
            mutualAuxCnx=vContains(networkB.nodes[id1].auxCnxs,id0);
        }
    }
    mutualAuxCnx=true;

    //overall flag
    bool consistent=netDualEquality*mutualNetCnx*mutualDualCnx*nbNetCnx*nbDualCnx*numAux*mutualAuxCnx;

    return consistent;
}

//Check linked networks have accurate descriptors
bool LinkedNetwork::checkDescriptorConsistency() {

    VecF<int> checkNodeA(networkA.nodeDistribution);
    VecF<int> checkNodeB(networkB.nodeDistribution);
    VecF< VecF<int> > checkEdgeA(networkA.edgeDistribution);
    VecF< VecF<int> > checkEdgeB(networkB.edgeDistribution);

    //Check node distribution
    bool nodeA, nodeB;
    checkNodeA=0;
    checkNodeB=0;
    for(int i=0; i<networkA.nodes.n; ++i){
        int n=networkA.nodes[i].netCnxs.n;
        ++checkNodeA[n];
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int n=networkB.nodes[i].netCnxs.n;
        ++checkNodeB[n];
    }
    nodeA=checkNodeA==networkA.nodeDistribution;
    nodeB=checkNodeB==networkB.nodeDistribution;


    //Check edge distribution
    bool edgeA=true, edgeB=true;
    for(int i=0; i<checkEdgeA.n; ++i) checkEdgeA[i]=0;
    for(int i=0; i<checkEdgeB.n; ++i) checkEdgeB[i]=0;
    for(int i=0; i<networkA.nodes.n; ++i){
        int m=networkA.nodes[i].netCnxs.n;
        for(int j=0; j<m; ++j){
            int n=networkA.nodes[networkA.nodes[i].netCnxs[j]].netCnxs.n;
            ++checkEdgeA[m][n];
        }
    }
    for(int i=0; i<networkB.nodes.n; ++i){
        int m=networkB.nodes[i].netCnxs.n;
        for(int j=0; j<m; ++j){
            int n=networkB.nodes[networkB.nodes[i].netCnxs[j]].netCnxs.n;
            ++checkEdgeB[m][n];
        }
    }
    for(int i=0; i<checkEdgeA.n; ++i){
        if(!(checkEdgeA[i]==networkA.edgeDistribution[i])){
            edgeA=false;
            break;
        }
    };
    for(int i=0; i<checkEdgeB.n; ++i){
        if(!(checkEdgeB[i]==networkB.edgeDistribution[i])){
//            for(int j=0; j<checkEdgeB[i].n; ++j){
//                if(checkEdgeB[i][j]!=networkB.edgeDistribution[i][j]) cout<<"+++ "<<i<<" "<<j<<" "<<checkEdgeB[i][j]<<endl;
//            }
            edgeB=false;
            break;
        }
    };

    //Overall flag
    bool consistent=nodeA*nodeB*edgeA*edgeB;

    return consistent;
}

//Check convexity of all angles
bool LinkedNetwork::checkConvexity() {

    bool convex;
    for(int i=0; i<networkA.nodes.n; ++i){
        convex=checkConvexity(i);
        if(!convex) break;
    }

    return convex;
}

//Check convexity by summing angles around node
bool LinkedNetwork::checkConvexity(int id) {

    double angleSum=0.0;
    if(networkA.geometryCode=="2DE") {
        //Coordinate vectors, coordination and pbc
        VecF<double> v0(2), v1(2), v2(2);
        v0[0] = crds[2 * id];
        v0[1] = crds[2 * id + 1];
        int cnd = networkA.nodes[id].netCnxs.n;
        int id1, id2;
        double pbx = networkA.pb[0], pby = networkA.pb[1];
        double pbrx = networkA.rpb[0], pbry = networkA.rpb[1];
        //Determine vectors to neighbours and sum angles
        for (int i = 0; i < cnd; ++i) {
            int j = (i + 1) % cnd;
            id1 = networkA.nodes[id].netCnxs[i];
            id2 = networkA.nodes[id].netCnxs[j];
            //Periodic vectors to adjacent neighbours
            v1[0] = crds[2 * id1];
            v1[1] = crds[2 * id1 + 1];
            v2[0] = crds[2 * id2];
            v2[1] = crds[2 * id2 + 1];
            v1 -= v0;
            v2 -= v0;
            v1[0] -= pbx * nearbyint(v1[0] * pbrx);
            v1[1] -= pby * nearbyint(v1[1] * pbry);
            v2[0] -= pbx * nearbyint(v2[0] * pbrx);
            v2[1] -= pby * nearbyint(v2[1] * pbry);
            //Angle from dot product
            double n1, n2;
            angleSum += vAngle(v1, v2, n1, n2);
        }
    }
    else if (networkA.geometryCode=="2DS"){
        //Project neighbour vectors onto tangent plane of sphere, then sum angles
        VecF<double> v0(3), v1(3), v2(3);
        v0[0] = crds[3 * id];
        v0[1] = crds[3 * id + 1];
        v0[2] = crds[3 * id + 2];
        double n0=vNorm(v0);
        VecF<double> normal=v0/n0;
        int cnd = networkA.nodes[id].netCnxs.n;
        int id1, id2;
        for (int i = 0; i < cnd; ++i) {
            int j = (i + 1) % cnd;
            id1 = networkA.nodes[id].netCnxs[i];
            id2 = networkA.nodes[id].netCnxs[j];
            v1[0] = crds[3 * id1];
            v1[1] = crds[3 * id1 + 1];
            v1[2] = crds[3 * id1 + 2];
            v2[0] = crds[3 * id2];
            v2[1] = crds[3 * id2 + 1];
            v2[2] = crds[3 * id2 + 2];
            v1 -= v0;
            v2 -= v0;
            v1-=normal*vSum(v1*normal);
            v2-=normal*vSum(v2*normal);
            //Angle from dot product
            double n1, n2;
            angleSum += vAngle(v1, v2, n1, n2);
        }
    }

    if (fabs(angleSum - 2 * M_PI) < 1e-12) return true;
    else return false;
}

//Calculate bond length and angle mean and standard deviation
VecF<double> LinkedNetwork::getOptimisationGeometry(Network network, VecF<double> &lenHist, VecF<double> &angHist) {

    //Calculate for current configuration
    double x=0.0,xSq=0.0,y=0.0,ySq=0.0; //len and angle
    int xN=0,yN=0; //count
    int cnd;
    VecF<double> v0(2),v1(2),v2(2);
    double pbx=network.pb[0],pby=network.pb[1];
    double pbrx=network.rpb[0],pbry=network.rpb[1];
    double lenBinWidth = 4.0/10000.0;
    double angBinWidth = 2*M_PI/10000.0;
    double bin;
    for(int i=0; i<network.nodes.n; ++i){
        cnd=network.nodes[i].netCnxs.n;
        v0[0]=crds[2*i];
        v0[1]=crds[2*i+1];
        for(int j=0; j<cnd; ++j){
            int id1=network.nodes[i].netCnxs[j];
            int id2=network.nodes[i].netCnxs[(j+1)%cnd];
            v1[0]=crds[2*id1];
            v1[1]=crds[2*id1+1];
            v2[0]=crds[2*id2];
            v2[1]=crds[2*id2+1];
            v1-=v0;
            v2-=v0;
            v1[0]-=pbx*nearbyint(v1[0]*pbrx);
            v1[1]-=pby*nearbyint(v1[1]*pbry);
            v2[0]-=pbx*nearbyint(v2[0]*pbrx);
            v2[1]-=pby*nearbyint(v2[1]*pbry);
            double n1,n2;
            double theta=vAngle(v1,v2,n1,n2);
            //Edge lengths avoiding double counting
            if(i<id1){
                x+=n1;
                xSq+=n1*n1;
                xN+=1;
                bin = floor(n1/lenBinWidth);
                if(bin<lenHist.n) lenHist[bin] += 1.0;
            }
            //Angles
            y+=theta;
            ySq+=theta*theta;
            yN+=1;
            bin = floor(theta/angBinWidth);
            if (bin<angHist.n) angHist[bin] += 1.0;
        }
    }

    //Return current configuration
    VecF<double> optGeom(8);
    optGeom[0]=x;
    optGeom[1]=xSq;
    optGeom[2]=x/xN;
    optGeom[3]=sqrt(xSq/xN-optGeom[2]*optGeom[2]);
    optGeom[4]=y;
    optGeom[5]=ySq;
    optGeom[6]=y/yN;
    optGeom[7]=ySq/yN-optGeom[6]*optGeom[6];
    if(optGeom[7]<0.0) optGeom[7]=0.0;
    else optGeom[7]=sqrt(optGeom[7]);

    return optGeom;
}

//Calculate bond length and angle mean and standard deviation
VecF<double> LinkedNetwork::getOptimisationGeometryTD(VecF<double> &lenHist, VecF<double> &angHist) {

    //Calculate for current configuration
    double x=0.0,xSq=0.0,y=0.0,ySq=0.0; //len and angle
    int xN=0,yN=0; //count
    int cnd;
    VecF<double> v0(2),v1(2),v2(2);
    double pbx=networkT.pb[0],pby=networkT.pb[1];
    double pbrx=networkT.rpb[0],pbry=networkT.rpb[1];
    double lenBinWidth = 4.0/10000.0;
    double angBinWidth = 2*M_PI/10000.0;
    double bin;
    for(int i=0; i<networkT.nodes.n; ++i){
        cnd=networkT.nodes[i].netCnxs.n;
        v0[0]=crds[2*i];
        v0[1]=crds[2*i+1];
        for(int j=0; j<cnd; ++j){
            int id1=networkT.nodes[i].netCnxs[j];
            int id2=networkT.nodes[i].netCnxs[(j+1)%cnd];
            v1[0]=crds[2*id1];
            v1[1]=crds[2*id1+1];
            v2[0]=crds[2*id2];
            v2[1]=crds[2*id2+1];
            v1-=v0;
            v2-=v0;
            v1[0]-=pbx*nearbyint(v1[0]*pbrx);
            v1[1]-=pby*nearbyint(v1[1]*pbry);
            v2[0]-=pbx*nearbyint(v2[0]*pbrx);
            v2[1]-=pby*nearbyint(v2[1]*pbry);
            double n1,n2;
            double theta=vAngle(v1,v2,n1,n2);
            //Edge lengths avoiding double counting
            if(i<id1){
                x+=n1;
                xSq+=n1*n1;
                xN+=1;
                bin = floor(n1/lenBinWidth);
                if(bin<lenHist.n) lenHist[bin] += 1.0;
            }
            //Angles
            y+=theta;
            ySq+=theta*theta;
            yN+=1;
            bin = floor(theta/angBinWidth);
            if (bin<angHist.n) angHist[bin] += 1.0;
        }
    }

    //Return current configuration
    VecF<double> optGeom(8);
    optGeom[0]=x;
    optGeom[1]=xSq;
    optGeom[2]=x/xN;
    optGeom[3]=sqrt(xSq/xN-optGeom[2]*optGeom[2]);
    optGeom[4]=y;
    optGeom[5]=ySq;
    optGeom[6]=y/yN;
    optGeom[7]=ySq/yN-optGeom[6]*optGeom[6];
    if(optGeom[7]<0.0) optGeom[7]=0.0;
    else optGeom[7]=sqrt(optGeom[7]);

    return optGeom;
}

//Get sum of areas and squared areas for each ring size
void LinkedNetwork::getRingAreas(VecF<double> &areaSum, VecF<double> &areaSqSum) {

    //Loop over rings, recentre and apply shoelace formula
    areaSum = 0.0;
    areaSqSum = 0.0;
    double pbx=networkA.pb[0],pby=networkA.pb[1];
    double pbrx=networkA.rpb[0],pbry=networkA.rpb[1];
    for(int i=0; i<networkB.nodes.n; ++i){
        VecR<int> ids = networkB.nodes[i].dualCnxs;
        int ringSize = ids.n;
        VecF<double> xCrds(ringSize), yCrds(ringSize);
        for(int j=0; j<ringSize; ++j){
            xCrds[j] = crds[2*ids[j]];
            yCrds[j] = crds[2*ids[j]+1];
        }
        xCrds = xCrds-xCrds[0];
        yCrds = yCrds-yCrds[0];
        for(int j=1; j<ringSize; ++j){
            xCrds[j]-=pbx*nearbyint(xCrds[j]*pbrx);
            yCrds[j]-=pby*nearbyint(yCrds[j]*pbry);
        }
        double a = 0.0;
        for(int j=0; j<ringSize-1; ++j) a += xCrds[j]*yCrds[j+1];
        for(int j=0; j<ringSize-1; ++j) a -= xCrds[j+1]*yCrds[j];
        a = 0.5*fabs(a+xCrds[ringSize-1]*yCrds[0]-xCrds[0]*yCrds[ringSize-1]);
        areaSum[ringSize] += a;
        areaSqSum[ringSize] += a*a;
    }
}

//Wrap coordinates of lattice A if periodic
void LinkedNetwork::wrapCoordinates() {

    if (networkT.geometryCode=="2DEtr"){
        HLJ2DP potModel(networkT.pb[0],networkT.pb[1]);
        potModel.wrap(crds);
    }

    if(networkA.geometryCode=="2DE"){
        HI2DP potModel(networkA.pb[0],networkA.pb[1]);
        potModel.wrap(crds);
    }

}

//Write xyz file format of networks
void LinkedNetwork::writeXYZ(string prefix) {
    networkA.writeXYZ(prefix+"_A","O");
    networkB.writeXYZ(prefix+"_B","N");

//    for (int i=0;i<networkT.nodes.n;++i) {
//        cout << "Atom "<<i<<" for ref : " << networkT.nodes[i].crd[0] << " " << networkT.nodes[i].crd[1] << endl;
//    }
    networkT.writeXYZ(prefix+"_T", "Si");
}
void LinkedNetwork::writeBilayer(string prefix, float lj_cutoff) {
    networkA.writeBilayerA(prefix, lj_cutoff);
    networkB.writeBilayerB(prefix);
}
//Write networks in format that can be loaded and visualised
void LinkedNetwork::write(string prefix) {
    networkA.write(prefix+"_A");
    networkB.write(prefix+"_B");
    networkT.write(prefix+"_T");



}


