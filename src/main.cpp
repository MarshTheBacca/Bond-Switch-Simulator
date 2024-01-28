#include <iostream>
#include <sstream>
#include "output_file.h"
#include "linked_network.h"
#include "lammps_object.h"
#include <sys/stat.h>
#include <unistd.h>
#include "input_data.h"
#include <iomanip>
#include <ctime>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
using namespace std;

// You were changing variable names

int main()
{
    // Create a file sink and a console sink with different names for clarity
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("./netmc.log", true);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

    // Combine the sinks into a multi-sink logger
    auto logger = std::make_shared<spdlog::logger>("multi_sink", spdlog::sinks_init_list{file_sink, console_sink});
    spdlog::set_level(spdlog::level::info);
    spdlog::register_logger(logger);
    logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");

    // Log some messages
    logger->info("Simulation Started");
    logger->info("Network Monte Carlo");
    logger->info("Written by David Ormrod Morely, Edited by Oliver Whitaker and Marshall Hunt, Wilson Group, 2024");

    // Read input file
    InputData inputData("./netmc.inpt", logger);

    OutputFile lammpsBondAngle(inputData.inputFolder + "/PARM_Si.lammps");
    lammpsBondAngle.write("bond_style harmonic");
    lammpsBondAngle.write("bond_coeff 1 " + to_string(inputData.harmonicBondForceConstant) + " 1.000");
    lammpsBondAngle.write("angle_style cosine/squared");
    lammpsBondAngle.write("angle_coeff 1 " + to_string(inputData.harmonicAngleForceConstant) + " 120");

    // Initialise network
    logger->info("Initialising network");
    logger->info("Number of rings: {}", inputData.numRings);
    logger->info("Min ring size: {}", inputData.minRingSize);
    logger->info("Max ring size: {}", inputData.maxRingSize);
    logger->info("Min node coordination: {}", inputData.minCoordination);
    logger->info("Max node coordination: {}", inputData.maxCoordination);
    logger->info("Monte Carlo move type: {}", inputData.moveType);
    logger->info("Monte carlo initial log temperature: {}", inputData.startTemperature);
    logger->info("Monte carlo final log temperature: {}", inputData.endTemperature);
    logger->info("Monte carlo log temperature increment: {}", inputData.temperatureIncrement);
    logger->info("Monte carlo steps per increment: {}", inputData.stepsPerTemperature);
    logger->info("Equilibrium steps: {}", inputData.initialThermalisationSteps);

    LinkedNetwork network(inputData.inputFolder, inputData.inputFilePrefix,
                          inputData.outputFolder,
                          inputData.minCoordination, inputData.maxCoordination,
                          inputData.minRingSize, inputData.maxRingSize,
                          inputData.isSimpleGrapheneEnabled,
                          inputData.isTriangleRaftEnabled, inputData.isBilayerEnabled,
                          inputData.isTersoffGrapheneEnabled,
                          inputData.isBNEnabled,
                          inputData.isRestartUsingLammpsObjectsEnabled, logger);

    int minCnd = network.minACnxs;
    int maxCnd = network.maxACnxs;
    int minRingSize = network.minBCnxs;
    int maxRingSize = network.maxBCnxs;

    // Find Fixed Rings
    network.pushPrefix(inputData.inputFolder, inputData.outputFolder);
    network.findFixedRings(inputData.isFixRingsEnabled, inputData.inputFolder + "/" + inputData.fixedRingsFile, logger);

    logger->info("Initialising potential model...");
    network.initialisePotentialModel(network.networkA, inputData.harmonicAngleForceConstant,
                                     inputData.harmonicBondForceConstant,
                                     inputData.harmonicGeometryConstraint,
                                     inputData.isMaintainConvexityEnabled, logger);

    logger->info("Initialising Geometry Optimisation...");
    network.initialiseGeometryOpt(inputData.monteCarloLocalMaxIterations,
                                  inputData.tauBacktrackingParameter,
                                  inputData.tolerance,
                                  inputData.localRegionSize);

    logger->info("Initialising Monte Carlo...");
    network.initialiseMonteCarlo(network.networkA, pow(10, inputData.startTemperature), logger, inputData.randomSeed);

    network.isOpenMP = inputData.isOpenMPIEnabled;
    network.MCWeighting = inputData.randomOrWeighted;
    network.isSimpleGraphene = inputData.isSimpleGrapheneEnabled;
    network.isTersoffGraphene = inputData.isTersoffGrapheneEnabled;
    network.isTriangleRaft = inputData.isTriangleRaftEnabled;
    network.isBilayer = inputData.isBilayerEnabled;
    network.isBN = inputData.isBNEnabled;
    network.spiralRadius = inputData.spiralRadius;
    if (inputData.isSpiralEnabled)
        network.makerFixed();
    network.MC_Routine = inputData.selectedMinimisationProtocol;

    logger->info("Network initialised");

    // Initialise output files
    logger->info("Initialising analysis output files...");
    int status;

    status = mkdir(inputData.outputFolder.c_str(), 0777);

    string prefixOut = inputData.outputFolder + "/" + inputData.outputFilePrefix;
    OutputFile outEnergyStats(prefixOut + "_e_compare.out");
    OutputFile outRingStats(prefixOut + "_ringstats.out");
    OutputFile outCorr(prefixOut + "_correlations.out");
    OutputFile outEnergy(prefixOut + "_energy.out");
    OutputFile outEntropy(prefixOut + "_entropy.out");
    OutputFile outTemperature(prefixOut + "_temperature.out");
    OutputFile outGeometry(prefixOut + "_geometry.out");
    OutputFile outEmatrix(prefixOut + "_ematrix.out");
    OutputFile outGeomHist(prefixOut + "_geomhist.out");
    OutputFile outAreas(prefixOut + "_areas.out");
    OutputFile outClusterA(prefixOut + "_cluster_a.out");
    OutputFile outClusterB(prefixOut + "_cluster_b.out");
    OutputFile outCndStats(prefixOut + "_cndstats.out");
    outGeometry.initVariables(6, 4, 60, 20);
    outAreas.initVariables(6, 4, 60, 30);
    outEmatrix.initVariables(1, 4, 60, int(log10(inputData.numRings * 12)) + 2);
    outGeomHist.initVariables(6, 4, 60, 20);
    outClusterB.initVariables(1, 4, 60, 10);
    logger->info("Ring statistics file created");
    logger->info("Correlations file created");
    logger->info("Energy file created");
    logger->info("Entropy file created");
    logger->info("Temperature file created");
    logger->info("Geometry file created");
    logger->info("Geometry histogram file created");
    logger->info("Edge distribution file created");
    logger->info("Coordination statistics file created");
    logger->info("Files initialised");

    bool lammps = true;

    // Initialise total analysis variables - only update in main simulation not equilibration
    VecF<double> lenHist(10000), angHist(10000);
    lenHist = 0.0;
    angHist = 0.0;

    // Run monte carlo
    logger->info("Running Monte Carlo...");
    logger->info("Initial energy: {}", network.mc.getEnergy());
    int accepted = 0, optIterations = 0;
    VecF<int> optCodes(5);
    optCodes = 0;
    int trackFreq = 100;
    VecF<int> moveStatus;

    // Run monte carlo thermalisation
    logger->info("Running Monte Carlo thermalisation");
    double energy = network.mc.getEnergy();
    logger->info("Initial energy: {}", energy);

    double SimpleGrapheneEnergy = 0.0;
    double TersoffGrapheneEnergy = 0.0;
    double TriangleRaftEnergy = 0.0;
    double BilayerEnergy = 0.0;
    double BNEnergy = 0.0;
    if (inputData.isSimpleGrapheneEnabled)
    {
        SimpleGrapheneEnergy = network.SimpleGraphene.GlobalPotentialEnergy();
    }
    if (inputData.isTersoffGrapheneEnabled)
    {
        TersoffGrapheneEnergy = network.TersoffGraphene.GlobalPotentialEnergy();
    }
    if (inputData.isTriangleRaftEnabled)
    {
        TriangleRaftEnergy = network.Triangle_Raft.GlobalPotentialEnergy();
    }
    if (inputData.isBilayerEnabled)
    {
        BilayerEnergy = network.Bilayer.GlobalPotentialEnergy();
    }
    if (inputData.isBNEnabled)
    {
        BNEnergy = network.BN.GlobalPotentialEnergy();
    }

    logger->info("Calculated energies for running systems");
    double mcT = pow(10, inputData.thermalisationTemperature);
    network.mc.setTemperature(mcT);
    if (inputData.isSpiralEnabled)
    {
        logger->info("Pre-run : attempting MC steps on all atoms within radius");
        bool disallowed_node;
        for (int i = 0; i < network.networkA.nodes.n; ++i)
        {
            if (network.rFixed[i] > 1000 / inputData.spiralRadius)
            {
                disallowed_node = false;
                for (int j = 0; j < network.networkB.nodes[network.fixedRings[0]].dualCnxs.n; ++j)
                {
                    if (i == network.networkB.nodes[network.fixedRings[0]].dualCnxs[j])
                        disallowed_node = true;
                }
                if (!disallowed_node)
                {
                    logger->info("Node {} withing radius at distance {}", i, network.rFixed[i]);
                    moveStatus = network.SpiralmonteCarloSwitchMoveLAMMPS(i, SimpleGrapheneEnergy, TersoffGrapheneEnergy,
                                                                          TriangleRaftEnergy, BilayerEnergy, BNEnergy, 0, logger);
                    VecF<double> ringStats = network.getNodeDistribution("B");
                    double r = network.getAssortativity("B");
                    double aEst = network.getAboavWeaireEstimate("B");
                    VecF<double> aw = network.getAboavWeaire("B");
                    VecF<double> s = network.getEntropy("B");
                    VecF<double> corr(6);
                    VecF<double> a(maxRingSize + 1), aSq(maxRingSize + 1);
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

                    VecF<double> emptyL, emptyA; // dummy histograms
                    VecF<double> geomStats = network.getOptimisationGeometry(network.networkA, emptyL, emptyA);
                    VecF<VecF<int>> edgeDist = network.getEdgeDistribution("B");
                    VecF<double> cndStats = network.getNodeDistribution("A");
                    outRingStats.writeRowVector(ringStats);
                    outCorr.writeRowVector(corr);
                    outEnergy.write(energy);
                    outEntropy.writeRowVector(s);
                    outTemperature.write(mcT);
                    outGeometry.writeRowVector(geomStats);
                    outAreas.writeRowVector(a);
                    outAreas.writeRowVector(aSq);

                    outEnergyStats.writeRowVector(eCompare);

                    for (int j = 0; j < edgeDist.n; ++j)
                        outEmatrix.writeRowVector(edgeDist[j]);
                    outCndStats.writeRowVector(cndStats);
                }
            }
        }
    }

    for (int i = 1; i <= inputData.initialThermalisationSteps; ++i)
    {
        logger->info("Equilibriating: {}  -- Energies: {} {} {} {} {}",
                     i, SimpleGrapheneEnergy, TersoffGrapheneEnergy, TriangleRaftEnergy, BilayerEnergy, BNEnergy);
        if (lammps)
        {
            moveStatus = network.monteCarloSwitchMoveLAMMPS(SimpleGrapheneEnergy, TersoffGrapheneEnergy,
                                                            TriangleRaftEnergy, BilayerEnergy, BNEnergy, 0, logger);
        }
        else
            moveStatus = network.monteCarloMixMove(energy);

        accepted += moveStatus[0];
        optCodes[moveStatus[1]] += 1;
        optIterations += moveStatus[2];
        if (i % inputData.analysisWriteFrequency == 0)
        {
            VecF<double> ringStats = network.getNodeDistribution("B");
            double r = network.getAssortativity("B");
            double aEst = network.getAboavWeaireEstimate("B");
            VecF<double> aw = network.getAboavWeaire("B");
            VecF<double> s = network.getEntropy("B");
            VecF<double> corr(6);
            VecF<double> a(maxRingSize + 1), aSq(maxRingSize + 1);
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

            VecF<double> emptyL, emptyA; // dummy histograms
            VecF<double> geomStats = network.getOptimisationGeometry(network.networkA, emptyL, emptyA);
            VecF<VecF<int>> edgeDist = network.getEdgeDistribution("B");
            VecF<double> cndStats = network.getNodeDistribution("A");
            outRingStats.writeRowVector(ringStats);
            outCorr.writeRowVector(corr);
            outEnergy.write(energy);
            outEntropy.writeRowVector(s);
            outTemperature.write(mcT);
            outGeometry.writeRowVector(geomStats);
            outAreas.writeRowVector(a);
            outAreas.writeRowVector(aSq);

            outEnergyStats.writeRowVector(eCompare);

            for (int j = 0; j < edgeDist.n; ++j)
                outEmatrix.writeRowVector(edgeDist[j]);
            outCndStats.writeRowVector(cndStats);
        }
        if (i % inputData.structureWriteFrequency == 0 && inputData.isWriteSamplingStructuresEnabled == 1)
        {
            network.syncCoordinates();
            std::string structureFilePath = prefixOut + "_therm_" + to_string(i);
            logger->info("Writing structure: {} to file {}", i, structureFilePath);
            network.writeXYZ(structureFilePath);
        }
    }
    logger->info("Monte Carlo equilibration complete");

    // Perform monte carlo simulation
    logger->info("Running Monte Carlo simulation");
    int nT = (inputData.endTemperature - inputData.startTemperature) / inputData.temperatureIncrement;

    for (int t = 0; t <= nT; ++t)
    {

        mcT = pow(10, inputData.startTemperature + t * inputData.temperatureIncrement);
        network.mc.setTemperature(mcT);
        logger->info("Temperature: {}", mcT);
        for (int i = 1; i <= inputData.stepsPerTemperature; ++i)
        {
            logger->info("Running: {}  -- Energies: {} {} {} {} {}",
                         i, SimpleGrapheneEnergy, TersoffGrapheneEnergy, TriangleRaftEnergy, BilayerEnergy, BNEnergy);
            if (lammps)
            {
                moveStatus = network.monteCarloSwitchMoveLAMMPS(SimpleGrapheneEnergy, TersoffGrapheneEnergy,
                                                                TriangleRaftEnergy, BilayerEnergy,
                                                                BNEnergy, 0, logger);
            }
            else
                moveStatus = network.monteCarloMixMove(energy);
            accepted += moveStatus[0];
            optCodes[moveStatus[1]] += 1;
            optIterations += moveStatus[2];
            if (i % inputData.analysisWriteFrequency == 0)
            {

                VecF<double> ringStats = network.getNodeDistribution("B");
                double r = network.getAssortativity("B");
                double aEst = network.getAboavWeaireEstimate("B");
                VecF<double> aw = network.getAboavWeaire("B");
                VecF<double> s = network.getEntropy("B");
                VecF<double> corr(6);
                VecF<double> a(maxRingSize + 1), aSq(maxRingSize + 1);
                double rr = network.getAssortativity("A");
                corr[0] = r;
                corr[1] = aEst;
                corr[2] = aw[0];
                corr[3] = aw[1];
                corr[4] = aw[2];
                corr[5] = rr;
                VecF<double> geomStats = network.getOptimisationGeometry(network.networkA, lenHist, angHist);
                VecF<VecF<int>> edgeDist = network.getEdgeDistribution("B");
                VecF<double> cndStats = network.getNodeDistribution("A");
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
                for (int j = 0; j < edgeDist.n; ++j)
                    outEmatrix.writeRowVector(edgeDist[j]);
                outCndStats.writeRowVector(cndStats);
            }
            if (i % inputData.structureWriteFrequency == 0 && inputData.isWriteSamplingStructuresEnabled == 1)
            {
                network.syncCoordinates();
                network.writeXYZ(prefixOut + "_t" + to_string(t) + "_" + to_string(i));
            }
        }
    }
    logger->info("Monte Carlo simulation complete");

    // Write total analysis
    for (int i = 0; i < 10000; ++i)
    {
        VecF<double> hist(4);
        hist[0] = i * 4.0 / 10000.0;
        hist[1] = lenHist[i];
        hist[2] = i * 2 * M_PI / 10000.0;
        hist[3] = angHist[i];
        outGeomHist.writeRowVector(hist);
    }

    if (inputData.isSimpleGrapheneEnabled)
    {
        logger->info("Writing Simple Graphene Results");
        network.SimpleGraphene.write_data(0);
        network.SimpleGraphene.write_restart(0);
    }
    if (inputData.isTriangleRaftEnabled)
    {
        logger->info("Writing Triangle Raft Results");
        network.Triangle_Raft.write_data(1);
        network.Triangle_Raft.write_restart(1);
    }
    if (inputData.isBNEnabled)
    {
        logger->info("Writing BN Results");
        network.BN.write_data(4);
        network.BN.write_restart(4);
    }

    // Check network
    logger->info("Simulation diagnostics");
    logger->info("Checking consistency and convexity");
    bool consistent = network.checkConsistency();
    bool convex = network.checkConvexity();
    logger->info("Network consistent:", consistent);
    logger->info("Rings convex:", convex);
    logger->info("Monte Carlo acceptance:", (double)accepted / inputData.stepsPerTemperature);
    logger->info("Geometry optimisation codes:");
    logger->info("Converged: ", optCodes[0]);
    logger->info("Converged (zero force): ", optCodes[1]);
    logger->info("Unconverged: ", optCodes[2]);
    logger->info("Failed (overlapping): ", optCodes[3]);
    logger->info("Failed (initially non-convex): ", optCodes[4]);
    logger->info("Geometry optimisation average iterations:", optIterations / vSum(optCodes));
    logger->info("Diagnostics complete");

    // Write files
    logger->info("Writing files");
    network.wrapCoordinates();
    network.syncCoordinates();
    network.write(prefixOut);
    network.writeXYZ(prefixOut);
    logger->info(prefixOut);
    logger->info("Writing network files");
    logger->info("Files written");

    // Close files
    // logFile.datetime("Simulation complete at: ");

    spdlog::shutdown();

    return 0;
}
