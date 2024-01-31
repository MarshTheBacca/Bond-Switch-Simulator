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
#include <unistd.h>

int main(int argc, char *argv[])
{
    // Set start time so we can calculate the total run time
    auto start = std::chrono::high_resolution_clock::now();

    // Create a file sink and a console sink with different names for clarity
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("./netmc.log", true);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

    // Combine the sinks into a multi-sink logger
    auto logger = std::make_shared<spdlog::logger>("multi_sink", spdlog::sinks_init_list{file_sink, console_sink});
    spdlog::register_logger(logger);
    logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");

    // Set the default log level to info
    logger->set_level(spdlog::level::info);

    // Check command line arguments for --debug flag
    int opt;
    while ((opt = getopt(argc, argv, "d")) != -1)
    {
        if (opt == 'd')
        {
            logger->set_level(spdlog::level::debug);
            logger->debug("Debug messages enabled");
            break;
        }
    }

    logger->info("Network Monte Carlo");
    logger->info("Written by David Ormrod Morely, Edited by Oliver Whitaker and Marshall Hunt, Wilson Group, 2024");

    // Read input file
    InputData inputData("./netmc.inpt", logger);

    OutputFile lammpsBondAngle(inputData.inputFolder + "/PARM_Si.lammps");
    lammpsBondAngle.write("bond_style harmonic");
    lammpsBondAngle.write("bond_coeff 1 " + std::to_string(inputData.harmonicBondForceConstant) + " 1.000");
    lammpsBondAngle.write("angle_style cosine/squared");
    lammpsBondAngle.write("angle_coeff 1 " + std::to_string(inputData.harmonicAngleForceConstant) + " 120");

    // Initialise network
    logger->info("Initialising network");
    logger->debug("Number of rings: {}", inputData.numRings);
    logger->debug("Min ring size: {}", inputData.minRingSize);
    logger->debug("Max ring size: {}", inputData.maxRingSize);
    logger->debug("Min node coordination: {}", inputData.minCoordination);
    logger->debug("Max node coordination: {}", inputData.maxCoordination);
    logger->debug("Monte Carlo move type: {}", inputData.moveType);
    logger->debug("Monte carlo initial log temperature: {}", inputData.startTemperature);
    logger->debug("Monte carlo final log temperature: {}", inputData.endTemperature);
    logger->debug("Monte carlo log temperature increment: {}", inputData.temperatureIncrement);
    logger->debug("Monte carlo steps per increment: {}", inputData.stepsPerTemperature);
    logger->debug("Equilibrium steps: {}", inputData.initialThermalisationSteps);

    LinkedNetwork network(inputData.inputFolder, inputData.inputFilePrefix,
                          inputData.outputFolder,
                          inputData.minCoordination, inputData.maxCoordination,
                          inputData.minRingSize, inputData.maxRingSize,
                          inputData.isSimpleGrapheneEnabled,
                          inputData.isTriangleRaftEnabled, inputData.isBilayerEnabled,
                          inputData.isTersoffGrapheneEnabled,
                          inputData.isBNEnabled,
                          inputData.isRestartUsingLAMMPSObjectsEnabled, logger);

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

    network.isOpenMPIEnabled = inputData.isOpenMPIEnabled;
    network.MCWeighting = inputData.randomOrWeighted;
    network.isSimpleGrapheneEnabled = inputData.isSimpleGrapheneEnabled;
    network.isTersoffGrapheneEnabled = inputData.isTersoffGrapheneEnabled;
    network.isTriangleRaftEnabled = inputData.isTriangleRaftEnabled;
    network.isBilayerEnabled = inputData.isBilayerEnabled;
    network.isBNEnabled = inputData.isBNEnabled;
    network.spiralRadius = inputData.spiralRadius;
    if (inputData.isSpiralEnabled)
        network.makerFixed();
    network.MC_Routine = inputData.selectedMinimisationProtocol;

    logger->info("Network initialised!");

    // Initialise output files
    logger->info("Initialising analysis output files...");
    int status;

    // Create output folder
    status = mkdir(inputData.outputFolder.c_str(), 0777);
    if (status == -1)
    {
        logger->error("Error creating output folder");
        return 1;
    }

    std::string prefixOut = inputData.outputFolder + "/" + inputData.outputFilePrefix;
    OutputFile outEnergyStats(prefixOut + "_e_compare.out");
    logger->debug("Energy file created");
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

    logger->info("Entropy file created");
    logger->info("Temperature file created");
    logger->info("Geometry file created");
    logger->info("Geometry histogram file created");
    logger->info("Edge distribution file created");
    logger->info("Coordination statistics file created");
    logger->info("Files initialised");

    bool lammps = true;

    // Initialise total analysis variables - only update in main simulation not equilibration
    VecF<double> lenHist(10000);
    VecF<double> angHist(10000);
    lenHist = 0.0;
    angHist = 0.0;

    // Run monte carlo
    logger->info("Running Monte Carlo...");
    logger->info("Initial energy: {}", network.mc.getEnergy());
    int numAcceptedMoves = 0;
    int optIterations = 0;
    VecF<int> optCodes(5);
    optCodes = 0;
    VecF<int> moveStatus;

    // Run monte carlo thermalisation
    logger->info("Running Monte Carlo thermalisation...");
    double energy = network.mc.getEnergy();
    logger->info("Initial energy: {}", energy);
    try
    {

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
                                                                              TriangleRaftEnergy, BilayerEnergy, BNEnergy, logger);
                        VecF<double> ringStats = network.getNodeDistribution("B");
                        double r = network.getAssortativity("B");
                        double aEst = network.getAboavWeaireEstimate("B");
                        VecF<double> aw = network.getAboavWeaire("B");
                        VecF<double> s = network.getEntropy("B");
                        VecF<double> corr(6);
                        VecF<double> a(maxRingSize + 1);
                        VecF<double> aSq(maxRingSize + 1);
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

                        VecF<double> emptyL;
                        VecF<double> emptyA; // dummy histograms
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
                                                                TriangleRaftEnergy, BilayerEnergy, BNEnergy, logger);
            }
            else
                moveStatus = network.monteCarloMixMove(energy);

            if (moveStatus[0])
                numAcceptedMoves++;
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
                VecF<double> a(maxRingSize + 1);
                VecF<double> aSq(maxRingSize + 1);
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

                VecF<double> emptyL;
                VecF<double> emptyA; // dummy histograms
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
                std::string structureFilePath = prefixOut + "_therm_" + std::to_string(i);
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
                                                                    BNEnergy, logger);
                }
                else
                    moveStatus = network.monteCarloMixMove(energy);
                if (moveStatus[0])
                    numAcceptedMoves++;
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
                    VecF<double> a(maxRingSize + 1);
                    VecF<double> aSq(maxRingSize + 1);
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
                    network.writeXYZ(prefixOut + "_t" + std::to_string(t) + "_" + std::to_string(i));
                }
            }
        }
    }
    catch (std::exception &e)
    {
        logger->error("Exception: {}", e.what());
        return 1;
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
        network.SimpleGraphene.write_restart("Si");
    }
    if (inputData.isTriangleRaftEnabled)
    {
        logger->info("Writing Triangle Raft Results");
        network.Triangle_Raft.write_data(1);
        network.Triangle_Raft.write_restart("Si2O3");
    }
    if (inputData.isBNEnabled)
    {
        logger->info("Writing BN Results");
        network.BN.write_data(4);
        network.BN.write_restart("BN");
    }

    // Check network
    logger->debug("Diagnosing simulation...");
    logger->debug("Checking consistency and convexity");
    bool consistent = network.checkConsistency();
    bool convex = network.checkConvexity();
    logger->debug("Network consistent: {}", consistent);
    logger->debug("Rings convex: {}", convex);
    logger->debug("Monte Carlo acceptance: {}", (double)numAcceptedMoves / inputData.stepsPerTemperature);
    logger->debug("Geometry optimisation codes: {}");
    logger->debug("Converged: {}", optCodes[0]);
    logger->debug("Converged (zero force): {}", optCodes[1]);
    logger->debug("Unconverged: {}", optCodes[2]);
    logger->debug("Failed (overlapping): {}", optCodes[3]);
    logger->debug("Failed (initially non-convex): {}", optCodes[4]);
    logger->debug("Geometry optimisation average iterations: {}", optIterations / vSum(optCodes));

    // Write files
    logger->info("Writing files...");
    network.wrapCoordinates();
    network.syncCoordinates();
    network.write(prefixOut);
    network.writeXYZ(prefixOut);
    logger->info(prefixOut);
    logger->debug("Files written");

    // Log time taken
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) / 1000.0;
    logger->info("Total run time: {} s", duration.count());
    spdlog::shutdown();

    return 0;
}
