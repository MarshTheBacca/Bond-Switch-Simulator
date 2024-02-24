#include "input_data.h"
#include "lammps_object.h"
#include "linked_network.h"
#include "output_file.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

LoggerPtr initialiseLogger() {
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("./netmc.log", true);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

    // Combine the sinks into a multi-sink logger
    auto logger = std::make_shared<spdlog::logger>("multi_sink", spdlog::sinks_init_list{file_sink, console_sink});
    spdlog::register_logger(logger);
    logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");

    // Set the default log level to info
    logger->set_level(spdlog::level::info);
    return logger;
}

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    LoggerPtr logger;
    try {
        logger = initialiseLogger();
    } catch (std::exception &e) {
        std::cerr << "Exception while initialising logger: " << e.what() << std::endl;
        return 1;
    }
    try {
        // Set start time so we can calculate the total run time

        // Create a file sink and a console sink with different names for clarity

        // Check command line arguments for --debug flag
        int opt;
        while ((opt = getopt(argc, argv, "d")) != -1) {
            if (opt == 'd') {
                logger->set_level(spdlog::level::debug);
                logger->debug("Debug messages enabled");
                break;
            }
        }

        logger->info("Network Monte Carlo");
        logger->info("Written by David Ormrod Morley, Edited by Oliver Whitaker and "
                     "Marshall Hunt, Wilson Group, 2024");

        // Read input file
        InputData inputData("./netmc.inpt", logger);

        // Create output folder
        if (!std::filesystem::create_directory(inputData.outputFolder)) {
            logger->error("Error creating output folder {}", inputData.outputFolder);
            return 1;
        }

        // Initialise linkedNetwork
        logger->info("Initialising linkedNetwork...");

        LinkedNetwork linkedNetwork;
        if (inputData.isFromScratchEnabled) {
            logger->info("Creating linkedNetwork from scratch...");
            logger->info("numRings: {}, minCoordination: {}, maxCoordination: {}, "
                         "minRingSize: {}, maxRingSize: {}",
                         inputData.numRings, inputData.minCoordination,
                         inputData.maxCoordination, inputData.minRingSize,
                         inputData.maxRingSize);
            // Generate hexagonal linkedNetwork from scratch
            linkedNetwork = LinkedNetwork(inputData.numRings, logger);
        } else {
            logger->info("Loading linkedNetwork from files...");
            linkedNetwork = LinkedNetwork(inputData, logger);
        }

        logger->info("Network initialised!");
        logger->info("Initial energy: {:.3f} Hartrees", linkedNetwork.energy);

        // Initialise output files
        logger->info("Initialising analysis output files...");

        std::string prefixOut = inputData.outputFolder + "/" + inputData.outputFilePrefix;
        OutputFile outRingStats(prefixOut + "_ringstats.out");
        OutputFile outCorr(prefixOut + "_correlations.out");
        OutputFile energyFile(prefixOut + "_energy.out");
        OutputFile outEntropy(prefixOut + "_entropy.out");
        OutputFile temperatureFile(prefixOut + "_temperature.out");
        OutputFile outEmatrix(prefixOut + "_ematrix.out");
        OutputFile outAreas(prefixOut + "_areas.out");
        OutputFile outCndStats(prefixOut + "_cndstats.out");

        // Run monte carlo thermalisation
        logger->info("Running Thermalisation...");
        double expTemperature = pow(10, inputData.thermalisationTemperature);
        linkedNetwork.mc.setTemperature(expTemperature);
        double timeSpentSwitching = 0;
        for (int i = 1; i <= inputData.initialThermalisationSteps; ++i) {
            std::chrono::high_resolution_clock::time_point startSwitch = std::chrono::high_resolution_clock::now();
            linkedNetwork.monteCarloSwitchMoveLAMMPS();
            std::chrono::high_resolution_clock::time_point endSwitch = std::chrono::high_resolution_clock::now();
            timeSpentSwitching += std::chrono::duration_cast<std::chrono::milliseconds>(endSwitch - startSwitch).count();
            if (i % inputData.analysisWriteInterval == 0) {
                std::vector<double> corr(6);
                double aboavWeaireEstimate = linkedNetwork.networkB.getAboavWeaireEstimate();
                std::vector<double> ringStats = linkedNetwork.networkB.getNodeDistribution();
                double networkBAssortativity = linkedNetwork.networkB.getAssortativity();
                std::tie(corr[2], corr[3], corr[4]) = linkedNetwork.networkB.getAboavWeaireParams();
                std::vector<double> entropy = linkedNetwork.networkB.getEntropy();
                std::vector<double> a(inputData.maxRingSize + 1);
                std::vector<double> aSq(inputData.maxRingSize + 1);
                double networkAAssortativity = linkedNetwork.networkA.getAssortativity();
                corr[0] = networkBAssortativity;
                corr[1] = aboavWeaireEstimate;
                corr[5] = networkAAssortativity;

                std::vector<double> emptyL;
                std::vector<double> emptyA; // dummy histograms
                std::vector<std::vector<int>> edgeDist = linkedNetwork.networkB.edgeDistribution;
                std::vector<double> cndStats = linkedNetwork.networkA.getNodeDistribution();
                outRingStats.writeVector(ringStats);
                outCorr.writeVector(corr);
                energyFile.writeValues(linkedNetwork.energy);
                outEntropy.writeVector(entropy);
                temperatureFile.writeValues(expTemperature);
                outAreas.writeVector(a);
                outAreas.writeVector(aSq);

                for (int j = 0; j < edgeDist.size(); ++j)
                    outEmatrix.writeVector(edgeDist[j]);
                outCndStats.writeVector(cndStats);
            }
        }
        logger->info("Thermalisation complete");

        // Perform monte carlo simulation
        logger->info("Annealing...");

        int numTemperatureSteps = std::abs(std::floor((inputData.endTemperature - inputData.startTemperature) / inputData.temperatureIncrement));
        logger->info("Number of temperature steps: {}", numTemperatureSteps);
        for (int i = 0; i < numTemperatureSteps; ++i) {
            expTemperature = pow(10, inputData.startTemperature + i * inputData.temperatureIncrement);
            linkedNetwork.mc.setTemperature(expTemperature);
            logger->info("Temperature: {:.2f}", expTemperature);
            for (int k = 1; k <= inputData.stepsPerTemperature; ++k) {
                std::chrono::high_resolution_clock::time_point startSwitch = std::chrono::high_resolution_clock::now();
                linkedNetwork.monteCarloSwitchMoveLAMMPS();
                std::chrono::high_resolution_clock::time_point endSwitch = std::chrono::high_resolution_clock::now();
                timeSpentSwitching += std::chrono::duration_cast<std::chrono::milliseconds>(endSwitch - startSwitch).count();
                if (k % inputData.analysisWriteInterval == 0) {
                    std::vector<double> corr(6);
                    std::vector<double> ringStats = linkedNetwork.networkB.getNodeDistribution();
                    double networkBAssortativity = linkedNetwork.networkB.getAssortativity();
                    double aboavWeaireEstimate = linkedNetwork.networkB.getAboavWeaireEstimate();
                    std::tie(corr[2], corr[3], corr[4]) = linkedNetwork.networkB.getAboavWeaireParams();
                    std::vector<double> entropy = linkedNetwork.networkB.getEntropy();
                    std::vector<double> a(inputData.maxRingSize + 1);
                    std::vector<double> aSq(inputData.maxRingSize + 1);
                    double networkAAssortativity = linkedNetwork.networkA.getAssortativity();
                    corr[0] = networkBAssortativity;
                    corr[1] = aboavWeaireEstimate;
                    corr[5] = networkAAssortativity;
                    std::vector<std::vector<int>> edgeDist = linkedNetwork.networkB.edgeDistribution;
                    std::vector<double> cndStats = linkedNetwork.networkA.getNodeDistribution();

                    energyFile.writeValues(linkedNetwork.energy);
                    temperatureFile.writeValues(expTemperature);

                    outRingStats.writeVector(ringStats);
                    outCorr.writeVector(corr);

                    outEntropy.writeVector(entropy);
                    outAreas.writeVector(a);
                    outAreas.writeVector(aSq);
                    for (int j = 0; j < edgeDist.size(); ++j)
                        outEmatrix.writeVector(edgeDist[j]);
                    outCndStats.writeVector(cndStats);
                }
            }
        }
        linkedNetwork.lammpsNetwork.stopMovie();
        logger->info("Annealing complete");

        // Check linkedNetwork
        logger->debug("Diagnosing simulation...");
        bool consistent = linkedNetwork.checkConsistency();
        logger->debug("Network consistent: {}", consistent ? "true" : "false");

        // Write files
        logger->info("Writing files...");
        linkedNetwork.write(prefixOut);
        logger->info("");
        logger->info("Number of attempted switches: {}", linkedNetwork.numSwitches);
        logger->info("Number of accepted switches: {}", linkedNetwork.numAcceptedSwitches);
        logger->info("Number of failed switches due to angle: {}", linkedNetwork.failedAngleChecks);
        logger->info("Number of failed switches due to bond length: {}", linkedNetwork.failedBondLengthChecks);
        logger->info("Number of failed switches due to energy: {}", linkedNetwork.failedEnergyChecks);
        logger->info("");
        logger->info("Monte Carlo acceptance: {:.3f}", (double)linkedNetwork.numAcceptedSwitches / linkedNetwork.numSwitches);
        if (linkedNetwork.checkAllClockwiseNeighbours()) {
            logger->info("All rings have clockwise neighbours");
        } else {
            logger->info("Not all rings have clockwise neighbours");
        }

        // Log time taken
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) / 1000.0;
        logger->info("Total run time: {:.3f} s", duration.count());
        logger->info("Time spent switching: {:.3f} s", timeSpentSwitching / 1000.0);
        spdlog::shutdown();
    } catch (std::exception &e) {
        logger->error("Exception: {}", e.what());
        logger->flush();
        spdlog::shutdown();
        return 1;
    }
    return 0;
}
