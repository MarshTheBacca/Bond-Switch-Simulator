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

LoggerPtr initialiseLogger(int argc, char *argv[]) {
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
    while ((opt = getopt(argc, argv, "d")) != -1) {
        if (opt == 'd') {
            logger->set_level(spdlog::level::debug);
            logger->debug("Debug messages enabled");
            break;
        }
    }

    return logger;
}
/**
 * @brief Attempts to switch the network at each temperature given in expTemperature
 */
void runSimulation(const std::vector<double> &expTemperatures, LinkedNetwork &linkedNetwork,
                   OutputFile &allStatsFile, const int &writeInterval, LoggerPtr logger) {
    for (int i = 1; i <= expTemperatures.size(); ++i) {
        linkedNetwork.metropolisCondition.setTemperature(expTemperatures[i]);
        logger->debug("Temperature: {:.2f}", expTemperatures[i]);
        linkedNetwork.monteCarloSwitchMoveLAMMPS();
        if (i % writeInterval == 0) {
            allStatsFile.writeValues(linkedNetwork.numSwitches, expTemperatures[i], linkedNetwork.energy,
                                     linkedNetwork.networkA.getEntropy(), linkedNetwork.networkA.getAssortativity(),
                                     linkedNetwork.networkA.getAboavWeaire(), linkedNetwork.getRingSizes(), linkedNetwork.getRingAreas());
        }
    }
}

int main(int argc, char *argv[]) {
    // Set start time so we can calculate the total run time
    auto start = std::chrono::high_resolution_clock::now();
    LoggerPtr logger;
    try {
        logger = initialiseLogger(argc, argv);
    } catch (std::exception &e) {
        std::cerr << "Exception while initialising logger: " << e.what() << std::endl;
        return 1;
    }
    try {
        logger->info("Network Monte Carlo");
        logger->info("Written by Marshall Hunt (Part II), Wilson Group, 2024");

        // Read input file
        InputData inputData("./netmc.inpt", logger);

        // Create output folder
        if (!std::filesystem::create_directory(inputData.outputFolder)) {
            logger->error("Error creating output folder {}", inputData.outputFolder);
            return 1;
        }

        // Initialise linkedNetwork
        logger->debug("Initialising linkedNetwork...");

        LinkedNetwork linkedNetwork;
        if (inputData.isFromScratchEnabled) {
            logger->debug("Creating linkedNetwork from scratch...");
            logger->debug("numRings: {}, minCoordination: {}, maxCoordination: {}, "
                          "minRingSize: {}, maxRingSize: {}",
                          inputData.numRings, inputData.minCoordination,
                          inputData.maxCoordination, inputData.minRingSize,
                          inputData.maxRingSize);
            // Generate hexagonal linkedNetwork from scratch
            linkedNetwork = LinkedNetwork(inputData.numRings, logger);
        } else {
            logger->debug("Loading linkedNetwork from files...");
            linkedNetwork = LinkedNetwork(inputData, logger);
        }

        logger->debug("Network initialised!");
        logger->info("Initial energy: {:.3f} Hartrees", linkedNetwork.energy);

        // Initialise output files
        logger->debug("Initialising analysis output file...");

        std::string prefixOut = inputData.outputFolder + "/" + inputData.outputFilePrefix;
        OutputFile allStatsFile(prefixOut + "_all_stats.csv");
        allStatsFile.writeDatetime("Written by LAMMPS-NetMC (Marshall Hunt, Wilson Group, 2024)");
        allStatsFile.writeLine("The data is structured as follows: Each value is comma separated, with inner vectors having their elements separated by semi-colons");
        std::string line = "For the ring size distribution, the sizes range from " + std::to_string(inputData.minRingSize) + " to " + std::to_string(inputData.maxRingSize);
        allStatsFile.writeLine(line);
        allStatsFile.writeLine("Step, Temperature, Energy, Entropy, Assortativity, Aboave Weaire, Ring Size Distribution (vector), Ring Areas (vector)");

        // Run monte carlo thermalisation
        std::vector<double>
            thermalisationTemperatures(inputData.initialThermalisationSteps, pow(10, inputData.thermalisationTemperature));
        logger->info("Thermalising...");
        runSimulation(thermalisationTemperatures, linkedNetwork, allStatsFile, inputData.analysisWriteInterval, logger);

        // Run monte carlo annealing

        std::vector<double> annealingTemperatures;
        int numTemperatureSteps = std::floor((inputData.endTemperature - inputData.startTemperature) / inputData.temperatureIncrement);
        annealingTemperatures.reserve((numTemperatureSteps + 1) * inputData.stepsPerTemperature);
        for (int i = 0; i <= numTemperatureSteps; ++i) {
            double temperature = inputData.startTemperature + i * inputData.temperatureIncrement;
            for (int j = 0; j < inputData.stepsPerTemperature; ++j) {
                annealingTemperatures.push_back(pow(10, temperature));
            }
        }
        logger->info("Annealing...");
        runSimulation(annealingTemperatures, linkedNetwork, allStatsFile, inputData.analysisWriteInterval, logger);
        logger->info("Simulation complete!");
        // Check linkedNetwork
        linkedNetwork.lammpsNetwork.stopMovie();
        logger->debug("Diagnosing simulation...");
        bool consistent = linkedNetwork.checkConsistency();
        logger->debug("Network consistent: {}", consistent ? "true" : "false");

        // Write files
        logger->debug("Writing files...");
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
        spdlog::shutdown();
    } catch (std::exception &e) {
        logger->error("Exception: {}", e.what());
        logger->flush();
        spdlog::shutdown();
        return 1;
    } catch (...) {
        logger->error("Unknown exception");
        logger->flush();
        spdlog::shutdown();
        return 1;
    }
    return 0;
}
