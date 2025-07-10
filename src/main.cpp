#include "input_data.h"
#include "linked_network.h"
#include "output_file.h"
#include "random_number_generator.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
#include "stats.h"
#include <atomic>
#include <ctime>
#include <expected>
#include <filesystem>
#include <iostream>
#include <signal.h>
#include <sys/stat.h>
#include <unistd.h>

// Global exit flag to cleanly exit when we use ctrl+c
static std::atomic<bool> exitFlag(false);
// Set start time so we can calculate the total run time
const auto start = std::chrono::high_resolution_clock::now();

/**
 * @brief Signal handler to set the exit flag to true when we use ctrl+c
 * @param sig The signal
 */
inline void exitFlagger([[maybe_unused]] const int sig) { exitFlag = true; }

void baseCleanup(const LoggerPtr &logger) {
  logger->flush();
  spdlog::shutdown();
  std::filesystem::remove("./log.lammps");
}

/**
 * @brief Cleans up the simulation by writing the network files and stopping the
 * LAMMPS movie
 * @param linkedNetwork The linked network to clean up
 */
void cleanup(LinkedNetwork &linkedNetwork, OutputFile &allStatsFile) {
  linkedNetwork.lammpsManager.stopMovie();
  linkedNetwork.write();
  linkedNetwork.lammpsManager.writeData();
  std::filesystem::remove("./log.lammps");
  allStatsFile.writeFooter(linkedNetwork.stats,
                           linkedNetwork.checkConsistency(), start);
  baseCleanup(linkedNetwork.logger);
}

/**
 * @brief Initialises the logger by creating a file sink and a console sink
 */
LoggerPtr initialiseLogger(const int argc, char *const *const argv) {
  // Create a file sink and a console sink with different names for clarity
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      std::filesystem::path("./output_files") / "bond_switch_simulator.log",
      true);
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

  // Combine the sinks into a multi-sink logger
  auto logger = std::make_shared<spdlog::logger>(
      "multi_sink", spdlog::sinks_init_list{file_sink, console_sink});
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
    }
  }
  return logger;
}

/**
 * @brief Attempts to switch the network at each temperature given in
 * expTemperature
 * @param expTemperatures The temperatures to switch at in raw form
 * @param linkedNetwork The linked network to switch
 * @param allStatsFile The file to write the statistics to
 * @param writeInterval The interval to write the statistics
 * @param logger The logger to log to
 */
void runSimulation(const std::vector<double> &expTemperatures,
                   LinkedNetwork &linkedNetwork, OutputFile &allStatsFile,
                   const int writeInterval, const LoggerPtr &logger) {
  if (expTemperatures.empty()) {
    logger->warn("No temperatures given, simulation not run");
    return;
  }
  double completion = 0.0;
  for (size_t i = 1; i <= expTemperatures.size(); ++i) {
    if (exitFlag) {
      logger->warn("Caught SIGINT, exiting...");
      cleanup(linkedNetwork, allStatsFile);
      exit(0);
    }
    linkedNetwork.performBondSwitch(expTemperatures[i - 1]);
    if (i % writeInterval == 0) {
      linkedNetwork.networkB.refreshStatistics();
      allStatsFile.writeValues(linkedNetwork.stats.getSwitches(),
                               expTemperatures[i - 1], linkedNetwork.energy,
                               linkedNetwork.networkB.entropy,
                               linkedNetwork.networkB.pearsonsCoeff,
                               linkedNetwork.networkA.getAboavWeaire(),
                               linkedNetwork.networkB.nodeSizes);
    }
    double currentCompletion =
        std::floor(static_cast<double>(i) /
                   static_cast<double>(expTemperatures.size()) / 0.1);
    if (currentCompletion > completion) {
      completion = currentCompletion;
      logger->info("{:.0f}% Complete", completion * 10);
    }
  }
}

std::expected<InputData, std::string> loadInputData(const LoggerPtr &logger) {
  // Read input file
  InputData inputData(
      std::filesystem::path("./input_files") / "bss_parameters.txt", logger);
  // Check if output folder already exists
  if (std::filesystem::exists("./output_files")) {
    logger->warn("Output folder already exists, files will be overwritten!");
  } else {
    // Try to create output folder
    if (!std::filesystem::create_directory("./output_files")) {
      return std::unexpected(
          std::format("Failed to create output directory: {}",
                      std::filesystem::current_path().string()));
    }
  }
  // Initialise Random Number Generator
  RandomNumberGenerator::initialize(inputData.randomSeed);
  return inputData;
}

std::expected<LinkedNetwork, std::string>
initialiseLinkedNetwork(const LoggerPtr &logger, const InputData &inputData) {
  // Initialise linkedNetwork
  logger->debug("Initialising linkedNetwork...");
  try {
    logger->debug("Loading linkedNetwork from files...");
    return LinkedNetwork(inputData, logger);
  } catch (const std::exception &e) {
    return std::unexpected(
        std::format("Failed to load linked network: {}", e.what()));
  }
}

std::expected<OutputFile, std::string>
initialiseOutputFile(const std::filesystem::path &path) {
  try {
    auto allStatsFile = OutputFile(path);
    allStatsFile.writeDatetime("Written by Bond-Switch-Simulator by Marshall "
                               "Hunt, Wilson Group, 2024");
    allStatsFile.writeLine(
        "The data is structured as follows: Each value is comma separated, "
        "with inner vectors having their elements separated by semi-colons");
    allStatsFile.writeLine(
        "Step, Temperature, Energy, Entropy, Pearson's Coefficient, Aboave "
        "Weaire, Ring Size Distribution (vector), Ring Areas (vector)");
    return allStatsFile;
  } catch (const std::exception &e) {
    return std::unexpected(
        std::format("Failed to initialise output file: {}", e.what()));
  }
}

void summarise(const LoggerPtr &logger, const Stats &stats,
               const bool &networkConsistent) {
  logger->info("");
  logger->info("Number of attempted switches: {}", stats.getSwitches());
  logger->info("Number of accepted switches: {}", stats.getAcceptedSwitches());
  logger->info("Number of failed switches due to angle: {}",
               stats.getFailedAngleChecks());
  logger->info("Number of failed switches due to bond length: {}",
               stats.getFailedBondLengthChecks());
  logger->info("Number of failed switches due to energy: {}",
               stats.getFailedEnergyChecks());
  logger->info("");
  logger->info("Monte Carlo acceptance: {:.3f}",
               (double)stats.getAcceptedSwitches() / stats.getSwitches());
  logger->info("Network consistent: {}", networkConsistent ? "true" : "false");
  logger->info("");
}

int main(const int argc, char *const *const argv) {
  // Set up signal handler to cleanly exit when we use ctrl+c
  signal(SIGINT, exitFlagger);
  LoggerPtr logger;
  try {
    logger = initialiseLogger(argc, argv);
  } catch (std::exception &e) {
    std::cerr << "Exception while initialising logger: " << e.what()
              << std::endl;
    return 1;
  }
  logger->info("Bond Switch Simulator");
  logger->info("Written by Marshall Hunt (Part II), Wilson Group, 2024");
  // Load input data
  logger->debug("Loading input data...");
  std::expected<InputData, std::string> inputDataResult = loadInputData(logger);
  if (!inputDataResult) {
    logger->error("Failed to load input data: {}", inputDataResult.error());
    baseCleanup(logger);
    return 1;
  }
  InputData inputData = std::move(inputDataResult.value());
  logger->info("Input data loaded successfully");

  logger->debug("Initialising linkedNetwork...");
  std::expected<LinkedNetwork, std::string> linkedNetworkResult =
      initialiseLinkedNetwork(logger, inputData);
  if (!linkedNetworkResult) {
    logger->error("Failed to initialise linked network: {}",
                  linkedNetworkResult.error());
    baseCleanup(logger);
    return 1;
  }
  LinkedNetwork linkedNetwork = std::move(linkedNetworkResult.value());
  logger->debug("Network initialised successfully");
  logger->info("Initial energy: {:.3f} Hartrees", linkedNetwork.energy);

  logger->debug("Initialising output file...");
  std::expected<OutputFile, std::string> outputFileResult =
      initialiseOutputFile(std::filesystem::path("./output_files") /
                           "bss_stats.csv");
  if (!outputFileResult) {
    logger->error("Failed to initialise output file: {}",
                  outputFileResult.error());
    baseCleanup(logger);
    return 1;
  }
  OutputFile allStatsFile = std::move(outputFileResult.value());
  logger->debug("Output file initialised successfully");

  try {

    // Run monte carlo thermalisation
    std::vector<double> thermalisationTemperatures(
        inputData.thermalisationSteps,
        pow(10, inputData.thermalisationTemperature));
    logger->info("Thermalising...");
    runSimulation(thermalisationTemperatures, linkedNetwork, allStatsFile,
                  inputData.analysisWriteInterval, logger);

    // Run monte carlo annealing
    std::vector<double> annealingTemperatures;
    annealingTemperatures.reserve(inputData.annealingSteps);
    double temperatureIncrement = (inputData.annealingEndTemperature -
                                   inputData.annealingStartTemperature) /
                                  (inputData.annealingSteps - 1);
    for (int i = 0; i < inputData.annealingSteps; ++i) {
      double temperature =
          inputData.annealingStartTemperature + i * temperatureIncrement;
      annealingTemperatures.push_back(pow(10, temperature));
    }

    logger->info("Annealing...");
    runSimulation(annealingTemperatures, linkedNetwork, allStatsFile,
                  inputData.analysisWriteInterval, logger);
    logger->info("Simulation complete!");
    linkedNetwork.lammpsManager.stopMovie();

    logger->debug("Writing final network files...");
    linkedNetwork.write();
    linkedNetwork.lammpsManager.writeData();
    bool networkConsistent = linkedNetwork.checkConsistency();
    summarise(logger, linkedNetwork.stats, networkConsistent);
    allStatsFile.writeFooter(linkedNetwork.stats, networkConsistent, start);

    // Log time taken
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start) /
        1000.0;
    logger->info("Total run time: {:.3f} s", duration.count());
    logger->info("Average time per step: {:.3f} us",
                 duration.count() / linkedNetwork.stats.getSwitches() * 1000.0);
    std::filesystem::remove("./log.lammps");
    logger->flush();
    spdlog::shutdown();
  } catch (std::exception &e) {
    logger->error("Exception: {}", e.what());
    cleanup(linkedNetwork, allStatsFile);
    return 1;
  } catch (...) {
    logger->error("Unknown exception");
    cleanup(linkedNetwork, allStatsFile);
    return 1;
  }
  return 0;
}
