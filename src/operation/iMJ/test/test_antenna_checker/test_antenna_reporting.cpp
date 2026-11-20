// ***************************************************************************************
// Copyright (c) 2023-2025 Peng Cheng Laboratory
// Copyright (c) 2023-2025 Institute of Computing Technology, Chinese Academy of Sciences
// Copyright (c) 2023-2025 Beijing Institute of Open Source Chip
// ***************************************************************************************
#include "IdbDesign.h"
#include "IdbLayout.h"
#include "Logger.hpp"
#include "MJInterface.hpp"
#include "Utility.hpp"
#include "idm.h"

class AntennaReportingTest
{
 public:
  void run()
  {
    setenv("ECC_LOGGER_THROW_ON_ERROR", "1", 1);
    loadDesign();
    MJI.initMJ({{"-temp_directory_path", _output.string()}});
    testCheckOnly();
    testClockNetsAreChecked();
    testDiffusionThreshold();
    testCumRoutingPlusCut();
    testZeroDiffusionOutputIsNotProtected();
    testLayerQualifiedPinAreas();
    testCumulativeAreaSurvivesGatelessLayer();
    MJI.destroyMJ();
    dmInst->reset();
  }

 private:
  std::filesystem::path _data = std::filesystem::path(__FILE__).parent_path() / "data";
  std::filesystem::path _output = std::filesystem::temp_directory_path() / ("test_imj_antenna_" + std::to_string(getpid()));

  void loadDesign(const std::string& def = "design.def.in")
  {
    dmInst->reset();
    require(dmInst->readLef(std::vector<std::string>{(_data / "tech.lef").string()}, true), "read technology LEF");
    require(dmInst->readLef(std::vector<std::string>{(_data / "cell.lef").string()}), "read cell LEF");
    require(dmInst->readLef(std::vector<std::string>{(_data / "cell_layer_qualified.lef").string()}), "read layer-qualified cell LEF");
    require(dmInst->readDef((_data / def).string()), "read routed DEF");
  }

  void require(bool condition, const std::string& message)
  {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  void testCheckOnly()
  {
    idb::IdbDesign* design = dmInst->get_idb_design();
    idb::IdbNet* net = design->get_net_list()->find_net("violating");
    require(dmInst->saveDef((_output / "before.def").string()), "save input design");
    MJI.checkAntenna({{"-step_json", (_output / "antennaCheck.step.json").string()}});
    require(dmInst->get_idb_design() == design, "checker replaced the design");
    require(design->get_net_list()->find_net("violating") == net, "checker replaced a net");
    require(dmInst->saveDef((_output / "after.def").string()), "save checked design");
    require(readText(_output / "before.def") == readText(_output / "after.def"), "checker modified routed DEF");

    std::string report = readText(_output / "antenna_checker/antenna_check.rpt");
    require(report.find("Total Violations: 10\n") != std::string::npos, "wrong total violation count");
    require(report.find("Net: clean ") == std::string::npos, "clean component reported");
    require(report.find("Net: missing ") == std::string::npos, "missing gate data used as a gate");
    require(report.find("Net: power ") == std::string::npos, "power net checked as a signal");
    std::string step_json = readText(_output / "antennaCheck.step.json");
    require(step_json.find("\"schema_version\": 1") != std::string::npos, "missing antenna step schema");
    require(step_json.find("\"step\": \"antenna_check\"") != std::string::npos, "missing antenna step name");
    require(step_json.find("\"violation_count\": 10") != std::string::npos, "wrong structured violation count");
    require(step_json.find("\"net\": \"violating\"") != std::string::npos, "missing structured violation");
    require(report.find("Gate pin: driver/Y") == std::string::npos, "diffusion-only output reported as a gate");
    require(report.find("Gate pin: isolated/A") == std::string::npos, "disconnected gate selected");

    checkViolation(report, "M1", "PAR", 9.0, 2.0, 1.0);
    checkViolation(report, "M1", "CAR", 9.0, 3.0, 1.0);
    checkViolation(report, "M1", "PSR", 2.0, 1.0, 1.0);
    checkViolation(report, "M1", "CSR", 2.0, 1.5, 1.0);
    checkViolation(report, "V1", "CutPAR", 0.04, 0.01, 1.0);
    checkViolation(report, "V1", "CutCAR", 0.04, 0.02, 1.0);
    checkViolation(report, "M2", "PAR", 3.5, 2.0, 2.0);
    checkViolation(report, "M2", "CAR", 8.0, 7.0, 2.0);
    checkViolation(report, "M2", "PSR", 0.8, 0.5, 2.0);
    checkViolation(report, "M2", "CSR", 1.8, 1.5, 2.0);
    require(getViolation(report, "M1", "PAR").find("Metal area: 9 | Side area: 2 | Cut area: 0") != std::string::npos, "M1 area context");
    require(getViolation(report, "V1", "CutPAR").find("Cut area: 0.04") != std::string::npos, "cut area context");
    require(getViolation(report, "M2", "PAR").find("Metal area: 7 | Side area: 1.6 | Cut area: 0") != std::string::npos, "M2 area context");

    MJI.checkAntenna({});
    require(report == readText(_output / "antenna_checker/antenna_check.rpt"), "repeated check changed its result");
  }

  std::string readText(const std::filesystem::path& path)
  {
    std::ifstream* stream = MJUTIL.getInputFileStream(path.string());
    std::string text((std::istreambuf_iterator<char>(*stream)), std::istreambuf_iterator<char>());
    MJUTIL.closeFileStream(stream);
    return text;
  }

  void checkViolation(const std::string& report, const std::string& layer, const std::string& type, double ratio, double threshold, double gate_area)
  {
    std::string violation = getViolation(report, layer, type);
    require(imj::Utility::equalDoubleByError(getNumber(violation, "Ratio: "), ratio, 1e-9), "ratio mismatch: " + layer + " " + type);
    require(imj::Utility::equalDoubleByError(getNumber(violation, "Threshold: "), threshold, 1e-9), "threshold mismatch: " + layer + " " + type);
    require(imj::Utility::equalDoubleByError(getNumber(violation, "Gate area: "), gate_area, 1e-9), "gate area mismatch: " + layer + " " + type);
    require(violation.find("Gate pin: early/A\n") != std::string::npos, "representative gate mismatch: " + layer + " " + type);
  }

  std::string getViolation(const std::string& report, const std::string& layer, const std::string& type)
  {
    return getViolationOf(report, "violating", layer, type);
  }

  std::string getViolationOf(const std::string& report, const std::string& net, const std::string& layer, const std::string& type)
  {
    std::string header = "Net: " + net + " | Layer: " + layer + " | Type: " + type + "\n";
    size_t start = report.find(header);
    require(start != std::string::npos, "missing violation: " + net + " " + layer + " " + type);
    size_t end = report.find("\n\n", start);
    return report.substr(start, end - start);
  }

  double getNumber(const std::string& text, const std::string& key)
  {
    size_t start = text.find(key);
    require(start != std::string::npos, "missing numeric field: " + key);
    return std::stod(text.substr(start + key.size()));
  }

  void testClockNetsAreChecked()
  {
    // Same routed net, but declared USE CLOCK: it must be analysed, not silently skipped.
    idb::IdbNet* net = dmInst->get_idb_design()->get_net_list()->find_net("violating");
    net->set_connect_type(idb::IdbConnectType::kClock);
    MJI.checkAntenna({{"-step_json", (_output / "clock.step.json").string()}});
    std::string report = readText(_output / "antenna_checker/antenna_check.rpt");
    require(report.find("Total Violations: 10\n") != std::string::npos, "clock net was not antenna-checked");
    checkViolation(report, "M1", "PAR", 9.0, 2.0, 1.0);
    require(readText(_output / "clock.step.json").find("\"signal_net_count\": 3") != std::string::npos,
            "clock net missing from the analysed net count");
    net->set_connect_type(idb::IdbConnectType::kSignal);
  }

  void testDiffusionThreshold()
  {
    idb::IdbLayerRouting* layer = dynamic_cast<idb::IdbLayerRouting*>(dmInst->get_idb_layout()->get_layers()->find_layer("M1"));
    require(layer != nullptr, "missing M1");
    layer->set_antenna_diff_area_ratio_pwl({{0.0, 3.0}, {4.0, 11.0}});
    idb::IdbPin* pin = dmInst->get_idb_design()->get_instance_list()->find_instance("driver")->get_pin("Y");
    pin->get_term()->set_antenna_diff_area(2.0);
    MJI.checkAntenna({});
    std::string report = readText(_output / "antenna_checker/antenna_check.rpt");
    checkViolation(report, "M1", "DiffPAR", 9.0, 7.0, 1.0);
    require(getViolation(report, "M1", "DiffPAR").find("Diffusion area: 2") != std::string::npos, "diffusion area context");
    require(report.find("Layer: M1 | Type: PAR\n") == std::string::npos, "plain threshold used instead of diffusion PWL");

    // Equality with the threshold is not a violation.
    layer->set_antenna_diff_area_ratio_pwl({{0.0, 5.0}, {4.0, 13.0}});
    MJI.checkAntenna({});
    report = readText(_output / "antenna_checker/antenna_check.rpt");
    require(report.find("Layer: M1 | Type: DiffPAR\n") == std::string::npos, "threshold equality reported as a violation");
    require(report.find("Total Violations: 9\n") != std::string::npos, "threshold equality changed unrelated checks");
  }

  void testCumRoutingPlusCut()
  {
    // Enabling ANTENNACUMROUTINGPLUSCUT on M2 must add BOTH prior routing (M1=9) and
    // cut (V1=0.04) cumulative areas into M2 CAR: (9 + 0.04 + 7) / 2 = 8.02 (not 3.52).
    loadDesign();
    auto* m2 = dynamic_cast<idb::IdbLayerRouting*>(dmInst->get_idb_layout()->get_layers()->find_layer("M2"));
    auto* v1 = dynamic_cast<idb::IdbLayerCut*>(dmInst->get_idb_layout()->get_layers()->find_layer("V1"));
    require(m2 != nullptr && v1 != nullptr, "missing M2 or V1");
    m2->set_antenna_cum_routing_plus_cut(true);
    v1->set_antenna_cum_routing_plus_cut(true);
    MJI.checkAntenna({});
    std::string report = readText(_output / "antenna_checker/antenna_check.rpt");
    checkViolation(report, "V1", "CutCAR", 9.04, 0.02, 1.0);
    checkViolation(report, "M2", "CAR", 8.02, 7.0, 2.0);
  }

  void testZeroDiffusionOutputIsNotProtected()
  {
    // The driver's Y pin is an OUTPUT with "ANTENNADIFFAREA 0": no real diffusion, so the
    // plain threshold (2) applies, not the scalar diode threshold (5).
    loadDesign();
    idb::IdbLayerRouting* layer = dynamic_cast<idb::IdbLayerRouting*>(dmInst->get_idb_layout()->get_layers()->find_layer("M1"));
    require(layer != nullptr, "missing M1");
    layer->set_antenna_diff_area_ratio(5.0);
    MJI.checkAntenna({});
    std::string report = readText(_output / "antenna_checker/antenna_check.rpt");
    checkViolation(report, "M1", "PAR", 9.0, 2.0, 1.0);
    require(report.find("Layer: M1 | Type: DiffPAR\n") == std::string::npos, "zero-diffusion output selected the diode threshold");

    // Positive control: real diffusion area switches to the diode threshold.
    idb::IdbPin* pin = dmInst->get_idb_design()->get_instance_list()->find_instance("driver")->get_pin("Y");
    pin->get_term()->set_antenna_diff_area(2.0);
    MJI.checkAntenna({});
    report = readText(_output / "antenna_checker/antenna_check.rpt");
    checkViolation(report, "M1", "DiffPAR", 9.0, 5.0, 1.0);
  }

  void testLayerQualifiedPinAreas()
  {
    // GATE_Q/A declares "ANTENNAGATEAREA 1 LAYER M1", "ANTENNAGATEAREA 4 LAYER M2" and
    // "ANTENNADIFFAREA 3 LAYER M2" but the net only reaches the pin on M1. Only the M1 gate (1)
    // applies: PAR = 9 / 1 > 2. Summing every layer (gate 5) would hide the violation, and the
    // M2-only diffusion must not be credited.
    loadDesign("design_layer_gate.def.in");
    MJI.checkAntenna({});
    std::string report = readText(_output / "antenna_checker/antenna_check.rpt");
    std::string par = getViolationOf(report, "layered", "M1", "PAR");
    require(imj::Utility::equalDoubleByError(getNumber(par, "Ratio: "), 9.0, 1e-9), "M1 PAR ratio");
    require(imj::Utility::equalDoubleByError(getNumber(par, "Gate area: "), 1.0, 1e-9), "gate area summed across layers");
    require(par.find("Diffusion area: 0") != std::string::npos, "diffusion of another layer credited");
  }

  void testCumulativeAreaSurvivesGatelessLayer()
  {
    // The only gate sits on M2. The M1 wire (area 9) and the V1 via have no gate yet when they are
    // processed, but they still count in the cumulative M2 ratio: (9 + 7) / 1 = 16 > 7.
    loadDesign("design_late_gate.def.in");
    MJI.checkAntenna({});
    std::string report = readText(_output / "antenna_checker/antenna_check.rpt");
    std::string car = getViolationOf(report, "island", "M2", "CAR");
    require(imj::Utility::equalDoubleByError(getNumber(car, "Ratio: "), 16.0, 1e-9), "cumulative area lost before the gate joined");
    require(imj::Utility::equalDoubleByError(getNumber(car, "Threshold: "), 7.0, 1e-9), "M2 CAR threshold");
    // The per-layer ratio is unaffected by the gateless M1 island.
    require(imj::Utility::equalDoubleByError(getNumber(getViolationOf(report, "island", "M2", "PAR"), "Ratio: "), 7.0, 1e-9),
            "M2 PAR changed");
  }
};

int main()
{
  try {
    AntennaReportingTest().run();
  } catch (const std::exception& error) {
    MJLOG.warn(imj::Loc::current(), "Antenna checking regression failed: ", error.what());
    return 1;
  }
  return 0;
}
