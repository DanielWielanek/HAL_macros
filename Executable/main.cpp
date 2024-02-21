/*
 * main.cpp
 *
 *  Created on: 21 lut 2024
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */


#include <AnalysisManager.h>
#include <DataFormat.h>
#include <DbgTask.h>
#include <OTFMcEvent.h>
#include <OTFReader.h>
#include <OTFSource.h>
#include <PropertyMonitorXY.h>
#include <TH2.h>
#include <TrackAna.h>

typedef Hal::DataFieldID::Track::EBasic trackField;
int main() {
  auto run = new Hal::AnalysisManager();
  TH2D* h  = new TH2D("a", "a", 1, -0.5, 0.5, 1, 0, 1);
  h->SetBinContent(1, 1, 1);
  run->SetSource(new HalOTF::Source(100000));
  auto reader = new HalOTF::Reader();
  reader->SetSpiecies(*h, 211, 400);
  run->AddTask(reader);
  run->SetOutput("dataHS.root");
  auto dbg = new HalDbg::Task();
  dbg->SetReportStep(10000);
  dbg->SetLogFile("log.txt");
  run->AddTask(dbg);
  Hal::TrackAna* ana = new Hal::TrackAna();
  ana->SetFormat(new HalOTF::McEvent());
  Hal::TrackFieldMonitorXY kinMon(trackField::kEta, {100, -1, 1}, trackField::kPt, {100, 0, 2});
  ana->AddCutMonitor(kinMon);
  run->AddTask(ana);
  run->Init();
  run->Run();
  return 0;
}
