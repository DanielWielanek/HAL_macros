/*
 * fit3d.C
 *
 *  Created on: 5 cze 2024
 *      Author: daniel
 */

#ifndef __CLING__
#include "Header.h"
#include <TF3.h>
#endif

#include "randcf.C"

Double_t cff(Double_t* x, Double_t* p) {
  double scale = 5.06842372;
  double gaus  = 0;
  for (int i = 0; i < 3; i++)
    gaus += x[i] * x[i] * scale * scale * p[i + 2] * p[i + 2];
  return p[0] * (1. + p[1] * TMath::Exp(-gaus));
}


void fit3d() {
  auto cf  = GetRandomFunc(2.0, 2.5, 3.0);
  TF3* fit = new TF3("func", cff, 0, 1, 0, 1, 0, 1, 5);
  fit->SetParLimits(0, 0.5, 1.5);
  fit->SetParLimits(1, 0, 1);
  fit->SetParLimits(2, 1, 5);
  fit->SetParLimits(3, 1, 5);
  fit->SetParLimits(4, 1, 5);
  // clasical fit
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
  auto h3 = cf->GetHist(kFALSE);
  h3->Fit(fit);

  // hal fit
  auto hal = new Hal::CorrFit3DCF_Gauss();
  hal->SetMinimizedFunc(Hal::CorrFit::EMinFunc::kChi2);  // minimize "standard chi2"
  hal->SetMinimizer(Hal::CorrFit::EMinAlgo::kMinuitMigrad);
  hal->SetParLimits(hal->NormID(), 0.5, 1.5);
  hal->SetParLimits(hal->LambdaID(), 0, 1);
  hal->SetParLimits(hal->RoutID(), 1, 5);
  hal->SetParLimits(hal->RsideID(), 1, 5);
  hal->SetParLimits(hal->RlongID(), 1, 5);
  cf->Fit(hal);

  cout << "HAL chi2\t" << hal->GetChiSquare() << endl;
  cout << "HAL ndf\t" << hal->GetNDF() << endl;
  hal->Draw();

  auto hal2 = new Hal::CorrFit3DCF_Gauss();
  hal2->SetLineColor(kBlue);
  hal2->SetLineStyle(7);
  hal2->FixParameter(hal2->RoutID(), fit->GetParameter(2));
  hal2->FixParameter(hal2->RsideID(), fit->GetParameter(3));
  hal2->FixParameter(hal2->RlongID(), fit->GetParameter(4));
  hal2->FixParameter(hal2->LambdaID(), fit->GetParameter(1));
  hal2->FixParameter(hal2->NormID(), fit->GetParameter(0));
  cf->FitDummy(hal2);
  for (int i = 0; i < 5; i++)
    cout << i << " " << fit->GetParameter(i) << " " << hal2->NormID() << endl;

  hal2->Draw("same");
}
