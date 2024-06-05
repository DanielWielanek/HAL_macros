/*
 * randcf.C
 *
 *  Created on: 5 cze 2024
 *      Author: daniel
 */
#ifndef __CLING__
#include "Header.h"
#endif

Hal::Femto3DCF* GetRandomFunc(double r_out, double r_side, double r_long) {
  int bins     = 25;
  double maxi  = 1.0;
  double mini  = 0.0;
  double pairs = 10000;
  auto cf      = new Hal::Femto3DCF("cf", bins, mini, maxi, bins, mini, maxi, bins, mini, maxi, Hal::Femto::EKinematics::kLCMS);
  auto num     = (TH3D*) cf->GetNum();
  auto den     = (TH3D*) cf->GetDen();
  r_out        = Hal::Femto::FmToGeV(r_out);
  r_side       = Hal::Femto::FmToGeV(r_side);
  r_long       = Hal::Femto::FmToGeV(r_long);

  for (int i = 1; i <= 100; i++) {
    double q_out = num->GetXaxis()->GetBinCenter(i);
    for (int j = 1; j <= 100; j++) {
      double q_side = num->GetYaxis()->GetBinCenter(j);
      for (int k = 1; k <= 100; k++) {
        double q_long  = num->GetZaxis()->GetBinCenter(k);
        double pairNum = gRandom->Gaus(pairs, TMath::Sqrt(pairs));
        double pairDen = gRandom->Gaus(pairs, TMath::Sqrt(pairs));
        Double_t cf =
          1.0
          + TMath::Exp(-(q_out * q_out * r_out * r_out + q_side * q_side * r_side * r_side + q_long * q_long * r_long * r_long));
        pairNum *= cf;
        num->SetBinContent(i, j, k, pairNum);
        num->SetBinError(i, j, k, TMath::Sqrt(pairNum));
        den->SetBinContent(i, j, k, pairDen);
        den->SetBinError(i, j, k, TMath::Sqrt(pairDen));
      }
    }
  }
  return cf;
}
