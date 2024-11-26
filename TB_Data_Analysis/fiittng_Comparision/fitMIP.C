#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

double langaufun(double *x, double *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[5]*TMath::Gaus(x[0],par[4],par[3]) +  par[2] * step * sum * invsq2pi / par[3]);
}

void fitMIP(){

  TFile *fin = TFile::Open("100GeV_Pion_without_abs_100K_out.root");
  //TCanvas *c = (TCanvas*)fin->Get("E0_L2");
  TH1F *h = (TH1F*)fin->Get("Layer 2 E0");

  /*TF1 *ffit = new TF1("f",langaufun,5,90,4);
  ffit->SetParameters(1,8,h->Integral(),1.5);
  ffit->SetParNames("Width","MP","Area","GSigma");
  */
  
  TF1 *ffit = new TF1("f",langaufun,10,60,6);
  ffit->SetParameters(1,8,h->Integral(),1.5,0,h->Integral());
  ffit->SetParNames("Width","MP","Area","GSigma", "Gmean", "Ginteg");

  
  /*
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  */
  
  h->Fit(ffit,"R");   // fit within specified range, use ParLimits, do not plot

  TCanvas *c1 = new TCanvas("c2","c",600,600);
  h->Draw();
  
}
