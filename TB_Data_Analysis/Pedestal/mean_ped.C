void mean_ped()
{
	TChain *Tout = new TChain("Events");

	for (int i = 0; i < 1; i++) {
            Tout->Add("/home/pravesh/Desktop/Cross/Ped_Run1695495152/input/Run1695495152.root");
            }

        TFile *fout = new TFile("/home/pravesh/Desktop/Cross/Ped_Run1695495152/Run1695495152_mean_pedestal.root","RECREATE");
	//TFile *fped = new TFile("/home/pravesh/Desktop/Without_absorber/Pedestal/Mean_Ped.root");
	
	
	//Defining tuples and addressing branches
        int N_HGC=444;
        
        Int_t           nHGC;
        UChar_t         HGC_halfrocChannel[N_HGC];
        UShort_t        HGC_adc[N_HGC];
        UChar_t         HGC_layer[N_HGC];
        UChar_t         HGC_econdeRx[N_HGC];
        Float_t         HGC_x[N_HGC];
        Float_t         HGC_y[N_HGC];
        
	Tout->SetBranchAddress("nHGC", &nHGC);
	Tout->SetBranchAddress("HGC_halfrocChannel", &HGC_halfrocChannel);
	Tout->SetBranchAddress("HGC_adc", &HGC_adc);
  Tout->SetBranchAddress("HGC_layer", &HGC_layer);     
  Tout->SetBranchAddress("HGC_econdeRx", &HGC_econdeRx); 
  Tout->SetBranchAddress("HGC_x", &HGC_x); 
  Tout->SetBranchAddress("HGC_y", &HGC_y);
  	
        int CHANNEL_MAX = 234;
        
        
        //Defining Histograms for Each channel for ADC
        
        TH1F *ADC_Hist_L1[CHANNEL_MAX],*ADC_Hist_L2[CHANNEL_MAX];
        char name[100];
   	char title[100];
	for(int u=0; u<CHANNEL_MAX; ++u){
		
		int bins=1000;
		int xmin=-20;
		int xmax=980;
		
      		sprintf(name,"layer1_ch_%i",u);
      		sprintf(title,"layer1_ch_%i",u);
      		ADC_Hist_L1[u] = new TH1F(name,title,bins,xmin,xmax);
		ADC_Hist_L1[u]->GetXaxis()->SetTitle("ADC");
		ADC_Hist_L1[u]->GetYaxis()->SetTitle("Events");
		
		sprintf(name,"layer2_ch_%i",u);
      		sprintf(title,"layer2_ch_%i",u);
      		ADC_Hist_L2[u] = new TH1F(name,title,bins,xmin,xmax);
		ADC_Hist_L2[u]->GetXaxis()->SetTitle("ADC");
		ADC_Hist_L2[u]->GetYaxis()->SetTitle("Events");
		
        	}	
       
	
	int nevt;
        nevt=Tout->GetEntries();
        
        //nevt=100;
        float n_ADC=0;
        int N_channel;
        
        
        
        for (int i = 0; i < nevt; i++)    //Event Loop Starts
        
        
	{
	Tout->GetEntry(i);
	for(int ij=0; ij<nHGC;++ij){
			if (HGC_adc[ij] == 0) continue;
                        if(HGC_x[ij]==-1 && HGC_y[ij]==-1) continue;
			if(HGC_layer[ij]==1){
				N_channel = 39*HGC_econdeRx[ij] + HGC_halfrocChannel[ij];
				n_ADC = HGC_adc[ij];
				ADC_Hist_L1[N_channel]->Fill(n_ADC);
				}
			
			
			if(HGC_layer[ij]==2){
				N_channel = 39*HGC_econdeRx[ij] + HGC_halfrocChannel[ij];
				n_ADC = HGC_adc[ij];
				ADC_Hist_L2[N_channel]->Fill(n_ADC);
				}
			}
				
	}
	
	std::ofstream outputFile("Mean_ADC.txt");
        
        outputFile <<"Channel	"<<"L1_Mean_ADC	"<<"L2_Mean_ADC	"<<std::endl;
	
	//Defining Fitting Function & Fitting Histograms of channels to get mean from gaussian and standard deviation from histogram not from gaussian

	TF1 *fitFunc = new TF1("fitFunc", "gaus", -1, -1);
	
	
	float        	 Ch[CHANNEL_MAX];
        float    	 m_L1[CHANNEL_MAX];
        float        	 m_L2[CHANNEL_MAX];
        float       	 sd_L1[CHANNEL_MAX];
        float       	 sd_L2[CHANNEL_MAX];
	
	for(int jk=0; jk<CHANNEL_MAX;++jk){
	
	ADC_Hist_L1[jk]->Fit(fitFunc, "Rq");
	gStyle->SetOptFit(1111);
	float mean1 = fitFunc->GetParameter(1);
    	float std1 = ADC_Hist_L1[jk]->GetStdDev();
    	
    	Long64_t nEntries = ADC_Hist_L1[jk]->GetEntries();
    	
    	if (nEntries==0){
    		mean1=0;
    		std1=0;
    		}
    	
	ADC_Hist_L1[jk]->Write();
	
	
	ADC_Hist_L2[jk]->Fit(fitFunc, "Rq");
	gStyle->SetOptFit(1111);
	float mean2 = fitFunc->GetParameter(1);
    	float std2 = ADC_Hist_L2[jk]->GetStdDev();
    	
    	nEntries = ADC_Hist_L2[jk]->GetEntries();
    	
    	if (nEntries==0){
    		mean2=0;
    		std2=0;
    		}

	ADC_Hist_L2[jk]->Write();
	
	Ch[jk]=jk;
	m_L1[jk] = mean1;
	m_L2[jk] = mean2;
	sd_L1[jk] = std1;
	sd_L2[jk] = std2;
	
	outputFile <<jk<<"	"<<mean1<<"	"<<mean2<<std::endl;
	}
	
	
	//Plotting Graphs of mean and standard deviation for each layer
	
	TGraph* mean_plot_l1 = new TGraph(CHANNEL_MAX, Ch, m_L1);
	mean_plot_l1->SetTitle("Mean vs Channel (Layer 1)");
	mean_plot_l1->GetXaxis()->SetTitle("Global Channel ID");
    	mean_plot_l1->GetYaxis()->SetTitle("Pedestal Mean (ADC)");
	mean_plot_l1->Write("Mean_L1");
    	
        TGraph* std_plot_l1 = new TGraph(CHANNEL_MAX, Ch, sd_L1);
	std_plot_l1->SetTitle("Std vs Channel (Layer 1)");
	std_plot_l1->GetXaxis()->SetTitle("Global Channel ID");
    	std_plot_l1->GetYaxis()->SetTitle("Pedestal Standard Deviation (ADC)");
	std_plot_l1->Write("std_L1");
  
        
        TGraph* mean_plot_l2 = new TGraph(CHANNEL_MAX, Ch, m_L2);
	mean_plot_l2->SetTitle("Mean vs Channel (Layer 2)");
	mean_plot_l2->GetXaxis()->SetTitle("Global Channel ID");
    	mean_plot_l2->GetYaxis()->SetTitle("Pedestal Mean (ADC)");
	mean_plot_l2->Write("Mean_L2");
    	
        TGraph* std_plot_l2 = new TGraph(CHANNEL_MAX, Ch, sd_L2);
	std_plot_l2->SetTitle("Std vs Channel (Layer 2)");
	std_plot_l2->GetXaxis()->SetTitle("Global Channel ID");
    	std_plot_l2->GetYaxis()->SetTitle("Pedestal Standard Deviation (ADC)");
	std_plot_l2->Write("std_L2");
    	 
}

