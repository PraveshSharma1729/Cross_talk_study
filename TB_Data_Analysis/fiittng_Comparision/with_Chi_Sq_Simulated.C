void with_Chi_Sq_Simulated()
{


	std::ofstream outputFile("Chi_Sq.txt");
	outputFile <<"Scaling_Factor	"<<"Chi_Sq	"<<endl;
	
for(int alp=0; alp<20;alp++){	
	
      

    // Open the input ROOT file
    TFile* file = TFile::Open("HGCTBSimOutput_MC_PiE100_33mBeam.root");

    // Get the TTree from the file
    TTree* Tout = (TTree*)file->Get("SimTree");

    // Create output file
    TFile* fout = new TFile("/home/pravesh/Desktop/Cross_Talk/Simulated_Output/100GeV_pi_without_abs_100K_out.root", "RECREATE");

    // Create histogram to fill
    TH1F *Cell_L1[19],*Cell_L2[19];
    char name[100];
    char title[100];
    char Energy[100];
    float xmin=0.0;
    float xmax=100.0;
    int bins=100000;
    
    for(int u=0;u<19;u++){
    	sprintf(name,"Layer 1 E%i",u);
    	sprintf(title,"Layer 1 E%i ",u);
    	Cell_L1[u]= new TH1F(name, title, 1400,-20,120);
    	Cell_L1[u]->GetXaxis()->SetTitle("ADC");
	Cell_L1[u]->GetYaxis()->SetTitle("Events");
	
	
	
	if(u==0){
	sprintf(name,"Layer 2 E%i",u);
    	sprintf(title,"Layer 2 E%i",u);
    	Cell_L2[u]= new TH1F(name, title, 100,-20,80);
    	Cell_L2[u]->GetXaxis()->SetTitle("ADC");
	Cell_L2[u]->GetYaxis()->SetTitle("Events");
	}
	
	else{
	sprintf(name,"Layer 2 E%i",u);
    	sprintf(title,"Layer 2 E%i",u);
    	Cell_L2[u]= new TH1F(name, title, 20,-10,10);
    	Cell_L2[u]->GetXaxis()->SetTitle("ADC");
	Cell_L2[u]->GetYaxis()->SetTitle("Events");
	}
	
    }
    
   
   
   
   

    // Variables to hold branch data
    UInt_t *nLayers = nullptr;
    vector<float> *E1Layer = nullptr;
    vector<float> *E7Layer = nullptr;
    vector<float> *E19Layer = nullptr;
    
    vector<vector<float>> Ener_l2;
    Ener_l2.clear();    
    for (int i = 0; i < 19; ++i) {
        // Create a new vector of doubles and add it to the vectorOfVectors
        vector<float> newVector;
        newVector.clear();
        Ener_l2.push_back(newVector);
    	}
    
    // Set branch address
    TBranch *b_nLayers;
    TBranch *b_E1Layer;
    TBranch *b_E7Layer;
    TBranch *b_E19Layer;
    
    
    Tout->SetBranchAddress("nLayers", &nLayers, &b_nLayers);
    Tout->SetBranchAddress("E1Layer", &E1Layer, &b_E1Layer);
    Tout->SetBranchAddress("E7Layer", &E7Layer, &b_E7Layer);
    Tout->SetBranchAddress("E19Layer", &E19Layer, &b_E19Layer);

    // Loop over events
    int nevt = Tout->GetEntries();
    
    // Define the Gaussian function
    float MIP = 17.5+alp*0.1;
    //float mean=0.0118;
    float mean=0.;
    float sigma=1.493;
        
    //TF1 *gaussian = new TF1("gaussian", "gaus(0)", -5.0, 5.0);
    //gaussian->SetParameters(1.0 / (sigma * TMath::Sqrt(2.0 * TMath::Pi())), mean, sigma);

    TRandom3 *gaussian= new TRandom3();
    
    
    
    TF1 *fitFunc = new TF1("fitFunc", "pol1", 20, 100);
    
    //nevt=100;
    for (int i = 0; i < nevt; ++i)                             //Event Loop Starts
    {
        // Get event
        Tout->GetEntry(i);
        // Process each element in the vector
        float value1_l1 = MIP*(*E1Layer)[0];
        float value7_l1 = MIP*(*E7Layer)[0];
        float value19_l1 = MIP*(*E19Layer)[0];
        
        
        //Defining Value to be filled in histogram as multiple of MIP
        float value1_l2 = MIP*(*E1Layer)[1];
        float value7_l2 = MIP*(*E7Layer)[1];
        float value19_l2 = MIP*(*E19Layer)[1];
    	
    	
    	//Defining Vectors for storing energy value in each cell
    	
    	vector<float> Ener_l1;
    	Ener_l1.clear();
    	
    	
    	
    	vector<float> gvalue;
    	gvalue.clear();
    	
    	
    	
    	//Energy in E0 Cell layer 1 
    	float E0_Energy_l1 = value1_l1;
    	Ener_l1.push_back(E0_Energy_l1);
    	
    	//Energy in 1st Ring Cells layer 1 
    	float E_1st_ring_Cell_l1 = (value7_l1-value1_l1)/6.0;
    	
    	for (int k=0; k<6;k++){
    	Ener_l1.push_back(E_1st_ring_Cell_l1);
    	}
    	//Energy in 2nd Ring Cells layer 1 
    	float E_2nd_ring_Cell_l1 = (value19_l1-value7_l1)/12.0;
    	
    	for (int k=0; k<12;k++){
    	Ener_l1.push_back(E_2nd_ring_Cell_l1);
    	}
    	
    	
    	//Filling Layer 1 Histograms with value
    	for(int u=0; u<19; u++){
        	Cell_L1[u]->Fill(Ener_l1[u]);
        }
        
        
        
        
        
        
        
        
        //Taking random values from Gaussian Distribution
    	for(int u=0;u<19;u++){
	  //float val = gaussian->GetRandom();
	  float val = gaussian->Gaus(mean,sigma);
	  
    	gvalue.push_back(val);
    	}
    	
    	
    	
    	
    	
        
        //Energy in E0 Cell layer 2 
    	//float E0_Energy_l2 = value1_l2-gvalue[0];
	float E0_Energy_l2 = value1_l2+gvalue[0];
    	Ener_l2[0].push_back(E0_Energy_l2);
    	
    	
    	
    	//Energy in 1st Ring Cells layer 2 
    	float E_1st_ring_Cell_l2 = (value7_l2-value1_l2)/6.0;
    	
    	for (int k=1; k<7;k++){
	  //Ener_l2[k].push_back(E_1st_ring_Cell_l2-gvalue[k]);
	  Ener_l2[k].push_back(gvalue[k]);
    	
    	//cout<<"U  "<<k<<"  val  "<<gvalue[k]<<endl;
    	
    	}
    	//Energy in 2nd Ring Cells layer 2
    	float E_2nd_ring_Cell_l2 = (value19_l2-value7_l2)/12.0;
    	
    	for (int k=7; k<19;k++){
	  //Ener_l2[k].push_back(E_2nd_ring_Cell_l2-gvalue[k]);
	  Ener_l2[k].push_back(gvalue[k]);
    	
    	}
    
    	
    	//Filling Layer 2 Histograms with value
    	for(int u=0; u<19; u++){
        	Cell_L2[u]->Fill(Ener_l2[u][i]);
      
        }
        
    	}                                                  //Event Loop Ends
    
    
    	
    
    for(int u=0;u<19;u++){
	
	Cell_L1[u]->Write();}
    
    
  
        
     

    // Write histogram to output file    
    for(int u=0;u<19;u++){
	
	Cell_L2[u]->Write();
	
	}
	
	

	
	
	
	
	
	
	
	
	
    // Clean up

    fout->Close();
    file->Close();

    // Delete the file objects
    delete fout;
    delete file;
    
    
    
    
    
    
    
     	// Open the ROOT files
    TFile* inFile1 = TFile::Open("/home/pravesh/Desktop/Cross_Talk/Run1695564190/Run1695564190.root");
    TFile* inFile2 = TFile::Open("/home/pravesh/Desktop/Cross_Talk/Simulated_Output/100GeV_pi_without_abs_100K_out.root");

    // Retrieve the histograms from each file
    TH1F* hist1 = (TH1F*)inFile1->Get("E0_L2");
    TH1F* hist2 = (TH1F*)inFile2->Get("Layer 2 E0");
    
     // Set statbox options
    gStyle->SetOptStat(1111); // Set to display all stats info

    // Create canvas
    //TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    
    // Normalize histograms
    hist1->Scale(1.0 / hist1->Integral());
    hist2->Scale(1.0 / hist2->Integral());
    
    
    
    // Get the number of bins
    int nBins = hist1->GetNbinsX();
    float chiSquare=0.0;
    int k=0;
    for(int i=0; i<nBins; i++){
    double observed = hist1->GetBinContent(i);
    double expected = hist2->GetBinContent(i);
    if (observed > 0){
    chiSquare += pow((observed - expected), 2) / observed ;
    k++;
    }
    }
    
    
    
    outputFile <<MIP<<"	"<<chiSquare<<"	"<<endl;
    
    
    
    
    
    
    inFile1->Close();
    inFile2->Close();
    
    delete inFile1;
    delete inFile2;
    
    }
    
    
    
    
    
   
    
    
    
    
    
    
    
}
