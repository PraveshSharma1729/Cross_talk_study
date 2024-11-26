void cross_talk()
{
	TChain *Tout = new TChain("Events");

	for(int i=0;i<3;i++){
	Tout->Add(Form("/home/pravesh/Desktop/Cross/Run1695564190/input/Run1695564190_Run1695564190_Link1_File000000000%i_NANO.root",i));
	}

        TFile *fout = new TFile("/home/pravesh/Desktop/Cross/Run1695564190/Run1695564190.root","RECREATE");
	TFile *fped = new TFile("/home/pravesh/Desktop/Cross/Ped_Run1695495152/Run1695495152_mean_pedestal.root");
	//Defining tuples and addressing branches
        int N_HGC=444;
        
        Int_t           nHGC;
        UChar_t         HGC_halfrocChannel[N_HGC];
        UShort_t        HGC_adc[N_HGC];
        UChar_t         HGC_layer[N_HGC];
        UChar_t         HGC_econdeRx[N_HGC];
        Float_t         HGC_x[N_HGC];
        Float_t         HGC_y[N_HGC];
        UInt_t          hgcMetadata_trigTime;
  	
  	 // List of branches
        TBranch        *b_nHGC;
        TBranch        *b_hgcMetadata_trigTime;
        TBranch        *b_HGC_halfrocChannel;
        TBranch        *b_HGC_adc;
        TBranch        *b_HGC_layer;
        TBranch        *b_HGC_econdeRx;
        TBranch        *b_HGC_x;
        TBranch        *b_HGC_y;

        Tout->SetBranchAddress("nHGC", &nHGC, &b_nHGC);
        Tout->SetBranchAddress("hgcMetadata_trigTime", &hgcMetadata_trigTime, &b_hgcMetadata_trigTime);
        Tout->SetBranchAddress("HGC_halfrocChannel", HGC_halfrocChannel, &b_HGC_halfrocChannel);
        Tout->SetBranchAddress("HGC_adc", HGC_adc, &b_HGC_adc);
        Tout->SetBranchAddress("HGC_econdeRx", HGC_econdeRx, &b_HGC_econdeRx);
        Tout->SetBranchAddress("HGC_x", HGC_x, &b_HGC_x);
        Tout->SetBranchAddress("HGC_y", HGC_y, &b_HGC_y);
        Tout->SetBranchAddress("HGC_layer", HGC_layer, &b_HGC_layer);

  	
        int CHANNEL_MAX = 234;
        
        
        //Defining Histograms for Each channel for ADC
        
        TH1F *ADC_Hist_L1[CHANNEL_MAX],*ADC_Hist_L2[CHANNEL_MAX];
        char name[100];
   	char title[100];
	for(int u=0; u<CHANNEL_MAX; ++u){
		
		int bins=100;
		int xmin=-20;
		int xmax=80;
		
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
       
       
       //Defining 2-D histograms for Trig Time vs ADC counts and also Defing ProfileX plots of those Channels
	TH2F *h_trigTime_l1[CHANNEL_MAX], *h_trigTime_l2[CHANNEL_MAX];
	TProfile *profileX_L1[CHANNEL_MAX],*profileX_L2[CHANNEL_MAX];
        for(int i=0; i<CHANNEL_MAX; ++i){
               sprintf(name,"TrigTime plot ch%i L1",i);
               sprintf(title,"ADC vs Tigger Time Layer1_ch%i",i);
               h_trigTime_l1[i] = new TH2F(name,title,32,104,136,120,-20,100);
               h_trigTime_l1[i]->GetXaxis()->SetTitle("Trig Time");
               h_trigTime_l1[i]->GetYaxis()->SetTitle("HGC_adc");  
               
               
               sprintf(name,"TrigTime plot ch%i L2",i);
               sprintf(title,"ADC vs Tigger Time Layer2_ch%i",i);
               h_trigTime_l2[i] = new TH2F(name,title,32,104,136,120,-20,100);
               h_trigTime_l2[i]->GetXaxis()->SetTitle("Trig Time");
               h_trigTime_l2[i]->GetYaxis()->SetTitle("HGC_adc");  
               
               sprintf(name,"profileX_Ch_%i_L1",i);
               sprintf(title,"profileX_Ch_%i_L1",i);
               profileX_L1[i] = new TProfile(name,title,15,105,120,-20,100,"p"); 
               
               
               sprintf(name,"profileX_Ch_%i_L2",i);
               sprintf(title,"profileX_Ch_%i_L2",i);
               profileX_L2[i] = new TProfile(name,title,15,105,120,-20,100,"p"); 
        }  
       
       //Defining Seed channel Histogram
       TH1F *Seed_L1, *Seed_L2;
       Seed_L1 = new TH1F("Seed L1","Seed Channel L1",234,0,234);
       Seed_L1->GetXaxis()->SetTitle("Channel No.");
       Seed_L1->GetYaxis()->SetTitle("Events");
       
       
       
       
       Seed_L2 = new TH1F("Seed L2","Seed Channel L2",234,0,234);
       Seed_L2->GetXaxis()->SetTitle("Channel No.");
       Seed_L2->GetYaxis()->SetTitle("Events");
       
       
       
       
       //Getting the mean pedestal values
	TGraphErrors *gr_l1 = (TGraphErrors*)fped->Get("Mean_L1");
	TGraphErrors *gr_l2 = (TGraphErrors*)fped->Get("Mean_L2");
	TGraphErrors *sd_l1 = (TGraphErrors*)fped->Get("std_L1");
	TGraphErrors *sd_l2 = (TGraphErrors*)fped->Get("std_L2");
	
        float ped_l1[CHANNEL_MAX],ped_l2[CHANNEL_MAX],std_l1[CHANNEL_MAX],std_l2[CHANNEL_MAX];
        for(int i=0; i< CHANNEL_MAX; ++i){
                ped_l1[i] = gr_l1->GetPointY(i);
                ped_l2[i] = gr_l2->GetPointY(i);
                std_l1[i] = sd_l1->GetPointY(i);
                std_l2[i] = sd_l2->GetPointY(i);
             
        }
        
        
        
     
	
	int nevt;
        nevt=Tout->GetEntries();
        //nevt=10000;
        
        vector<int>Highest_ADC_channel_L1;
        vector<int>Highest_ADC_channel_L2;
	
	Highest_ADC_channel_L1.clear();
        Highest_ADC_channel_L2.clear();
        
        float n_ADC=0;
        int N_channel;
        
        
        vector<float> Com1_L1;
        Com1_L1.clear();
        
        vector<float> Com2_L1;
        Com2_L1.clear();
        
        vector<float> Com3_L1;
        Com3_L1.clear();
        
        vector<float> Com1_L2;
        Com1_L2.clear();
        
        vector<float> Com2_L2;
        Com2_L2.clear();
        
        vector<float> Com3_L2;
        Com3_L2.clear();
     
        
        vector<vector<float>> trigg_time_l2;
        trigg_time_l2.clear();
        
        vector<vector<float>>Ch_ADC_L1;
        vector<vector<float>>Ch_ADC_L2;
        Ch_ADC_L1.clear();
        Ch_ADC_L2.clear();
        
        for (int i = 0; i < CHANNEL_MAX; ++i) {
        // Create a new vector of doubles and add it to the vectorOfVectors
        vector<float> newVector;
        newVector.clear();
        Ch_ADC_L1.push_back(newVector);
        Ch_ADC_L2.push_back(newVector);
        
        
        
        trigg_time_l2.push_back(newVector);
    	}
        
        int t=0;
        for (int i = 0; i < nevt; i++)    //Event Loop Starts
        
	{
	Tout->GetEntry(i);
	
	
	/*if(i>=231340 && i<=231344){
				
	cout<<"E  "<<i-231340<<endl;}*/
	
	//Common Mode Calculations per HGCROC
	float CM1_L1=0.0;
	float CM2_L1=0.0;
	float CM3_L1=0.0;
	float CM1_L2=0.0;
	float CM2_L2=0.0;
	float CM3_L2=0.0;
	
	int count1_L1 = 0;
        int count2_L1 = 0;
        int count3_L1 = 0;
        int count1_L2 = 0;
        int count2_L2 = 0;
        int count3_L2 = 0;
	
	
	for(int ik=0; ik<nHGC;++ik){
			//if (HGC_adc[ik] == 0) continue;
                        if(HGC_x[ik]==-1 && HGC_y[ik]==-1) continue;
			if(HGC_layer[ik]==1){
				N_channel = 39*HGC_econdeRx[ik] + HGC_halfrocChannel[ik];
				if(N_channel == 50 ||N_channel == 148 ||N_channel == 152) continue;
				if(HGC_econdeRx[ik]==0 || HGC_econdeRx[ik]==1){
					if(fabs(HGC_adc[ik]-ped_l1[N_channel])<3.0*std_l1[N_channel]){
					CM1_L1 = CM1_L1+HGC_adc[ik]-ped_l1[N_channel];
					count1_L1 = count1_L1+1;
					}
				}
				
                                if(HGC_econdeRx[ik]==2 || HGC_econdeRx[ik]==3){
                                	if(fabs(HGC_adc[ik]-ped_l1[N_channel])<3.0*std_l1[N_channel]){
                                	CM2_L1 = CM2_L1+HGC_adc[ik]-ped_l1[N_channel];
					count2_L1 = count2_L1+1;
                                	}
                                }
                                
                                if(HGC_econdeRx[ik]==4 || HGC_econdeRx[ik]==5){
                                	if(fabs(HGC_adc[ik]-ped_l1[N_channel])<3.0*std_l1[N_channel]){
                                	CM3_L1 = CM3_L1+HGC_adc[ik]-ped_l1[N_channel];
					count3_L1 = count3_L1+1;
                                	}
                                }
                                
                      	}
			
			
			if(HGC_layer[ik]==2){
				N_channel = 39*HGC_econdeRx[ik] + HGC_halfrocChannel[ik];
				
				if(HGC_econdeRx[ik]==0 || HGC_econdeRx[ik]==1){
					if(fabs(HGC_adc[ik]-ped_l2[N_channel])<3.0*std_l2[N_channel]){
					CM1_L2 = CM1_L2+HGC_adc[ik]-ped_l2[N_channel];
					count1_L2 = count1_L2+1;
					}
				}
				
                                if(HGC_econdeRx[ik]==2 || HGC_econdeRx[ik]==3){
                                	if(fabs(HGC_adc[ik]-ped_l2[N_channel])<3.0*std_l2[N_channel]){
                                	CM2_L2 = CM2_L2+HGC_adc[ik]-ped_l2[N_channel];
					count2_L2 = count2_L2+1;
                                	}
                                }
                                
                                if(HGC_econdeRx[ik]==4 || HGC_econdeRx[ik]==5){
                                	if(fabs(HGC_adc[ik]-ped_l2[N_channel])<3.0*std_l2[N_channel]){
                                	CM3_L2 = CM3_L2+HGC_adc[ik]-ped_l2[N_channel];
					count3_L2 = count3_L2+1;
                                	}
                               }
                                
			}
			
			}
			
		CM1_L1 = CM1_L1/count1_L1;
		CM2_L1 = CM2_L1/count2_L1;
		CM3_L1 = CM3_L1/count3_L1;
		CM1_L2 = CM1_L2/count1_L2;
		CM2_L2 = CM2_L2/count2_L2;
		CM3_L2 = CM3_L2/count3_L2;
		
		
		
		Com1_L1.push_back(CM1_L1);
		Com2_L1.push_back(CM2_L1);
		Com3_L1.push_back(CM3_L1);
		Com1_L2.push_back(CM1_L2);
		Com2_L2.push_back(CM2_L2);
		Com3_L2.push_back(CM3_L2);
	
	
	
	
	//Filling the histograms with pedestal and CM subtracted ADC and finding seed Channel
	
	
	float Highest_ADC_L1=-999;
	float Highest_ADC_L2=-999;
	int seed1=-1;
	int seed2=-1;
	
	
	for(int ij=0; ij<nHGC;++ij){
			//if (HGC_adc[ij] == 0) continue;
                        if(HGC_x[ij]==-1 && HGC_y[ij]==-1) continue;
                        
			if(HGC_layer[ij]==1){
			
				N_channel = 39*HGC_econdeRx[ij] + HGC_halfrocChannel[ij];
				if(N_channel == 50 ||N_channel == 148 ||N_channel == 152) continue;
				if(HGC_econdeRx[ij]==0 || HGC_econdeRx[ij]==1){
				n_ADC = HGC_adc[ij]-ped_l1[N_channel]-CM1_L1;
				}
				
                                if(HGC_econdeRx[ij]==2 || HGC_econdeRx[ij]==3){
                                n_ADC = HGC_adc[ij]-ped_l1[N_channel]-CM2_L1;
                                }
                                
                                if(HGC_econdeRx[ij]==4 || HGC_econdeRx[ij]==5){
                                n_ADC = HGC_adc[ij]-ped_l1[N_channel]-CM3_L1;
                                }
                         
                                Ch_ADC_L1[N_channel].push_back(n_ADC);   // Filling vectors of each channel of Layer 1 alongside, for further calculation 
				
				
				
				if (Highest_ADC_L1<n_ADC){
				Highest_ADC_L1=n_ADC;
				seed1=N_channel;
				}
				
				
				if(hgcMetadata_trigTime>=118 && hgcMetadata_trigTime<=124){
				
				h_trigTime_l1[N_channel]->Fill(hgcMetadata_trigTime,n_ADC);
				
				profileX_L1[N_channel]->Fill(hgcMetadata_trigTime,n_ADC);
				ADC_Hist_L1[N_channel]->Fill(n_ADC);
				}
				
				
				
				
				}
			
			
			if(HGC_layer[ij]==2){
				N_channel = 39*HGC_econdeRx[ij] + HGC_halfrocChannel[ij];
				
				if(HGC_econdeRx[ij]==0 || HGC_econdeRx[ij]==1){
				n_ADC = HGC_adc[ij]-ped_l2[N_channel]-CM1_L2;
				}
				
                                if(HGC_econdeRx[ij]==2 || HGC_econdeRx[ij]==3){
                                n_ADC = HGC_adc[ij]-ped_l2[N_channel]-CM2_L2;
                                }
                                
                                if(HGC_econdeRx[ij]==4 || HGC_econdeRx[ij]==5){
                                n_ADC = HGC_adc[ij]-ped_l2[N_channel]-CM3_L2;
                                }
                                
                                 Ch_ADC_L2[N_channel].push_back(n_ADC);  // Filling vectors of each channel of Layer 2 alongside, for further calculation 
                                
				trigg_time_l2[N_channel].push_back(hgcMetadata_trigTime);
				
			
				
				
				
				if (Highest_ADC_L2<n_ADC){
				Highest_ADC_L2=n_ADC;
				seed2=N_channel;
				}
				
				
				
				if(hgcMetadata_trigTime>=115 && hgcMetadata_trigTime<=120){
				
				h_trigTime_l2[N_channel]->Fill(hgcMetadata_trigTime,n_ADC);
				
				profileX_L2[N_channel]->Fill(hgcMetadata_trigTime,n_ADC);
				ADC_Hist_L2[N_channel]->Fill(n_ADC);
				
				}
				
				
				
			}
			} //HGC Cells loop closed
			
			
			
			if(seed1>-1){
			Highest_ADC_channel_L1.push_back(seed1);   // Vectors are storing highest ADC channel of each Event
			} 
			if(seed2>-1){
			Highest_ADC_channel_L2.push_back(seed2);
			}
			
			Seed_L1->Fill(seed1);
			Seed_L2->Fill(seed2);
				
	
	}     //Event Loop Closed
	
	
	std::ofstream outputFile("Mean_ADC_Ch%i.txt");
        
        outputFile <<"Channel	"<<"L1_Mean_ADC	"<<"L2_Mean_ADC	"<<std::endl;
	
	//Writing Channel and Seed iteration Histogram in root file
	fout->cd();
	for(int jk=0; jk<CHANNEL_MAX;++jk){
	
    	
	//ADC_Hist_L1[jk]->Write();
	//ADC_Hist_L2[jk]->Write();
	
	float mean_L1 = ADC_Hist_L1[jk]->GetMean();
	float mean_L2 = ADC_Hist_L2[jk]->GetMean();
	
	
	//h_trigTime_l1[jk]->Write();
	//h_trigTime_l2[jk]->Write();
	
	//profileX_L1[jk]->Write();
	//profileX_L2[jk]->Write();
	
	
	outputFile <<jk<<"	"<<mean_L1<<"	"<<mean_L2<<std::endl;
	}
	
	outputFile.close();
	
	//Seed_L1->Write();
	//Seed_L2->Write();
	
	//Finding Max Hit Channel in each layer
	int maxbin1 = Seed_L1->GetMaximumBin();
	int max_bin_value1 = Seed_L1->GetBinContent(maxbin1);
	int max_e_cell_L1 = Seed_L1->GetBinLowEdge(maxbin1);
	
	cout<<"no. of channel L1 = "<<max_e_cell_L1<<endl;
	cout<<"no. of events in L1 = "<<max_bin_value1<<endl;
	
	
	int maxbin2 = Seed_L2->GetMaximumBin();
	int max_bin_value2 = Seed_L2->GetBinContent(maxbin2);
	int max_e_cell_L2 = Seed_L2->GetBinLowEdge(maxbin2);
	
	cout<<"no. of channel L2 = "<<max_e_cell_L2<<endl;
	cout<<"no. of events in L2 = "<<max_bin_value2<<endl;
	
	
	
	
	
	
	
	
	
	//Finding Nearby Channels to Seed Channel in Layer 1 and Layer 2
	
	vector<int> HGC_channel_l1;
	vector<double> hgcXN_l1;
	vector<double> hgcYN_l1;


	vector<int> HGC_channel_l2;
	vector<double> hgcXN_l2;
	vector<double> hgcYN_l2;
	
	
	
		HGC_channel_l1.clear();
		hgcXN_l1.clear();
		hgcYN_l1.clear();
		
		
		HGC_channel_l2.clear();
		hgcXN_l2.clear();
		hgcYN_l2.clear();
		
		vector<int>nearby_L1;
		vector<int>nearby_L2;
		
		nearby_L1.clear();
		nearby_L2.clear();
		
		for(int jn=0; jn<1000; jn++){
		if(HGC_channel_l1.size()<198 || HGC_channel_l2.size()<198){
		Tout->GetEntry(jn);    //choosing entry with all 198 non zero values
		
		
		
		
		HGC_channel_l1.clear();
		hgcXN_l1.clear();
		hgcYN_l1.clear();
		
		
		HGC_channel_l2.clear();
		hgcXN_l2.clear();
		hgcYN_l2.clear();
		
		
	for(int ij=0; ij<nHGC;++ij){
			//if (HGC_adc[ij] == 0) continue;
                        if(HGC_x[ij]==-1 && HGC_y[ij]==-1) continue;
			if(HGC_layer[ij]==1){
				N_channel = 39*HGC_econdeRx[ij] + HGC_halfrocChannel[ij];
				//cout<<"Test ch  "<<N_channel<<endl;
				HGC_channel_l1.push_back(N_channel);
				hgcXN_l1.push_back(HGC_x[ij]); 
				hgcYN_l1.push_back(HGC_y[ij]);
				
				}
			
			
			if(HGC_layer[ij]==2){
				N_channel = 39*HGC_econdeRx[ij] + HGC_halfrocChannel[ij];
				HGC_channel_l2.push_back(N_channel);
				hgcXN_l2.push_back(HGC_x[ij]); 
				hgcYN_l2.push_back(HGC_y[ij]);
			        
				}
			
			}
			
			
                nearby_L1.clear();
		nearby_L2.clear();
		            
    		// Finding the position of the max cell in Layer 1
    		auto it1 = find(HGC_channel_l1.begin(), HGC_channel_l1.end(), max_e_cell_L1);

        	int pos_l1 = distance(HGC_channel_l1.begin(), it1);
        	
        	cout<<"Event   "<<jn<<endl;
        	cout<<"Active cells Layer 1  "<<HGC_channel_l1.size()<<endl;
        	
			
	for (int j = 0; j < HGC_channel_l1.size(); j++)
		{
			float s1 = sqrt((hgcXN_l1[j] - hgcXN_l1[pos_l1])*(hgcXN_l1[j] - hgcXN_l1[pos_l1]) + (hgcYN_l1[j] - hgcYN_l1[pos_l1])*(hgcYN_l1[j] - hgcYN_l1[pos_l1]));
			if (s1 == 0) continue;
			if (s1 <= 1.21)
			{
                  		nearby_L1.push_back(HGC_channel_l1[j]);
                  		
			}
			
			}
			
			///
	for (int j = 0; j < HGC_channel_l1.size(); j++)
		{
			float s2 = sqrt((hgcXN_l1[j] - hgcXN_l1[pos_l1])*(hgcXN_l1[j] - hgcXN_l1[pos_l1]) + (hgcYN_l1[j] - hgcYN_l1[pos_l1])*(hgcYN_l1[j] - hgcYN_l1[pos_l1]));
			if (s2 == 0) continue;
			if (s2 > 1.21 && s2 < 2.5){
                  		nearby_L1.push_back(HGC_channel_l1[j]);
                  		
			
			
			}
			
			
			}
			
			
			///
			
	for (int j = 0; j < HGC_channel_l1.size(); j++)
		{
			float s2 = sqrt((hgcXN_l1[j] - hgcXN_l1[pos_l1])*(hgcXN_l1[j] - hgcXN_l1[pos_l1]) + (hgcYN_l1[j] - hgcYN_l1[pos_l1])*(hgcYN_l1[j] - hgcYN_l1[pos_l1]));
			if (s2 == 0) continue;
			if (s2 > 2.5 && s2 < 3.7){
                  		nearby_L1.push_back(HGC_channel_l1[j]);
                  		
			
			
			}
			
			
			}	
			
			
			
			
			
			// Finding the position of the max cell in Layer 2
    		auto it2 = find(HGC_channel_l2.begin(), HGC_channel_l2.end(), max_e_cell_L2);

        	int pos_l2 = distance(HGC_channel_l2.begin(), it2);
        	
        	cout<<"Active cells Layer 2  "<<HGC_channel_l2.size()<<endl;
        	
		
		
		
			
	for (int j = 0; j < HGC_channel_l2.size(); j++)
		{
			float s1 = sqrt((hgcXN_l2[j] - hgcXN_l2[pos_l2])*(hgcXN_l2[j] - hgcXN_l2[pos_l2]) + (hgcYN_l2[j] - hgcYN_l2[pos_l2])*(hgcYN_l2[j] - hgcYN_l2[pos_l2]));
			if (s1 == 0) continue;
			if (s1 <= 1.21){
                  		nearby_L2.push_back(HGC_channel_l2[j]);
                  		
			}
			
			}
			
			
	for (int j = 0; j < HGC_channel_l2.size(); j++)
		{
			float s2 = sqrt((hgcXN_l2[j] - hgcXN_l2[pos_l2])*(hgcXN_l2[j] - hgcXN_l2[pos_l2]) + (hgcYN_l2[j] - hgcYN_l2[pos_l2])*(hgcYN_l2[j] - hgcYN_l2[pos_l2]));
			if (s2 == 0) continue;
			if (s2 > 1.21 && s2 < 2.5){
                  		nearby_L2.push_back(HGC_channel_l2[j]);
                  		
			}
			
			}
		
		
		
		
		for (int j = 0; j < HGC_channel_l2.size(); j++)
		{
			float s2 = sqrt((hgcXN_l2[j] - hgcXN_l2[pos_l2])*(hgcXN_l2[j] - hgcXN_l2[pos_l2]) + (hgcYN_l2[j] - hgcYN_l2[pos_l2])*(hgcYN_l2[j] - hgcYN_l2[pos_l2]));
			if (s2 == 0) continue;
			if (s2 > 2.5 && s2 < 3.7){
                  		nearby_L2.push_back(HGC_channel_l2[j]);
                  		
			}
			
			}
				
		
		}
		
		}///
		
			
			cout<<"Max channel L1 = "<<max_e_cell_L1<<endl;
			int k1=nearby_L1.size();
			for(int k=0; k<nearby_L1.size();k++){
			cout<<"nearby channel l1 = "<<nearby_L1[k]<<endl;
			
			}
			
			cout<<"Max channel L2 = "<<max_e_cell_L2<<endl;
			int k2=nearby_L2.size();
			for(int k=0; k<nearby_L2.size();k++){
			cout<<"nearby channel l2 = "<<nearby_L2[k]<<endl;
			}
			
	
	//Modify the part of the code according to the number of nearby cells for eg. in case of calibration cell nearby Cell no. would increase by 1
	
	//Defining Histograms for nearby Channels
	TH1F *Cells_L1[39], *Cells_L2[39];
	for(int u=0;u<39;u++){
    	sprintf(name,"E%i_L1",u);
    	sprintf(title,"E%i_L1",u);
    	if(u==0){
    	Cells_L1[u]= new TH1F(name, title, 100,-20.0,80);
    	Cells_L1[u]->GetXaxis()->SetTitle("ADC");
	Cells_L1[u]->GetYaxis()->SetTitle("Events");
    	}
    	
    	else{
    	Cells_L1[u]= new TH1F(name, title, 20,-10.0,10);
    	Cells_L1[u]->GetXaxis()->SetTitle("ADC");
	Cells_L1[u]->GetYaxis()->SetTitle("Events");}
	
	}
	
	for(int u=0;u<39;u++){
	
	sprintf(name,"E%i_L2",u);
    	sprintf(title,"E%i_L2",u);
    	if(u==0){
    	Cells_L2[u]= new TH1F(name, title, 100,-20,80);
    	Cells_L2[u]->GetXaxis()->SetTitle("ADC");
	Cells_L2[u]->GetYaxis()->SetTitle("Events");
    	}
    	
    	else{
    	Cells_L2[u]= new TH1F(name, title, 20,-10.0,10);
    	Cells_L2[u]->GetXaxis()->SetTitle("ADC");
	Cells_L2[u]->GetYaxis()->SetTitle("Events");}
        }
	
	
	
	
	
	
	
	//Defining Histogram for Ratio of nearby channel ADC with seed ADC
	
	TH1F *Ratio_E_by_E0_L1[38],*Ratio_E_by_E0_L2[38];
	for(int u=0;u<38;u++){
    	sprintf(name,"E%i_by_E0_L1",u+1);
    	sprintf(title,"E%i_by_E0_L1",u+1);
    	Ratio_E_by_E0_L1[u]= new TH1F(name, title, 200,-1.0,1.0);
    	Ratio_E_by_E0_L1[u]->GetXaxis()->SetTitle("ADC");
	Ratio_E_by_E0_L1[u]->GetYaxis()->SetTitle("Events");
	}
	
	
	
	for(int u=0;u<38;u++){
	sprintf(name,"E%i_by_E0_L2",u+1);
    	sprintf(title,"E%i_by_E0_L2",u+1);
    	Ratio_E_by_E0_L2[u]= new TH1F(name, title, 200,-1.0,1.0);
    	Ratio_E_by_E0_L2[u]->GetXaxis()->SetTitle("ADC");
	Ratio_E_by_E0_L2[u]->GetYaxis()->SetTitle("Events");
	Ratio_E_by_E0_L2[u]->SetMaximum(45000);
	
	}
	
	
	
	
	
	
	
	
	
	// Defining 1D Histograms of Ratio of E0/1st_Ring_including_E0 for Layer 1
        TH1F *ring_hist1_L1;
    	sprintf(name,"E0_by_1st_ring_&_Seed_L1");
    	sprintf(title,"E0_by_1st_ring_&_Seed_L1");
    	ring_hist1_L1= new TH1F(name, title, 200,0,3.0);
    	ring_hist1_L1->GetXaxis()->SetTitle("E0/1st_ring_&_Seed");
	ring_hist1_L1->GetYaxis()->SetTitle("Events");
        
        
        // Defining 1D Histograms of Ratio of E0/1st_and_2nd_Ring_including_E0 for Layer 1
        TH1F *ring_hist2_L1;
    	sprintf(name,"E0_by_1st_2nd_ring_&_Seed_L1");
    	sprintf(title,"E0_by_1st_2nd_ring_&_Seed_L1");
    	ring_hist2_L1= new TH1F(name, title, 200,0,3.0);
    	ring_hist2_L1->GetXaxis()->SetTitle("E0/1st_2nd_ring_&_Seed");
	ring_hist2_L1->GetYaxis()->SetTitle("Events");
	
	
	// Defining 1D Histograms of Ratio of E0/1st_Ring_including_E0 for Layer 2
        TH1F *ring_hist1_L2;
    	sprintf(name,"E0_by_1st_ring_&_Seed_L2");
    	sprintf(title,"E0_by_1st_ring_&_Seed_L2");
    	ring_hist1_L2= new TH1F(name, title, 200,0,3.0);
    	ring_hist1_L2->GetXaxis()->SetTitle("E0/1st_ring_&_Seed");
	ring_hist1_L2->GetYaxis()->SetTitle("Events");
        
        
        // Defining 1D Histograms of Ratio of E0/1st_and_2nd_Ring_including_E0 for Layer 2
        TH1F *ring_hist2_L2;
    	sprintf(name,"E0_by_1st_2nd_ring_&_Seed_L2");
    	sprintf(title,"E0_by_1st_2nd_ring_&_Seed_L2");
    	ring_hist2_L2= new TH1F(name, title, 200,0,3.0);
    	ring_hist2_L2->GetXaxis()->SetTitle("E0/1st_2nd_ring_&_Seed");
	ring_hist2_L2->GetYaxis()->SetTitle("Events");
	
	
	
	
	
	// Defining 2D histograms for E_nearby vs E0 for Layer 1
    	TH2D *hist_E_vsE0_L1[38];
    	char Energy[100];
    	for(int u=0;u<38;u++){
    	sprintf(name,"hist_E%ivsE0_L1",u+1);
    	sprintf(title,"E%i vs E0 Layer 1",u+1);
    	sprintf(Energy,"E%i",u+1);
    	hist_E_vsE0_L1[u]= new TH2D(name, title, 30,-20,100,  120,-10,10);
    	hist_E_vsE0_L1[u]->GetXaxis()->SetTitle("E0");
	hist_E_vsE0_L1[u]->GetYaxis()->SetTitle(Energy);
        } 
        
        
        
        
        // Defining ProfileX plots of nearby cells Layer 1
        TProfile *profx_E_vsE0_L1[38];
        Double_t binEdges[31];
	for (int i = 0; i <= 20; ++i) {
    	binEdges[i] = i * 1.0; // Separation of 1.0
	}
	for (int i = 21; i <= 24; ++i) {
    	binEdges[i] = 20+(i-20.0) * 5.0; // Separation of 5.0
	}
	for (int i = 25; i <= 30; ++i) {
    	binEdges[i] = 40+(i-24) * 10.0; // Separation of 10.0
	}
    	Int_t nBins = sizeof(binEdges) / sizeof(Double_t) - 1;
    	for(int u=0;u<38;u++){
    	sprintf(name,"E%ivsE0_L1",u+1);
    	sprintf(title,"E%i vs E0 Layer 1",u+1);
    	sprintf(Energy,"E%i",u+1);
    	profx_E_vsE0_L1[u]= new TProfile(name, title, 30,binEdges);
    	profx_E_vsE0_L1[u]->GetXaxis()->SetTitle("E0");
	profx_E_vsE0_L1[u]->GetYaxis()->SetTitle(Energy);
        } 
        
        
        
        // Defining 2D histograms for E_nearby vs E0 for Layer 2
    	TH2D *hist_E_vsE0_L2[38];
    	for(int u=0;u<38;u++){
    	sprintf(name,"hist_E%ivsE0_L2",u+1);
    	sprintf(title,"E%i vs E0 Layer 2",u+1);
    	sprintf(Energy,"E%i",u+1);
    	hist_E_vsE0_L2[u]= new TH2D(name, title, 30,-20,100,  120,-10,10);
    	hist_E_vsE0_L2[u]->GetXaxis()->SetTitle("E0");
	hist_E_vsE0_L2[u]->GetYaxis()->SetTitle(Energy);
        } 
        
        
        
        
        
        
        
        //Defining 2D Histogram for CM vs E0 for Layer 2
	TH2D *hist_CM_vsE0_L1[3];
    	for(int u=0;u<3;u++){
    	sprintf(name,"hist_CM%ivsE0_L1",u+1);
    	sprintf(title,"CM%i vs E0 Layer 1",u+1);
    	sprintf(Energy,"CM%i",u+1);
    	hist_CM_vsE0_L1[u]= new TH2D(name, title, 100,-20,80,  120,-10,10);
    	hist_CM_vsE0_L1[u]->GetXaxis()->SetTitle("E0");
	hist_CM_vsE0_L1[u]->GetYaxis()->SetTitle(Energy);
        } 
        
        
        
        
        
        
        TH2D *hist_CM_vsE0_L2[3];
    	for(int u=0;u<3;u++){
    	sprintf(name,"hist_CM%ivsE0_L2",u+1);
    	sprintf(title,"CM%i vs E0 Layer 2",u+1);
    	sprintf(Energy,"CM%i",u+1);
    	hist_CM_vsE0_L2[u]= new TH2D(name, title, 100,-20,80,  120,-10,10);
    	hist_CM_vsE0_L2[u]->GetXaxis()->SetTitle("E0");
	hist_CM_vsE0_L2[u]->GetYaxis()->SetTitle(Energy);
        } 
        
        
        
        
        
        
        TProfile *profx_CM_vsE0_L2[3];
    	for(int u=0;u<3;u++){
    	sprintf(name,"prof_CM%ivsE0_L2",u+1);
    	sprintf(title,"CM%i vs E0 Layer 2",u+1);
    	sprintf(Energy,"CM%i",u+1);
    	profx_CM_vsE0_L2[u]= new TProfile(name, title, 30,binEdges);
    	profx_CM_vsE0_L2[u]->GetXaxis()->SetTitle("E0");
	profx_CM_vsE0_L2[u]->GetYaxis()->SetTitle(Energy);
        } 
        
        
        
        
        
        
        
        
        
	
	
	
	
	// Defining ProfileX plots of nearby cells Layer 2
        TProfile *profx_E_vsE0_L2[38];
    	for(int u=0;u<38;u++){
    	sprintf(name,"E%ivsE0_L2",u+1);
    	sprintf(title,"E%i vs E0 Layer 2",u+1);
    	sprintf(Energy,"E%i",u+1);
    	profx_E_vsE0_L2[u]= new TProfile(name, title, 30,binEdges);
    	profx_E_vsE0_L2[u]->GetXaxis()->SetTitle("E0");
	profx_E_vsE0_L2[u]->GetYaxis()->SetTitle(Energy);
        } 
        
        
        
        
        
        //Defining Efficiency Histogram of Layer 1
        TH1F *x1_L1,*x2_L1,*x1_L2,*x2_L2;
        x1_L1= new TH1F("x1_L1", "x1 Layer 1", 200,-1.0,1.0);
    	x1_L1->GetXaxis()->SetTitle("Ratio");
	x1_L1->GetYaxis()->SetTitle("Events");
        
        
        x2_L1= new TH1F("x2_L1", "x2 Layer 1", 200,-1.0,1.0);
    	x2_L1->GetXaxis()->SetTitle("Ratio");
	x2_L1->GetYaxis()->SetTitle("Events");
	
	
	//Defining Cross Talk Factor Histogram of Layer 2
        x1_L2= new TH1F("x1_L2", "x1 Layer 2", 2000,-1.0,1.0);
    	x1_L2->GetXaxis()->SetTitle("Ratio");
	x1_L2->GetYaxis()->SetTitle("Events");
        
        
        x2_L2= new TH1F("x2_L2", "x2 Layer 2", 2000,-1.0,1.0);
    	x2_L2->GetXaxis()->SetTitle("Ratio");
	x2_L2->GetYaxis()->SetTitle("Events");
        
        
        
        
        
        
        //Defining Efficiency Histogram of Layer 2
        TH1F *n01_L2,*n02_L2;
        n01_L2= new TH1F("n1_L2", "n1 Layer 2", 2000,-1.0,1.0);
    	n01_L2->GetXaxis()->SetTitle("Efficiency");
	n01_L2->GetYaxis()->SetTitle("Events");
        
        n02_L2= new TH1F("n2_L2", "n2 Layer 2", 2000,-1.0,1.0);
    	n02_L2->GetXaxis()->SetTitle("Efficiency");
	n02_L2->GetYaxis()->SetTitle("Events");
       
	 //Defining 2D Histogram of Layer 2 for n2 vs n1
	TH2D *hist_n2_vs_n1_L2;
    	hist_n2_vs_n1_L2= new TH2D("n2_vs_n1_L2", "n2 vs n1 Layer 2", 40,-0.2,0.2,  40,-0.2,0.2);
    	hist_n2_vs_n1_L2->GetXaxis()->SetTitle("n1");
	hist_n2_vs_n1_L2->GetYaxis()->SetTitle("n2");
        
        
        
        
        
        
        //Defining profileX Histogram of Layer 2 for n2 vs n1
        TProfile *profx_n2_vs_n1_L2;
    	profx_n2_vs_n1_L2= new TProfile("profileX_n2_vs_n1_L2", "n2 vs n1 Layer 2 PofileX", 40,-0.2,0.2,-0.2,0.2);
    	profx_n2_vs_n1_L2->GetXaxis()->SetTitle("n1");
	profx_n2_vs_n1_L2->GetYaxis()->SetTitle("n2");
	 
	 
	 
	//Defining Histogrm of x0 Layer 2 
	TH1F *x0_L2;
        x0_L2= new TH1F("x0_L2", "x0_L2", 500,0,5.0);
    	x0_L2->GetXaxis()->SetTitle("x0");
	x0_L2->GetYaxis()->SetTitle("Events"); 
	 
	 
	 
	 //Defining 2D Histogram of Layer 2 for x1 vs E0 & x2 vs E0
	TH2D *hist_x1_vs_E0_L2,*hist_x2_vs_E0_L2;
    	hist_x1_vs_E0_L2= new TH2D("x1_vs_E0_L2", "x1 vs E0 Layer 2", 100,-20,80,  40,-0.2,0.2);
    	hist_x1_vs_E0_L2->GetXaxis()->SetTitle("E0");
	hist_x1_vs_E0_L2->GetYaxis()->SetTitle("x1");
	
	hist_x2_vs_E0_L2= new TH2D("x2_vs_E0_L2", "x2 vs E0 Layer 2", 100,-20,80,  40,-0.2,0.2);
    	hist_x2_vs_E0_L2->GetXaxis()->SetTitle("E0");
	hist_x2_vs_E0_L2->GetYaxis()->SetTitle("x2");
        
        
        
        
        
        
        //Defining profileX Histogram of Layer 2 for x1 vs E0 & x1 vs E0
        TProfile *profx_x1_vs_E0_L2,*profx_x2_vs_E0_L2;
    	profx_x1_vs_E0_L2= new TProfile("profx_x1_vs_E0_L2", "x1 vs E0 Layer 2 PofileX", 100,-20,80,-0.2,0.2);
    	profx_x1_vs_E0_L2->GetXaxis()->SetTitle("E0");
	profx_x1_vs_E0_L2->GetYaxis()->SetTitle("x1"); 
	
	profx_x2_vs_E0_L2= new TProfile("profx_x2_vs_E0_L2", "x2 vs E0 Layer 2 PofileX", 100,-20,80,-0.2,0.2);
    	profx_x2_vs_E0_L2->GetXaxis()->SetTitle("E0");
	profx_x2_vs_E0_L2->GetYaxis()->SetTitle("x2"); 
	 
	
	
	
	
	
	
	
	
	
	
	
	
	
	 //Defining Efficiency Histogram of Layer 2 per cell
        TH1F *x_Ring1_L2[6],*x_Ring2_L2[13];
	
	//Defining Cross Talk Factor Histogram of Layer 2 per cell
	for(int u =0;u<6;u++){
	sprintf(name,"x1_%i_L2",u+1);
    	sprintf(title,"x1_%i_L2",u+1);
        x_Ring1_L2[u]= new TH1F(name, title, 2000,-1.0,1.0);
    	x_Ring1_L2[u]->GetXaxis()->SetTitle("Ratio");
	x_Ring1_L2[u]->GetYaxis()->SetTitle("Events");
	}
        
        for(int u =0;u<13;u++){
        sprintf(name,"x2_%i_L2",u+1);
    	sprintf(title,"x2_%i_L2",u+1);
        x_Ring2_L2[u]= new TH1F(name, title, 2000,-1.0,1.0);
    	x_Ring2_L2[u]->GetXaxis()->SetTitle("Ratio");
	x_Ring2_L2[u]->GetYaxis()->SetTitle("Events");
        }
        
        
        
        
        
        //Defining Efficiency Histogram of Layer 2 per cell
        TH1F *n_Ring1_L2[6],*n_Ring2_L2[13];
        
	for(int u =0;u<6;u++){
	sprintf(name,"n1_%i_L2",u+1);
    	sprintf(title,"n1_%i_L2",u+1);
        n_Ring1_L2[u]= new TH1F(name, title, 2000,-1.0,1.0);
    	n_Ring1_L2[u]->GetXaxis()->SetTitle("Efficiency");
	n_Ring1_L2[u]->GetYaxis()->SetTitle("Events");
	}
        
        for(int u =0;u<13;u++){
        sprintf(name,"n2_%i_L2",u+1);
    	sprintf(title,"n2_%i_L2",u+1);
        n_Ring2_L2[u]= new TH1F(name, title, 2000,-1.0,1.0);
    	n_Ring2_L2[u]->GetXaxis()->SetTitle("Efficiency");
	n_Ring2_L2[u]->GetYaxis()->SetTitle("Events");
        }
	
	
	
	
	
	
	 
	 
	 //Defining Histograms for Each channel for ADC for internal calculation after event selection cut
        
        TH1F *ADC_H_L1[CHANNEL_MAX],*ADC_H_L2[CHANNEL_MAX];
	for(int u=0; u<CHANNEL_MAX; ++u){
		
		int bins=120;
		int xmin=-20;
		int xmax=100;
		
      		sprintf(name,"layer1__ch_%i",u);
      		sprintf(title,"layer1__ch_%i",u);
      		ADC_H_L1[u] = new TH1F(name,title,bins,xmin,xmax);
		ADC_H_L1[u]->GetXaxis()->SetTitle("ADC");
		ADC_H_L1[u]->GetYaxis()->SetTitle("Events");
		
		sprintf(name,"layer2__ch_%i",u);
      		sprintf(title,"layer2__ch_%i",u);
      		ADC_H_L2[u] = new TH1F(name,title,bins,xmin,xmax);
		ADC_H_L2[u]->GetXaxis()->SetTitle("ADC");
		ADC_H_L2[u]->GetYaxis()->SetTitle("Events");
		
        	}
	
	
	
	
	
	
	
	float E_mean_L1[39];
        float sq_sum_L1[39];

        float corr_L1[39];
        for(int i=0;i<39;i++){
        corr_L1[i]=0;
        sq_sum_L1[i]=0;
        }
        
        
        
        
        float E_mean_L2[39];
        float sq_sum_L2[39];

        float corr_L2[39];
        for(int i=0;i<39;i++){
        corr_L2[i]=0;
        sq_sum_L2[i]=0;
        }
      
        
        int t1=0;
        
        
        std::ofstream outputFile1("Mean_ADC_selected_Events.txt");
        
        outputFile1 <<"Channel	"<<"L1_Mean_ADC	"<<"L2_Mean_ADC	"<<std::endl;
        
	for (int i =0; i<nevt; ++i){
	
		Tout->GetEntry(i);
		
		/*if(i>=231320 && i<=231500){
		cout<<"Event  vector "<<i-231320<<"  ADC  "<<Ch_ADC_L2[41][i]<<endl;}*/
	
		if (Highest_ADC_channel_L1[i] == max_e_cell_L1){
		if(hgcMetadata_trigTime>=118 && hgcMetadata_trigTime<=124){
		
		
		
		hist_CM_vsE0_L1[0]->Fill(Ch_ADC_L1[max_e_cell_L1][i],Com1_L1[i]);
		hist_CM_vsE0_L1[1]->Fill(Ch_ADC_L1[max_e_cell_L1][i],Com2_L1[i]);
		hist_CM_vsE0_L1[2]->Fill(Ch_ADC_L1[max_e_cell_L1][i],Com3_L1[i]);
		
		//Filling the Histogram for calculating mean ADC layer 1
		
		for(int jk=0; jk<CHANNEL_MAX;++jk){
		if(Ch_ADC_L1[jk].size()==0)continue;
		ADC_H_L1[jk]->Fill(Ch_ADC_L1[jk][i]);
		
		}
		
		
		
	
		Cells_L1[0]->Fill(Ch_ADC_L1[max_e_cell_L1][i]);
		
		
		
	
		for(int u=0;u<38;u++){
        		Cells_L1[u+1]->Fill(Ch_ADC_L1[nearby_L1[u]][i]);
        		
        		
        		
        
        		
        		Ratio_E_by_E0_L1[u]->Fill(Ch_ADC_L1[nearby_L1[u]][i]/Ch_ADC_L1[max_e_cell_L1][i]);
        		
        		profx_E_vsE0_L1[u]->Fill(Ch_ADC_L1[max_e_cell_L1][i], Ch_ADC_L1[nearby_L1[u]][i]);
        		
        		hist_E_vsE0_L1[u]->Fill(Ch_ADC_L1[max_e_cell_L1][i], Ch_ADC_L1[nearby_L1[u]][i]);
        		
		}
		
		
		//Filling the Histograms for ring 1 and Ring 2 of layer 1
		
		
		float sum_ring1_L1=0;
        	float sum_ring2_L1=0;
        	float n1_L1=0;
        	float n2_L1=0;
        	for(int k=0; k<6;k++){
        		sum_ring1_L1 = sum_ring1_L1+Ch_ADC_L1[nearby_L1[k]][i];
			}
		ring_hist1_L1->Fill(Ch_ADC_L1[max_e_cell_L1][i]/(Ch_ADC_L1[max_e_cell_L1][i]+sum_ring1_L1));
		n1_L1=sum_ring1_L1/(6.0*Ch_ADC_L1[max_e_cell_L1][i]);
		
		
		
		for(int k=6; k<18;k++){
        		sum_ring2_L1 = sum_ring2_L1+Ch_ADC_L1[nearby_L1[k]][i];
			}
		ring_hist2_L1->Fill(Ch_ADC_L1[max_e_cell_L1][i]/(Ch_ADC_L1[max_e_cell_L1][i]+sum_ring1_L1+sum_ring2_L1));
		n2_L1=sum_ring2_L1/(12.0*Ch_ADC_L1[max_e_cell_L1][i]);
		
		
		
		if (Ch_ADC_L1[max_e_cell_L1][i]>=20){
		
		x1_L1->Fill(n1_L1/(1+6.0*n1_L1+12.0*n2_L1));
		x2_L1->Fill(n2_L1/(1+6.0*n1_L1+12.0*n2_L1));
		
		}
		
		}
	  }
	
		if (Highest_ADC_channel_L2[i] == max_e_cell_L2){
		
	
		
		if(trigg_time_l2[max_e_cell_L2][i]>=115 && trigg_time_l2[max_e_cell_L2][i]<=120){
		
		
		profx_CM_vsE0_L2[0]->Fill(Ch_ADC_L2[max_e_cell_L2][i],Com1_L2[i]);
		profx_CM_vsE0_L2[1]->Fill(Ch_ADC_L2[max_e_cell_L2][i],Com2_L2[i]);
		profx_CM_vsE0_L2[2]->Fill(Ch_ADC_L2[max_e_cell_L2][i],Com3_L2[i]);
		
		
		hist_CM_vsE0_L2[0]->Fill(Ch_ADC_L2[max_e_cell_L2][i],Com1_L2[i]);
		hist_CM_vsE0_L2[1]->Fill(Ch_ADC_L2[max_e_cell_L2][i],Com2_L2[i]);
		hist_CM_vsE0_L2[2]->Fill(Ch_ADC_L2[max_e_cell_L2][i],Com3_L2[i]);
		
		
		
	
		//Filling the Histogram for calculating mean ADC layer 2
		for(int jk=0; jk<CHANNEL_MAX;++jk){
		if(Ch_ADC_L2[jk].size()==0)continue;
		ADC_H_L2[jk]->Fill(Ch_ADC_L2[jk][i]);
		
		}
	
		Cells_L2[0]->Fill(Ch_ADC_L2[max_e_cell_L2][i]);
		
		if (Ch_ADC_L2[max_e_cell_L2][i]>=-20){
		for(int u=0;u<38;u++){
       			Cells_L2[u+1]->Fill(Ch_ADC_L2[nearby_L2[u]][i]);
       			
       			
       			
			Ratio_E_by_E0_L2[u]->Fill(Ch_ADC_L2[nearby_L2[u]][i]/Ch_ADC_L2[max_e_cell_L2][i]);
			
			profx_E_vsE0_L2[u]->Fill(Ch_ADC_L2[max_e_cell_L2][i], Ch_ADC_L2[nearby_L2[u]][i]);
			
			hist_E_vsE0_L2[u]->Fill(Ch_ADC_L2[max_e_cell_L2][i], Ch_ADC_L2[nearby_L2[u]][i]);
			
			}
		
		}
		
		//Filling the Histograms for ring 1 and Ring 2 of layer 2
		
		
		float sum_ring1_L2=0;
        	float sum_ring2_L2=0;
        	float n1_L2=0;
        	float n2_L2=0;
        	float x01_L2=0;
        	float x02_L2=0;
        	float n_r1[6];
        	float n_r2[13];
        	float sum_r1=0;
        	float sum_r2=0;
        	//float x0_L2=0;
        	for(int k=0; k<6;k++){
        	//if(k==2)continue;
        		sum_ring1_L2 = sum_ring1_L2+Ch_ADC_L2[nearby_L2[k]][i];
        		
        		//n_Ring1_L2[k]->Fill(Ch_ADC_L2[nearby_L2[k]][i]/Ch_ADC_L2[max_e_cell_L2][i]);
			}
		
		n1_L2=sum_ring1_L2/(6.0*Ch_ADC_L2[max_e_cell_L2][i]);
		
		
		
		for(int k=6; k<19;k++){
        		sum_ring2_L2 = sum_ring2_L2+Ch_ADC_L2[nearby_L2[k]][i];
        		
        		//n_Ring2_L2[k-6]->Fill(Ch_ADC_L2[nearby_L2[k]][i]/Ch_ADC_L2[max_e_cell_L2][i]);
			}
			
		
		n2_L2=sum_ring2_L2/(13.0*Ch_ADC_L2[max_e_cell_L2][i]);
		
		
		//x0_L2 = 1/(1+6.0*n1_L2+13.0*n2_L2);
		x01_L2=n1_L2/(1+6.0*n1_L2+13.0*n2_L2);
		x02_L2=n2_L2/(1+6.0*n1_L2+13.0*n2_L2);
		
		
		hist_x1_vs_E0_L2->Fill(Ch_ADC_L2[max_e_cell_L2][i],x01_L2);
		hist_x2_vs_E0_L2->Fill(Ch_ADC_L2[max_e_cell_L2][i],x02_L2);
		
		profx_x1_vs_E0_L2->Fill(Ch_ADC_L2[max_e_cell_L2][i],x01_L2);
		profx_x2_vs_E0_L2->Fill(Ch_ADC_L2[max_e_cell_L2][i],x02_L2);
		
		if (Ch_ADC_L2[max_e_cell_L2][i]>=20){
		ring_hist1_L2->Fill(Ch_ADC_L2[max_e_cell_L2][i]/(Ch_ADC_L2[max_e_cell_L2][i]+sum_ring1_L2));
		ring_hist2_L2->Fill(Ch_ADC_L2[max_e_cell_L2][i]/(Ch_ADC_L2[max_e_cell_L2][i]+sum_ring2_L2));
		x1_L2->Fill(x01_L2);
		x2_L2->Fill(x02_L2);
		
		//if(x01_L2<0){
		n01_L2->Fill(n1_L2);
		n02_L2->Fill(n2_L2);
		
		x0_L2->Fill(1/(1+6.0*n1_L2+13.0*n2_L2));
		
		
		hist_n2_vs_n1_L2->Fill(n1_L2,n2_L2);
		
		profx_n2_vs_n1_L2->Fill(n1_L2,n2_L2);
		
		
		for(int k=0; k<6;k++){
        	//if(k==2)continue;
        		n_r1[k]=Ch_ADC_L2[nearby_L2[k]][i]/Ch_ADC_L2[max_e_cell_L2][i];
        		n_Ring1_L2[k]->Fill(n_r1[k]);
        		sum_r1=sum_r1+n_r1[k];
			}
			
		
		
		for(int k=6; k<19;k++){
			n_r2[k-6]=Ch_ADC_L2[nearby_L2[k]][i]/Ch_ADC_L2[max_e_cell_L2][i];
        		n_Ring2_L2[k-6]->Fill(n_r2[k-6]);
        		sum_r2=sum_r2+n_r2[k-6];
			}
			
		for (int k =0; k<13;k++){
		if(k<6){
		
		x_Ring1_L2[k]->Fill(n_r1[k]/(1+sum_r1+sum_r2));
		}
	
		
		x_Ring2_L2[k]->Fill(n_r2[k]/(1+sum_r1+sum_r2));
		}
		//}
		
		}
		}
	}
	
	}  //End of selected event loop
	
	
	
	
	for(int kl=0; kl<CHANNEL_MAX;++kl){
		
		float mean01=ADC_H_L1[kl]->GetMean();
		float mean02=ADC_H_L2[kl]->GetMean();
		
		outputFile1 <<kl<<"	"<<mean01<<"	"<<mean02<<std::endl;
		}
	
	
	outputFile1.close();
	
	
	
	
	TF1 *fitFunc = new TF1("fitFunc", "pol1", 20, 100);
        
	//Writing Nearby Channel and Ratio Histogram in Root File
	for(int u=0;u<39;u++){
	
	//Cells_L1[u]->Write();
	E_mean_L1[u]=Cells_L1[u]->GetMean();
	
	
	if (u<38){
	
		
		//Ratio_E_by_E0_L1[u]->Write();
		profx_E_vsE0_L1[u]->Fit(fitFunc, "Rq");
		//gStyle->SetOptFit(1111);
		
		//profx_E_vsE0_L1[u]->Write("PROFILE");
		//gStyle->SetOptFit(16);
		//hist_E_vsE0_L1[u]->Write();
		}
	
	}
	
	float x[20];
	
	gStyle->SetOptStat(111111);
	for(int u=0;u<39;u++){
	
	//Cells_L2[u]->Write();
	E_mean_L2[u]=Cells_L2[u]->GetMean();
	x[u] = u;
    	
    	
	
	if (u<38){
		
		
		//Ratio_E_by_E0_L2[u]->Write();
		profx_E_vsE0_L2[u]->Fit(fitFunc, "Rq");
		//gStyle->SetOptFit(1111);
		
		//profx_E_vsE0_L2[u]->Write("PROFILE");
		//gStyle->SetOptFit(16);
		//hist_E_vsE0_L2[u]->Write();
		}
	}
	
	
	
	TGraph* graph = new TGraph(19, x, E_mean_L2);
	graph->SetTitle("Mean ADC vs Cell No.");
	graph->GetXaxis()->SetTitle("Cell No.");
    	graph->GetYaxis()->SetTitle("Mean ADC");
    	graph->SetMarkerStyle(20); 
	//graph->Write("Mean ADC");
	
	
	

	
	
	
	
	for (int i =0; i<nevt; ++i){
	
		Tout->GetEntry(i);
	
		if (Highest_ADC_channel_L1[i] == max_e_cell_L1){
		if(hgcMetadata_trigTime>=118 && hgcMetadata_trigTime<=124){
		
        	if (Ch_ADC_L1[max_e_cell_L1][i]>=20){
        	corr_L1[0]=corr_L1[0]+(Ch_ADC_L1[max_e_cell_L1][i]-E_mean_L1[0])*(Ch_ADC_L1[max_e_cell_L1][i]-E_mean_L1[0]);
	
		sq_sum_L1[0] = sq_sum_L1[0] + (Ch_ADC_L1[max_e_cell_L1][i]-E_mean_L1[0])*(Ch_ADC_L1[max_e_cell_L1][i]-E_mean_L1[0]);
		
		
		for(int u=0;u<38;u++){
        		
        		
        		corr_L1[u+1]=corr_L1[u+1]+(Ch_ADC_L1[max_e_cell_L1][i]-E_mean_L1[0])*(Ch_ADC_L1[nearby_L1[u]][i]-E_mean_L1[u+1]);
        
        		sq_sum_L1[u+1] = sq_sum_L1[u+1] + (Ch_ADC_L1[nearby_L1[u]][i]-E_mean_L1[u+1])*(Ch_ADC_L1[nearby_L1[u]][i]-E_mean_L1[u+1]);
		
		
		    }
		   }
		  }
		}
		
		
		
		
		
		if (Highest_ADC_channel_L2[i] == max_e_cell_L2){
		
		if(trigg_time_l2[max_e_cell_L2][i]>=115 && trigg_time_l2[max_e_cell_L2][i]<=120){
		
        	if (Ch_ADC_L2[max_e_cell_L2][i]>=20){
        	
        	corr_L2[0]=corr_L2[0]+(Ch_ADC_L2[max_e_cell_L2][i]-E_mean_L2[0])*(Ch_ADC_L2[max_e_cell_L2][i]-E_mean_L2[0]);
        	
        	sq_sum_L2[0] = sq_sum_L2[0] + (Ch_ADC_L2[max_e_cell_L2][i]-E_mean_L2[0])*(Ch_ADC_L2[max_e_cell_L2][i]-E_mean_L2[0]);
		
		
		   
		   
		   
		   
		for(int u=0;u<38;u++){
       			
       			corr_L2[u+1]=corr_L2[u+1]+(Ch_ADC_L2[max_e_cell_L2][i]-E_mean_L2[0])*(Ch_ADC_L2[nearby_L2[u]][i]-E_mean_L2[u+1]);
        
        		sq_sum_L2[u+1] = sq_sum_L2[u+1] + (Ch_ADC_L2[nearby_L2[u]][i]-E_mean_L2[u+1])*(Ch_ADC_L2[nearby_L2[u]][i]-E_mean_L2[u+1]);
       			
       			
       				}
       			}
	 		 }
		   }
		
		
	}
		
		
		
	std::ofstream outputFile2("Correlation_Coefficients.txt");
        
		
		
		
	
	outputFile2 <<"Layer	"<<"Cell_name	"<<"correlation_with_E0	"<<"Global_ID	"<<endl;
	for(int u=0;u<39;u++){
        corr_L1[u]=corr_L1[u]/pow(sq_sum_L1[0]*sq_sum_L1[u],0.5);
        if(u==0){
        outputFile2 <<"1"<<"	"<<u<<"	"<<corr_L1[u]<<"	"<<max_e_cell_L1<<endl;
        }
        else{
        outputFile2 <<"1"<<"	"<<u<<"	"<<corr_L1[u]<<"	"<<nearby_L1[u-1]<<endl;
        }
        }
        
        
        
        
        
        
        
        
        for(int u=0;u<39;u++){
        corr_L2[u]=corr_L2[u]/pow(sq_sum_L2[0]*sq_sum_L2[u],0.5);
      	if(u==0){
        outputFile2 <<"2"<<"	"<<u<<"	"<<corr_L2[u]<<"	"<<max_e_cell_L2<<endl;
        }
        else{
        outputFile2 <<"2"<<"	"<<u<<"	"<<corr_L2[u]<<"	"<<nearby_L2[u-1]<<endl;
        }
        }
        
        outputFile2.close();
        
        
	
	
	//Writing 1st Ring and 2nd Ring including seed histogram in Root file
	
	//ring_hist1_L1->Write();
	//ring_hist2_L1->Write();
	
	//ring_hist1_L2->Write();
	//ring_hist2_L2->Write();
	
	
	//hist_CM_vsE0_L1[0]->Write();
	//hist_CM_vsE0_L1[1]->Write();
	//hist_CM_vsE0_L1[2]->Write();
	
	//hist_CM_vsE0_L2[0]->Write();
	//hist_CM_vsE0_L2[1]->Write();
	//hist_CM_vsE0_L2[2]->Write();
	
	
	//profx_CM_vsE0_L2[0]->Write("PROFILE");
	//profx_CM_vsE0_L2[1]->Write("PROFILE");
	//profx_CM_vsE0_L2[2]->Write("PROFILE");
	
	
	for (int k =0; k<13;k++){
	if(k<6){
	n_Ring1_L2[k]->Write();
	x_Ring1_L2[k]->Write();
	}
	
	n_Ring2_L2[k]->Write();
	x_Ring2_L2[k]->Write();
	}
	
	
	
	
	//Writing Histograms of efficiency
	
	//x1_L1->Write();
	//x2_L1->Write();
	
	//x1_L2->Write();
	//x2_L2->Write();
	
	//n01_L2->Write();
	//n02_L2->Write();
	
	
	
	
	//hist_n2_vs_n1_L2->Write();
	
	//profx_n2_vs_n1_L2->Write();
	
	
	//x0_L2->Write();
	
	
	//hist_x1_vs_E0_L2->Write();
	//hist_x2_vs_E0_L2->Write();
	
	//profx_x1_vs_E0_L2->Write();
	//profx_x2_vs_E0_L2->Write();
	//Creating E1/E0, E2/E0, E3/E0, E4/E0, E5/E0, E6/E0 on same canvas for Layer 1
	
	TCanvas *canvas_L1 = new TCanvas("canvas_L1", "Overlayed Histograms", 800, 600);
	
	Ratio_E_by_E0_L1[0]->SetLineColor(kRed);
    	Ratio_E_by_E0_L1[0]->Draw("same");
    	Ratio_E_by_E0_L1[1]->SetLineColor(kBlue);
    	Ratio_E_by_E0_L1[1]->Draw("same");
    	Ratio_E_by_E0_L1[2]->SetLineColor(kGreen);
    	Ratio_E_by_E0_L1[2]->Draw("same");
    	Ratio_E_by_E0_L1[3]->SetLineColor(kOrange);
    	Ratio_E_by_E0_L1[3]->Draw("same");
    	Ratio_E_by_E0_L1[4]->SetLineColor(kViolet);
    	Ratio_E_by_E0_L1[4]->Draw("same");
    	Ratio_E_by_E0_L1[5]->SetLineColor(kBlack);
    	Ratio_E_by_E0_L1[5]->Draw("same");
    	
    	TLegend *legend1 = new TLegend(0.78, 0.2, 0.9, 0.65);
	legend1->AddEntry(Ratio_E_by_E0_L1[0], "E1/E0", "l");
    	legend1->AddEntry(Ratio_E_by_E0_L1[1], "E2/E0", "l");
    	legend1->AddEntry(Ratio_E_by_E0_L1[2], "E3/E0", "l");
    	legend1->AddEntry(Ratio_E_by_E0_L1[3], "E4/E0", "l");
    	legend1->AddEntry(Ratio_E_by_E0_L1[4], "E5/E0", "l");
    	legend1->AddEntry(Ratio_E_by_E0_L1[5], "E6/E0", "l");
    	
    	legend1->Draw();
	
	//canvas_L1->Write();
	delete canvas_L1;
	
	
	
	
	
	
	//Creating E1/E0, E2/E0, E3/E0, E4/E0, E5/E0, E6/E0 on same canvas for Layer 2
	
	TCanvas *canvas_L2 = new TCanvas("canvas_L2", "Overlayed Histograms", 800, 600);
	
	Ratio_E_by_E0_L2[0]->SetLineColor(kRed);
    	Ratio_E_by_E0_L2[0]->Draw("same");
    	Ratio_E_by_E0_L2[1]->SetLineColor(kBlue);
    	Ratio_E_by_E0_L2[1]->Draw("same");
    	Ratio_E_by_E0_L2[2]->SetLineColor(kGreen);
    	Ratio_E_by_E0_L2[2]->Draw("same");
    	Ratio_E_by_E0_L2[3]->SetLineColor(kOrange);
    	Ratio_E_by_E0_L2[3]->Draw("same");
    	Ratio_E_by_E0_L2[4]->SetLineColor(kViolet);
    	Ratio_E_by_E0_L2[4]->Draw("same");
    	Ratio_E_by_E0_L2[5]->SetLineColor(kBlack);
    	Ratio_E_by_E0_L2[5]->Draw("same");
    	
    	TLegend *legend2 = new TLegend(0.78, 0.2, 0.9, 0.65);
	legend2->AddEntry(Ratio_E_by_E0_L2[0], "E1/E0", "l");
    	legend2->AddEntry(Ratio_E_by_E0_L2[1], "E2/E0", "l");
    	legend2->AddEntry(Ratio_E_by_E0_L2[2], "E3/E0", "l");
    	legend2->AddEntry(Ratio_E_by_E0_L2[3], "E4/E0", "l");
    	legend2->AddEntry(Ratio_E_by_E0_L2[4], "E5/E0", "l");
    	legend2->AddEntry(Ratio_E_by_E0_L2[5], "E6/E0", "l");
    	
    	legend2->Draw();
	
	//canvas_L2->Write();
	delete canvas_L2;
	
	fout->Close();
    	 

}
