void ordering(TString name, TString file){
	TFile* F = TFile::Open(name);
	TTree* T = (TTree*) F->Get("T");

	vector<int>* Cells = 0;
	vector<int>* Channel = 0;
	vector<double>* CellTime = 0;

	vector<int> CellsAll;
	vector<int> ChannelAll;
	vector<double> CellTimeAll;

	T->SetBranchAddress("Cells", &Cells);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("CellTime", &CellTime);

	int eventsPerCycle = 100;
	int Ni = (T->GetEntries() / eventsPerCycle);
	if(T->GetEntries()%eventsPerCycle == 0) Ni -= 1;
	int nHalf = 0;
	vector<int> tempCells;
	vector<int> tempChannel;
	vector<double> tempCellTime;

	TFile* N = TFile::Open(file + ".root", "RECREATE");
	TTree* Tnew = new TTree("Tnew", "signals");
	N->cd();
	int sCells = 0, sChannel = 0, sFiber = 0;
	double sCellTime = 0;

	Tnew->Branch("Cells", &sCells);
	Tnew->Branch("Channel", &sChannel);
	Tnew->Branch("Fiber", &sFiber);
	Tnew->Branch("CellTime", &sCellTime);

	for(int i = 0; i < Ni; i++){
		for(int k = i * eventsPerCycle; k < eventsPerCycle * (i + 1); k++){
			T->GetEntry(k);
			int n = Cells->size();
			if(n > 0){
				for(int j = 0; j < n; j++){
					CellsAll.push_back(Cells->at(j));
					ChannelAll.push_back(Channel->at(j));
					CellTimeAll.push_back(CellTime->at(j));
				}
			}
		}
		Cells->clear();
		Channel->clear();
		CellTime->clear();

		int n = CellsAll.size();
		double CellTimeL[n];
		int I[n];

		for(int i = 0; i < n; i++){
			CellTimeL[i] = CellTimeAll.at(i);
		}

		CellTimeAll.clear();
		TMath::Sort(n, CellTimeL, I, false);
		
		int tempnHalf = n/2;
		if(i == 0){
			for(int j = 0; j < tempnHalf; j++){
				sCells = CellsAll.at(I[j]);
				sChannel = ChannelAll.at(I[j]);
				sCellTime = CellTimeAll[I[j]];
				sFiber = sChannel/2;
				sChannel = sChannel%2;
				Tnew->Fill();
			}

			for(int j = tempnHalf; j < n; j++){
				tempCells.push_back(CellsAll.at(I[j]));
				tempChannel.push_back(ChannelAll.at(I[j]));
				tempCellTime.push_back(CellTimeL[I[j]]);
			}
			nHalf = n - tempnHalf;
		}

		else{
			for(int j = 0; j < tempnHalf; j++){
				tempCells.push_back(CellsAll.at(I[j]));
				tempChannel.push_back(ChannelAll.at(I[j]));
				tempCellTime.push_back(CellTimeL[I[j]]);
			}

			int ntemp = tempCellTime.size();
			double tempCellTimeL[ntemp];
			int tempI[ntemp];
			for(int j = 0; j < ntemp; j++){
				tempCellTimeL[j] = tempCellTime.at(j);
			}
			tempCellTime.clear();

			TMath::Sort(ntemp, tempCellTimeL, tempI, false);

			for(int j = 0; j < ntemp; j++){
				sCells = tempCells.at(tempI[j]);
				sChannel = tempChannel.at(tempI[j]);
				sCellTime = tempCellTimeL[tempI[j]];
				sFiber = sChannel/2;
				sChannel = sChannel%2;
				Tnew->Fill();
			}
			tempCells.clear();
			tempChannel.clear();

			nHalf = n - tempnHalf;
			for(int j = tempnHalf; j < n; j++){
				tempCells.push_back(CellsAll.at(j));
				tempChannel.push_back(ChannelAll.at(j));
				tempCellTime.push_back(CellTimeL[j]);
			}
		}
		CellsAll.clear();
		ChannelAll.clear();
	}

	for(int i = Ni*eventsPerCycle; i < T->GetEntries(); ++i){
		T->GetEntry(i);
		int n = Cells->size();
		if(n > 0){
			for(int j = 0; j < n; j++){
				CellsAll.push_back(Cells->at(j));
				ChannelAll.push_back(Channel->at(j));
				CellTimeAll.push_back(CellTime->at(j));
			}
		}
	}
	Cells->clear();
	Channel->clear();
	CellTime->clear();

	int n = CellsAll.size();
	double CellTimeL[n];
	int I[n];

	for(int i = 0; i < n; i++){
		CellTimeL[i] = CellTimeAll.at(i);
	}

	CellTimeAll.clear();
	TMath::Sort(n, CellTimeL, I, false);
	
	for(int j = 0; j < n; j++){
		tempCells.push_back(CellsAll.at(I[j]));
		tempChannel.push_back(ChannelAll.at(I[j]));
		tempCellTime.push_back(CellTimeL[I[j]]);
	}

	CellsAll.clear();
	ChannelAll.clear();

	int ntemp = tempCellTime.size();
	double tempCellTimeL[ntemp];
	int tempI[ntemp];
	for(int j = 0; j < ntemp; j++){
		tempCellTimeL[j] = tempCellTime.at(j);
	}
	tempCellTime.clear();

	TMath::Sort(ntemp, tempCellTimeL, tempI, false);

	for(int j = 0; j < ntemp; j++){
		sCells = tempCells.at(tempI[j]);
		sChannel = tempChannel.at(tempI[j]);
		sCellTime = tempCellTimeL[tempI[j]];
		sFiber = sChannel/2;
		sChannel = sChannel%2;
		Tnew->Fill();
	}

	tempCells.clear();
	tempChannel.clear();


	N->cd();
	Tnew->Write("T", TObject::kOverwrite);
	N->Close();
	F->Close();
}		

void preprocessing(TString file){
	TFile* F = TFile::Open(file + ".root");
	TTree* T = (TTree*) F->Get("T");

	double activationTime[84][2668];
	for(int i = 0; i < 84; i++){
		for(int j = 0; j < 2668; j++){
			activationTime[i][j] = -30; // ns
		}
	}

	vector<TTree*> Tnew;
	int Cell = 0, Channel = 0, Fiber = 0;
	double Time = 0;

	T->SetBranchAddress("Cells", &Cell);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Fiber", &Fiber);
	T->SetBranchAddress("CellTime", &Time);

	TFile* N = TFile::Open(file + ".root", "RECREATE");

	for(int i = 0; i < 42; i++){
		Tnew.push_back(new TTree(TString::Format("T%02d", i), TString::Format("Ch%02d", i)));
		Tnew.at(i)->Branch("Cell", &Cell);
		Tnew.at(i)->Branch("Channel", &Channel);
		Tnew.at(i)->Branch("Time", &Time);
	}


	int Ntot = T->GetEntries();
	for(int i = 0; i < Ntot; i++){
		T->GetEntry(i);
		if(Time - activationTime[Fiber * 2 + Channel][Cell] > 20){
			activationTime[Fiber * 2 + Channel][Cell] = Time;
			Tnew.at(Fiber)->Fill();
		}
	}

	F->Close();
	N->cd();
	for(int i = 0; i < 42; i++){
		Tnew.at(i)->Write(TString::Format("Fiber%02d", i), TObject::kOverwrite);
	}
	N->Close();
}

double expo(double t, double offset){
	if(t > offset) return TMath::Exp(-(t - offset) / 15.);
	else return 0;
}

void createSignals(double threashold, TString file){
	TFile* F = TFile::Open(file + ".root", "UPDATE");
	TTree* T;
	TTree* Twaves = new TTree("waves", "signals");
	
	double deltaT = 200; // ns
	double pitch  = 0.1; // ns

	vector<double> signal[2];
	vector<double> signalT[2];
	signal[0].resize(int(deltaT / pitch) + 1);
	signal[1].resize(int(deltaT / pitch) + 1);
	for(int i = 0; i < int(deltaT / pitch) + 1; i++){
		for(int j = 0; j < 2; j++) signalT[j].push_back(-200 + pitch * i);
	}
	vector<double> activationTime[2];

	int Channel = 0, Fiber = 0;
	double signalTime, Amplitude = 0;

	TGraph G;
	
	Twaves->Branch("Signal", &G);
	Twaves->Branch("Channel", &Channel);
	Twaves->Branch("Fiber", &Fiber);
	Twaves->Branch("Amplitude", &Amplitude);
	Twaves->Branch("Time", &signalTime);

	int thCheck[2] = {0}, timeCheck[2] = {- 30};

	double globalTime[2] = {0}; // ns, time of the Last Point
	double Time = 0;

//	TF1* f = new TF1("f", "[0]*(x<[1])*TMath::Exp(-(x-[1])*(x-[1])/2/[2]) + [0]*(x>[1])*TMath::Exp(-(x-[1])/[3])");
	TF1* f = new TF1("f", "[0]*(x<[1])*TMath::Exp(-(x-[1])*(x-[1])/2/[2]) + [0]*(x>[1])*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2] + [3]*(x-[1])))");
	f->SetParLimits(0, 10, 1000);
	f->SetParLimits(2, 1, 10);
	f->SetParLimits(3, 10, 20);
	
	for(int i = 0; i < 42; i++){
		thCheck[0] = 0;
		thCheck[1] = 0;
		Fiber = i;
		T = (TTree*) F->Get(TString::Format("Fiber%02d", i));
		T->SetBranchAddress("Time", &Time);
		T->SetBranchAddress("Channel", &Channel);
		for(int j = 0; j < T->GetEntries(); j++){
			T->GetEntry(j);
			while(globalTime[Channel] < Time){
				signal[Channel].erase(signal[Channel].begin());
				signalT[Channel].erase(signalT[Channel].begin());
				double amplitude = 0;
				int n = activationTime[Channel].size();
				int actual = 0;
				for(int k = 0; k < n; k++){
					amplitude += expo(globalTime[Channel] , activationTime[Channel].at(actual));
					if(globalTime[Channel] - activationTime[Channel].at(actual) > 30){ 
						activationTime[Channel].erase(activationTime[Channel].begin() + actual);
						actual --;
					}
					actual ++;
				}
				signal[Channel].push_back(amplitude);
				globalTime[Channel] += pitch;
				signalT[Channel].push_back(globalTime[Channel] );
				if(signal[Channel].at(int(deltaT/pitch/2)) > threashold && thCheck[Channel] == 0){
					thCheck[Channel] = 1;
					timeCheck[Channel]= signalT[Channel].at(int(deltaT/pitch/2));
					signalTime = globalTime[Channel] - deltaT/2;
					G = TGraph(signal[Channel].size());
					for(int l = 0; l < signal[Channel].size(); l++){
						G.SetPoint(l, signalT[Channel].at(l), signal[Channel].at(l));
					}
					f->SetParameter(0, 10);
					f->SetParameter(1, 5 + signalT[Channel].at(int(deltaT/pitch/2)));
					f->SetParLimits(1, signalT[Channel].at(int(deltaT/pitch/2)),20 + signalT[Channel].at(int(deltaT/pitch/2)));
					f->SetParameter(2, 5);
					f->SetParameter(3, 15);
					f->SetRange(signalT[Channel].at(int(deltaT/pitch/2)), signalT[Channel].at(int(deltaT/pitch/2)) + 16.5);
					for(int g = 0; g < 10; g++) G.Fit("f", "0", "", signalT[Channel].at(int(deltaT/pitch/2)), signalT[Channel].at(int(deltaT/pitch/2)) + 16.5);
					Amplitude = G.GetFunction("f")->GetParameter(0);
					Twaves->Fill();
				}
				else if(signal[Channel].at(int(deltaT/pitch/2)) < threashold && thCheck[Channel] == 1 && signalT[Channel].at(int(deltaT/pitch/2)) - timeCheck[Channel]> 20){
					thCheck[Channel] = 0;
				}
			}
			activationTime[Channel].push_back(Time);
		}

		for(int j = 0; j < int(deltaT/pitch/2); j++){
			signal[Channel].erase(signal[Channel].begin());
			signalT[Channel].erase(signalT[Channel].begin());
			double amplitude = 0;
			int n = activationTime[Channel].size();
			int actual = 0;
			for(int k = 0; k < n; k++){
				amplitude += expo(globalTime[Channel] , activationTime[Channel].at(actual));
				if(globalTime[Channel] - activationTime[Channel].at(actual) > 30){ 
					activationTime[Channel].erase(activationTime[Channel].begin() + actual);
					actual --;
				}
				actual ++;
			}
			signal[Channel].push_back(amplitude);
			globalTime[Channel] += pitch;
			signalT[Channel].push_back(globalTime[Channel] );
			if(signal[Channel].at(int(deltaT/pitch/2)) > threashold && thCheck[Channel] == 0){
				thCheck[Channel] = 1;
				timeCheck[Channel]= signalT[Channel].at(int(deltaT/pitch/2));
				signalTime = globalTime[Channel] - deltaT/2;
				G = TGraph(signal[Channel].size());
				for(int l = 0; l < signal[Channel].size(); l++){
					G.SetPoint(l, signalT[Channel].at(l), signal[Channel].at(l));
				}
				f->SetParameter(0, 10);
				f->SetParameter(1, 5 + signalT[Channel].at(int(deltaT/pitch/2)));
				f->SetParLimits(1, signalT[Channel].at(int(deltaT/pitch/2)),20 + signalT[Channel].at(int(deltaT/pitch/2)));
				f->SetParameter(2, 5);
				f->SetParameter(3, 15);
				f->SetRange(signalT[Channel].at(int(deltaT/pitch/2)), signalT[Channel].at(int(deltaT/pitch/2)) + 16.5);
				for(int g = 0; g < 10; g++) G.Fit("f", "0", "", signalT[Channel].at(int(deltaT/pitch/2)), signalT[Channel].at(int(deltaT/pitch/2)) + 16.5);
				Amplitude = G.GetFunction("f")->GetParameter(0);
				Twaves->Fill();
			}
			else if(signal[Channel].at(int(deltaT/pitch/2)) < threashold && thCheck[Channel] == 1 && signalT[Channel].at(int(deltaT/pitch/2)) - timeCheck[Channel]> 20){
				thCheck[Channel] = 0;
			}
		}
		signal[0].clear();
		signal[0].resize(int(deltaT / pitch) + 1);
		globalTime[0] = 0;
		timeCheck[0]= -30;
		signal[1].clear();
		signal[1].resize(int(deltaT / pitch) + 1);
		globalTime[1] = 0;
		timeCheck[1]= -30;
	}
	
	Twaves->Write("waves", TObject::kOverwrite);
	F->Close();
}

void reordering(TString file){
	TFile* F = TFile::Open(file + ".root", "UPDATE");
	TTree* T = (TTree*) F->Get("waves");

	int Channel = 0;
	int Fiber = 0;
	double Amplitude = 0;
	double Time = 0;

	vector<int> ChannelAll;
	vector<int> FiberAll;
	vector<double> AmplitudeAll;
	vector<double> TimeAll;

	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Fiber", &Fiber);
	T->SetBranchAddress("Amplitude", &Amplitude);
	T->SetBranchAddress("Time", &Time);

	int eventsPerCycle = 1e6;
	int Ni = (T->GetEntries() / eventsPerCycle);
	if(T->GetEntries()%eventsPerCycle == 0) Ni -= 1;
	int nHalf = 0;

	vector<int> tempChannel;
	vector<int> tempFiber;
	vector<double> tempAmplitude;
	vector<double> tempTime;

	TTree* Tnew = new TTree("Tnew", "reordered");

	int sChannel = 0, sFiber = 0;
	double sAmplitude = 0, sTime = 0;

	Tnew->Branch("Channel", &sChannel);
	Tnew->Branch("Fiber", &sFiber);
	Tnew->Branch("Amplitude", &sAmplitude);
	Tnew->Branch("Time", &sTime);

	for(int i = 0; i < Ni; i++){
		for(int k = i * eventsPerCycle; k < eventsPerCycle * (i + 1); k++){
			T->GetEntry(k);
			ChannelAll.push_back(Channel);
			FiberAll.push_back(Fiber);
			AmplitudeAll.push_back(Amplitude);
			TimeAll.push_back(Time);
		}

		int n = ChannelAll.size();
		double TimeL[n];
		int I[n];

		for(int i = 0; i < n; i++){
			TimeL[i] = TimeAll.at(i);
		}

		TimeAll.clear();
		TMath::Sort(n, TimeL, I, false);
		
		int tempnHalf = n/2;
		if(i == 0){
			for(int j = 0; j < tempnHalf; j++){
				sChannel = ChannelAll.at(I[j]);
				sFiber = FiberAll.at(I[j]);
				sAmplitude = AmplitudeAll[I[j]];
				sTime = TimeAll[I[j]];
				Tnew->Fill();
			}

			for(int j = tempnHalf; j < n; j++){
				tempChannel.push_back(ChannelAll.at(I[j]));
				tempFiber.push_back(FiberAll.at(I[j]));
				tempAmplitude.push_back(AmplitudeAll.at(I[j]));
				tempTime.push_back(TimeL[I[j]]);
			}
			nHalf = n - tempnHalf;
		}

		else{
			for(int j = 0; j < tempnHalf; j++){
				tempChannel.push_back(ChannelAll.at(I[j]));
				tempFiber.push_back(FiberAll.at(I[j]));
				tempAmplitude.push_back(AmplitudeAll.at(I[j]));
				tempTime.push_back(TimeL[I[j]]);
			}

			int ntemp = tempTime.size();
			double tempTimeL[ntemp];
			int tempI[ntemp];
			for(int j = 0; j < ntemp; j++){
				tempTimeL[j] = tempTime.at(j);
			}
			tempTime.clear();

			TMath::Sort(ntemp, tempTimeL, tempI, false);

			for(int j = 0; j < ntemp; j++){
				sChannel = tempChannel.at(tempI[j]);
				sFiber = tempFiber.at(tempI[j]);
				sAmplitude = tempAmplitude.at(tempI[j]);
				sTime = tempTimeL[tempI[j]];
				Tnew->Fill();
			}
			tempChannel.clear();
			tempFiber.clear();
			tempAmplitude.clear();

			nHalf = n - tempnHalf;
			for(int j = tempnHalf; j < n; j++){
				tempChannel.push_back(ChannelAll.at(j));
				tempFiber.push_back(FiberAll.at(j));
				tempAmplitude.push_back(AmplitudeAll.at(j));
				tempTime.push_back(TimeL[j]);
			}
		}
		ChannelAll.clear();
		FiberAll.clear();
		AmplitudeAll.clear();
	}

	for(int i = Ni*eventsPerCycle; i < T->GetEntries(); ++i){
		T->GetEntry(i);
		ChannelAll.push_back(Channel);
		FiberAll.push_back(Fiber);
		AmplitudeAll.push_back(Amplitude);
		TimeAll.push_back(Time);
	}

	int n = ChannelAll.size();
	double TimeL[n];
	int I[n];

	for(int i = 0; i < n; i++){
		TimeL[i] = TimeAll.at(i);
	}

	TimeAll.clear();
	TMath::Sort(n, TimeL, I, false);
	
	for(int j = 0; j < n; j++){
		tempChannel.push_back(ChannelAll.at(I[j]));
		tempFiber.push_back(FiberAll.at(I[j]));
		tempAmplitude.push_back(AmplitudeAll.at(I[j]));
		tempTime.push_back(TimeL[I[j]]);
	}

	ChannelAll.clear();
	FiberAll.clear();
	AmplitudeAll.clear();

	int ntemp = tempTime.size();
	double tempTimeL[ntemp];
	int tempI[ntemp];
	for(int j = 0; j < ntemp; j++){
		tempTimeL[j] = tempTime.at(j);
	}
	tempTime.clear();

	TMath::Sort(ntemp, tempTimeL, tempI, false);

	for(int j = 0; j < ntemp; j++){
		sChannel = tempChannel.at(tempI[j]);
		sFiber = tempFiber.at(tempI[j]);
		sAmplitude = tempAmplitude.at(tempI[j]);
		sTime = tempTimeL[tempI[j]];
		Tnew->Fill();
	}

	tempChannel.clear();
	tempFiber.clear();
	tempAmplitude.clear();

	F->cd();
	Tnew->Write("rwaves", TObject::kOverwrite);
	F->Close();
}		

void processing(double threashold, TString file){
	TFile* F = TFile::Open(file + ".root");
	TTree* T = (TTree*) F->Get("rwaves");

	double Amplitude = 0, Time = 0;
	double activationTime[42] = {-30}, preTime[42] = {-30}, deltaTime = 0;
	int Channel = 0, Fiber = 0, NCounts[42] = {0}, preChannel[42] = {0};

	T->SetBranchAddress("Amplitude", &Amplitude);
	T->SetBranchAddress("Time", &Time);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Fiber", &Fiber);

	int n_entries = T->GetEntries();
	for(int i = 0; i < n_entries; i++){
		T->GetEntry(i);
		if(preTime[Fiber] < 0 && Amplitude > threashold){
			preTime[Fiber] = Time;
			preChannel[Fiber] = Channel;
		}
		else if(Time - preTime[Fiber] < 20 && preChannel[Fiber] != Channel && Amplitude > threashold){
			activationTime[Fiber] = Time;
			NCounts[Fiber]++;
		}
		else if(Time - activationTime[Fiber] > 20 && Amplitude > threashold){
			preTime[Fiber] = Time;
			preChannel[Fiber] = Channel;
		}

		if(deltaTime < Time) deltaTime = Time;
	}

	ofstream myfile;
	myfile.open(file + ".txt");
	for(int i = 0; i < 42; i++){
		myfile << NCounts[i] << endl;
	}
	myfile << 0.22 << endl;
	myfile << deltaTime << endl;
	myfile.close();
	F->Close();
}


void signals(){
	ordering("data.root", "out");
	preprocessing("out");
	createSignals(0.5, "out");
	reordering("out");
	processing(0.5, "out");
}

void signals(double threashold){
	ordering("data.root", "out");
	preprocessing("out");
	createSignals(threashold, "out");
	reordering("out");
	processing(threashold, "out");
}

void signals(TString name){
	ordering(name, "out");
	preprocessing("out");
	createSignals(0.5, "out");
	reordering("out");
	processing(0.5, "out");
}


void signals(TString name, TString out){
	ordering(name, out);
	preprocessing(out);
	createSignals(0.5, out);
	reordering(out);
	processing(0.5, out);
}

void signals(TString name, double threashold, TString out){
	ordering(name, out);
	preprocessing(out);
	createSignals(threashold, out);
	reordering(out);
	processing(threashold, out);
}


