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

	for(int i = 0; i < 84; i++){
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
			Tnew.at(Fiber * 2 + Channel)->Fill();
		}
	}

	F->Close();
	N->cd();
	for(int i = 0; i < 84; i++){
		Tnew.at(i)->Write(TString::Format("Fiber%02dCh%d", i/2, i%2), TObject::kOverwrite);
	}
	N->Close();
}

double Gp[5] = {8.88913032e-01,  6.46867235e+01,  4.16687779e-01, 212.49027149107357, 1.5};

double C(double x,double a,double b){
	return b * sqrt(TMath::Pi()/2) * exp(b*b/2/a/a - x/a)*(TMath::Erf((a*x - b*b)/(a*b*sqrt(2))) + 1);
}

double Signal(double x){
	return Gp[3]*(C(x - Gp[4], Gp[0], Gp[2]) - C(x - Gp[4], Gp[0]*Gp[1]/(Gp[0]+Gp[1]), Gp[2]));
}

void createSignals(double threashold, TString file){
	TFile* F = TFile::Open(file + ".root", "UPDATE");
	TTree* T;
	TTree* Twaves[2];
	Twaves[0] = new TTree("waves0", "signals");
	Twaves[1] = new TTree("waves1", "signals");
	
	double deltaT = 200; // ns
	double pitch  = 0.05; // ns

	vector<double> signal;
	vector<double> signalT;
	signal.resize(int(deltaT / pitch) + 1);
	for(int i = 0; i < int(deltaT / pitch) + 1; i++) signalT.push_back(-200 + pitch * i);
	vector<double> activationTime;

	int Channel = 0, Fiber = 0;
	double signalTime, Amplitude = 0;
	
	Twaves[0]->Branch("Channel", &Channel);
	Twaves[0]->Branch("Fiber", &Fiber);
	Twaves[0]->Branch("Amplitude", &Amplitude);
	Twaves[0]->Branch("Time", &signalTime);

	Twaves[1]->Branch("Channel", &Channel);
	Twaves[1]->Branch("Fiber", &Fiber);
	Twaves[1]->Branch("Amplitude", &Amplitude);
	Twaves[1]->Branch("Time", &signalTime);

	int thCheck = 0, timeCheck = - 30;

	double globalTime = 0; // ns, time of the Last Point
	double Time = 0;

	for(int i = 0; i < 84; i++){
		thCheck = 0;
		Fiber = i/2;
		T = (TTree*) F->Get(TString::Format("Fiber%02dCh%d", i/2, i%2));
		T->SetBranchAddress("Time", &Time);
		T->SetBranchAddress("Channel", &Channel);
		int n_entries = T->GetEntries();
		int last = -1;
		for(int j = 0; j < n_entries; j++){
			if(int(100. * j / n_entries)%10 == 0 && int(100. * j / n_entries)/10 > last){
				TString output = TString::Format("Ch. %02d, Fiber %02d, advancement: [", i, i / 2);
				last = int(100. * j / n_entries)/10;
				for(int i = 0; i < 10; i++){
					if (i < last) output = output + "=";
					else if (i == last) output = output + ">";
					else output = output + " ";
				}
				output = output + "]";
				cout << output << endl;
			}
			T->GetEntry(j);
			int g = -1; // to skip dead time
			bool check = false;
			while(globalTime < Time){
				if(check) g++;
				if(Time - globalTime > 100 && g == -1){
					check = true;
				}
				if(g + 1 == int(deltaT/pitch/2)){
					globalTime = Time - 10;
					check = false;
					g = -1;
				}
				signal.erase(signal.begin());
				signalT.erase(signalT.begin());
				double amplitude = 0;
				int n = activationTime.size();
				int actual = 0;
				for(int k = 0; k < n; k++){
					amplitude += Signal(globalTime - activationTime.at(actual));
					if(globalTime - activationTime.at(actual) > 40){ 
						activationTime.erase(activationTime.begin() + actual);
						actual --;
					}
					actual ++;
				}
				signal.push_back(amplitude);
				globalTime += pitch;
				signalT.push_back(globalTime );
				if(signal.at(int(deltaT/pitch/2)) > threashold && thCheck == 0){
					thCheck = 1;
					timeCheck= signalT.at(int(deltaT/pitch/2));
					signalTime = globalTime - deltaT/2;
					int l = 0;
					double DeltaT = 0;
					while(true){
						if(signal.at(int(deltaT/pitch/2) + l) < threashold){
							DeltaT = signalT.at(int(deltaT/pitch/2) + l) - signalT.at(int(deltaT/pitch/2));
							break;
						}
						l++;
					}
					Amplitude = TMath::MaxElement(int(DeltaT/pitch), &signal.at(int(deltaT/pitch/2)));
					Twaves[Channel]->Fill();
				}
				else if(signal.at(int(deltaT/pitch/2)) < threashold && thCheck == 1){
// && signalT.at(int(deltaT/pitch/2)) - timeCheck> 10
					thCheck = 0;
				}
			}
			activationTime.push_back(Time);
		}

		for(int j = 0; j < int(deltaT/pitch/2); j++){
			signal.erase(signal.begin());
			signalT.erase(signalT.begin());
			double amplitude = 0;
			int n = activationTime.size();
			int actual = 0;
			for(int k = 0; k < n; k++){
				amplitude += Signal(globalTime - activationTime.at(actual));
				if(globalTime - activationTime.at(actual) > 40){ 
					activationTime.erase(activationTime.begin() + actual);
					actual --;
				}
				actual ++;
			}
			signal.push_back(amplitude);
			globalTime += pitch;
			signalT.push_back(globalTime );
			if(signal.at(int(deltaT/pitch/2)) > threashold && thCheck == 0){
				thCheck = 1;
				timeCheck= signalT.at(int(deltaT/pitch/2));
				signalTime = globalTime - deltaT/2;
				int l = 0;
				double DeltaT = 0;
				while(true){
					if(signal.at(int(deltaT/pitch/2) + l) < threashold){
						DeltaT = signalT.at(int(deltaT/pitch/2) + l) - signalT.at(int(deltaT/pitch/2));
						break;
					}
					l++;
				}
				Amplitude = TMath::MaxElement(int(DeltaT/pitch), &signal.at(int(deltaT/pitch/2)));
				Twaves[Channel]->Fill();
			}
			else if(signal.at(int(deltaT/pitch/2)) < threashold && thCheck == 1){
// && signalT.at(int(deltaT/pitch/2)) - timeCheck> 10
				thCheck = 0;
			}
		}
		signal.clear();
		signal.resize(int(deltaT / pitch) + 1);
		globalTime = 0;
		timeCheck= -30;
	}

	Twaves[0]->Write("waves0", TObject::kOverwrite);
	Twaves[1]->Write("waves1", TObject::kOverwrite);
	F->Close();
}

void reordering(TString file){
	TFile* F = TFile::Open(file + ".root", "UPDATE");
	TTree* T[2];
	T[0] = (TTree*) F->Get("waves0");
	T[1] = (TTree*) F->Get("waves1");

	int Channel = 0;
	int Fiber = 0;
	double Amplitude = 0;
	double Time = 0;

	int Channel0 = 0;
	int Fiber0 = 0;
	double Amplitude0 = 0;
	double Time0 = 0;

	int Channel1 = 0;
	int Fiber1 = 0;
	double Amplitude1 = 0;
	double Time1 = 0;

	T[0]->SetBranchAddress("Channel", &Channel0);
	T[0]->SetBranchAddress("Fiber", &Fiber0);
	T[0]->SetBranchAddress("Amplitude", &Amplitude0);
	T[0]->SetBranchAddress("Time", &Time0);

	T[1]->SetBranchAddress("Channel", &Channel1);
	T[1]->SetBranchAddress("Fiber", &Fiber1);
	T[1]->SetBranchAddress("Amplitude", &Amplitude1);
	T[1]->SetBranchAddress("Time", &Time1);

	TTree* Tnew = new TTree("new", "reordered");
	Tnew->Branch("Channel", &Channel);
	Tnew->Branch("Fiber", &Fiber);
	Tnew->Branch("Amplitude", &Amplitude);
	Tnew->Branch("Time", &Time);

	int i0 = 0, i1 = 0;
	while(true){
		if(T[0]->GetEntry(i0) == 0 && T[1]->GetEntry(i1) == 0) break;

		else if(T[0]->GetEntry(i0) == 0){
			T[1]->GetEntry(i1);
			Channel = Channel1;
			Fiber = Fiber1;
			Amplitude = Amplitude1;
			Time = Time1;
			Tnew->Fill();
			i1++;
		}

		else if(T[1]->GetEntry(i1) == 0){
			T[0]->GetEntry(i0);
			Channel = Channel0;
			Fiber = Fiber0;
			Amplitude = Amplitude0;
			Time = Time0;
			Tnew->Fill();
			i0++;
		}

		else{
			T[0]->GetEntry(i0);
			T[1]->GetEntry(i1);

			if(Time0 < Time1){
				Channel = Channel0;
				Fiber = Fiber0;
				Amplitude = Amplitude0;
				Time = Time0;
				Tnew->Fill();
				i0++;
			}

			else{
				Channel = Channel1;
				Fiber = Fiber1;
				Amplitude = Amplitude1;
				Time = Time1;
				Tnew->Fill();
				i1++;
			}
		}
	}

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


