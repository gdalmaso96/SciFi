void rateCheck(TString name){
	TFile* F = TFile::Open(name + ".root");
	TTree* T = (TTree*) F->Get("T");

	int NgammaOut = 0;
	double GunTime = 0;
	vector<int> *Channel = 0;

	T->SetBranchAddress("NgammaOut", &NgammaOut);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("GunTime", &GunTime);

	int Counts[42] = {0};
	int n_entries = T->GetEntries();

	for(int i = 0; i < n_entries; i++){
		T->GetEntry(i);

		if(NgammaOut > 1 && Channel->size() > 0) Counts[Channel->at(0)/2]++;
	}

	ofstream myfile;
	myfile.open(name + "check.txt");
	for(int i = 0; i < 42; i++){
		myfile << Counts[i] << endl;
	}
	myfile << 0.22 << endl;
	myfile << GunTime << endl;
	myfile.close();
}
