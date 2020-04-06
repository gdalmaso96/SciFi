#include "./signalslite.C"
#include "./readgio.cpp"

void scores(){
	vector<std::string> name;
	name.push_back("scores/025mmpos50polish10e6");
	name.push_back("scores/025mmpos50polish56e6");
	name.push_back("scores/025mmpos50polish31e7");
	name.push_back("scores/025mmpos50polish18e8");
	name.push_back("scores/025mmpos50polish10e9");

	name.push_back("scores/025mmpos50ground10e6");
	name.push_back("scores/025mmpos50ground56e6");
	name.push_back("scores/025mmpos50ground31e7");
	name.push_back("scores/025mmpos50ground18e8");
	name.push_back("scores/025mmpos50ground10e9");

	name.push_back("scores/025mmpos25polish10e6");
	name.push_back("scores/025mmpos25polish56e6");
	name.push_back("scores/025mmpos25polish31e7");
	name.push_back("scores/025mmpos25polish18e8");
	name.push_back("scores/025mmpos25polish10e9");

	name.push_back("scores/025mmpos25ground10e6");
	name.push_back("scores/025mmpos25ground56e6");
	name.push_back("scores/025mmpos25ground31e7");
	name.push_back("scores/025mmpos25ground18e8");
	name.push_back("scores/025mmpos25ground10e9");

	vector<double> Rate;
	Rate.push_back(1.0e6);
	Rate.push_back(5.6e6);
	Rate.push_back(3.1e7);
	Rate.push_back(1.8e8);
	Rate.push_back(1.0e9);

	TFile* N = TFile::Open("results/results.root", "UPDATE");
	vector<TString> name_graph;
	name_graph.push_back("025mmpos50polish");
	name_graph.push_back("025mmpos50ground");
	name_graph.push_back("025mmpos25polish");
	name_graph.push_back("025mmpos25ground");

	for(int i = 0; i < 4; i++){
		TGraphErrors* G = new TGraphErrors(5);

		for(int j = 0; j < 5; j++){
			signals(name.at(i*5 + j) + ".root", name.at(i*5 + j) + "out");
			processing(2, name.at(i*5 + j) + "out");
			
			double R[2] = {0};
			readSciFi(name.at(i*5 + j) + "out.txt", R);

			G->SetPoint(j, Rate.at(j), R[0]/Rate.at(j));
		}

		N->cd();
		G->Write(name_graph.at(i), TObject::kOverwrite);
	}

	N->Close();
}
