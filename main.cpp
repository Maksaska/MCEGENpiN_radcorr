#include "header.h"

using namespace std;

auto c1 = new TCanvas("c1", "Histogram", 1280, 1080);   
TH2* h1 = new TH2F("h1", "Histogram (W,Q^{2})_2d", 186, 1.08, 2, 440, 0, 12);
TH1F* h3 = new TH1F("h3", "Histogram W", 186, 1.08, 2);
TH1F* h4 = new TH1F("h4", "Histogram Q^{2}", 440, 0, 12);
TH2* h5 = new TH2F("h5", "Histogram (#phi,cos(#theta^{*}))", 180, 0, 2*M_PI, 100, -1 , 1);

int main(int argc, char* argv[])
{	
	auto start = std::chrono::high_resolution_clock::now();
	int counter_(500), inter(1); 

	input_check(argc, argv); 
	if(seed_ == 0){seed_ = time(NULL);}
	srand(seed_);
	if(rad_corr){counter_ = 5; inter = 3;}
	
	for(int k = 0; k < N; k += inter)
	{
		generate_particle(k);
		if(k % counter_ == 0)
		{
			auto finish_ = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_ = double(N - k)*(finish_ - start)/double(k); 
			cout << "Stand by ... " << k*100/N << "%   Time remaining: " << floor(elapsed_.count()/3600) << " h " << floor((elapsed_.count() - 3600*floor(elapsed_.count()/3600))/60) << " min " << elapsed_.count() - 60*floor(elapsed_.count()/60) << " s             \r" << flush;
		}
	}

	if(histogram)
	{
		c1 -> Divide(2 , 2);
		
		c1->cd(1);
		c1->cd(1)->SetLogz();
		h1->Draw("COL");

		c1->cd(2);
		c1->cd(2)->SetLogz();
		
		h5->GetYaxis()->SetTitle("cos(#theta^{*})");
		h5->GetYaxis()->SetTitleOffset(1.3);
		h5->GetYaxis()->CenterTitle(true);
		h5->GetXaxis()->SetTitle("#phi , grad");
		h5->GetXaxis()->CenterTitle(true);		

		h5->Draw("COL");

		c1->cd(3);
		h3->Draw();

		c1->cd(4);
		c1->cd(4)->SetLogy();
		h4->Draw();

		c1 -> Print("Histogram.jpeg");	
	}	
	
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Elapsed time: " << floor(elapsed.count()/3600) << " h " << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min " << elapsed.count() - 60*floor(elapsed.count()/60) << " s\n";
	
	return 0;
}
