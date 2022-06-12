#include "header.h"

using namespace std;

auto c1 = new TCanvas("c1", "Histogram", 1280, 1080);
TH2* h1 = new TH2F("h1", "Histogram (W,Q^{2})_2d", 9600, 0.0, 40.0, 3200, 0.0, 40.0);
TH1F* h3 = new TH1F("h3", "Histogram W", 9600, 0.0, 40.0);
TH1F* h4 = new TH1F("h4", "Histogram Q^{2}", 3200, 0.0, 40.0);
TH2* h5 = new TH2F("h5", "Histogram (#phi,cos(#theta^{*}))", 180, 0, 2*M_PI, 100, -1 , 1);

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();
    auto finish_ = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ = (finish_ - start);

    int counter_(500), inter(1), counter_events(0), limit(0);

    input_check(argc, argv);
    if(seed_ == 0){seed_ = time(NULL);}
    srand(seed_);
    if(rad_corr){counter_ = 5; inter = 3;}

    ofstream File;
    string PATH;
    int id_dot = path.find_last_of('.');

    if(truncate_out)
    {
        for(int counter = 0; counter < ceil(double(N)/10000.0); counter++)
        {
            limit = (counter == ceil(double(N)/10000.0) - 1) ? N - counter*10000 : 10000;

            for(int k = 0; k < limit; k += inter)
            {
                generate_particle(k);
            }

            PATH = path;
            PATH.insert(id_dot,std::to_string(counter));

            File.open(PATH,fstream::in | fstream::out | fstream::app);

            for(auto i:LUND_OUTPUT)
            {
                for(auto j:i)
                {
                    File << j << "\t";
                }
                File << endl;
            }

            finish_ = std::chrono::high_resolution_clock::now();
            elapsed_ = double(ceil(double(N)/10000.0) - 1 - counter)*(finish_ - start)/double(counter);
            cout << "Stand by ..." << ceil(counter*10000.0/(ceil(double(N)/10000.0) - 1))/100.0 << "%   Time remaining: " << floor(elapsed_.count()/3600) << " h " << floor((elapsed_.count() - 3600*floor(elapsed_.count()/3600))/60) << " min " << elapsed_.count() - 60*floor(elapsed_.count()/60) << " s             \r" << flush;

            File.close(); LUND_OUTPUT.clear();
        }
    }else
    {
        for(int k = 0; k < N; k += inter)
        {
            generate_particle(k);
            if(k % counter_ == 0)
            {
                finish_ = std::chrono::high_resolution_clock::now();
                elapsed_ = double(N - k)*(finish_ - start)/double(k);
                if(!method) cout << "Acceptance Rate: " << floor(10000*(double(Accepted_Number)/double(Total_Number)))/100 << "%\tStand by ..." << k*100/N << "%   Time remaining: " << floor(elapsed_.count()/3600) << " h " << floor((elapsed_.count() - 3600*floor(elapsed_.count()/3600))/60) << " min " << elapsed_.count() - 60*floor(elapsed_.count()/60) << " s             \r" << flush;
                else cout << "Stand by ..." << k*100/N << "%   Time remaining: " << floor(elapsed_.count()/3600) << " h " << floor((elapsed_.count() - 3600*floor(elapsed_.count()/3600))/60) << " min " << elapsed_.count() - 60*floor(elapsed_.count()/60) << " s             \r" << flush;
            }
        }

        File.open(path,fstream::in | fstream::out | fstream::app);

        for(int i = 0; i < LUND_OUTPUT.size(); i++)
        {
            for(auto j:LUND_OUTPUT[i])
            {
                File << j << "\t";
            }
            File << endl;

            if(i % counter_ == 0) cout << "Recording LUND ... " << ceil((i+1)*10000.0/LUND_OUTPUT.size())/100.0 << "%                                                                                  \r" << flush;
        }

        File.close();
    }

    cout << "Recording complete                                            " << endl;

    if(histogram)
    {
        c1 -> Divide(2 , 2);

        h1->GetYaxis()->SetTitle("Q^{2}, GeV^{2}");
        h1->GetYaxis()->SetTitleOffset(1.3);
        h1->GetYaxis()->CenterTitle(true);
        h1->GetXaxis()->SetTitle("W, GeV");
        h1->GetXaxis()->CenterTitle(true);

        c1->cd(1);
        c1->cd(1)->SetLogz();
        h1->SetAxisRange(W_min-0.1, W_max+0.1, "X");
        h1->SetAxisRange(Q2_min-0.5, Q2_max+0.5, "Y");
        h1->Draw("COL");

        c1->cd(2);
        c1->cd(2)->SetLogz();

        h5->GetYaxis()->SetTitle("cos(#theta^{*})");
        h5->GetYaxis()->SetTitleOffset(1.3);
        h5->GetYaxis()->CenterTitle(true);
        h5->GetXaxis()->SetTitle("#phi");
        h5->GetXaxis()->CenterTitle(true);

        h5->Draw("COL");

        h3->GetXaxis()->SetTitle("W, GeV");
        h3->GetXaxis()->CenterTitle(true);

        c1->cd(3);
        h3->SetAxisRange(W_min-0.1, W_max+0.1);
        h3->Draw();

        h4->GetXaxis()->SetTitle("Q^{2}, GeV^{2}");
        h4->GetXaxis()->CenterTitle(true);

        c1->cd(4);
        h4->SetAxisRange(Q2_min-0.5, Q2_max+0.5);
        c1->cd(4)->SetLogy();
        h4->Draw();

        c1 -> Print("Histogram.jpeg");
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Elapsed time: " << floor(elapsed.count()/3600) << " h " << floor((elapsed.count() - 3600*floor(elapsed.count()/3600))/60) << " min " << elapsed.count() - 60*floor(elapsed.count()/60) << " s\n";

    data.clear(); data_interp.clear(); values_rad.clear(), LUND_OUTPUT.clear();

    return 0;
}
